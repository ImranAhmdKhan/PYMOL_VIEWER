#!/usr/bin/env python3
"""
PyMOL Publication-Quality Protein–Ligand Complex Viewer
========================================================
Renders multiple protein–ligand complexes (PDB / MAE / MAEGZ) as a
publication-ready subplot grid image using PyMOL's headless API.

Requirements
------------
    pip install PyQt5 pymol-open-source Pillow matplotlib

Usage
-----
    python pymol_viewer.py
"""

import os
import sys
import gzip
import math
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

# ---------------------------------------------------------------------------
# PyQt5 imports
# ---------------------------------------------------------------------------
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QPushButton, QLabel, QComboBox, QSpinBox, QDoubleSpinBox,
    QFileDialog, QListWidget, QListWidgetItem, QGroupBox, QCheckBox,
    QTabWidget, QScrollArea, QFrame, QSplitter, QMessageBox,
    QColorDialog, QProgressBar, QStatusBar, QLineEdit, QSlider,
    QAbstractItemView, QSizePolicy, QFormLayout,
    QDialog, QDialogButtonBox, QAction, QMenuBar,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QSize, QSettings
from PyQt5.QtGui import QPixmap, QImage, QColor, QFont, QIcon

# ---------------------------------------------------------------------------
# Third-party imaging
# ---------------------------------------------------------------------------
try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:
    raise SystemExit("Pillow is required: pip install Pillow")

# ---------------------------------------------------------------------------
# Optional matplotlib for subplot compositing (fallback if PIL fonts fail)
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch


# ===========================================================================
# Data model
# ===========================================================================

PROTEIN_REPRS = ["cartoon", "ribbon", "surface", "lines", "sticks", "sphere"]
LIGAND_REPRS = ["sticks", "ball_and_stick", "sphere", "lines"]
COLOR_SCHEMES = [
    "chain",      # util.cbc – different color per chain
    "element",    # CPK coloring
    "spectrum",   # rainbow along residue index
    "bfactor",    # blue→white→red by B-factor
    "white",
    "cyan",
    "green",
    "yellow",
    "orange",
    "red",
    "magenta",
    "blue",
    "grey",
]
RAY_MODES = {
    "Default (Colour + Outline)": 1,
    "Cartoon (No Outline)": 0,
    "Black Outline": 2,
    "Quantised Colour": 3,
}
EXPORT_FORMATS = ["PNG", "TIFF", "PDF"]

# ---------------------------------------------------------------------------
# App-level settings keys
# ---------------------------------------------------------------------------
SETTINGS_ORG  = "PyMOLViewer"
SETTINGS_APP  = "PyMOLViewer"
KEY_PYMOL_PATH = "pymol_scripts_path"
KEY_VMD_PATH   = "vmd_exe_path"

_DEFAULT_PYMOL_PATH = r"C:\Users\<username>\AppData\Local\Schrodinger\PyMOL2\Scripts"
_DEFAULT_VMD_PATH   = r"C:\Program Files\University of Illinois\VMD\vmd.exe"


def _load_app_settings() -> QSettings:
    return QSettings(SETTINGS_ORG, SETTINGS_APP)


def _setup_pymol_path():
    """Inject the configured PyMOL scripts directory into sys.path so that
    ``import pymol2`` resolves to the Schrödinger (or any custom) build."""
    s = _load_app_settings()
    path = s.value(KEY_PYMOL_PATH, _DEFAULT_PYMOL_PATH)
    if path and os.path.isdir(path) and path not in sys.path:
        sys.path.insert(0, path)
        # Also add the parent directory (sometimes needed for Schrödinger)
        parent = str(Path(path).parent)
        if parent not in sys.path:
            sys.path.insert(1, parent)


def _resolve_vmd_exe() -> str:
    """Return the VMD executable path from settings.

    On Windows the user may have supplied a .lnk shortcut; we fall back to
    resolving via the shell (``start`` command) in that case.
    """
    s = _load_app_settings()
    return s.value(KEY_VMD_PATH, _DEFAULT_VMD_PATH)


class ComplexEntry:
    """All settings for one protein–ligand complex."""

    def __init__(self, filepath: str):
        self.filepath: str = filepath
        self.filename: str = Path(filepath).name
        self.name: str = Path(filepath).stem.replace(" ", "_").replace("-", "_")
        self.label: str = Path(filepath).stem          # subplot caption

        # Representations
        self.protein_repr: str = "cartoon"
        self.ligand_repr: str = "sticks"
        self.ligand_selection: str = "organic"         # PyMOL selection

        # Colours
        self.protein_color: str = "chain"
        self.ligand_color: str = "element"

        # Surface
        self.show_surface: bool = False
        self.surface_transparency: float = 0.4

        # H-bonds
        self.show_hbonds: bool = True

        # Ray settings
        self.bg_color: str = "white"
        self.ray_shadows: bool = True

        # Docking analysis / annotation
        self.show_residue_labels: bool = True   # label residues within radius
        self.label_radius: float = 5.0          # Å from ligand to show labels
        self.label_color: str = "black"         # label text colour in PyMOL
        self.show_distances: bool = False        # draw extra distance lines
        self.distance_cutoff: float = 4.0       # Å cutoff for contact distances
        self.zoom_extra: float = 0.0            # extra zoom-out (0 = default)


# ===========================================================================
# Rendering worker (runs in a background QThread)
# ===========================================================================

class RenderWorker(QThread):
    progress = pyqtSignal(int)      # 0-100
    status   = pyqtSignal(str)
    finished = pyqtSignal(str)      # path to output file
    error    = pyqtSignal(str)

    def __init__(
        self,
        complexes: List[ComplexEntry],
        layout: Tuple[int, int],
        output_path: str,
        settings: dict,
        parent=None,
    ):
        super().__init__(parent)
        self.complexes    = complexes
        self.layout       = layout          # (rows, cols)
        self.output_path  = output_path
        self.settings     = settings

    # ------------------------------------------------------------------
    def run(self):
        try:
            self._render()
        except Exception as exc:
            import traceback
            self.error.emit(traceback.format_exc())

    # ------------------------------------------------------------------
    def _render(self):
        _setup_pymol_path()   # ensure configured PyMOL is on sys.path
        try:
            import pymol2
        except ImportError:
            self.error.emit(
                "PyMOL is not installed.\n"
                "Install it with:  pip install pymol-open-source\n"
                "or from the official Schrödinger bundle."
            )
            return

        rows, cols     = self.layout
        dpi            = self.settings.get("dpi", 300)
        img_w          = self.settings.get("img_width", 800)
        img_h          = self.settings.get("img_height", 600)
        ray_mode       = self.settings.get("ray_mode", 1)
        fast_preview   = self.settings.get("fast_preview", False)

        rendered: List[Tuple[Image.Image, str]] = []

        for idx, cx in enumerate(self.complexes):
            self.status.emit(f"Rendering [{idx + 1}/{len(self.complexes)}]  {cx.filename} …")
            self.progress.emit(int(idx / len(self.complexes) * 80))

            with pymol2.PyMOL() as pymol:
                cmd = pymol.cmd

                # ── global quality settings ──────────────────────────────
                cmd.set("antialias",             2 if not fast_preview else 1)
                cmd.set("ray_shadows",           "1" if cx.ray_shadows else "0")
                cmd.set("ambient",               0.35)
                cmd.set("specular",              0.6)
                cmd.set("shininess",             50)
                cmd.set("ray_trace_mode",        0 if fast_preview else ray_mode)
                cmd.set("ray_opaque_background", 1)
                cmd.set("cartoon_fancy_helices", 1)
                cmd.set("cartoon_fancy_sheets",  1)
                cmd.set("cartoon_loop_radius",   0.3)
                cmd.set("stick_radius",          0.14)
                cmd.set("stick_h_scale",         1.0)
                cmd.set("line_width",            2.0)
                cmd.set("sphere_scale",          0.25)
                cmd.set("depth_cue",             0 if fast_preview else 1)

                # ── load file ────────────────────────────────────────────
                filepath = cx.filepath
                # MAEGZ: decompress to a temp .mae so PyMOL can load it
                _tmp_mae = None
                if filepath.lower().endswith(".maegz"):
                    _tmp_mae = tempfile.NamedTemporaryFile(
                        suffix=".mae", delete=False
                    )
                    with gzip.open(filepath, "rb") as gz:
                        shutil.copyfileobj(gz, _tmp_mae)
                    _tmp_mae.close()
                    filepath = _tmp_mae.name

                try:
                    cmd.load(filepath, cx.name)
                except Exception as exc:
                    self.error.emit(f"Cannot load '{cx.filename}': {exc}")
                    return
                finally:
                    if _tmp_mae and os.path.exists(_tmp_mae.name):
                        os.unlink(_tmp_mae.name)

                # ── remove solvent / ions ─────────────────────────────────
                cmd.remove("solvent")
                cmd.remove("inorganic and not metals")

                # ── selections ──────────────────────────────────────────
                protein_sel = f"({cx.name}) and polymer"
                ligand_sel  = f"({cx.name}) and ({cx.ligand_selection})"

                cmd.hide("everything", cx.name)
                cmd.bg_color(cx.bg_color)

                # ── protein ─────────────────────────────────────────────
                cmd.show(cx.protein_repr, protein_sel)
                self._apply_color(cmd, cx.protein_color, protein_sel)

                # ── ligand ──────────────────────────────────────────────
                if cx.ligand_repr == "ball_and_stick":
                    cmd.show("sticks", ligand_sel)
                    cmd.show("spheres", ligand_sel)
                    cmd.set("sphere_scale", 0.25, ligand_sel)
                else:
                    cmd.show(cx.ligand_repr, ligand_sel)
                self._apply_color(cmd, cx.ligand_color, ligand_sel)

                # ── surface ─────────────────────────────────────────────
                if cx.show_surface:
                    cmd.show("surface", protein_sel)
                    cmd.set("transparency", cx.surface_transparency, protein_sel)

                # ── H-bonds ─────────────────────────────────────────────
                if cx.show_hbonds:
                    try:
                        hb_name = f"hbonds_{idx}"
                        cmd.dist(hb_name, protein_sel, ligand_sel,
                                 mode=2, cutoff=3.5, label=0)
                        cmd.hide("labels", hb_name)
                        cmd.color("yellow", hb_name)
                        cmd.set("dash_width",  2.5, hb_name)
                        cmd.set("dash_length", 0.3, hb_name)
                        cmd.set("dash_gap",    0.3, hb_name)
                    except Exception:
                        pass

                # ── contact distances (docking analysis) ─────────────────
                if cx.show_distances:
                    try:
                        dist_name = f"contacts_{idx}"
                        cmd.dist(dist_name, ligand_sel, protein_sel,
                                 mode=0, cutoff=cx.distance_cutoff, label=1)
                        cmd.color("cyan", dist_name)
                        cmd.set("dash_width",  1.5, dist_name)
                        cmd.set("label_size",  10,  dist_name)
                        cmd.set("label_color", "cyan", dist_name)
                    except Exception:
                        pass

                # ── residue labels around ligand ─────────────────────────
                if cx.show_residue_labels:
                    try:
                        nearby_sel = (
                            f"({cx.name}) and polymer and "
                            f"(byres (all within {cx.label_radius} of ({ligand_sel})))"
                        )
                        cmd.set("label_size",  12)
                        cmd.set("label_color", cx.label_color)
                        cmd.set("label_font_id", 7)   # bold
                        cmd.label(
                            f"{nearby_sel} and name CA",
                            r'"%s%s" % (resn, resi)',
                        )
                    except Exception:
                        pass

                # ── camera – zoom on ligand ──────────────────────────────
                try:
                    cmd.orient(ligand_sel)
                    zoom_buf = 6 + cx.zoom_extra
                    cmd.zoom(ligand_sel, zoom_buf)
                except Exception:
                    cmd.zoom(cx.name)

                # ── render ───────────────────────────────────────────────
                with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tf:
                    tmp_png = tf.name

                if fast_preview:
                    cmd.png(tmp_png, width=img_w, height=img_h, dpi=72, ray=0)
                else:
                    cmd.ray(img_w, img_h)
                    cmd.png(tmp_png, dpi=dpi)

                img = Image.open(tmp_png).copy()
                os.unlink(tmp_png)

            rendered.append((img, cx.label))

        # ── composite ────────────────────────────────────────────────────
        self.status.emit("Compositing subplot image …")
        self.progress.emit(88)

        composite = self._compose_subplot(
            rendered,
            rows, cols,
            dpi           = dpi,
            title         = self.settings.get("title", ""),
            label_size    = self.settings.get("label_size", 28),
            caption_size  = self.settings.get("caption_size", 22),
            bg_color      = self.settings.get("subplot_bg", "white"),
            border        = self.settings.get("border_px", 30),
            show_panels   = self.settings.get("show_panel_letters", True),
        )

        fmt = Path(self.output_path).suffix.lower().lstrip(".")
        if fmt == "pdf":
            # PIL cannot save multi-page; use matplotlib for PDF
            self._save_as_pdf(composite, self.output_path, dpi)
        else:
            composite.save(self.output_path, dpi=(dpi, dpi))

        self.progress.emit(100)
        self.finished.emit(self.output_path)

    # ------------------------------------------------------------------
    @staticmethod
    def _apply_color(cmd, scheme: str, selection: str):
        if scheme == "element":
            cmd.color("gray80",  f"{selection} and elem C")
            cmd.color("red",     f"{selection} and elem O")
            cmd.color("blue",    f"{selection} and elem N")
            cmd.color("yellow",  f"{selection} and elem S")
            cmd.color("white",   f"{selection} and elem H")
        elif scheme == "chain":
            cmd.util.cbc(selection)
        elif scheme == "spectrum":
            cmd.spectrum("resi", "rainbow", selection)
        elif scheme == "bfactor":
            cmd.spectrum("b", "blue_white_red", selection)
        else:
            cmd.color(scheme, selection)

    # ------------------------------------------------------------------
    @staticmethod
    def _load_font(size: int, bold: bool = False):
        candidates = [
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
            "/Library/Fonts/Arial Bold.ttf",
            "C:/Windows/Fonts/arialbd.ttf",
        ]
        for path in candidates:
            if os.path.exists(path):
                try:
                    from PIL import ImageFont
                    return ImageFont.truetype(path, size)
                except Exception:
                    pass
        from PIL import ImageFont
        return ImageFont.load_default()

    # ------------------------------------------------------------------
    def _compose_subplot(
        self,
        images_labels: List[Tuple[Image.Image, str]],
        rows: int, cols: int,
        dpi: int = 300,
        title: str = "",
        label_size: int = 28,
        caption_size: int = 22,
        bg_color: str = "white",
        border: int = 30,
        show_panels: bool = True,
    ) -> Image.Image:

        if not images_labels:
            return Image.new("RGB", (400, 400), bg_color)

        iw, ih = images_labels[0][0].size

        font_panel   = self._load_font(label_size, bold=True)
        font_caption = self._load_font(caption_size)
        font_title   = self._load_font(int(label_size * 1.4), bold=True)

        # Estimate text heights
        caption_h = (caption_size + 10) if any(lbl for _, lbl in images_labels) else 0
        title_h   = (int(label_size * 1.4) + 20) if title else 0

        total_w = cols * iw + (cols + 1) * border
        total_h = rows * (ih + caption_h) + (rows + 1) * border + title_h

        canvas = Image.new("RGB", (total_w, total_h), bg_color)
        draw   = ImageDraw.Draw(canvas)

        # ── title ──────────────────────────────────────────────────────
        if title:
            draw.text(
                (total_w // 2, border // 2),
                title, fill="black", font=font_title, anchor="mt",
            )

        panel_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        for idx, (img, caption) in enumerate(images_labels):
            if idx >= rows * cols:
                break
            row = idx // cols
            col = idx % cols

            x = border + col * (iw + border)
            y = title_h + border + row * (ih + caption_h + border)

            canvas.paste(img, (x, y))

            # Panel letter (A, B, C …)
            if show_panels:
                letter = panel_letters[idx] if idx < len(panel_letters) else str(idx + 1)
                # Shadow for readability on any background
                draw.text((x + 12, y + 8),  letter, fill="black",  font=font_panel)
                draw.text((x + 10, y + 6),  letter, fill="white",  font=font_panel)

            # Caption below image
            if caption and caption_h:
                draw.text(
                    (x + iw // 2, y + ih + 4),
                    caption, fill="black", font=font_caption, anchor="mt",
                )

        return canvas

    # ------------------------------------------------------------------
    @staticmethod
    def _save_as_pdf(img: Image.Image, path: str, dpi: int):
        fig, ax = plt.subplots(
            figsize=(img.width / dpi, img.height / dpi), dpi=dpi
        )
        ax.imshow(img)
        ax.axis("off")
        fig.tight_layout(pad=0)
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0)
        plt.close(fig)


# ===========================================================================
# Per-complex settings panel
# ===========================================================================

class ComplexSettingsPanel(QWidget):
    """Edit settings for one ComplexEntry."""

    changed = pyqtSignal()

    def __init__(self, entry: ComplexEntry, parent=None):
        super().__init__(parent)
        self.entry = entry
        self._build_ui()

    def _build_ui(self):
        layout = QFormLayout(self)
        layout.setLabelAlignment(Qt.AlignRight)

        # Label / caption
        self.w_label = QLineEdit(self.entry.label)
        self.w_label.textChanged.connect(self._sync)
        layout.addRow("Caption:", self.w_label)

        # Ligand selection
        self.w_ligand_sel = QLineEdit(self.entry.ligand_selection)
        self.w_ligand_sel.setToolTip(
            "PyMOL selection for the ligand, e.g. 'organic', 'resn LIG', '/chain A//"
        )
        self.w_ligand_sel.textChanged.connect(self._sync)
        layout.addRow("Ligand selection:", self.w_ligand_sel)

        # Protein representation
        self.w_prot_repr = QComboBox()
        self.w_prot_repr.addItems(PROTEIN_REPRS)
        self.w_prot_repr.setCurrentText(self.entry.protein_repr)
        self.w_prot_repr.currentTextChanged.connect(self._sync)
        layout.addRow("Protein repr:", self.w_prot_repr)

        # Ligand representation
        self.w_lig_repr = QComboBox()
        self.w_lig_repr.addItems(LIGAND_REPRS)
        self.w_lig_repr.setCurrentText(self.entry.ligand_repr)
        self.w_lig_repr.currentTextChanged.connect(self._sync)
        layout.addRow("Ligand repr:", self.w_lig_repr)

        # Protein colour
        self.w_prot_color = QComboBox()
        self.w_prot_color.addItems(COLOR_SCHEMES)
        self.w_prot_color.setCurrentText(self.entry.protein_color)
        self.w_prot_color.currentTextChanged.connect(self._sync)
        layout.addRow("Protein color:", self.w_prot_color)

        # Ligand color
        self.w_lig_color = QComboBox()
        self.w_lig_color.addItems(COLOR_SCHEMES)
        self.w_lig_color.setCurrentText(self.entry.ligand_color)
        self.w_lig_color.currentTextChanged.connect(self._sync)
        layout.addRow("Ligand color:", self.w_lig_color)

        # Surface
        self.w_surface = QCheckBox("Show protein surface")
        self.w_surface.setChecked(self.entry.show_surface)
        self.w_surface.toggled.connect(self._sync)
        layout.addRow("", self.w_surface)

        self.w_surf_trans = QDoubleSpinBox()
        self.w_surf_trans.setRange(0.0, 1.0)
        self.w_surf_trans.setSingleStep(0.05)
        self.w_surf_trans.setValue(self.entry.surface_transparency)
        self.w_surf_trans.valueChanged.connect(self._sync)
        layout.addRow("Surface transparency:", self.w_surf_trans)

        # H-bonds
        self.w_hbonds = QCheckBox("Show H-bonds")
        self.w_hbonds.setChecked(self.entry.show_hbonds)
        self.w_hbonds.toggled.connect(self._sync)
        layout.addRow("", self.w_hbonds)

        # Background colour
        self.w_bg = QComboBox()
        self.w_bg.addItems(["white", "black", "grey"])
        self.w_bg.setCurrentText(self.entry.bg_color)
        self.w_bg.currentTextChanged.connect(self._sync)
        layout.addRow("Background:", self.w_bg)

        # Ray shadows
        self.w_shadows = QCheckBox("Ray shadows")
        self.w_shadows.setChecked(self.entry.ray_shadows)
        self.w_shadows.toggled.connect(self._sync)
        layout.addRow("", self.w_shadows)

        # ── Docking analysis ────────────────────────────────────────────
        sep = QLabel("─── Docking Analysis ───")
        sep.setAlignment(Qt.AlignCenter)
        layout.addRow(sep)

        # Residue labels
        self.w_res_labels = QCheckBox("Label binding-site residues")
        self.w_res_labels.setChecked(self.entry.show_residue_labels)
        self.w_res_labels.toggled.connect(self._sync)
        layout.addRow("", self.w_res_labels)

        self.w_label_radius = QDoubleSpinBox()
        self.w_label_radius.setRange(1.0, 20.0)
        self.w_label_radius.setSingleStep(0.5)
        self.w_label_radius.setSuffix(" Å")
        self.w_label_radius.setValue(self.entry.label_radius)
        self.w_label_radius.valueChanged.connect(self._sync)
        layout.addRow("Label radius:", self.w_label_radius)

        self.w_label_color = QComboBox()
        self.w_label_color.addItems(["black", "white", "yellow", "cyan",
                                     "magenta", "green", "red", "blue"])
        self.w_label_color.setCurrentText(self.entry.label_color)
        self.w_label_color.currentTextChanged.connect(self._sync)
        layout.addRow("Label color:", self.w_label_color)

        # Contact distances
        self.w_distances = QCheckBox("Show contact distances")
        self.w_distances.setChecked(self.entry.show_distances)
        self.w_distances.toggled.connect(self._sync)
        layout.addRow("", self.w_distances)

        self.w_dist_cutoff = QDoubleSpinBox()
        self.w_dist_cutoff.setRange(1.0, 10.0)
        self.w_dist_cutoff.setSingleStep(0.5)
        self.w_dist_cutoff.setSuffix(" Å")
        self.w_dist_cutoff.setValue(self.entry.distance_cutoff)
        self.w_dist_cutoff.valueChanged.connect(self._sync)
        layout.addRow("Distance cutoff:", self.w_dist_cutoff)

        # Zoom extra
        self.w_zoom_extra = QDoubleSpinBox()
        self.w_zoom_extra.setRange(-5.0, 30.0)
        self.w_zoom_extra.setSingleStep(1.0)
        self.w_zoom_extra.setSuffix(" Å")
        self.w_zoom_extra.setValue(self.entry.zoom_extra)
        self.w_zoom_extra.valueChanged.connect(self._sync)
        layout.addRow("Extra zoom-out:", self.w_zoom_extra)

    def _sync(self):
        self.entry.label                 = self.w_label.text()
        self.entry.ligand_selection      = self.w_ligand_sel.text() or "organic"
        self.entry.protein_repr          = self.w_prot_repr.currentText()
        self.entry.ligand_repr           = self.w_lig_repr.currentText()
        self.entry.protein_color         = self.w_prot_color.currentText()
        self.entry.ligand_color          = self.w_lig_color.currentText()
        self.entry.show_surface          = self.w_surface.isChecked()
        self.entry.surface_transparency  = self.w_surf_trans.value()
        self.entry.show_hbonds           = self.w_hbonds.isChecked()
        self.entry.bg_color              = self.w_bg.currentText()
        self.entry.ray_shadows           = self.w_shadows.isChecked()
        self.entry.show_residue_labels   = self.w_res_labels.isChecked()
        self.entry.label_radius          = self.w_label_radius.value()
        self.entry.label_color           = self.w_label_color.currentText()
        self.entry.show_distances        = self.w_distances.isChecked()
        self.entry.distance_cutoff       = self.w_dist_cutoff.value()
        self.entry.zoom_extra            = self.w_zoom_extra.value()
        self.changed.emit()


# ===========================================================================
# Settings dialog
# ===========================================================================

class SettingsDialog(QDialog):
    """Configure paths to PyMOL (Schrödinger build) and VMD."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Application Settings")
        self.setMinimumWidth(560)

        s = _load_app_settings()

        layout = QVBoxLayout(self)

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignRight)

        # PyMOL scripts path
        pymol_row = QHBoxLayout()
        self.w_pymol = QLineEdit(s.value(KEY_PYMOL_PATH, _DEFAULT_PYMOL_PATH))
        self.w_pymol.setPlaceholderText(r"e.g. C:\…\PyMOL2\Scripts")
        btn_pymol = QPushButton("Browse …")
        btn_pymol.clicked.connect(self._browse_pymol)
        pymol_row.addWidget(self.w_pymol)
        pymol_row.addWidget(btn_pymol)
        form.addRow("PyMOL Scripts path:", pymol_row)

        pymol_hint = QLabel(
            "Directory that contains pymol2.py / pymol package.\n"
            r"Example: C:\Users\<username>\AppData\Local\Schrodinger\PyMOL2\Scripts"
        )
        pymol_hint.setStyleSheet("color: #aaa; font-size: 10px;")
        form.addRow("", pymol_hint)

        # VMD executable path
        vmd_row = QHBoxLayout()
        self.w_vmd = QLineEdit(s.value(KEY_VMD_PATH, _DEFAULT_VMD_PATH))
        self.w_vmd.setPlaceholderText(r'e.g. C:\Program Files\University of Illinois\VMD\vmd.exe')
        btn_vmd = QPushButton("Browse …")
        btn_vmd.clicked.connect(self._browse_vmd)
        vmd_row.addWidget(self.w_vmd)
        vmd_row.addWidget(btn_vmd)
        form.addRow("VMD executable:", vmd_row)

        vmd_hint = QLabel(
            "Full path to vmd.exe or a .lnk shortcut.\n"
            r"Example: C:\Users\<username>\Desktop\VMD 2.0.lnk"
        )
        vmd_hint.setStyleSheet("color: #aaa; font-size: 10px;")
        form.addRow("", vmd_hint)

        layout.addLayout(form)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._save_and_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _browse_pymol(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select PyMOL Scripts directory",
            self.w_pymol.text() or "C:\\"
        )
        if path:
            self.w_pymol.setText(path)

    def _browse_vmd(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select VMD executable",
            self.w_vmd.text() or "C:\\",
            "Executables / Shortcuts (*.exe *.lnk);;All (*)"
        )
        if path:
            self.w_vmd.setText(path)

    def _save_and_accept(self):
        s = _load_app_settings()
        s.setValue(KEY_PYMOL_PATH, self.w_pymol.text().strip())
        s.setValue(KEY_VMD_PATH,   self.w_vmd.text().strip())
        self.accept()


# ===========================================================================
# Main window
# ===========================================================================

class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("PyMOL / VMD – Publication-Quality Docking Viewer")
        self.resize(1280, 820)

        # Inject configured PyMOL path before any pymol2 import attempt
        _setup_pymol_path()

        self.complexes: List[ComplexEntry] = []
        self._worker: Optional[RenderWorker] = None

        self._build_menu()
        self._build_ui()
        self._update_button_states()

    # -----------------------------------------------------------------------
    def _build_menu(self):
        menubar = self.menuBar()

        file_menu = menubar.addMenu("File")
        act_add = QAction("Add Files …", self)
        act_add.setShortcut("Ctrl+O")
        act_add.triggered.connect(self._add_files)
        file_menu.addAction(act_add)

        act_clear = QAction("Clear All", self)
        act_clear.triggered.connect(self._clear_all)
        file_menu.addAction(act_clear)

        file_menu.addSeparator()
        act_quit = QAction("Quit", self)
        act_quit.setShortcut("Ctrl+Q")
        act_quit.triggered.connect(self.close)
        file_menu.addAction(act_quit)

        tools_menu = menubar.addMenu("Tools")
        act_open_pymol = QAction("Open Selected in PyMOL (Interactive) …", self)
        act_open_pymol.triggered.connect(self._open_in_pymol_interactive)
        tools_menu.addAction(act_open_pymol)

        act_open_vmd = QAction("Open Selected in VMD …", self)
        act_open_vmd.triggered.connect(self._open_in_vmd)
        tools_menu.addAction(act_open_vmd)

        tools_menu.addSeparator()
        act_settings = QAction("Settings …", self)
        act_settings.setShortcut("Ctrl+,")
        act_settings.triggered.connect(self._open_settings)
        tools_menu.addAction(act_settings)

    # -----------------------------------------------------------------------
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(6)

        # ── left panel ──────────────────────────────────────────────────
        left = QWidget()
        left.setFixedWidth(420)
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)

        left_layout.addWidget(self._build_file_panel())
        left_layout.addWidget(self._build_tabs())

        root.addWidget(left)

        # ── right panel (preview + progress) ────────────────────────────
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)

        right_layout.addWidget(self._build_preview_area())
        right_layout.addWidget(self._build_action_bar())

        root.addWidget(right, stretch=1)

        # status bar
        self.statusBar().showMessage("Ready")

    # ------------------------------------------------------------------
    def _build_file_panel(self) -> QGroupBox:
        grp = QGroupBox("Structures")
        lay = QVBoxLayout(grp)

        btn_row = QHBoxLayout()
        self.btn_add = QPushButton("➕ Add Files …")
        self.btn_add.clicked.connect(self._add_files)
        self.btn_remove = QPushButton("➖ Remove")
        self.btn_remove.clicked.connect(self._remove_selected)
        self.btn_clear = QPushButton("🗑 Clear All")
        self.btn_clear.clicked.connect(self._clear_all)
        btn_row.addWidget(self.btn_add)
        btn_row.addWidget(self.btn_remove)
        btn_row.addWidget(self.btn_clear)
        lay.addLayout(btn_row)

        self.file_list = QListWidget()
        self.file_list.setSelectionMode(QAbstractItemView.SingleSelection)
        self.file_list.setMinimumHeight(110)
        self.file_list.currentRowChanged.connect(self._on_file_selected)
        lay.addWidget(self.file_list)

        return grp

    # ------------------------------------------------------------------
    def _build_tabs(self) -> QTabWidget:
        tabs = QTabWidget()

        # ── Tab 1: Complex Settings ─────────────────────────────────────
        self.settings_scroll = QScrollArea()
        self.settings_scroll.setWidgetResizable(True)
        self.settings_placeholder = QLabel(
            "Add files and select one to edit its settings."
        )
        self.settings_placeholder.setAlignment(Qt.AlignCenter)
        self.settings_scroll.setWidget(self.settings_placeholder)
        tabs.addTab(self.settings_scroll, "Complex Settings")

        # ── Tab 2: Layout & Export ──────────────────────────────────────
        layout_widget = QWidget()
        lf = QFormLayout(layout_widget)
        lf.setLabelAlignment(Qt.AlignRight)

        # subplot grid
        self.w_rows = QSpinBox(); self.w_rows.setRange(1, 10); self.w_rows.setValue(1)
        self.w_cols = QSpinBox(); self.w_cols.setRange(1, 10); self.w_cols.setValue(2)
        grid_row = QHBoxLayout()
        grid_row.addWidget(QLabel("Rows:")); grid_row.addWidget(self.w_rows)
        grid_row.addWidget(QLabel("Cols:")); grid_row.addWidget(self.w_cols)
        lf.addRow("Grid:", grid_row)

        # auto-layout button
        self.btn_auto_layout = QPushButton("Auto-fit grid")
        self.btn_auto_layout.clicked.connect(self._auto_layout)
        lf.addRow("", self.btn_auto_layout)

        # figure title
        self.w_title = QLineEdit()
        self.w_title.setPlaceholderText("Optional figure title …")
        lf.addRow("Figure title:", self.w_title)

        # panel letters
        self.w_panel_letters = QCheckBox("Show panel letters (A, B, C …)")
        self.w_panel_letters.setChecked(True)
        lf.addRow("", self.w_panel_letters)

        # subplot background
        self.w_subplot_bg = QComboBox()
        self.w_subplot_bg.addItems(["white", "black", "#eeeeee"])
        lf.addRow("Canvas background:", self.w_subplot_bg)

        # border
        self.w_border = QSpinBox(); self.w_border.setRange(0, 200); self.w_border.setValue(30)
        lf.addRow("Border (px):", self.w_border)

        # image size per panel
        self.w_img_w = QSpinBox(); self.w_img_w.setRange(200, 4000); self.w_img_w.setValue(800)
        self.w_img_h = QSpinBox(); self.w_img_h.setRange(200, 4000); self.w_img_h.setValue(600)
        sz_row = QHBoxLayout()
        sz_row.addWidget(QLabel("W:")); sz_row.addWidget(self.w_img_w)
        sz_row.addWidget(QLabel("H:")); sz_row.addWidget(self.w_img_h)
        lf.addRow("Panel size (px):", sz_row)

        # DPI
        self.w_dpi = QSpinBox(); self.w_dpi.setRange(72, 1200); self.w_dpi.setValue(300)
        lf.addRow("DPI:", self.w_dpi)

        # Ray-trace mode
        self.w_ray_mode = QComboBox()
        self.w_ray_mode.addItems(list(RAY_MODES.keys()))
        self.w_ray_mode.setCurrentIndex(0)
        lf.addRow("Ray-trace mode:", self.w_ray_mode)

        # label sizes
        self.w_label_size  = QSpinBox(); self.w_label_size.setRange(8, 72); self.w_label_size.setValue(28)
        self.w_caption_size = QSpinBox(); self.w_caption_size.setRange(8, 72); self.w_caption_size.setValue(22)
        lf.addRow("Panel letter size:", self.w_label_size)
        lf.addRow("Caption font size:", self.w_caption_size)

        # output path
        out_row = QHBoxLayout()
        self.w_out_path = QLineEdit()
        self.w_out_path.setPlaceholderText("Output file path …")
        btn_browse = QPushButton("Browse …")
        btn_browse.clicked.connect(self._browse_output)
        out_row.addWidget(self.w_out_path); out_row.addWidget(btn_browse)
        lf.addRow("Output file:", out_row)

        # format
        self.w_fmt = QComboBox()
        self.w_fmt.addItems(EXPORT_FORMATS)
        self.w_fmt.currentTextChanged.connect(self._sync_output_extension)
        lf.addRow("Format:", self.w_fmt)

        tabs.addTab(layout_widget, "Layout & Export")

        return tabs

    # ------------------------------------------------------------------
    def _build_preview_area(self) -> QGroupBox:
        grp = QGroupBox("Preview")
        lay = QVBoxLayout(grp)

        self.preview_label = QLabel("No preview yet.\nClick 'Quick Preview' to generate a low-resolution preview.")
        self.preview_label.setAlignment(Qt.AlignCenter)
        self.preview_label.setStyleSheet("background: #1a1a2e; color: #aaaaaa;")
        self.preview_label.setMinimumSize(500, 350)
        self.preview_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.preview_label)
        lay.addWidget(scroll, stretch=1)

        return grp

    # ------------------------------------------------------------------
    def _build_action_bar(self) -> QWidget:
        bar = QWidget()
        lay = QHBoxLayout(bar)
        lay.setContentsMargins(0, 0, 0, 0)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(False)
        lay.addWidget(self.progress_bar, stretch=1)

        self.btn_open_pymol = QPushButton("🔬 Open in PyMOL")
        self.btn_open_pymol.setToolTip("Launch selected file(s) in PyMOL interactive GUI")
        self.btn_open_pymol.clicked.connect(self._open_in_pymol_interactive)

        self.btn_open_vmd = QPushButton("🧬 Open in VMD")
        self.btn_open_vmd.setToolTip("Launch selected file(s) in VMD")
        self.btn_open_vmd.clicked.connect(self._open_in_vmd)

        self.btn_preview = QPushButton("🔍 Quick Preview")
        self.btn_preview.setToolTip("Low-resolution preview (fast, no ray-tracing)")
        self.btn_preview.clicked.connect(self._run_preview)

        self.btn_render = QPushButton("🖼 Render & Save")
        self.btn_render.setToolTip("Full publication-quality render")
        self.btn_render.setStyleSheet(
            "QPushButton { background: #2a6abb; color: white; font-weight: bold; padding: 6px 16px; }"
            "QPushButton:disabled { background: #555; }"
        )
        self.btn_render.clicked.connect(self._run_render)

        self.btn_cancel = QPushButton("⛔ Cancel")
        self.btn_cancel.setVisible(False)
        self.btn_cancel.clicked.connect(self._cancel_render)

        lay.addWidget(self.btn_open_pymol)
        lay.addWidget(self.btn_open_vmd)
        lay.addWidget(self.btn_preview)
        lay.addWidget(self.btn_render)
        lay.addWidget(self.btn_cancel)

        return bar

    # -----------------------------------------------------------------------
    # File management
    # -----------------------------------------------------------------------

    def _add_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Open Structure Files", "",
            "Structure Files (*.pdb *.mae *.maegz);;PDB (*.pdb);;MAE (*.mae);;MAEGZ (*.maegz);;All (*)"
        )
        for p in paths:
            if any(cx.filepath == p for cx in self.complexes):
                continue
            cx = ComplexEntry(p)
            self.complexes.append(cx)
            item = QListWidgetItem(cx.filename)
            item.setToolTip(p)
            self.file_list.addItem(item)
        self._update_button_states()
        self._auto_layout()

    def _remove_selected(self):
        row = self.file_list.currentRow()
        if row < 0:
            return
        self.file_list.takeItem(row)
        del self.complexes[row]
        self._update_button_states()
        self._auto_layout()

    def _clear_all(self):
        self.complexes.clear()
        self.file_list.clear()
        self._update_button_states()

    def _on_file_selected(self, row: int):
        if row < 0 or row >= len(self.complexes):
            self.settings_scroll.setWidget(self.settings_placeholder)
            return
        panel = ComplexSettingsPanel(self.complexes[row])
        panel.changed.connect(lambda: self._update_list_label(row))
        self.settings_scroll.setWidget(panel)

    def _update_list_label(self, row: int):
        if 0 <= row < self.file_list.count():
            self.file_list.item(row).setText(
                f"{self.complexes[row].filename}  [{self.complexes[row].label}]"
            )

    def _update_button_states(self):
        has = len(self.complexes) > 0
        self.btn_preview.setEnabled(has)
        self.btn_render.setEnabled(has)
        self.btn_remove.setEnabled(has)
        self.btn_clear.setEnabled(has)
        self.btn_open_pymol.setEnabled(has)
        self.btn_open_vmd.setEnabled(has)

    # -----------------------------------------------------------------------
    # Layout helpers
    # -----------------------------------------------------------------------

    def _auto_layout(self):
        """Square-ish grid that fits all loaded complexes."""
        n = len(self.complexes)
        if n == 0:
            return
        import math
        cols = math.ceil(math.sqrt(n))
        rows = math.ceil(n / cols)
        self.w_cols.setValue(cols)
        self.w_rows.setValue(rows)

    def _browse_output(self):
        fmt  = self.w_fmt.currentText()
        ext  = fmt.lower()
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Output Image", f"subplot.{ext}",
            f"{fmt} Image (*.{ext});;All (*)"
        )
        if path:
            self.w_out_path.setText(path)

    def _sync_output_extension(self, fmt: str):
        path = self.w_out_path.text()
        if path:
            p = Path(path)
            self.w_out_path.setText(str(p.with_suffix("." + fmt.lower())))

    # -----------------------------------------------------------------------
    # Render helpers
    # -----------------------------------------------------------------------

    def _collect_settings(self, fast: bool) -> dict:
        return {
            "dpi":               self.w_dpi.value(),
            "img_width":         self.w_img_w.value(),
            "img_height":        self.w_img_h.value(),
            "ray_mode":          RAY_MODES[self.w_ray_mode.currentText()],
            "title":             self.w_title.text().strip(),
            "label_size":        self.w_label_size.value(),
            "caption_size":      self.w_caption_size.value(),
            "subplot_bg":        self.w_subplot_bg.currentText(),
            "border_px":         self.w_border.value(),
            "show_panel_letters": self.w_panel_letters.isChecked(),
            "fast_preview":      fast,
        }

    def _run_preview(self):
        if not self.complexes:
            return
        out = os.path.join(tempfile.gettempdir(), "_pymol_preview.png")
        # Use lower resolution for speed
        settings = self._collect_settings(fast=True)
        settings["dpi"]       = 72
        settings["img_width"] = 400
        settings["img_height"] = 300
        self._start_worker(out, settings)

    def _run_render(self):
        if not self.complexes:
            return
        out = self.w_out_path.text().strip()
        if not out:
            QMessageBox.warning(self, "No output path",
                                "Please specify an output file path in the Layout & Export tab.")
            return
        self._start_worker(out, self._collect_settings(fast=False))

    def _start_worker(self, out_path: str, settings: dict):
        if self._worker and self._worker.isRunning():
            return

        layout = (self.w_rows.value(), self.w_cols.value())
        self._worker = RenderWorker(
            self.complexes, layout, out_path, settings, parent=self
        )
        self._worker.progress.connect(self._on_progress)
        self._worker.status.connect(self.statusBar().showMessage)
        self._worker.finished.connect(self._on_finished)
        self._worker.error.connect(self._on_error)

        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(True)
        self.btn_cancel.setVisible(True)
        self.btn_preview.setEnabled(False)
        self.btn_render.setEnabled(False)

        self._worker.start()

    def _cancel_render(self):
        if self._worker and self._worker.isRunning():
            self._worker.terminate()
            self._worker.wait()
        self._reset_ui()
        self.statusBar().showMessage("Cancelled.")

    # -----------------------------------------------------------------------
    # External viewer launchers
    # -----------------------------------------------------------------------

    def _files_to_open(self) -> List[str]:
        """Return paths of all loaded complexes (all files when none selected)."""
        row = self.file_list.currentRow()
        if row >= 0 and row < len(self.complexes):
            return [self.complexes[row].filepath]
        return [cx.filepath for cx in self.complexes]

    def _open_in_pymol_interactive(self):
        """Launch PyMOL GUI with the selected (or all) structure files."""
        files = self._files_to_open()
        if not files:
            return

        # Try to find the pymol executable next to the configured scripts path
        s = _load_app_settings()
        scripts_path = s.value(KEY_PYMOL_PATH, _DEFAULT_PYMOL_PATH)
        pymol_exe = None

        candidates: List[str] = []
        if scripts_path:
            # Schrödinger layout: Scripts/pymol.exe  or  Scripts/pymol
            candidates += [
                os.path.join(scripts_path, "pymol.exe"),
                os.path.join(scripts_path, "pymol"),
                os.path.join(str(Path(scripts_path).parent), "pymol.exe"),
                os.path.join(str(Path(scripts_path).parent), "pymol"),
            ]
        candidates += ["pymol"]  # system PATH fallback

        for c in candidates:
            if os.path.isfile(c) or shutil.which(c):
                pymol_exe = c
                break

        if not pymol_exe:
            QMessageBox.warning(
                self, "PyMOL not found",
                "Could not locate a PyMOL executable.\n"
                "Please set the correct path in Tools → Settings."
            )
            return

        try:
            self._launch_external([pymol_exe] + files)
            self.statusBar().showMessage(f"Launched PyMOL with {len(files)} file(s).")
        except Exception as exc:
            QMessageBox.critical(self, "Launch Error", str(exc))

    def _open_in_vmd(self):
        """Launch VMD with the selected (or all) structure files."""
        files = self._files_to_open()
        if not files:
            return

        vmd_path = _resolve_vmd_exe()

        # Resolve .lnk shortcuts on Windows via the shell
        if vmd_path.lower().endswith(".lnk"):
            self._launch_vmd_lnk(vmd_path, files)
            return

        if not (os.path.isfile(vmd_path) or shutil.which(vmd_path)):
            QMessageBox.warning(
                self, "VMD not found",
                f"VMD executable not found at:\n{vmd_path}\n\n"
                "Please set the correct path in Tools → Settings."
            )
            return

        # Decompress any .maegz files to temporary .mae files for VMD
        vmd_files, tmp_paths = self._prepare_vmd_files(files)
        try:
            self._launch_external([vmd_path] + vmd_files)
            self.statusBar().showMessage(f"Launched VMD with {len(vmd_files)} file(s).")
        except Exception as exc:
            QMessageBox.critical(self, "Launch Error", str(exc))
        finally:
            # Schedule cleanup of temp files after a short delay
            self._schedule_tmp_cleanup(tmp_paths)

    @staticmethod
    def _prepare_vmd_files(files: List[str]) -> Tuple[List[str], List[str]]:
        """Decompress .maegz files to temporary .mae files for VMD.

        Returns a tuple of (file_paths_for_vmd, temp_paths_to_cleanup).
        """
        out: List[str] = []
        tmp_paths: List[str] = []
        for f in files:
            if f.lower().endswith(".maegz"):
                tmp = tempfile.NamedTemporaryFile(suffix=".mae", delete=False)
                with gzip.open(f, "rb") as gz:
                    shutil.copyfileobj(gz, tmp)
                tmp.close()
                out.append(tmp.name)
                tmp_paths.append(tmp.name)
            else:
                out.append(f)
        return out, tmp_paths

    @staticmethod
    def _schedule_tmp_cleanup(paths: List[str]):
        """Remove temporary files after a short delay (30 s) to allow VMD to load them."""
        if not paths:
            return
        from PyQt5.QtCore import QTimer

        def _cleanup():
            for p in paths:
                try:
                    os.unlink(p)
                except OSError:
                    pass

        QTimer.singleShot(30_000, _cleanup)

    @staticmethod
    def _launch_external(cmd: List[str]):
        """Start an external process without blocking the GUI."""
        subprocess.Popen(cmd, close_fds=True)

    @staticmethod
    def _launch_vmd_lnk(lnk_path: str, files: List[str]):
        """Open a Windows .lnk shortcut via the shell with file arguments."""
        # 'start' can resolve .lnk but doesn't pass extra args reliably;
        # open via explorer with the shortcut then inform the user.
        try:
            subprocess.Popen(["cmd", "/c", "start", "", lnk_path], shell=False)
        except Exception:
            subprocess.Popen(["explorer", lnk_path])

    # -----------------------------------------------------------------------
    # Settings
    # -----------------------------------------------------------------------

    def _open_settings(self):
        dlg = SettingsDialog(self)
        if dlg.exec_() == QDialog.Accepted:
            _setup_pymol_path()
            self.statusBar().showMessage("Settings saved.")

    def _on_progress(self, val: int):
        self.progress_bar.setValue(val)

    def _on_finished(self, path: str):
        self._reset_ui()
        self.statusBar().showMessage(f"Done → {path}")

        # Show preview in the right panel
        if path.lower().endswith(".png"):
            pix = QPixmap(path)
        else:
            # For TIFF / PDF convert first panel to QPixmap via PIL
            try:
                img = Image.open(path)
                img.thumbnail((1200, 900))
                img = img.convert("RGB")
                data = img.tobytes("raw", "RGB")
                qimg = QImage(data, img.width, img.height, QImage.Format_RGB888)
                pix = QPixmap.fromImage(qimg)
            except Exception:
                pix = None

        if pix:
            scaled = pix.scaled(
                self.preview_label.size(),
                Qt.KeepAspectRatio,
                Qt.SmoothTransformation,
            )
            self.preview_label.setPixmap(scaled)

        msg = QMessageBox(self)
        msg.setWindowTitle("Render complete")
        msg.setText(f"Image saved to:\n{path}")
        msg.setIcon(QMessageBox.Information)
        msg.exec_()

    def _on_error(self, text: str):
        self._reset_ui()
        self.statusBar().showMessage("Error during rendering.")
        QMessageBox.critical(self, "Render Error", text)

    def _reset_ui(self):
        self.progress_bar.setVisible(False)
        self.btn_cancel.setVisible(False)
        self._update_button_states()


# ===========================================================================
# Entry point
# ===========================================================================

def main():
    # Inject PyMOL path early so imports work before the window is shown
    _setup_pymol_path()

    app = QApplication(sys.argv)
    app.setApplicationName("PyMOL Viewer")
    app.setStyle("Fusion")

    # Dark-ish palette
    palette = app.palette()
    palette.setColor(palette.Window,          QColor(45, 45, 48))
    palette.setColor(palette.WindowText,      QColor(220, 220, 220))
    palette.setColor(palette.Base,            QColor(30, 30, 30))
    palette.setColor(palette.AlternateBase,   QColor(53, 53, 53))
    palette.setColor(palette.ToolTipBase,     QColor(255, 255, 220))
    palette.setColor(palette.ToolTipText,     QColor(0, 0, 0))
    palette.setColor(palette.Text,            QColor(220, 220, 220))
    palette.setColor(palette.Button,          QColor(60, 60, 60))
    palette.setColor(palette.ButtonText,      QColor(220, 220, 220))
    palette.setColor(palette.Highlight,       QColor(42, 106, 187))
    palette.setColor(palette.HighlightedText, QColor(255, 255, 255))
    app.setPalette(palette)

    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
