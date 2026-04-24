# PYMOL_VIEWER

A PyQt5 GUI for generating **publication-quality images** of multiple
protein–ligand docking complexes arranged in a **subplot grid**, using PyMOL
as the rendering back-end, with optional interactive launch via PyMOL or VMD.

---

## Features

| Feature | Detail |
|---|---|
| **File formats** | `.pdb`, `.mae`, `.maegz` |
| **Multiple complexes** | Load any number; displayed in an N × M grid |
| **Per-complex settings** | Protein/ligand representation, colour scheme, surface, H-bonds, background |
| **Docking analysis** | Binding-site residue labels, contact distance measurements, H-bond display |
| **Subplot compositing** | Automatic panel letters (A, B, C …), optional captions and figure title |
| **Export formats** | PNG, TIFF, PDF |
| **DPI control** | 72 – 1200 DPI (300 DPI default for publications) |
| **Ray-trace modes** | Default colour + outline / Cartoon / Black outline / Quantised |
| **Quick Preview** | Low-resolution, no ray-tracing – instant feedback |
| **Open in PyMOL** | Launch files in PyMOL interactive GUI (Schrödinger or open-source) |
| **Open in VMD** | Launch files in VMD for further analysis |

---

## Requirements

```
Python ≥ 3.8
PyMOL (open-source or Schrödinger bundle)
PyQt5
Pillow
matplotlib
```

### Install dependencies

```bash
pip install PyQt5 Pillow matplotlib

# PyMOL (open-source edition):
pip install pymol-open-source
```

> **MAE / MAEGZ support**: the open-source PyMOL build reads `.mae` files
> natively.  `.maegz` files are automatically decompressed before loading.
> The Schrödinger commercial PyMOL bundle provides the fullest MAE support.

---

## Configuring PyMOL and VMD paths

Open **Tools → Settings** (or press `Ctrl+,`) and set:

| Setting | Default / Example |
|---|---|
| **PyMOL Scripts path** | `C:\Users\<username>\AppData\Local\Schrodinger\PyMOL2\Scripts` |
| **VMD executable** | `C:\Users\<username>\Desktop\VMD 2.0.lnk` or the full path to `vmd.exe` |

These settings are saved between sessions.  The PyMOL scripts directory is
automatically added to `sys.path` so that `import pymol2` resolves to your
chosen build.

---

## Usage

```bash
python pymol_viewer.py
```

### Workflow

1. **Add Files** – click *➕ Add Files …* and select one or more `.pdb`,
   `.mae`, or `.maegz` files.
2. **Complex Settings** (left tab) – for each selected file, customise:
   - Ligand PyMOL selection (default `organic`; use `resn LIG` or chain
     notation for specific residues)
   - Protein / ligand representations and colour schemes
   - Optional surface with transparency
   - H-bond display
   - **Docking analysis**: residue labels (radius & colour), contact distance
     lines (cutoff), extra zoom-out
3. **Layout & Export** (right tab) – set the grid dimensions (or click
   *Auto-fit grid*), figure title, panel letters, DPI, output path and format.
4. **Quick Preview** – generates a fast, low-resolution composite to check
   the layout before the final render.
5. **Render & Save** – performs full ray-traced rendering and writes the
   output file.
6. **Open in PyMOL / VMD** – launch the loaded file(s) directly in the
   corresponding interactive viewer for further analysis.

---

## Docking analysis features

| Setting | Description |
|---|---|
| **Label binding-site residues** | Displays `ResnResi` labels on Cα atoms within the specified radius of the ligand |
| **Label radius (Å)** | Distance from the ligand used to select residues to label (default 5 Å) |
| **Label color** | Text colour of residue labels |
| **Show contact distances** | Draws dashed lines between ligand and protein atoms within the cutoff |
| **Distance cutoff (Å)** | Maximum distance for contact lines (default 4 Å) |
| **Extra zoom-out (Å)** | Increase the zoom buffer beyond the default 6 Å to show more context |

---

## Ligand selection syntax

PyMOL selection language is used for the *Ligand selection* field:

| Example | Selects |
|---|---|
| `organic` | All non-polymer organic molecules (default) |
| `resn LIG` | Residue named LIG |
| `resn ATP+ADP` | Residues ATP or ADP |
| `/complex//A/LIG` | Chain A, residue LIG in object *complex* |

---

## Output example

```
subplot.png  (300 DPI, 2 × 2 grid, white background)
├── A  1ABC – cartoon protein + sticks ligand + residue labels
├── B  2XYZ – surface + sticks ligand + H-bonds
├── C  3DEF – ribbon + ball-and-stick ligand + contact distances
└── D  4GHI – cartoon + sticks, H-bonds + residue labels
```

---

## File format notes

| Extension | Notes |
|---|---|
| `.pdb` | Standard Protein Data Bank format – universal support |
| `.mae` | Schrödinger Maestro format – loaded directly by PyMOL |
| `.maegz` | Gzip-compressed MAE – auto-decompressed before loading |
