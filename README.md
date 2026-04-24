# PYMOL_VIEWER

A PyQt5 GUI for generating **publication-quality images** of multiple
protein–ligand complexes arranged in a **subplot grid**, using PyMOL as the
rendering back-end.

---

## Features

| Feature | Detail |
|---|---|
| **File formats** | `.pdb`, `.mae`, `.maegz` |
| **Multiple complexes** | Load any number; displayed in an N × M grid |
| **Per-complex settings** | Protein/ligand representation, colour scheme, surface, H-bonds, background |
| **Subplot compositing** | Automatic panel letters (A, B, C …), optional captions and figure title |
| **Export formats** | PNG, TIFF, PDF |
| **DPI control** | 72 – 1200 DPI (300 DPI default for publications) |
| **Ray-trace modes** | Default colour + outline / Cartoon / Black outline / Quantised |
| **Quick Preview** | Low-resolution, no ray-tracing – instant feedback |

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
3. **Layout & Export** (right tab) – set the grid dimensions (or click
   *Auto-fit grid*), figure title, panel letters, DPI, output path and format.
4. **Quick Preview** – generates a fast, low-resolution composite to check
   the layout before the final render.
5. **Render & Save** – performs full ray-traced rendering and writes the
   output file.

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
├── A  1ABC – cartoon protein + sticks ligand
├── B  2XYZ – surface + sticks ligand
├── C  3DEF – ribbon + ball-and-stick ligand
└── D  4GHI – cartoon + sticks, H-bonds shown
```

---

## File format notes

| Extension | Notes |
|---|---|
| `.pdb` | Standard Protein Data Bank format – universal support |
| `.mae` | Schrödinger Maestro format – loaded directly by PyMOL |
| `.maegz` | Gzip-compressed MAE – auto-decompressed before loading |
