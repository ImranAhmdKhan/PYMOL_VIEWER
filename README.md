# Molecular Viewer & Analytics

A pure **R / Shiny** dashboard for protein–ligand docking analysis and
publication-quality visualisation. No Python, no PyMOL, no external binaries
required — just R.

---

## Features

| Feature | Detail |
|---|---|
| **File formats** | `.pdb` |
| **Multiple complexes** | Upload any number of PDB files simultaneously |
| **H-bond analysis** | Donor/acceptor N & O atom proximity (configurable cutoff) |
| **Close-contact count** | All heavy-atom contacts within a configurable distance |
| **Binding-site residues** | Protein Cα atoms within a configurable radius of the ligand |
| **B-factor statistics** | Per-complex mean ± SD B-factor |
| **7 plot templates** | Bar chart, scatter, frequency chart, radar, heatmap, property space, B-factor |
| **Interactive plots** | Plotly-powered hover tooltips |
| **Export formats** | PNG (with DPI control), PDF, SVG |
| **CSV / RDS export** | Raw analysis data for further processing |
| **Structure detail tab** | Atom table, binding-site residue list, per-complex summary |

---

## Requirements

- **R ≥ 4.1**

All R packages are installed automatically when `app.R` starts, or you can
install them manually:

```r
Rscript install_packages.R
```

### R packages used

| Package | Purpose |
|---|---|
| `shiny` | Web application framework |
| `shinydashboard` | Dashboard layout |
| `shinyWidgets` | Enhanced UI widgets |
| `colourpicker` | Colour picker widget |
| `bio3d` | PDB file I/O, distance calculations |
| `ggplot2` | Publication-quality static plots |
| `ggrepel` | Non-overlapping text labels |
| `plotly` | Interactive plots |
| `DT` | Interactive data tables |
| `dplyr` / `tidyr` | Data wrangling |
| `scales` | Colour palettes |
| `fmsb` | Radar / spider charts |

---

## Usage

### Launch the app

```r
# From an R console:
shiny::runApp("app.R")

# Or from the terminal:
Rscript -e "shiny::runApp('app.R')"
```

The app opens in your default browser at `http://127.0.0.1:<port>`.

---

## Workflow

### 1 — Load Structures tab

1. Click **Browse…** and select one or more `.pdb` files.
2. Optionally enter comma-separated **labels** (used as complex names in plots).
3. Adjust analysis parameters:
   - **Contact cutoff (Å)** – maximum distance for a close contact (default 4 Å).
   - **H-bond cutoff (Å)** – maximum N/O distance for an H-bond (default 3.5 Å).
   - **Binding-site radius (Å)** – Cα atoms within this distance of any ligand atom are counted as binding-site residues (default 5 Å).
4. Click **🔬 Analyse Structures**.
5. The summary table updates with per-complex metrics.
6. Use **📥 Download CSV** or **💾 Save Analysis (.rds)** to export results.

### 2 — Analysis & Plots tab

1. Select a **Plot template** from the drop-down.
2. Click **▶ Render Plot** to generate the static ggplot.
3. The interactive Plotly version updates automatically.
4. Set **Export format**, **DPI**, and dimensions, then click **💾 Save Plot**.

### 3 — Structure Detail tab

Select a complex to view:
- A text summary (residue count, chain count, B-factor, etc.)
- The binding-site residue list
- A scrollable table of the first 500 atoms from the PDB file

### 4 — Settings tab

Adjust default cutoffs and plot aesthetics. View installed package versions.

---

## Plot Templates

| Template | Description |
|---|---|
| **H-Bond Bar Chart** | Bar chart of H-bond counts per complex |
| **H-Bonds vs. Contacts Scatter** | Scatter plot; bubble size = binding-site size |
| **Binding-Site Residue Frequency** | Horizontal bar chart of the most frequent binding-site residues |
| **Interaction Radar Chart** | Spider/radar chart of normalised interaction metrics |
| **Fingerprint Heatmap** | Presence/absence heatmap of binding-site residues across complexes |
| **Ligand Property Space** | 2-D plot: heavy atoms vs. H-bonds, coloured by contact count |
| **B-Factor Comparison** | Mean ± SD B-factor bar chart per complex |

---

## Exported CSV columns

| Column | Description |
|---|---|
| `name` | Complex label |
| `filename` | Source PDB file name |
| `hbond_count` | Number of H-bonds (N/O proximity) |
| `contact_count` | Number of close contacts within `contact_cutoff` |
| `binding_site_size` | Number of Cα atoms within `site_radius` of the ligand |
| `binding_site_residues` | Semicolon-separated residue labels |
| `ligand_atoms` | Ligand heavy-atom count |
| `contact_cutoff` | Distance cutoff used (Å) |
| `site_radius` | Binding-site radius used (Å) |
| `bfactor_mean` | Mean B-factor of protein ATOM records |
| `bfactor_sd` | Standard deviation of B-factor |
| `n_residues` | Number of Cα atoms (proxy for residue count) |
| `n_chains` | Number of unique chain IDs |

---

## Notes

- **Ligand detection**: atoms with `HETATM` record type that are *not* water (`HOH`) are treated as ligand atoms.
- **H-bond approximation**: donor/acceptor pairs are identified as N and O atoms (and F) within `hbond_cutoff` of each other. This is a geometric approximation; no hydrogen placement or angle criterion is applied.
- **Contact count**: counts pairs of ligand heavy atoms and protein Cα atoms within `contact_cutoff`.
