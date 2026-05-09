# CBG Modelling

Computational modelling for the **Chemical Brain Gate (CBG)** project.
CBG is a system where neurons expressing a BBB-opening protein (under an activity-dependent promoter, e.g. c-Fos) locally disrupt tight junctions, allowing a systemically administered drug to enter the brain regionally.
The core value proposition: sub-therapeutic systemic drug doses achieve therapeutic local concentrations only at sites of pathological hyperactivity.

All models are parameterised against published literature values. No experimental work is conducted here.

---

## Quick start

```bash
# 1. Fetch source data (real two-photon vasculature + fallbacks)
python fetch_data.py

# 2. (First time only) Build Panel A data arrays
#    Open and run all cells of CBG_panel_A_dev.ipynb
#    → writes panel_a_data.npz and panel_a_provenance.json

# 3. Generate the combined figure
#    Open CBG_combined_figure.ipynb and run all cells (Cells 1–7)
#    → writes CBG_outputs/CBG_figure.pdf + .png

# 4. Verify before committing
python verify_figure.py   # must exit 0
```

---

## Repository structure

```
.
├── CBG_combined_figure.ipynb   # Main notebook — generates the publication figure
├── CBG_panel_A_dev.ipynb       # Panel A pipeline (vessel segmentation, PDE)
├── fetch_data.py               # Fetches real data with multi-tier fallbacks
├── verify_figure.py            # Automated pixel- and value-level verification
├── panel_a_data.npz            # Intermediate data produced by Panel A notebook
├── panel_a_provenance.json     # Source record for Panel A data
├── AGENTS.md                   # Agent instructions and model parameters
├── data/
│   ├── manifest.json           # Record of fetched vs synthetic sources
│   ├── atlas/                  # Allen CCF v3 annotation data
│   ├── enhancers/              # Allen ISH expression images
│   ├── focal/                  # Focal epilepsy model metadata
│   ├── raw/                    # MiniVess two-photon NIfTI volumes
│   └── vasculature/            # Processed vascular data
└── CBG_outputs/
    ├── CBG_figure.pdf          # Publication figure (PDF)
    ├── CBG_figure.png          # Publication figure (PNG, 600 dpi)
    └── CBG_figure_provenance.json
```

---

## Figure layout

| Region | Content |
|---|---|
| Left column, top | Panel A — IV infusion only: drug confined to vessels |
| Left column, bottom | Panel B — IV + CBG opening: drug leaks into active-region tissue |
| Right top | 4 × 5 small multiples: untreated cFos / untreated drug / treated cFos / treated drug across 5 time windows (W1–W5, 0–60 min) |
| Right bottom | Expression-weighted mean cFos and mean interstitial drug time courses |

---

## Model parameters

All parameters fixed across all tasks unless otherwise stated.

| Parameter | Symbol | Value | Unit | Source |
|---|---|---|---|---|
| ECF diffusion coefficient | D_tissue | 3×10⁻⁴ | cm²/min | Vendel et al. 2019, Fluids Barriers CNS |
| Intact BBB permeability | P_intact | 1×10⁻⁴ | cm/min | Bickel 2022, Pharmaceutics |
| Opened BBB permeability | P_open | 5×10⁻³ | cm/min | ~50× intact; FUS-BBBO literature |
| Brain elimination rate | k_e | 0.05 | /min | Vendel et al. 2019 |
| ECF volume fraction | φ | 0.20 | — | Nicholson & Hrabetova 2017 |

**Assumed (free) parameters** (labelled in provenance JSON):
`EC50`, `opening_gain`, `source_gain`, `interstitial_gain`, `open_rise`, `open_decay`, `supp_amp_gain`

---

## Data sources

| Data | Source | Licence |
|---|---|---|
| Two-photon vasculature | MiniVess, Poon et al. 2023 *Scientific Data* (EBRAINS) | CC BY-NC-SA |
| Allen CCF v3 annotation | Allen Institute for Brain Science | CC BY 4.0 |
| ISH expression images | Allen Brain Atlas API (Prox1, Wfs1) | Allen Institute terms |

`fetch_data.py` tries each source in order and falls back to synthetic data automatically.
The provenance JSON records which tier was actually used.

---

## Verification

`verify_figure.py` runs 20 checks (pixel-level and value-level) and exits 0 only if all pass.
**Do not commit until `verify_figure.py` exits 0.**

Commit messages must include the verification summary, e.g.:
```
CBG figure v2: 20/20 checks pass
```
