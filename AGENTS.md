# CBG Modelling Agent — Persistent Context

## What this repository is

Computational modelling for the **Chemical Brain Gate (CBG)** project.
CBG is a system where neurons expressing a BBB-opening protein (under an
activity-dependent promoter, e.g. c-Fos) locally disrupt tight junctions,
allowing a systemically administered drug to enter the brain regionally.
The core value proposition: sub-therapeutic systemic drug doses achieve
therapeutic local concentrations only at sites of pathological hyperactivity.

This repo contains multi-compartment PDE-based models of BBB opening dynamics,
drug delivery, and machine-learning-based inference. No experimental work is
conducted here — all models are parameterised against published literature values.

---

## Model parameters (fixed across all tasks)

| Parameter | Symbol | Value | Unit | Source |
|---|---|---|---|---|
| ECF diffusion coefficient | D_tissue | 3e-4 | cm²/min | Vendel et al. 2019, Fluids Barriers CNS |
| Intact BBB permeability | P_intact | 1e-4 | cm/min | Bickel 2022, Pharmaceutics |
| Opened BBB permeability | P_open | 5e-3 | cm/min | ~50x intact; FUS-BBBO literature |
| Brain elimination rate | k_e | 0.05 | /min | Vendel et al. 2019 |
| ECF volume fraction | phi | 0.20 | — | Nicholson & Hrabětová 2017 |
| Protein diffusion (large) | D_protein | 1e-5 | cm²/min | Nance et al. 2014 |
| BBB opening half-width | sigma | 0.5 | mm | Estimated; FUS-BBBO analogy |

**Benchmark values to validate against:**
- Osmotic BBB opening amplification for 153 kDa antibody: ~9.5x
  (Chu et al. 2019, J Control Release — 33.9 vs 3.56 %ID/g)
- Small molecule amplification expected: 10–100x (higher than antibody)
- Physiological c-Fos induction (exploration, fear conditioning): 2–5x baseline
- Seizure-level c-Fos induction: 20–50x baseline (use 30x as conservative estimate)
- Physiological activity → c-Fos → BBB opener: must NOT reach opening threshold
- Seizure activity → c-Fos → BBB opener: must exceed opening threshold

---

## Code conventions

- Language: Python 3.10+
- Required packages: numpy, scipy, matplotlib, ipywidgets, jupyter
- PDE solver: scipy.integrate.solve_ivp (method='Radau') or finite difference
- Figure style:
  - Font: Arial or DejaVu Sans, minimum 7pt
  - Layout: Nature single-column standard (18 cm wide)
  - Palette: colourblind-safe (viridis or RdBu_r for heatmaps;
    blue/red for two-condition comparisons)
  - No gridlines on spatial panels
  - Bold panel labels (A, B, C) upper-left of each panel
  - Legends inside panels, no box
- Output files: save PDF + PNG (300 dpi) for all figures
- Notebooks: all parameter values defined in a single PARAMETERS cell;
  notebook parameter values must exactly match values used in saved figures

---

## Quality checks (apply to every figure-generating task)

These checks apply universally. Task files may add task-specific checks.

- [ ] Peak concentration ratio (CBG / systemic-only) is 10–100x for small molecules
- [ ] Off-target regions remain below therapeutic threshold
- [ ] Notebook parameter values match figure parameter values exactly
- [ ] All axis labels present with units
- [ ] No panel labels or tick labels truncated
- [ ] Parameter table in notebook cites a source for every value

---

## Output file naming

| Output type | Naming convention |
|---|---|
| Publication figure | `CBG_<descriptor>_figure.pdf` + `.png` |
| Jupyter notebook | `CBG_<descriptor>_explainer.ipynb` |
| Intermediate data | `CBG_<descriptor>_data.npz` |

---

## Consistency rule

**The notebook and the figure are always generated from the same parameter set
in the same script run.** Never update one without regenerating the other.
