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
| ECF volume fraction | phi | 0.20 | — | Nicholson & Hrabetova 2017 |
| Protein diffusion (large) | D_protein | 1e-5 | cm²/min | Nance et al. 2014 |
| c-Fos lag (physiol. activity) | fold_phys | 2–5 | x baseline | Bhatt et al. 2020, Nat Neurosci |
| c-Fos induction (seizure) | fold_sz | 20–50 | x baseline | Qiu et al. 2022, Science |
| BBB opener protein t½ | t_half_prot | 45 | min | estimated from CCL2 literature |

**Benchmarks:**
- Osmotic BBB amplification for 153 kDa antibody: ~9.5x (Chu et al. 2019)
- Small molecule amplification expected: 10–100x

---

## Code conventions

- Language: Python 3.10+
- Required packages: numpy, scipy, matplotlib, ipywidgets, jupyter
- PDE solver: scipy.integrate.solve_ivp (method='Radau') or finite difference
- **Figure style — Nature Communications standard:**
  - Font: Arial or Helvetica, 7pt axis labels, 8pt panel labels
  - Layout: single-column 85 mm wide or double-column 170 mm wide
  - Panel labels: bold uppercase A, B, C top-left, 8pt
  - Line weights: 0.75pt axes, 1.0pt data lines
  - Colormaps: viridis, inferno, RdBu_r (perceptually uniform only)
  - No legend box; no gridlines unless data requires
  - No chartjunk — remove all non-essential ink
  - Colour-blind safe (simulate deuteranopia before accepting)
  - Scale bars on spatial panels; no axis tick labels on heatmaps unless essential
- Output: PDF + PNG (600 dpi) for all figures

---

## Universal quality checks

- [ ] CBG/systemic peak ratio 10–100x for small molecules
- [ ] Off-target regions below therapeutic threshold
- [ ] Notebook parameters exactly match figure parameters
- [ ] All axes labelled with units; scale bars on spatial panels
- [ ] No labels truncated

---

## MANDATORY CRITIC LAYER — run before accepting any figure

This is not optional. Every figure must pass all four tests before it is saved.
If any test fails, state which test failed and why, revise, and rerun.

### Test 1: Single-claim test
For each panel, write one sentence:
"This panel shows that [X] because [Y], which supports [Z]."
The sentence must be falsifiable. If it contains "illustrates," "demonstrates,"
or "suggests" as weasel words with no quantitative backing, the panel fails.
Revise until the claim is specific and defensible.

### Test 2: Grounded parameter test
For every threshold, boundary, contour line, or colour scale break, state the
source. If the source is "chosen to make the figure look convincing," either:
(a) replace with a literature-derived value, OR
(b) label explicitly in the figure as "assumed" or "free parameter."
Unlabelled asserted thresholds are a Nature-tier rejection criterion.

### Test 3: Biological realism test
Spatial panels must reflect real biology, not mathematical convenience.
Specifically:
- Drug distributions must follow vessel geometry, not smooth Gaussians
- Anatomical boundaries (grey/white matter, laminar structure) must produce
  visible concentration differences
- Vascular patterns must show branching and tortuosity, not uniform grids
If any panel looks like a textbook diagram rather than data, it fails this test.

### Test 4: Mock reviewer test
Write two sentences as a hostile but fair Nature Communications reviewer.
Acceptable reviewer comment: "The model produces convincing results consistent with
the literature but the authors should clarify the basis for the P_open value."
Failing reviewer comment: "The parameters appear chosen to support the
conclusion; no independent validation is shown. The figures are not clear."
If your mock reviewer writes the failing version, revise before accepting.

### Test 5: Story coherence test
Remove all axis labels and captions mentally. Can a field-literate reader
still follow the argument from panel A to C? Each panel must add something
the previous panel cannot show. If panels A and B make the same point at
different scales, one of them is redundant.

---

## Output file naming

| Output type | Naming convention |
|---|---|
| Publication figure | CBG_<descriptor>_figure.pdf + .png |
| Jupyter notebook | CBG_<descriptor>_explainer.ipynb |
| Intermediate data | CBG_<descriptor>_data.npz |

---

## Consistency rule

Notebook and figure always generated from the same parameter set in the same run.
