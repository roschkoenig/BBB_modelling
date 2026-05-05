# Task: Generate CBG Pilot Figure + Explainer Notebook

Generate two outputs:
1. `CBG_pilot_figure.pdf` + `CBG_pilot_figure.png`
2. `CBG_model_explainer.ipynb`

Use parameters and conventions from AGENTS.md.
Iterate until all quality checks in AGENTS.md and the task-specific checks
below all pass. Print PASS/FAIL for each check before finishing.

---

## Outputs

### Figure: 3 panels side by side, 18 cm wide x 7 cm tall

---

#### Panel A — Dose Amplification

**What it shows:** CBG amplifies local drug concentration at the seizure focus
without increasing systemic dose.

**Geometry:** 1D brain cross-section, x = 0–8 mm.
BBB open zone: gaussian profile centred at x = 4 mm, half-width = sigma.

**PDE:**
  dC/dt = D_tissue * d²C/dx² + P(x) * (C_blood - C_tissue) - k_e * C

  where P(x) = P_open inside opening zone, P_intact elsewhere.
  Use parabolic PDE framework (Charemis et al. 2026).
  Run to steady state (t = 60 min). Use C_blood = 1.0 (normalised).

**Plot two curves:**
- Blue: systemic only (P_intact everywhere)
- Red: CBG + systemic (P_open in focal zone)

**Add horizontal dashed lines:**
- Therapeutic threshold: 3x the systemic-only peak
- Side-effect threshold: 8x the systemic-only peak

---

#### Panel B — Safe Operating Window

**What it shows:** There is a designable threshold separating seizure-level
from physiological c-Fos activation, determined by enhancer expression density.

**2D heatmap. Axes (both log scale):**
- X: enhancer expression density — 10 to 10,000 neurons/mm³
- Y: c-Fos fold-change above baseline — 1x to 100x

**Compute for each point:**
  C_opener = density * fold_change * alpha
  Gate opens if C_opener > threshold_open

Set alpha and threshold_open so that:
- (1000/mm³, 30x) → RED (opens)
- (1000/mm³, 3x)  → BLUE (closed)
- (500/mm³, 2x)   → BLUE (closed)

**Colour:** RED = opens, BLUE = closed.
White contour at threshold boundary.

**Mark with symbols:**
- ★  Seizure: (1000/mm³, 30x) — must be RED
- ●  Learning: (1000/mm³, 3x) — must be BLUE
- ◆  Normal waking: (500/mm³, 2x) — must be BLUE

---

#### Panel C — Multi-Focal Network Targeting

**What it shows:** Activity-dependent BBB opening scales with seizure severity
at each focus; sub-therapeutic systemic dose reaches therapeutic levels only
at active foci.

**Geometry:** 2D cortical slice, 10 mm x 10 mm.

**Two foci:**
- Focus 1: position (3, 5) mm, c-Fos = 40x → opening radius r1
- Focus 2: position (7, 5) mm, c-Fos = 15x → opening radius r2 < r1
- Derive radii from Panel B threshold model:
  r = radius at which C_opener falls to threshold_open

**Solve 2D diffusion PDE for drug concentration across the slice.**

**Show as filled contour heatmap** with:
- Dashed white circles: BBB opening zone boundaries
- Solid contour line: therapeutic threshold
- Second contour (if reached): side-effect threshold

**Add inset:** 1D cross-section along y = 5 mm with both thresholds marked.

Set C_blood so that systemic-only gives < 10% of therapeutic threshold,
confirming that drug reaches therapeutic levels only within opening zones.

---

### Notebook: CBG_model_explainer.ipynb

Structure with these cells in order:

**Cell 1 — Markdown:** Title + 3–4 sentence plain-language description of CBG
and what this notebook models. Audience: non-specialist grant reviewer.

**Cell 2 — Markdown:** Physical model overview. Explain the PDE in plain
language. Include equation in LaTeX. One paragraph per panel. Include a simple
matplotlib schematic of the 1D geometry (Panel A) and 2D grid (Panel C).

**Cell 3 — Code:** Imports. Pin package versions in comments.
Include: numpy, scipy, matplotlib, ipywidgets, IPython.display.

**Cell 4 — Markdown + Code:** TUNABLE PARAMETERS.
Header: `## Tunable Parameters — change these to explore the model`

Define all parameters as named constants with inline source citations
(copy values from AGENTS.md parameter table).

Add ipywidgets sliders for these five parameters:
- P_open:          range 1e-4 to 1e-2,  step 1e-4
- sigma_mm:        range 0.1  to 2.0,   step 0.1
- C_blood:         range 0.1  to 2.0,   step 0.1
- alpha:           range 0.01 to 1.0,   step 0.01
- c_fos_seizure:   range 10   to 100,   step 5

Wire sliders to regenerate Panel A and Panel B live on change.
Add a plain-English sentence below each slider explaining what it controls.

**Cell 5 — Markdown + Code:** Panel A in isolation.
Below plot, print:
- Peak amplification ratio: Xx (benchmark: 10–100x)
- Therapeutic threshold reached within opening zone: YES/NO
- Side-effect threshold exceeded anywhere: YES/NO

**Cell 6 — Markdown + Code:** Panel B in isolation.
Below plot, print:
- Seizure (30x, 1000/mm³) in RED: YES/NO
- Learning (3x, 1000/mm³) in BLUE: YES/NO
- Normal (2x, 500/mm³) in BLUE: YES/NO

**Cell 7 — Markdown + Code:** Panel C in isolation.
Below plot, print:
- Focus 1 centre concentration: X (above therapeutic: YES/NO)
- Focus 2 centre concentration: X (above therapeutic: YES/NO)
- Off-target max concentration: X (below therapeutic: YES/NO)

**Cell 8 — Code:** Assemble and save the final figure (all 3 panels).
Print: `Figure saved: CBG_pilot_figure.pdf and CBG_pilot_figure.png`

**Cell 9 — Markdown:** Parameter sources table.
Columns: Parameter | Value used | Source | Assumption/caveat
One row per parameter. Flag any estimated values.

**Cell 10 — Markdown:** Limitations and next steps.
5–7 bullet points covering: what the model does not capture (tortuosity,
perivascular flow, protein binding, inflammation-driven permeability changes),
what experimental data would constrain each free parameter, and what the
Phase 2 model extension would include (multi-site, pharmacodynamic coupling,
patient-specific geometry from MRI).

---

## Task-specific quality checks

Run after generation. Print PASS/FAIL. Revise and regenerate if any FAIL.
Also run all universal checks from AGENTS.md.

- [ ] CHECK A1: Peak ratio (CBG / systemic) is 10–100x → else adjust P_open
- [ ] CHECK A2: CBG curve exceeds therapeutic threshold inside zone,
               falls below it outside → else adjust sigma or P_open
- [ ] CHECK A3: CBG curve stays below side-effect threshold everywhere
               → else lower C_blood
- [ ] CHECK B1: (1000/mm³, 30x) → RED; (1000/mm³, 3x) → BLUE
               → else adjust alpha
- [ ] CHECK B2: Boundary cleanly separates seizure from physiological regimes
               at some density value → else revise threshold model
- [ ] CHECK C1: Focus 1 centre exceeds therapeutic threshold
               → else adjust opening radius or C_blood
- [ ] CHECK C2: Off-target regions remain below therapeutic threshold
               → else lower C_blood
- [ ] CHECK N1: Notebook has exactly 10 cells; sliders functional;
               parameter table complete
- [ ] CHECK N2: Notebook Cell 4 parameter values exactly match values
               used to generate the saved figure

---

## Caption

Generate an 80–120 word figure caption alongside the figure.
State what each panel shows, name the key result in each, cite parameter sources.
Do not use the word "novel."
Write this into Cell 8 as a markdown sub-cell and print it as standalone output.
