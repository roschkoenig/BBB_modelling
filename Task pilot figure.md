# Task: Generate CBG Pilot Figure + Explainer Notebook

Generate two outputs:
1. `CBG_pilot_figure.pdf` + `CBG_pilot_figure.png`
2. `CBG_model_explainer.ipynb`

Use parameters and conventions from AGENTS.md.
Run ALL critic layer tests (Tests 1–5 in AGENTS.md) after generating.
State pass/fail explicitly. Revise and regenerate if any fail.

Target aesthetic: Nature Communications. Every visual element must be either
derived from the model with cited parameters, or explicitly labelled as assumed.

---

## Figure layout

Total: 170 mm wide × 150 mm tall.

**Top row** (55 mm tall): Panels A, B, C side by side.
  A: 35% width. B: 35%. C: 30%.

**Bottom row** (85 mm tall): Panel D — full width.
  Internally: thin time-course strip (15 mm) above a 4×5 spatial grid (65 mm).

Panel labels A B C D bold 8pt top-left of each panel.

---

## Panel A: Vascular-constrained drug delivery

**Claim:** "Drug entry without BBB opening is spatially uniform and vessel-proximal;
CBG creates discrete high-concentration hotspots at vessels within the
expressing neuronal population."

Three sub-panels stacked vertically (labels i ii iii). 2×2 mm tissue patch.
Shared viridis colourbar. 200 µm scale bar on sub-panel i.

**Vascular tree:** Fractal branching, Murray's law, target density 400–600
vessels/mm² (Tsai et al. 2009 PNAS), ±10° tortuosity per 50 µm segment.

**2D PDE (finite difference dx=10 µm, steady state):**
  dC/dt = D * ∇²C - k_e * C + P(x,y) * (C_blood - C) * delta_vessel(x,y)
  delta_vessel: Gaussian kernel sigma=10 µm on vessel walls.

Sub-panel i: P = P_intact everywhere. Steady-state heatmap.
Sub-panel ii: BBB-opener expression mask only (circular patch r=0.7 mm,
  Gaussian edge). Colour = expression density. Vessel tree as white lines.
Sub-panel iii: P = P_open where (mask > threshold) AND (within 30 µm of vessel).
  Hotspots must appear at vessels inside mask only.

---

## Panel B: Anatomical specificity — single time point reference

**Claim:** "CBG activation with hippocampus-specific enhancers produces bilateral
drug accumulation in hippocampal grey matter that respects white matter boundaries."

Coronal mouse brain at Bregma -2.0 mm. Anatomical polygons: cortex, CA1/CA3/DG
(bilateral), corpus callosum, thalamus, lateral ventricles, fimbria.

Two-compartment diffusion: grey matter D=3e-4, white matter D=5e-5 cm²/min
(Vorisek & Sykova 1997). BBB open in bilateral CA1 and DG only.
Steady-state drug concentration heatmap overlaid on outlines.
White matter must be visibly darker than hippocampal grey.
Add CA1, DG, CC labels (6pt italic). 1 mm scale bar.

---

## Panel C: Single-event time course — the feedback cycle

**Claim:** "A single seizure event triggers the full CBG negative feedback cycle:
pathological activity drives BBB opener accumulation, which gates drug entry,
which terminates the seizure — all within a ~60 min window."

This panel is the KEY to reading Panel D. It shows the variables that become
the spatial heatmaps in Panel D.
Time axis 0–90 min. Two stacked sub-panels sharing x-axis.

Top: c-Fos fold-change A(t). Red = seizure, blue dashed = physiological.
  Dashed horizontal at A_thresh ("BBB-opening threshold").
  Shaded region showing BBB opener P(t) accumulation window.

Bottom: Drug concentration C(t). Red = CBG active. Grey flat = no CBG.
  Dashed lines: therapeutic threshold (green, "assumed") and
  side-effect threshold (orange, "assumed").

**ODE system:**
  dA/dt = k_A * S(t) - lambda_A * A
  dP/dt = k_P * max(A - A_thresh, 0) - lambda_P * P
  dC/dt = k_C * max(P - P_thresh, 0) * (C_blood - C) - lambda_C * C
  dS/dt = -k_suppress * C * S    ← drug terminates seizure activity

Parameters: k_A=2.0, lambda_A=0.05, A_thresh=3.0 (Tullai 2007; Bhatt 2020),
k_P=0.1, lambda_P=0.015 (estimated from CCL2), k_C=0.5, lambda_C=0.05
(Vendel 2019), P_thresh=0.3 (LABEL AS ASSUMED), k_suppress=0.8 (LABEL AS ASSUMED).

Seizure: S₀=15 for t∈[10,40] min, subject to suppression. Linear axes only.

Mark the 5 time points used in Panel D as vertical dashed lines on Panel C
(e.g. t = 0, 15, 30, 45, 70 min) so the reader can map Panel D back to Panel C.

---

## Panel D: Spatial evolution — 4 rows × 5 time snapshots

**Claim:** "Without CBG, seizure-driven c-Fos expression spreads through the
hippocampus while drug concentration remains zero throughout. With CBG, drug
entry coincides spatially and temporally with peak c-Fos expression and
accelerates seizure resolution, visible as localised drug hotspots that
track and then outlast the activity."

**Layout:**
Full panel width. 4 rows × 5 columns of spatial heatmaps.
Each small panel: same coronal hippocampal geometry as Panel B.
Each small panel approx 28 mm wide × 14 mm tall.

**The 5 columns are 5 time snapshots** of the hippocampal section as the
seizure evolves. Use the same 5 time points marked in Panel C:
  t₁ = 0 min    (pre-seizure baseline)
  t₂ = 15 min   (early seizure, c-Fos rising)
  t₃ = 30 min   (seizure peak)
  t₄ = 45 min   (resolution phase, with vs without CBG diverges)
  t₅ = 70 min   (post-seizure)

Time point labels as column headers above the grid: "t=0", "t=15 min", etc.

**The 4 rows are:**
  Row 1: c-Fos spatial distribution — NO CBG condition
  Row 2: Drug concentration — NO CBG condition
  Row 3: c-Fos spatial distribution — CBG active condition
  Row 4: Drug concentration — CBG active condition

Row labels on left margin (7pt): "c-Fos (no CBG)", "Drug (no CBG)",
  "c-Fos (CBG)", "Drug (CBG)"

**How to generate the spatial heatmaps:**

The hippocampal c-Fos spatial distribution at each time point is derived from
the ODE solution. c-Fos expression A(t) from the ODE defines the AMPLITUDE of
the spatial c-Fos pattern. The spatial PATTERN is a 2D Gaussian seeded at the
seizure initiation site (left CA1) and spreading across the hippocampus with a
diffusion-like spread rate of 0.3 mm/min (consistent with hippocampal seizure
propagation speeds; Trevelyan et al. 2006, J Neurosci).

c-Fos spatial field at time t:
  F(x,y,t) = A(t) * G(x,y,t)
  where G is a 2D spreading Gaussian:
    G(x,y,t) = exp(-r²/ (2 * (sigma_0 + v_spread * t)²))
  sigma_0 = 0.3 mm (initial focus), v_spread = 0.05 mm/min

Drug concentration spatial field at time t:
  Use the same coronal geometry as Panel B (grey/white matter compartments).
  Drug source at each time point = k_C * max(P(t) - P_thresh, 0) * C_blood
  applied only within regions where F(x,y,t) > F_thresh (BBB opens where
  neurons are sufficiently active AND expressing the transgene).
  Solve 2D diffusion at each time point (use steady-state approximation within
  each snapshot — drug diffusion is fast relative to ODE timescale).
  NO CBG condition: drug source = 0 everywhere at all times.

**Colourscales:**
  Rows 1 & 3 (c-Fos): hot colourmap (black-red-yellow). Shared scale
    0 to max(A(t)) across all panels in rows 1 and 3.
  Rows 2 & 4 (drug): viridis. Shared scale 0 to max_C across all panels
    in rows 2 and 4 (CBG condition sets the max; no-CBG will be uniformly dark).
  Single colourbar for rows 1&3, placed right of row 1.
  Single colourbar for rows 2&4, placed right of row 2.

**Key visual results that must be visible:**
1. Row 1 (c-Fos, no CBG): spatial spread of activity from t₂ to t₃,
   then persistence through t₄ and t₅ (seizure runs full course uninterrupted).
2. Row 2 (drug, no CBG): uniformly dark (near-zero) at all time points.
   The visual message is silence — nothing is entering the brain.
3. Row 3 (c-Fos, CBG): same spatial spread at t₂ and t₃ as Row 1, but
   at t₄ the activity is visibly less intense — seizure is terminating early
   because drug has entered. By t₅ the field is darker than in Row 1.
4. Row 4 (drug, CBG): near-zero at t₁ and t₂ (drug not yet entered),
   bright hotspots appear at t₃ (drug entering at peak BBB-opener window),
   hotspots persist at t₄ then begin to decay at t₅.
   Crucially: drug hotspots are spatially co-localised with the c-Fos pattern
   in Row 3, not diffuse — they track the pathological activity.

**Formatting:**
- No tick labels on any small panel.
- No axes on small panels — clean image, just the heatmap and the anatomical outlines.
- Anatomical outlines as thin white lines (0.3pt) within each panel.
- Column separator: 1.5 mm gap. Row separator: 2 mm gap.
- Thicker separator (3 mm) between rows 2 and 3 (condition break).
- Small horizontal bar between rows 2 and 3 with labels "No CBG" and "CBG active".

---

## Critic layer (mandatory, run before saving)

```
CRITIC ASSESSMENT:
Test 1 (Single-claim):
  Panel A: [sentence] — PASS/FAIL
  Panel B: [sentence] — PASS/FAIL
  Panel C: [sentence] — PASS/FAIL
  Panel D: [sentence] — PASS/FAIL

Test 2 (Grounded parameters):
  [Each threshold/boundary → source or ASSUMED label]

Test 3 (Biological realism):
  [Vascular branching / seizure spatial spread / tissue compartments] — PASS/FAIL

Test 4 (Mock reviewer — two sentences):
  [Write them] — PASS/FAIL

Test 5 (Story coherence — A→B→C→D):
  [Progressive argument?] — PASS/FAIL

OVERALL: PASS / FAIL
```

Revise and regenerate any failing panels before saving.

---

## Notebook: CBG_model_explainer.ipynb

Cell 1: Title + 3-sentence description.
Cell 2: Physical model overview — LaTeX for A, B, C ODEs. Explain seizure
  spatial propagation model used in Panel D. Explain feedback term explicitly.
Cell 3: Imports + pinned versions.
Cell 4: TUNABLE PARAMETERS with ipywidgets sliders:
  P_open, A_thresh, P_thresh, k_suppress, v_spread (seizure propagation speed).
  Wire to regenerate Panel C time course and Panel D grid live.
Cell 5: Panel A code.
Cell 6: Panel B code.
Cell 7: Panel C code. Verify 5 time point markers match Panel D columns.
Cell 8: Panel D code. Verify spatial snapshots match Panel C time course.
Cell 9: Assemble and save. Print CRITIC ASSESSMENT block.
Cell 10: Parameter table (value | source | ASSUMED flag).
Cell 11: Limitations and next steps.

---

## Final checks

- [ ] A1: CBG hotspot 10–100x systemic-only
- [ ] A2: Off-target vessels at baseline
- [ ] B1: White matter darker than hippocampal grey
- [ ] C1: Seizure condition drug enters therapeutic range
- [ ] C2: Physiological condition stays below threshold
- [ ] C3: 5 time points marked as vertical lines in Panel C
- [ ] D1: Panel D column times match Panel C markers exactly
- [ ] D2: Row 2 (drug, no CBG) uniformly dark at all time points
- [ ] D3: Row 3 (c-Fos, CBG) shows visibly earlier termination vs Row 1 at t₄ and t₅
- [ ] D4: Row 4 (drug, CBG) hotspots are spatially co-localised with Row 3 c-Fos at t₃
- [ ] D5: Shared colourscale within rows 1&3; shared within rows 2&4
- [ ] N2: P_thresh and k_suppress labelled ASSUMED in figure and notebook
- [ ] CRITIC: All 5 tests pass
