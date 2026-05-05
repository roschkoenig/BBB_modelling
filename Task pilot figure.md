# Task: Generate CBG Pilot Figure + Explainer Notebook

## STEP 0 — Run fetch_data.py first

Before writing any figure code:
  python fetch_data.py

Check data/manifest.json to confirm what was fetched vs synthetic.
Note in the critic assessment which panels use real vs generated data.
If fetch_data.py fails on a source, retry up to 3 times before accepting
the synthetic fallback. Do not skip this step.

---

## Outputs

1. `CBG_pilot_figure.pdf` + `CBG_pilot_figure.png`
2. `CBG_model_explainer.ipynb`

Target aesthetic: Nature Communications. Every threshold, boundary, and
contour must be either literature-derived or explicitly labelled "assumed".
Run the full critic layer from AGENTS.md before saving. Revise until all
5 tests pass.

---

## Figure layout

Total: 170 mm wide × 150 mm tall (double-column, two-row).

Top row (55 mm): Panels A, B, C side by side.
  A: 35% width.  B: 35%.  C: 30%.

Bottom row (85 mm): Panel D — full width.
  Internally: thin time-course header strip (15 mm) + 4×5 spatial grid (65 mm).

Panel labels A B C D bold 8pt top-left. No figure title. No outer box.

---

## Panel A: Vascular-constrained drug delivery

**Claim:** "Drug entry without BBB opening is vessel-proximal and spatially
uniform; CBG-mediated opening creates discrete high-concentration hotspots
restricted to vessels within the expressing neuronal population."

### Geometry — use real data if available

Load `data/vasculature/minivess_slice.tif` if present.
Extract a 2×2 mm sub-region from the image centred on a capillary-dense patch.
Use the actual vessel geometry from the image as the spatial domain.

If MiniVess is absent, load `data/vasculature/vessel_graph.json`.
If both absent, generate a synthetic fractal tree (Murray's law,
r_d = r_p / 2^(1/3), tortuosity ±10° per 50 µm, target density
400–600 vessels/mm², Tsai et al. 2009 PNAS).

Whichever source is used: extract vessel centrelines as a binary mask
on a 10 µm grid. This mask is used for the line-source PDE.

### Layout: three sub-panels stacked vertically (labels i ii iii)

Shared viridis colourbar on right side.
200 µm scale bar bottom-left of sub-panel i.
Sub-panel labels i ii iii: italic, 7pt, top-left of each sub-panel.
Separate colourbar for sub-panel ii (expression density, orange).

### Sub-panel i — Systemic only (intact BBB)

Solve 2D PDE to steady state (finite difference, dx=10 µm, 60 min):
  dC/dt = D_tissue * ∇²C - k_e * C + P_intact * (C_blood - C) * vessel_mask(x,y)

vessel_mask: 1 on vessel pixels (Gaussian-smoothed, sigma=10 µm), 0 elsewhere.
Parameters from AGENTS.md: D_tissue=3e-4, k_e=0.05, P_intact=1e-4, C_blood=1.0.

Expected: concentration peaks at vessels (~0.01 normalised), falls to ~5%
within 200 µm. Off-vessel tissue near-zero.

### Sub-panel ii — Expression mask

Do NOT show drug concentration here.
If `data/enhancers/Prox1_ISH_*.jpg` is present:
  Load image, apply Gaussian blur (sigma=20px), threshold at 60th percentile,
  crop to the same 2×2 mm spatial domain as sub-panel i.
  Colour = expression intensity (orange colourmap, 0–1 normalised).
If absent: draw a patchy, non-Gaussian expression domain:
  Use 3–5 overlapping blobs of random size (r = 0.1–0.3 mm) with random
  centres within the patch. Sum them with random weights. Do NOT use a
  single centred Gaussian — this was the failure in the previous figure.
Overlay vessel centrelines as thin white lines (0.3pt).

### Sub-panel iii — CBG active

P(x,y) = P_open where (expression_mask > 0.5) AND (vessel_mask > 0.1)
P(x,y) = P_intact elsewhere.
Solve same PDE as sub-panel i. Show steady-state C heatmap (viridis).

Key visual: high-concentration hotspots appear ONLY at vessel pixels that
are also within the expression zone. Vessels outside the expression zone
remain at systemic-only baseline.

### Verification checks for Panel A
- Hotspot concentration / systemic-only background ratio: must be 10–100x.
  If outside range adjust P_open.
- Off-mask vessel concentration must match sub-panel i (no elevation).
- Expression mask must NOT be a single smooth Gaussian — if it is, regenerate
  with the patchy multi-blob approach above.

---

## Panel B: Anatomical specificity — coronal hippocampal section

**Claim:** "CBG activation with hippocampus-specific enhancer expression
produces bilateral drug accumulation in hippocampal grey matter that respects
white matter boundaries, at a systemically sub-therapeutic dose."

### Geometry — use real data if available

Strategy A: Load `data/atlas/P56_Annotation_downsample2.nii.gz`.
  Extract coronal slice at z ≈ 332 (Bregma -2.0 mm in CCF v3 at 25 µm/vox).
  Use structure IDs from `data/atlas/structure_ids.json` to generate
  binary masks for: CA1 (382), CA3 (463), DG (726), corpus callosum (776),
  thalamus (549), cortex (315), lateral ventricle (73).
  Convert masks to polygon outlines using skimage.measure.find_contours.

Strategy B: If NIfTI absent, load
  `data/atlas/hippocampus_polygons_Bregma-2mm.json`.
  Use the polygon coordinates directly to fill anatomical regions.

The brain outline must be shaped like a CORONAL section (butterfly shape
with hippocampi flanking central structures), NOT a dorsal view oval.
If the previous figure generated a dorsal-view shape, this was wrong.
A coronal section at Bregma -2.0 mm has approximate dimensions:
  width ~8 mm, height ~6 mm, widest at the hippocampal level.

### Tissue compartment properties

Grey matter (cortex, CA1, CA3, DG, thalamus):  D = 3e-4 cm²/min
White matter (corpus callosum, fimbria):         D = 5e-5 cm²/min
CSF (lateral ventricles):                        D = 1e-2 cm²/min
Source: Vorisek & Sykova 1997; white matter D is 10x lower than grey.

BBB is open ONLY in bilateral CA1 and DG.
Drug enters from a uniform source within those regions (k_source in those voxels).
Solve 2D steady-state diffusion with spatially varying D.

### Visual requirements

- Anatomical region boundaries: thin dark lines 0.5pt.
- Colourmap: viridis for drug concentration.
- White matter must be visibly DARKER than hippocampal grey matter —
  if it is not, check that the D values are applied correctly.
- Bilateral symmetry should be apparent but with ±5% noise added to
  the drug concentration to avoid a perfectly mirrored appearance.
- Add subregion labels CA1, DG, CC in 6pt italic.
- 1 mm scale bar bottom-left.
- Drug concentration outside BBB-open zones should be near-zero
  (below 10% of the CA1 peak).

---

## Panel C: Single-event time course — the feedback cycle

**Claim:** "A single seizure event triggers the full CBG feedback cycle:
pathological activity drives c-Fos induction, BBB opener accumulation,
drug entry, and seizure termination — all within a ~60 min window;
physiological activity does not cross the opening threshold."

Two stacked sub-panels sharing x-axis, time 0–90 min.
Mark the 5 Panel D time points as vertical dashed lines: t=0, 15, 30, 45, 70 min.
Label them "D:t1" through "D:t5" in 6pt above the top sub-panel.

### Top sub-panel: c-Fos fold-change A(t)

Two traces:
  Red solid:    Seizure + CBG condition
  Blue dashed:  Physiological condition
Shaded orange region showing window where P(t) > P_thresh (BBB open).
Horizontal dashed line at A_thresh=3.0 labelled "BBB-opening threshold (A_thresh=3.0)".
Y-axis: "c-Fos fold-change (A)", linear scale 0–35.

### Bottom sub-panel: Drug concentration C(t)

Two traces:
  Red solid:   CBG active
  Grey flat:   No CBG (C=0 everywhere)
Horizontal dashed lines:
  Green at therapeutic threshold — label "Therapeutic threshold (assumed)"
  Orange at side-effect threshold — label "Side-effect threshold (assumed)"
Y-axis: "Drug conc. (C, norm.)", linear scale.

### ODE system

dA/dt = k_A * S(t) - lambda_A * A
dP/dt = k_P * max(A - A_thresh, 0) - lambda_P * P
dC/dt = k_C * max(P - P_thresh, 0) * (C_blood - C) - lambda_C * C
dS/dt = -k_suppress * C * S    ← drug terminates seizure via this feedback term

Parameters (all to appear in notebook Cell 4 with citations):
  k_A = 2.0 /min          Tullai et al. 2007, Mol Cell Biol
  lambda_A = 0.05 /min    c-Fos mRNA t½ ~14 min; Tullai et al. 2007
  A_thresh = 3.0          physiological ceiling; Bhatt et al. 2020, Nat Neurosci
  k_P = 0.1 /min          transcription-translation lag ~20–30 min (estimated)
  lambda_P = 0.015 /min   protein t½ ~45 min (estimated, CCL2 analogue)
  P_thresh = 0.3          ← FREE PARAMETER — label "assumed" in figure
  k_C = 0.5               effective permeability when BBB open
  lambda_C = 0.05 /min    brain elimination; Vendel et al. 2019
  k_suppress = 0.8        ← FREE PARAMETER — label "assumed" in figure
  C_blood = 1.0           normalised systemic dose

Seizure input:       S(t) = 15 for t ∈ [10, 40] min, subject to suppression.
Physiological input: S(t) = 3  for t ∈ [10, 25] min.

Linear axes only. No log scale.

### Verification checks for Panel C
- Seizure condition: C(t) must cross therapeutic threshold.
- Physiological condition: C(t) must stay below therapeutic threshold.
- Both conditions: C returns to near-zero by t=90 min.
- 5 vertical dashed lines must align exactly with Panel D time points.

---

## Panel D: Spatial evolution — 4 rows × 5 time snapshots

**Claim:** "Without CBG, seizure-driven c-Fos spreads bilaterally through the
hippocampus while drug concentration remains zero. With CBG, drug entry
coincides with peak c-Fos and accelerates seizure resolution — visible as
drug hotspots that spatially track and then outlast the c-Fos activity."

### Layout

Full panel width. 4 rows × 5 columns of hippocampal spatial heatmaps.
Each small panel: same coronal geometry as Panel B (use same polygon/mask).
Approximate size per panel: 28 mm wide × 14 mm tall.

The 5 COLUMNS are 5 time snapshots (same times as Panel C vertical lines):
  t₁ = 0 min     (baseline — pre-seizure)
  t₂ = 15 min    (early seizure, c-Fos rising)
  t₃ = 30 min    (seizure peak)
  t₄ = 45 min    (resolution, with vs without CBG diverges)
  t₅ = 70 min    (post-seizure)

Column headers above grid (8pt): "t=0", "t=15 min", "t=30 min", "t=45 min", "t=70 min".

The 4 ROWS are:
  Row 1: c-Fos spatial field — NO CBG condition
  Row 2: Drug concentration — NO CBG condition
  Row 3: c-Fos spatial field — CBG condition
  Row 4: Drug concentration — CBG condition

Row labels on left margin (7pt, rotated 90°):
  "c-Fos (no CBG)", "Drug (no CBG)", "c-Fos (CBG)", "Drug (CBG)"

### Spatial field generation

The c-Fos spatial field at time t is:
  F(x,y,t) = A(t) * G(x,y,t)

where G is a spreading Gaussian seeded at left CA1 (seizure onset):
  G(x,y,t) = exp(-r² / (2 * sigma(t)²))
  sigma(t) = sigma_0 + v_spread * t
  sigma_0 = 0.3 mm, v_spread = 0.05 mm/min
  (seizure propagation speed ~0.3 mm/min in hippocampus;
   Trevelyan et al. 2006, J Neurosci — use this citation)

IMPORTANT: The seizure is bilateral. Mirror the Gaussian to right CA1 with a
delay of 5 min (seizure propagates via commissural fibres). The bilateral
spreading makes the figure visually distinct from a unilateral model and is
anatomically accurate for hippocampal seizures.

IMPORTANT: Restrict F(x,y,t) to grey matter only. Values in white matter
(corpus callosum) and ventricles must be set to zero. This was wrong in the
previous figure — white matter should be dark in the c-Fos rows.

For the NO CBG condition: A(t) from the ODE with k_C=0 (drug never enters,
seizure runs its natural course without pharmacological termination).
For the CBG condition: A(t) from the full ODE with k_suppress active.

Drug concentration field D(x,y,t):
  Drug source at time t = k_C * max(P(t) - P_thresh, 0) * C_blood
  applied only within hippocampal grey matter voxels (CA1 + DG bilaterally).
  Solve 2D diffusion at each time point using the same compartmental D
  values as Panel B (grey/white matter diffusion contrast).
  NO CBG: drug source = 0 everywhere at all times → all-dark panels.

### Colour scales

Rows 1 & 3 (c-Fos): 'hot' colourmap (black → red → yellow → white).
  Shared colourscale: 0 to max(A(t)) across ALL panels in rows 1 AND 3.
  Single colourbar labelled "c-Fos (A)" placed right of column 5, rows 1–3.

Rows 2 & 4 (drug): viridis colourmap.
  Shared colourscale: 0 to max(C_CBG) across ALL panels in rows 2 AND 4.
  CBG row sets the max; no-CBG panels will appear uniformly dark (correct).
  Single colourbar labelled "Drug conc. (norm.)" placed right of column 5, rows 2–4.

### Key visual results that MUST be visible

1. Row 1 (c-Fos, no CBG): Activity spreads from left CA1 at t₂, becomes
   bilateral by t₃, persists at t₄ and t₅ (no termination without drug).
   White matter (corpus callosum) remains dark throughout.

2. Row 2 (drug, no CBG): Uniformly near-zero (all dark) at ALL time points.
   This is intentional — it demonstrates the failure of systemic delivery.
   Do NOT show a colourbar that makes this look like something is happening.

3. Row 3 (c-Fos, CBG): Identical to Row 1 at t₂ and t₃. At t₄ the signal
   is VISIBLY LESS INTENSE than Row 1 — the seizure is terminating because
   drug has entered. At t₅ the field should be clearly darker than Row 1.
   If rows 1 and 3 look identical at t₄ and t₅, increase k_suppress.

4. Row 4 (drug, CBG): Near-zero at t₁ and t₂. Bilateral hotspots appear
   in hippocampal grey matter at t₃ (peak BBB-opening window). Hotspots
   are spatially co-localised with the c-Fos pattern in Row 3 (NOT diffuse).
   At t₄ hotspots persist; at t₅ they begin decaying.
   If the drug fills the whole brain outline uniformly, the diffusion
   model is not applying the tissue compartment mask correctly — fix it.

### Formatting

- No tick labels on any small panel.
- No axes on small panels: just the heatmap + thin white anatomical outlines (0.3pt).
- Column gap: 1.5 mm. Row gap: 2 mm.
- Thicker gap (3 mm) + thin horizontal rule between rows 2 and 3.
- Labels "No CBG" and "CBG active" centred in the gap, 7pt bold.
- Colourbar height spans rows 1–3 (c-Fos) and rows 2–4 (drug) respectively.

---

## Critic layer (mandatory before saving)

After generating, write the full assessment block:

```
CRITIC ASSESSMENT:
Test 1 (Single-claim):
  Panel A: [one falsifiable sentence] — PASS/FAIL
  Panel B: [one falsifiable sentence] — PASS/FAIL
  Panel C: [one falsifiable sentence] — PASS/FAIL
  Panel D: [one falsifiable sentence] — PASS/FAIL

Test 2 (Grounded parameters):
  A_thresh=3.0: Bhatt et al. 2020
  P_thresh=0.3: FREE PARAMETER — labelled ASSUMED ✓/✗
  k_suppress=0.8: FREE PARAMETER — labelled ASSUMED ✓/✗
  P_open/P_intact ratio ~50x: FUS-BBBO literature ✓/✗
  v_spread=0.05 mm/min: Trevelyan et al. 2006 ✓/✗
  [list any others]

Test 3 (Biological realism):
  Panel A: real vasculature used (MiniVess) or synthetic? State which.
  Panel A sub-ii: patchy multi-blob mask (not Gaussian)? PASS/FAIL
  Panel B: coronal butterfly shape (not dorsal oval)? PASS/FAIL
  Panel B: white matter darker than grey? PASS/FAIL
  Panel D: white matter dark in c-Fos rows? PASS/FAIL
  Panel D: bilateral hippocampal seizure spread visible? PASS/FAIL
  Panel D: drug restricted to hippocampus (not whole brain)? PASS/FAIL

Test 4 (Mock reviewer — write as hostile but fair):
  [two sentences]

Test 5 (Story coherence A→B→C→D):
  [one paragraph]

DATA SOURCES USED:
  Vasculature: [MiniVess real / VesselGraph / synthetic]
  Atlas: [Allen CCF NIfTI / polygon fallback]
  Enhancers: [Prox1 ISH real / synthetic mask]

OVERALL: PASS / FAIL (list panels to revise)
```

If OVERALL is FAIL: revise failing panels, regenerate, re-run assessment.
The most common failure modes from the previous iteration were:
  - Panel B drawn as a dorsal view instead of coronal (fix: use polygon coords)
  - Rows 1 and 3 of Panel D indistinguishable (fix: increase k_suppress)
  - Drug in Panel D Row 4 filling whole brain (fix: apply tissue mask)
  - Expression mask in Panel A sub-ii a single Gaussian (fix: multi-blob)
Address these explicitly.

---

## Notebook: CBG_model_explainer.ipynb

Cell 1: Title + 3-sentence plain-language description.
Cell 2: Physical model. LaTeX for all ODEs. Explain the dS/dt feedback term
  explicitly: "The drug suppresses seizure activity at rate k_suppress,
  creating a negative feedback loop that terminates the seizure."
  Explain the bilateral seizure propagation model used in Panel D.
Cell 3: Imports + pinned versions + call to fetch_data.py if data/ absent.
Cell 4: TUNABLE PARAMETERS with ipywidgets sliders for:
  P_open (1e-4 to 1e-2), A_thresh (1 to 10), P_thresh (0.1 to 1.0),
  k_suppress (0 to 3), v_spread (0.01 to 0.2 mm/min).
  Wire to regenerate Panel C time course and Panel D grid live.
  Inline comment for each: value used | source or ASSUMED.
Cell 5–8: Individual panels with inline checks.
Cell 9: Assemble and save figure. Print CRITIC ASSESSMENT block.
Cell 10: Parameter table (value | source | ASSUMED flag for free parameters).
Cell 11: Limitations: tortuosity, perivascular flow, protein binding,
  inflammation-driven permeability, single seizure type, no pharmacokinetic
  model for drug distribution. Next steps for Phase 2 model.

---

## Final checks before saving

- [ ] A1: Real vasculature loaded (MiniVess) OR synthetic noted in assessment
- [ ] A2: Expression mask is patchy/multi-blob, NOT a single Gaussian
- [ ] A3: CBG hotspot 10–100x systemic-only
- [ ] A4: Off-mask vessels at systemic-only level
- [ ] B1: Coronal butterfly shape, NOT dorsal oval
- [ ] B2: White matter visibly darker than hippocampal grey matter
- [ ] B3: Drug bilateral within ±10%
- [ ] C1: Seizure C(t) crosses therapeutic threshold
- [ ] C2: Physiological C(t) stays below threshold
- [ ] C3: 5 vertical dashed lines present and match Panel D times exactly
- [ ] D1: Row 1 and Row 3 visibly diverge at t₄ and t₅
- [ ] D2: Row 2 uniformly dark (drug never enters without CBG)
- [ ] D3: Row 4 drug hotspots co-localised with Row 3 c-Fos at t₃
- [ ] D4: Drug in Row 4 restricted to hippocampus, not whole brain
- [ ] D5: White matter dark in all c-Fos rows
- [ ] D6: Bilateral hippocampal spread visible in both c-Fos rows
- [ ] D7: Shared colourscale within rows 1&3; shared within rows 2&4
- [ ] N1: P_thresh and k_suppress labelled ASSUMED in figure and notebook
- [ ] N2: fetch_data.py called in Cell 3 if data/ absent
- [ ] CRITIC: All 5 tests pass
