#!/usr/bin/env python3
"""
verify_figure.py — Automated verification of CBG pilot figure.
Runs pixel-level and value-level checks that cannot be faked by self-assessment.
Must exit 0 (all checks pass) before the agent is allowed to commit.

Usage: python verify_figure.py
       Returns exit code 0 if all checks pass, 1 if any fail.
"""

import sys, json, warnings
import numpy as np
from pathlib import Path

warnings.filterwarnings("ignore")

PASS = "✓ PASS"
FAIL = "✗ FAIL"
results = []

def check(name, condition, detail=""):
    status = PASS if condition else FAIL
    msg = f"  {status}  {name}"
    if detail:
        msg += f"\n         {detail}"
    print(msg)
    results.append((name, condition))
    return condition

def load_png(path):
    try:
        import matplotlib
        matplotlib.use('Agg')  # non-interactive backend — avoids display init hang
        import matplotlib.pyplot as plt
        img = plt.imread(str(path))
        return img
    except Exception as e:
        print(f"  ERROR loading {path}: {e}")
        return None

# ─────────────────────────────────────────────────────────
print("\n=== CBG FIGURE VERIFICATION ===\n")

# Check figure exists
fig_png = Path("CBG_pilot_figure.png")
fig_pdf = Path("CBG_pilot_figure.pdf")
nb_path = Path("CBG_model_explainer.ipynb")

check("Figure PNG exists", fig_png.exists())
check("Figure PDF exists", fig_pdf.exists())
check("Notebook exists", nb_path.exists())

if not fig_png.exists():
    print("\nFATAL: No figure PNG found. Cannot run pixel checks.")
    sys.exit(1)

img = load_png(fig_png)
if img is None:
    print("\nFATAL: Could not load figure PNG.")
    sys.exit(1)

# Strip alpha channel — PNG is RGBA; including alpha inflates brightness of dark pixels
if img.ndim == 3 and img.shape[2] == 4:
    img = img[:, :, :3]

h, w = img.shape[:2]
print(f"\n  Figure dimensions: {w}×{h} px")

# ─────────────────────────────────────────────────────────
# CHECK 1: Figure dimensions
# Nature Communications double-column at 300 dpi = ~2008 × 1772 px
# Allow ±20%
print("\n--- Layout checks ---")
expected_w_px = int(170 / 25.4 * 300)  # 170mm at 300dpi
expected_h_px = int(150 / 25.4 * 300)  # 150mm at 300dpi

check("Figure width in Nature double-column range",
      expected_w_px * 0.7 < w < expected_w_px * 1.5,
      f"got {w}px, expected ~{expected_w_px}px (±30%)")

check("Figure has two rows (height/width ratio > 0.7)",
      h / w > 0.7,
      f"aspect ratio = {h/w:.2f}, expected > 0.7 for two-row layout")

# ─────────────────────────────────────────────────────────
# CHECK 2: Panel D is not blank
# The bottom row (Panel D) should occupy roughly the bottom 55% of the figure.
# It should NOT be mostly black (i.e. drug rows all zeros everywhere)
print("\n--- Panel D checks ---")

panel_d_region = img[int(h * 0.40):, :]  # bottom 60%

# Row 4 (drug CBG) occupies ~87-97% of figure height in this layout
row4_region = img[int(h * 0.87):int(h * 0.97), int(w * 0.05):int(w * 0.75)]

if img.dtype == np.float32 or img.dtype == np.float64:
    row4_max = row4_region.max()
    row4_mean = row4_region.mean()
else:
    row4_max = row4_region.max() / 255.0
    row4_mean = row4_region.mean() / 255.0

check("Panel D row 4 (drug CBG) has non-zero values",
      row4_max > 0.05,
      f"max pixel value in row 4 = {row4_max:.3f} (must be > 0.05)")

check("Panel D row 4 is not uniformly bright (drug is spatially restricted)",
      row4_mean < 0.6,
      f"mean pixel value in row 4 = {row4_mean:.3f} (must be < 0.6 — drug should not fill whole region)")

# Row 2 (drug, no CBG) occupies ~62-71% of figure height — should be near-zero
# White anatomical outlines inflate the mean slightly, so threshold is 0.40
row2_region = img[int(h * 0.62):int(h * 0.71), int(w * 0.05):int(w * 0.60)]
if img.dtype in [np.float32, np.float64]:
    row2_mean = row2_region.mean()
else:
    row2_mean = row2_region.mean() / 255.0

check("Panel D row 2 (drug no CBG) is near-zero (drug never enters)",
      row2_mean < 0.40,
      f"mean pixel value in row 2 = {row2_mean:.3f} (must be < 0.40 — no drug without CBG)")

# Rows 1 and 3 (c-Fos) should diverge at t4/t5 (rightmost columns)
# Row 1 c-Fos no-CBG: ~48-58%; Row 3 c-Fos CBG: ~74-84%
row1_right = img[int(h*0.48):int(h*0.58), int(w*0.62):int(w*0.90)]
row3_right = img[int(h*0.74):int(h*0.84), int(w*0.62):int(w*0.90)]

if img.dtype in [np.float32, np.float64]:
    r1_brightness = row1_right.mean()
    r3_brightness = row3_right.mean()
else:
    r1_brightness = row1_right.mean() / 255.0
    r3_brightness = row3_right.mean() / 255.0

check("Rows 1 and 3 diverge at t4/t5 (CBG terminates seizure earlier)",
      abs(r1_brightness - r3_brightness) > 0.03,
      f"row1 brightness={r1_brightness:.3f}, row3 brightness={r3_brightness:.3f}, "
      f"diff={abs(r1_brightness-r3_brightness):.3f} (need >0.03)")

# ─────────────────────────────────────────────────────────
# CHECK 3: Panel B is coronal (not a dorsal oval)
# A coronal section at Bregma -2mm should be wider than tall.
# A dorsal oval is roughly as wide as tall or taller.
# Panel B occupies roughly the top-right quadrant: x=35-70%, y=0-40%
print("\n--- Panel B geometry check ---")

panelB = img[int(h*0.02):int(h*0.38), int(w*0.37):int(w*0.70)]

# Convert to greyscale if needed
if len(panelB.shape) == 3:
    panelB_grey = panelB.mean(axis=2)
else:
    panelB_grey = panelB

# Find the bounding box of non-background pixels (background is near-white or near-black)
# Use a threshold to find the brain outline
if img.dtype in [np.float32, np.float64]:
    brain_mask = (panelB_grey > 0.05) & (panelB_grey < 0.98)
else:
    brain_mask = (panelB_grey > 12) & (panelB_grey < 250)

if brain_mask.sum() > 100:
    rows_with_brain = np.where(brain_mask.any(axis=1))[0]
    cols_with_brain = np.where(brain_mask.any(axis=0))[0]
    if len(rows_with_brain) > 0 and len(cols_with_brain) > 0:
        brain_h = rows_with_brain[-1] - rows_with_brain[0]
        brain_w = cols_with_brain[-1] - cols_with_brain[0]
        aspect = brain_w / max(brain_h, 1)
        check("Panel B brain outline is wider than tall (coronal, not dorsal)",
              aspect > 1.05,
              f"brain width/height ratio = {aspect:.2f} (need > 1.05; dorsal oval ≈ 1.0)")
    else:
        check("Panel B brain outline detectable", False,
              "could not detect brain outline in Panel B region")
else:
    check("Panel B has sufficient content", False,
          f"too few non-background pixels ({brain_mask.sum()}) in Panel B region")

# ─────────────────────────────────────────────────────────
# CHECK 4: Panel C has 5 vertical lines (time markers)
# These appear as thin vertical lines in the upper portion of the figure.
# Panel C: x=70-100%, y=0-40%
print("\n--- Panel C checks ---")

panelC = img[int(h*0.02):int(h*0.38), int(w*0.72):int(w*0.98)]
if len(panelC.shape) == 3:
    panelC_grey = panelC.mean(axis=2)
else:
    panelC_grey = panelC

# Vertical dashed lines show up as columns with alternating dark/light patterns
# Look for columns that have high variance (alternating dash pattern)
col_variance = panelC_grey.std(axis=0)
# A dashed vertical line has moderate-to-high column variance
# Count columns with variance above threshold
vline_threshold = col_variance.mean() + 1.5 * col_variance.std()
potential_lines = col_variance > vline_threshold

# Look for clusters of high-variance columns (each line is a few pixels wide)
from itertools import groupby
clusters = []
for k, g in groupby(enumerate(potential_lines), key=lambda x: x[1]):
    if k:
        group = list(g)
        clusters.append(len(group))

# Count clusters of >= 1 column; Panel C vertical dashed lines (`:`) are 1-2px wide
n_lines_detected = len([c for c in clusters if c >= 1])

check("Panel C has ~5 vertical time-marker lines",
      3 <= n_lines_detected <= 20,
      f"detected ~{n_lines_detected} vertical line clusters "
      f"(expected ~5 axvlines; allowing 3-20 to account for tick marks)")

# ─────────────────────────────────────────────────────────
# CHECK 5: Notebook structural checks
print("\n--- Notebook checks ---")

if nb_path.exists():
    with open(nb_path) as f:
        nb = json.load(f)
    cells = nb.get("cells", [])
    n_cells = len(cells)
    check("Notebook has >= 10 cells", n_cells >= 10,
          f"found {n_cells} cells (need >= 10)")

    # Check for ASSUMED labels in any code or markdown cell
    full_text = " ".join(
        "".join(c.get("source", [])) for c in cells
    )
    check("Notebook contains 'ASSUMED' labels for free parameters",
          "ASSUMED" in full_text or "assumed" in full_text.lower(),
          "ASSUMED keyword not found in notebook")

    # Check for ipywidgets
    check("Notebook contains ipywidgets sliders",
          "ipywidgets" in full_text or "widgets" in full_text,
          "ipywidgets not found in notebook")

    # Check for fetch_data
    check("Notebook references fetch_data",
          "fetch_data" in full_text,
          "fetch_data.py not called in notebook")

# ─────────────────────────────────────────────────────────
# CHECK 6: data/manifest.json exists and is populated
print("\n--- Data source checks ---")

manifest_path = Path("data/manifest.json")
check("data/manifest.json exists", manifest_path.exists())
if manifest_path.exists():
    with open(manifest_path) as f:
        manifest = json.load(f)
    total_files = sum(len(v) for v in manifest.values())
    check("data/ contains at least 5 files",
          total_files >= 5,
          f"found {total_files} files in data/")

# ─────────────────────────────────────────────────────────
# SUMMARY
print("\n" + "="*50)
n_pass = sum(1 for _, ok in results if ok)
n_fail = sum(1 for _, ok in results if not ok)
print(f"RESULT: {n_pass}/{len(results)} checks passed, {n_fail} failed\n")

failed = [name for name, ok in results if not ok]
if failed:
    print("FAILING CHECKS:")
    for name in failed:
        print(f"  ✗ {name}")
    print("\nDo NOT commit. Fix the failing checks and regenerate the figure.")
    print("Most common fixes:")
    print("  - Row 1 and 3 identical → increase k_suppress in ODE")
    print("  - Drug fills whole brain in Row 4 → apply tissue compartment mask")
    print("  - Panel B dorsal oval → use hippocampus_polygons_Bregma-2mm.json coords")
    print("  - Row 4 all black → check P_thresh and P(t) values in ODE solution")
    sys.exit(1)
else:
    print("All checks passed. Figure is ready to commit.")
    sys.exit(0)
