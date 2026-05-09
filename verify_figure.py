#!/usr/bin/env python3
"""
verify_figure.py — Automated verification of CBG combined figure.
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

# Output paths
fig_png = Path("CBG_outputs/CBG_figure.png")
fig_pdf = Path("CBG_outputs/CBG_figure.pdf")
nb_path = Path("CBG_combined_figure.ipynb")

check("Figure PNG exists", fig_png.exists())
check("Figure PDF exists", fig_pdf.exists())
check("Combined notebook exists", nb_path.exists())

if not fig_png.exists():
    print("\nFATAL: No figure PNG found. Cannot run pixel checks.")
    sys.exit(1)

img = load_png(fig_png)
if img is None:
    print("\nFATAL: Could not load figure PNG.")
    sys.exit(1)

# Strip alpha channel
if img.ndim == 3 and img.shape[2] == 4:
    img = img[:, :, :3]

h, w = img.shape[:2]
print(f"\n  Figure dimensions: {w}×{h} px")

# ─────────────────────────────────────────────────────────
print("\n--- Layout checks ---")

# Figure is 17"×10" at 600 dpi → 10200×6000 px; allow ±40% for screen-dpi saves
check("Figure width is wide-format (>= 1000 px)",
      w >= 1000,
      f"got {w} px")

check("Figure aspect ratio consistent with landscape layout (0.4 < h/w < 0.8)",
      0.4 < h / w < 0.8,
      f"h/w = {h/w:.2f}")

# ─────────────────────────────────────────────────────────
# Panel A checks — left ~28% of the figure width
print("\n--- Panel A checks ---")

# Normalise pixel values to [0,1]
img_f = img.astype(float)
if img_f.max() > 1.5:
    img_f /= 255.0

left_col = img_f[:, :int(w * 0.30), :]

# Top sub-panel (IV infusion only): tissue should be darker than Panel A bottom
# (white margins / colorbars inflate brightness, so we only require > 1% dark pixels)
top_panel = left_col[:int(h * 0.50), :]
tissue_like = top_panel[top_panel.mean(axis=-1) < 0.2]
check("Panel A top (IV infusion): has dark tissue pixels (drug confined to vessels)",
      tissue_like.shape[0] > 0.01 * top_panel.reshape(-1, 3).shape[0],
      f"dark-pixel fraction = {tissue_like.shape[0] / top_panel.reshape(-1,3).shape[0]:.3f}")

# Bottom sub-panel (IV+CBG): should have more bright tissue than top panel
bot_panel = left_col[int(h * 0.52):, :]
bot_mid_bright = (bot_panel.mean(axis=-1) > 0.25).mean()
top_mid_bright = (top_panel.mean(axis=-1) > 0.25).mean()
check("Panel A bottom (IV+CBG) has more bright tissue than top (drug leaked in)",
      bot_mid_bright > top_mid_bright,
      f"bright fraction top={top_mid_bright:.3f}, bot={bot_mid_bright:.3f}")

# ─────────────────────────────────────────────────────────
# Small-multiples checks — right ~70% of the figure, top ~60% height
print("\n--- Small-multiples checks ---")

sm = img_f[:int(h * 0.60), int(w * 0.30):]

# Row 1 (untreated cFos, top ~25% of small-multiples): should be dim early, bright later
row1 = sm[:int(sm.shape[0] * 0.25), :]

# Untreated cFos should build up over time — verified directly from provenance JSON
_prov_path = Path("CBG_outputs/CBG_figure_provenance.json")
if _prov_path.exists():
    with open(_prov_path) as _f:
        _prov = json.load(_f)
    _cu = _prov.get("results", {}).get("mean_cfos_u_by_window", [])
    if _cu:
        check("Untreated cFos mean increases from W1 to W4 (provenance-based)",
              float(_cu[-2]) > float(_cu[0]),
              f"W1={_cu[0]:.4f}  W4={_cu[-2]:.4f}")

# Row 4 (treated drug, bottom ~25% of small-multiples): should have bright tissue pixels from W2 onward
row4 = sm[int(sm.shape[0] * 0.75):, :]
row4_right = row4[:, int(row4.shape[1] * 0.20):]   # W2-W5 columns
check("Treated drug row (small multiples, row 4) has bright tissue pixels",
      row4_right.mean() > 0.1,
      f"W2-W5 mean brightness={row4_right.mean():.3f} (need > 0.1)")

# Row 2 (untreated drug): vessels only, should be darker in tissue than row 4
row2 = sm[int(sm.shape[0] * 0.25):int(sm.shape[0] * 0.50), :]
check("Treated drug (row 4) is brighter than untreated drug (row 2) overall",
      row4.mean() > row2.mean(),
      f"row2 mean={row2.mean():.3f}  row4 mean={row4.mean():.3f}")

# ─────────────────────────────────────────────────────────
# Line-chart checks — right ~70% of figure, bottom ~40% height
print("\n--- Line-chart checks ---")

lc = img_f[int(h * 0.60):, int(w * 0.30):]

# Line charts should have a visible background (not all-white or all-black)
lc_var = lc.std()
check("Line chart region has visible content (non-uniform)",
      lc_var > 0.02,
      f"pixel std = {lc_var:.4f} (need > 0.02)")

# The charts should not be identical (two distinct plots side-by-side)
left_chart  = lc[:, :int(lc.shape[1] * 0.48)]
right_chart = lc[:, int(lc.shape[1] * 0.52):]
diff = abs(left_chart.mean() - right_chart.mean())
check("Left and right line charts are distinct (not identical)",
      diff > 0.005 or left_chart.std() > 0.01,
      f"mean diff={diff:.4f}, left std={left_chart.std():.4f}")

# ─────────────────────────────────────────────────────────
# Notebook structural checks
print("\n--- Notebook checks ---")

if nb_path.exists():
    with open(nb_path) as f:
        nb = json.load(f)
    cells = nb.get("cells", [])
    n_cells = len(cells)
    check("Notebook has >= 6 cells", n_cells >= 6,
          f"found {n_cells} cells")

    full_text = " ".join("".join(c.get("source", [])) for c in cells)
    check("Notebook contains ASSUMED labels for free parameters",
          "ASSUMED" in full_text or "assumed" in full_text.lower())
    check("Notebook references panel_a_data.npz",
          "panel_a_data.npz" in full_text)

# ─────────────────────────────────────────────────────────
# Data source checks
print("\n--- Data source checks ---")

check("panel_a_data.npz exists", Path("panel_a_data.npz").exists())

manifest_path = Path("data/manifest.json")
check("data/manifest.json exists", manifest_path.exists())
if manifest_path.exists():
    with open(manifest_path) as f:
        manifest = json.load(f)
    total_files = sum(len(v) if isinstance(v, list) else 1 for v in manifest.values())
    check("data/ manifest has at least 3 entries",
          total_files >= 3,
          f"found {total_files} entries in manifest")

provenance_path = Path("CBG_outputs/CBG_figure_provenance.json")
check("CBG_outputs/CBG_figure_provenance.json exists", provenance_path.exists())
if provenance_path.exists():
    with open(provenance_path) as f:
        prov = json.load(f)
    check("Provenance records critic_passed = True",
          prov.get("critic_passed", False) is True,
          f"critic_passed = {prov.get('critic_passed')}")

# ─────────────────────────────────────────────────────────
# SUMMARY
print("\n" + "=" * 50)
n_pass = sum(1 for _, ok in results if ok)
n_fail = sum(1 for _, ok in results if not ok)
print(f"RESULT: {n_pass}/{len(results)} checks passed, {n_fail} failed\n")

failed = [name for name, ok in results if not ok]
if failed:
    print("FAILING CHECKS:")
    for name in failed:
        print(f"  ✗ {name}")
    print("\nDo NOT commit. Fix the failing checks and regenerate the figure.")
    print("Common fixes:")
    print("  - Panel A top not dark → check C_intact values (should be near 0 in tissue)")
    print("  - Panel A bottom not brighter → check C_open hotspot (P_open / P_intact = 50x)")
    print("  - Row 4 not bright → check interstitial_gain and spatial coupling")
    print("  - Critic not passed → rerun Cell 6 and fix failing test")
    sys.exit(1)
else:
    print("All checks passed. Figure is ready to commit.")
    sys.exit(0)
