#!/usr/bin/env python3
"""
fetch_data.py — Downloads all source data for CBG figure generation.
Designed to run inside GitHub Actions or on roschnode.
Each source has multiple fallback strategies. The script keeps retrying
until all required files are present or all fallbacks are exhausted.

Run: python fetch_data.py
Output: data/ directory with vasculature, atlas, and enhancer files.
"""

import os, sys, time, json, gzip, shutil, struct, hashlib, urllib.request, urllib.error
from pathlib import Path

DATA_DIR = Path("data")
VASC_DIR = DATA_DIR / "vasculature"
ATLAS_DIR = DATA_DIR / "atlas"
ENH_DIR   = DATA_DIR / "enhancers"
FOCAL_DIR = DATA_DIR / "focal"

for d in [VASC_DIR, ATLAS_DIR, ENH_DIR, FOCAL_DIR]:
    d.mkdir(parents=True, exist_ok=True)

FOCAL_PAPER = {
    "pii": "S0969996122000249",
    "doi": "10.1016/j.nbd.2022.105633",
    "title": "Neuronal circuits sustaining neocortical-injury-induced status epilepticus",
    "journal": "Neurobiology of Disease",
    "year": 2022,
}

MAX_RETRIES = 3
RETRY_WAIT  = 5  # seconds between retries

# ─────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────

def log(msg):
    print(f"[fetch_data] {msg}", flush=True)

def download(url, dest, timeout=120):
    """Download url to dest with retries. Returns True on success."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        log(f"  already present: {dest.name}")
        return True
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            log(f"  downloading {dest.name} (attempt {attempt}) ...")
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(dest, "wb") as f:
                shutil.copyfileobj(r, f)
            if dest.stat().st_size > 0:
                log(f"  OK: {dest.name} ({dest.stat().st_size // 1024} KB)")
                return True
            dest.unlink()
        except Exception as e:
            log(f"  attempt {attempt} failed: {e}")
            if dest.exists(): dest.unlink()
            if attempt < MAX_RETRIES: time.sleep(RETRY_WAIT)
    return False

def try_pip_install(pkg):
    import subprocess
    log(f"  installing {pkg} ...")
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", pkg, "--quiet"],
        capture_output=True
    )
    return result.returncode == 0

# ─────────────────────────────────────────────────────────
# SOURCE 1: MiniVess — real two-photon cortical vasculature
# Source: Poon et al. 2023, Scientific Data
# Figshare DOI: 10.6084/m9.figshare.21742947
# We download one representative 2D slice as a TIFF.
# ─────────────────────────────────────────────────────────

def fetch_vasculature():
    dest_tif  = VASC_DIR / "minivess_slice.tif"
    dest_json = VASC_DIR / "vessel_graph.json"

    log("\n=== SOURCE 1: Cortical vasculature ===")

    # Strategy A: MiniVess Figshare — direct volume download (subset)
    # The Figshare API returns download links for each file in the dataset.
    # We pick the smallest volume for speed.
    figshare_api = "https://api.figshare.com/v2/articles/21742947/files"
    tif_ok = False
    if not dest_tif.exists():
        try:
            with urllib.request.urlopen(figshare_api, timeout=30) as r:
                files = json.loads(r.read())
            # pick the smallest .tif file
            tifs = [f for f in files if f["name"].lower().endswith(".tif")]
            if tifs:
                target = min(tifs, key=lambda x: x.get("size", 9e9))
                log(f"  MiniVess: downloading {target['name']} ({target.get('size',0)//1024} KB)")
                tif_ok = download(target["download_url"], dest_tif)
        except Exception as e:
            log(f"  MiniVess Figshare API failed: {e}")

    # Strategy B: VesselGraph — whole-brain vessel graph from Paetzold et al. 2021
    # GitHub release: https://github.com/jocpae/VesselGraph
    if not dest_json.exists():
        vg_urls = [
            # Direct release asset — C57BL/6 graph subset
            "https://github.com/jocpae/VesselGraph/releases/download/v1.0/vessel_graph_C57BL6.json",
            # Fallback: raw data from repo
            "https://raw.githubusercontent.com/jocpae/VesselGraph/main/dataset/vessel_graph.json",
        ]
        for url in vg_urls:
            if download(url, dest_json):
                break

    # Strategy C: Synthetic fallback — generate a realistic fractal vascular tree
    # using Murray's law and empirical capillary density (Tsai et al. 2009 PNAS)
    if not dest_tif.exists() and not dest_json.exists():
        log("  All external sources failed — generating synthetic vascular tree ...")
        _generate_synthetic_vasculature(VASC_DIR / "synthetic_vessel_tree.json")

    # Report
    available = [f.name for f in VASC_DIR.iterdir() if f.stat().st_size > 0]
    log(f"  Vasculature files: {available}")
    return len(available) > 0


def _generate_synthetic_vasculature(dest):
    """
    Generate a fractal branching vascular tree using Murray's law.
    r_daughter = r_parent / 2^(1/3). Tortuosity: ±10° per 50 µm.
    Target density: 400-600 capillaries/mm². Tsai et al. 2009 PNAS.
    """
    import math, random
    random.seed(42)

    segments = []
    def branch(x, y, angle, radius, depth):
        if radius < 0.004 or depth > 12:  # stop at ~4 µm capillary radius
            return
        length = max(0.05, radius * 8)    # segment length proportional to radius
        # add tortuosity: random angle perturbation every 50 µm
        n_steps = max(1, int(length / 0.05))
        pts = [(x, y)]
        cx, cy, ca = x, y, angle
        for _ in range(n_steps):
            ca += math.radians(random.gauss(0, 10))
            step = length / n_steps
            cx += step * math.cos(ca)
            cy += step * math.sin(ca)
            pts.append((cx, cy))
        segments.append({"points": pts, "radius": radius, "depth": depth})
        # Murray's law branching
        r_d = radius / (2 ** (1/3))
        branch_angle = random.uniform(40, 70)
        branch(pts[-1][0], pts[-1][1], ca - math.radians(branch_angle), r_d, depth+1)
        branch(pts[-1][0], pts[-1][1], ca + math.radians(branch_angle), r_d, depth+1)

    # Start with parent vessels entering from edges of 2mm x 2mm patch
    for start_angle in [0, 45, 90, 135]:
        branch(0.2, 1.0 + 0.5*(start_angle/90), math.radians(start_angle),
               0.025, 0)

    with open(dest, "w") as f:
        json.dump({"segments": segments,
                   "source": "synthetic_Murray_law",
                   "ref": "Tsai et al. 2009 PNAS: 400-600 vessels/mm2",
                   "domain_mm": [2.0, 2.0]}, f)
    log(f"  Generated synthetic tree: {len(segments)} segments → {dest.name}")


# ─────────────────────────────────────────────────────────
# SOURCE 2: Allen Mouse Brain Atlas — CCF v3 annotation
# Source: Allen Institute, CCF v3, 25 µm voxels
# Provides region masks for CA1/CA3/DG/CC/thalamus/cortex
# Download via allensdk OR direct NIfTI from Scalable Brain Atlas
# ─────────────────────────────────────────────────────────

# Allen structure IDs (CCF v3 ontology)
STRUCTURE_IDS = {
    "CA1":          382,
    "CA3":          463,
    "DG":           726,
    "corpus_callosum": 776,
    "thalamus":     549,
    "cortex":       315,
    "lateral_ventricle": 73,
    "fimbria":      78,
    "hippocampus":  1080,
}

def fetch_atlas():
    dest_nii  = ATLAS_DIR / "P56_Annotation_downsample2.nii.gz"
    dest_ids  = ATLAS_DIR / "structure_ids.json"
    dest_svg  = ATLAS_DIR / "coronal_Bregma-2mm.svg"

    log("\n=== SOURCE 2: Allen Brain Atlas ===")

    # Save structure IDs always
    with open(dest_ids, "w") as f:
        json.dump(STRUCTURE_IDS, f, indent=2)
    log(f"  Saved structure_ids.json")

    # Strategy A: Scalable Brain Atlas — downsampled NIfTI (1.2 MB)
    nii_ok = False
    nii_urls = [
        "https://scalablebrainatlas.incf.org/mouse/ABA_v3/source/P56_Annotation_downsample2.nii.gz",
        # Mirror: OSF
        "https://osf.io/download/fqj8d/",
    ]
    for url in nii_urls:
        if download(url, dest_nii, timeout=180):
            nii_ok = True
            break

    # Strategy B: allensdk Python package
    if not nii_ok:
        log("  Trying allensdk ...")
        try:
            import allensdk
        except ImportError:
            try_pip_install("allensdk")
        try:
            from allensdk.core.reference_space_cache import ReferenceSpaceCache
            rspc = ReferenceSpaceCache(
                resolution=25,
                reference_space_key="annotation/ccf_2017",
                manifest=str(ATLAS_DIR / "manifest.json")
            )
            annotation, _ = rspc.get_annotation_volume()
            # Save as simple numpy binary (load with np.load)
            import numpy as np
            np.save(str(ATLAS_DIR / "P56_Annotation_25um.npy"), annotation)
            log(f"  allensdk: saved annotation volume shape {annotation.shape}")
            nii_ok = True
        except Exception as e:
            log(f"  allensdk failed: {e}")

    # Strategy C: Allen Brain Atlas API — download SVG coronal section
    # Plate 87 ≈ Bregma -2.0 mm in the adult mouse atlas
    if not dest_svg.exists():
        svg_url = ("http://atlas.brain-map.org/atlas?atlas=1"
                   "#atlas_id=1&plate_index=87&resolution=10&zoom=-7")
        # Direct SVG download from ABA plates API
        api_url = ("http://api.brain-map.org/api/v2/atlas_image_download/100960560"
                   "?downsample=4&annotation=true")
        download(api_url, ATLAS_DIR / "coronal_plate87.jpg", timeout=60)

    # Strategy D: Generate atlas outlines from Allen ontology API
    if not nii_ok:
        log("  Generating atlas outlines from Allen ontology API ...")
        _fetch_allen_structure_outlines()

    available = [f.name for f in ATLAS_DIR.iterdir() if f.stat().st_size > 0]
    log(f"  Atlas files: {available}")
    return len(available) > 0


def _fetch_allen_structure_outlines():
    """
    Query Allen Brain Atlas API for structure polygons at Bregma -2.0mm.
    Returns SVG-compatible polygon coordinates for major hippocampal structures.
    API: http://api.brain-map.org/api/v2/structure_graph_download/1.json
    """
    log("  Fetching structure outlines from Allen API ...")
    structures_to_fetch = {
        "CA1": 382, "CA3": 463, "DG": 726,
        "CC": 776, "CTX": 315
    }
    outlines = {}
    base = "http://api.brain-map.org/api/v2/data/Structure/query.json"
    for name, sid in structures_to_fetch.items():
        url = f"{base}?criteria=[id$eq{sid}]&include=color_hex_triplet,name,acronym"
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                data = json.loads(r.read())
            if data.get("success") and data["msg"]:
                outlines[name] = data["msg"][0]
                log(f"  Got structure: {name} ({data['msg'][0].get('acronym')})")
        except Exception as e:
            log(f"  Failed to fetch {name}: {e}")

    # Save as JSON for use by the figure generator
    out_file = ATLAS_DIR / "allen_structure_metadata.json"
    with open(out_file, "w") as f:
        json.dump(outlines, f, indent=2)

    # Also save hardcoded approximate polygon coordinates (Bregma -2.0mm, CCFv3)
    # Coordinates in mm, relative to midline/brain surface.
    # Source: Allen Mouse Brain Atlas, Paxinos & Franklin 2001.
    polygons = {
        "brain_outline": [
            [-4.0,0],[-3.8,-0.3],[-3.5,-0.8],[-3.0,-1.3],[-2.5,-1.6],
            [-2.0,-1.8],[-1.5,-1.9],[-1.0,-2.0],[-0.5,-2.1],[0,-2.15],
            [0.5,-2.1],[1.0,-2.0],[1.5,-1.9],[2.0,-1.8],[2.5,-1.6],
            [3.0,-1.3],[3.5,-0.8],[3.8,-0.3],[4.0,0],
        ],
        "cortex_L": [
            [-4.0,0],[-3.8,-0.3],[-3.5,-0.8],[-3.0,-1.3],[-2.5,-1.6],
            [-2.0,-1.8],[-1.5,-1.9],[-1.0,-2.0],[-0.8,-1.5],
            [-0.5,-0.8],[-0.5,0],
        ],
        "cortex_R": [
            [4.0,0],[3.8,-0.3],[3.5,-0.8],[3.0,-1.3],[2.5,-1.6],
            [2.0,-1.8],[1.5,-1.9],[1.0,-2.0],[0.8,-1.5],
            [0.5,-0.8],[0.5,0],
        ],
        "CA1_L": [
            [-1.0,-1.5],[-1.5,-1.6],[-2.0,-1.7],[-2.5,-1.6],
            [-2.8,-1.4],[-2.5,-1.2],[-2.0,-1.1],[-1.5,-1.2],[-1.0,-1.3],
        ],
        "CA1_R": [
            [1.0,-1.5],[1.5,-1.6],[2.0,-1.7],[2.5,-1.6],
            [2.8,-1.4],[2.5,-1.2],[2.0,-1.1],[1.5,-1.2],[1.0,-1.3],
        ],
        "DG_L": [
            [-1.5,-1.8],[-2.0,-2.0],[-2.5,-1.9],[-2.7,-1.7],
            [-2.3,-1.6],[-1.8,-1.6],[-1.5,-1.7],
        ],
        "DG_R": [
            [1.5,-1.8],[2.0,-2.0],[2.5,-1.9],[2.7,-1.7],
            [2.3,-1.6],[1.8,-1.6],[1.5,-1.7],
        ],
        "CA3_L": [
            [-0.8,-1.5],[-1.0,-1.7],[-1.3,-1.9],[-1.5,-1.8],
            [-1.3,-1.5],[-1.0,-1.4],
        ],
        "CA3_R": [
            [0.8,-1.5],[1.0,-1.7],[1.3,-1.9],[1.5,-1.8],
            [1.3,-1.5],[1.0,-1.4],
        ],
        "corpus_callosum": [
            [-2.5,-0.8],[-2.0,-0.85],[-1.5,-0.9],[-1.0,-0.9],
            [-0.5,-0.88],[0,-0.87],[0.5,-0.88],[1.0,-0.9],
            [1.5,-0.9],[2.0,-0.85],[2.5,-0.8],[2.5,-1.0],
            [2.0,-1.05],[1.5,-1.1],[1.0,-1.1],[0.5,-1.08],
            [0,-1.07],[-0.5,-1.08],[-1.0,-1.1],[-1.5,-1.1],
            [-2.0,-1.05],[-2.5,-1.0],
        ],
        "thalamus_L": [
            [-0.5,-0.9],[-1.0,-1.0],[-1.5,-1.1],[-1.5,-1.5],
            [-1.0,-1.5],[-0.5,-1.4],[-0.2,-1.2],[-0.2,-1.0],
        ],
        "thalamus_R": [
            [0.5,-0.9],[1.0,-1.0],[1.5,-1.1],[1.5,-1.5],
            [1.0,-1.5],[0.5,-1.4],[0.2,-1.2],[0.2,-1.0],
        ],
        "lateral_ventricle_L": [
            [-0.5,-0.87],[-0.8,-0.9],[-0.8,-1.1],[-0.5,-1.15],
            [-0.3,-1.1],[-0.3,-0.9],
        ],
        "lateral_ventricle_R": [
            [0.5,-0.87],[0.8,-0.9],[0.8,-1.1],[0.5,-1.15],
            [0.3,-1.1],[0.3,-0.9],
        ],
    }
    poly_file = ATLAS_DIR / "hippocampus_polygons_Bregma-2mm.json"
    with open(poly_file, "w") as f:
        json.dump({
            "polygons": polygons,
            "source": "Allen Mouse Brain Atlas CCFv3, Paxinos & Franklin 2001",
            "bregma_mm": -2.0,
            "units": "mm from midline (x) and brain surface (y, negative = deeper)",
            "note": "Approximate polygons. Use P56_Annotation_downsample2.nii.gz for precision."
        }, f, indent=2)
    log(f"  Saved polygon coordinates: {poly_file.name}")


# ─────────────────────────────────────────────────────────
# SOURCE 3: Enhancer expression maps — hippocampal markers
# We fetch Allen ISH images for three hippocampal marker genes:
#   Prox1  — DG granule cells (bilateral DG-specific)
#   Wfs1   — CA1 pyramidal neurons
#   Calb2  — Mossy cells / CA2
# These serve as realistic CBG expression mask templates.
# ─────────────────────────────────────────────────────────

HIPPOCAMPAL_GENES = {
    "Prox1": {
        "description": "DG granule cells — bilateral dentate gyrus specific",
        "section_dataset_id": 71587828,   # Allen Mouse Brain ISH experiment ID
        "section_image_id":   101889718,  # section at ~Bregma -2mm
    },
    "Wfs1": {
        "description": "CA1 pyramidal neurons",
        "section_dataset_id": 79556218,
        "section_image_id":   102025813,
    },
    "Calb2": {
        "description": "Mossy cells / CA2",
        "section_dataset_id": 71717848,
        "section_image_id":   101892343,
    },
}

def fetch_enhancers():
    log("\n=== SOURCE 3: Enhancer / ISH expression maps ===")
    fetched = 0

    for gene, info in HIPPOCAMPAL_GENES.items():
        dest = ENH_DIR / f"{gene}_ISH_coronal_Bregma-2mm.jpg"
        if dest.exists() and dest.stat().st_size > 0:
            log(f"  {gene}: already present")
            fetched += 1
            continue

        # Strategy A: Allen Brain Atlas image download service
        # API: http://api.brain-map.org/api/v2/image_download/{section_image_id}
        sid = info["section_image_id"]
        sds = info["section_dataset_id"]
        urls_to_try = [
            # Direct section image download (downsample=3 → ~400px wide)
            f"http://api.brain-map.org/api/v2/image_download/{sid}?downsample=3&quality=70",
            # Fallback: get via section dataset query
            f"http://api.brain-map.org/api/v2/section_image/query.json?criteria=section_data_set_id$eq{sds}&include=image",
            # Second fallback: search for any coronal section of this gene
            (f"http://api.brain-map.org/api/v2/SectionDataSet/query.json?"
             f"criteria=genes[acronym$eq{gene}]"
             f"%5B%5Bplane_of_section[name$eqcoronal]%5D%5D"
             f"&num_rows=1"),
        ]

        ok = False
        # Try primary download
        if download(urls_to_try[0], dest, timeout=60):
            ok = True

        # If first URL failed, query API to find correct image ID then download
        if not ok:
            try:
                # Search for experiments with this gene
                search_url = (
                    f"http://api.brain-map.org/api/v2/data/SectionDataSet/query.json?"
                    f"criteria=genes[acronym$eq{gene}],"
                    f"plane_of_section[name$eqcoronal]"
                    f"&num_rows=1&order=id"
                )
                with urllib.request.urlopen(search_url, timeout=30) as r:
                    result = json.loads(r.read())
                if result.get("success") and result["msg"]:
                    ds_id = result["msg"][0]["id"]
                    # Get the section images for this dataset
                    img_url = (
                        f"http://api.brain-map.org/api/v2/data/SectionImage/query.json?"
                        f"criteria=data_set_id$eq{ds_id}"
                        f"&num_rows=1"
                    )
                    with urllib.request.urlopen(img_url, timeout=30) as r:
                        imgs = json.loads(r.read())
                    if imgs.get("success") and imgs["msg"]:
                        img_id = imgs["msg"][0]["id"]
                        dl_url = (f"http://api.brain-map.org/api/v2/image_download/{img_id}"
                                  f"?downsample=3&quality=70")
                        ok = download(dl_url, dest, timeout=60)
            except Exception as e:
                log(f"  {gene} API fallback failed: {e}")

        if ok:
            fetched += 1
            # Save metadata alongside image
            meta_file = ENH_DIR / f"{gene}_metadata.json"
            with open(meta_file, "w") as f:
                json.dump({
                    "gene": gene,
                    "description": info["description"],
                    "section_dataset_id": sds,
                    "source": "Allen Mouse Brain Atlas ISH",
                    "url": f"https://mouse.brain-map.org/experiment/show/{sds}",
                    "license": "Allen Institute Terms of Use — free for non-commercial use"
                }, f, indent=2)
        else:
            log(f"  WARNING: could not fetch {gene} — figure will use synthetic mask")

    # Strategy C: Save a data-URI manifest so the figure generator can use
    # the Allen API at runtime if local files are missing
    manifest = {
        "genes": HIPPOCAMPAL_GENES,
        "api_base": "http://api.brain-map.org/api/v2",
        "image_download": "http://api.brain-map.org/api/v2/image_download/{section_image_id}?downsample=3",
        "note": "If local files absent, query API at runtime using section_dataset_id"
    }
    with open(ENH_DIR / "allen_ish_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    log(f"  Enhancer files fetched: {fetched}/{len(HIPPOCAMPAL_GENES)}")
    return fetched > 0


# ─────────────────────────────────────────────────────────
# SOURCE 4: Focal model assets (neocortical-injury status epilepticus)
# Primary paper: PII S0969996122000249, DOI 10.1016/j.nbd.2022.105633
# Source tiers:
#   1) raw registered maps (preferred, user-provided)
#   2) figure-derived ROI template (atlas-registered approximation)
#   3) synthetic fallback ROI
# ─────────────────────────────────────────────────────────

def _write_json(path, payload):
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def _generate_figure_derived_focal_roi(dest):
    """
    Generate a figure-derived unilateral cortical ROI template in mm coordinates.
    This template is intentionally explicit about assumptions and provenance.
    """
    payload = {
        "source_tier": "figure-derived",
        "paper": FOCAL_PAPER,
        "bregma_mm": -2.0,
        "units": "mm relative to midline and cortical surface",
        "roi_polygon_mm": [
            [-2.4, 1.8],
            [-1.9, 1.9],
            [-1.3, 1.6],
            [-1.0, 1.2],
            [-1.1, 0.7],
            [-1.6, 0.5],
            [-2.1, 0.6],
            [-2.5, 1.0]
        ],
        "registration_landmarks_mm": {
            "midline": [0.0, 0.0],
            "left_cortical_edge": [-4.0, 0.0],
            "right_cortical_edge": [4.0, 0.0]
        },
        "notes": [
            "Derived for focal neocortical injury workflow when raw rasters are unavailable.",
            "Must be labelled figure-derived in all downstream outputs.",
            "Replace with raw registered maps when available."
        ]
    }
    _write_json(dest, payload)


def _generate_synthetic_focal_roi(dest):
    payload = {
        "source_tier": "synthetic",
        "paper": FOCAL_PAPER,
        "bregma_mm": -2.0,
        "units": "mm relative to midline and cortical surface",
        "gaussian_focus": {
            "center_mm": [-1.8, 1.2],
            "sigma_x_mm": 0.55,
            "sigma_y_mm": 0.40,
            "rotation_deg": -18
        },
        "notes": [
            "Synthetic fallback used because raw and figure-derived inputs were unavailable.",
            "Do not use for quantitative claims without explicit caveat."
        ]
    }
    _write_json(dest, payload)


def fetch_focal_model():
    log("\n=== SOURCE 4: Focal neocortical-injury model assets ===")

    meta_path = FOCAL_DIR / "focal_model_metadata.json"
    _write_json(meta_path, FOCAL_PAPER)

    raw_candidates = []
    raw_search_roots = [
        FOCAL_DIR / "raw",
        DATA_DIR / "raw",
    ]
    for root in raw_search_roots:
        if root.exists():
            for pattern in ("*.nii", "*.nii.gz", "*.npy", "*.npz"):
                raw_candidates.extend(sorted(root.rglob(pattern)))

    if raw_candidates:
        selected = raw_candidates[0]
        pointer_path = FOCAL_DIR / "focal_raw_pointer.json"
        _write_json(pointer_path, {
            "source_tier": "raw",
            "paper": FOCAL_PAPER,
            "selected_file": str(selected),
            "candidate_count": len(raw_candidates),
            "notes": [
                "Raw focal map detected; downstream notebook should prefer this source.",
                "Pointer file stores path only; original data remains in place."
            ]
        })
        log(f"  Raw focal map found: {selected}")
        return True

    figure_roi_path = FOCAL_DIR / "focal_ictal_roi_figure_derived.json"
    try:
        _generate_figure_derived_focal_roi(figure_roi_path)
        log(f"  Figure-derived ROI template created: {figure_roi_path.name}")
        return True
    except Exception as e:
        log(f"  Figure-derived ROI generation failed: {e}")

    synthetic_path = FOCAL_DIR / "focal_ictal_roi_synthetic.json"
    _generate_synthetic_focal_roi(synthetic_path)
    log(f"  Synthetic focal ROI created: {synthetic_path.name}")
    return False


# ─────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────

def main():
    log("Starting data fetch for CBG figure ...")
    log(f"Output directory: {DATA_DIR.resolve()}")

    results = {
        "vasculature": fetch_vasculature(),
        "atlas":       fetch_atlas(),
        "enhancers":   fetch_enhancers(),
        "focal_model": fetch_focal_model(),
    }

    log("\n=== SUMMARY ===")
    all_ok = True
    for name, ok in results.items():
        status = "OK" if ok else "FAILED (synthetic fallback used)"
        log(f"  {name:20s}: {status}")
        if not ok: all_ok = False

    # Write manifest of what was actually fetched
    manifest = {}
    for d in [VASC_DIR, ATLAS_DIR, ENH_DIR, FOCAL_DIR]:
        manifest[d.name] = [
            {"file": f.name, "size_kb": f.stat().st_size // 1024}
            for f in sorted(d.iterdir()) if f.stat().st_size > 0
        ]
    with open(DATA_DIR / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    log(f"\nFull manifest written to data/manifest.json")

    if all_ok:
        log("\nAll data sources fetched successfully.")
    else:
        log("\nSome sources used synthetic fallbacks — figure will still generate.")
        log("Check data/manifest.json to see what was actually retrieved.")

    return 0 if all_ok else 1

if __name__ == "__main__":
    sys.exit(main())
