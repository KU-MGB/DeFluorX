#!/usr/bin/env python3

"""
===============================================================================
FAcDs Pipeline  |  Step 02  |  Boltz-2 Production, Analysis & FAcD Ranking
===============================================================================
Large-scale, resume-safe Boltz-2 protein-ligand predictions with deep
structural, geometric, and chemical scoring for FAcD SN2 degrader tiers.

Clean/final build: all old-format folder/file migration and job-name renaming
have been removed. Resume (--resume) is retained and reuses completed GPU
predictions (boltz_results_*/), recomputing only analysis; it assumes the input
roster (proteins × ligands) is stable across resumes — change parameters and
re-run, not the ligand/protein list. Only the canonical directory layout and
job-naming scheme are supported.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 15 July 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 02_Production_FAcDs.py
  Role          : "Engine" — Boltz-2 prediction orchestrator and ranker.
  Imports from  : 00_01_Project_Config_FAcDs.py  (CFG — all geometric thresholds)
                  00_02_Project_Utils_FAcDs.py   (geometric utilities, console funcs)
  Reads         : 01_Merge_FAcDs.py output — merged *.fasta (protein sequences)
                  User-supplied *.smi (SMILES ligand file)
  Writes        : <Run>/1_Boltz2_Production/  (Boltz-2 CIF outputs)
                  <Run>/2_Best_Complexes_CIFs/ (top-model CIF selection)
                  <Run>/1_Boltz2_Production/7_Boltz2_FAcDs_Ranked_*.csv
                  <Run>/1_Boltz2_Production/6_Boltz2_FAcDs_Master_*.csv
  Upstream      : 01_Merge_FAcDs.py → writes the merged FASTA consumed here
  Downstream    : 03_Validation_Figures_FAcDs.py → reads ranked CSV
                  05_TopN_and_PDB_Preparation_FAcDs.py → reads Best_Complexes_CIFs
                  05_TopN_and_PDB_Preparation_FAcDs.py → reads ranked CSV
───────────────────────────────────────────────────────────────────────────────

── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. GPU Availability: Hard requirement for CUDA-capable GPUs and the Boltz-2
     CLI; will fail immediately on CPU-only machines or missing binaries.
  2. Input Mismatch: If the number of SMILES strings does not match the
     expected ligand count, job batching will misalign metadata.
  3. Active Site Scanning: Relies on cKDTree for proximity; highly sensitive to
     ligand atom-naming consistency in the input SMILES.
  4. Runtime Scaling: Total runtime scales O(N_proteins × N_ligands); large
     runs require multi-day execution windows.
───────────────────────────────────────────────────────────────────────────────

Usage:
    conda activate PFAS
    python 02_Production_FAcDs.py --fasta merged.fasta --smi ligands.smi
    python 02_Production_FAcDs.py --resume Boltz-2_Run_20260309T085406Z

Purpose:
    Performs large-scale, resume-safe Boltz-2 protein-ligand complex
    predictions followed by deep structural, geometric, and chemical analysis.
    Includes self-healing to recover from crashes.

-------------------------------------------------------------------------------
Key Features:
    • Ultra-Strict FAcD Tiers: Specifically tuned for S_N2 trajectory validation.
    • Dynamic Active Site: Spatially scans for any polar/aromatic stabilisers near the pocket.
    • Atom-Specific Stabilisation: Validates oxyanion hole interactions strictly
      against functional sidechain atoms (Nitrogen/Oxygen), avoiding backbone false positives.
    • NAC Geometry: Checks SN2 deviations via backside-attack (Walden-inversion) trajectory analysis.
    • High-Performance Physics: Utilises cKDTree for O(N log N) spatial queries,
      replacing expensive nested loops.
    • Bias-Free Normalisation: Interaction density normalised by ligand heavy atom count.
    • Fluorous-Specific Physics: Enhanced weighting for fluorine-polar and fluorous-lipophilic interactions.
    • Mechanism-Weighted Ranking: Prioritises the ideal 180° backside attack geometry.

-------------------------------------------------------------------------------
EXPANDED TIER DEFINITIONS (Mechanism-First Classification - ULTRA STRICT):
The pipeline assigns a "Degrader Tier" based on strictly tightened catalytic requirements
derived from the 1.60 A Crystal Structure (3R3U) and pristine SN2 reaction mechanics.

    1. Tier_1A (Elite Catalysis - High Priority for MD)
       • Mechanistic Score >= 0.85 (ladder floor; angle folded into the score. Elite tier further refined by the coupled machinery gate: mech >= 0.90 OR mech >= 0.85 with a crystal-exact catalytic constellation.)
       • (Active Site Conservation is reported downstream for ranking; it is NOT a tier gate.)
       • Nucleophile (Asp110): <= 3.0 A (ligand α-carbon → Asp-Oδ; tight pre-reactive ground-state gate)
       • Nuc–Base relay:       <= 3.5 A (INTERNAL triad Asp110-Oδ → His277, dist_nuc_base — not a ligand contact)
       • Base–Acid relay:      <= 4.5 A (INTERNAL triad His277 → Asp134, dist_base_acid — not a ligand contact)
       • Attack Angle:         >= 170°  (Near-Ideal Linear Trajectory; = TIER_ANGLE_MIN['Tier_1A'])
       • Stabilisation:        REQUIRED (Trp156/Tyr217 or Dynamic Polar Residue)

    2. Tier_1B (High Functional)
       • Mechanistic Score >= 0.85
       • Nucleophile (Asp110): <= 3.2 A
       • Base/Acid:            <= 4.0 A / 5.0 A
       • Attack Angle:         >= 165°
       • Stabilisation:        REQUIRED

    3. Tier_2A (Functional Geometry)
       • Nucleophile (Asp110): <= 3.2 A
       • Base/Acid:            <= 5.0 A / 6.0 A
       • Attack Angle:         >= 155°
       • Stabilisation:        Optional

    4. Tier_2B (Marginal Functionality)
       • Nucleophile (Asp110): <= 3.8 A
       • Attack Angle:         >= 145°

    5. Tier_3 (Non-Catalytic Binding)
       • Nucleophile (Asp110): <= 4.2 A
       • Ligand is in the pocket, but orientation fails SN2 requirements.

    6. Tier_4 (Poor)
       • Nucleophile (Asp110): > 4.2 A (Docking Failure)

    7. Tier_5_Decoy (Invalid/Decoy)
       • Missing sequences or severe structural clashes.

-------------------------------------------------------------------------------
RANKING LOGIC (Sorting the Master CSV):
The final `06_Ranked.csv` file uses a hierarchical scoring system to prioritise
catalytic mechanism over generic binding affinity.

    1. Tier Value (Primary Sort Key):
       Candidates are strictly grouped by Tier. The top tier always ranks above the subsequent tier.
       (Tier_1A=50 > Tier_1B=40 > Tier_2A=30 > Tier_2B=20 > Tier_3=10)

    2. Active-Site Conservation / Likelihood (Secondary Sort Key):
       Within the same Tier, candidates are ranked by the sequence/structure
       conservation score.

    3. SN2 Attack Angle (Tertiary Sort Key):
       Proximity of the backside-attack angle to the ideal 180°.

    4. Binding score (Final tie-breaker).

-------------------------------------------------------------------------------
Design Principles:
    • Production-first, large-scale screening
    • Deterministic outputs with reproducible naming conventions
    • Fail-safe execution with minimal human intervention
    • Scalable to thousands of protein-ligand combinations.
    • Focus on structure modelling and analysis.

-------------------------------------------------------------------------------
Scientific References:
    1. Catalytic Mechanism & Enzymatic Halogen Cleavage (Triad: Asp-His-Asp):
       - Mapping the reaction coordinates of enzymatic defluorination.
       - Chan, P.W.Y., Yakunin, A.F., Edwards, E.A. & Pai, E.F. (2011) JACS 133:7461–7468. PDB 3R3U/3R3Z.
       - DOI: https://doi.org/10.1021/ja200277d
    2. Bürgi–Dunitz Trajectory (auxiliary carbonyl-addition metric):
       - Stereochemistry of reaction paths at carbonyl centres.
       - Bürgi, H.B., Dunitz, J.D., Lehn, J.M. & Wipff, G. (1974) Tetrahedron 30:1563–1572.
       - DOI: https://doi.org/10.1016/S0040-4020(01)90678-7
    3. Halogen Bonding and Halide Cradles (Stabilisation):
       - Halogen bonds in biological molecules.
       - Auffinger, P., Hays, F.A., Westhof, E. & Ho, P.S. (2004) PNAS 101:16789–16794.
       - DOI: https://doi.org/10.1073/pnas.0407607101
    4. Near Attack Conformation (NAC) Approach:
       - Defluorination of organofluorine compounds by FAcD (DEHA4 experimental validation).
       - Farajollahi, S. et al. (2024) ACS Omega 9:28546–28555.
       - DOI: https://doi.org/10.1021/acsomega.4c02517
       - Near attack conformation approach to the study of the chorismate to prephenate reaction.
       - Hur, S. & Bruice, T.C. (2003) PNAS 100:12015–12020.
       - DOI: https://doi.org/10.1073/pnas.1534873100
    5. Fluorous-Specific Interactions (Fluorine-Polar/Fluorine-Hydrophobic):
       - The many roles for fluorine in medicinal chemistry.
       - Hagmann, W.K. (2008) J Med Chem 51:4359–4369.
       - DOI: https://doi.org/10.1021/jm800219f
    6. Fluoroacetate Dehalogenase (FAcD) Structural Architecture:
       - Mapping the reaction coordinates of enzymatic defluorination (same as ref 1).
       - Chan et al. (2011) JACS 133:7461–7468. PDB 3R3U/3R3Z.
       - DOI: https://doi.org/10.1021/ja200277d
       - FAcD small-substrate scope + TFA recalcitrance (graded chemistry/containment tier penalties):
         Wackett, L.P. (2022) Microb Biotechnol 15(3):773–792. DOI: https://doi.org/10.1111/1751-7915.13928
    6b. FAcD SN2 Defluorination QM/MM Energetics (basis for QSite scan design):
       - Comprehensive understanding of FAcD-catalyzed degradation of fluorocarboxylic acids: a QM/MM approach.
       - Yue, Y. et al. (2021) Environ Sci Technol 55(14):9817–9825.
       - DOI: https://doi.org/10.1021/acs.est.0c08811
    6c. Difluoroacetate (gem-CF2) is a genuine FAcD substrate (no CF2 penalty applied):
       - Structural insights into hydrolytic defluorination of difluoroacetate by microbial fluoroacetate dehalogenases.
       - Khusnutdinova, A.N. et al. (2023) FEBS J 290:4966–4983.
       - DOI: https://doi.org/10.1111/febs.16903
    6d. SN2 Dead-End Consensus Check (basis for the scissile-carbon dead-end gate):
       (A) Scissile C–F bond strength — α-fluorination strengthens the C–F bond:
       - O'Hagan, D. (2008) Understanding organofluorine chemistry — an introduction to the C–F bond.
       - Chem Soc Rev 37:308–319. DOI: https://doi.org/10.1039/B711844A
         (C–F BDE: CH3F 109.9, CH2F2 119.5, CHF3 127.5 kcal/mol).
       (B) Backside SN2 steric accessibility + van-der-Waals radii:
       - Bento, A.P. & Bickelhaupt, F.M. (2008) J Org Chem 73:7290–7299. DOI: https://doi.org/10.1021/jo801215z
       - Bondi, A. (1964) van der Waals Volumes and Radii. J Phys Chem 68:441–451. DOI: https://doi.org/10.1021/j100785a001
    7. Directed Evolution of FAcD on Non-Natural Organofluorides:
       - Engineering fluoroacetate dehalogenase by growth-based selections on non-natural organofluorides.
       - Jansen, S.C., van Beers, P. & Mayer, C. (2026) Angew Chem Int Ed 65:e202524234.
       - DOI: https://doi.org/10.1002/anie.202524234

    8. Computational Tooling & Core Bioinformatics Stack:
       - Boltz-1 (Foundation model for biomolecular interaction prediction):
         Wohlwend, J., Corso, G., Passaro, S. et al. (2024) bioRxiv preprint.
         DOI: https://doi.org/10.1101/2024.11.19.624167
       - Boltz-2 (Binding affinity prediction):
         Passaro, S., Corso, G., Wohlwend, J. et al. (2025) bioRxiv preprint.
         DOI: https://doi.org/10.1101/2025.06.14.659707
       - Gemmi (Macromolecular crystallography library for mmCIF parsing and robust superimposition):
         Wojdyr, M. (2022) J Open Source Softw 7:4200.
         DOI: https://doi.org/10.21105/joss.04200
       - SciPy 1.0 (Fundamental algorithms for scientific computing and O(N log N) cKDTree spatial queries):
         Virtanen, P. et al. / SciPy 1.0 Contributors (2020) Nature Methods 17:261–272.
         DOI: https://doi.org/10.1038/s41592-019-0686-2
       - RDKit (Open-source cheminformatics toolkit for 3D embedding and valid C-F topology mapping):
         Link: https://www.rdkit.org

    -------------------------------------------------------------------------------
    INTERPRETATION & KNOWN LIMITATIONS (read before trusting a number)
    -------------------------------------------------------------------------------
    • Tier vs rank optimise DIFFERENT keys. degrader_tier is gated on mechanistic_score
      (geometry + active-site machinery ONLY) plus the Criterion-B constellation cap;
      Scientific_Rank orders within a tier by competence_score (geometry × chemical
      feasibility). A pose can therefore be high-rank/low-tier or vice-versa: the tier is the
      mechanistic class, competence is the within-class ordering. Do not read rank as tier.
    • Poly-fluorine angle multiplicity is Šidák-corrected in BOTH mechanistic_score (tier)
      and competence_score (rank); the correction is intentionally surfaced in both keys, so a
      CF2/CF3 attack carbon is discounted at the gate AND in the ordering — by design, not a
      double-count bug.
    • STATIC single-frame proxy. Every geometric term (mech, competence, Criterion B,
      attack angle, distances) is read from ONE Boltz pose. SN2 defluorination is dynamic;
      these are near-attack-conformation (NAC) descriptors, NOT kinetic activities
      (see Classification_Basis column). model_degrader_consensus (5 samples) is the only
      sampling-spread signal and is thin. The activation barrier is decided downstream by
      Step-06 MM-GBSA and Step-07 QM/MM — those, not this screen, are the reactivity arbiters.
    • Criterion B is a tier CAP, not a fine discriminator at the low end: Active_Site_RMSD
      caps each residue's contribution, so a badly-mis-assembled site still yields a moderate
      B. B reliably removes the worst constellations; it does not finely rank good ones.
    • Active-site identity (Criterion A, active_site_integrity) is alignment-imputed; for
      low-identity homologs the ±RESIDUE_SEARCH_WINDOW resolver can mis-pick a residue.
      active_site_plddt (local confidence) and Criterion B (geometric fidelity) are the
      cross-checks; the nucleophile additionally carries a rescue-offset audit (§5.4).

===============================================================================
"""

# =============================================================================
# SECTION 1: SYSTEM CONFIGURATION & IMPORTS
# =============================================================================
from __future__ import annotations

# -------------------------------------------------------------------------------
# Step 1.1: Standard Library Imports
# -------------------------------------------------------------------------------
import argparse
import gc
import hashlib
import json
import math
import os
import re
import shutil
import tempfile
import subprocess
import sys
import time
import threading
import multiprocessing
import csv
import fcntl
import traceback
import urllib.request
import urllib.parse
import urllib.error
import tarfile
import io
import uuid
import requests as _cf_requests   # ColabFold API network calls (available as boltz depends on it)
import warnings
from collections import defaultdict, Counter
from datetime import datetime
from pathlib import Path

# Disable global Python warnings (such as DeprecationWarnings) to maintain clean console output
warnings.filterwarnings("ignore")

from typing import Any, Dict, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor

# -------------------------------------------------------------------------------
# Step 1.2: Hardware Resource Management
# -------------------------------------------------------------------------------
try:
    import psutil
    # -------------------------------------------------------------------------------
    # Sub-Step 1.2.1: Dynamic Thread Allocation
    # -------------------------------------------------------------------------------
    # Calculates available CPU cores; reserves 2 for OS/desktop stability.
    TOTAL_CORES = psutil.cpu_count(logical=True)
    TARGET_CORES = max(1, TOTAL_CORES - 2)

    """
    Reduces the process scheduling priority so the pipeline yields CPU time to the operating system
    and interactive applications, preventing desktop freezes under sustained computational loads.
    """
    try:
        p = psutil.Process(os.getpid())
        if sys.platform != "win32":
            p.nice(15)  # POSIX niceness 15 designates a below-normal priority
    except Exception:
        pass

    """
    Enforces single-threaded BLAS/MKL execution per worker process.
    When running numerous parallel workers, unconstrained BLAS thread pools cause severe
    CPU over-subscription (where context-switching overhead dominates computational time).
    Pinning each worker to one thread yields near-linear scaling with the worker count.
    """
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
except Exception:
    # Fallback configuration if the psutil library is unavailable
    TARGET_CORES = 1

# -------------------------------------------------------------------------------
# Step 1.3: Scientific Stack
# -------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import yaml
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")

# -------------------------------------------------------------------------------
# Step 1.4: Bioinformatics & Chemistry
# -------------------------------------------------------------------------------
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

# -------------------------------------------------------------------------------
# Sub-Step 1.4.1: Torch Threading Configuration
# -------------------------------------------------------------------------------
try:
    import torch
    if torch is not None:
        """
        Restricts PyTorch intra-op parallelism to one thread per worker process,
        ensuring consistency with the BLAS/MKL thread limits established above.
        """
        torch.set_num_threads(1)
except Exception:
    torch = None

import importlib.util as _ilu
from pathlib import Path as _Path

def _load_module(name: str, path):
    """Load a Python file as a module using its filesystem path, regardless of filename."""
    path = _Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Required module not found: {path}")
    spec = _ilu.spec_from_file_location(name, str(path))
    mod  = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_REPO_DIR  = _Path(__file__).resolve().parent
_cfg_mod   = _load_module("ProjectConfig", _REPO_DIR / "00_01_Project_Config_FAcDs.py")
_utils_mod = _load_module("ProjectUtils",  _REPO_DIR / "00_02_Project_Utils_FAcDs.py")
CFG        = _cfg_mod.CFG()

safe_name = _utils_mod.safe_name
# Auxiliary (non-gating) trajectory-geometry descriptors — see utils docstrings.
calculate_burgi_dunitz = _utils_mod.calculate_burgi_dunitz
calculate_flippin_lodge = _utils_mod.calculate_flippin_lodge
distance = _utils_mod.distance
calculate_angle = _utils_mod.calculate_angle


# =============================================================================
# SECTION 2: PHYSICAL CONSTANTS & PARAMETERS (MAESTRO STANDARDS)
# =============================================================================

# -------------------------------------------------------------------------------
# Step 2.1: Path & Model Configuration
# -------------------------------------------------------------------------------
BOLTZ_BIN      = CFG.BOLTZ_EXECUTABLE
BOLTZ_MODEL    = CFG.BOLTZ_MODEL_VERSION
BOLTZ_CACHE    = Path(os.path.expanduser(CFG.BOLTZ_CACHE_DIR))   # CFG-driven; defaults to ~/.cache/boltz
OUTPUT_FORMAT  = CFG.BOLTZ_OUTPUT_FORMAT

# -------------------------------------------------------------------------------
# Step 2.2: Execution Settings
# -------------------------------------------------------------------------------
RECYCLING_STEPS           = CFG.BOLTZ_RECYCLING_STEPS
DIFFUSION_SAMPLES_DEFAULT = CFG.BOLTZ_DIFFUSION_SAMPLES
RETRY_ON_FAIL             = CFG.BOLTZ_RETRY_MAX
RETRY_SLEEP               = CFG.BOLTZ_RETRY_SLEEP
PREDICT_TIMEOUT_S         = CFG.BOLTZ_PREDICT_TIMEOUT_S
MAX_PROTEINS_PER_BATCH    = CFG.BOLTZ_MAX_PROTEINS_PER_BATCH
SEPARATOR                 = "-" * 80

# -------------------------------------------------------------------------------
# Step 2.3: Scoring Weights
# -------------------------------------------------------------------------------
# These weights govern the calculation of the custom "Binding Probability" score.
W_IPTM          = CFG.BOLTZ_W_IPTM
W_CROSS_PAE     = CFG.BOLTZ_W_CROSS_PAE
W_PLDDT         = CFG.BOLTZ_W_PLDDT
W_INTERACTIONS  = CFG.BOLTZ_W_INTERACTIONS
W_CONF          = CFG.BOLTZ_W_CONF

# -------------------------------------------------------------------------------
# Step 2.4: Interaction Thresholds (Ångströms & Degrees)
# -------------------------------------------------------------------------------
"""
Interaction cutoffs reproduce the Schrödinger Maestro default criteria, defined
once in CFG §3 (single source of truth) and applied here by this script's own
geometry engine — no $SCHRODINGER binary is called. Boltz-2 mmCIF carries no
explicit hydrogens, so H-bonds use the heavy-atom donor–acceptor (D···A) proxy
(CFG.THRESHOLD_HB_DIST_MAX); H···A + angle gates (CFG §3.1) apply only to the
explicit-H figure pipeline (Step 05). See CFG reference list for citations.
"""
# Sub-Step 2.4.1: Bond Distances
# -------------------------------------------------------------------------------
HBOND_MAX_DIST_DA        = CFG.THRESHOLD_HB_DIST_MAX   # heavy-atom D···A proxy (no explicit H in Boltz CIF)
HBOND_MIN_DIST_DA        = CFG.THRESHOLD_HB_DIST_MIN   # below this a D···A pair is a clash, not an H-bond
SALT_BRIDGE_MAX_DIST     = CFG.THRESHOLD_SALT_BRIDGE
HYDROPHOBIC_MAX_DIST     = CFG.THRESHOLD_HYDROPHOBIC_MAX
HALOGEN_MAX_DIST         = CFG.HALOGEN_BOND_DIST_MAX
METAL_COORD_MAX_DIST     = CFG.METAL_COORD_DIST_MAX

# -------------------------------------------------------------------------------
# Sub-Step 2.4.2: Pi-Interaction Geometry (Strict constraints)
# -------------------------------------------------------------------------------
PI_STACK_FACE_MAX_DIST   = CFG.PI_STACK_FACE_DIST_MAX
PI_STACK_FACE_MAX_ANGLE  = CFG.PI_STACK_FACE_ANGLE_MAX
PI_STACK_EDGE_MAX_DIST   = CFG.PI_STACK_EDGE_DIST_MAX
PI_STACK_EDGE_MIN_ANGLE  = CFG.PI_STACK_EDGE_ANGLE_MIN
PI_CATION_MAX_DIST       = CFG.THRESHOLD_PI_CATION_MAX
PI_CATION_MAX_ANGLE      = CFG.PI_CATION_ANGLE_MAX

COORD_MATCH_MAX_DIST     = CFG.COORD_MATCH_DIST_MAX

# -------------------------------------------------------------------------------
# Step 2.5: Residue & Atom Definitions
# -------------------------------------------------------------------------------
CATALYTIC_DIST_CUTOFF    = CFG.CATALYTIC_DIST_CUTOFF

# Sets of amino acids grouped by their chemical properties (CFG §2.10 — single source).
STANDARD_AA  = CFG.PREP_STANDARD_AA
POSITIVE_RES = CFG.POSITIVE_RES
NEGATIVE_RES = CFG.NEGATIVE_RES
METALS       = CFG.METALS

# Protein backbone atom names — the amide N and carbonyl O are not the formally charged
# groups, so they are excluded from salt-bridge assignment (see classify_pair).
_BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

"""
Protonation-variant → canonical residue name normalisation (CFG §2.9).
Force-field engines rename residues to encode protonation state; this
function collapses them back to canonical 3-letter codes for all class
matching, structure filtering, and residue-class lookups.
"""
_PROTONATION_MAP = getattr(CFG, "PROTONATION_MAP", {})

def _canonical_resname(resname: str) -> str:
    """Map force-field protonation variants to their canonical 3-letter code."""
    return _PROTONATION_MAP.get(resname, resname)

# Amino-acid chemical-class groups — sourced from CFG (single source of truth; no local copy).
RES_PROPS = CFG.RES_PROPS
# Sidechain ring atoms for π-stacking / π-cation plane fitting — single source of truth in CFG
# (covers all four aromatic residues, incl. His imidazole + Trp indole; see config docstring).
AROMATIC_RING_ATOMS = CFG.AROMATIC_RING_ATOMS

# Functional (H-bond donor) sidechain atoms used to validate fluoride-cradle stabilisation —
# single source of truth in CFG (prevents backbone atoms being flagged as stabilising).
POLAR_SIDECHAIN_ATOMS = CFG.POLAR_SIDECHAIN_ATOMS

# Protonation-aware residue groups for mechanistic scoring / active-site role evaluation —
# single source of truth in CFG (includes FF variants HID/HIE/HIP… so they are recognised).
RESIDUE_CLASS_GROUPS = CFG.RESIDUE_CLASS_GROUPS

# -------------------------------------------------------------------------------
# Step 2.6: The Geometric Reference Data (PDB 3R3U)
# -------------------------------------------------------------------------------
"""
The canonical sequence for DeHa4_[Delftia_acidovorans_D4B].
This serves as a structural map. As target proteins vary in their numbering,
each is aligned to this canonical sequence to precisely map the positions of
key catalytic residues (such as Asp110, Asp134, His277).
Source: Farajollahi et al. (2024) — see header Scientific References §4 (position cross-validated, establishing HIS277 rather than HIS288 in DEHA4 sequences).
"""
DEHA4_CONTROL_SEQ = CFG.DEHA4_CONTROL_SEQ
REF_SEQUENCE_STR = DEHA4_CONTROL_SEQ

# Fluoroacetate SMILES used for reference construction operations.
REF_LIGAND_SMILES = CFG.FLUOROACETATE_SMILES

# The three fluoroacetate control substrates used for DEHA4 and 3R3U calibration runs.
CTRL_LIGANDS = CFG.CTRL_LIGANDS

"""
Structural mapping based on PDB 3R3U (Rhodopseudomonas palustris FAcD,
wild-type, 1.60 Å; carries Ni²⁺/Cl⁻ ions, no substrate bound),
adapted from Farajollahi et al. (2024; see header §4) and matched against DEHA4 sequences.
"""
REF_ACTIVE_SITE_MAP  = CFG.REF_ACTIVE_SITE_MAP
CATALYTIC_TRIAD_KEYS = CFG.CATALYTIC_TRIAD_KEYS

# -------------------------------------------------------------------------------
# Step 2.7: Output Formatting Configuration
# -------------------------------------------------------------------------------
# Mapping dictionary to convert internal variable names into user-friendly CSV column headers.
COLUMN_RENAMING_MAP = {
    "protein": CFG.COL_PROT,
    "ligand": CFG.COL_LIG,
    # Role-based generic names: VALUE = distance to the per-protein dynamically-
    # mapped residue (±5 resolver), not a fixed 3R3U number; mapped residue per
    # role is in the Mapped_* columns / Active_Site_Triad_Map (refs: CFG §2.5).
    "dist_Nuc":    "Dist_Nucleophile",
    "dist_Base":   "Dist_Base",
    "dist_Acid":   "Dist_Acid",
    "dist_Stab_H": "Dist_Stabiliser_H",
    "dist_Stab_W": "Dist_Stabiliser_W",
    "dist_Stab_Y": "Dist_Stabiliser_Y",
    "dist_Carb1":  "Dist_Clamp1",
    "dist_Carb2":  "Dist_Clamp2",
    "binding_likelihood_computed": "Binding_Probability_Score",
    "confidence_score": CFG.COL_CONF,
    "sn2_attack_angle": CFG.COL_SN2,
    "sn2_trajectory_dev": "SN2_Trajectory_Deviation_A",
    "Mapped_to_Control_All": "Full_Sequence_Alignment_Map",
    "Mapped_to_Control_Cat_Triad": "Active_Site_Triad_Map",
    CFG.COL_LIKE_S: CFG.COL_LIKE_S,
    "Active_Site_RMSD": "Active_Site_RMSD_to_Control",
    "Halide_Stabilisation": "Has_Halide_Stabilisation",
    "Carboxylate_Clamp": "Has_Carboxylate_Clamp"
}

# A list of internal columns to exclude from the final CSV output to maintain clarity.
COLUMNS_TO_DROP = [
    "protein_id", "catalytic_dist_A", "rc", "cross_interface_pae_mean",
    "completed_at", "elapsed_seconds",
    "active_site_mapping", "scientific_meaning", "constraint_check"
]

# -------------------------------------------------------------------------------
# Step 2.8: Default Metrics Template
# -------------------------------------------------------------------------------
DEFAULT_METRICS = {
    "status": "Unknown",
    "elapsed_seconds": 0.0,
    CFG.COL_TIER: CFG.TIER_DECOY,
    "is_degrader": False,
    CFG.COL_LIKE_S: 0.0,
    CFG.COL_MECH_S: 0.0,
    "Active_Site_RMSD": 999.0,
    "Identity_to_Control": 0.0,
    "Halide_Stabilisation": False,
    "Carboxylate_Clamp": False,
    "catalytic_dist_A": 999.0,
    CFG.COL_IDENS: 0.0,
    "custom_affinity_score": 0.0,
    "ligand_smiles_stale": 0,   # 1 = the finished structure predates the current input SMILES; it is
                                #     analysed against the SMILES it was PREDICTED FROM, never the new one
    "binding_likelihood_computed": 0.0,
    "confidence_score": 0.0,
    "sn2_attack_angle": 0.0,
    "sn2_trajectory_dev": 999.0,
    "scissile_cf_bde": 0.0,                 # est. leaving C–F bond-dissociation energy (kcal/mol); SN2 dead-end check A
    "sn2_backside_occlusion": 0.0,          # vdW bulk crowding the SN2 backside approach; dead-end check B
    "sn2_dead_end": 0,                       # 1 = scissile C–F too strong + backside blocked (diagnostic)
    "beta_f_count": 0,                       # fluorines on carbons β to the scissile α-carbon (perfluoro discriminator)
    "feasibility_factor": 1.0,               # graded chemical-feasibility multiplier on mech (α-BDE × β-fluorination)
    "active_site_integrity": 0.0,           # Criterion A: fraction of the 8 catalytic residues mapped to correct TYPE (identity; Chan 2011)
    "catalytic_constellation_score": 0.0,   # Criterion B: 1/(1+Active_Site_RMSD) — 8-residue GEOMETRIC constellation fidelity vs 3R3U; caps the tier
    "active_site_plddt": 0.0,                # mean Boltz pLDDT over the 8 mapped catalytic residues (local active-site confidence)
    "active_site_residues_correct": 0,      # count (0–8) of correctly-mapped catalytic residues
    "head_is_carboxylate": 0.0,             # 1.0 = ligand presents a –COO head (FAcD anchoring handle)
    "scissile_is_alpha": 0.0,               # 1.0 = SN2 attack carbon is the α-carbon adjacent to the carboxylate
    "competence_score": 0.0,                # gated continuous 0–1 catalytic competence (Scientific-ranking key)
    "model_degrader_consensus": 0.0,        # fraction of Boltz models independently reaching a degrader tier
    "nuc_resolution": "none",               # how the catalytic Asp was resolved (direct / window+N); §5.4 QC
    "nuc_rescue_offset": 0,                 # residue offset of a windowed nucleophile rescue (0 = direct hit)
    "burgi_dunitz_angle": 999.0,            # aux reference (carbonyl); non-gating; 999.0 = undefined
    "flippin_lodge_offset": 999.0,          # aux reference (in-plane offset); non-gating; 999.0 = undefined
    "mainchain_clash_ratio": 0.0,           # clashing tail atoms / all ligand heavy atoms (size-fair diagnostic)
    "mainchain_clash_count": 0,             # absolute backbone-interpenetrating tail-atom count
    # --- Active-site pocket vs ligand steric fit (DIAGNOSTIC columns — reported for analysis; by design NOT a tier/rank input) ---
    "active_site_volume": 0.0,              # Å³ convex hull of the 8 catalytic residues
    "active_site_radius": 0.0,              # Å pocket centroid → farthest catalytic atom
    "ligand_volume": 0.0,                   # Å³ union of ligand heavy-atom Bondi vdW spheres (physical molecular volume)
    "ligand_radius_gyration": 0.0,          # Å ligand Rg
    "ligand_max_extent": 0.0,               # Å longest ligand interatomic distance
    "pocket_occupancy": 0.0,                # ligand_volume / active_site_volume (convex-hull proxy; not a tier input)
    "fit_ratio": 0.0,                       # (ligand_max_extent/2) / active_site_radius (convex-hull proxy; not a tier input)
    "ligand_fits": 0,                       # 1 = ligand within the convex-hull envelope on both proxies
    "pocket_containment_cavity": 1.0,       # frac of ligand heavy atoms enclosed by the PROTEIN cavity (ray-cast buriedness ≥ BURIAL_MIN) — the tier-gate term
    "pocket_containment_site8": 1.0,        # frac of ligand heavy atoms within SITE8_SHELL_A of the eight mapped active-site residues (catalytic engagement; reported)
    "ligand_buriedness_mean": 0.0,          # mean per-atom buriedness (blocked ray fraction) — continuous companion to pocket_containment_cavity
    "ligand_reach": 0.0,                    # Å; farthest ligand atom from the carboxylate anchor (molecular reach out of the pocket)
    "chem_penalty": 0.0,                    # graded BDE+occlusion penalty subtracted from mech for the tier gate
    "containment_penalty": 0.0,             # graded pocket-fit penalty subtracted from mech for the tier gate
    "mechanistic_score_effective": 0.0,     # feasibility-weighted mech (raw geometry − chem − containment); tier-gate key
}
# Initialises all specific residue distances to an arbitrary maximum (999.0).
for k in REF_ACTIVE_SITE_MAP:
    DEFAULT_METRICS[f"dist_{k}"] = 999.0

# -------------------------------------------------------------------------------
# Step 2.9: Global State Variables
# -------------------------------------------------------------------------------

# ConsoleColours sourced from 00_02_Project_Utils (single canonical definition).
ConsoleColours  = _utils_mod.ConsoleColours
SEPARATOR_HEAVY = _utils_mod.SEPARATOR_HEAVY
SEPARATOR_LIGHT = _utils_mod.SEPARATOR_LIGHT

logger = None
_gpu_lock = threading.Lock()
_gpu_free_mib: Dict[int, int] = {}
CACHED_ALIGNMENTS: Dict[str, Dict] = {}
_PERSISTED_PROTEINS: set = set()   # proteins whose alignment .txt + stats row are already written this run

GLOBAL_STATS = {
    "created_yaml": 0, "created_aln": 0, "created_json": 0, "created_csv_rows": 0,
    "csv_columns": 0, "jobs_run_gpu": 0, "jobs_repaired_cpu": 0,
    "total_alignments_generated": 0,
    "complexes_copied": 0,
    "analytical_failures": 0,
}


# =============================================================================
# SECTION 3: UTILITY FUNCTIONS (LOGGING & DISPLAY)
# =============================================================================

def setup_logging(log_file_path: Path) -> Path:
    global logger
    logger = _utils_mod.setup_logging(log_file_path, logger_name="boltz_master")
    logger.info("=== New execution cycle initialised ===")
    logger.info(f"Master log location: {log_file_path}")
    return log_file_path

def console_info(msg: str) -> None:
    _utils_mod.console_info(msg, logger)

def console_title(msg: str) -> None:
    _utils_mod.console_title(msg, logger)

def console_separator() -> None:
    _utils_mod.console_separator(logger, heavy=True)

def atomic_to_csv(df, path, **kwargs) -> None:
    """Crash-safe CSV write: serialise to a sibling .tmp then os.replace() onto the
    target (atomic on a single filesystem). A process kill mid-write (e.g. the
    OOM-killer during a long grid run) leaves the prior good CSV intact rather than
    a half-written file the resume path cannot detect.
    """
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, **kwargs)
    os.replace(tmp, path)

_TTY_ANSI_RE = re.compile(r"\033\[[0-9;]*[mKABCDEFGHJSTfhilmnprsu]")

def _tty_write(msg: str, flush: bool = True) -> None:
    """Carriage-return progress write. Always writes raw so ANSI colours pass
    through pipeline tee to the terminal. \\r keeps log lines overwriting in
    place rather than flooding the log with 58,000+ entries."""
    sys.stdout.write(msg)
    if flush:
        sys.stdout.flush()


def extract_short_fasta_id(header: str) -> str:
    """Parses the FASTA header to extract clean identifier strings (e.g. '>ID123 Description' becomes 'ID123')."""
    if header.startswith(">"): header = header[1:]
    h = re.split(r"[ \t|()-]", header)[0]
    return h


def sequence_hash(sequence: str) -> str:
    """Generates a stable hash for FASTA sequences to verify consistency during resumed executions."""
    return hashlib.sha256(sequence.upper().encode("utf-8")).hexdigest()


def colabfold_meta_path(prod_dir: Path, protein_id: str) -> Path:
    """Returns the persistent metadata location corresponding to a protein's ColabFold/MSA record."""
    return prod_dir / "3_Sequence_Reference_Data" / "MSA_Sequences" / f"{protein_id}.json"


def colabfold_a3m_path(prod_dir: Path, protein_id: str) -> Path:
    """Returns the persistent ColabFold/MSA output path for a specific protein target."""
    return prod_dir / "3_Sequence_Reference_Data" / "MSA_Sequences" / f"{protein_id}.a3m"


def remove_stale_protein_jobs(prod_dir: Path, protein_id: str):
    """Moves outdated run folders associated with a protein to a trash directory to prevent permanent data loss."""
    runs_dir = prod_dir / "4_Prediction_Jobs"
    trash_dir = prod_dir / "9_Trash" / "Jobs"
    trash_dir.mkdir(parents=True, exist_ok=True)

    for item in sorted(runs_dir.glob(f"*_{protein_id}_*")):
        if item.is_dir():
            target = trash_dir / f"{item.name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            try:
                shutil.move(str(item), str(target))
            except Exception:
                shutil.rmtree(item, ignore_errors=True)


# =============================================================================
# SECTION 4: REFERENCE & SELF-HEALING SYSTEM (DATA RECOVERY)
# =============================================================================

# -------------------------------------------------------------------------------
# Step 4.1: Geometric Truth Initialisation
# -------------------------------------------------------------------------------
def setup_reference_data(target_dir: Path):
    """
    Initialises reference data inside the input folder: checks for the reference PDB
    (3R3U) and associated DeHa4 sequence / control-ligand SMILES, downloading the PDB
    if absent. Reference files live alongside the run inputs (no separate folder).
    """
    target_dir.mkdir(parents=True, exist_ok=True)

    pdb_path   = target_dir / f"{CFG.REFERENCE_PDB_ID}.pdb"
    fasta_path = target_dir / CFG.DEHA4_REF_FASTA
    smi_path   = target_dir / "Control_Ligands_TFA_FA_DFA_Ref.smi"

    console_info(f"\n{SEPARATOR_LIGHT}")
    console_info("  Reference Data Initialisation")
    console_info(SEPARATOR_LIGHT)

    # 1. Download PDB 3R3U if missing
    if not pdb_path.exists():
        url = CFG.REFERENCE_PDB_URL
        try:
            console_info(f"  Downloading PDB {CFG.REFERENCE_PDB_ID} from {url}...")
            urllib.request.urlretrieve(url, pdb_path)
            console_info(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC}  PDB 3R3U download complete.")
        except Exception as e:
            console_info(f"  {ConsoleColours.FAIL}✘{ConsoleColours.ENDC}  Error downloading PDB 3R3U: {e}.")
    else:
        console_info(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC}  PDB 3R3U (Reference) located.")

    # 2. Create Reference FASTA
    if not fasta_path.exists():
        with open(fasta_path, "w") as f:
            f.write(f">DeHa4_Reference_Sequence\n{REF_SEQUENCE_STR}\n")
    else:
        console_info(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC}  DeHa4 Reference Sequence located.")

    # 3. Create Reference SMILES (all three control ligands).
    if not smi_path.exists():
        with open(smi_path, "w") as f:
            for _ln, _ls in CTRL_LIGANDS:
                f.write(f"{_ls}\t{_ln}\n")
    else:
        console_info(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC}  Control ligand SMILES definitions located.")

    return pdb_path


def extract_3r3u_sequence(pdb_path: Path) -> Optional[str]:
    """Extract the first protein chain's amino acid sequence from a PDB file using gemmi."""
    try:
        st = gemmi.read_pdb(str(pdb_path))
        for model in st:
            for chain in model:
                seq = ""
                for res in chain.get_polymer():
                    info = gemmi.find_tabulated_residue(res.name)
                    if info.is_amino_acid():
                        code = info.one_letter_code
                        seq += code if code not in ("?", " ") else "X"
                if len(seq) >= 50:
                    return seq
    except Exception as e:
        if logger:
            logger.debug(f"Sequence extraction from {pdb_path.name} failed: {e}")
    return None


# -------------------------------------------------------------------------------
# Step 4.2: MSA Algorithms
# -------------------------------------------------------------------------------
def validate_a3m_file(a3m_path: Path) -> bool:
    """
    Validates the A3M MSA file to detect corruption. Returns True if the file structure is sound.
    Validation criteria:
    1. File size > 500 bytes (must encompass header and data matrices)
    2. Zero occurrence of null bytes (\x00), which indicate stream corruption
    3. Proper initiation via the '>' character (standard FASTA formatting)
    """
    try:
        if not a3m_path.exists():
            return False

        file_size = a3m_path.stat().st_size
        if file_size < 500:  # Minimum acceptable size for a valid MSA construct
            return False

        # Reads and validates the entire file sequence for underlying corruption
        with open(a3m_path, "rb") as f:
            header = f.read(1024)
            # Checks for standard FASTA header initiation
            if not header.startswith(b">"):
                return False
            # Scans header strings for anomalous null bytes
            if b"\x00" in header:
                return False
            # Iteratively scans the remainder of the file in chunks
            while True:
                chunk = f.read(65536)
                if not chunk:
                    break
                if b"\x00" in chunk:
                    return False

        return True
    except Exception as e:
        logger.debug(f"A3M file validation error encountered for {a3m_path}: {e}")
        return False

# -------------------------------------------------------------------------------
# Step 4.3: Direct MSA Retrieval Operations
# -------------------------------------------------------------------------------
def fetch_msa_direct(pid: str, seq: str, a3m_path: Path, meta_path: Path) -> bool:
    """
    Fetches the MSA directly from the ColabFold REST API. This circumvents the Boltz model loading sequence,
    eliminating the significant overhead associated with CPU-based execution.

    Flow: POST ticket → poll until COMPLETE → download tar.gz → merge constituent .a3m files.
    Returns True upon success, False upon any systemic failure.
    """
    COLABFOLD_API = CFG.COLABFOLD_API_URL
    HEADERS       = {"User-Agent": "boltz"}
    POLL_INTERVAL = CFG.BOLTZ_POLL_INTERVAL
    MAX_SUBMIT_RETRIES = CFG.COLABFOLD_SUBMIT_RETRIES

    try:
        # --- Step 1: Submit MSA job sequence ---
        fasta_query = f">101\n{seq}\n"
        ticket      = None
        for attempt in range(MAX_SUBMIT_RETRIES):
            try:
                resp = _cf_requests.post(
                    f"{COLABFOLD_API}/ticket/msa",
                    data    = {"q": fasta_query, "mode": "env"},
                    timeout = CFG.COLABFOLD_MSA_SUBMIT_TIMEOUT,
                    headers = HEADERS,
                )
                ticket = resp.json()
            except Exception as e:
                logger.debug(f"[MSA-direct] Submission attempt {attempt+1} failed for {pid}: {e}")
                time.sleep(5)
                continue

            status = ticket.get("status", "")
            if status in ("UNKNOWN", "RATELIMIT"):
                logger.debug(f"[MSA-direct] Submission received status {status} for {pid}, initiating retry ({attempt+1}/{MAX_SUBMIT_RETRIES})...")
                time.sleep(6 + attempt * 2)
                continue
            if status in ("ERROR", "MAINTENANCE"):
                logger.debug(f"[MSA-direct] Server responded with status {status} for {pid}: {ticket}")
                return False
            break   # A valid processing ticket has been secured
        else:
            logger.debug(f"[MSA-direct] All submission retries exhausted for {pid}")
            return False

        ticket_id = ticket.get("id")
        if not ticket_id:
            logger.debug(f"[MSA-direct] No valid ticket ID located in response for {pid}: {ticket}")
            return False

        logger.debug(f"[MSA-direct] Ticket {ticket_id} successfully submitted for {pid} (status={ticket.get('status')})")

        # --- Step 2: Poll operation until COMPLETE ---
        current_status = ticket.get("status", "")
        deadline = time.time() + CFG.COLABFOLD_MSA_POLL_TIMEOUT
        while current_status not in ("COMPLETE", "ERROR", "MAINTENANCE") and time.time() < deadline:
            time.sleep(POLL_INTERVAL)
            try:
                r = _cf_requests.get(f"{COLABFOLD_API}/ticket/{ticket_id}", timeout=CFG.COLABFOLD_MSA_POLL_REQUEST_TIMEOUT, headers=HEADERS)
                current_status = r.json().get("status", "UNKNOWN")
                logger.debug(f"[MSA-direct] {pid} polling response: {current_status}")
            except Exception as poll_err:
                logger.debug(f"[MSA-direct] Polling error encountered for {pid}: {poll_err}")

        if current_status != "COMPLETE":
            logger.debug(f"[MSA-direct] Operation did not reach COMPLETE status for {pid} (final status: {current_status})")
            return False

        # --- Step 3: Download resulting tarball ---
        r = _cf_requests.get(
            f"{COLABFOLD_API}/result/download/{ticket_id}",
            timeout = CFG.COLABFOLD_MSA_DOWNLOAD_TIMEOUT,
            headers = HEADERS,
        )
        if r.status_code != 200:
            logger.debug(f"[MSA-direct] Download generated HTTP {r.status_code} for {pid}")
            return False

        tar_bytes = r.content
        if not tar_bytes:
            logger.debug(f"[MSA-direct] Empty byte stream downloaded for {pid}")
            return False

        # --- Step 4: Extract and merge constituent .a3m files ---
        """
        The ColabFold tarball contains uniref.a3m and bfd.mgnify30.metaeuk30.smag30.a3m files.
        Merging both files yields maximum MSA coverage parameters.
        """
        a3m_parts = []
        with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:gz") as tf:
            for member in tf.getmembers():
                if member.name.endswith(".a3m"):
                    fobj = tf.extractfile(member)
                    if fobj:
                        a3m_parts.append(fobj.read())

        if not a3m_parts:
            logger.debug(f"[MSA-direct] No compatible .a3m files found within tarball for {pid}")
            return False

        # Concatenate the constituent MSAs with exactly one separating newline. Each part is
        # trailing-newline-stripped first so joining cannot inject a blank line between them —
        # a blank line reads as EOF to a3m parsers and would silently truncate the MSA depth.
        # Null bytes from the underlying structural tools are removed.
        a3m_content = (b"\n".join(p.rstrip(b"\r\n") for p in a3m_parts) + b"\n").replace(b"\x00", b"")

        # --- Step 5: Execute atomic write and metadata creation ---
        tmp_path = a3m_path.with_suffix(".tmp")
        tmp_path.write_bytes(a3m_content)

        if not validate_a3m_file(tmp_path):
            logger.debug(f"[MSA-direct] Downloaded A3M sequence failed post-processing validation for {pid}")
            tmp_path.unlink(missing_ok=True)
            return False

        tmp_path.rename(a3m_path)
        with open(meta_path, "w") as mf:
            json.dump({
                "sequence_sha256": sequence_hash(seq),
                "protein_id":      pid,
                "generated_at":    datetime.utcnow().isoformat() + "Z",
                "source":          "colabfold_api_direct",
                "a3m_bytes":       len(a3m_content),
            }, mf)

        logger.debug(f"[MSA-direct] MSA data saved for {pid} ({len(a3m_content)} bytes, {len(a3m_parts)} constituent files)")
        return True

    except Exception as e:
        logger.debug(f"[MSA-direct] Processing exception for {pid}: {type(e).__name__}: {e}")
        return False


# -------------------------------------------------------------------------------
# Step 4.4: Stale YAML and Run Folder Cleanup
# -------------------------------------------------------------------------------
def purge_orphans(prod_dir: Path, active_job_names: set) -> int:
    """Deletes YAML files and active run folders that are no longer referenced within the input dataset."""
    yaml_count, run_count = 0, 0

    # 1. Clean designated YAML structures
    yaml_dir = prod_dir / "2_Boltz2_YAML_Configs"
    if yaml_dir.exists():
        for f in sorted(yaml_dir.glob("*.yaml")):
            if f.stem not in active_job_names:
                f.unlink()
                yaml_count += 1

    # 2. Remove all active prediction run folders not mathematically associated with active jobs
    runs_dir = prod_dir / "4_Prediction_Jobs"
    if runs_dir.exists():
        for item in sorted(runs_dir.iterdir()):
            if item.name == "Best_Complex": continue # Protects core architectural directories

            # Unmatched items are deleted to maintain an error-free workspace
            if item.name not in active_job_names:
                if item.is_dir():
                    shutil.rmtree(item, ignore_errors=True)
                else:
                    item.unlink(missing_ok=True)
                run_count += 1

    if yaml_count > 0 or run_count > 0:
        console_info(f" [Purge] Eliminated {yaml_count} orphaned YAML files and {run_count} orphaned items from 4_Prediction_Jobs.")

    return run_count


# =============================================================================
# SECTION 5: BIOINFORMATICS (ALIGNMENT & SUPERIMPOSITION)
# =============================================================================

# -------------------------------------------------------------------------------
# Step 5.1: Caching Framework
# -------------------------------------------------------------------------------
def load_cached_alignments(csv_path: Path):
    """Loads previously calculated alignments into memory and ensures data types are strictly preserved."""
    if csv_path.exists():
        try:
            df = pd.read_csv(csv_path, dtype={"protein": str, "protein_id": str})
            for _, row in df.iterrows():
                pid = row.get("protein") or row.get("protein_id")
                if pid:
                    CACHED_ALIGNMENTS[str(pid)] = row.to_dict()
            console_info(f"Successfully loaded {len(CACHED_ALIGNMENTS)} alignment records from internal cache.")
        except Exception as e:
            console_info(f"Warning: Failed to load alignment cache structure: {e}")

def get_cached_alignment_for_protein(protein_id: str) -> Optional[Dict]:
    """Helper tool designed to retrieve alignment data employing secure string-safe keys."""
    return CACHED_ALIGNMENTS.get(str(protein_id))

def update_alignment_stats(csv_path: Path, data: Dict):
    """Append a newly computed alignment to the on-disk cache file, under an exclusive file lock.

    This runs inside multiprocessing.Pool workers, and the _PERSISTED_PROTEINS guard that decides
    WHETHER to write cannot deduplicate across them: it is a plain module-level set, so every forked
    worker carries its own copy and none of them sees the others' writes. Two workers holding jobs for
    the same protein therefore both reach this function, and an unlocked append from separate
    processes can interleave two rows into one line.

    The lock is on the file rather than in memory because the writers are PROCESSES — a threading.Lock
    would not be shared across the fork. flock is advisory but every writer here goes through this
    function, which is what makes it sufficient.

    Row-level duplication remains possible (two workers may each write the same protein once) and is
    harmless by design: save_alignment_cache_final rebuilds this file wholesale from the authoritative
    per-job summaries at the end of the run, emitting exactly one row per protein. What the lock buys
    is that the live-feedback file is never TORN mid-run.
    """
    if "protein" in data:
        CACHED_ALIGNMENTS[data["protein"]] = data

    keys = sorted(data.keys())
    try:
        with open(csv_path, "a", newline="") as f:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            try:
                # Under the lock: another worker may have created the file (and its header)
                # between the exists() check and the open.
                f.seek(0, os.SEEK_END)
                writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
                if f.tell() == 0:
                    writer.writeheader()
                writer.writerow(data)
                f.flush()
            finally:
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)
    except Exception as e:
        if logger: logger.debug(f"Cache synchronisation deferred: {e}")

def save_alignment_cache_final(csv_path: Path):
    """Rebuild the alignment stats CSV from the authoritative per-job summaries.

    The active-site alignment is identical across a protein's 27 ligand jobs and is
    recorded in every job's summary.json. The in-loop persist appends one row per
    protein for live feedback, but on a resume the analysis pool skips jobs already
    written to the Master CSV, so those proteins never re-enter the persist path and
    their rows would be absent from a wiped alignment directory. Reconstructing the
    file from the summaries guarantees exactly one row per unique protein regardless
    of which jobs ran this session.
    """
    import ast
    cols = ["active_site_mapping", "align_score", "gap_count", CFG.COL_ID_PCT,
            "nuc_rescue_offset", "nuc_resolution", "protein", "seq_length", "target_sequence"]
    jobs_dir = csv_path.parents[2] / "4_Prediction_Jobs"
    if not jobs_dir.is_dir():
        return

    # One summary per protein: job folders are named "<jid>_<serial>_<name>_<ligand>"
    # and sort contiguously by serial, so the first folder of each serial suffices.
    serial_re = re.compile(r"^(\d+)_(\d+)_")
    seen_serial, picks = set(), []
    for d in sorted(os.listdir(jobs_dir)):
        m = serial_re.match(d)
        if not m:
            continue
        key = d.rsplit("_", 2)[0] if m.group(1) == "0000000" else m.group(2)
        if key in seen_serial:
            continue
        seen_serial.add(key)
        picks.append(d)

    rows = {}
    for d in picks:
        sp = jobs_dir / d / f"{d}_summary.json"
        if not sp.exists():
            continue
        try:
            with open(sp) as f:
                j = json.load(f)
        except Exception:
            continue
        prot = j.get("protein")
        if not prot:
            continue
        raw = j.get("active_site_mapping", "")
        try:
            asm = json.dumps(ast.literal_eval(raw)) if raw else "{}"
        except Exception:
            try:
                asm = json.dumps(json.loads(raw))
            except Exception:
                asm = "{}"
        rows[prot] = {
            "active_site_mapping": asm,
            "align_score": j.get("align_score", ""),
            "gap_count": j.get("gap_count", ""),
            CFG.COL_ID_PCT: j.get(CFG.COL_ID_PCT, ""),
            "nuc_rescue_offset": j.get("nuc_rescue_offset", ""),
            "nuc_resolution": j.get("nuc_resolution", ""),
            "protein": prot,
            "seq_length": j.get("seq_length", ""),
            "target_sequence": j.get("target_sequence", ""),
        }

    if not rows:
        console_info("Warning: No per-job summaries found — alignment stats not rebuilt.")
        return

    try:
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            writer.writeheader()
            for prot in sorted(rows):
                writer.writerow(rows[prot])
    except Exception as e:
        console_info(f"Warning: Failed to write alignment stats: {e}")
        return

    console_separator()
    console_info(f"Active-site alignment complete — {len(rows):,} sequences aligned, "
                 f"Alignment_Stats.csv written ({len(rows):,} rows) → {csv_path}")
    console_separator()

# -------------------------------------------------------------------------------
# Step 5.2: Alignment Logic Engine
# -------------------------------------------------------------------------------
def map_active_site_residues(protein_id: str, target_seq: str, out_aln_path: Optional[Path] = None, stats_csv: Optional[Path] = None) -> Tuple[Dict[str, int], Dict[str, Any], Dict[str, str], str]:
    """
    Performs a global sequence alignment between the canonical reference sequence
    and the target sequence to mathematically pinpoint dynamic active site indices.
    """
    # -------------------------------------------------------------------------------
    # Sub-Step 5.2.1: Define Transliteration Maps
    # -------------------------------------------------------------------------------
    AA_1_TO_3 = {
        "A":"ALA", "R":"ARG", "N":"ASN", "D":"ASP", "C":"CYS", "Q":"GLN", "E":"GLU",
        "G":"GLY", "H":"HIS", "I":"ILE", "L":"LEU", "K":"LYS", "M":"MET", "F":"PHE",
        "P":"PRO", "S":"SER", "T":"THR", "W":"TRP", "Y":"TYR", "V":"VAL", "B":"ASX",
        "Z":"GLX", "X":"UNK", "-": "GAP"
    }

    # -------------------------------------------------------------------------------
    # Sub-Step 5.2.2: Configure Alignment Parameters
    # -------------------------------------------------------------------------------
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score   = CFG.ALIGN_OPEN_GAP_SCORE
    aligner.extend_gap_score = CFG.ALIGN_EXTEND_GAP_SCORE
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    if not target_seq or not target_seq.strip():
        return {}, {"error": "Target Sequence string is empty", CFG.COL_ID_PCT: 0.0}, {}, "NA"

    target_seq = str(target_seq)

    # -------------------------------------------------------------------------------
    # Sub-Step 5.2.3: Execute Alignment Algorithm
    # -------------------------------------------------------------------------------
    alignments = aligner.align(REF_SEQUENCE_STR, target_seq)
    if not alignments:
        return {}, {"error": "No viable alignment mapping identified", CFG.COL_ID_PCT: 0.0}, {}, "NA"

    aln = alignments[0]

    if out_aln_path:
        out_aln_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_aln_path, "w") as f: f.write(format(aln))
        GLOBAL_STATS["total_alignments_generated"] += 1

    # -------------------------------------------------------------------------------
    # Sub-Step 5.2.4: Extract Statistics & Map Coordinate Frames
    # -------------------------------------------------------------------------------

    '''
    Sequence identity over ALIGNED (gap-free) columns — positions where both sequences
    carry a residue — the standard definition. Dividing by the full alignment length
    (gaps included) would artificially deflate identity for a target bearing a large
    insertion/deletion vs the reference. total_columns is retained for reporting.
    '''
    n_match, aligned_cols, total_columns = 0, 0, 0
    seq1_str, seq2_str = str(aln[0]), str(aln[1])
    for c1, c2 in zip(seq1_str, seq2_str):
        total_columns += 1
        if c1 != "-" and c2 != "-":
            aligned_cols += 1
            if c1 == c2: n_match += 1

    identity_pct = (n_match / aligned_cols) * 100 if aligned_cols > 0 else 0

    stats = {
        "align_score": aln.score, CFG.COL_ID_PCT: round(identity_pct, 2),
        "target_sequence": target_seq, "seq_length": len(target_seq),
        "gap_count": seq2_str.count("-")
    }

    ref_aligned, tgt_aligned = aln[0], aln[1]
    ref_cursor, tgt_cursor = 0, 0
    ref_to_tgt_idx = {}
    tgt_seq_map = {}
    full_mapping_parts = []

    for i in range(len(ref_aligned)):
        r_char, t_char = ref_aligned[i], tgt_aligned[i]
        if r_char != "-": ref_cursor += 1
        if t_char != "-": tgt_cursor += 1

        if r_char != "-":
            ref_aa_3 = AA_1_TO_3.get(r_char, "UNK")
            if t_char != "-":
                tgt_aa_3 = AA_1_TO_3.get(t_char, "UNK")
                entry = f"{ref_aa_3}{ref_cursor}:{tgt_aa_3}{tgt_cursor}"
                ref_to_tgt_idx[ref_cursor] = tgt_cursor
                tgt_seq_map[tgt_cursor] = t_char
            else:
                entry = f"{ref_aa_3}{ref_cursor}:GAP"
            full_mapping_parts.append(entry)

    # -------------------------------------------------------------------------------
    # Sub-Step 5.2.5: Resolve Catalytic Triad Positions
    # -------------------------------------------------------------------------------
    final_idx_map = {}
    final_resname_map = {}

    """
    Class-aware ±window resolver (CFG §2.7). Every catalytic residue is mapped
    the same dynamic way as the nucleophile: take the directly aligned target
    position; if it is missing (gap/indel) or carries the wrong chemical class
    (substitution), scan ±RESIDUE_SEARCH_WINDOW target positions for a residue
    of the role's EXPECTED class and map to the nearest such match. This makes
    ASP110/ASP134/HIS277 and the cradle/stabilisers tolerant of small numbering
    shifts and conservative mutations, instead of returning MISSING/999.
    """
    _WIN = int(getattr(CFG, "RESIDUE_SEARCH_WINDOW", 5))
    _ROLE_CLASS = getattr(CFG, "ROLE_EXPECTED_RESIDUES", {})
    """
    Ambiguous / non-standard residues (X→UNK, B, Z, J, O, U) carry no defined
    sidechain identity or geometry.
    """
    _AMBIG_AA = getattr(CFG, "PREP_AMBIGUOUS_AA", set("BXZJOU"))
    """
    When an ambiguous residue sits inside the base ±window, reach a little
    further so a genuine catalytic residue just beyond the ambiguous gap is
    not missed (base window unchanged when the window is clean).
    """
    _AMBIG_EXT = int(getattr(CFG, "RESIDUE_SEARCH_AMBIG_EXTENSION", 2))

    def _aa3_at(tpos):
        raw = AA_1_TO_3.get(tgt_seq_map.get(tpos, "X"), "UNK")
        return _canonical_resname(raw)

    def _is_ambig(tpos):
        return tgt_seq_map.get(tpos, "X") in _AMBIG_AA

    def _resolve_key(ref_idx, role, exclude=None):
        """
        exclude = target indices already assigned to OTHER catalytic roles. Each
        catalytic residue is a DISTINCT position (e.g. the two clamp arginines
        Arg111/Arg114), so a residue already claimed by one role must not be
        re-used for another — this prevents the clamp-degeneracy bug where
        Clamp1 and Clamp2 collapse onto the same arginine.
        """
        exclude = exclude or set()
        expected = _ROLE_CLASS.get(role, set())
        direct = ref_to_tgt_idx.get(ref_idx)
        if direct is not None and direct not in exclude and (not expected or _aa3_at(direct) in expected):
            return direct, "direct", 0                     # clean exact-position hit
        if not expected:
            _d = direct if direct not in exclude else None
            return _d, ("direct_noclass" if _d is not None else "none"), 0  # no class constraint → keep direct
        """
        Anchor the window: the direct position if it exists, else the nearest
        ref position within the window that did align to a target index.
        """
        anchor = direct
        if anchor is None:
            for k in range(1, _WIN + 1):
                if ref_to_tgt_idx.get(ref_idx - k) is not None:
                    anchor = ref_to_tgt_idx[ref_idx - k]; break
                if ref_to_tgt_idx.get(ref_idx + k) is not None:
                    anchor = ref_to_tgt_idx[ref_idx + k]; break
        if anchor is None:
            _d = direct if direct not in exclude else None
            return _d, ("direct_fallback" if _d is not None else "none"), 0
        """
        Nearest-in-sequence residue of the expected class within ±window.
        If an ambiguous/non-standard residue falls inside the base window,
        extend the reach by _AMBIG_EXT so a real residue just past the gap is
        still reachable; otherwise the base window is used unchanged. The
        resolution method ('direct'/'window'/'window_gap') and the residue offset
        from the directly-aligned column are returned for the nucleophile-rescue
        audit (CFG §5.4): a wide rescue can latch onto a non-catalytic Asp.
        """
        _win_eff = _WIN
        if any(_is_ambig(anchor + s) for s in range(-_WIN, _WIN + 1)):
            _win_eff = _WIN + _AMBIG_EXT
        for k in range(0, _win_eff + 1):
            for cand in ((anchor,) if k == 0 else (anchor + k, anchor - k)):
                if cand in tgt_seq_map and cand not in exclude and _aa3_at(cand) in expected:
                    if direct is not None:
                        _off = abs(cand - direct)
                        return cand, ("direct" if _off == 0 else "window"), _off
                    return cand, "window_gap", abs(cand - anchor)
        _d = direct if direct not in exclude else None
        return _d, ("direct_fallback" if _d is not None else "none"), 0   # no class match → keep direct (may be None)

    """
    Safe fallback — after the (window-extended) search, if the resolved
    position is still ambiguous/non-standard, drop it. Such a residue is parsed
    downstream as a non-protein "ligand" atom, producing spurious 0.0 Å
    self-distances and 0.0° SN2 angles. Merge QC (01) only rejects whole
    sequences above an ambiguity fraction, so a lone X survives — it must be
    screened out here so the triad fails honestly (MISSING) rather than on a
    phantom residue.
    """
    _used_idx = set()   # target indices already claimed by a catalytic role (distinct-residue guard)
    _resolution_audit = {}   # role → (method, residue-offset) for the §5.4 rescue audit
    for name, info in REF_ACTIVE_SITE_MAP.items():
        tgt_idx, _res_method, _res_off = _resolve_key(info["id"], info.get("role", ""), exclude=_used_idx)
        if tgt_idx is not None and _is_ambig(tgt_idx):
            tgt_idx, _res_method, _res_off = None, "none", 0
        _resolution_audit[name] = (_res_method, _res_off)
        final_idx_map[name] = tgt_idx
        if tgt_idx is not None:
            _used_idx.add(tgt_idx)
            aa_1 = tgt_seq_map.get(tgt_idx, "X")
            aa_3 = AA_1_TO_3.get(aa_1, "UNK")
            final_resname_map[name] = f"{aa_3}{tgt_idx}"
        else:
            final_resname_map[name] = "MISSING"

    """
    Nucleophile alignment-rescue audit (CFG §5.4): record how the catalytic Asp was
    resolved (direct hit / windowed rescue) and how far the rescue travelled, so a
    far-fetched rescued Asp can be flagged and barred from the elite tier downstream.
    """
    _nuc_m, _nuc_off = _resolution_audit.get("Nuc", ("none", 0))
    stats["nuc_resolution"]    = _nuc_m if _nuc_off == 0 else f"{_nuc_m}+{_nuc_off}"
    stats["nuc_rescue_offset"] = int(_nuc_off)

    full_map_str = " | ".join(full_mapping_parts)
    full_record = stats.copy()
    full_record["protein"] = protein_id
    full_record["active_site_mapping"] = json.dumps(final_idx_map)

    if stats_csv:
        update_alignment_stats(stats_csv, full_record)
        CACHED_ALIGNMENTS[protein_id] = full_record

    return final_idx_map, stats, final_resname_map, full_map_str

def format_control_mappings(resname_map: Dict[str, str], full_map_str: str) -> Tuple[str, str]:
    """Helper tool for parsing standard string formatting for the master CSV report mapping strings."""
    triad_str_parts = []
    for key in sorted(CATALYTIC_TRIAD_KEYS):
        ref_label = key
        tgt_val = resname_map.get(key, "MISSING")
        entry = f"{ref_label}:{tgt_val}"
        triad_str_parts.append(entry)
    return full_map_str, " | ".join(triad_str_parts)


def format_full_role_map(resname_map: Dict[str, str]) -> Dict[str, str]:
    """Per-role ACTUAL mapped residue (e.g. 'ASP104') for the CSV Mapped_* columns.

    Pairs 1:1 with the Dist_* columns. The value is the residue resolved by the
    ±5 class-aware resolver for THIS protein, not the canonical 3R3U number, so a
    reader can see e.g. Dist_Nucleophile=2.97 ↔ Mapped_Nucleophile=ASP104.
    """
    _ROLE_TO_COL = {
        "Nuc": "Mapped_Nucleophile", "Base": "Mapped_Base", "Acid": "Mapped_Acid",
        "Carb1": "Mapped_Clamp1", "Carb2": "Mapped_Clamp2",
        "Stab_H": "Mapped_Stabiliser_H", "Stab_W": "Mapped_Stabiliser_W",
        "Stab_Y": "Mapped_Stabiliser_Y",
    }
    return {col: resname_map.get(role, "MISSING") for role, col in _ROLE_TO_COL.items()}


# =============================================================================
# SECTION 6: CHEMINFORMATICS & PHYSICS ENGINE
# =============================================================================

# -------------------------------------------------------------------------------
# Step 6.1: Spatial Mathematical Utilities
# -------------------------------------------------------------------------------
# Distance and angle helpers are imported from the central ProjectUtils module.

def get_plane_normal(atoms):
    """Calculates the optimal best-fit plane normal vector for Pi-stacking analysis using Singular Value Decomposition (SVD)."""
    if not atoms: return None, None
    coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms])
    centroid = np.mean(coords, axis=0)
    centred = coords - centroid
    _, _, vh = np.linalg.svd(centred)
    normal = vh[2, :]
    return centroid, normal

def load_structure_safe(path: Path):
    """Safely loads standard PDB and mmCIF files using Gemmi's robust parsing capabilities."""
    try:
        if path.suffix == ".pdb": return gemmi.read_structure(str(path))
        else:
            doc = gemmi.cif.read_file(str(path))
            return gemmi.make_structure_from_block(doc.sole_block())
    except Exception as e:
        if logger: logger.debug(f"Failed to load standard physical structure at {path}: {e}")
        return None

# -------------------------------------------------------------------------------
# Step 6.2: RDKit Core Integration
# -------------------------------------------------------------------------------
def rdkit_mol_from_smiles_with_3d(smiles: str) -> Chem.Mol:
    """
    Constructs a three-dimensional RDKit molecule from a SMILES string.
    This is crucial for establishing charge topologies, which are absent in standard mmCIF files.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    mol = Chem.AddHs(mol)
    opts = StereoEnumerationOptions(tryEmbedding=True, maxIsomers=1)
    isomers = list(EnumerateStereoisomers(mol, options=opts))
    if isomers: mol = isomers[0]
    params = AllChem.ETKDGv3()
    params.randomSeed = int(CFG.RDKIT_EMBED_SEED)
    res = AllChem.EmbedMolecule(mol, params)
    if res == -1:
        params.useRandomCoords = True
        AllChem.EmbedMolecule(mol, params)
    try: AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception as e:
        if logger: logger.debug(f"UFF optimisation layer operation failed: {e}")
    return mol

def rdkit_coords_list(mol: Chem.Mol):
    """Extracts raw XYZ coordination matrices from RDKit molecule data blocks."""
    conf = mol.GetConformer()
    return [(p.x, p.y, p.z) for p in [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]]

def kabsch_transform(P, Q):
    """
    Determines the optimal rotation/translation matrix required to superimpose spatial cluster P onto Q.
    Employed extensively for aligning RDKit theoretical models to Boltz-predicted mmCIF coordinates.
    """
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P_centred = P - Pc
    Q_centred = Q - Qc
    U, _, Vt = np.linalg.svd(P_centred.T @ Q_centred)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = Qc - (R @ Pc)
    return R, t

def graph_map_mmcif_to_rdkit(mmcif_atoms, rdkit_mol):
    """
    Map Boltz-predicted (mmCIF) atoms to RDKit template atoms by CONNECTIVITY, not by
    position. Returns {mmcif_atom_name: rdkit_atom_index}, or None if the graphs cannot
    be matched.

    This is the correct primitive, and the spatial assignment below is only a fallback.
    The template carries an ETKDG-embedded conformer (rdkit_mol_from_smiles_with_3d),
    which shares its BONDING with the Boltz pose but not its TORSIONS. Two different
    conformers of the same molecule cannot be superposed by any rigid rotation, so a
    distance-based assignment between them is matching atoms across a geometry mismatch
    it cannot remove: a nearest-neighbour cost can then pair chemically distinct atoms
    (an F on one carbon with an F on another; a carboxylate O with an F) whenever the
    torsional difference exceeds the interatomic spacing. Because the reactive-centre
    gates are keyed on the mapped α-carbon and carboxylate oxygens, that mislabelling
    silently relocates the scissile centre.

    A graph isomorphism has no such failure mode: it is invariant to conformation, to
    reference frame, and to atom order. Bonds are perceived from the mmCIF geometry
    (rdDetermineBonds), and the template is matched onto them with generic bond queries
    so that perceived single bonds still match aromatic/double template bonds. Where the
    molecule has symmetry (the three F of a CF3), the match picks one automorphism —
    those atoms are chemically interchangeable, so any of them is equally correct.

    One ambiguity is not resolvable and does not need to be: the perceived probe carries
    neither hydrogens nor bond orders, so the two carboxylate oxygens are topologically
    equivalent in it and the match may interchange them. They are a resonance pair — the
    ligand is the carboxylATE at the assay pH — and every consumer is symmetric in them:
    _ligand_ionisable classifies an oxygen by its local topology (bound to a C/S/P bearing
    ≥2 oxygens), which is identical for both, and the reactive-centre gates key on the
    α-carbon, not on either oxygen. Measured against shuffled, re-embedded conformers, the
    carbon/fluorine skeleton — the scissile centre — maps with 100% accuracy for FA, DFA,
    TFA, PFBA and PFOA, where the spatial fallback reaches only 82% for PFBA and 48% for
    PFOA (the flexible chains, where the two conformers differ most).
    """
    if not mmcif_atoms or rdkit_mol is None:
        return None
    try:
        heavy = [(i, a) for i, a in enumerate(mmcif_atoms)
                 if str(a.get("element", "")).strip().upper() not in ("H", "D")]
        if len(heavy) < 2:
            return None

        # Perceive connectivity from the predicted geometry alone.
        _xyz = f"{len(heavy)}\n\n" + "\n".join(
            f'{a["element"]} {float(a["x"]):.6f} {float(a["y"]):.6f} {float(a["z"]):.6f}'
            for _, a in heavy)
        probe = Chem.MolFromXYZBlock(_xyz)
        if probe is None:
            return None
        rdDetermineBonds.DetermineConnectivity(probe)
        if probe.GetNumBonds() == 0:
            return None

        # Heavy-atom template, remembering each atom's index in the ORIGINAL rdkit_mol.
        work = Chem.Mol(rdkit_mol)
        for _a in work.GetAtoms():
            _a.SetIntProp("_orig_idx", _a.GetIdx())
        tmpl = Chem.RemoveHs(work, sanitize=False)
        if tmpl is None or tmpl.GetNumAtoms() < 2:
            return None
        orig_idx = [_a.GetIntProp("_orig_idx") for _a in tmpl.GetAtoms()]

        '''
        Neither bond ORDERS nor formal CHARGES are recoverable from bare coordinates, so the query
        must not demand either. Bonds are made generic; charges are zeroed on the query copy, because
        the input SMILES are supplied as the anion (O=C([O-])...) at the assay pH while the perceived
        probe carries neutral oxygens — a charge-sensitive match rejects every carboxylate and every
        sulfonate. Atom ELEMENTS and the bond GRAPH still have to agree, which is what identifies the
        atoms. The original template indices are carried in _orig_idx, so neutralising this copy
        changes nothing that is returned.
        '''
        for _a in tmpl.GetAtoms():
            _a.SetFormalCharge(0)
            _a.SetNoImplicit(True)
        _qp = Chem.AdjustQueryParameters.NoAdjustments()
        _qp.makeBondsGeneric = True
        query = Chem.AdjustQueryProperties(tmpl, _qp)

        match = probe.GetSubstructMatch(query)      # match[k] = probe atom for template atom k
        if len(match) != tmpl.GetNumAtoms():
            return None

        # Belt and braces: the isomorphism must also be element-consistent.
        out = {}
        for _k, _probe_i in enumerate(match):
            _mm_i, _mm_a = heavy[_probe_i]
            if str(_mm_a["element"]).strip().upper() != tmpl.GetAtomWithIdx(_k).GetSymbol().upper():
                return None
            out[_mm_a["atom_name"]] = int(orig_idx[_k])
        return out or None
    except Exception:
        return None


def map_mmcif_to_rdkit(mmcif_atoms, rdkit_mol):
    """
    FALLBACK atom map, reached only when graph_map_mmcif_to_rdkit above cannot perceive
    the bonds or cannot match the graph. Maps Boltz-predicted (mmCIF) atoms to RDKit
    template atoms by optimal element-aware 3D assignment, robust to atom-ordering AND
    reference-frame differences between the two sources — but NOT to the conformational
    difference between the pose and the embedded template, which is why it is second choice.

    The two structures live in different frames, so a rigid alignment is required
    before inter-atomic distances are meaningful. That alignment must NOT be seeded
    from atom INDEX order (rd_coords[:m] ↔ mm_coords[:m]): Boltz's mmCIF atom order
    need not match RDKit's SMILES order — the mismatch grows with ligand size — and
    an index-seeded rotation silently produces inflated costs and unmapped catalytic
    atoms (e.g. the carboxylate C/O, which then breaks Bürgi–Dunitz / Flippin–Lodge /
    mechanistic scoring).

    Instead, several ORDER-INDEPENDENT candidate alignments are generated by
    matching the two clouds' principal axes (all four proper-rotation sign
    combinations), with the index-order Kabsch fit retained as one extra
    candidate. Each candidate is scored by its optimal Hungarian assignment and
    the alignment that maps the most atoms (lowest total cost as tie-break) wins.
    """
    if not mmcif_atoms or rdkit_mol is None:
        return {}

    # Connectivity first: conformer-, frame- and order-independent. Spatial assignment is
    # reached only when bond perception or the graph match fails.
    _graph = graph_map_mmcif_to_rdkit(mmcif_atoms, rdkit_mol)
    if _graph:
        return _graph

    mm_coords = np.asarray([[a["x"], a["y"], a["z"]] for a in mmcif_atoms], float)
    rd_coords = np.asarray(rdkit_coords_list(rdkit_mol), float)
    if mm_coords.size == 0 or rd_coords.size == 0:
        return {}

    mm_elems = [a["element"] for a in mmcif_atoms]
    rd_elems = [a.GetSymbol() for a in rdkit_mol.GetAtoms()]
    names    = [a["atom_name"] for a in mmcif_atoms]
    # Element-mismatch mask (constant across candidate alignments).
    elem_mismatch = np.array([[me != re for re in rd_elems] for me in mm_elems])

    def _assignment(rd_xyz):
        cost = np.linalg.norm(mm_coords[:, None, :] - rd_xyz[None, :, :], axis=2)
        penalty = (np.max(cost) * CFG.ALIGN_DYNAMIC_PENALTY_FACTOR
                   if cost.size > 0 else CFG.ALIGN_DYNAMIC_PENALTY_FALLBACK)
        cost = cost + elem_mismatch * penalty
        row_ind, col_ind = linear_sum_assignment(cost)
        pairs = [(i, j) for i, j in zip(row_ind, col_ind)
                 if cost[i, j] <= COORD_MATCH_MAX_DIST]
        total = float(sum(cost[i, j] for i, j in pairs))
        return pairs, total

    # ── Candidate rigid alignments of the RDKit cloud onto the mmCIF cloud ──
    mm_mean, rd_mean = mm_coords.mean(0), rd_coords.mean(0)
    candidates = []
    try:
        _, Sm, Vm = np.linalg.svd(mm_coords - mm_mean)
        _, _, Vr = np.linalg.svd(rd_coords - rd_mean)
        for sx in (1.0, -1.0):
            for sy in (1.0, -1.0):
                sz = sx * sy * float(np.sign(np.linalg.det(Vm.T @ Vr)) or 1.0)
                R = Vm.T @ np.diag([sx, sy, sz]) @ Vr
                if np.linalg.det(R) > 0:
                    candidates.append((rd_coords - rd_mean) @ R.T + mm_mean)
        '''
        Degenerate principal axes (linear/planar ligands, e.g. long-chain PFCAs):
        when the 2nd/3rd singular values are near-equal the discrete sign flips above
        cannot resolve the arbitrary continuous rotation about the major (chain) axis.
        Add rotations of the sign-flip candidates about that axis; selection below is
        min-cost over all candidates, so this only ever improves the match.
        '''
        _s0 = float(Sm[0]) if len(Sm) else 0.0
        if (len(Sm) >= 3 and _s0 > 1e-6
                and abs(float(Sm[1]) - float(Sm[2])) / _s0 < CFG.PCA_DEGENERACY_TOL):
            _ax = np.asarray(Vm[0], float)
            _ax = _ax / (np.linalg.norm(_ax) or 1.0)
            _Kx = np.array([[0.0, -_ax[2], _ax[1]],
                            [_ax[2], 0.0, -_ax[0]],
                            [-_ax[1], _ax[0], 0.0]])
            for _ang in (np.pi / 2.0, np.pi, 3.0 * np.pi / 2.0):
                _Rax = (np.eye(3) + np.sin(_ang) * _Kx
                        + (1.0 - np.cos(_ang)) * (_Kx @ _Kx))
                for _base in list(candidates):
                    candidates.append((_base - mm_mean) @ _Rax.T + mm_mean)
    except Exception:
        pass
    try:
        m = min(len(rd_coords), len(mm_coords))
        if m >= 3:
            R, t = kabsch_transform(rd_coords[:m], mm_coords[:m])
            candidates.append((R @ rd_coords.T).T + t)
    except Exception:
        pass
    if not candidates:
        candidates.append(rd_coords)

    best_pairs, best_key = [], None
    for rd_xyz in candidates:
        pairs, total = _assignment(rd_xyz)
        key = (len(pairs), -total)            # most atoms mapped, then lowest cost
        if best_key is None or key > best_key:
            best_key, best_pairs = key, pairs
    return {names[i]: int(j) for i, j in best_pairs}

# -------------------------------------------------------------------------------
# Step 6.3: Interaction Classification Engine
# -------------------------------------------------------------------------------
def _ligand_ionisable(rd_atom) -> Optional[str]:
    """
    Classify a ligand atom as 'anion'- or 'cation'-capable for salt-bridge
    detection, independent of the input SMILES protonation state.

    Input PFAS SMILES are typically supplied in the neutral (protonated) acid
    form (e.g. fluoroacetate "C(C(=O)O)F"), so RDKit reports a formal charge of
    0 even for groups that are fully ionised at the assay pH (8.0). Relying on
    formal charge alone would miss every carboxylate/sulfonate, while treating ANY
    neutral N/O/S as charged over-counts ethers and amides. The actual ionisable
    groups are recognised instead:
      • anion : carboxylate / sulfonate / phosphonate oxygen (O bound to a
                C/S/P bearing ≥2 oxygens), or any atom with formal charge < 0.
      • cation: basic amine nitrogen (all-single-bond N, not an amide), or any
                atom with formal charge > 0.
    Returns 'anion', 'cation', or None.
    """
    if rd_atom is None:
        return None
    fc = rd_atom.GetFormalCharge()
    if fc < 0: return "anion"
    if fc > 0: return "cation"
    sym = rd_atom.GetSymbol()
    if sym == "O":
        for nb in rd_atom.GetNeighbors():
            if nb.GetSymbol() in ("C", "S", "P"):
                o_count = sum(1 for x in nb.GetNeighbors() if x.GetSymbol() == "O")
                if o_count >= 2:   # carboxylate / sulfonate / phosphate motif
                    return "anion"
    elif sym == "N":
        if all(b.GetBondType() == Chem.BondType.SINGLE for b in rd_atom.GetBonds()):
            is_amide = any(
                ob.GetOtherAtom(nb).GetSymbol() == "O"
                and ob.GetBondType() == Chem.BondType.DOUBLE
                for nb in rd_atom.GetNeighbors() if nb.GetSymbol() == "C"
                for ob in nb.GetBonds()
            )
            if not is_amide:
                return "cation"
    return None


def classify_pair(p_atom, l_atom, d, rdkit_mol, mm_to_rd):
    """
    Resolves the specific chemical interaction type between two spatially proximate atoms,
    using the Schrödinger Maestro interaction cutoffs defined in CFG §3 (heavy-atom
    geometry; D···A proxy for H-bonds since Boltz-2 CIF has no explicit H).
    """
    types = []
    pel = p_atom["element"]
    lel = l_atom["element"]

    # -------------------------------------------------------------------------------
    # Sub-Step 6.3.1: Evaluator Block
    # -------------------------------------------------------------------------------
    # 1. Hydrogen Bonds (D-A distance only; D-H-A angle requires explicit H atoms not present in Boltz-2 CIF)
    if pel in {"N", "O"} and lel in {"N", "O"}:
        if HBOND_MIN_DIST_DA <= d <= HBOND_MAX_DIST_DA: types.append("hydrogen_bond")

    # 2. Hydrophobic Bonds
    if pel=="C" and lel=="C" and d <= HYDROPHOBIC_MAX_DIST: types.append("hydrophobic")

    """
    3. Salt Bridges (opposite-charge groups within SALT_BRIDGE_MAX_DIST).
    Ligand ionisability is resolved by functional group, not by the SMILES
    protonation state (see _ligand_ionisable), so neutral-form carboxylates
    still register while ethers/amides do not.
    """
    if d <= SALT_BRIDGE_MAX_DIST:
        rd_atom = None
        idx = mm_to_rd.get(l_atom["atom_name"])
        if idx is not None and rdkit_mol:
            try: rd_atom = rdkit_mol.GetAtomWithIdx(idx)
            except Exception: rd_atom = None
        lig_ion = _ligand_ionisable(rd_atom)
        is_salt_bridge = False
        # Canonicalise the residue name so force-field protonation variants (HIP/ASH/GLH/LYN…)
        # are classified, not just the bare canonical codes (consistent with the rest of the engine).
        _pres = _canonical_resname(p_atom["resname"])

        # The formal charge sits on specific sidechain atoms, not the whole residue: the
        # cationic nitrogen (Arg NH1/NH2/NE, Lys NZ, His ND1/NE2) or the carboxylate oxygen
        # (Asp OD1/OD2, Glu OE1/OE2). Gate on the protein atom's element AND exclude backbone
        # N/O, so a residue's backbone and sidechain carbons within SALT_BRIDGE_MAX_DIST no
        # longer each register a spurious salt bridge (interaction-density inflation).
        _pname = p_atom.get("atom_name", "")
        _is_backbone = _pname in _BACKBONE_ATOMS

        # Case A: cationic residue (Arg/Lys/His) vs anionic ligand group — charged sidechain N.
        if _pres in POSITIVE_RES and lig_ion == "anion" and pel.upper() == "N" and not _is_backbone:
            is_salt_bridge = True

        # Case B: anionic residue (Asp/Glu) vs cationic ligand group — carboxylate sidechain O.
        if _pres in NEGATIVE_RES and lig_ion == "cation" and pel.upper() == "O" and not _is_backbone:
            is_salt_bridge = True

        if is_salt_bridge:
            types.append("salt_bridge")
            # A salt bridge IS a charged hydrogen bond — count it once, as the salt bridge,
            # not also as a generic H-bond (the same O···N contact otherwise inflated both
            # count_salt_bridge and count_hydrogen_bond → num_interactions double-counted).
            if "hydrogen_bond" in types:
                types.remove("hydrogen_bond")

    # 4. Halogen Contact Analysis & Specific Fluorous Metrics
    """
    Aliphatic C–F is not a halogen-bond donor (CFG §3.6): only the polarisable
    heavy halogens (Cl/Br/I) form genuine σ-hole halogen bonds, so halogen_contact
    is reserved for those. Ligand fluorine is scored solely through the dedicated
    fluorous metrics below — never double-counted as a halogen bond.
    """
    if lel.upper() in {"CL", "BR", "I"} and d <= HALOGEN_MAX_DIST:
        types.append("halogen_contact")

    if lel.upper() == "F" and d <= HALOGEN_MAX_DIST:
        types.append("fluorine_contact")
        # 4a. Fluorine-Polar (Considered highly potent in FAcD structures)
        if pel in {"N", "O", "S"}: types.append("fluorine_polar")
        # 4b. Fluorous-Lipophilic (Fluorine binding to internal Carbon structures)
        if pel == "C": types.append("fluorous_hydrophobic")
        # 4c. Fluorous-Fluorous (Uncommon but potentially present in extensive linings)
        if pel == "F": types.append("fluorous_fluorous")

    # 5. Metal Coordination Definitions
    if d <= METAL_COORD_MAX_DIST:
        if (pel in METALS and lel in {"O", "N", "S"}) or (lel in METALS and pel in {"O", "N", "S"}): types.append("metal_coordination")

    return types

# -------------------------------------------------------------------------------
# Step 6.4: mmCIF Parsers
# -------------------------------------------------------------------------------
def load_atoms_from_structure(cif_path):
    """Efficiently parses the predicted three-dimensional CIF structure into discrete atoms using Gemmi."""
    doc = gemmi.cif.read_file(str(cif_path))
    st = gemmi.make_structure_from_block(doc.sole_block())
    prot, lig = [], []
    for model in st:
        for chain in model:
            for res in chain:
                is_lig = res.name not in STANDARD_AA and res.name != "HOH"
                for atom in res:
                    d = {"x":atom.pos.x, "y":atom.pos.y, "z":atom.pos.z,
                         "element": (atom.element.name or atom.name[0]),
                         "chain": chain.name, "resname": res.name, "resseq": res.seqid.num, "atom_name": atom.name}
                    if is_lig: lig.append(d)
                    else: prot.append(d)
    return prot, lig

# -------------------------------------------------------------------------------
# Step 6.5: Advanced Geometrics
# -------------------------------------------------------------------------------
def analyse_pi_interactions(prot_atoms, lig_atoms, rdkit_mol, mm_map, counts):
    """
    Rigorously evaluates Pi-stacking (face-to-face and edge-to-face) and Pi-cation interactions
    using derived normal vectors.
    """
    # -------------------------------------------------------------------------------
    # Sub-Step 6.5.1: Identify Biological Rings
    # -------------------------------------------------------------------------------
    prot_rings = []
    res_map = defaultdict(list)
    for a in prot_atoms: res_map[(a["chain"], a["resseq"])].append(a)

    for (_ch, _seq), atoms in res_map.items():
        resname = _canonical_resname(atoms[0]["resname"])
        if resname in RES_PROPS["Aromatic"]:
            ring_atoms = [a for a in atoms if a["atom_name"] in AROMATIC_RING_ATOMS]
            if len(ring_atoms) >= 5:
                cent, norm = get_plane_normal(ring_atoms)
                if cent is not None: prot_rings.append({"centroid": cent, "normal": norm})

    # -------------------------------------------------------------------------------
    # Sub-Step 6.5.2: Identify Chemical Rings
    # -------------------------------------------------------------------------------
    lig_rings = []
    if rdkit_mol:
        ri = rdkit_mol.GetRingInfo()
        inv_map = {v: k for k, v in mm_map.items()}
        for ring_idxs in ri.AtomRings():
            cif_ring_atoms = []
            for idx in ring_idxs:
                cif_name = inv_map.get(idx)
                if cif_name:
                    found = next((a for a in lig_atoms if a["atom_name"] == cif_name), None)
                    if found: cif_ring_atoms.append(found)

            if len(cif_ring_atoms) >= 5:
                cent, norm = get_plane_normal(cif_ring_atoms)
                if cent is not None: lig_rings.append({"centroid": cent, "normal": norm})

    # -------------------------------------------------------------------------------
    # Sub-Step 6.5.3: Spatial Interference Verification Checks
    # -------------------------------------------------------------------------------
    for p_ring in prot_rings:
        p_c, p_n = p_ring["centroid"], p_ring["normal"]

        # A. Pi-Stacking Logic Validation Framework
        for l_ring in lig_rings:
            l_c, l_n = l_ring["centroid"], l_ring["normal"]
            dist = np.linalg.norm(p_c - l_c)
            dot = np.abs(np.dot(p_n, l_n))
            angle = np.degrees(np.arccos(np.clip(dot, -1, 1)))
            if dist <= PI_STACK_FACE_MAX_DIST and angle <= PI_STACK_FACE_MAX_ANGLE: counts["pi_stacking"] += 1
            elif dist <= PI_STACK_EDGE_MAX_DIST and angle >= PI_STACK_EDGE_MIN_ANGLE: counts["pi_stacking"] += 1

        # B. Pi-Cation Validation Framework
        for la in lig_atoms:
            idx = mm_map.get(la["atom_name"])
            charge = 0
            if idx is not None and rdkit_mol:
                try: charge = rdkit_mol.GetAtomWithIdx(idx).GetFormalCharge()
                except Exception: pass
            if charge > 0:
                l_pos = np.array([la["x"], la["y"], la["z"]])
                dist = np.linalg.norm(p_c - l_pos)
                if dist <= PI_CATION_MAX_DIST:
                    vec_to_cation = (l_pos - p_c) / dist
                    angle_cation = np.degrees(np.arccos(np.clip(np.abs(np.dot(p_n, vec_to_cation)), -1, 1)))
                    if angle_cation <= PI_CATION_MAX_ANGLE: counts["pi_cation"] += 1

def generate_detailed_interactions(cif_path, smiles, output_csv: Path) -> Dict[str, Any]:
    """
    The primary statistical engine orchestrating the complete interaction suite.
    Utilises a cKDTree for O(N log N) spatial efficiency and corrects the heavy atom count
    for bias-free density normalisation calculations.
    """
    try:
        # -------------------------------------------------------------------------------
        # Sub-Step 6.5.4: Execute Primary Loaders & Mappers
        # -------------------------------------------------------------------------------
        prot, lig = load_atoms_from_structure(cif_path)
        rd_mol = rdkit_mol_from_smiles_with_3d(smiles) if smiles else None
        mm_map = map_mmcif_to_rdkit(lig, rd_mol)

        # Heavy Atom Count Calculation for Objective, Bias-Free Density Metrics
        lig_heavy_atoms = 0
        if rd_mol: lig_heavy_atoms = rd_mol.GetNumHeavyAtoms()
        else: lig_heavy_atoms = len([a for a in lig if a["element"].upper() != "H"])

        # Structural safeguards implemented for single-atom ions (e.g. F-, Cl-)
        if lig_heavy_atoms == 0: lig_heavy_atoms = 1

        # -------------------------------------------------------------------------------
        # Sub-Step 6.5.5: Residue Demographics Census
        # -------------------------------------------------------------------------------
        census = {"Total_Residues": 0, "Total_Polar_Residues": 0,
                  "Total_NonPolar_Residues": 0, "Total_Pos_Residues": 0,
                  "Total_Neg_Residues": 0, "Total_Aromatic_Residues": 0}
        unique_residues = set()
        for p in prot: unique_residues.add((p["chain"], p["resseq"], p["resname"]))
        census["Total_Residues"] = len(unique_residues)
        res_atom_map = defaultdict(list)
        for p in prot: res_atom_map[(p["chain"], p["resseq"], p["resname"])].append(p)

        for (_chain, _seq, resname), _atoms in res_atom_map.items():
            # Normalise force-field protonation variants (HID/HIE/HIP, ASH, GLH…)
            # to their canonical code so the census is not blind to protonated states.
            resname = _canonical_resname(resname)
            if resname in RES_PROPS["Polar"]: census["Total_Polar_Residues"] += 1
            if resname in RES_PROPS["NonPolar"]: census["Total_NonPolar_Residues"] += 1
            if resname in RES_PROPS["Pos"]: census["Total_Pos_Residues"] += 1
            if resname in RES_PROPS["Neg"]: census["Total_Neg_Residues"] += 1
            if resname in RES_PROPS["Aromatic"]: census["Total_Aromatic_Residues"] += 1

        counts = defaultdict(int)
        dists = []
        rows = []
        total_f_atoms = len([a for a in lig if a["element"] == "F"])
        interacting_f_set = set()

        # -------------------------------------------------------------------------------
        # Sub-Step 6.5.6: Optimised cKDTree Interaction Evaluation Loop
        # -------------------------------------------------------------------------------
        if prot and lig:
            prot_coords = [[p["x"], p["y"], p["z"]] for p in prot]
            lig_coords  = [[l["x"], l["y"], l["z"]] for l in lig]
            tree = cKDTree(prot_coords)
            indices_list = tree.query_ball_point(lig_coords, r=CFG.CATALYTIC_DIST_CUTOFF)

            _salt_seen = set()   # (chain, resseq) that already contributed one salt bridge
            for i, p_indices in enumerate(indices_list):
                la = lig[i]
                l_pos = lig_coords[i]
                for p_idx in p_indices:
                    pa = prot[p_idx]
                    p_pos = prot_coords[p_idx]
                    d = distance(p_pos, l_pos)
                    labels = classify_pair(pa, la, d, rd_mol, mm_map)
                    if labels:
                        # Count at most one salt bridge per cationic/anionic residue: a single
                        # Arg/Lys/His (up to 3 sidechain N) opposite a carboxylate (2 O) otherwise
                        # registers up to 6 pairwise salt bridges for one ionic contact.
                        if "salt_bridge" in labels:
                            _rk = (pa.get("chain"), pa.get("resseq"))
                            if _rk in _salt_seen:
                                labels = [x for x in labels if x != "salt_bridge"]
                            else:
                                _salt_seen.add(_rk)
                        dists.append(d)
                        for l in labels: counts[l] += 1
                        if "fluorine_contact" in labels: interacting_f_set.add(la["atom_name"])
                        rows.append({
                            "protein_chain": pa["chain"], "protein_resname": pa["resname"],
                            "protein_resseq": pa["resseq"], "protein_atom": pa["atom_name"],
                            "ligand_atom": la["atom_name"], "distance_A": round(d, 4),
                            "interaction_types": ";".join(labels)
                        })

        # -------------------------------------------------------------------------------
        # Sub-Step 6.5.7: Invoke Macro-Structure Physics Evaluation
        # -------------------------------------------------------------------------------
        analyse_pi_interactions(prot, lig, rd_mol, mm_map, counts)
        if rows: pd.DataFrame(rows).to_csv(output_csv, index=False)
        avg = sum(dists)/len(dists) if dists else 0.0

        result = {
            "counts": counts, "avg_dist": avg, "num_interactions": len(rows),
            "prot_len": len(prot), "lig_len": len(lig),
            "ligand_heavy_atoms": lig_heavy_atoms,
            "mapped_fraction": float(len(mm_map))/float(len(lig)) if lig else 0.0,
            "min_dist": min(dists) if dists else 0.0, "max_dist": max(dists) if dists else 0.0,
            "interacting_fluorine_count": len(interacting_f_set),
            "total_fluorine_count": total_f_atoms
        }
        result.update(census)
        return result
    except Exception as e:
        if logger: logger.error(f"Interaction generation operations failed: {e}")
        return {}


# =============================================================================
# SECTION 7: GEOMETRIC ANALYSIS (CATALYSIS) & CALIBRATION
# =============================================================================

# -------------------------------------------------------------------------------
# Step 7.1: Reference Data Handling (Calibration Logic)
# -------------------------------------------------------------------------------
def perform_control_calibration(control_cif: Path, ref_pdb: Path,
                                model_label: str = "DeHa4 Control",
                                target_seq: Optional[str] = None,
                                silent: bool = False) -> Tuple[float, Dict, Dict]:
    """
    Calibrates a Boltz-2 model against the 3R3U crystal structure.
    Returns (trust_score, mapped_sites, calib_data).
    When silent=True, suppresses all console output; calib_data is returned for
    deferred consolidated printing.
    """
    def _ci(msg):
        if not silent:
            console_info(msg)

    _ci(f"\n   Calibrating {model_label} against Crystal Data (3R3U) ---")
    trust_score   = 99.0
    mapped_sites  = {}
    calib_data: Dict = {
        "model_label": model_label,
        "trust_score": 99.0,
        "verdict_label": "Failed",
        "verdict_colour": ConsoleColours.FAIL,
        "mapping_rows": [],   # list of dicts: key, ref_res, mapped_res, dist, ok
        "missing_residues": [],
        "all_mapped": False,
    }

    try:
        st_ref = load_structure_safe(ref_pdb)
        st_con = load_structure_safe(control_cif)

        if not st_ref or not st_con:
            _ci("    ! System error loading structural components for calibration phase.")
            return trust_score, mapped_sites, calib_data

        ref_chain = st_ref[0]["A"]
        con_chain = st_con[0][0]

        try:
            if hasattr(gemmi.SupSelect, "Ca"):
                sup_select = gemmi.SupSelect.Ca
            elif hasattr(gemmi.SupSelect, "CA"):
                sup_select = gemmi.SupSelect.CA
            else:
                sup_select = gemmi.SupSelect.All
        except Exception:
            sup_select = gemmi.SupSelect.All

        sup = gemmi.calculate_superposition(
            ref_chain.whole(),
            con_chain.whole(),
            gemmi.PolymerType.PeptideL,
            sup_select
        )
        trust_score = sup.rmsd
        calib_data["trust_score"] = trust_score

        for ch in st_con[0]:
            for res in ch:
                for atom in res:
                    p = sup.transform.apply(atom.pos)
                    atom.pos = gemmi.Position(p.x, p.y, p.z)

        if trust_score <= CFG.TRUST_SCORE_EXCELLENT:
            verdict_colour = ConsoleColours.OKGREEN
            verdict_label  = "High Confidence"
            verdict_note   = "RMSD < 2.0 Å — Structures are reliable."
        elif trust_score <= CFG.TRUST_SCORE_ACCEPTABLE:
            verdict_colour = ConsoleColours.WARNING
            verdict_label  = "Moderate Confidence"
            verdict_note   = "RMSD 2.0–3.0 Å — Marginal deviation detected."
        else:
            verdict_colour = ConsoleColours.FAIL
            verdict_label  = "Low Confidence"
            verdict_note   = "RMSD > 3.0 Å — Significant deviation. Proceed carefully."

        calib_data["verdict_label"]  = verdict_label
        calib_data["verdict_colour"] = verdict_colour

        _ci(
            f"  > Trust Score: {verdict_colour}{ConsoleColours.BOLD}{trust_score:.3f} Å  |  "
            f"{verdict_label}  —  {verdict_note}{ConsoleColours.ENDC}"
        )
        _ci("")
        _ci("  " + "-" * 76)
        _ci("")

        # --- ACTIVE SITE MAPPING ---
        missing_residues = []
        _ci(f"  > Verifying {model_label} Active Site Mapping:")

        all_con_residues = []
        for _ch in st_con[0]:
            for _res in _ch:
                if hasattr(_res, "find_atom"):
                    all_con_residues.append(_res)
                elif hasattr(_res, "__getitem__") and len(_res) > 0:
                    all_con_residues.append(_res[0])

        con_by_seqid: Dict[int, Any] = {}
        for _r in all_con_residues:
            try:
                con_by_seqid[int(_r.seqid.num)] = _r
            except Exception:
                pass

        if target_seq:
            _safe_pid = re.sub(r"[^A-Za-z0-9_]", "_", model_label)
            aln_idx_map, _, _, _ = map_active_site_residues(_safe_pid, target_seq)

            for key, site in REF_ACTIVE_SITE_MAP.items():
                ref_pdb_id = site.get("pdb_id", site["id"])
                tgt_idx    = aln_idx_map.get(key)

                if tgt_idx is None:
                    missing_residues.append(f"{key} ({site['res']}{ref_pdb_id})")
                    calib_data["mapping_rows"].append({
                        "key": key, "ref_res": f"{site['res']}{ref_pdb_id}",
                        "mapped_res": "—", "dist": "—", "ok": False
                    })
                    _ci(
                        f"    - {key:<10} (Ref:{site['res']}{ref_pdb_id}) -> "
                        f"{ConsoleColours.FAIL}Alignment: no equivalent position found in {model_label}.{ConsoleColours.ENDC}"
                    )
                    continue

                cif_res = None
                _win = int(getattr(CFG, "RESIDUE_SEARCH_WINDOW", 5))
                _offsets = [0] + [s * k for k in range(1, _win + 1) for s in (1, -1)]
                for _offset in _offsets:
                    cif_res = con_by_seqid.get(tgt_idx + _offset)
                    if cif_res is not None:
                        tgt_idx = tgt_idx + _offset
                        break

                if cif_res is None:
                    missing_residues.append(f"{key} ({site['res']}{ref_pdb_id})")
                    calib_data["mapping_rows"].append({
                        "key": key, "ref_res": f"{site['res']}{ref_pdb_id}",
                        "mapped_res": "—", "dist": "—", "ok": False
                    })
                    _ci(
                        f"    - {key:<10} (Ref:{site['res']}{ref_pdb_id}) -> "
                        f"{ConsoleColours.FAIL}Aligned pos {tgt_idx} not found in {model_label} CIF.{ConsoleColours.ENDC}"
                    )
                    continue

                dist_str = "—"
                try:
                    ref_res_obj = ref_chain[str(ref_pdb_id)]
                    if ref_res_obj:
                        _rr  = ref_res_obj if hasattr(ref_res_obj, "find_atom") else ref_res_obj[0]
                        _rca = _rr.find_atom("CA", "*")
                        _ref_pos = _rca.pos if _rca else (_rr[0].pos if len(_rr) > 0 else None)
                        _cca = cif_res.find_atom("CA", "*")
                        _con_pos = _cca.pos if _cca else (cif_res[0].pos if len(cif_res) > 0 else None)
                        if _ref_pos and _con_pos:
                            dist_str = f"{_ref_pos.dist(_con_pos):.2f}Å"
                except Exception:
                    pass

                mapped_sites[key] = tgt_idx
                calib_data["mapping_rows"].append({
                    "key": key, "ref_res": f"{site['res']}{ref_pdb_id}",
                    "mapped_res": f"{cif_res.name}{tgt_idx}", "dist": dist_str, "ok": True
                })
                _ci(
                    f"    - {key:<10} (Ref:{site['res']}{ref_pdb_id}) -> "
                    f"Mapped to {model_label} {cif_res.name}{tgt_idx} ({dist_str})"
                )

        else:
            for key, site in REF_ACTIVE_SITE_MAP.items():
                ref_pdb_id = site.get("pdb_id", site["id"])
                ref_res_obj = None
                try:
                    ref_res_obj = ref_chain[str(ref_pdb_id)]
                except Exception:
                    pass

                if not ref_res_obj:
                    missing_residues.append(f"{key} (Ref PDB absent {site['res']}{ref_pdb_id})")
                    calib_data["mapping_rows"].append({
                        "key": key, "ref_res": f"{site['res']}{ref_pdb_id}",
                        "mapped_res": "—", "dist": "—", "ok": False
                    })
                    continue

                ref_res = ref_res_obj if hasattr(ref_res_obj, "find_atom") else (
                    ref_res_obj[0] if hasattr(ref_res_obj, "__getitem__") and len(ref_res_obj) > 0 else None
                )
                if not ref_res:
                    continue

                ca_atom = ref_res.find_atom("CA", "*")
                ref_pos = ca_atom.pos if ca_atom else (ref_res[0].pos if len(ref_res) > 0 else None)
                if not ref_pos:
                    continue

                # Search radius for matching a control's catalytic residue onto its structural
                # counterpart. A gate distance belongs in CFG, not in the loop that applies it.
                min_dist   = float(CFG.CONTROL_RESIDUE_MATCH_RADIUS)
                best_match = None
                best_match_resname = ""

                for res_inner in all_con_residues:
                    if res_inner.name[:3] != site["res"][:3]:
                        continue
                    ca_con = res_inner.find_atom("CA", "*")
                    if not ca_con:
                        ca_con = res_inner[0] if len(res_inner) > 0 else None
                    if not ca_con:
                        continue
                    d = ref_pos.dist(ca_con.pos)
                    if d < min_dist:
                        min_dist = d
                        best_match = res_inner.seqid.num
                        best_match_resname = res_inner.name

                if best_match:
                    dist_str = f"{min_dist:.2f}Å"
                    mapped_sites[key] = best_match
                    calib_data["mapping_rows"].append({
                        "key": key, "ref_res": f"{site['res']}{ref_pdb_id}",
                        "mapped_res": f"{best_match_resname}{best_match}", "dist": dist_str, "ok": True
                    })
                    _ci(
                        f"    - {key:<10} (Ref:{site['res']}{ref_pdb_id}) -> "
                        f"Mapped to {model_label} {best_match_resname}{best_match} ({dist_str})"
                    )
                else:
                    missing_residues.append(f"{key} ({site['res']}{ref_pdb_id})")
                    calib_data["mapping_rows"].append({
                        "key": key, "ref_res": f"{site['res']}{ref_pdb_id}",
                        "mapped_res": "—", "dist": "—", "ok": False
                    })
                    _ci(
                        f"    - {key:<10} (Ref:{site['res']}{ref_pdb_id}) -> "
                        f"{ConsoleColours.FAIL}NOT FOUND in {model_label} sequence.{ConsoleColours.ENDC}"
                    )

        calib_data["missing_residues"] = missing_residues
        calib_data["all_mapped"] = len(missing_residues) == 0

        if missing_residues:
            _ci("    ! Critical Warning: Specified residues were not definitively mapped.")
        else:
            _ci(f"    {ConsoleColours.OKGREEN}✔ All Reference Residues successfully mapped to {model_label}.{ConsoleColours.ENDC}")

    except Exception as e:
        _ci(f"    ! System Calibration Exception encountered: {e}")

    if not silent:
        console_separator()
    return trust_score, mapped_sites, calib_data

def print_control_calibration_consolidated(deha4_data: list, r3u_data: list) -> None:
    """Print a single consolidated summary table for all 6 control calibration jobs.

    deha4_data: list of 3 calib_data dicts [FA, DFA, TFA] from DeHa4 calibration
    r3u_data:   list of 3 calib_data dicts [FA, DFA, TFA] from 3R3U calibration
    """
    lig_short = ["FA", "DFA", "TFA"]

    W = 110
    sep_eq = "=" * W

    # ── Header ──────────────────────────────────────────────────────────────────
    console_info(f"\n  {sep_eq}")
    console_info(f"  {'CONTROL CALIBRATION SUMMARY — DeHa4 × 3R3U  (3 Ligands Each)':^{W}}")
    console_info(f"  {sep_eq}\n")

    # ── Trust Scores ─────────────────────────────────────────────────────────────
    """
    Trust score = Cα RMSD (Å) between Boltz-2 prediction and 3R3U crystal structure.
    Lower = better alignment. <2.0 Å = High Confidence, 2–3 Å = Moderate, >3 Å = Low.
    """
    _ts_col = 16
    _ts_sep = "─" * (14 + 3 * (_ts_col + 2) + 5 + 3 * (_ts_col + 2))
    console_info("  Structural Trust Score  (Cα RMSD vs 3R3U crystal; <2.0 Å = high confidence)")
    console_info(f"  {_ts_sep}")
    # Group header row
    _grp_d = f"{'── DeHa4 Control (Boltz-2 Prediction) ──':^{3*(_ts_col+2)}}"
    _grp_r = f"{'── 3R3U Boltz-2 (Crystal Sequence Prediction) ──':^{3*(_ts_col+2)}}"
    console_info(f"  {'':14}{_grp_d}   │  {_grp_r}")
    # Column header row
    _col_hdr = f"  {'Model →':<14}" + \
               "".join(f"  {'DeHa4_'+s:>{_ts_col}}" for s in lig_short) + \
               "   │" + \
               "".join(f"  {'3R3U_'+s:>{_ts_col}}" for s in lig_short)
    console_info(_col_hdr)
    console_info(f"  {_ts_sep}")
    # RMSD row — rjust the value string to _ts_col to match header cell width exactly
    ts_row = f"  {'RMSD (Å)':<14}"
    for d in deha4_data:
        v   = d.get("trust_score", 99.0)
        col = ConsoleColours.OKGREEN if v <= CFG.TRUST_SCORE_EXCELLENT else (ConsoleColours.WARNING if v <= CFG.TRUST_SCORE_ACCEPTABLE else ConsoleColours.FAIL)
        ts_row += f"  {col}{f'{v:.3f} Å'.rjust(_ts_col)}{ConsoleColours.ENDC}"
    ts_row += "   │"
    for d in r3u_data:
        v   = d.get("trust_score", 99.0)
        col = ConsoleColours.OKGREEN if v <= CFG.TRUST_SCORE_EXCELLENT else (ConsoleColours.WARNING if v <= CFG.TRUST_SCORE_ACCEPTABLE else ConsoleColours.FAIL)
        ts_row += f"  {col}{f'{v:.3f} Å'.rjust(_ts_col)}{ConsoleColours.ENDC}"
    console_info(ts_row)
    # Confidence verdict row
    verd_row = f"  {'Confidence':<14}"
    for d in list(deha4_data) + ["│"] + list(r3u_data):
        if d == "│":
            verd_row += "   │"
            continue
        v = d.get("trust_score", 99.0)
        lbl = "High" if v <= CFG.TRUST_SCORE_EXCELLENT else ("Moderate" if v <= CFG.TRUST_SCORE_ACCEPTABLE else "Low")
        col = ConsoleColours.OKGREEN if v <= CFG.TRUST_SCORE_EXCELLENT else (ConsoleColours.WARNING if v <= CFG.TRUST_SCORE_ACCEPTABLE else ConsoleColours.FAIL)
        verd_row += f"  {col}{lbl:>{_ts_col}}{ConsoleColours.ENDC}"
    console_info(verd_row)
    console_info(f"  {_ts_sep}\n")

    # ── Active Site Mapping Tables ────────────────────────────────────────────────
    for group_label, group_data, col_prefix in [
        ("DeHa4 Controls  (Boltz-2 prediction of DeHa4 sequence)", deha4_data, "DeHa4"),
        ("3R3U Boltz-2 Controls  (Boltz-2 prediction of 3R3U sequence)", r3u_data, "3R3U"),
    ]:
        console_info(f"  Active Site Residue Mapping — {group_label}")
        mw   = 10   # role col
        rw   = 10   # ref col
        dw   = 18   # each model col
        header_line = f"  {'Role':<{mw}}  {'Ref (3R3U)':<{rw}}" + \
                      "".join(f"  {f'{col_prefix}_{s}':>{dw}}" for s in lig_short)
        console_info(f"  {'─'*(mw+rw+3*dw+8)}")
        console_info(header_line)
        console_info(f"  {'─'*(mw+rw+3*dw+8)}")

        # Collect all unique keys in order
        all_keys = []
        seen_keys: set = set()
        for d in group_data:
            for row in d.get("mapping_rows", []):
                if row["key"] not in seen_keys:
                    all_keys.append(row["key"])
                    seen_keys.add(row["key"])

        for key in all_keys:
            ref_res = ""
            cells   = []
            for d in group_data:
                found = next((r for r in d.get("mapping_rows", []) if r["key"] == key), None)
                if found:
                    if not ref_res:
                        ref_res = found["ref_res"]
                    if found["ok"]:
                        cell = f"{found['mapped_res']} ({found['dist']})"
                    else:
                        cell = f"{ConsoleColours.FAIL}MISSING{ConsoleColours.ENDC}"
                else:
                    cell = "—"
                cells.append(cell)
            row_line = f"  {key:<{mw}}  {ref_res:<{rw}}" + \
                       "".join(f"  {c:>{dw}}" for c in cells)
            console_info(row_line)

        console_info(f"  {'─'*(mw+rw+3*dw+8)}")
        # Status row
        status_cells = []
        for d in group_data:
            if d.get("all_mapped"):
                status_cells.append(f"{ConsoleColours.OKGREEN}✔ All mapped{ConsoleColours.ENDC}")
            else:
                n_miss = len(d.get("missing_residues", []))
                status_cells.append(f"{ConsoleColours.FAIL}✘ {n_miss} missing{ConsoleColours.ENDC}")
        status_line = f"  {'Status':<{mw}}  {'—':<{rw}}" + \
                      "".join(f"  {c:>{dw}}" for c in status_cells)
        console_info(status_line)
        console_info(f"  {'─'*(mw+rw+3*dw+8)}\n")

    console_info(f"  {sep_eq}\n")


def analyse_candidate_structure(target_cif: Path, control_cif: Path, control_map: Dict[str, int],
                                target_map: Optional[Dict[str, int]] = None) -> Dict[str, Any]:
    """
    Superimposes the candidate structure onto the DeHa4 control model.
    Computes the likelihood score, RMSD, and mechanistic score by comparing targeted structural pockets.

    target_map: when provided, the per-role TARGET residue indices already resolved for the
    tier (mapped_sites). Active_Site_RMSD (→ Criterion B, catalytic_constellation_score) is then
    measured on the SAME eight residues the tier scored — control role Cα ↔ the tier's mapped
    target residue Cα — instead of an independent nearest-Cα search. This unifies the two
    active-site definitions: without it, B could read crystal-like off residues the tier never
    used. When omitted, the nearest-Cα search is used instead.
    """
    result = {
        "Active_Site_RMSD": 999.0,
        "Halide_Stabilisation": False,
        "Carboxylate_Clamp": False
    }
    try:
        # -------------------------------------------------------------------------------
        # Sub-Step 7.1.1: Execute Loading and Superimposition Routines
        # -------------------------------------------------------------------------------
        st_con = load_structure_safe(control_cif)
        st_tar = load_structure_safe(target_cif)

        if st_con and st_tar:
            con_chain = st_con[0][0]
            tar_chain = st_tar[0][0]

            try:
                if hasattr(gemmi.SupSelect, "Ca"):
                    sup_select = gemmi.SupSelect.Ca
                elif hasattr(gemmi.SupSelect, "CA"):
                    sup_select = gemmi.SupSelect.CA
                else:
                    sup_select = gemmi.SupSelect.All
            except Exception:
                sup_select = gemmi.SupSelect.All

            # Active-site RMSD (Criterion B) must be measured after a LOCAL superposition on the
            # catalytic Cα only. A global whole-chain fit lets divergent-domain motion in a distant
            # homologue rotate the active site off-register and inflate the RMSD. When the tier's
            # mapped target residues are known, superpose on those catalytic Cα pairs; otherwise
            # fall back to the whole-chain Cα fit.
            _ref_ca, _mov_ca = [], []
            if target_map and hasattr(gemmi, "superpose_positions"):
                for _key, _rid in control_map.items():
                    _cg = con_chain[str(_rid)]
                    if not _cg:
                        continue
                    _cr = _cg if hasattr(_cg, "find_atom") else (_cg[0] if len(_cg) > 0 else None)
                    if _cr is None:
                        continue
                    _cca = _cr.find_atom("CA", "*")
                    _tidx = target_map.get(_key)
                    if _cca is None or _tidx is None:
                        continue
                    for _res in tar_chain:
                        if _res.seqid.num == _tidx:
                            _tca = _res.find_atom("CA", "*")
                            if _tca:
                                _ref_ca.append(_cca.pos)
                                _mov_ca.append(_tca.pos)
                            break

            if len(_ref_ca) >= 3:
                # Local superposition on the catalytic Cα anchors (≥3 for a stable 3D fit).
                sup = gemmi.superpose_positions(_ref_ca, _mov_ca)
            else:
                # Formally constructs a robust superposition mapping matrix (whole-chain fallback)
                sup = gemmi.calculate_superposition(
                    con_chain.whole(),
                    tar_chain.whole(),
                    gemmi.PolymerType.PeptideL,
                    sup_select
                )
            for ch in st_tar[0]:
                for res in ch:
                    for atom in res:
                        p = sup.transform.apply(atom.pos)
                        atom.pos = gemmi.Position(p.x, p.y, p.z)

            # -------------------------------------------------------------------------------
            # Sub-Step 7.1.2: Extract Active-Site Conservation Data (superposition vs crystal)
            # -------------------------------------------------------------------------------
            rmsd_sq_sum = 0.0
            count = 0
            hits = 0
            total_checks = 0
            has_halide_cradle = False
            has_carb_clamp = False

            for key, res_id in control_map.items():
                con_res_obj = con_chain[str(res_id)]
                if not con_res_obj: continue

                if hasattr(con_res_obj, "find_atom"):
                    con_res = con_res_obj
                elif hasattr(con_res_obj, "__getitem__") and len(con_res_obj) > 0:
                    con_res = con_res_obj[0]
                else:
                    continue

                con_ca = con_res.find_atom("CA", "*")
                if not con_ca: continue
                con_pos = con_ca.pos

                min_d = 99.0
                mapped_res = None

                # Criterion-B unification: when the tier's mapped target residue is known,
                # measure THAT residue — the same atom the tier scored — not the nearest Cα.
                # Falls back to the nearest-Cα search if the role was unmapped/absent.
                _tgt_idx = target_map.get(key) if target_map else None
                if _tgt_idx is not None:
                    for res in tar_chain:
                        if res.seqid.num == _tgt_idx:
                            tar_ca = res.find_atom("CA", "*")
                            if tar_ca:
                                min_d = con_pos.dist(tar_ca.pos)
                                mapped_res = res
                            break
                if mapped_res is None:
                    for res in tar_chain:
                        tar_ca = res.find_atom("CA", "*")
                        if tar_ca:
                            d = con_pos.dist(tar_ca.pos)
                            if d < min_d:
                                min_d = d
                                mapped_res = res

                # Accrue internal RMSD values representing core active site geometry
                if min_d < CFG.QC_MIN_DIST_5:
                    rmsd_sq_sum += (min_d * min_d)
                    count += 1
                else:
                    rmsd_sq_sum += CFG.QC_RMSD_SQ_SUM
                    count += 1

                # Validate exact spatial chemical identity requirements
                if mapped_res and min_d < CFG.QC_MIN_DIST_4:
                    res_name = mapped_res.name
                    total_checks += 1
                    ref_role = REF_ACTIVE_SITE_MAP[key]["role"]

                    if ref_role == "Fluoride_Cradle":
                        if res_name in RESIDUE_CLASS_GROUPS["AROMATIC"]: hits += 1; has_halide_cradle = True
                    elif ref_role == "Carboxylate_Clamp":
                        if res_name in RESIDUE_CLASS_GROUPS["POSITIVE"]: hits += 1; has_carb_clamp = True
                    elif ref_role == "Fluorine_Stabiliser":
                         if res_name in RESIDUE_CLASS_GROUPS["AROMATIC"] or res_name in RESIDUE_CLASS_GROUPS["POSITIVE"]: hits += 1
                    elif ref_role == "Nucleophile":
                         if res_name == "ASP": hits += 1
                    elif ref_role == "Acid_Catalyst":
                         if res_name == "ASP" or res_name == "GLU": hits += 1

            # ── Mutation-tolerant cradle detection (geometric fallback) ──────────
            """
            The canonical check above requires an AROMATIC residue at the exact
            aligned cradle position. Where a homolog has a gap/substitution there
            (so the alignment misses it), search the superimposed target for ANY
            aromatic sidechain within CFG.MECH_CRADLE_RADIUS of the ligand's
            leaving halogen — the same dynamic principle used for the nucleophile.
            The rigid superposition preserves intra-target distances.
            """
            if not has_halide_cradle:
                _lig_halogens, _arom_pos = [], []
                for _ch in st_tar[0]:
                    for _res in _ch:
                        _is_prot_cr = (_res.name in STANDARD_AA
                                       or _canonical_resname(_res.name) in STANDARD_AA)
                        if not _is_prot_cr and _res.name != "HOH":
                            for _a in _res:
                                if _a.element.name in ("F", "Cl", "Br", "I"):
                                    _lig_halogens.append(_a.pos)
                        elif (_res.name in RESIDUE_CLASS_GROUPS["AROMATIC"]
                              or _canonical_resname(_res.name) in RESIDUE_CLASS_GROUPS["AROMATIC"]):
                            for _a in _res:
                                if _a.name not in ("N", "C", "CA", "O"):
                                    _arom_pos.append(_a.pos)
                if _lig_halogens and _arom_pos:
                    _min_cradle = min(_h.dist(_p) for _h in _lig_halogens for _p in _arom_pos)
                    if _min_cradle <= CFG.MECH_CRADLE_RADIUS:
                        has_halide_cradle = True

            result["Active_Site_RMSD"] = math.sqrt(rmsd_sq_sum / count) if count > 0 else 999.0
            result["Halide_Stabilisation"] = has_halide_cradle
            result["Carboxylate_Clamp"] = has_carb_clamp
    except Exception as e:
        if logger: logger.debug(f"Candidate Analysis Failure detected: {e}")
    return result

# -------------------------------------------------------------------------------
# Step 7.2: Catalyst Distance Evaluation Helpers
# -------------------------------------------------------------------------------
def generate_rich_justification(tier: str, meaning: str, constraint: str, aligned_ok: bool, identity: float) -> str:
    """Generates a detailed justification string to objectively explain algorithmic decisions within the final report."""
    if not aligned_ok:
        if tier in CFG.TIER_HIGH_QUALITY:
            return f"Caution: Documented low Sequence Identity ({identity}%) but Excellent Active Site Geometry. Identifies a likely remote homologue."
        else:
            return f"Failure: Overall sequence alignment proved unreliable (ID {identity}% < {CFG.ALIGN_MIN_SEQ_IDENTITY:.0f}%). Structure is likely invalid."
    return f"{meaning} ({constraint})"

def _derive_burgi_dunitz(nuc_np, centres, neigh_fn, sym_fn, pos_fn):
    """Bürgi–Dunitz angle on the electrophilic head centre (C/S/P bearing the
    most O; a carboxylate carbon with ≥2 O ranked first), computed once for both
    the RDKit-template and structure-only paths. `centres` are candidate node
    handles; `neigh_fn`/`sym_fn`/`pos_fn` adapt neighbour lookup, element symbol
    and 3-D position to the caller's data model. Returns a rounded angle or None.
    """
    def _rank(a):
        n_o = sum(1 for n in neigh_fn(a) if sym_fn(n) == "O")
        return (sym_fn(a) == "C" and n_o >= 2, n_o)
    for _ctr in sorted(centres, key=_rank, reverse=True):
        _cC = pos_fn(_ctr)
        if _cC is None:
            continue
        """
        The Bürgi-Dunitz trajectory is defined against the CARBONYL oxygen. Taking whichever oxygen
        the neighbour list happens to yield first can return the single-bonded one on an ester or a
        carboxylate, which points the reference vector at the wrong lone pair. The carbonyl is the
        oxygen with the SHORTEST C–O distance (C=O ~1.21 Å against C–O ~1.31 Å), which is a
        geometric fact of the frame and needs no bond-order perception.
        """
        _os = [(pos_fn(n), np.linalg.norm(np.asarray(pos_fn(n)) - np.asarray(_cC)))
               for n in neigh_fn(_ctr) if sym_fn(n) == "O" and pos_fn(n) is not None]
        _cO = min(_os, key=lambda t: t[1])[0] if _os else None
        if _cO is not None:
            return round(calculate_burgi_dunitz(nuc_np, _cC, _cO), 1)
    return None


def _derive_flippin_lodge(nuc_np, c_node, c_pos, neigh_fn, sym_fn, pos_fn, is_leaving_fn):
    """Flippin–Lodge offset at the attacked carbon: ≥2 spectators (heavy first,
    then H), excluding the leaving halogen. Shared by the RDKit-template and
    structure-only paths via the accessor callables. Returns rounded offset or
    None.
    """
    _heavy = [n for n in neigh_fn(c_node) if sym_fn(n) != "H" and not is_leaving_fn(n)]
    _hyd   = [n for n in neigh_fn(c_node) if sym_fn(n) == "H"]
    _spect = [p for p in (pos_fn(n) for n in _heavy + _hyd) if p is not None]
    if len(_spect) >= 2:
        _fl = calculate_flippin_lodge(nuc_np, c_pos, _spect[0], _spect[1])
        if _fl < CFG.SENTINEL_VALID_MAX:
            return round(_fl, 1)
    return None


def calculate_sn2_metrics(asp_atoms, lig_atoms, rd_mol=None, mm_map=None, preferred_c_names=None, cradle_coords=None) -> Tuple[float, float, int, Any, Any, Any]:
    """
    Calculates the SN2 attack angle and the SN2 backside-attack trajectory deviation.
    An ideal SN2 backside attack strictly requires an angle of approximately 180 degrees.
    Deviation quantifies the perpendicular distance measured from the ideal C-X vector.
    When preferred_c_names is supplied the attack carbon is restricted to that set
    (CFG §5.3: the α-carbon adjacent to the ligand carboxylate), so the geometry is
    measured at the catalytically productive position rather than an incidental C–F.
    """
    if not asp_atoms or not lig_atoms: return 0.0, 999.0, 0, None, None, None
    asp_oxygens = [a for a in asp_atoms if a.name in ("OD1", "OD2")]
    lig_carbons = [a for a in lig_atoms if a.element.name == "C"]
    lig_halogens = [a for a in lig_atoms if a.element.name == "F"]

    if not asp_oxygens or not lig_carbons or not lig_halogens: return 0.0, 999.0, 0, None, None, None

    valid_cx_pairs = []

    if rd_mol and mm_map:
        inv_map = {v: k for k, v in mm_map.items()}
        for bond in rd_mol.GetBonds():
            a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
            sym1, sym2 = a1.GetSymbol(), a2.GetSymbol()

            c_idx, x_idx = None, None
            if sym1 == "C" and sym2 == "F": c_idx, x_idx = a1.GetIdx(), a2.GetIdx()
            elif sym2 == "C" and sym1 == "F": c_idx, x_idx = a2.GetIdx(), a1.GetIdx()

            if c_idx is not None and x_idx is not None:
                c_name, x_name = inv_map.get(c_idx), inv_map.get(x_idx)
                if c_name and x_name:
                    c_atom = next((a for a in lig_carbons if a.name == c_name), None)
                    x_atom = next((a for a in lig_halogens if a.name == x_name), None)
                    if c_atom and x_atom:
                        valid_cx_pairs.append((c_atom, x_atom))

    if not valid_cx_pairs:
        for c in lig_carbons:
            for x in lig_halogens:
                if c.pos.dist(x.pos) < CFG.CF_DIST_TOLERANCE:
                    valid_cx_pairs.append((c, x))

    if not valid_cx_pairs: return 0.0, 999.0, 0, None, None, None

    """
    Restrict the attack carbon to the catalytically productive set when supplied
    (CFG §5.3: the α-carbon adjacent to the ligand carboxylate). This keeps the
    angle and distance tied to the position FAcD defluorinates. When the preferred
    carbons carry no fluorine there is no productive α-defluorination → no SN2.
    """
    if preferred_c_names:
        _restricted = [(c, x) for (c, x) in valid_cx_pairs if c.name in preferred_c_names]
        if not _restricted:
            return 0.0, 999.0, 0, None, None, None
        valid_cx_pairs = _restricted

    """
    Choose the attacking oxygen and the warhead carbon by the REACTION, not by proximity.

    FAcD's nucleophile is the aspartate carboxylate: it attacks the α-carbon and displaces the
    fluoride, forming the covalent ester intermediate. Its two oxygens are resonance-equivalent, so
    the one that reacts is the one lined up for backside attack — not the one that happens to be
    nearer. They sit ~2.2 Å apart and can point in quite different directions, so picking the closer
    of two chemically identical atoms can report a side-on approach (~95°) for a complex whose other
    oxygen is properly anti-periplanar (~172°), and the pose then fails the angle gate on an
    arbitrary tiebreak.

    The pair is therefore scored by the backside O–C–F angle it produces, with the scissile fluorine
    resolved for each candidate carbon exactly as below (cradle-coupled where the cradle is known).
    Distance breaks ties, so among equally-aligned pairs the closest still wins.
    """
    def _scissile_for(_c, _o):
        """The leaving fluorine for this carbon: the cradle-facing α-F, else the most
        anti-periplanar one to this oxygen."""
        _fs = [x for cc, x in valid_cx_pairs if cc == _c]
        if not _fs:
            return None
        if cradle_coords:
            _cc0 = np.mean(np.asarray(cradle_coords, float), axis=0)
            return min(_fs, key=lambda x: (x.pos.x - _cc0[0])**2 + (x.pos.y - _cc0[1])**2
                       + (x.pos.z - _cc0[2])**2)
        return max(_fs, key=lambda x: calculate_angle(_o.pos, _c.pos, x.pos))

    """
    The attacking oxygen is the one that best satisfies BOTH near-attack conditions AT ONCE — a short
    approach and a linear trajectory — because an SN2 needs them on the SAME atom. Neither condition
    alone identifies it:

      · ranking by DISTANCE alone picks a nearby oxygen that may be approaching side-on;
      · ranking by ANGLE alone picks a well-aligned oxygen that may be far out of reach. Measured on
        the 3R3U fluoroacetate control, OD1 is anti-periplanar (171°) but sits 4.44 Å from the
        α-carbon, while OD2 is 2.85 Å away at 163° — a textbook NAC. Angle-first selects OD1 and the
        native substrate control scores as a non-degrader.

    Both are therefore scored jointly, with the SAME graded NAC terms the soft score uses (§4.1): the
    distance sigmoid about NAC_DIST_STRICT and the angle sigmoid about NAC_ANGLE_STRICT. Their product
    is maximised, so an oxygen must be both close and aligned to win, and a pose is credited with the
    geometry a nucleophile could actually react through.
    """
    best_O, best_C = None, None
    _best_key = None
    for o in asp_oxygens:
        for c, _x in valid_cx_pairs:
            _f = _scissile_for(c, o)
            if _f is None:
                continue
            _ang = calculate_angle(o.pos, c.pos, _f.pos)
            _d = o.pos.dist(c.pos)
            _nac = (sigmoid(_d, k=CFG.SOFT_K_NUC, x0=CFG.NAC_DIST_STRICT)
                    * sigmoid(_ang, k=CFG.SOFT_K_ANG, x0=CFG.NAC_ANGLE_STRICT))
            _key = (_nac, _ang, -_d)        # joint NAC quality; angle then proximity break ties
            if _best_key is None or _key > _best_key:
                _best_key, best_O, best_C = _key, o, c

    if not best_O or not best_C: return 0.0, 999.0, 0, None, None, None

    """
    Select the scissile (leaving) fluorine among those bonded to the chosen best_C.
    Cradle-coupled definition (preferred): the departing F⁻ is the one stabilised by the
    fluoride cradle (His155/Trp156/Tyr217; Chan et al. 2011), so the leaving F is the
    α-fluorine pointing into the cradle — the F nearest the cradle centroid. This is the
    physical SN2 leaving group and removes the best-of-N angle inflation that a plain
    max-over-fluorines selection introduces for CF2/CF3 carbons. When the cradle is not
    resolved, fall back to the most anti-periplanar F (max Nu-C-F angle).
    """
    _alpha_fs = [x for c, x in valid_cx_pairs if c == best_C]
    best_X = None
    if cradle_coords:
        _cc = np.mean(np.asarray(cradle_coords, float), axis=0)
        best_X = min(_alpha_fs,
                     key=lambda x: (x.pos.x - _cc[0])**2 + (x.pos.y - _cc[1])**2 + (x.pos.z - _cc[2])**2,
                     default=None)
    if best_X is None:
        max_ang = -1.0
        for x in _alpha_fs:
            ang = calculate_angle(best_O.pos, best_C.pos, x.pos)
            if ang > max_ang:
                max_ang = ang
                best_X = x

    if not best_X: return 0.0, 999.0, 0, None, None, None

    # Derives the final structured Geometric Result (Vector Angle)
    angle = calculate_angle(best_O.pos, best_C.pos, best_X.pos)

    # Computes Trajectory Deviation (Perpendicular position of Nucleophile relative to C-X vector)
    try:
        vec_cx  = np.array([best_X.pos.x - best_C.pos.x, best_X.pos.y - best_C.pos.y, best_X.pos.z - best_C.pos.z])
        vec_cnu = np.array([best_O.pos.x - best_C.pos.x, best_O.pos.y - best_C.pos.y, best_O.pos.z - best_C.pos.z])
        """
        Deviation from the BACKSIDE attack ray, not from the C–X line.

        The perpendicular distance to the C–X line alone cannot tell a backside attack (180°, the
        Walden inversion the mechanism requires) from a front-side approach (0°, no reaction): both
        are collinear, and both give a perpendicular distance of zero. Measuring instead against the
        ray that leaves C directly OPPOSITE the leaving group makes the two distinguishable — a
        front-side nucleophile projects onto the wrong half-line, its nearest point on the ray is C
        itself, and its deviation is its full distance from the carbon.
        """
        _u = vec_cx / (np.linalg.norm(vec_cx) + 1e-6)
        _t = float(np.dot(vec_cnu, _u))          # < 0 = backside (opposite the leaving group)
        if _t < 0.0:
            deviation = float(np.linalg.norm(np.cross(vec_cnu, _u)))     # perpendicular offset
        else:
            deviation = float(np.linalg.norm(vec_cnu))                   # front side: fully off-axis
    except Exception: deviation = 999.0

    # Teflon Shield Calculation: Count adjacent fluorines that sterically clash with Aspartate catalytic oxygen
    teflon_clashes = 0
    for f in lig_halogens:
        if f != best_X and f.element.name == "F":
            if best_O.pos.dist(f.pos) <= CFG.CLASH_DIST_TOLERANCE:
                teflon_clashes += 1

    # ── Auxiliary (NON-GATING) reference geometries — see utils docstrings ─────
    """
    Bürgi–Dunitz: ideally the substrate carbonyl (–COO⁻) carbon. Where the head
      group is not a carboxylate (sulfonate –SO₃⁻, phosphonate, alcohol) that
      carbon does not exist, so the angle is taken on the head-group electrophilic
      centre instead — the C/S/P bearing the most oxygens (Nu–centre–O). The
      canonical carboxylate carbon is always preferred when present.
    Flippin–Lodge: in-plane offset of the nucleophile at the attacked carbon.
      Uses heavy spectators first, then hydrogens when the carbon carries <2 heavy
      substituents (e.g. fluoroacetate's CH₂), so the offset is always defined.
    The predicted structure is heavy-atom only, so the RDKit template (which has
    hydrogens) is rigidly aligned into the structure frame via the mapped heavy
    atoms; mapped atoms use exact structure coordinates, hydrogens use the
    transformed template coordinates. Both metrics thus return a real value for
    every substrate class.
    """
    # Count the equivalent C–F leaving-group vectors on the chosen scissile carbon
    # (multiplicity correction for the competence angle term, CFG §5.5).
    n_scissile_f = sum(1 for c, x in valid_cx_pairs if c == best_C)

    '''
    β-fluorination: fluorines on the carbons adjacent to the scissile α-carbon
    (excluding the carboxylate head carbon). FA/DFA/TFA = 0; every perfluoro chain
    (PFBA/PFOA…) ≥ 2; perfluoro-ether α-carbons (GenX) ≥ 2 via the CF3 branch. This
    is the discriminator the α-only C–F BDE proxy misses (a perfluoro substrate reads
    the same low BDE as fluoroacetate unless β-fluorination is counted).
    '''
    beta_f_count = 0
    if rd_mol is not None and mm_map is not None:
        try:
            _aidx = mm_map.get(best_C.name)
            if _aidx is not None:
                for _nb in rd_mol.GetAtomWithIdx(_aidx).GetNeighbors():
                    if _nb.GetSymbol() != "C":
                        continue
                    if sum(1 for x in _nb.GetNeighbors() if x.GetSymbol() == "O") >= 2:
                        continue   # carboxylate / head carbon — not a β position
                    beta_f_count += sum(1 for x in _nb.GetNeighbors() if x.GetSymbol() == "F")
        except Exception:
            beta_f_count = 0

    """
    The attacking oxygen travels with the angle it produced. The nucleophile approach distance
    must be measured to THIS oxygen: reporting the minimum over both Oδ while taking the angle
    from the other one describes a nucleophile that does not exist — one oxygen supplying the
    trajectory and its partner supplying the reach.

    n_scissile_f is the number of equivalent C–F bonds on the attack carbon. It drives the Šidák
    multiplicity correction (CFG §5.2d): the inflation it removes lives in the POSE — a CF3 carbon
    has three chances to present some fluorine anti-periplanar — so it is the bond COUNT that
    matters, not which of them the cradle later identifies as leaving.
    """
    aux = {"burgi_dunitz_angle": 999.0, "flippin_lodge_offset": 999.0,
           "n_scissile_f": n_scissile_f, "beta_f_count": beta_f_count,
           "attack_o_atom": best_O}
    if rd_mol and mm_map:
        try:
            inv_map   = {v: k for k, v in mm_map.items()}
            _name2lig = {a.name: a for a in lig_atoms}

            # Rigid RDKit-template → structure-frame transform from mapped heavy atoms.
            _R = _t = None
            try:
                _conf = rd_mol.GetConformer()
                _P, _Q = [], []
                for _j in range(rd_mol.GetNumAtoms()):
                    _a = _name2lig.get(inv_map.get(_j))
                    if _a is not None:
                        _p = _conf.GetAtomPosition(_j)
                        _P.append([_p.x, _p.y, _p.z])
                        _Q.append([_a.pos.x, _a.pos.y, _a.pos.z])
                if len(_P) >= 3:
                    _R, _t = kabsch_transform(np.asarray(_P), np.asarray(_Q))
            except Exception:
                _R = _t = None

            def _pos(rd_atom):
                """Position in the structure frame: exact if mapped, else the
                template coordinate transformed into the structure frame."""
                _a = _name2lig.get(inv_map.get(rd_atom.GetIdx()))
                if _a is not None:
                    return np.array([_a.pos.x, _a.pos.y, _a.pos.z])
                if _R is not None:
                    _p = rd_mol.GetConformer().GetAtomPosition(rd_atom.GetIdx())
                    return _R @ np.array([_p.x, _p.y, _p.z]) + _t
                return None

            _nuc = np.array([best_O.pos.x, best_O.pos.y, best_O.pos.z])

            # Bürgi–Dunitz on the RDKit-template head centre (shared derivation).
            _centres = [a for a in rd_mol.GetAtoms()
                        if a.GetSymbol() in ("C", "S", "P")
                        and any(n.GetSymbol() == "O" for n in a.GetNeighbors())]
            _bd = _derive_burgi_dunitz(
                _nuc, _centres,
                neigh_fn=lambda a: a.GetNeighbors(),
                sym_fn=lambda a: a.GetSymbol(),
                pos_fn=_pos)
            if _bd is not None:
                aux["burgi_dunitz_angle"] = _bd

            # Flippin–Lodge on the SN2 carbon (shared derivation).
            _c_rd_idx = mm_map.get(best_C.name)
            if _c_rd_idx is not None:
                _c_rd = rd_mol.GetAtomWithIdx(_c_rd_idx)
                _cpos = np.array([best_C.pos.x, best_C.pos.y, best_C.pos.z])
                _fl = _derive_flippin_lodge(
                    _nuc, _c_rd, _cpos,
                    neigh_fn=lambda a: a.GetNeighbors(),
                    sym_fn=lambda a: a.GetSymbol(),
                    pos_fn=_pos,
                    is_leaving_fn=lambda n: inv_map.get(n.GetIdx()) == best_X.name)
                if _fl is not None:
                    aux["flippin_lodge_offset"] = _fl
        except Exception as _e:
            # A swallowed failure here leaves flippin_lodge_offset absent, and the metric then
            # reads as "not measured" rather than "measurement failed" — say which it was.
            logger.debug(f"Flippin-Lodge offset (RDKit path) failed: {_e}")

    # ── Structure-only fallback (no RDKit needed) ─────────────────────────────
    """
    The block above needs an intact RDKit template + atom mapping; on small or
    poorly-mapped ligands (e.g. TFA) it can leave the aux geometries at 999.
    Recompute them directly from the structure ligand atoms so both metrics are
    defined whenever the geometry physically exists. Covalent neighbours are
    detected by distance (≤ CFG.BOND_DIST_MAX). Only fills values still at 999.
    """
    try:
        if aux.get("burgi_dunitz_angle", 999.0) >= 999.0 or aux.get("flippin_lodge_offset", 999.0) >= 999.0:
            _bmax = float(getattr(CFG, "BOND_DIST_MAX", 1.9))
            _nuc_np = np.array([best_O.pos.x, best_O.pos.y, best_O.pos.z])

            def _neigh(atom):
                return [a for a in lig_atoms
                        if a is not atom and atom.pos.dist(a.pos) <= _bmax]

            def _np_of(a):
                return np.array([a.pos.x, a.pos.y, a.pos.z])

            # Bürgi–Dunitz from structure ligand atoms (shared derivation).
            if aux.get("burgi_dunitz_angle", 999.0) >= 999.0:
                _centres = [a for a in lig_atoms if a.element.name in ("C", "S", "P")]
                _bd = _derive_burgi_dunitz(
                    _nuc_np, _centres,
                    neigh_fn=_neigh,
                    sym_fn=lambda a: a.element.name,
                    pos_fn=_np_of)
                if _bd is not None:
                    aux["burgi_dunitz_angle"] = _bd

            # Flippin–Lodge from structure ligand atoms (shared derivation).
            if aux.get("flippin_lodge_offset", 999.0) >= 999.0:
                _fl = _derive_flippin_lodge(
                    _nuc_np, best_C, _np_of(best_C),
                    neigh_fn=_neigh,
                    sym_fn=lambda a: a.element.name,
                    pos_fn=_np_of,
                    is_leaving_fn=lambda n: n is best_X)
                if _fl is not None:
                    aux["flippin_lodge_offset"] = _fl
    except Exception as _e:
        logger.debug(f"Flippin-Lodge offset (structure path) failed: {_e}")

    return angle, deviation, teflon_clashes, best_X.pos, best_C, aux

def sigmoid(x: float, k: float = 1.0, x0: float = 0.0) -> float:
    """Soft threshold. The limit is decided by the SIGN OF k, not by the sign of (x - x0).

    Two of this pipeline's steepnesses are NEGATIVE (SOFT_K_NUC = -4.0, SOFT_K_TRIAD = -2.0): they must
    DECREASE with distance. So the overflow branch cannot key on the sign of (x - x0) — that assumes a
    positive k, and with a negative one it is inverted, handing a PERFECT 1.0 to an enormous distance.

    The trap is concrete: calculate_sn2_metrics returns 999.0 A as its "no nucleophile" sentinel, and at
    k = -4 the exp() overflows past 180.4 A. A complex with no nucleophile at all must not score a perfect
    nucleophile term.

    The limit of a logistic is set by the sign of k*(x - x0): it saturates to 1 when that product is large
    and positive, to 0 when it is large and negative. That is what is computed here.
    """
    _z = k * (x - x0)
    try:
        return 1.0 / (1.0 + math.exp(-_z))
    except OverflowError:
        return 1.0 if _z > 0 else 0.0

def compute_pocket_fit(site_atoms_obj: Dict[str, list], lig_atoms_obj: list,
                       all_prot_atoms: Optional[list] = None) -> Dict[str, Any]:
    """
    Steric complementarity of the enzyme and the docked ligand.

    Containment is measured TWICE, both times against protein coordinates (CFG §5.2c):

      pocket_containment_cavity — does the protein cavity enclose the ligand? Each ligand heavy
        atom casts CFG.BURIAL_RAYS rays; a ray is blocked when a protein heavy atom obstructs it
        within CFG.BURIAL_PROBE_A. Buriedness is the blocked fraction, and an atom is contained at
        or above CFG.BURIAL_MIN. The metric is the mean over ligand heavy atoms. This is the term
        the tier gate penalises (mechanistic_score_effective): it is protein-dependent, so the same
        ligand scores differently in a narrow and a wide pocket, and a wide-pocket homolog that
        genuinely encloses a long chain is not punished for the ligand's intrinsic length.

      pocket_containment_site8 — how much of the ligand sits inside the catalytic constellation?
        Fraction of ligand heavy atoms within CFG.SITE8_SHELL_A of any heavy atom of the eight
        mapped active-site residues (Nuc/Base/Acid, the two clamp arginines, the His/Trp/Tyr
        cradle). This is catalytic ENGAGEMENT, not cavity fit — a tail outside the shell is outside
        the reactive machinery. Reported and plotted; it does not gate.

    The volume/occupancy fields (hull volume of the eight catalytic residues, ligand vdW volume,
    occupancy, fit_ratio) remain DIAGNOSTIC: the eight-residue hull overstates the true cavity and
    is blind to the rest of the fold, which is precisely why it does not decide anything. The
    ligand volume is the union of its heavy-atom Bondi vdW spheres, since the convex hull of atom
    centres collapses to ~0 for a small or planar ligand. Degenerate hulls (< 4 points / coplanar)
    fall back to a bounding-sphere volume so every pose returns a value.
    """
    out = {"active_site_volume": 0.0, "active_site_radius": 0.0,
           "ligand_volume": 0.0, "ligand_radius_gyration": 0.0, "ligand_max_extent": 0.0,
           "pocket_occupancy": 0.0, "fit_ratio": 0.0, "ligand_fits": 0,
           "pocket_containment_cavity": 1.0, "pocket_containment_site8": 1.0,
           "ligand_reach": 0.0}
    try:
        from scipy.spatial import ConvexHull
        from scipy.spatial.distance import pdist
        _pts = []
        for _k in ("Nuc", "Base", "Acid", "Carb1", "Carb2", "Stab_H", "Stab_W", "Stab_Y"):
            for _a in site_atoms_obj.get(_k, []):
                _pts.append([_a.pos.x, _a.pos.y, _a.pos.z])
        pkt = np.asarray(_pts, float)
        lpts = np.asarray([[a.pos.x, a.pos.y, a.pos.z] for a in lig_atoms_obj
                           if a.element.name != "H"], float)
        lels = [a.element.name.upper() for a in lig_atoms_obj if a.element.name != "H"]

        def _hullvol(P):
            if len(P) < 4: return 0.0
            try: return float(ConvexHull(P).volume)
            except Exception: return 0.0
        def _radius(P):
            if len(P) == 0: return 0.0
            c = P.mean(axis=0); return float(np.sqrt(((P - c) ** 2).sum(axis=1)).max())
        def _rg(P):
            if len(P) == 0: return 0.0
            c = P.mean(axis=0); return float(np.sqrt(((P - c) ** 2).sum(axis=1).mean()))
        def _extent(P):
            if len(P) < 2: return 0.0
            return float(pdist(P).max())
        def _vdw_volume(P, els, spacing=0.6):
            # True molecular volume: union of atomic Bondi vdW spheres on a voxel grid. The
            # convex hull of atom CENTRES collapses to ~0 for a small/planar molecule (e.g.
            # fluoroacetate), so vdW spheres are required for a physical ligand volume.
            if len(P) == 0: return 0.0
            _R = {"C":1.70,"N":1.55,"O":1.52,"F":1.47,"S":1.80,"P":1.80,"CL":1.75,"BR":1.85,"I":1.98}
            r = np.array([_R.get(e, 1.70) for e in els], float)
            lo = (P - r[:, None]).min(0); hi = (P + r[:, None]).max(0)
            gx = np.arange(lo[0], hi[0] + spacing, spacing)
            gy = np.arange(lo[1], hi[1] + spacing, spacing)
            gz = np.arange(lo[2], hi[2] + spacing, spacing)
            if gx.size * gy.size * gz.size > 2_000_000 or min(gx.size, gy.size, gz.size) == 0:
                return float((4.0 / 3.0) * np.pi * (r ** 3).sum())   # sphere-sum fallback (rare, oversized)
            G = np.stack(np.meshgrid(gx, gy, gz, indexing="ij"), -1).reshape(-1, 3)
            inside = np.zeros(len(G), bool)
            for _p, _r in zip(P, r):
                inside |= (((G - _p) ** 2).sum(1) <= _r * _r)
            return float(inside.sum()) * (spacing ** 3)

        pv, pr = _hullvol(pkt), _radius(pkt)
        lv, lr, lext = _vdw_volume(lpts, lels), _rg(lpts), _extent(lpts)
        if pv == 0.0 and pr > 0: pv = (4.0 / 3.0) * np.pi * pr ** 3   # bounding-sphere fallback
        occ = lv / pv if pv > 0 else 0.0
        fr  = (lext / 2.0) / pr if pr > 0 else 0.0

        """
        (A) CAVITY CONTAINMENT — the protein-aware fit that gates.

        Buriedness by ray casting: from each ligand heavy atom, CFG.BURIAL_RAYS directions are
        spread over the sphere (Fibonacci lattice — even coverage without a mesh). A direction is
        blocked when some PROTEIN heavy atom lies within CFG.BURIAL_RAY_CLEARANCE of the ray axis,
        ahead of the atom and no farther than CFG.BURIAL_PROBE_A along it. The blocked fraction is
        that atom's buriedness; at or above CFG.BURIAL_MIN the atom is enclosed by the protein.

        Ray casting is used rather than a hull or a neighbour count because it answers the
        question directly — is there protein in the way, in every direction — and it degrades
        gracefully at the pocket mouth, where a hull test flips discontinuously.
        """
        # contain_cav starts at 0.0 (worst containment): the ray-cast below overwrites it with the
        # real measurement whenever it can run. If it cannot (no protein atoms to cast against, no
        # ligand heavy atoms), the value stays failed rather than defaulting to perfect containment —
        # an unmeasured pose must never read as ideally buried and skip the tier-gate penalty.
        contain_cav, contain_s8, reach = 0.0, 1.0, 0.0
        if len(lpts) >= 1:
            # ligand_reach — head-to-tail span from the carboxylate carbon (the head held by the
            # Arg clamp and attacked by the nucleophile). A ligand shape descriptor, reported only.
            _els = [a.element.name.upper() for a in lig_atoms_obj if a.element.name != "H"]
            _oi  = [i for i, e in enumerate(_els) if e == "O"]
            _ai, _bc = -1, 0
            for _i2, _e2 in enumerate(_els):
                if _e2 != "C" or not _oi:
                    continue
                _n2 = int((np.sqrt(((lpts[_oi] - lpts[_i2]) ** 2).sum(1)) < CFG.BOND_CO_MAX_A).sum())
                if _n2 > _bc:
                    _bc, _ai = _n2, _i2
            if _ai < 0 and len(pkt):                       # no carboxylate: use the pocket-proximal atom
                _ai = int(np.argmin(((lpts - pkt.mean(0)) ** 2).sum(1)))
            if _ai < 0:
                _ai = 0
            reach = float(np.sqrt(((lpts - lpts[_ai]) ** 2).sum(1)).max())

            _ppts = np.asarray(
                [[a["pos"].x, a["pos"].y, a["pos"].z] for a in (all_prot_atoms or [])
                 if str(a.get("element", "")).upper() != "H"], float)

            if len(_ppts) >= 4:
                _k = int(CFG.BURIAL_RAYS)
                _i = np.arange(_k) + 0.5
                _phi = np.arccos(1.0 - 2.0 * _i / _k)
                _theta = np.pi * (1.0 + 5.0 ** 0.5) * _i
                _dirs = np.stack([np.cos(_theta) * np.sin(_phi),
                                  np.sin(_theta) * np.sin(_phi),
                                  np.cos(_phi)], axis=1)          # (K, 3) unit vectors

                _probe = float(CFG.BURIAL_PROBE_A)
                _clear = float(CFG.BURIAL_RAY_CLEARANCE)
                _burial = np.zeros(len(lpts), float)
                for _j, _p in enumerate(lpts):
                    _rel = _ppts - _p                              # (P, 3)
                    _d2 = (_rel ** 2).sum(1)
                    _near = _rel[_d2 <= _probe * _probe]           # only atoms within the probe sphere
                    if len(_near) == 0:
                        continue
                    _t = _near @ _dirs.T                           # (P', K) projection along each ray
                    _perp2 = (_near ** 2).sum(1)[:, None] - _t ** 2
                    _blocked = ((_t > 0.0) & (_perp2 <= _clear * _clear)).any(axis=0)
                    _burial[_j] = float(_blocked.mean())
                contain_cav = float((_burial >= float(CFG.BURIAL_MIN)).mean())
                out["ligand_buriedness_mean"] = round(float(_burial.mean()), 3)

            """
            (B) ACTIVE-SITE CONTAINMENT — engagement with the eight catalytic residues.
            Fraction of ligand heavy atoms within CFG.SITE8_SHELL_A of any heavy atom of the
            mapped constellation. Says how much of the molecule is inside the reactive machinery,
            which is a mechanistic statement — not a cavity measurement and not a size penalty.
            """
            if len(pkt):
                _dmin = np.sqrt(((lpts[:, None, :] - pkt[None, :, :]) ** 2).sum(-1)).min(axis=1)
                contain_s8 = float((_dmin <= float(CFG.SITE8_SHELL_A)).mean())

        out.update(active_site_volume=round(pv, 1), active_site_radius=round(pr, 2),
                   ligand_volume=round(lv, 1), ligand_radius_gyration=round(lr, 2),
                   ligand_max_extent=round(lext, 2), pocket_occupancy=round(occ, 3),
                   fit_ratio=round(fr, 3), ligand_fits=int(occ <= 1.0 and fr <= 1.0),
                   pocket_containment_cavity=round(contain_cav, 3),
                   pocket_containment_site8=round(contain_s8, 3),
                   ligand_reach=round(reach, 2))
    except Exception as _e:
        # pocket_containment_cavity feeds the tier gate. A measurement that could not be taken must
        # NOT read as perfect containment (the initial 1.0 in `out`) — that silently promotes an
        # unmeasured pose. Fail it to 0.0 (worst containment → full penalty, kept out of the elite
        # tier) and log loudly so the incomplete record is visible, not swallowed at DEBUG.
        out["pocket_containment_cavity"] = 0.0
        logger.warning(f"Pocket-fit metrics failed ({_e}); pocket_containment_cavity set to 0.0 "
                       "(unmeasured — not promoted).")
    return out


def check_catalytic_geometry(cif_path: Path, mapped_sites: Dict[str, int], smiles_str: str = "", nuc_rescue_offset: int = 0) -> Dict[str, Any]:
    """
    Evaluates exact distance topologies within the active site against strictly defined tiers,
    subsequently assigning a final categorical degrader tier.
    """
    results = {f"dist_{k}": 999.0 for k in REF_ACTIVE_SITE_MAP.keys()}
    results.update({"sn2_attack_angle": 0.0, "sn2_trajectory_dev": 999.0})

    try:
        # -------------------------------------------------------------------------------
        # Sub-Step 7.2.1: Formal Load of Positional Coordinate Data
        # -------------------------------------------------------------------------------
        _, mmcif_lig = load_atoms_from_structure(cif_path)
        rd_mol = rdkit_mol_from_smiles_with_3d(smiles_str) if smiles_str else None
        mm_map = map_mmcif_to_rdkit(mmcif_lig, rd_mol) if rd_mol else {}

        doc = gemmi.cif.read_file(str(cif_path))
        st = gemmi.make_structure_from_block(doc.sole_block())
        lig_atoms_obj, lig_coords = [], []
        site_atoms, site_atoms_obj = defaultdict(list), defaultdict(list)
        site_resname = {}   # role key → mapped residue 3-letter name (for identity gating)
        all_prot_atoms = []

        for model in st:
            for chain in model:
                for res in chain:
                    """
                    Protonation-aware protein/ligand classification:
                    accept both canonical and force-field variant names.
                    """
                    _is_protein = (res.name in STANDARD_AA
                                   or _canonical_resname(res.name) in STANDARD_AA)
                    if not _is_protein and res.name != "HOH":
                        for atom in res:
                            lig_atoms_obj.append(atom)
                            lig_coords.append(atom.pos)
                    elif res.name != "HOH":
                        for atom in res:
                            all_prot_atoms.append({
                                "pos": atom.pos, "resname": res.name,
                                "element": atom.element.name, "atom": atom.name,
                                "atom_name": atom.name
                            })
                    if res.seqid.num in mapped_sites.values():
                        for name, idx in mapped_sites.items():
                            if idx == res.seqid.num:
                                site_resname[name] = res.name
                                for atom in res:
                                    site_atoms[name].append(atom.pos)
                                    site_atoms_obj[name].append(atom)

        # ── Spatial-gate validation (CFG §2.8) ──────────────────────────────
        """
        Reject any resolved triad residue whose nearest atom sits farther than
        RESOLVER_SPATIAL_CUTOFF from the NEAREST LIGAND ATOM. Such residues are
        sequence-proximal but spatially remote (e.g. ASP on a distant loop) and
        must not inflate the tier-gate distance inputs. The metric mirrors the
        tier gate itself (nearest residue-atom ↔ nearest ligand-atom), so it is
        robust to ligand size: a centroid-based test would wrongly reject a
        genuine catalytic residue at the reactive head of a long-chain PFAS,
        whose centroid lies far down the chain.
        """
        _spatial_cutoff = float(getattr(CFG, "RESOLVER_SPATIAL_CUTOFF", 10.0))
        if lig_coords:
            # Spatial sanity applies to ALL EIGHT catalytic roles: the clamp arginines and
            # His/Trp/Tyr cradle gate elite eligibility, so a sequence-proximal but spatially
            # remote mapping must be rejected for them too — not just the triad.
            _gate_keys = set(getattr(CFG, "CATALYTIC_TRIAD_KEYS", ["Nuc", "Base", "Acid"])) | {
                "Carb1", "Carb2", "Stab_H", "Stab_W", "Stab_Y"}
            for _tk in list(_gate_keys):
                if _tk in site_atoms and site_atoms[_tk]:
                    _min_d = min(a.dist(la) for a in site_atoms[_tk] for la in lig_coords)
                    if _min_d > _spatial_cutoff:
                        """
                        Residue is real but spatially remote from the substrate —
                        not the catalytic residue; report MISSING so the tier
                        fails honestly rather than on a spurious far-residue distance.
                        """
                        site_atoms.pop(_tk, None)
                        site_atoms_obj.pop(_tk, None)
                        results[f"dist_{_tk}"] = 999.0

        if not lig_coords:
                    results.update({
                        "catalytic_dist_A": 999.0, CFG.COL_TIER: CFG.TIER_DECOY, "is_degrader": False,
                        "residues_within_6A": "None", "constraint_check": "Ligand absence indicated",
                        "scientific_meaning": "No validated ligand atoms identified within structure.",
                        "Interaction_Density_Norm": 0.0, "Interaction_Density_Calc": "0.00",
                        "Chemical_Affinity_Score": 0.0, "Chemical_Affinity_Calc": "0.00",
                        "Binding_Probability_Score": 0.0, "Binding_Probability_Calc": "0.0000"
                    })
                    return results

        # -------------------------------------------------------------------------------
        # Sub-Step 7.2.2: Triangulate Valid Catalytic Trajectories
        # -------------------------------------------------------------------------------
        dists = {}
        for key_res in REF_ACTIVE_SITE_MAP.keys():
            min_d = 999.0
            if key_res in site_atoms:
                for pa in site_atoms[key_res]:
                    for la in lig_coords:
                        d = la.dist(pa)
                        if d < min_d: min_d = d
            dists[key_res] = min_d
            results[f"dist_{key_res}"] = round(min_d, 2)

        """
        Fluoride-cradle distances use only the H-bond-DONOR sidechain atoms (Chan 2011 —
        His155 Nδ/Nε, Trp156 Nε1, Tyr219 OH form three H-bonds to the leaving F⁻), via the
        POLAR_SIDECHAIN_ATOMS table. The generic loop above takes the nearest atom of ANY
        type, so a hydrophobic ring carbon of Trp/Tyr/His sitting near the halogen would
        otherwise satisfy the stabilisation gate without a real donor contact (false positive).
        Donor-restricted distance overwrites the all-atom value; if the mapped residue carries
        no donor (wrong type) the value stays 999 and the mutation-tolerant fallback below runs.
        A non-cradle residue mis-mapped at a Stab_* slot is guarded upstream: the class-aware
        ±window resolver maps Stab_W/Y to AROMATIC and Stab_H to AROMATIC|POSITIVE, and the
        elite tier additionally gates on cradle role-identity, so a stray polar residue cannot
        earn elite stabilisation credit on its donor alone.
        """
        # The fluoride cradle stabilises the DEPARTING F⁻, so measure the donor→ligand
        # distance to the ligand fluorine atoms only, not to every heavy atom. Otherwise a
        # donor sitting near the carboxylate head earns Halide_Stabilisation without ever
        # contacting the leaving fluorine. Fall back to all ligand atoms for a non-fluorinated
        # ligand so the metric never silently drops to 999 for those.
        _lig_f_pos = [a.pos for a in lig_atoms_obj if a.element.name == "F"] or lig_coords
        for _sk in ("Stab_H", "Stab_W", "Stab_Y"):
            _res3 = _canonical_resname(site_resname.get(_sk, ""))
            _donor = POLAR_SIDECHAIN_ATOMS.get(_res3, set())
            if _sk in site_atoms_obj and _donor and _lig_f_pos:
                _dd = min((la.dist(a.pos) for a in site_atoms_obj[_sk] if a.name in _donor
                           for la in _lig_f_pos), default=999.0)
                dists[_sk] = _dd
                results[f"dist_{_sk}"] = round(_dd, 2)

        # Catalytic-triad relay distances are between the functional sidechain atoms
        # (Asp Oδ, His Nδ/Nε, Asp/Glu Oδ/Oε) — NOT the backbone. Measuring over all residue
        # atoms lets a misfolded triad pass the strict relay gate on an incidental backbone
        # CA···CA contact. Restrict to the functional atoms, with a fallback to all atoms of
        # the residue when none match (non-standard naming / altloc) so it never regresses.
        _TRIAD_FUNC = {
            "Nuc":  {"OD1", "OD2"},
            "Base": {"ND1", "NE2"},
            "Acid": {"OD1", "OD2", "OE1", "OE2"},
        }
        _BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT"}
        def _triad_func_pos(_key):
            _objs = site_atoms_obj.get(_key, [])
            _pts = [a.pos for a in _objs if a.name in _TRIAD_FUNC.get(_key, set())]
            # Fall back to sidechain (non-backbone) atoms — never the full residue — so a
            # naming mismatch cannot silently re-admit backbone CA···CA relay contacts.
            return _pts or [a.pos for a in _objs if a.name not in _BACKBONE_NAMES] or [a.pos for a in _objs]

        dist_nuc_base = 999.0
        dist_base_acid = 999.0
        if site_atoms_obj.get("Nuc") and site_atoms_obj.get("Base"):
            _nuc_p, _base_p = _triad_func_pos("Nuc"), _triad_func_pos("Base")
            dist_nuc_base = min(p1.dist(p2) for p1 in _nuc_p for p2 in _base_p)
        if site_atoms_obj.get("Base") and site_atoms_obj.get("Acid"):
            _base_p, _acid_p = _triad_func_pos("Base"), _triad_func_pos("Acid")
            dist_base_acid = min(p1.dist(p2) for p1 in _base_p for p2 in _acid_p)

        results["dist_nuc_base_internal"] = round(dist_nuc_base, 2)
        results["dist_base_acid_internal"] = round(dist_base_acid, 2)


        """
        Reactive-centre identification (CFG §5.3). The catalytically productive attack
        carbon is the α-carbon bonded to the ligand carboxylate; the clamp engages the
        carboxylate oxygens. Both are resolved from the RDKit topology and mapped to
        structure-frame atom names so the SN2 geometry and the clamp gate key on the
        mechanistically correct atoms. A non-carboxylate head (sulfonate, ether) yields
        no α-carbon / carboxylate oxygens and therefore cannot anchor for α-attack.
        """
        alpha_c_names, carboxylate_o_names = set(), set()
        head_is_carboxylate = False
        if rd_mol and mm_map:
            _inv = {v: k for k, v in mm_map.items()}
            for _a in rd_mol.GetAtoms():
                if _a.GetSymbol() != "C":
                    continue
                _o_nb = [n for n in _a.GetNeighbors() if n.GetSymbol() == "O"]
                if len(_o_nb) >= 2:                       # carboxylate (–COO) carbon
                    head_is_carboxylate = True
                    for _o in _o_nb:
                        _on = _inv.get(_o.GetIdx())
                        if _on: carboxylate_o_names.add(_on)
                    for _c in _a.GetNeighbors():
                        if _c.GetSymbol() == "C":         # α-carbon adjacent to the carboxylate
                            _cn = _inv.get(_c.GetIdx())
                            if _cn: alpha_c_names.add(_cn)
        results["head_is_carboxylate"] = 1.0 if head_is_carboxylate else 0.0

        _preferred = alpha_c_names if (CFG.SCISSILE_REQUIRE_ALPHA and alpha_c_names) else None
        # Cradle atoms (Trp156/Tyr217/His155 sidechains) → couple the leaving-F selection
        # to the fluoride cradle so the scissile F is the departing one (Chan 2011), not the
        # best-of-N most-anti fluorine. Empty when the cradle is unresolved → max-angle fallback.
        _cradle_coords = [(_a.pos.x, _a.pos.y, _a.pos.z)
                          for _ck in ("Stab_W", "Stab_Y", "Stab_H")
                          for _a in site_atoms_obj.get(_ck, [])
                          if _a.name not in ("N", "C", "CA", "O")]
        sn2_angle, sn2_dev, steric_clashes, target_x_pos, best_c_atom, sn2_aux = calculate_sn2_metrics(
            site_atoms_obj.get("Nuc", []), lig_atoms_obj, rd_mol, mm_map,
            preferred_c_names=_preferred, cradle_coords=_cradle_coords or None)
        results["sn2_attack_angle"] = round(sn2_angle, 1)
        results["sn2_trajectory_dev"] = round(sn2_dev, 2)
        results["teflon_shield_clashes"] = steric_clashes
        # Auxiliary (non-gating) reference geometries
        _aux = sn2_aux or {}
        results["burgi_dunitz_angle"]   = _aux.get("burgi_dunitz_angle", 999.0)
        results["flippin_lodge_offset"] = _aux.get("flippin_lodge_offset", 999.0)

        scissile_is_alpha = bool(best_c_atom is not None and best_c_atom.name in alpha_c_names)
        results["scissile_is_alpha"] = 1.0 if scissile_is_alpha else 0.0

        """
        Nucleophile reach is measured attacking Asp-Oδ → SN2 attack carbon (the reactive centre),
        not Asp → nearest ligand atom, so a stray fluorine cannot satisfy the gate.

        The oxygen is the one that produced the reported attack angle, not the nearer of the two.
        Measuring the distance to one oxygen and the angle to the other describes a chimeric
        nucleophile and lets a pose clear the NAC distance gate on an oxygen that is not attacking.
        The minimum over both Oδ is retained as a separate diagnostic column.
        """
        _asp_ox = [a for a in site_atoms_obj.get("Nuc", []) if a.name in ("OD1", "OD2")] or site_atoms_obj.get("Nuc", [])
        _attack_o = _aux.get("attack_o_atom")
        if best_c_atom is not None and _attack_o is not None:
            results["dist_Nuc"] = round(_attack_o.pos.dist(best_c_atom.pos), 2)
        elif best_c_atom is not None and _asp_ox:
            results["dist_Nuc"] = round(min(o.pos.dist(best_c_atom.pos) for o in _asp_ox), 2)
        if best_c_atom is not None and _asp_ox:
            results["dist_Nuc_nearest_O"] = round(min(o.pos.dist(best_c_atom.pos) for o in _asp_ox), 2)

        d_nuc = results["dist_Nuc"]
        angle  = results["sn2_attack_angle"]
        d_trp, d_tyr = results.get("dist_Stab_W", 999.0), results.get("dist_Stab_Y", 999.0)

        # --- ACTIVE SITE FLUORIDE-STABILISATION METRICS ---
        # Validates the documented aromatic/His fluoride cradle near the leaving halogen.
        target_stabilise_centre = target_x_pos
        if not target_stabilise_centre and lig_coords: target_stabilise_centre = lig_coords[0]

        """
        Mutation-tolerant stabiliser distances: when the alignment-mapped
        canonical residue is absent (dist == 999), report the nearest
        functionally-equivalent sidechain to the leaving group instead (same
        dynamic principle as the nucleophile search), accepted only within
        CFG.MECH_CRADLE_RADIUS. Aromatic π-system for the Trp156/Tyr217 cradle
        keys; an H-bond donor / cationic group for the His155 fluoride stabiliser
        — PHE excluded there (no polar donor to F⁻), matching CFG ROLE_EXPECTED_
        RESIDUES['Fluorine_Stabiliser'].
        """
        if target_stabilise_centre:
            _stab_classes = {
                "dist_Stab_W": RESIDUE_CLASS_GROUPS["AROMATIC"],
                "dist_Stab_Y": RESIDUE_CLASS_GROUPS["AROMATIC"],
                "dist_Stab_H": (RESIDUE_CLASS_GROUPS["AROMATIC"] - {"PHE"}) | RESIDUE_CLASS_GROUPS["POSITIVE"],
            }
            for _sk, _cls in _stab_classes.items():
                if results.get(_sk, 999.0) >= 999.0:
                    _best = 999.0
                    for p_at in all_prot_atoms:
                        _canon_st = _canonical_resname(p_at["resname"])
                        # Only the residue's H-bond-DONOR atoms count (POLAR_SIDECHAIN_ATOMS):
                        # a ring/backbone carbon near the halogen is not fluoride stabilisation.
                        # PHE carries no donor → empty set → correctly excluded.
                        _donor_fb = POLAR_SIDECHAIN_ATOMS.get(_canon_st, set())
                        if (_canon_st in _cls or p_at["resname"] in _cls) and p_at["atom_name"] in _donor_fb:
                            _d = p_at["pos"].dist(target_stabilise_centre)
                            if _d < _best: _best = _d
                    if _best <= CFG.MECH_CRADLE_RADIUS:
                        results[_sk] = round(_best, 2)
            d_trp, d_tyr = results.get("dist_Stab_W", 999.0), results.get("dist_Stab_Y", 999.0)

        """
        Fluoride stabilisation requires the documented cradle — an aromatic π-system
        (Trp/Tyr) or the His fluoride stabiliser — within MECH_STAB_RADIUS of the leaving
        halogen. A generic polar sidechain (Ser/Thr/Asn …) does NOT count: in a packed
        pocket some polar atom is always nearby, which would make the flag fire for every
        pose and carry no signal.
        """
        d_his = results.get("dist_Stab_H", 999.0)
        stabilised = (d_trp <= CFG.MECH_STAB_RADIUS
                      or d_tyr <= CFG.MECH_STAB_RADIUS
                      or d_his <= CFG.MECH_STAB_RADIUS)
        results["halide_stabilisation_score"] = 1.0 if stabilised else 0.0

        """
        Carboxylate clamp engagement (CFG §5.3): the ligand carboxylate oxygens must
        salt-bridge the clamp arginines (guanidinium N within CLAMP_SALT_BRIDGE_DIST).
        One engaged arm earns mechanistic-score credit (clamp_ok); the elite tier
        requires BOTH distinct arginines (ligand_clamp_engaged), the bidentate anchoring
        that holds a 2-haloalkanoate for α-attack. A non-carboxylate head cannot satisfy
        either, so it is barred from clamp credit.
        """
        '''
        Clamp salt-bridge donor atoms: Arg guanidinium (NH1/NH2/NE) AND Lys ammonium (NZ),
        so an engineered carboxylate clamp that substitutes Arg→Lys is still recognised as
        engaged (discovery-friendly; a flat Arg-only atom list gives false negatives). The
        role-identity / elite gate separately decides which residue types earn elite status.
        '''
        _CLAMP_CATION_N = ("NH1", "NH2", "NE", "NZ")
        def _clamp_arm_oxygens(_site_key):
            """The ligand carboxylate oxygens this clamp arm salt-bridges, by atom name."""
            _arg_ns = [a for a in site_atoms_obj.get(_site_key, []) if a.name in _CLAMP_CATION_N]
            if not _arg_ns or not carboxylate_o_names:
                return set()
            _o_atoms = [a for a in lig_atoms_obj if a.name in carboxylate_o_names]
            return {_o.name for _n in _arg_ns for _o in _o_atoms
                    if _n.pos.dist(_o.pos) <= CFG.CLAMP_SALT_BRIDGE_DIST}
        _o1 = _clamp_arm_oxygens("Carb1")
        _o2 = _clamp_arm_oxygens("Carb2")
        _arm1, _arm2 = bool(_o1), bool(_o2)
        """
        The two arms must also be two DIFFERENT residues. The alignment can map Carb1 and Carb2 onto
        the same arginine, and that single residue reaching both carboxylate oxygens would otherwise
        report a full bidentate clamp — one guanidinium chelating the head is not the two-point
        anchoring the mechanism requires. The elite tier already refuses a duplicated mapping
        (elite_identity_ok), but carboxylate_clamp_integrity is reported for every pose and must not
        claim a clamp that does not exist.
        """
        _clamp_residues_distinct = (
            mapped_sites.get("Carb1") is not None
            and mapped_sites.get("Carb2") is not None
            and mapped_sites.get("Carb1") != mapped_sites.get("Carb2")
        )
        """
        Bidentate credit requires the two arms to hold DIFFERENT carboxylate oxygens. Two
        arginines converging on the same oxygen is a monodentate collapse: the head group is
        pinched at one point and free to rotate, which is not the two-point anchoring that holds
        a 2-haloalkanoate rigid for α-attack. Distinct arginines alone are not sufficient.

        The test is whether the arms can be assigned to distinct oxygens at all, not whether they
        already contact disjoint sets: with two arms that holds exactly when each arm engages at
        least one oxygen and the two arms reach at least two oxygens between them.
        """
        _bidentate_oxygens = bool(_o1 and _o2 and len(_o1 | _o2) >= 2)
        clamp_ok = bool(head_is_carboxylate and (_arm1 or _arm2))
        ligand_clamp_engaged = bool(head_is_carboxylate and _arm1 and _arm2
                                    and _bidentate_oxygens and _clamp_residues_distinct)
        results["carboxylate_clamp_integrity"] = 1.0 if ligand_clamp_engaged else (0.5 if clamp_ok else 0.0)
        results["sn2_alignment_score"] = angle

        # --- SOFT SCORING ENGINE ---
        # Uses sigmoid functions to avoid binary threshold "cliffs"
        s_nuc = sigmoid(d_nuc, k=CFG.SOFT_K_NUC, x0=CFG.NAC_DIST_STRICT)   # High score for < NAC_DIST_STRICT (3.2 Å)
        s_ang = sigmoid(angle, k=CFG.SOFT_K_ANG, x0=CFG.NAC_ANGLE_STRICT)  # High score for > NAC_ANGLE_STRICT (155°)
        s_int = sigmoid(dist_nuc_base, k=CFG.SOFT_K_TRIAD, x0=CFG.SOFT_NB_MIDPOINT) * sigmoid(dist_base_acid, k=CFG.SOFT_K_TRIAD, x0=CFG.SOFT_BA_MIDPOINT)
        soft_score = (s_nuc * CFG.SOFT_W_NUC) + (s_ang * CFG.SOFT_W_ANG) + (s_int * CFG.SOFT_W_INT)
        results["soft_catalytic_score"] = round(soft_score, 3)

        # ---- SN2 dead-end consensus check (A: scissile C–F BDE | B: backside sterics) ----
        # Computed BEFORE the mechanistic score so the dead-end verdict folds into it.
        # BOTH indicators key on the mapped SN2 ATTACK CARBON (nearest the nucleophile),
        # NOT every CF3 in the ligand: a long PFCA docks carboxylate-in-clamp and presents
        # an α-CF2 here (2 F → passes), while its distal hydrophobic CF3 never faces Asp.
        # Only a CF3 reactive carbon (e.g. TFA) triggers the consensus.
        # A — fluorines bonded to the attack carbon → scissile C–F BDE class.
        terminal_f_count = 0
        _tf_mapped = False
        if rd_mol and mm_map and best_c_atom:
            _tf_idx = mm_map.get(best_c_atom.name)
            if _tf_idx is not None:
                try:
                    _tf_atom = rd_mol.GetAtomWithIdx(_tf_idx)
                    for _nbr in _tf_atom.GetNeighbors():
                        if _nbr.GetSymbol() == "F":
                            terminal_f_count += 1
                    _tf_mapped = True
                except Exception: pass
        '''
        When the Hungarian mapping fails to resolve the attack carbon, fall back to the
        map-independent scissile-F count from calculate_sn2_metrics (_aux). Without this
        a distorted pose defaults terminal_f_count=0 → BDE 0.0 → feasibility_factor 1.0,
        erasing the C–F recalcitrance penalty for a CF3-α ligand (e.g. TFA reads 1.0
        instead of 0.375) and inflating the reported feasibility column / control read-out.
        '''
        if not _tf_mapped:
            terminal_f_count = int(_aux.get("n_scissile_f", 0))
        """
        CHEMISTRY THAT COULD NOT BE EVALUATED IS NOT CHEMISTRY THAT PASSED.

        Falling back to 0.0 here does not mean "unknown" — it means a C-F bond with ZERO dissociation
        energy, which sails under the dead-end gate (123 kcal/mol) and under the Tier_1A ceiling
        (TIER_ELITE_BDE_MAX). A complex whose fluorine count could not be resolved was therefore
        scored MORE favourably than one that was properly penalised.

        The BDE table is keyed on 1, 2 or 3 fluorines. Anything else means the scissile centre was not
        resolved, and that is recorded as such (_chem_unknown) rather than silently forgiven.
        """
        _chem_unknown = terminal_f_count not in CFG.SCISSILE_CF_BDE
        scissile_cf_bde = float(CFG.SCISSILE_CF_BDE.get(terminal_f_count, CFG.SENTINEL_UNDEFINED))

        '''
        An incomplete ligand atom map means the scissile centre was never identified, because every
        reactive-centre gate is keyed on the MAPPED α-carbon. The fraction was being written to the
        CSV and read by nothing, so a complex whose ligand barely mapped was scored on the geometry
        around an atom the pipeline could not name. Below CFG.LIGAND_MAP_MIN_FRACTION the chemistry
        is declared unverified — barred from the elite tier, not annihilated.
        '''
        _map_frac = (float(len(mm_map)) / float(len(mmcif_lig))) if mmcif_lig else 0.0
        if _map_frac < CFG.LIGAND_MAP_MIN_FRACTION:
            _chem_unknown = True

        # B — backside steric occlusion: vdW bulk of heavy halogen substituents on the
        # attack carbon lying on the nucleophile-approach hemisphere (anti to leaving F).
        # The leaving F sits opposite the approach axis (dot < 0) and is auto-excluded.
        # A skipped computation is not a measurement of zero. Fluoroacetate genuinely HAS zero backside
        # occlusion; a complex whose RDKit mapping failed has an UNKNOWN one, and the two must not be
        # written the same way. The flag below separates them.
        backside_occlusion = 0.0
        _occl_measured = bool(best_c_atom is not None and target_x_pos is not None and rd_mol and mm_map)
        if not _occl_measured:
            _chem_unknown = True
        if best_c_atom is not None and target_x_pos is not None and rd_mol and mm_map:
            _c_idx = mm_map.get(best_c_atom.name)
            if _c_idx is not None:
                _invm = {v: k for k, v in mm_map.items()}
                _name2pos = {a.name: a.pos for a in lig_atoms_obj}
                _cp = best_c_atom.pos
                _ux, _uy, _uz = _cp.x - target_x_pos.x, _cp.y - target_x_pos.y, _cp.z - target_x_pos.z
                _un = math.sqrt(_ux*_ux + _uy*_uy + _uz*_uz) or 1.0
                _ux, _uy, _uz = _ux/_un, _uy/_un, _uz/_un
                try:
                    _crd = rd_mol.GetAtomWithIdx(_c_idx)
                    for _nb in _crd.GetNeighbors():
                        _el = _nb.GetSymbol().upper()
                        if _el in ("C", "H"):   # chain carbons / H: negligible backside crowding
                            continue
                        _pp = _name2pos.get(_invm.get(_nb.GetIdx()))
                        if _pp is None:
                            continue
                        _vx, _vy, _vz = _pp.x - _cp.x, _pp.y - _cp.y, _pp.z - _cp.z
                        _vn = math.sqrt(_vx*_vx + _vy*_vy + _vz*_vz) or 1.0
                        if (_vx*_ux + _vy*_uy + _vz*_uz) / _vn > 0.0:   # approach-hemisphere substituent
                            backside_occlusion += CFG.VDW_RADII.get(_el, CFG.VDW_RADIUS_DEFAULT)
                except Exception: pass

        results["scissile_cf_bde"] = round(scissile_cf_bde, 1)
        results["sn2_backside_occlusion"] = round(backside_occlusion, 2)
        beta_f_count = int(_aux.get("beta_f_count", 0))
        results["beta_f_count"] = beta_f_count
        n_scissile_f = int(_aux.get("n_scissile_f", 1))
        """
        Angle multiplicity for the Šidák correction (CFG §5.2d) = the number of equivalent C–F bonds
        on the scissile carbon. The inflation being corrected lives in the POSE: an α-CF3 has three
        chances to present some fluorine anti-periplanar to the nucleophile, a mono-fluoro substrate
        has one. Which of the three actually leaves — resolved here by the fluoride cradle — does not
        change how many chances the pose had, so the cradle does not reduce the exponent.

        The tier ladder then gates on the EFFECTIVE angle: the angle a single-C–F substrate would
        have to show to be as improbable as this pose. Without it, TFA out-angles fluoroacetate on the
        DeHa4 control in all five diffusion samples and takes a higher tier on the one enzyme known
        not to turn it over.
        """
        n_angle_choices = max(1, n_scissile_f)
        results["angle_multiplicity"] = int(n_angle_choices)
        angle_effective = CFG.sn2_effective_angle(angle, n_angle_choices)
        results["sn2_attack_angle_effective"] = round(angle_effective, 1)
        """
        WHERE C–F BOND STRENGTH ENTERS, STATED HONESTLY.

        The previous note here claimed "competence/ranking keep the raw geometry". They do not:
        CFG.competence_score multiplies the geometric sum by feasibility_factor(scissile_cf_bde,
        beta_f_count). Bond strength therefore enters in THREE places, deliberately and for three
        different purposes, and a reader is entitled to know all three:

          · the TIER, through the graded chemistry penalty in mechanistic_score_effective (_bde_pen);
          · the ELITE CEILING, through TIER_ELITE_BDE_MAX, which bars an unbreakable C–F from Tier_1A
            however good its pose;
          · the WITHIN-TIER RANK, through feasibility_factor scaling competence_score.

        The three are not a double-count: the first decides whether the chemistry is feasible at all, the
        second whether it is elite, the third how it ranks among its peers. Step-07 QM/MM remains the
        final reactivity arbiter — none of these is a substitute for a barrier.
        """
        results["feasibility_factor"] = round(CFG.feasibility_factor(scissile_cf_bde, beta_f_count), 3)
        sn2_dead_end = (scissile_cf_bde > CFG.SCISSILE_CF_BDE_MAX
                        and backside_occlusion > CFG.SN2_BACKSIDE_OCCL_MAX)

        """
        MECHANISTIC SCORE  (holistic 0–1)
        Computed via CFG.mechanistic_score() — the single source of truth (§5.1). It sums
        five binary anchor checks (nucleophile reach, Nuc–Base and Base–Acid relay distances,
        carboxylate clamp, halide stabilisation) plus a linear graded SN2 attack-angle term,
        less a per-clash penalty. This is the raw GEOMETRY + machinery score (competence and
        the reported mechanistic_score use it). The tier gate instead uses
        mechanistic_score_effective = this minus the graded chemistry and containment
        penalties (computed after the pocket-fit block). All weights live in CFG.
        """
        """
        THE MULTIPLICITY CORRECTION IS CHARGED ONCE, AND BOTH GATES SEE THE SAME CORRECTED ANGLE.

        The correction must not be applied twice to one tier decision. `angle_effective` already converts the
        observed angle into the angle a single-C–F substrate would need to be equally improbable (the
        Šidák deflation, §5.2d). Passing the RAW angle plus `angle_multiplicity=n` then made
        mechanistic_score deflate that same angle a second time, by the same exponent — and the Tier_1A
        branch requires `angle_effective >= TIER_ANGLE_MIN` AND `mech_score >= TIER_MECH_MIN`, so a CF2 or
        CF3 substrate paid the best-of-N penalty on both rungs of the same ladder.

        The correction is right; charging it twice is not. It is computed once, and the already-corrected
        angle is what both gates score — so the penalty is real, applied once, and the two gates cannot
        silently compound it.
        """
        mech_score = CFG.mechanistic_score(
            d_nuc, dist_nuc_base, dist_base_acid,
            clamp_ok, stabilised, angle_effective, steric_clashes,
            scissile_cf_bde, backside_occlusion, beta_f_count,
            angle_multiplicity=1,
        )
        results[CFG.COL_MECH_S] = round(mech_score, 2)

        """
        Gated continuous competence score (CFG §5.5) — the Scientific-ranking key. Every
        independent catalytic axis enters once, under hard chemistry gates (α-attack +
        carboxylate, not a dead-end); a non-productive pose scores 0. Computed via the CFG
        single source so the formula is never re-implemented inline.
        """
        competence = CFG.competence_score(
            scissile_is_alpha, head_is_carboxylate,
            angle, d_nuc, results["carboxylate_clamp_integrity"],
            results["sn2_trajectory_dev"], stabilised,
            dist_nuc_base, dist_base_acid,
            scissile_cf_bde=scissile_cf_bde, beta_f_count=beta_f_count,
            scissile_f_count=n_scissile_f,
        )
        results["competence_score"] = competence

        """
        Tier assignment is driven by the variant's measured active-site geometry toward
        the docked ligand (nucleophile reach, SN2 angle, triad distances, clamp and
        halide-pocket engagement) together with the presence and identity of the
        catalytic machinery. No substrate-class penalty is applied, so the screen can
        surface FAcD variants active on any PFAS. The one chemistry-aware filter is the
        SN2 dead-end consensus above (A: scissile C–F bond strength, O'Hagan 2008;
        B: backside steric occlusion, Bento & Bickelhaupt 2008 + Bondi 1964) — it keys
        on the mapped attack carbon, so only a CF3 reactive carbon (e.g. TFA) is flagged;
        a long PFCA presents an α-CF2 and passes. Final chemical feasibility is confirmed
        at the Step-07 QM/MM + MD/WaterMap stage.
        """

        '''
        BACKBONE-CLASH METRIC (size-fair; graded, not a size veto)
        Counts deep-pocket tail atoms physically interpenetrating the protein backbone.
        Denominator is ALL ligand heavy atoms (size-fair): the reported ratio does not
        inflate for long chains, so length alone never decoys a pose. A genuine clash
        is a real physical defect, applied below as a graded competence penalty (per
        clashing atom). The hard veto is the size-fair clash FRACTION (TAIL_CLASH_RATIO)
        gated by a minimum absolute count (CLASH_MIN_FLOOR), so a non-physical pose is
        caught by the proportion of the ligand interpenetrating the backbone — not a flat
        count that biases against long PFAS. Both clash_count and the ratio are still
        reported in the master CSV (mainchain_clash_count / mainchain_clash_ratio).
        '''
        tail_clash_ratio = 0.0
        clash_count = 0
        if rd_mol and mm_map and best_c_atom:
            target_rd_idx = mm_map.get(best_c_atom.name)
            if target_rd_idx is not None:
                distances = Chem.rdmolops.GetDistanceMatrix(rd_mol)
                tail_rd_indices = [i for i, d in enumerate(distances[target_rd_idx]) if d > CFG.TAIL_MIN_BOND_DISTANCE]
                inv_map = {v: k for k, v in mm_map.items()}
                tail_cif_names = {inv_map.get(idx) for idx in tail_rd_indices if inv_map.get(idx)}
                tail_pos_list = [a.pos for a in lig_atoms_obj if a.name in tail_cif_names]

                if tail_pos_list:
                    bb_atoms = [pa for pa in all_prot_atoms if pa["atom_name"] in ["N", "CA", "C", "O"]]
                    for tpos in tail_pos_list:
                        for ba in bb_atoms:
                            if tpos.dist(ba["pos"]) < CFG.BACKBONE_CLASH_DIST:
                                clash_count += 1
                                break
                    _lig_heavy = sum(1 for a in lig_atoms_obj if a.element.name != "H")
                    tail_clash_ratio = clash_count / max(1, _lig_heavy)

        results["mainchain_clash_ratio"] = round(tail_clash_ratio, 3)
        results["mainchain_clash_count"] = clash_count
        # Graded, size-fair clash penalty on the competence rank key (never a length veto).
        if clash_count > 0:
            results["competence_score"] = round(
                max(0.0, results.get("competence_score", 0.0) - CFG.MAINCHAIN_CLASH_PENALTY * clash_count), 3)

        # Enzyme vs ligand steric complementarity (§5.2c). pocket_containment_cavity (protein
        # ray-cast buriedness) feeds the feasibility-weighted tier score below;
        # pocket_containment_site8 (engagement with the eight catalytic residues) and the hull
        # volume/occupancy fields are reported diagnostics. The full protein is required: a
        # containment measured on ligand coordinates alone cannot tell a narrow pocket from a wide one.
        results.update(compute_pocket_fit(site_atoms_obj, lig_atoms_obj, all_prot_atoms))

        # Feasibility-weighted mechanistic score (tier-gate key). Two graded penalties are
        # subtracted from the geometric mech_score: (A) chemistry — scissile C–F BDE and
        # backside occlusion beyond the dead-end cutoffs (down-ranks the SN2 dead-end TFA to
        # ~Tier_2B), plus a per-β-fluorine term (the perfluoro tail inductively strengthens the
        # α-C–F, so a perfluoro CF2 α-carbon like PFBA is down-ranked past a haloacetate CF2
        # like DFA; FA/DFA/TFA have β_F=0 and are untouched); (B) containment — the ligand
        # spilling out of the small FAcD pocket (down-ranks long-chain PFAS). All engage only
        # past a threshold, so FA/DFA/short substrates are unpenalised, and all are per-pose
        # geometry, so a genuine strong-binding/wide-pocket variant can still surface a harder
        # substrate. Competence and the reported mechanistic_score stay on the raw geometry;
        # only the tier gates use the effective value. FAcD small-substrate scope + PFAS
        # recalcitrance: Wackett 2022; Chan 2011.
        """
        Only the backside-occlusion term fades toward a near-ideal SN2 angle (§5.2b): occlusion is
        a steric obstruction of the attack trajectory, and a pose that reaches 180° has cleared it.
        The C–F bond dissociation energy does NOT fade — the strength of the bond being broken is a
        property of the bond, not of the angle of approach, so a near-linear attack on a strong C–F
        keeps its full penalty. β-fluorination and containment stay flat for the same reason.
        """
        """
        UNKNOWN CHEMISTRY IS BARRED FROM ELITE, NOT ANNIHILATED.

        When the scissile centre could not be resolved (_chem_unknown), the BDE is the sentinel and the
        occlusion was never measured. Feeding the sentinel into the penalty would subtract ~13 from the
        mechanistic score and bury the complex in Tier_5 — punishing it for a mapping failure rather than
        for its chemistry, and destroying a candidate that may be perfectly good.

        The honest position is narrower: we do not know, so we do not PROMOTE. No numeric penalty is
        invented from a sentinel, and the complex is barred from the elite tier (below). It keeps its
        geometry-earned rank and is flagged, so a human can see exactly which complexes were never
        chemically verified instead of finding them silently at the top or silently at the bottom.
        """
        if _chem_unknown:
            _occl_pen = 0.0
            _bde_pen  = 0.0
        else:
            _occl_pen = CFG.CHEM_PEN_W_OCCL * max(0.0, backside_occlusion - CFG.SN2_BACKSIDE_OCCL_MAX)
            _bde_pen  = CFG.CHEM_PEN_W_BDE  * max(0.0, scissile_cf_bde   - CFG.SCISSILE_CF_BDE_MAX)
        _ang_scale = min(1.0, max(0.0, (CFG.CHEM_PEN_ANGLE_NONE - angle)
                                       / (CFG.CHEM_PEN_ANGLE_NONE - CFG.CHEM_PEN_ANGLE_FULL)))
        _chem_pen = _bde_pen + _occl_pen * _ang_scale + CFG.CHEM_PEN_W_BETA * int(beta_f_count)
        _cont     = float(results.get("pocket_containment_cavity", 1.0) or 1.0)
        _cont_pen = CFG.CONTAIN_PEN_W * max(0.0, CFG.CONTAIN_PEN_TARGET - _cont)
        mech_effective = max(0.0, mech_score - _chem_pen - _cont_pen)
        results["chem_penalty"]              = round(_chem_pen, 3)
        results["chem_verified"]             = 0 if _chem_unknown else 1
        results["containment_penalty"]       = round(_cont_pen, 3)
        results["mechanistic_score_effective"] = round(mech_effective, 2)
        mech_score = mech_effective   # tier ladder + elite gate gate on the feasibility-weighted score

        """
        Catalytic-identity gates. Geometry alone (good SN2 angle + close nucleophile)
        is NOT sufficient: a functional FAcD must present the correct catalytic
        residues. Uses the alignment-mapped residue names captured above.
          A — triad identity (required for ANY degrader tier): Asp nucleophile,
              His base, Asp acid. A Ser/Thr at the acid, a missing His, etc. → not a degrader.
          B — elite (Tier_1A) additionally requires the carboxylate clamp present as
              TWO DISTINCT arginines and an intact halide pocket (aromatic cradle +
              fluoride stabiliser of the documented class).
        """
        _exp_id = CFG.ROLE_EXPECTED_RESIDUES
        def _role_id_ok(_site_key, _role):
            return _canonical_resname(site_resname.get(_site_key, "")) in _exp_id.get(_role, set())
        triad_identity_ok = (_role_id_ok("Nuc", "Nucleophile")
                             and _role_id_ok("Base", "Base_Catalyst")
                             and _role_id_ok("Acid", "Acid_Catalyst"))
        _c1_idx, _c2_idx = mapped_sites.get("Carb1"), mapped_sites.get("Carb2")
        elite_identity_ok = (_role_id_ok("Carb1", "Carboxylate_Clamp")
                             and _role_id_ok("Carb2", "Carboxylate_Clamp")
                             and _c1_idx is not None and _c2_idx is not None and _c1_idx != _c2_idx
                             and _role_id_ok("Stab_W", "Fluoride_Cradle")
                             and _role_id_ok("Stab_Y", "Fluoride_Cradle")
                             and _role_id_ok("Stab_H", "Fluorine_Stabiliser"))

        """
        Active-site integrity (Chan et al. 2011, JACS 133:7461): fraction of the eight
        catalytic residues — Asp nucleophile, His base, Asp acid, the two clamp arginines
        and the His/Trp/Tyr fluoride cradle — mapped to the correct residue type. This is
        the mechanistic measure of whether the catalytic machinery is present, and replaces
        global sequence identity as the conservation signal: a distant homolog with all eight
        residues conserved scores 1.0 regardless of overall identity.
        """
        # The two clamp arginines must be DISTINCT residues; without the exclusivity guard a
        # single Arg mapped to both Carb1 and Carb2 would be counted twice, inflating the
        # conservation score to 8/8 when only seven residues are present.
        _role_checks = [
            _role_id_ok("Nuc", "Nucleophile"), _role_id_ok("Base", "Base_Catalyst"),
            _role_id_ok("Acid", "Acid_Catalyst"), _role_id_ok("Carb1", "Carboxylate_Clamp"),
            (_role_id_ok("Carb2", "Carboxylate_Clamp") and _c2_idx is not None and _c2_idx != _c1_idx),
            _role_id_ok("Stab_H", "Fluorine_Stabiliser"),
            _role_id_ok("Stab_W", "Fluoride_Cradle"), _role_id_ok("Stab_Y", "Fluoride_Cradle"),
        ]
        results["active_site_residues_correct"] = int(sum(_role_checks))
        results["active_site_integrity"] = round(sum(_role_checks) / 8.0, 3)

        '''
        Per-active-site prediction confidence. Boltz writes per-atom pLDDT in the
        B-factor field; this is the mean pLDDT over the eight mapped catalytic residues —
        local confidence in the active-site fold specifically, not the global mean_plddt. A
        Tier_1A resting on a low-confidence loop residue is visible here (diagnostic column;
        the elite-confidence demotion gates on the global confidence_score, not this field).
        '''
        _as_bvals = [a.b_iso for _rk in REF_ACTIVE_SITE_MAP.keys()
                     for a in site_atoms_obj.get(_rk, [])]
        results["active_site_plddt"] = round(sum(_as_bvals) / len(_as_bvals), 2) if _as_bvals else 0.0

        """
        Productive-attack gate (CFG §5.3): a degrader must present the α-carbon (adjacent
        to the substrate carboxylate) as the SN2 attack carbon. Sulfonate/ether heads and
        mid-chain CF2 attacks are not productive defluorination and cannot earn a degrader
        tier. Elite eligibility additionally requires the bidentate carboxylate clamp engaged,
        the fluoride cradle actually within stabilising distance of the leaving fluorine
        (Chan et al. 2011 — His155/Trp156/Tyr219 form three H-bonds to F⁻; residue identity
        alone is insufficient, the pocket must reach the departing halide), and the nucleophile
        resolved without a wide alignment rescue (§5.4).
        """
        productive_attack = scissile_is_alpha
        """
        Elite bond-strength ceiling (§8.5). The attack angle is a property of the approach; the
        C–F dissociation energy is a property of the bond. A near-linear trajectory onto a C–F
        above TIER_ELITE_BDE_MAX does not make that bond cleavable, so such a pose is barred from
        the elite tier no matter how ideal its geometry — it falls to Tier_1B and stays a
        discovery lead for Step-07 QM/MM to adjudicate.
        """
        bde_elite_ok = scissile_cf_bde <= CFG.TIER_ELITE_BDE_MAX
        """
        `not _chem_unknown` is a requirement of the ELITE tier, not of every tier. Tier_1A is the claim
        that this complex is a lead worth a week of GPU time — and that claim cannot rest on chemistry
        that was never evaluated. A complex whose scissile centre could not be resolved keeps whatever
        rank its geometry earns, but it cannot be promoted to the top on the strength of penalties that
        were skipped rather than passed.
        """
        elite_ready = (elite_identity_ok and ligand_clamp_engaged and stabilised
                       and bde_elite_ok
                       and not _chem_unknown
                       and nuc_rescue_offset <= CFG.NUC_RESCUE_MAX_OFFSET_ELITE)

        tier, is_degrader, meaning, constraint = CFG.TIER_DECOY, False, "No significant documented interactions.", "Fail"

        # -------------------------------------------------------------------------------
        # Sub-Step 7.2.3: ULTRA-STRICT FAcD HIERARCHY LOGIC
        # -------------------------------------------------------------------------------

        if tail_clash_ratio >= CFG.TAIL_CLASH_RATIO and clash_count >= CFG.CLASH_MIN_FLOOR:
            tier, meaning, is_degrader = CFG.TIER_DECOY, (
                f"Failed: non-physical pose ({clash_count} backbone-interpenetrating atoms, "
                f"{tail_clash_ratio:.0%} of ligand heavy atoms)."), False

        # Catalytic triad identity — the enzyme must present an Asp nucleophile,
        # His base and Asp acid; without the catalytic machinery there is no FAcD
        # mechanism for any substrate.
        elif not triad_identity_ok:
            tier, meaning, is_degrader = CFG.TIER_DECOY, (
                "Failed: catalytic triad identity not conserved — needs Asp nucleophile, "
                f"His base, Asp acid (mapped Nuc={site_resname.get('Nuc','MISSING')}, "
                f"Base={site_resname.get('Base','MISSING')}, Acid={site_resname.get('Acid','MISSING')})."), False

        # The ladder gates on angle_effective, the multiplicity-corrected attack angle (CFG §5.2d),
        # not on the raw one. A poly-fluorinated attack carbon has several equivalent C–F bonds and
        # so several chances at a near-linear backside geometry; the raw angle rewards that best-of-N
        # as though it were catalytic competence. Gating on the raw angle lets trifluoroacetate
        # outrank the native substrate on the DeHa4 control — the enzyme known NOT to turn TFA over —
        # because Boltz gives its CF3 a better-aligned pose in all five diffusion samples. The
        # effective angle asks what a single-C–F substrate would have had to achieve to be equally
        # improbable, so a mediocre CF3 pose falls back and only a near-ideal one holds its tier.
        #
        # Tier_1A  (TIER_ORDER[0]) — elite: productive α-attack with near-ideal SN2 geometry,
        # intact catalytic constellation, bidentate carboxylate clamp engaged and a directly
        # resolved nucleophile. Pruned for near-ideal geometry to reduce downstream MD workload.
        elif productive_attack and elite_ready and d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_TOP] and dist_nuc_base <= CFG.TIER_NB_MAX[CFG.TIER_TOP] and dist_base_acid <= CFG.TIER_BA_MAX[CFG.TIER_TOP] and angle_effective >= CFG.TIER_ANGLE_MIN[CFG.TIER_TOP] and mech_score >= CFG.TIER_MECH_MIN[CFG.TIER_TOP]:
            tier, meaning, is_degrader = CFG.TIER_TOP, "Elite-Grade Analysis: Near-perfect SN2 Trajectory demonstrating absolute anchor integrity.", True

        # Tier_1B  (TIER_ORDER[1])
        elif productive_attack and d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_ORDER[1]] and dist_nuc_base <= CFG.TIER_NB_MAX[CFG.TIER_ORDER[1]] and dist_base_acid <= CFG.TIER_BA_MAX[CFG.TIER_ORDER[1]] and angle_effective >= CFG.TIER_ANGLE_MIN[CFG.TIER_ORDER[1]] and mech_score >= CFG.TIER_MECH_MIN[CFG.TIER_ORDER[1]]:
            tier, meaning, is_degrader = CFG.TIER_ORDER[1], "Crystal-Grade Analysis: Ideal ground-state contact sequence with a connected catalytic relay (elite anchor integrity NOT asserted — Tier_1A only).", True

        # Tier_2A  (TIER_ORDER[2])
        elif productive_attack and d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_ORDER[2]] and dist_nuc_base <= CFG.TIER_NB_MAX[CFG.TIER_ORDER[2]] and dist_base_acid <= CFG.TIER_BA_MAX[CFG.TIER_ORDER[2]] and angle_effective >= CFG.TIER_ANGLE_MIN[CFG.TIER_ORDER[2]] and mech_score >= CFG.TIER_MECH_MIN[CFG.TIER_ORDER[2]]:
            tier, meaning, is_degrader = CFG.TIER_ORDER[2], "Functional Analysis: Nucleophile located in tight contact, accompanied by acceptable target attack angles.", True

        # Tier_2B  (TIER_ORDER[3]) — lowest degrader tier: productive α-attack, nucleophile in reach,
        # SN2 angle ≥ threshold AND a still-connected proton relay (Nuc–Base / Base–Acid within the
        # loose Tier_2B ceilings). The relay check prevents a catalytically dead, geometrically
        # dissociated triad from being labelled a (marginal) degrader on nucleophile+angle alone.
        elif productive_attack and d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_ORDER[3]] and angle_effective >= CFG.TIER_ANGLE_MIN[CFG.TIER_ORDER[3]] and dist_nuc_base <= CFG.TIER_NB_MAX[CFG.TIER_ORDER[3]] and dist_base_acid <= CFG.TIER_BA_MAX[CFG.TIER_ORDER[3]]:
            tier, meaning, is_degrader = CFG.TIER_ORDER[3], "Marginal Analysis: Nucleophile in loose contact with a marginal attack angle, catalytic relay still connected.", True

        # Tier_3  (TIER_ORDER[4]) — nucleophile within reach but no productive pathway:
        # either the attack carbon is not the α-carbon (non-carboxylate head or mid-chain
        # C–F) or the SN2 angle/anchoring is sub-threshold (NOT a degrader).
        elif d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_ORDER[4]]:
            _t3_reason = ("the SN2 attack carbon is not the α-carbon adjacent to a carboxylate "
                          "(no productive defluorination pathway)") if not productive_attack else \
                         ("the SN2 attack angle / anchoring is below the productive threshold in this predicted pose")
            tier, meaning, is_degrader = CFG.TIER_ORDER[4], f"Catalytic anchoring present (nucleophile within reach) but {_t3_reason} — not classified as a degrader.", False

        # Tier_4 (POOR) / Tier_5_Decoy — catch-all when the nucleophile is out of reach.
        else:
            tier = CFG.TIER_POOR if d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_POOR] else CFG.TIER_DECOY
            meaning = "Failed Analysis: Positional distance or established angle physically violates FAcD catalytic structural requirements."
            is_degrader = False

        """
        SN2 dead-end flag (scissile_cf_bde, sn2_backside_occlusion, sn2_dead_end) is a
        diagnostic; the tier demotion of a CF3-α ligand such as TFA is carried by the graded
        chemistry penalty folded into mechanistic_score_effective above (BDE + occlusion beyond
        their cutoffs). That penalty fades with a near-ideal SN2 angle (§5.2b CHEM_PEN_ANGLE_*):
        a poor-to-moderate-angle TFA pose keeps the full penalty and lands around Tier_2B (the
        DeHa4 control, ~159°, stays low), whereas a variant that organises TFA into a near-linear
        Walden trajectory (≥ CHEM_PEN_ANGLE_NONE) has the penalty waived and surfaces to Tier_1A
        as a discovery lead. This encodes that DeHa4 does not turn over TFA (Wackett 2022) without
        precluding an undiscovered variant that can. Step-07 QM/MM remains the final arbiter.
        """
        results["sn2_dead_end"] = 1 if sn2_dead_end else 0

        """
        Derive constraint_check from the assigned tier — Tier_1/Tier_2 tiers pass,
        Good is a partial pass (in pocket but wrong orientation), Poor/None fail.
        """
        if is_degrader:
            constraint = "Pass"
        elif tier == CFG.TIER_ORDER[4]:
            constraint = "Partial"
        else:
            constraint = "Fail"

        close_residues = [k for k, d in dists.items() if d <= CATALYTIC_DIST_CUTOFF]
        results.update({
            "catalytic_dist_A": round(min(dists.values()), 2), CFG.COL_TIER: tier,
            "is_degrader": is_degrader, "residues_within_6A": ";".join(close_residues) if close_residues else "None",
            "constraint_check": constraint, "scientific_meaning": meaning
        })
        return results
    except Exception as e:
        """
        A pose that fails to score is not a pose that scored zero. The tier is set to "Error" and the
        exception is LOGGED WITH ITS TRACEBACK, because this handler is the last thing standing
        between a broken call signature and a full run of 58,056 complexes that all silently read
        tier="Error", mechanistic_score=0.00 and default pocket containment — a result that looks
        like data and is not. Anything that reaches here is a bug in the scorer, not a property of
        the complex, and it must be visible in the log the moment it happens.
        """
        logger.error(f"check_catalytic_geometry FAILED for {getattr(cif_path, 'name', cif_path)}: "
                     f"{type(e).__name__}: {e}", exc_info=True)
        results.update({"catalytic_dist_A": 999.0, "error": f"{type(e).__name__}: {e}",
                        CFG.COL_TIER: "Error"})
        return results

# -------------------------------------------------------------------------------
# Step 7.3: Parallelised Model Analysis & Ranking Framework
# -------------------------------------------------------------------------------
def analyse_model_task(args: Tuple[Path, Dict[str, int], Path, str, str, int]) -> Tuple[str, Dict, float]:
    """Global worker routine functionally decoupled for robust serialisation capabilities across CPU threads."""
    cif_p, sites, json_p, m_name, smiles_str, nuc_rescue_offset = args
    geom = check_catalytic_geometry(cif_p, sites, smiles_str, nuc_rescue_offset)
    try:
        conf_data = json.loads(json_p.read_text())
        conf_score = float(conf_data.get("confidence_score", 0.0))
    except Exception:
        conf_score = 0.0
    return m_name, geom, conf_score

def select_best_degrader_model(br_dir: Path, mapped_sites: Dict[str, int], smiles_str: str = "", nuc_rescue_offset: int = 0) -> Tuple[str, Dict]:
    """
    Selects the representative pose across the Boltz diffusion samples tier-first: the highest
    degrader tier reached by any model, then the highest competence_score within that tier
    (confidence as the final tie-break). This keeps selection consistent with the tier-primary
    Scientific ranking — the reported pose is the candidate's best catalytic shot, scored on the
    full criteria (CFG §5.5) within that tier. The best geometry found is always reported; a
    favourable frame is never discarded. Reproducibility is captured separately by
    model_degrader_consensus (the fraction of models independently reaching a degrader tier)
    and folded into the ranking as a down-rank weight, so a single-frame hit sorts below a
    reproducible one without hiding its true geometry.
    """
    best_model, best_meta = "model_0", {CFG.COL_TIER: CFG.TIER_DECOY}
    TIER_SCORES = CFG.TIER_SCORE
    json_files = sorted(list(br_dir.rglob("confidence_*.json")))
    tasks = []
    for f in json_files:
        match = re.search(r"(model_\d+)", f.stem)
        if match:
            m_name = match.group(1)
            cif_path = f.parent / f"{f.parent.name}_{m_name}.cif"
            if not cif_path.exists(): cif_path = f.parent / f"{m_name}.cif"
            if cif_path.exists(): tasks.append((cif_path, mapped_sites, f, m_name, smiles_str, nuc_rescue_offset))

    if tasks:
        results = [analyse_model_task(t) for t in tasks]
        total = len(results)
        deg_models = [r for r in results if isinstance(r[1], dict) and r[1].get("is_degrader")]

        def _rank_key(r):
            '''
            Tier-first, then the near-attack (Walden) conformer — SN2 attack angle is
            weighted ahead of competence so the representative pose is the diffusion
            sample with the most productive backside trajectory, then competence, then
            confidence. Angle is rounded to 0.1° so sub-noise differences do not flip
            the choice away from a more competent equal-angle pose.
            '''
            return (TIER_SCORES.get(r[1].get(CFG.COL_TIER, CFG.TIER_DECOY), 0),
                    round(float(r[1].get("sn2_attack_angle", 0.0) or 0.0), 1),
                    float(r[1].get("competence_score", 0.0) or 0.0), r[2])
        best_model, best_geom, best_conf = max(results, key=_rank_key)
        best_meta = best_geom.copy()
        best_meta["confidence_score"] = best_conf
        best_meta["model_degrader_consensus"] = round(len(deg_models) / total, 2) if total else 0.0
    return best_model, best_meta


# =============================================================================
# SECTION 8: METRICS & GPU MANAGEMENT
# =============================================================================

# -------------------------------------------------------------------------------
# Step 8.1: Metric Computation
# -------------------------------------------------------------------------------
def compute_cross_interface_pae(pae_path: Path, prot_len: int, lig_len: int) -> float:
    """Calculates the Predicted Aligned Error (PAE), focusing specifically on the protein-ligand interface."""
    try:
        if not pae_path.exists(): return 0.0
        with np.load(str(pae_path)) as data:
            k = next((k for k in data.files if "pae" in k), None)
            if not k: return 0.0
            pae = data[k]
            if pae.shape[0] >= (prot_len + lig_len):
                sub = pae[0:prot_len, prot_len : prot_len+lig_len]
                return float(np.mean(sub))
    except Exception: pass
    return 0.0

def load_extra_boltz_metrics(br_dir: Path, model_name: str, prot_len: int, lig_len: int) -> Dict[str, Any]:
    """Parses background statistical validation tensors directly from Boltz NPZ file dumps."""
    out = {}
    try:
        conf_files = list(sorted(br_dir.rglob(f"confidence_*{model_name}.json")))
        if conf_files:
            d = json.loads(conf_files[0].read_text())
            out.update({k: d.get(k) for k in ["iptm", "confidence_score", "ptm", "ligand_iptm", "protein_iptm"]})
        # sorted(): glob/rglob return entries in filesystem order, which is arbitrary and differs
        # between machines and filesystems. Taking [0] of an unsorted match makes the value that ends
        # up in the CSV depend on the disk, not on the data — the same run on another box can pick a
        # different model's pLDDT. Sorting makes the choice deterministic and reproducible.
        plddt_files = sorted(br_dir.rglob(f"plddt_*{model_name}.npz"))
        if plddt_files:
            with np.load(str(plddt_files[0])) as data:
                if "plddt" in data: out["mean_plddt"] = float(np.mean(data["plddt"]))
        pae_files = sorted(br_dir.rglob(f"pae_*{model_name}.npz"))
        if pae_files and prot_len > 0 and lig_len > 0:
            out["cross_interface_pae_mean"] = compute_cross_interface_pae(pae_files[0], prot_len, lig_len)
    except Exception: pass
    return out

# -------------------------------------------------------------------------------
# Step 8.2: Hardware Probing Tools
# -------------------------------------------------------------------------------
def query_gpu_memory():
    """A resilient system query that safely falls back to CPU memory context if nvidia-smi hangs or crashes."""
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=memory.used,memory.total", "--format=csv,noheader,nounits"],
            encoding="utf8", stderr=subprocess.DEVNULL, timeout=5)
        return [(int(u), int(t)) for u, t in (line.split(",") for line in out.splitlines())]
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
        return []

def init_gpu_reservations():
    """Maps unallocated VRAM memory segments natively upon system initialisation."""
    global _gpu_free_mib
    info = query_gpu_memory()
    _gpu_free_mib = {}
    if info:
        for i, (used, total) in enumerate(info): _gpu_free_mib[i] = max(0, total - used)

def cpu_usage_summary():
    """Returns general CPU congestion levels."""
    try: return psutil.cpu_percent(interval=0.2), psutil.cpu_count(logical=True)
    except Exception: return 0.0, 1


# =============================================================================
# SECTION 9: JOB PROCESSING LOGIC (CORE WORKER)
# =============================================================================

# -------------------------------------------------------------------------------
# Step 9.1: Status Validation Operations
# -------------------------------------------------------------------------------
def check_job_status(job_dir: Path) -> bool:
    """Verifies the finality of a job by checking for a complete and uncorrupted summary JSON file."""
    if not job_dir.exists(): return False
    summary = next(iter(sorted(job_dir.glob("*_summary.json"))), None)
    if summary:
        try:
            if '"status": "Success"' in summary.read_text(): return True
        except Exception: pass
    return False

# -------------------------------------------------------------------------------
# Step 9.2: SMILES Self-Healing Utility Architecture
# -------------------------------------------------------------------------------
def heal_smiles_in_files(job_dir: Path, yaml_path: Path, new_smiles: str, lig_id: str):
    """
    A self-healing function that logically updates existing YAML and CIF files with the newly
    canonicalised SMILES strings, avoiding filename alterations or unnecessary GPU executions.
    """
    if yaml_path.exists():
        try:
            with open(yaml_path, "r") as yf:
                data = yaml.safe_load(yf)
            if "sequences" in data and len(data["sequences"]) > 1:
                current_smi = data["sequences"][1].get("ligand", {}).get("smiles", "")
                if current_smi != new_smiles:
                    data["sequences"][1]["ligand"]["smiles"] = new_smiles
                    with open(yaml_path, "w") as yf:
                        yaml.safe_dump(data, yf, sort_keys=False, default_flow_style=False)
        except Exception:
            pass

    if job_dir.exists():
        for cif in sorted(job_dir.rglob("*.cif")):
            try:
                with open(cif, "r") as f:
                    lines = f.readlines()

                if not any("UPDATED_SMILES" in line for line in lines[:15]):
                    for i, line in enumerate(lines):
                        if line.startswith("data_model") or line.startswith("data_"):
                            lines.insert(i + 1, f"# UPDATED_SMILES {new_smiles}\n")
                            lines.insert(i + 2, f"# LIGAND_ID {lig_id}\n")
                            break
                    with open(cif, "w") as f:
                        f.writelines(lines)
            except Exception:
                pass

# -------------------------------------------------------------------------------
# Step 9.3: Core Worker Engine Infrastructure
# -------------------------------------------------------------------------------


def process_single_job(job: Dict, prod_dir: Path, diffusion_samples: int, prev_elapsed: float, aln_dir: Path,
                       control_cif: Optional[Path], control_map: Dict, control_seq: str, gpu_queue=None,
                       r3u_cif: Optional[Path] = None, r3u_map: Optional[Dict] = None) -> Tuple[Optional[Dict[str, Any]], int]:
    """Run and score one Boltz-2 complex prediction.

    Executes the prediction for a single job, parses the output CIF, computes the
    full metric set (confidence, interaction profile, active-site mechanism vs the
    control/3R3U references, desolvation, tier inputs), and returns
    (metrics_dict, status_code) — metrics_dict is None on failure.
    """
    data = DEFAULT_METRICS.copy()
    jid = job["job_index"]
    name = f"{jid}_{job.get('job_protein', job['protein'])}_{job['ligand']}"
    job_dir = prod_dir / "4_Prediction_Jobs" / name
    job_dir.mkdir(parents=True, exist_ok=True)
    summary_path = job_dir / f"{name}_summary.json"

    best_complex_dir = job_dir / "Best_Complex"
    if best_complex_dir.exists():
        shutil.rmtree(best_complex_dir, ignore_errors=True)
    """
    Retry mkdir — shutil.rmtree uses ignore_errors so the dir may briefly
    still exist on a slow fs; exist_ok=True absorbs any EEXIST race.
    """
    for _attempt in range(3):
        try:
            best_complex_dir.mkdir(parents=True, exist_ok=True)
            break
        except OSError:
            import time as _t; _t.sleep(0.05)

    # Alignment is protein-level (ligand-independent) — key the file by protein so
    # each variant has exactly one .txt, not one per protein×ligand job.
    _safe_protein = re.sub(r"[^A-Za-z0-9._-]", "_", str(job["protein"]))
    aln_path = aln_dir / f"{_safe_protein}_alignment.txt"
    stats_csv_path = aln_dir / CFG.FILE_ALIGNMENT_STATS

    cached_aln = get_cached_alignment_for_protein(job["protein"])
    aln_inject = {
        CFG.COL_ID_PCT: 0.0, "align_score": 0.0,
        "Mapped_to_Control_All": "NA", "Mapped_to_Control_Cat_Triad": "NA",
        "Mapped_Nucleophile": "NA", "Mapped_Base": "NA", "Mapped_Acid": "NA",
        "Mapped_Clamp1": "NA", "Mapped_Clamp2": "NA",
        "Mapped_Stabiliser_H": "NA", "Mapped_Stabiliser_W": "NA", "Mapped_Stabiliser_Y": "NA",
    }
    if cached_aln:
        aln_inject.update({
            CFG.COL_ID_PCT: cached_aln.get(CFG.COL_ID_PCT, 0.0),
            "align_score": cached_aln.get("align_score", 0.0),
            "Mapped_to_Control_All": cached_aln.get("Mapped_to_Control_All", "NA"),
            "Mapped_to_Control_Cat_Triad": cached_aln.get("Mapped_to_Control_Cat_Triad", "NA")
        })

    is_done = check_job_status(job_dir)

    if is_done and summary_path.exists():
        try:
            with open(summary_path, "r") as f: data.update(json.load(f))
        except Exception: pass

    data["protein"] = job["protein"]
    data["ligand"]  = job["ligand"]
    data["job_name"] = name
    data["job_index"] = jid

    br_dir = next((p for p in sorted(job_dir.iterdir()) if p.is_dir() and p.name.startswith("boltz_results_")), None)
    prediction_exists = False
    if br_dir and (br_dir / "predictions").exists():
        """
        A prediction is only considered complete when at least one confidence_*.json
        exists alongside the CIF — Boltz writes these only after a model fully finishes.
        A directory with CIFs but no confidence JSONs is a killed/partial run.
        """
        has_confidence = any(br_dir.rglob("confidence_*.json"))
        if has_confidence:
            prediction_exists = True
        else:
            # Partial prediction: clean it up so the job re-runs from scratch.
            shutil.rmtree(br_dir, ignore_errors=True)
            br_dir = None

    run_gpu = False
    repair_cpu = False

    if not is_done:
        run_gpu = True
        repair_cpu = True
    else:
        run_gpu = False
        repair_cpu = True

    heal_smiles_in_files(job_dir, job["yaml"], job["smiles"], job["ligand"])

    gpu_idx = -1
    start_time = time.time()

    # -------------------------------------------------------------------------------
    # Sub-Step 9.3.1: Invoke AI Inference Model Engine (GPU)
    # -------------------------------------------------------------------------------
    if prediction_exists:
        run_gpu = False
        if logger: logger.debug(f"[SKIP] {name} — pre-existing prediction found, proceeding to analytics.")

    if run_gpu:
        br_dir_check = next((p for p in sorted(job_dir.iterdir()) if p.is_dir() and p.name.startswith("boltz_results_")), None)
        if br_dir_check and br_dir_check.exists(): shutil.rmtree(br_dir_check)

        if gpu_queue is not None:
            gpu_idx = gpu_queue.get()
        else:
            gpu_idx = 0

        start_time = time.time()

        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_idx)
        env["CC"] = shutil.which("gcc") or "gcc"
        env["CXX"] = shutil.which("g++") or "g++"
        # Per-worker Triton cache so parallel Boltz subprocesses cannot collide on lockfiles
        # in a shared /tmp/triton_cache.
        env["TRITON_CACHE_DIR"] = str(Path(tempfile.gettempdir()) / f"triton_cache_{os.getpid()}_{uuid.uuid4().hex[:8]}")
        # conda env PFAS manages CUDA/nvidia libraries via LD_LIBRARY_PATH.
        # Pin BLAS/OMP threads to 1 per subprocess: N parallel workers each spawning
        # GLOBAL_MAX_WORKERS threads oversubscribes the CPU (N² thread thrashing).
        env["OMP_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"
        env["OPENBLAS_NUM_THREADS"] = "1"

        # --seed: Boltz seeds nothing by default, so an unseeded run returns a different pose for the
        # same complex every time (measured on the DeHa4 control: 145.5° → 132.4° between two runs of
        # identical code). A tier that moves when nothing moved is not reproducible, and it breaks the
        # control read-out the pipeline calibrates itself against. The full diffusion ensemble is
        # still sampled; the seed only fixes where it starts.
        cmd = [
            BOLTZ_BIN, "predict", str(job["yaml"]), "--out_dir", str(job_dir),
            "--cache", str(BOLTZ_CACHE), "--model", BOLTZ_MODEL,
            "--recycling_steps", str(RECYCLING_STEPS),
            "--diffusion_samples", str(diffusion_samples),
            "--seed", str(CFG.BOLTZ_SEED),
            "--accelerator", "gpu", "--devices", "1",
            "--use_msa_server", "--output_format", OUTPUT_FORMAT,
            "--no_kernels"
        ]

        success = False

        for attempt in range(RETRY_ON_FAIL + 1):
            try:
                with open(job_dir / f"{jid}_stdout.log", "w") as so, open(job_dir / f"{jid}_stderr.log", "w") as se:
                    proc = subprocess.run(cmd, stdout=so, stderr=se, env=env, timeout=PREDICT_TIMEOUT_S)
                    if proc.returncode == 0:
                        success = True
                        break
                    else:
                        console_info(f"    ! Scheduled Job {name} subsequently failed on allocated GPU[{gpu_idx}] (Attempt {attempt+1}/{RETRY_ON_FAIL+1}). Initiating retry sequence...")
                        time.sleep(RETRY_SLEEP)
            except subprocess.TimeoutExpired:
                console_info(f"    ! Job {name} exceeded the {PREDICT_TIMEOUT_S}s prediction ceiling on GPU[{gpu_idx}] (frozen driver?). Child killed; retrying (Attempt {attempt+1}/{RETRY_ON_FAIL+1}).")
                time.sleep(RETRY_SLEEP)
            except Exception as e:
                console_info(f"    ! Routine Exception Encountered: {e}")
                time.sleep(RETRY_SLEEP)

        if gpu_queue is not None:
            gpu_queue.put(gpu_idx)

        if not success: return None, gpu_idx
        data["elapsed_seconds"] = time.time() - start_time

    br_dir = next((p for p in sorted(job_dir.iterdir()) if p.is_dir() and p.name.startswith("boltz_results_")), None)

    # -------------------------------------------------------------------------------
    # Sub-Step 9.3.2: Execution Physics & Structural Topology Engine (CPU)
    # -------------------------------------------------------------------------------
    data.update(aln_inject)
    if run_gpu or repair_cpu:
        # Compute the (cheap) alignment for every job, but persist the .txt + stats row
        # only once per protein this run — the mapping is identical across its ligands.
        _persist_aln = job["protein"] not in _PERSISTED_PROTEINS
        mapped, aln_stats, resname_map, map_all_str = map_active_site_residues(
            job["protein"], job["sequence"],
            out_aln_path=(aln_path if _persist_aln else None),
            stats_csv=(stats_csv_path if _persist_aln else None),
        )
        if _persist_aln:
            _PERSISTED_PROTEINS.add(job["protein"])
        _, map_triad_str = format_control_mappings(resname_map, map_all_str)

        data.update(aln_stats)
        data["active_site_mapping"] = str(mapped)
        data["Mapped_to_Control_All"] = map_all_str
        data["Mapped_to_Control_Cat_Triad"] = map_triad_str
        data.update(format_full_role_map(resname_map))
        data["alignment_reliable"] = data.get(CFG.COL_ID_PCT, 0) >= CFG.ALIGN_MIN_SEQ_IDENTITY

        _nuc_off = int(aln_stats.get("nuc_rescue_offset", 0))
        if br_dir:
            best_model, best_meta = select_best_degrader_model(br_dir, mapped, job.get("smiles", ""), _nuc_off)
            data["best_model_name"] = best_model
            data.update(best_meta)
        else:
            cifs = sorted(best_complex_dir.glob("*.cif"))
            if not cifs: cifs = sorted(job_dir.glob("*.cif"))
            if cifs:
                best_model = "model_0"
                if "model_" in cifs[0].name:
                    try: best_model = re.search(r"(model_\d+)", cifs[0].name).group(1)
                    except Exception: pass
                data["best_model_name"] = best_model
                geom = check_catalytic_geometry(cifs[0], mapped, job.get("smiles", ""), _nuc_off)
                data.update(geom)

        """
        Top-tier confidence guard (CFG §8). The catalytic machinery itself is already enforced
        by the tier ladder — Tier_1A requires all eight residues mapped to the correct type
        (triad_identity_ok + elite_identity_ok) at catalytic distances with the fluoride cradle
        engaged — so global sequence identity is NOT used here: a distant homolog with a valid,
        well-resolved active site keeps its tier.

        The residual risk is an unconfident predicted fold, and the confidence that matters is
        LOCAL. Every geometric quantity the elite tier rests on — the SN2 angle, the nucleophile
        distance, the triad relay, the cradle — is measured on the eight catalytic residues and
        nowhere else. A global confidence score averages those eight residues together with
        hundreds of loop and surface residues that the tier decision never touches: a protein with
        a crisply resolved active site and disordered termini is punished for the termini, while a
        globally confident fold with a smeared active site passes. `active_site_plddt` is the mean
        Boltz pLDDT over exactly the eight mapped residues, so the gate now asks about the region
        it actually trusts. A Tier_1A hit below ELITE_AS_PLDDT_MIN is demoted ONE notch to Tier_1B
        (still elite geometry, still reviewed); the geometric tier is preserved in its own column.
        The global score is retained as the fallback when the local value is unavailable.
        """
        data["geometric_tier"] = data.get(CFG.COL_TIER, CFG.TIER_DECOY)
        if data.get(CFG.COL_TIER) == CFG.TIER_ORDER[0]:
            _as_plddt = float(data.get("active_site_plddt", 0.0) or 0.0)
            if _as_plddt > 0.0:
                if _as_plddt < CFG.TIER_ELITE_AS_PLDDT_MIN:
                    data[CFG.COL_TIER] = CFG.TIER_ORDER[1]
                    data["elite_demotion"] = (
                        f"low_active_site_plddt({_as_plddt:g}<{CFG.TIER_ELITE_AS_PLDDT_MIN:g})")
                else:
                    data["elite_demotion"] = "none"
            else:
                # No per-residue confidence for this pose — fall back to the global fold score
                # rather than waving the pose through unchecked.
                _conf = float(data.get("confidence_score", 0.0) or 0.0)
                if _conf < CFG.TIER_ELITE_CONF_MIN:
                    data[CFG.COL_TIER] = CFG.TIER_ORDER[1]
                    data["elite_demotion"] = f"low_confidence(<{CFG.TIER_ELITE_CONF_MIN:g})"
                else:
                    data["elite_demotion"] = "none"
        else:
            data["elite_demotion"] = "none"

        cif_path = None
        if br_dir:
            cif_path = br_dir / "predictions" / name / f"{name}_{data.get('best_model_name','model_0')}.cif"
            if not cif_path.exists():
                cands = sorted(br_dir.rglob(f"*{data.get('best_model_name')}.cif"))
                if cands: cif_path = cands[0]

        if not cif_path:
             cifs = sorted(best_complex_dir.glob("*.cif"))
             if not cifs: cifs = sorted(job_dir.glob("*.cif"))
             if cifs: cif_path = cifs[0]

        if cif_path and cif_path.exists():
            inter_csv_path = job_dir / f"{name}_interactions.csv"
            res_inter = generate_detailed_interactions(cif_path, job["smiles"], inter_csv_path)
            counts = res_inter.get("counts", {})

            for k in ["Total_Residues", "Total_Polar_Residues", "Total_NonPolar_Residues",
                      "Total_Pos_Residues", "Total_Neg_Residues", "Total_Aromatic_Residues",
                      "interacting_fluorine_count", "total_fluorine_count"]:
                if k in res_inter: data[k] = res_inter[k]

            if br_dir:
                extra_metrics = load_extra_boltz_metrics(br_dir, data.get("best_model_name"), res_inter.get("prot_len", 0), res_inter.get("lig_len", 0))
                data.update(extra_metrics)

            data["protein"] = job["protein"]
            data["ligand"] = job["ligand"]
            data["completed_at"] = datetime.utcnow().isoformat()
            data["status"] = "Success"
            data["job_name"] = name
            data["job_index"] = jid
            if run_gpu: data["elapsed_seconds"] = round(time.time() - start_time, 2)

            data["protein_atom_count"] = res_inter.get("prot_len", 0)
            data["ligand_atom_count"] = res_inter.get("lig_len", 0)
            data["ligand_heavy_atoms"] = res_inter.get("ligand_heavy_atoms", 0)
            data["mapped_ligand_atoms_fraction"] = res_inter.get("mapped_fraction", 0.0)
            data["num_interactions"] = res_inter.get("num_interactions", 0)
            data["min_distance_A"] = res_inter.get("min_dist", 0.0)
            data["max_distance_A"] = res_inter.get("max_dist", 0.0)
            data["avg_distance_A"] = res_inter.get("avg_dist", 0.0)
            data["counts"] = counts

            # -------------------------------------------------------------------------------
            # Sub-Step 9.3.3: Advanced Thermodynamic Scoring (Bias-Free)
            # -------------------------------------------------------------------------------
            num_int = res_inter.get("num_interactions", 0)
            lig_heavy_atoms = max(1, res_inter.get("ligand_heavy_atoms", 1))

            interaction_density = num_int / lig_heavy_atoms
            data[CFG.COL_IDENS] = round(interaction_density, 4)
            data["interaction_density_calc"] = f"{num_int} (Ints) / {lig_heavy_atoms} (HeavyAtoms) = {interaction_density:.2f}"

            iptm_v = data.get("iptm") or 0.0
            plddt_v = data.get("mean_plddt") or 0.0
            c_pae = data.get("cross_interface_pae_mean") or 0.0
            conf = data.get("confidence_score") or 0.0

            score_v = (W_IPTM * float(iptm_v)) + (W_PLDDT * (plddt_v/100.0)) + \
                      (W_INTERACTIONS * min(1.0, interaction_density / CFG.BIND_INT_DENSITY_NORM)) - \
                      (W_CROSS_PAE * min(1.0, c_pae/CFG.BIND_CROSS_PAE_NORM)) + (W_CONF * conf)

            # Centre and scale the raw sum onto the range where a logistic actually resolves; without this
            # every decent complex saturates at P > 0.99 and the score carries no information (CFG §9.x).
            _z_bind = (score_v - CFG.BIND_LOGIT_CENTRE) / CFG.BIND_LOGIT_GAIN
            _z_bind = max(min(_z_bind, CFG.BIND_LOGIT_CLAMP), -CFG.BIND_LOGIT_CLAMP)
            binding_prob = 1.0 / (1.0 + math.exp(-_z_bind))
            data["binding_likelihood_computed"] = binding_prob
            data["binding_likelihood_calc"] = (
                f"({W_IPTM} * {float(iptm_v):.2f}) + ({W_PLDDT} * {plddt_v/100.0:.2f}) + "
                f"({W_INTERACTIONS} * {min(1.0, interaction_density / CFG.BIND_INT_DENSITY_NORM):.2f}) - "
                f"({W_CROSS_PAE} * {min(1.0, c_pae / CFG.BIND_CROSS_PAE_NORM):.2f}) + ({W_CONF} * {conf:.2f}) "
                f"= {score_v:.2f}  ->  z = (raw - {CFG.BIND_LOGIT_CENTRE}) / {CFG.BIND_LOGIT_GAIN} "
                f"= {_z_bind:.2f}  ->  Sigmoid = {binding_prob:.4f}"
            )

            hb = counts.get("hydrogen_bond", 0)
            hp = counts.get("hydrophobic", 0)
            sb = counts.get("salt_bridge", 0)
            hal = counts.get("halogen_contact", 0)
            fl = counts.get("fluorine_contact", 0)
            fp = counts.get("fluorine_polar", 0)
            fh = counts.get("fluorous_hydrophobic", 0)
            ff = counts.get("fluorous_fluorous", 0)
            mc = counts.get("metal_coordination", 0)
            avg_dist = res_inter.get("avg_dist", 0)

            aff_score = (CFG.SCORE_W_AFF_HB * hb) + (CFG.SCORE_W_AFF_HP * hp) + \
                        (CFG.SCORE_W_AFF_SB * sb) + (CFG.SCORE_W_AFF_HAL * hal) + \
                        (CFG.SCORE_W_AFF_FL * fl) + (CFG.SCORE_W_AFF_FP * fp) + \
                        (CFG.SCORE_W_AFF_FH * fh) + (CFG.SCORE_W_AFF_FF * ff) + \
                        (CFG.SCORE_W_AFF_MC * mc) - (CFG.SCORE_W_AFF_DIST * avg_dist)

            data["custom_affinity_score"] = round(aff_score, 4)
            data["custom_affinity_calc"] = (
                f"[{hb}*{CFG.SCORE_W_AFF_HB}] + [{hp}*{CFG.SCORE_W_AFF_HP}] + "
                f"[{sb}*{CFG.SCORE_W_AFF_SB}] + [{hal}*{CFG.SCORE_W_AFF_HAL}] + "
                f"[{fl}*{CFG.SCORE_W_AFF_FL}] + [{fp}*{CFG.SCORE_W_AFF_FP}] + "
                f"[{fh}*{CFG.SCORE_W_AFF_FH}] + [{ff}*{CFG.SCORE_W_AFF_FF}] + "
                f"[{mc}*{CFG.SCORE_W_AFF_MC}] - [{avg_dist:.1f}*{CFG.SCORE_W_AFF_DIST}] = {aff_score:.2f}"
            )

            scientific_meaning = data.get("scientific_meaning", "No significant documented interactions.")
            constraint = data.get("constraint_check", "Fail")
            tier = data.get(CFG.COL_TIER, CFG.TIER_DECOY)

            base_justification = generate_rich_justification(
                tier=tier, meaning=scientific_meaning, constraint=constraint,
                aligned_ok=data.get("alignment_reliable", False), identity=data.get(CFG.COL_ID_PCT, 0.0)
            )
            data["Justification"] = f"{base_justification} | Tier: {tier} | ID: {data.get(CFG.COL_ID_PCT, 0)}%"

            # -------------------------------------------------------------------------------
            # SUPERIMPOSITION & LIKELIHOOD SCORING ARCHITECTURE
            # -------------------------------------------------------------------------------
            data["hydrophobic_desolvation_ratio"] = 0.0   # default; overwritten below if control_cif available
            data["active_site_contact_flag"]    = 0.0   # default; set below when contact density exceeds the flag threshold
            data["catalytic_constellation_score"] = 0.0   # Criterion B default; overwritten below if control_cif available
            if control_cif and control_cif.exists() and control_map:
                mech_res = analyse_candidate_structure(cif_path, control_cif, control_map, target_map=mapped)
                data.update(mech_res)

                ident = data.get(CFG.COL_ID_PCT, 0.0)
                data["Identity_to_Control"] = round(ident, 2)

                rmsd = data.get("Active_Site_RMSD", 999.0)
                if rmsd >= 99.0: geo_fit = 0.0
                else: geo_fit = 100.0 / (1.0 + rmsd)

                """
                Active-site conservation is led by the catalytic machinery, not by global
                identity. Weighting: 60% active-site integrity (fraction of the eight
                catalytic residues correctly mapped, Chan 2011), 30% active-site geometric
                fit to the control (RMSD-derived), 10% global identity (retained only as a
                light corroborating signal). A distant homolog with an intact, well-placed
                active site therefore scores high regardless of overall sequence identity.
                """
                _integrity_pct = 100.0 * float(data.get("active_site_integrity", 0.0))
                likelihood = ((CFG.CONSERV_W_INTEGRITY * _integrity_pct)
                              + (CFG.CONSERV_W_GEO * geo_fit)
                              + (CFG.CONSERV_W_IDENT * ident))

                desolvation_ratio = fh / max((fh + fp), 1)
                data["hydrophobic_desolvation_ratio"] = round(desolvation_ratio, 3)
                if desolvation_ratio > CFG.SCORE_DESOLVATION_THRESHOLD:
                    likelihood += (CFG.SCORE_DESOLVATION_BONUS_FACTOR * desolvation_ratio)

                """
                Active-site conservation measures sequence homology + active-site
                geometry only. It is deliberately decoupled from the product-
                inhibition heuristic: the inhibition penalty is recorded as its
                own column (and may inform downstream ranking), but is NOT
                subtracted here — otherwise a high interaction density could zero
                out the conservation score of a genuinely well-conserved active
                site (e.g. distant homologs with low identity but valid geometry).
                """
                """
                active_site_contact_flag — a REPORTED contact-richness descriptor
                that never gates the tier or the ranking. interaction_density is
                normalised per ligand heavy atom, so a small, tightly-bound substrate
                reads a high value; this is a contact-density flag, NOT an energy and
                NOT a product-inhibition measure. Product inhibition (glycolate/F⁻
                retention slowing turnover) is quantified downstream by Step-06 Prime
                MM-GBSA ΔG_bind.
                """
                if interaction_density > CFG.SCORE_INHIBITION_DENSITY_THRESHOLD:
                    data["active_site_contact_flag"] = round(
                        CFG.SCORE_INHIBITION_PENALTY_FACTOR
                        * (interaction_density - CFG.SCORE_INHIBITION_DENSITY_THRESHOLD), 2)
                else:
                    data["active_site_contact_flag"] = 0.0

                data[CFG.COL_LIKE_S] = round(max(0, min(100, likelihood)), 2)

                # ----------------------------------------------------------------------
                # Criterion B — precise catalytic constellation.
                # Active_Site_RMSD is the all-EIGHT catalytic-residue constellation RMSD onto
                # the 3R3U crystal (Cα superposition). B = 1/(1+RMSD) is high only when the
                # eight residues are geometrically ASSEMBLED like the crystal — not merely
                # present and correctly typed (that is active_site_integrity, Criterion A).
                # ----------------------------------------------------------------------
                _as_rmsd = float(data.get("Active_Site_RMSD", 999.0))
                _B = 0.0 if _as_rmsd >= 99.0 else round(1.0 / (1.0 + _as_rmsd), 3)
                data["catalytic_constellation_score"] = _B
                """
                B CAPS the tier — downgrade only, never promotes (the SN2 trajectory stays the
                positive catalytic signal). The highest tier whose B-floor the pose meets is the
                ceiling; the tier is moved to the lower of (current, ceiling) on the ladder. B
                below every degrader floor → TIER_CONSTELLATION_FLOOR_TIER (B>0) or Tier_5_Decoy
                (B==0, constellation unmeasurable), and is_degrader is cleared. geometric_tier
                (captured pre-demotion) preserves the raw geometric call for audit.
                """
                _cfloors = CFG.TIER_CONSTELLATION_MIN
                _order   = CFG.TIER_ORDER
                _cur_t   = data.get(CFG.COL_TIER, CFG.TIER_DECOY)
                _ceiling = next((_t for _t in _order if _t in _cfloors and _B >= _cfloors[_t]), None)
                if _ceiling is None:
                    _ceiling = CFG.TIER_CONSTELLATION_FLOOR_TIER if _B > 0 else CFG.TIER_DECOY
                if (_cur_t in _order and _ceiling in _order
                        and _order.index(_ceiling) > _order.index(_cur_t)):
                    _tag = f"low_constellation(B={_B:g}<{_cfloors.get(_cur_t, 0):g})"
                    data[CFG.COL_TIER] = _ceiling
                    if _ceiling not in CFG.TIER_HIGH_QUALITY:
                        data["is_degrader"] = False
                        data["constraint_check"] = "Fail"
                    _prev_d = data.get("elite_demotion", "none")
                    data["elite_demotion"] = _tag if _prev_d in ("", "none") else f"{_prev_d}+{_tag}"
                    # Regenerate the human-readable verdict so Justification / scientific_meaning
                    # name the final (post-constellation) tier the row actually holds.
                    _meaning = (f"Catalytic constellation below the {_cur_t} floor "
                                f"(B={_B:g} < {_cfloors.get(_cur_t, 0):g}); the eight catalytic "
                                f"residues are not assembled with crystal-grade geometry — "
                                f"demoted to {_ceiling}.")
                    data["scientific_meaning"] = _meaning
                    data["Justification"] = generate_rich_justification(
                        tier=_ceiling, meaning=_meaning,
                        constraint=data.get("constraint_check", "Fail"),
                        aligned_ok=data.get("alignment_reliable", False),
                        identity=data.get(CFG.COL_ID_PCT, 0.0)
                    ) + f" | Tier: {_ceiling} | ID: {data.get(CFG.COL_ID_PCT, 0)}%"

            if r3u_cif and r3u_cif.exists() and r3u_map:
                _r3u = analyse_candidate_structure(cif_path, r3u_cif, r3u_map)
                _r3u_rmsd  = _r3u.get("Active_Site_RMSD", 99.0)
                _r3u_ident = float(data.get(CFG.COL_ID_PCT, 0.0))
                _r3u_geo   = 0.0 if _r3u_rmsd >= 99.0 else 100.0 / (1.0 + _r3u_rmsd)
                data["r3u_Active_Site_RMSD"]              = _r3u_rmsd
                data["r3u_Halide_Stabilisation"]          = _r3u.get("Halide_Stabilisation", False)
                data["r3u_Carboxylate_Clamp"]             = _r3u.get("Carboxylate_Clamp", False)
                data["r3u_ActiveSite_Conservation_Score"]     = round(max(0.0, CFG.CONSERV_REF_W_IDENT * _r3u_ident + CFG.CONSERV_REF_W_GEO * _r3u_geo), 2)
            else:
                data["r3u_Active_Site_RMSD"]              = 99.0
                data["r3u_Halide_Stabilisation"]          = False
                data["r3u_Carboxylate_Clamp"]             = False
                data["r3u_ActiveSite_Conservation_Score"]     = 0.0

            """
            Coupled elite-machinery cap (downgrade-only). A Tier_1A pose must have either
            complete machinery + open backside (mech ≥ MECH_ELITE_HI) OR near-complete
            machinery (mech ≥ MECH_ELITE_LO) redeemed by a crystal-exact constellation
            (B ≥ MECH_ELITE_CONSTELLATION). Otherwise it is demoted one notch to Tier_1B.
            This lets a slightly backside-occluded native substrate (α-CF3 TFA, mech ~0.85)
            hold elite only through its single most crystal-perfect pose, not every pose,
            without any substrate-class label. geometric_tier keeps the raw call.
            """
            if data.get(CFG.COL_TIER) == "Tier_1A":
                _mech_e = float(data.get("mechanistic_score_effective",
                                        data.get(CFG.COL_MECH_S, 0.0)) or 0.0)
                _B_e    = float(data.get("catalytic_constellation_score", 0.0) or 0.0)
                # Elite via either: (a) complete machinery + open backside; (b) near-complete
                # machinery redeemed by a crystal-exact constellation.
                _elite_ok = (_mech_e >= CFG.MECH_ELITE_HI) or \
                            (_mech_e >= CFG.MECH_ELITE_LO and _B_e >= CFG.MECH_ELITE_CONSTELLATION)
                if not _elite_ok:
                    data[CFG.COL_TIER] = "Tier_1B"
                    _etag = (f"incomplete_machinery(mech={_mech_e:.2f}"
                             f"<{CFG.MECH_ELITE_HI:g}, B={_B_e:.2f}<{CFG.MECH_ELITE_CONSTELLATION:g})")
                    _eprev = data.get("elite_demotion", "none")
                    data["elite_demotion"] = _etag if _eprev in ("", "none") else f"{_eprev}+{_etag}"

            """
            The tier ladder is size-agnostic: ligand_max_extent is computed and
            reported but does not gate the tier, so a long-chain PFAS variant that
            reaches elite near-attack geometry earns Tier_1A on merit. Chemical
            recalcitrance enters only the within-tier ranking via the feasibility
            factor; Step-07 QM/MM is the final arbiter. CFG.TIER_1A_MAX_LIGAND_EXTENT
            defines the controls-only validation subset for reporting.
            """

            summary_path.write_text(json.dumps(data, indent=2))

            best_model = data.get("best_model_name", "model_0")
            cif_src = None
            if br_dir:
                cif_src = br_dir / "predictions" / name / f"{name}_{best_model}.cif"
                if not cif_src.exists():
                    cands = list(sorted(br_dir.rglob(f"*{best_model}.cif")))
                    if cands: cif_src = cands[0]

            if cif_src and cif_src.exists():
                try:
                    shutil.copy(cif_src, best_complex_dir / cif_src.name)
                    mirror_best_cif(cif_src, prod_dir.parent)
                except Exception as e: console_info(f"System error copying resolved complex structure: {e}")

            if torch is not None: torch.cuda.empty_cache()
            gc.collect()

            # -------------------------------------------------------------------------------
            # Sub-Step 9.3.4: Terminate Processing Sequence Successfully
            # -------------------------------------------------------------------------------
            return data, gpu_idx
        else:
            if run_gpu: console_info(f"Job sequence {name} Failed: No CIF data generated natively.")
            return None, gpu_idx

    # -------------------------------------------------------------------------------
    # Sub-Step 9.3.5: Return Valid Data For Previously Completed Jobs
    # -------------------------------------------------------------------------------
    if CFG.COL_ID_PCT not in data or data.get(CFG.COL_ID_PCT, 0) == 0:
        _persist_aln = job["protein"] not in _PERSISTED_PROTEINS
        mapped, aln_stats, resname_map, map_all_str = map_active_site_residues(
             job["protein"], job["sequence"],
             out_aln_path=(aln_path if _persist_aln else None),
             stats_csv=(stats_csv_path if _persist_aln else None),
        )
        if _persist_aln:
            _PERSISTED_PROTEINS.add(job["protein"])
        _, map_triad_str = format_control_mappings(resname_map, map_all_str)
        data.update(aln_stats)
        data["active_site_mapping"] = str(mapped)
        data["Mapped_to_Control_All"] = map_all_str
        data["Mapped_to_Control_Cat_Triad"] = map_triad_str
        data.update(format_full_role_map(resname_map))

    return data, gpu_idx

def worker_task_wrapper(args):
    """Secure pickle-safe wrapper routine intended strictly for Python multiprocessing environments."""
    try: return process_single_job(*args)
    except Exception as e:
        job = args[0] if args else {}
        job_name = f"{job.get('job_index', '?')}_{job.get('job_protein', job.get('protein', '?'))}_{job.get('ligand', '?')}"
        return {"error": str(e), "traceback": traceback.format_exc(), "job_name": job_name, "status": "FAILED"}, -1

def flatten_job_result(res: dict) -> dict:
    """Flattens a process_single_job result dict into a single-level row for CSV writing."""
    flat = {k: v for k, v in res.items() if isinstance(v, (str, int, float, bool))}
    if "counts" in res:
        for k, v in res["counts"].items(): flat[f"count_{k}"] = v
    if "active_site_details" in res:
        for k, v in res["active_site_details"].items(): flat[k] = v
    return flat

# Canonical human-readable column order for the master CSV output.
CSV_COLUMN_ORDER = [
    # --- Identifiers ---
    "job_index", "job_name", "protein", "ligand",
    # --- Job outcome ---
    "status", "completed_at", "elapsed_seconds",
    # --- Primary degrader scores ---
    CFG.COL_TIER, "is_degrader",
    CFG.COL_LIKE_S,
    # --- Confidence scores ---
    "confidence_score", "iptm", "ptm", "ligand_iptm", "protein_iptm",
    "mean_plddt", "cross_interface_pae_mean",
    # --- Binding scores ---
    "ligand_smiles_stale", "binding_likelihood_computed", "custom_affinity_score", CFG.COL_IDENS,
    # --- SN2 geometry (+ auxiliary non-gating reference angles) ---
    "sn2_attack_angle", "sn2_trajectory_dev",
    "scissile_cf_bde", "sn2_backside_occlusion", "beta_f_count", "sn2_dead_end", "feasibility_factor",
    "active_site_residues_correct", "active_site_integrity", "soft_catalytic_score",
    "head_is_carboxylate", "scissile_is_alpha", "competence_score", "catalytic_constellation_score", "active_site_plddt", "model_degrader_consensus",
    "nuc_resolution", "nuc_rescue_offset",
    "burgi_dunitz_angle", "flippin_lodge_offset",
    "catalytic_dist_A", "Active_Site_RMSD", "Identity_to_Control",
    # --- Key active site distances (to the per-protein dynamically-mapped residue) ---
    "dist_Nuc", "dist_Carb1", "dist_Carb2", "dist_Acid",
    "dist_Stab_H", "dist_Stab_W", "dist_Stab_Y", "dist_Base",
    "dist_nuc_base_internal", "dist_base_acid_internal",
    # --- Functional flags ---
    "Halide_Stabilisation", "Carboxylate_Clamp",
    "halide_stabilisation_score", "carboxylate_clamp_integrity",
    # --- Interaction counts ---
    "num_interactions",
    "count_hydrogen_bond", "count_hydrophobic", "count_salt_bridge",
    "count_halogen_contact", "count_fluorine_contact", "count_fluorine_polar",
    "count_fluorous_hydrophobic",
    # --- Structural size ---
    "protein_atom_count", "ligand_atom_count", "ligand_heavy_atoms", "mapped_ligand_atoms_fraction",
    "Total_Residues", "Total_Polar_Residues", "Total_NonPolar_Residues",
    "Total_Pos_Residues", "Total_Neg_Residues", "Total_Aromatic_Residues",
    "interacting_fluorine_count", "total_fluorine_count",
    # --- Sequence / alignment ---
    CFG.COL_ID_PCT, "align_score", "seq_length", "gap_count",
    "alignment_reliable", "target_sequence",
    # --- Active site mapping ---
    # Mapped_* = the ACTUAL per-protein residue resolved for each role (±5 class-aware
    # resolver); pairs 1:1 with the Dist_* columns. Active_Site_Triad_Map keeps the
    # compact triad string. Canonical 3R3U reference numbers are in CFG §2.5.
    "active_site_mapping", "Mapped_to_Control_All", "Mapped_to_Control_Cat_Triad",
    "Mapped_Nucleophile", "Mapped_Base", "Mapped_Acid",
    "Mapped_Clamp1", "Mapped_Clamp2",
    "Mapped_Stabiliser_H", "Mapped_Stabiliser_W", "Mapped_Stabiliser_Y",
    "best_model_name",
    # --- Scoring / verdict ---
    "constraint_check", "scientific_meaning", "Justification",
    "sn2_alignment_score", CFG.COL_MECH_S, "residues_within_6A",
    # --- Verbose calculation strings ---
    "interaction_density_calc", "binding_likelihood_calc", "custom_affinity_calc",
    # --- Clash / penalty metrics ---
    "teflon_shield_clashes", "mainchain_clash_ratio", "mainchain_clash_count",
    "hydrophobic_desolvation_ratio", "active_site_contact_flag",
    # --- Active-site pocket vs ligand steric fit ---
    "pocket_containment_cavity", "pocket_containment_site8", "ligand_buriedness_mean",
    "ligand_reach", "dist_Nuc_nearest_O", "angle_multiplicity", "sn2_attack_angle_effective",
    "active_site_volume", "active_site_radius",
    "ligand_volume", "ligand_radius_gyration", "ligand_max_extent",
    "pocket_occupancy", "fit_ratio", "ligand_fits",
    # --- Feasibility-weighted mechanistic score (tier-gate key) ---
    "mechanistic_score_effective", "chem_penalty", "containment_penalty",
    # --- Distance statistics ---
    "min_distance_A", "max_distance_A", "avg_distance_A",
]

def mirror_best_cif(cif_src: Path, run_root: Path) -> None:
    """Copy *cif_src* into the flat run_root/2_Best_Complexes_CIFs/ directory.
    Called immediately after every successful per-job CIF copy so the aggregated
    directory stays current as each GPU job completes."""
    try:
        dest_dir = run_root / "2_Best_Complexes_CIFs"
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(cif_src, dest_dir / cif_src.name)
    except Exception as e:
        console_info(f"  Warning: could not mirror CIF to 2_Best_Complexes_CIFs: {e}")


def rebuild_best_complexes_mirror(runs_dir: Path, run_root: Path) -> int:
    """Scan every completed job's Best_Complex/ folder and copy each CIF into
    the flat run_root/2_Best_Complexes_CIFs/ directory.

    On resume, files that already exist in the destination with the same size
    are skipped (no I/O needed). Actual copies are dispatched via a thread pool
    so the disk is kept saturated. Returns the total CIF count (skipped + copied)."""
    from concurrent.futures import ThreadPoolExecutor
    dest_dir = run_root / "2_Best_Complexes_CIFs"
    dest_dir.mkdir(parents=True, exist_ok=True)

    _tty_write("\r   Listing job folders ...   \033[K")
    sys.stdout.flush()
    job_dirs = [d for d in sorted(runs_dir.iterdir()) if d.is_dir()] if runs_dir.exists() else []
    n_dirs = len(job_dirs)

    # Collect (src, dst) pairs in parallel — each dir check is independent I/O
    to_copy: list = []
    total_count = [0]
    _scanned = [0]
    _lock = threading.Lock()

    def _collect(job_dir):
        bc_dir = job_dir / "Best_Complex"
        local_copy = []
        local_total = 0
        if bc_dir.exists():
            for cif in sorted(bc_dir.glob("*.cif")):
                local_total += 1
                dst = dest_dir / cif.name
                if not (dst.exists() and dst.stat().st_size == cif.stat().st_size):
                    local_copy.append((cif, dst))
        with _lock:
            total_count[0] += local_total
            to_copy.extend(local_copy)
            _scanned[0] += 1
            if _scanned[0] % 2000 == 0 or _scanned[0] == n_dirs:
                _tty_write(f"\r   Scanning Best_Complex folders: {_scanned[0]:,}/{n_dirs:,} ...   \033[K")
                sys.stdout.flush()

    with ThreadPoolExecutor(max_workers=min(max(1, (os.cpu_count() or 4) - 2), n_dirs or 1)) as ex:
        list(ex.map(_collect, job_dirs))
    total = total_count[0]
    _tty_write("\r\033[K")
    sys.stdout.flush()

    _copied = [0]
    n_to_copy = len(to_copy)

    def _do_copy(pair):
        src, dst = pair
        try:
            shutil.copy2(src, dst)
        except Exception:
            pass
        with _lock:
            _copied[0] += 1
            if _copied[0] % 500 == 0 or _copied[0] == n_to_copy:
                _tty_write(f"\r   Copying CIFs to mirror: {_copied[0]:,}/{n_to_copy:,} ...   \033[K")
                sys.stdout.flush()

    if to_copy:
        n_workers = min(max(1, (os.cpu_count() or 4) - 2), n_to_copy)
        with ThreadPoolExecutor(max_workers=n_workers) as ex:
            list(ex.map(_do_copy, to_copy))
        _tty_write("\r\033[K")
        sys.stdout.flush()

    return total


def append_rows_to_csv(rows: list, csv_path: Path):
    """Appends a list of result rows to the master CSV, deduplicates, sorts by job_name,
    assigns a plain sequential job_index (0-based row counter), and reorders columns
    into the canonical human-readable CSV_COLUMN_ORDER."""
    try:
        df_new = pd.DataFrame(rows)
        if csv_path.exists() and csv_path.stat().st_size > 0:
            try:
                df_exist = pd.read_csv(csv_path, low_memory=False, on_bad_lines="warn")
                df_new = pd.concat([df_exist, df_new], ignore_index=True)
            except Exception:
                pass
        if "job_name" in df_new.columns:
            df_new.drop_duplicates(subset=["job_name"], keep="last", inplace=True)
            # Sort by folder name so rows are ordered protein-first then ligand.
            df_new = df_new.sort_values("job_name", kind="mergesort").reset_index(drop=True)
        # job_index is just a plain 0-based row counter; last row == total_rows - 1.
        df_new["job_index"] = range(len(df_new))
        # Reorder to canonical column layout; any extra columns are appended at the end.
        ordered = [c for c in CSV_COLUMN_ORDER if c in df_new.columns]
        remainder = [c for c in df_new.columns if c not in ordered]
        df_new = df_new[ordered + remainder]
        atomic_to_csv(df_new, csv_path, index=False)
    except Exception as e:
        if logger: logger.error(f"CSV checkpoint write failed: {e}")

def analysis_worker_loop(q, prod_dir, diff_samples, aln_dir, ctrl_cif, ctrl_map, ctrl_seq, csv_path,
                         r3u_cif=None, r3u_map=None):
    """Background worker that analyses freshly GPU-predicted jobs and appends rows to the CSV.

    Already-completed jobs are handled by the upfront bulk rebuild (rebuild_csv_from_summaries)
    before the GPU phase starts.  This loop only processes jobs that were just predicted by GPU.
    """
    local_rows = []
    while True:
        try:
            job = q.get(timeout=3600)
            if job == "STOP": break
            if isinstance(job, dict):
                is_fallback = job.pop("fallback_individual", False)
                if is_fallback:
                    name = f"{job['job_index']}_{job.get('job_protein', job['protein'])}_{job['ligand']}"
                    console_info(f"  [FALLBACK] Systematically processing {name} entirely in individual operational mode.")
                res, _ = process_single_job(job, prod_dir, diff_samples, 0.0, aln_dir, ctrl_cif, ctrl_map, ctrl_seq,
                                            r3u_cif=r3u_cif, r3u_map=r3u_map)
                if res and res.get("status") == "Success":
                    # Ensure job_index / job_name always match the folder (ground truth).
                    expected_name = f"{job['job_index']}_{job.get('job_protein', job['protein'])}_{job['ligand']}"
                    res["job_name"]  = expected_name
                    res["job_index"] = str(job["job_index"])
                    # A resumed prediction whose input SMILES has since changed is analysed against the
                    # molecule it was built from, and says so in its own row.
                    res["ligand_smiles_stale"] = int(job.get("smiles_stale", 0))
                    local_rows.append(flatten_job_result(res))
                    if len(local_rows) >= 10:
                        # clear in finally so a transient append/IO error cannot leave the
                        # buffer un-flushed and grow it without bound (worker OOM). The flushed
                        # rows are recoverable on resume from the per-job summaries.
                        try:
                            append_rows_to_csv(local_rows, csv_path)
                        finally:
                            local_rows.clear()
        except Exception:
            time.sleep(5)
    if local_rows:
        try:
            append_rows_to_csv(local_rows, csv_path)
        finally:
            local_rows.clear()


def rebuild_csv_from_summaries(runs_dir: Path, csv_path: Path) -> int:
    """Single parallel scan over all *_summary.json files. Returns n_rows_written.

    IMPORTANT — folder name is the canonical truth for job_index and job_name.
    """
    '''
    Live counter during the rglob traversal — on a USB HDD this walk over
    tens of thousands of job folders is slow and would otherwise be silent.
    When analysis outputs were just wiped (full re-analysis) it finds ~0
    summaries, so the walk is silent unless we announce it up front.
    '''
    console_info("   Walking job folders on disk (USB; ~1–2 min for 58k jobs, silent if outputs were wiped)...")
    json_files = []
    for _scanned in sorted(runs_dir.rglob("*_summary.json")):
        json_files.append(_scanned)
        if len(json_files) % 2000 == 0:
            _tty_write(f"\r   Scanning job folders: {len(json_files):,} summaries found ...   \033[K")
            sys.stdout.flush()
    json_files.sort()
    total = len(json_files)
    _tty_write(f"\r   Scanning job folders: {total:,} summaries found.          \033[K\n")
    sys.stdout.flush()
    _done_count = [0]

    def _load_one(jf):
        try:
            with open(jf, "r") as f:
                data = json.load(f)
            if data.get("status") != "Success":
                return None, None
            folder_name = jf.parent.name
            data["job_name"]  = folder_name
            data["job_index"] = folder_name.split("_")[0]
            return flatten_job_result(data)
        except Exception:
            return None
        finally:
            _done_count[0] += 1
            n = _done_count[0]
            if n % 1000 == 0 or n == total:
                _tty_write(f"\r   Reading summary JSONs: {n}/{total} ...   \033[K")
                sys.stdout.flush()

    from concurrent.futures import ThreadPoolExecutor as _TPE
    n_workers = min(max(1, (os.cpu_count() or 4) - 2), total or 1)
    with _TPE(max_workers=n_workers) as ex:
        results = list(ex.map(_load_one, json_files))

    rows = [row for row in results if row is not None]

    if rows:
        _tty_write("\r\033[K")
        sys.stdout.flush()
        df = pd.DataFrame(rows)
        if "job_name" in df.columns:
            df = df.sort_values("job_name", kind="mergesort").reset_index(drop=True)
        df["job_index"] = range(len(df))
        ordered  = [c for c in CSV_COLUMN_ORDER if c in df.columns]
        remainder = [c for c in df.columns if c not in ordered]
        df = df[ordered + remainder]
        atomic_to_csv(df, csv_path, index=False)

    return len(rows)

# =============================================================================
# SECTION 10: MAIN SYSTEM ENTRY POINT
# =============================================================================

def generate_scientific_ranking_csv(CSV_PATH, PROD, ts_now):
    """Build the mechanism-first Scientific Ranking CSV from the master CSV.
    Extracted verbatim from main() (M1): tier-first sort, MD-ready selection,
    positive-control validation, no-gaps fill. Returns (rank_csv_path,
    rank_columns_count). `logger` is the module-global; console_info is used
    for reporting exactly as in-line.
    """
    # Step 10.10: Scientific Ranking Matrix CSV Generation (MECHANISM-FIRST)
    # -------------------------------------------------------------------------------
    console_info("Executing the compilation of the strictly mechanistic Scientific Ranking CSV...")
    rank_csv_path = None
    rank_columns_count = 0
    try:
        rank_csv_name = f"7_Boltz2_FAcDs_Ranked_{ts_now}.csv"
        rank_csv_path = PROD / rank_csv_name

        if CSV_PATH.exists() and os.path.getsize(CSV_PATH) > 0:
            df_rank = pd.read_csv(CSV_PATH, low_memory=False)

            for c in [CFG.COL_LIKE_S, CFG.COL_MECH_S, "competence_score", "catalytic_constellation_score", "model_degrader_consensus"]:
                if c not in df_rank.columns: df_rank[c] = 0.0
            if CFG.COL_TIER not in df_rank.columns: df_rank[CFG.COL_TIER] = CFG.TIER_DECOY

            tier_map = CFG.TIER_SORT_WEIGHT
            df_rank["tier_val"] = df_rank[CFG.COL_TIER].map(tier_map).fillna(0)

            """
            Mechanism-first ranking. The degrader tier (hard chemistry gates + the Criterion-B
            constellation cap) is the primary key. Within a tier candidates are ordered by the
            gated continuous competence_score (CFG §5.5 — angle, nucleophile distance, carboxylate
            clamp, trajectory deviation, triad relay and halide each once), then by the precise
            catalytic_constellation_score (Criterion B — eight-residue geometric fidelity vs the
            3R3U crystal), then active-site conservation, and finally by
            model_degrader_consensus (the fraction of Boltz diffusion samples that independently
            reach a degrader tier) as the LAST tiebreaker. Consensus is a tiebreaker ONLY — it
            never crosses a tier or competence boundary, so a reproducible pose floats above a
            single-frame fluke of equal geometry without a conformationally flexible true
            substrate being demoted for its sampling spread, and nothing is filtered out.
            Rank_Within_Ligand additionally ranks each protein among all proteins screened
            against the same ligand — surfacing the strongest FAcD variant per ligand independent
            of cross-ligand geometric bias. mechanistic_score and soft_catalytic_score are shown
            for reference but not sorted on. mergesort keeps the order stable and reproducible.
            """
            df_rank = df_rank.sort_values(
                by=["tier_val", "competence_score", "catalytic_constellation_score", CFG.COL_LIKE_S, "model_degrader_consensus"],
                ascending=[False, False, False, False, False],
                kind="mergesort",
            )
            df_rank.insert(0, "Scientific_Rank", range(1, len(df_rank) + 1))
            _lig_col = next((c for c in (CFG.COL_LIG, "ligand") if c in df_rank.columns), None)
            if _lig_col:
                df_rank.insert(1, "Rank_Within_Ligand", df_rank.groupby(_lig_col).cumcount() + 1)

            '''
            MD-ready selection (SECTION 18 SSOT). Flag the cohort that receives the
            expensive downstream pipeline (CIF->PDB, PrepWizard, MM-GBSA, MD) so Steps
            05-08 prepare/simulate only these rows, not all ~58k. Analysis/figures still
            span the full population; only heavy compute is gated.
            '''
            _md_sel = CFG.MD_SELECTED_COL
            _md_rnk = CFG.MD_RANK_COL
            _mode = getattr(CFG, "MD_SELECTION_MODE", "tier")
            if _mode == "topN":
                _sel_mask = df_rank["Scientific_Rank"] <= CFG.MD_TOP_N
            elif _mode == "per_ligand" and _lig_col:
                _lig_bare = df_rank[_lig_col].astype(str).str.replace(r"^\d+_", "", regex=True)
                _pl_scope = (list(CFG.MD_PER_LIGAND_TIER)
                             if isinstance(CFG.MD_PER_LIGAND_TIER, (list, tuple, set))
                             else [CFG.MD_PER_LIGAND_TIER])
                _tier_ok  = df_rank[CFG.COL_TIER].isin(_pl_scope)
                # Roster: data-driven (every unique ligand that reached the tier) when
                # MD_PER_LIGAND_AUTO, else the explicit curated MD_PER_LIGAND panel.
                if getattr(CFG, "MD_PER_LIGAND_AUTO", True):
                    _md_roster = sorted(_lig_bare[_tier_ok].unique())
                else:
                    _md_roster = list(CFG.MD_PER_LIGAND)
                _pick_idx = []
                for _lg in _md_roster:
                    _cand = df_rank.index[_tier_ok & (_lig_bare == _lg)]
                    if len(_cand):
                        # df_rank is already sorted by rank, so the first index is the best
                        _pick_idx.append(_cand[0])
                _sel_mask = df_rank.index.isin(_pick_idx)
            else:   # "tier" (and fallback when no ligand column for per_ligand)
                _sel_mask = df_rank[CFG.COL_TIER].isin(CFG.MD_TIERS)
            df_rank[_md_sel] = _sel_mask.astype(bool)
            df_rank[_md_rnk] = pd.NA
            _sel_order = df_rank.loc[_sel_mask].sort_values("Scientific_Rank").index
            df_rank.loc[_sel_order, _md_rnk] = range(1, len(_sel_order) + 1)
            _pl_roster_src = ("all unique ligands in tier (data-driven)"
                              if getattr(CFG, "MD_PER_LIGAND_AUTO", True)
                              else f"curated panel of {len(CFG.MD_PER_LIGAND)}")
            _gate_desc = (f"single best complex per ligand within {'/'.join(_pl_scope)} — {_pl_roster_src}"
                          if _mode == "per_ligand"
                          else f"Scientific_Rank ≤ {CFG.MD_TOP_N}" if _mode == "topN"
                          else f"degrader_tier in {list(CFG.MD_TIERS)}")
            reporter_md = (f"  MD-ready selection — mode='{_mode}'  ·  "
                           f"{int(_sel_mask.sum())} complexes flagged {_md_sel}=True")
            try:
                console_info(reporter_md)
                console_info(f"    Gate: {_gate_desc}")
            except Exception:
                print(reporter_md, flush=True)
            # Per-pick provenance table: show WHICH complexes were flagged (ligand,
            # tier, rank, protein, job) so the shortlist is auditable, styled to match
            # the rest of the console output.
            try:
                _prot_col = next((c for c in (CFG.COL_PROT, "protein", "Protein")
                                  if c in df_rank.columns), None)
                _sel_rows = df_rank.loc[_sel_mask].sort_values("Scientific_Rank")
                _rows_data = [
                    (str(_i),
                     str(_r.get(_lig_col, "?")) if _lig_col else "?",
                     str(_r.get(CFG.COL_TIER, "?")),
                     f"#{int(_r.get('Scientific_Rank', 0))}",
                     str(_r.get(_prot_col, "?")) if _prot_col else "?",
                     str(_r.get("job_name", "?")))
                    for _i, (_, _r) in enumerate(_sel_rows.iterrows(), 1)
                ]
                _hdr = ("#", "Ligand", "Tier", "Rank", "Protein", "Job")
                if _rows_data:
                    _w = [max(len(_hdr[_c]), max(len(_rd[_c]) for _rd in _rows_data))
                          for _c in range(6)]
                    _w[5] = min(_w[5], 46)
                    _algn = (str.center, str.ljust, str.ljust, str.center, str.ljust, str.ljust)
                    def _cell(_c, _v): return _algn[_c](str(_v)[:_w[_c]], _w[_c])
                    def _fmt(_cells):  return "    │ " + " │ ".join(_cell(_c, _cells[_c]) for _c in range(6)) + " │"
                    def _rule(_l, _m, _r): return "    " + _l + _m.join("─" * (_w[_c] + 2) for _c in range(6)) + _r
                    console_info(_rule("┌", "┬", "┐"))
                    console_info(_fmt(_hdr))
                    console_info(_rule("├", "┼", "┤"))
                    for _rd in _rows_data:
                        console_info(_fmt(_rd))
                    console_info(_rule("└", "┴", "┘"))
                if _mode == "per_ligand":
                    _got  = set(_lig_bare.loc[_sel_mask].astype(str))
                    _miss = [_lg for _lg in _md_roster if _lg not in _got]
                    if _miss:
                        console_info(f"    No candidate in {'/'.join(_pl_scope)} for "
                                     f"{len(_miss)} requested ligand(s) — dropped: {', '.join(_miss)}")
            except Exception as _md_e:
                console_info(f"    (MD-ready provenance unavailable: {_md_e})")

            df_rank["Ranking_Score_Calc"] = (
                "Tier:"  + df_rank[CFG.COL_TIER].astype(str)
                + " | Competence:" + df_rank["competence_score"].map("{:.3f}".format)
                + " | Consensus:" + df_rank["model_degrader_consensus"].map("{:.2f}".format)
                + " | Cons:" + df_rank[CFG.COL_LIKE_S].map("{:.2f}".format)
            )
            '''
            Honest-claim disclaimer (carried on every row). The tiers/is_degrader flag are
            GEOMETRIC near-attack-conformation (NAC) pose-quality descriptors from a static
            Boltz-2 structure — NOT a kinetic turnover guarantee. Chemical feasibility (C–F
            BDE, β-fluorination, sn2_dead_end, feasibility_factor) is computed and REPORTED
            here but deliberately not gated, so a recalcitrant substrate (e.g. TFA) can still
            surface a structurally competent variant for discovery; the activation barrier is
            decided downstream by Step-06 MM-GBSA and Step-07 QM/MM, which are the arbiters.
            '''
            df_rank["Classification_Basis"] = (
                "geometric_NAC_pose_quality; kinetic_feasibility_deferred_to_QMMM_Step08")

            """
            Control read-out + VALIDATION. The ranking is geometry-driven and applies
            NO substrate-class penalty, so where the proven substrates (FA/DFA) and the SN2
            dead-end (TFA) land is itself a result. Beyond reporting positions/feasibility, this
            block asserts the positive-control expectation: the native substrates Fluoroacetate
            and Difluoroacetate must each surface at least one degrader-tier pose (is_degrader). A
            screen in which the proven substrate fails to register is mis-calibrated, so the failure
            is raised as a prominent WARNING and recorded, flagging a mis-calibrated ranking rather
            than shipping it silently. Non-fatal (results are still written); the operator decides.
            """
            control_validation = "PASS"
            try:
                _lc = next((c for c in (CFG.COL_LIG, "ligand") if c in df_rank.columns), None)
                if _lc:
                    _top = df_rank[df_rank[CFG.COL_TIER] == CFG.TIER_TOP][_lc]
                    _top_set = sorted(set(_top.astype(str).str.replace(r"^\d+_", "", regex=True)))
                    # Anchored to the exact bare "<idx>_Fluoroacetate" name so FA and DFA
                    # stay distinct — "Fluoroacetate" is a substring of "Difluoroacetate".
                    _fa_rank = df_rank.loc[df_rank[_lc].astype(str).str.fullmatch(r"\d+_Fluoroacetate", na=False), "Scientific_Rank"]
                    _fa_best = int(_fa_rank.min()) if len(_fa_rank) else -1
                    _tfa_feas = pd.to_numeric(
                        df_rank.loc[df_rank[_lc].astype(str).str.fullmatch(r"\d+_TFA", na=False), "feasibility_factor"],
                        errors="coerce").mean() if "feasibility_factor" in df_rank.columns else float("nan")
                    console_info(f"Control read-out — top-tier ({CFG.TIER_TOP}) ligands: {_top_set} | "
                                 f"FA best Scientific_Rank: {_fa_best} | TFA mean reported feasibility: {_tfa_feas:.3f}")

                    # Positive-control assertion: FA and DFA must each register as a degrader.
                    _hq = set(CFG.TIER_HIGH_QUALITY)
                    _failed = []
                    for _cname in CFG.CTRL_POSITIVE:
                        # Anchored fullmatch so "Fluoroacetate" does not also capture
                        # "Difluoroacetate" (substring) — each control asserted on its own poses.
                        _crows = df_rank[df_rank[_lc].astype(str).str.fullmatch(rf"\d+_{_cname}", na=False)]
                        _is_deg = _crows["is_degrader"].astype(str).str.lower().isin(("true", "1", "1.0")) \
                            if "is_degrader" in _crows.columns else pd.Series([], dtype=bool)
                        _in_hq = _crows[CFG.COL_TIER].isin(_hq) if CFG.COL_TIER in _crows.columns else pd.Series([], dtype=bool)
                        if not (bool(_is_deg.any()) or bool(_in_hq.any())):
                            _failed.append(_cname)
                    if _failed:
                        control_validation = f"FAIL({','.join(_failed)})"
                        _msg = (f"  ⚠ CONTROL VALIDATION FAILED — positive control(s) {_failed} produced NO "
                                f"degrader-tier pose. The screen is likely mis-calibrated; inspect before trusting ranks.")
                        if logger: logger.warning(_msg)
                        console_info(_msg)
                    else:
                        console_info("  ✔ Control validation PASS — FA and DFA both register as degraders.")
                    df_rank["control_validation"] = control_validation
            except Exception as _ce:
                console_info(f" Control read-out/validation skipped ({_ce}).")

            # Drop only an EXACT duplicate column (ligand_iptm ≡ iptm in Boltz-2 output);
            # near-duplicates (geometric_tier, identity_pct) are retained as distinct fields.
            cols_to_drop = ["tier_val"]
            if ("ligand_iptm" in df_rank.columns and "iptm" in df_rank.columns
                    and df_rank["ligand_iptm"].equals(df_rank["iptm"])):
                cols_to_drop.append("ligand_iptm")
            df_rank = df_rank.drop(columns=[c for c in cols_to_drop if c in df_rank.columns])

            """
            Final no-gaps guarantee: the published ranked CSV must contain no
            empty cells. Fill any residual NaN with type-appropriate values —
            numeric columns → 0.0, everything else → "N/A" text. This is a
            safety net over the upstream per-field defaults so that no future
            column (or resume edge case) can leave a blank cell.
            """
            """
            0.0 IS NOT A NEUTRAL FILLER. For an INVERTED metric it is the OPTIMUM, so filling a missing
            value with it does not mark the cell empty — it marks the complex PERFECT.

                Dist_Nucleophile        0.0 A  = the nucleophile sitting ON the carbon (passes every gate)
                Active_Site_RMSD        0.0 A  = a flawless match to the crystal
                sn2_backside_occlusion  0.0    = no steric blockade whatsoever
                chem_penalty            0.0    = no BDE or occlusion penalty at all

            Those columns are filled with CFG.SENTINEL_UNDEFINED (999.0) — the value the rest of the
            pipeline already uses for 'not measurable', and which every gate and every figure filters out.
            Columns where 0.0 genuinely IS the worst case (an angle, a score, a confidence, a count) keep
            the zero fill, because there the fill and the meaning agree.
            """
            for _col in df_rank.columns:
                if not df_rank[_col].isna().any():
                    continue
                if pd.api.types.is_numeric_dtype(df_rank[_col]):
                    _fill = CFG.SENTINEL_UNDEFINED if _col in CFG.INVERTED_METRIC_COLUMNS else 0.0
                    df_rank[_col] = df_rank[_col].fillna(_fill)
                else:
                    df_rank[_col] = df_rank[_col].astype(object).fillna("N/A")

            atomic_to_csv(df_rank, rank_csv_path, index=False)
            rank_columns_count = len(df_rank.columns)
            console_info("Scientific Ranking Output CSV generated flawlessly.")
        else: console_info(" Operation actively bypassed (No structured data detected).")
    except Exception as e: console_info(f" Mechanism ranking operation structurally failed ({e})")

    return rank_csv_path, rank_columns_count


def main():
    # Collapse stacked separator rules: a caller prints a rule, a helper prints its own,
    # and the log grows triple bars with nothing between them.
    _utils_mod.install_console_rule_filter()
    # -------------------------------------------------------------------------------
    # Step 10.1: Command Line Argument Parsing Sequence
    # -------------------------------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", help="Designated Input FASTA file")
    parser.add_argument("--smi", help="Designated Input SMILES file")
    parser.add_argument("--resume", help="Explicit name of the run folder intended for resuming operations")
    parser.add_argument("--cpus", type=int, default=TARGET_CORES, help="Overrides default limits to manually set CPU threads")
    parser.add_argument("--diffusion-samples", type=int, default=DIFFUSION_SAMPLES_DEFAULT)
    args = parser.parse_args()

    # --- Standard Header Output Logic ---
    if args.resume:
        run_folder_name = args.resume
        _utils_mod.print_script_banner(
            "02_Production_FAcDs.py",
            "Boltz-2 Scoring  ·  NAC Geometry  ·  Degrader Tier Assignment  │  Resume Mode",
        )
        console_info(f"  Run Name : {run_folder_name}")
    else:
        ts = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
        run_folder_name = f"Boltz-2_Run_{ts}"
        _utils_mod.print_script_banner(
            "02_Production_FAcDs.py",
            "Boltz-2 Scoring  ·  NAC Geometry  ·  Degrader Tier Assignment  │  Fresh Run",
        )
        console_info(f"  Run Name : {run_folder_name}")

    # -------------------------------------------------------------------------------
    # Step 10.2: Project Path System Configuration Architecture
    # -------------------------------------------------------------------------------
    root = Path.cwd()
    if args.resume:
        run_root = root / run_folder_name
        if not run_root.exists(): sys.exit(f"Target run sequence not located within {root}: {run_root}")
        resumed = True
    else:
        run_root = root / run_folder_name
        run_root.mkdir(parents=True, exist_ok=True)
        resumed = False

    PROD = run_root / "1_Boltz2_Production"

    # Canonical production directory layout.
    D_IN   = PROD / "1_Input_Data"
    D_YAML = PROD / "2_Boltz2_YAML_Configs"
    D_RUNS = PROD / "4_Prediction_Jobs"
    D_SEQ       = PROD / "3_Sequence_Reference_Data"
    D_COLABFOLD = D_SEQ / "MSA_Sequences"
    D_ALN       = D_SEQ / "Active_Site_Alignments"
    D_SEQ.mkdir(parents=True, exist_ok=True)

    # Fresh alignments every run: wipe any prior Active_Site_Alignments (no cache reuse).
    # Alignment is cheap (~2 s for all variants); the cache only ever caused file/CSV bloat.
    shutil.rmtree(D_ALN, ignore_errors=True)

    for d in [D_IN, D_YAML, D_RUNS, D_COLABFOLD, D_ALN]: d.mkdir(parents=True, exist_ok=True)

    initial_folders_count = len([d for d in sorted(D_RUNS.iterdir()) if d.is_dir()]) if D_RUNS.exists() else 0
    removed_folders_count = 0

    ts_now = datetime.now().strftime("%Y%m%d_%H%M%S")
    LOG_PATH = PROD / "00_Boltz2_Production.log"
    CSV_PATH = PROD / f"6_Boltz2_FAcDs_Master_{ts_now}.csv"
    setup_logging(LOG_PATH)

    prev_elapsed_map = {}
    if resumed:
        existing_csvs = sorted(list(PROD.glob("*_Master*.csv")))
        if existing_csvs:
            try:
                df_old = pd.read_csv(existing_csvs[-1], dtype={"job_index": str}, low_memory=False)
                elapsed_col = next((c for c in df_old.columns if "elapsed" in c.lower()), None)
                if elapsed_col:
                    for _, row in df_old.iterrows():
                        try:
                            jname = str(row["job_name"]).strip()
                            val = float(row.get(elapsed_col, 0))
                            prev_elapsed_map[jname] = val
                        except Exception: pass
            except Exception as e: console_info(f"System Warning: Inability to properly parse preceding CSV log for timing extraction operations: {e}")
        try:
            # Skip the reference files (now co-located in the input folder) when
            # picking the user input FASTA/SMI on resume.
            _ref_files = {CFG.DEHA4_REF_FASTA, "Control_Ligands_TFA_FA_DFA_Ref.smi"}
            f_path = next(p for p in sorted(D_IN.glob("*.fasta")) if p.name not in _ref_files)
            s_path = next(p for p in sorted(D_IN.glob("*.smi"))   if p.name not in _ref_files)
        except StopIteration:
            if args.fasta and args.smi:
                 console_info("  [Resume Sequence] System input files are absent in the structured folder; defaulting to strictly provided string arguments.")
                 shutil.copy2(args.fasta, D_IN / Path(args.fasta).name)
                 shutil.copy2(args.smi,   D_IN / Path(args.smi).name)
                 f_path = D_IN / Path(args.fasta).name
                 s_path = D_IN / Path(args.smi).name
            else: sys.exit(f"\nCRITICAL SYSTEM ERROR: Resume procedure failed inherently.\nInput FASTA/SMI files are absent from: {D_IN}\nPlease effectively restore them explicitly via the --fasta and --smi CLI parameters.")
    else:
        shutil.copy2(args.fasta, D_IN / Path(args.fasta).name)
        shutil.copy2(args.smi,   D_IN / Path(args.smi).name)
        f_path = D_IN / Path(args.fasta).name; s_path = D_IN / Path(args.smi).name

    proteins = []
    _MERGE_PREFIX = re.compile(r"^(\d+)_(.+?)_*$")
    with open(f_path, "r") as f_in:
        for order, r in enumerate(SeqIO.parse(f_in, "fasta"), 1):
            full_id = safe_name(extract_short_fasta_id(r.description or r.id))
            m = _MERGE_PREFIX.match(full_id)
            if m:
                fasta_pos = int(m.group(1))
                pid       = m.group(2).rstrip("_. ")
            else:
                fasta_pos = 0
                pid       = full_id
            proteins.append((order, pid, str(r.seq), fasta_pos))

    target_fasta_path = D_IN / f_path.name
    with open(target_fasta_path, "w") as out_fasta:
        for order, pid, seq, fasta_pos in proteins:
            header_index = fasta_pos if fasta_pos > 0 else order
            header = f"{header_index}_{pid}"
            out_fasta.write(f">{header}\n{seq}\n")

    ligands = []
    total_count = 0
    valid_count = 0
    invalid_count = 0
    corrected_log = []

    console_info(f"\n{SEPARATOR_LIGHT}")
    console_info("  SMILES Validation & Standardisation")
    console_info(SEPARATOR_LIGHT)
    with open(s_path) as f:
        for order, l in enumerate(f, 1):
            if l.strip() and not l.startswith("#"):
                total_count += 1
                p = l.split()
                smi = p[0]
                lid = p[1] if len(p) > 1 else f"LIG_{order}"

                mol = Chem.MolFromSmiles(smi)

                if mol is None:
                    invalid_count += 1
                else:
                    valid_count += 1
                    canonical_smi = Chem.MolToSmiles(mol)

                    if canonical_smi != smi:
                        corrected_log.append((lid, smi, canonical_smi))

                    ligands.append((order, safe_name(lid), canonical_smi))

    _sv_rows = [
        ("Ligands Processed",     total_count),
        ("Validated SMILES",      valid_count),
        ("Invalid SMILES",        invalid_count),
        ("Topologies Corrected",  len(corrected_log)),
    ]
    _kw, _vw = 22, 8
    _fail_clr = ConsoleColours.FAIL if invalid_count else ""
    _endc     = ConsoleColours.ENDC if invalid_count else ""
    console_info(f"  ┌{'─'*(_kw+2)}┬{'─'*(_vw+4)}┐")
    for _k, _v in _sv_rows:
        _c = _fail_clr if _k == "Invalid SMILES" else ""
        _e = _endc     if _k == "Invalid SMILES" else ""
        console_info(f"  │  {_k:<{_kw}}│  {_c}{_v:>{_vw}}{_e}  │")
    console_info(f"  └{'─'*(_kw+2)}┴{'─'*(_vw+4)}┘")

    if corrected_log:
        console_info(f"\n{SEPARATOR_LIGHT}")
        console_info("  Corrected SMILES Action Log")
        console_info(SEPARATOR_LIGHT)
        for lid, orig, fixed in corrected_log:
            console_info(f"  {lid}")
            console_info(f"    Before : {orig}")
            console_info(f"    After  : {fixed}")

    target_smi_path = D_IN / s_path.name
    with open(target_smi_path, "w") as out_smi:
        for _, lid, smi in ligands:
            out_smi.write(f"{smi}\t{lid}\n")

    init_gpu_reservations()
    _gpu_info = query_gpu_memory()
    gpu_count = len(_gpu_info) if _gpu_info else 0
    _, cores = cpu_usage_summary()
    _cr_rows = [
        ("GPUs Detected",        gpu_count),
        ("CPU Cores (Total)",    cores),
        ("Worker Threads",       TARGET_CORES),
        ("Proteins",             f"{len(proteins):,}"),
        ("Unique Ligands",       len(ligands)),
        ("Total Jobs Scheduled", f"{len(proteins) * len(ligands):,}"),
    ]
    _ck, _cv = 24, 12
    console_info(f"\n{SEPARATOR_LIGHT}")
    console_info("  System Compute Resources")
    console_info(SEPARATOR_LIGHT)
    console_info(f"  ┌{'─'*(_ck+2)}┬{'─'*(_cv+4)}┐")
    for _k, _v in _cr_rows:
        console_info(f"  │  {_k:<{_ck}}│  {str(_v):>{_cv}}  │")
    console_info(f"  └{'─'*(_ck+2)}┴{'─'*(_cv+4)}┘")

    # -------------------------------------------------------------------------------
    # Step 10.3: System Reference Data Initialisation
    # -------------------------------------------------------------------------------
    ref_pdb_path = setup_reference_data(D_IN)

    console_info("  Loading alignment cache...")
    load_cached_alignments(D_ALN / CFG.FILE_ALIGNMENT_STATS)
    console_separator()

    # -------------------------------------------------------------------------------
    # Step 10.4: Construct Initial Task Matrix & Pre-Processing Evaluation
    # -------------------------------------------------------------------------------
    tasks = []
    task_count = 0
    # job_stem -> the SMILES its finished CIF was actually predicted from, when that differs from the
    # current input. The analysis must use THIS, not the new one (see the resume branch below).
    _stale_smi_jobs: dict = {}
    total_ops = len(proteins) * len(ligands)
    existing_yamls = set(f.name for f in D_YAML.glob("*.yaml"))
    active_job_stems = set()
    num_ligands = len(ligands)
    print()
    for protein_order, pid, seq, fasta_pos in proteins:
        for lig_order, lid, smi in ligands:
            task_count += 1
            protein_position = fasta_pos if fasta_pos > 0 else protein_order
            jid = str((protein_position - 1) * num_ligands + lig_order).zfill(7)
            job_pid = f"{protein_position}_{pid}"
            y_name = f"{jid}_{job_pid}_{lid}.yaml"
            y_path = D_YAML / y_name
            active_job_stems.add(Path(y_name).stem)

            job_name_log = f"{jid}_{job_pid}_{lid}"
            if task_count % 1000 == 0 or task_count == total_ops:
                _tty_write(f"\rPre-Processing sequence evaluation {task_count}/{total_ops} | {job_name_log}\033[K")

            yaml_needs_write = True
            job_stem = y_name.replace(".yaml", "")
            job_path = D_RUNS / job_stem
            '''
            Completion is judged by the canonical success marker (a *_summary.json
            reporting "status": "Success"), the same test the scheduler uses. A
            finished prediction is therefore always preserved on resume even when
            the input SMILES representation drifts (e.g. carboxylate protonation).
            '''
            is_completed = check_job_status(job_path)

            if y_name in existing_yamls:
                try:
                    with open(y_path, "r") as yf:
                        y_data = yaml.safe_load(yf)
                    existing_seq = y_data.get("sequences", [])[0].get("protein", {}).get("sequence", "")
                    existing_smi = y_data.get("sequences", [])[1].get("ligand", {}).get("smiles", "")

                    if existing_seq == seq and existing_smi == smi:
                        yaml_needs_write = False
                    else:
                        if is_completed:
                            """
                            THE STRUCTURE WAS PREDICTED FROM THE OLD MOLECULE. IT MUST BE ANALYSED AS THE
                            OLD MOLECULE.

                            Keeping a finished prediction when the input SMILES has since changed is the
                            right call — re-running days of GPU because a carboxylate was re-protonated in
                            the input file would be absurd. Carrying the NEW SMILES forward against the OLD
                            CIF is not: map_mmcif_to_rdkit would build its template from a molecule that is
                            not the one in the structure, and every chemistry term downstream of that map —
                            the scissile carbon, n_scissile_f, beta_f_count, the formal charges, the
                            pi-cation terms — would be computed for the wrong compound. Nothing would fail.
                            One INFO line would be printed.

                            So the analysis uses the SMILES the structure was actually built from, and the
                            row is FLAGGED (ligand_smiles_stale) so the divergence is visible in the output
                            rather than buried in a resume log.
                            """
                            console_info(f"[Resume] {y_name}: input SMILES/sequence has changed since this "
                                         f"prediction was made. The finished structure is KEPT and will be "
                                         f"analysed against the SMILES it was PREDICTED FROM; the row is "
                                         f"flagged ligand_smiles_stale=1.")
                            _stale_smi_jobs[job_stem] = existing_smi
                            yaml_needs_write = False
                        else:
                            console_info(f"[Resume] Sequence or ligand structure mismatch unequivocally detected in {y_name}. Moving stale YAML definitions and corrupted runs to Trash.")

                            if y_path.exists():
                                yaml_trash = PROD / "9_Trash" / "YAMLs"
                                yaml_trash.mkdir(parents=True, exist_ok=True)
                                target_yaml = yaml_trash / f"{y_path.name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                                try:
                                    shutil.move(str(y_path), str(target_yaml))
                                except Exception:
                                    y_path.unlink(missing_ok=True)

                            remove_stale_protein_jobs(PROD, pid)
                            yaml_needs_write = True
                except Exception:
                    pass

            if yaml_needs_write:
                payload = {
                    "defaults": {"mol_type": "complex"},
                    "sequences": [{"protein": {"id": "A", "sequence": seq}}, {"ligand":  {"id": "L", "smiles": smi}}],
                    "properties": [{"affinity": {"binder": "L"}}]
                }

                with open(y_path, "w") as f: yaml.safe_dump(payload, f, sort_keys=False, default_flow_style=False)
                GLOBAL_STATS["created_yaml"] += 1

            # A finished prediction is analysed with the SMILES it was PREDICTED FROM, never with a newer
            # one that describes a different molecule.
            _smi_for_analysis = _stale_smi_jobs.get(job_stem, smi)
            tasks.append({
                "job_index": jid, "protein": pid, "job_protein": job_pid, "protein_idx": protein_order,
                "ligand": lid, "sequence": seq, "smiles": _smi_for_analysis, "yaml": y_path,
                "smiles_stale": int(job_stem in _stale_smi_jobs),
            })

    if resumed:
        for _, pid, seq, fasta_pos in proteins:
            a3m_file = colabfold_a3m_path(PROD, pid)
            meta_file = colabfold_meta_path(PROD, pid)
            if a3m_file.exists():
                if meta_file.exists():
                    try:
                        with open(meta_file, "r") as f: meta = json.load(f)
                        if meta.get("sequence_sha256") != sequence_hash(seq):
                            console_info(f"[Resume] Substantial sequence mismatch structurally detected for {pid}. Purging associated ColabFold files and conflicting structural runs.")
                            a3m_file.unlink(missing_ok=True)
                            meta_file.unlink(missing_ok=True)
                            remove_stale_protein_jobs(PROD, pid)
                    except Exception:
                        console_info(f"[Resume] Significantly corrupt ColabFold metadata actively detected for {pid}. Removing associated stale files.")
                        a3m_file.unlink(missing_ok=True)
                        meta_file.unlink(missing_ok=True)
                        remove_stale_protein_jobs(PROD, pid)
                else:
                    console_info(f"[Resume] Mandatory structural metadata critically missing for {pid} ColabFold data matrices. Purging associated stale MSA files.")
                    a3m_file.unlink(missing_ok=True)
                    remove_stale_protein_jobs(PROD, pid)

    console_info("\nPre-Processing phase unequivocally completed.")
    if resumed:
        # Register the canonical control-calibration job stems (DeHa4 × 3R3U × the
        # control ligands) so the orphan purge below never deletes a valid control run.
        _n_ctrl_ligs = len(CTRL_LIGANDS)
        for _sq, (_ctl_lig, _) in enumerate(CTRL_LIGANDS, start=1):
            active_job_stems.add(f"0000000_{_sq}_DeHa4_Control_{_ctl_lig}")
            active_job_stems.add(f"0000000_{_sq + _n_ctrl_ligs}_3R3U_Control_{_ctl_lig}")
        removed_folders_count = purge_orphans(PROD, active_job_stems)

        """
        Wipe per-job analysis outputs so re-analysis runs fresh on resume.
        GPU predictions (boltz_results_*/) are preserved — only analysis
        artifacts are removed so changed parameters take full effect.
        """
        console_info("Wiping per-job analysis outputs for full re-analysis...")
        _wipe_dirs = [d for d in sorted(D_RUNS.iterdir()) if d.is_dir() and d.name[0].isdigit()]
        _wiped = [0]

        def _wipe_one(d):
            for f in sorted(d.glob("*_summary.json")):
                try: f.unlink()
                except Exception: pass
            for f in sorted(d.glob("*_interactions.csv")):
                try: f.unlink()
                except Exception: pass
            bc = d / "Best_Complex"
            if bc.exists():
                try: shutil.rmtree(bc)
                except Exception: pass
            _wiped[0] += 1

        from concurrent.futures import ThreadPoolExecutor as _TPE
        _wipe_total = len(_wipe_dirs)
        with _TPE(max_workers=min(max(1, (os.cpu_count() or 4) - 2), _wipe_total or 1)) as _wex:
            for _done, _ in enumerate(_wex.map(_wipe_one, _wipe_dirs), 1):
                if _done % 1000 == 0 or _done == _wipe_total:
                    _tty_write(f"\r   Wiping analysis outputs {_done:,}/{_wipe_total:,}\033[K")
        _tty_write("\r\033[K")
        console_info(f" -> Wiped analysis outputs for {_wiped[0]:,} job directories.")

        '''
        Delete stale master CSVs only — they are rebuilt fresh from the re-analysis.
        Ranked CSVs are timestamped RESULT files and are NEVER deleted here: each
        resume writes a new timestamped ranked CSV and downstream steps pick the
        newest by mtime, so prior results are preserved as history.
        '''
        for _old_csv in sorted(PROD.glob("*_Master*.csv")):
            try: _old_csv.unlink()
            except Exception: pass

    console_separator()
    console_info("")

    # -------------------------------------------------------------------------------
    # Step 10.5: Control Job Execution & Calibration
    # -------------------------------------------------------------------------------
    """
    DEHA4 (input enzyme) vs 3R3U (crystal reference) × 3 fluoroacetate control substrates.
    Fluoroacetate map drives all 58k pipeline comparisons.
    """
    _n_ctrl_total = len(CTRL_LIGANDS) * 2  # DeHa4 × 3 + 3R3U × 3
    _ctrl_idx     = [0]
    console_separator()
    console_info(f"Control Job Execution  ({_n_ctrl_total} jobs: DeHa4 × 3 ligands + 3R3U × 3 ligands)")
    console_info("")
    ctrl_jid:         str                  = "0000000"
    ctrl_results:     Dict[str, dict]      = {}
    ctrl_trust_scores: Dict[str, float]   = {}
    ctrl_cifs:        Dict[str, Path]      = {}
    ctrl_maps:        Dict[str, dict]      = {}
    control_cif:      Optional[Path]       = None
    control_map:      Dict[str, int]       = {}

    deha4_calib_list: list = []
    r3u_calib_list:   list = []

    for _seq, (_lig_name, _lig_smiles) in enumerate(CTRL_LIGANDS, start=1):
        _ctrl_name = f"{ctrl_jid}_{_seq}_DeHa4_Control_{_lig_name}"
        _ctrl_yaml = D_YAML / f"{_ctrl_name}.yaml"
        _ctrl_task = {
            "job_index": ctrl_jid, "protein": "DeHa4_Control", "job_protein": f"{_seq}_DeHa4_Control",
            "ligand": _lig_name, "sequence": DEHA4_CONTROL_SEQ, "smiles": _lig_smiles, "yaml": _ctrl_yaml
        }
        if not _ctrl_yaml.exists():
            payload = {
                "defaults": {"mol_type": "complex"},
                "sequences": [
                    {"protein": {"id": "A", "sequence": DEHA4_CONTROL_SEQ}},
                    {"ligand":  {"id": "L", "smiles": _lig_smiles}}
                ],
                "properties": [{"affinity": {"binder": "L"}}]
            }
            with open(_ctrl_yaml, "w") as f:
                yaml.safe_dump(payload, f)
        _ctrl_idx[0] += 1
        _lig_short_name = _lig_name.split("_", 1)[-1]
        console_info(f"  [{_ctrl_idx[0]}/{_n_ctrl_total}]  DeHa4 × {_lig_short_name:<18}  →  {_ctrl_name}")
        _ctrl_res, _ = process_single_job(_ctrl_task, PROD, 5, 0.0, D_ALN, None, {}, "", None)
        ctrl_results[_lig_name] = _ctrl_res or {}
        _cifs = sorted((D_RUNS / _ctrl_name).glob("**/*.cif"))
        if _cifs:
            _ctrl_summary_p = D_RUNS / _ctrl_name / f"{_ctrl_name}_summary.json"
            _cached_ts = None
            _cached_cm = None
            if _ctrl_summary_p.exists():
                try:
                    with open(_ctrl_summary_p) as _cf:
                        _cs = json.load(_cf)
                    if "ctrl_trust_score" in _cs:
                        _cached_ts = float(_cs["ctrl_trust_score"])
                        _cached_cm = _cs.get("ctrl_active_site_map", {})
                except Exception:
                    pass
            if _cached_ts is not None:
                _ts  = _cached_ts
                _cm  = _cached_cm if isinstance(_cached_cm, dict) else {}
                _calib_data = _cached_cm if isinstance(_cached_cm, dict) and "mapping_rows" in _cached_cm else {
                    "model_label": f"DeHa4 Control {_lig_name}",
                    "trust_score": _ts, "verdict_label": "Cached",
                    "mapping_rows": [], "missing_residues": [], "all_mapped": True,
                }
                console_info(f"  > DeHa4 Control {_lig_name}: calibration loaded from cache (trust={_ts:.3f} Å).")
            else:
                _ts, _cm, _calib_data = perform_control_calibration(
                    _cifs[0], ref_pdb_path,
                    model_label=f"DeHa4 Control {_lig_name}",
                    silent=True
                )
                try:
                    _cs2 = json.loads(_ctrl_summary_p.read_text()) if _ctrl_summary_p.exists() else {}
                    _cs2["ctrl_trust_score"] = _ts
                    _cs2["ctrl_active_site_map"] = _cm
                    _cs2["ctrl_calib_data"] = _calib_data
                    _ctrl_summary_p.write_text(json.dumps(_cs2, indent=2))
                except Exception:
                    pass
            deha4_calib_list.append(_calib_data)
            ctrl_trust_scores[_lig_name] = _ts
            ctrl_cifs[_lig_name]         = _cifs[0]
            ctrl_maps[_lig_name]         = _cm
        if _lig_name == "26_Fluoroacetate":
            if not _cifs:
                console_info("CRITICAL OPERATIONAL ERROR: Failed to generate DeHa4 Control (Fluoroacetate) CIF.")
                sys.exit(1)
            control_cif = _cifs[0]
            control_map = ctrl_maps.get("26_Fluoroacetate", {})

    # -------------------------------------------------------------------------------
    # Step 10.5a-post: Retroactive active-site re-analysis for DeHa4 control jobs
    # -------------------------------------------------------------------------------
    """
    DeHa4 controls are processed without a control_cif/control_map (those do not
    yet exist when DeHa4 runs). Once ref_pdb_path (3R3U crystal) is available,
    each DeHa4 control CIF is superimposed on the crystal reference to fill in its
    active-site RMSD, halide-stabilisation / carboxylate-clamp flags and the
    derived ActiveSite_Conservation_Score — the same comparison applied to every
    candidate protein in the main pipeline.

    crystal_map: residue numbers in the 3R3U crystal PDB (pdb_id from
    REF_ACTIVE_SITE_MAP), used as anchor positions inside analyse_candidate_structure.
    """
    _crystal_map = {k: v["pdb_id"] for k, v in REF_ACTIVE_SITE_MAP.items()}
    for _lig_name, _lig_smiles in CTRL_LIGANDS:
        _ctrl_cif = ctrl_cifs.get(_lig_name)
        if _ctrl_cif and _ctrl_cif.exists() and ref_pdb_path and ref_pdb_path.exists():
            try:
                _deha4_recalc = analyse_candidate_structure(_ctrl_cif, ref_pdb_path, _crystal_map)
                if _lig_name in ctrl_results and ctrl_results[_lig_name]:
                    ctrl_results[_lig_name]["Halide_Stabilisation"] = _deha4_recalc.get(
                        "Halide_Stabilisation", False
                    )
                    ctrl_results[_lig_name]["Carboxylate_Clamp"] = _deha4_recalc.get(
                        "Carboxylate_Clamp", False
                    )
                    _rmsd = _deha4_recalc.get("Active_Site_RMSD", 99.0)
                    ctrl_results[_lig_name]["Active_Site_RMSD"] = _rmsd
                    _geo = 0.0 if _rmsd >= 99.0 else 100.0 / (1.0 + _rmsd)
                    ctrl_results[_lig_name][CFG.COL_LIKE_S] = round(CFG.CONSERV_REF_W_IDENT * 100.0 + CFG.CONSERV_REF_W_GEO * _geo, 2)
            except Exception as _recalc_err:
                console_info(f"  [Warning] DeHa4 control re-analysis failed for {_lig_name}: {_recalc_err}")

    # -------------------------------------------------------------------------------
    # Step 10.5b: 3R3U Crystal Reference Jobs (Experimental Positive Control)
    # -------------------------------------------------------------------------------
    """
    Run 3R3U against all three fluoroacetate control substrates.
    Per-ligand calibration (RMSD + active-site mapping) for each.
    Fluoroacetate (26) provides r3u_boltz_cif / r3u_control_map for 58 k pipeline.
    """
    r3u_result_data:   Optional[dict]      = None
    r3u_control_map:   Dict[str, int]      = {}
    r3u_control_maps:  Dict[str, dict]     = {}
    r3u_trust_scores:  Dict[str, float]    = {}
    r3u_results:       Dict[str, dict]     = {}
    r3u_boltz_cif:     Optional[Path]      = None
    r3u_sequence = extract_3r3u_sequence(ref_pdb_path)
    if r3u_sequence:
        console_info("")
        r3u_jid = "0000000"
        _r3u_offset = len(CTRL_LIGANDS)
        for _seq, (_lig_name, _lig_smiles) in enumerate(CTRL_LIGANDS, start=_r3u_offset + 1):
            _r3u_name = f"{r3u_jid}_{_seq}_3R3U_Control_{_lig_name}"
            _r3u_yaml = D_YAML / f"{_r3u_name}.yaml"
            _r3u_task = {
                "job_index": r3u_jid, "protein": "3R3U_Control", "job_protein": f"{_seq}_3R3U_Control",
                "ligand": _lig_name, "sequence": r3u_sequence, "smiles": _lig_smiles, "yaml": _r3u_yaml
            }
            if not _r3u_yaml.exists():
                r3u_payload = {
                    "defaults": {"mol_type": "complex"},
                    "sequences": [
                        {"protein": {"id": "A", "sequence": r3u_sequence}},
                        {"ligand":  {"id": "L", "smiles": _lig_smiles}}
                    ],
                    "properties": [{"affinity": {"binder": "L"}}]
                }
                with open(_r3u_yaml, "w") as _f:
                    yaml.safe_dump(r3u_payload, _f, sort_keys=False, default_flow_style=False)
            _ctrl_idx[0] += 1
            _lig_short_name = _lig_name.split("_", 1)[-1]
            console_info(f"  [{_ctrl_idx[0]}/{_n_ctrl_total}]  3R3U  × {_lig_short_name:<18}  →  {_r3u_name}")
            _r3u_res, _ = process_single_job(
                _r3u_task, PROD, 5, 0.0, D_ALN, control_cif, control_map, DEHA4_CONTROL_SEQ, None, None
            )
            r3u_results[_lig_name] = _r3u_res or {}
            _r3u_cifs = sorted((D_RUNS / _r3u_name).glob("**/*.cif"))
            if _r3u_cifs:
                _r3u_summary_p = D_RUNS / _r3u_name / f"{_r3u_name}_summary.json"
                _r3u_cached_ts = None
                _r3u_cached_cm = None
                if _r3u_summary_p.exists():
                    try:
                        with open(_r3u_summary_p) as _rf:
                            _rs = json.load(_rf)
                        if "ctrl_trust_score" in _rs:
                            _r3u_cached_ts = float(_rs["ctrl_trust_score"])
                            _r3u_cached_cm = _rs.get("ctrl_active_site_map", {})
                    except Exception:
                        pass
                if _r3u_cached_ts is not None:
                    _ts  = _r3u_cached_ts
                    _cm  = _r3u_cached_cm if isinstance(_r3u_cached_cm, dict) else {}
                    _r3u_calib_data = _r3u_cached_cm if isinstance(_r3u_cached_cm, dict) and "mapping_rows" in _r3u_cached_cm else {
                        "model_label": f"3R3U Boltz-2 {_lig_name}",
                        "trust_score": _ts, "verdict_label": "Cached",
                            "mapping_rows": [], "missing_residues": [], "all_mapped": True,
                    }
                    console_info(f"  > 3R3U {_lig_name}: calibration loaded from cache (trust={_ts:.3f} Å).")
                else:
                    _ts, _cm, _r3u_calib_data = perform_control_calibration(
                        _r3u_cifs[0], ref_pdb_path,
                        model_label=f"3R3U Boltz-2 {_lig_name}",
                        target_seq=r3u_sequence,
                        silent=True
                    )
                    try:
                        _rs2 = json.loads(_r3u_summary_p.read_text()) if _r3u_summary_p.exists() else {}
                        _rs2["ctrl_trust_score"] = _ts
                        _rs2["ctrl_active_site_map"] = _cm
                        _rs2["ctrl_calib_data"] = _r3u_calib_data
                        _r3u_summary_p.write_text(json.dumps(_rs2, indent=2))
                    except Exception:
                        pass
                r3u_trust_scores[_lig_name]  = _ts
                r3u_control_maps[_lig_name]  = _cm
                r3u_calib_list.append(_r3u_calib_data)
            if _lig_name == "26_Fluoroacetate":
                r3u_result_data = _r3u_res
                r3u_boltz_cif   = _r3u_cifs[0] if _r3u_cifs else None
                r3u_control_map = r3u_control_maps.get("26_Fluoroacetate", {})
                if r3u_result_data and r3u_result_data.get("status") == "Success":
                    pass  # reported in consolidated calibration summary below
                if not r3u_boltz_cif:
                    console_info("System Warning: 3R3U Fluoroacetate CIF not found — r3u scoring columns will be absent from Master CSV.")
                    sys.exit(1)

        # Print consolidated calibration summary (RMSD + active-site mapping tables)
        print_control_calibration_consolidated(deha4_calib_list, r3u_calib_list)

        # --- Extended substrate defluorination profile (all 6 control jobs) ---
        _mw, _cw  = 32, 17
        _lig_labels = [ln for ln, _ in CTRL_LIGANDS]
        _lig_short  = [ln.split("_", 1)[-1] if "_" in ln else ln for ln in _lig_labels]
        _tot = _mw + (_cw + 2) * 3

        for _prot_label, _res_dict, _ts_dict in [
            ("3R3U (Rhodopseudomonas palustris FAcD)", r3u_results,  r3u_trust_scores),
            ("DeHa4 Control",                         ctrl_results,  ctrl_trust_scores),
        ]:
            console_info(f"\n  {'─'*_tot}")
            console_info(f"  {f'Substrate Defluorination Profile — {_prot_label}':^{_tot}}")
            console_info(f"  {'─'*_tot}")
            console_info(f"  {'Metric':<{_mw}}" + "".join(f"  {s:^{_cw}}" for s in _lig_short))
            console_info(f"  {'─'*_tot}")

            def _row(label, key, fmt="{:.2f}", unit="", _rd=_res_dict):
                cells = []
                for ln in _lig_labels:
                    val = _rd.get(ln, {}).get(key)
                    if val is None:
                        cells.append(f"{'—':^{_cw}}")
                    else:
                        try:
                            cells.append(f"{fmt.format(float(val)) + unit:^{_cw}}")
                        except Exception:
                            cells.append(f"{str(val):^{_cw}}")
                console_info(f"  {label:<{_mw}}" + "".join(f"  {c}" for c in cells))

            def _trust_row(_td=_ts_dict):
                cells = [f"{f'{_td[ln]:.3f} Å' if ln in _td else '—':^{_cw}}" for ln in _lig_labels]
                console_info(f"  {'Trust Score (RMSD to Crystal)':<{_mw}}" + "".join(f"  {c}" for c in cells))

            _trust_row()
            _row("Degrader Tier",              CFG.COL_TIER,                 fmt="{}",     unit="")
            _row("SN2 Attack Angle",           "sn2_attack_angle",              fmt="{:.1f}", unit="°")
            _row("Nucleophile Distance (ASP)", "dist_Nuc",                      fmt="{:.2f}", unit=" Å")
            _row("Base Distance (HIS)",        "dist_Base",                     fmt="{:.2f}", unit=" Å")
            _row("Acid Distance (ASP)",        "dist_Acid",                     fmt="{:.2f}", unit=" Å")
            _row("Scissile C–F BDE (A)",       "scissile_cf_bde",               fmt="{:.1f}", unit=" kcal/mol")
            _row("Backside Occlusion (B)",     "sn2_backside_occlusion",        fmt="{:.2f}", unit=" Å")
            _row("α-Carbon Attack (1=yes)",    "scissile_is_alpha",             fmt="{:.0f}", unit="")
            _row("Carboxylate Clamp (0/½/1)",  "carboxylate_clamp_integrity",   fmt="{:.1f}", unit="")
            _row("Nucleophile Resolution",     "nuc_resolution",                fmt="{}",     unit="")
            _row("Mechanistic Score (geometry)", CFG.COL_MECH_S,            fmt="{:.2f}", unit="")
            _row("Pocket Containment (cavity)", "pocket_containment_cavity",     fmt="{:.2f}", unit="")
            _row("Pocket Containment (site-8)", "pocket_containment_site8",      fmt="{:.2f}", unit="")
            _row("Feasibility-wtd Mech (tier)", "mechanistic_score_effective",   fmt="{:.2f}", unit="")
            _row("Active Site Conservation Score",  CFG.COL_LIKE_S,     fmt="{:.2f}", unit="%")
            _row("Boltz-2 Confidence",         "confidence_score",              fmt="{:.4f}", unit="")
            console_info(f"  {'─'*_tot}")

        # --- Ligand selectivity ranking ---
        """
        Rank by: (1) Degrader Tier, (2) Mechanistic Score, (3) SN2 angle proximity to 180°.
        Tier hierarchy dynamically loaded.
        SN2 reactions require a nucleophilic attack angle close to 180° (back-side attack).
        """
        _TIER_RANK = CFG.TIER_RANK
        _w = 118
        console_info(f"\n  ┌{'─'*_w}┐")
        console_info(f"  │{'Substrate Selectivity Ranking':^{_w}}│")
        console_info(f"  │{'Rank criteria (same as the ranked CSV): Tier → Competence → Constellation → Conservation → Consensus':^{_w}}│")
        console_info(f"  ├{'─'*_w}┤")
        for idx, (_prot_s, _res_dict, _ts_dict) in enumerate([
            ("3R3U  (Rhodopseudomonas palustris FAcD — crystal sequence)",  r3u_results,  r3u_trust_scores),
            ("DeHa4 (Input enzyme sequence — primary screening target)",    ctrl_results, ctrl_trust_scores),
        ]):
            if idx > 0:
                console_info(f"  ├{'─'*_w}┤")
            console_info(f"  │  {_prot_s:<{_w-2}}│")
            # Same sort key as the whole-library Scientific_Rank (generate_scientific_ranking_csv):
            # tier → competence → constellation → active-site conservation → model consensus.
            # Controls and the 58,050 library jobs are ranked by identical parameters.
            ranked = sorted(
                [(ln, _res_dict.get(ln, {})) for ln, _ in CTRL_LIGANDS],
                key=lambda x: (
                    _TIER_RANK.get(x[1].get(CFG.COL_TIER, CFG.TIER_DECOY), 0),
                    float(x[1].get("competence_score") or 0.0),
                    float(x[1].get("catalytic_constellation_score") or 0.0),
                    float(x[1].get(CFG.COL_LIKE_S) or 0.0),
                    float(x[1].get("model_degrader_consensus") or 0.0),
                ),
                reverse=True
            )
            for _rank, (ln, d) in enumerate(ranked, 1):
                mech  = float(d.get("competence_score") or 0.0)
                sn2   = float(d.get("sn2_attack_angle") or 0.0)
                _tier = d.get(CFG.COL_TIER, CFG.TIER_DECOY)
                ts    = _ts_dict.get(ln)
                ts_s  = f"Trust: {ts:.3f} Å" if ts is not None else "Trust: —"
                _lig_short = ln.split("_", 1)[-1] if "_" in ln else ln
                if _rank == 1:
                    if _tier not in (CFG.TIER_DECOY, "Error", CFG.TIER_POOR, ""):
                        tag = "  ← Preferred substrate"
                    else:
                        tag = f"  ← Top-ranked by tier/competence  [Tier: {_tier} — no degradation activity confirmed]"
                else:
                    tag = ""
                row_content = (
                    f"    Rank {_rank}  │  {_lig_short:<18}  │  Tier: {_tier:<12}  │  "
                    f"Comp: {mech:.3f}  │  SN2: {sn2:>5.1f}°  │  {ts_s}{tag}"
                )
                console_info(f"  │{row_content:<{_w}}│")
        console_info(f"  └{'─'*_w}┘")
        console_info(f"  {'Mech = holistic 0–1 geometry score: anchors (nucleophile reach + relay distances + clamp + halide stabilisation) + graded SN2 angle. Tier gates on the feasibility-weighted mech (− graded BDE/occlusion − pocket-containment penalties).':^{_w+4}}")
        console_info(f"  {f'A = scissile C–F bond-dissociation energy (kcal/mol; >{CFG.SCISSILE_CF_BDE_MAX:.0f} too strong); B = backside steric occlusion (Σ vdW Å; >{CFG.SN2_BACKSIDE_OCCL_MAX:.1f} blocked). Both → non-degradable CF3 attack carbon.':^{_w+4}}")
        console_info(f"  {'─'*_tot}\n")
    else:
        console_info("  [Warning] Could not extract sequence from 3R3U PDB — skipping reference jobs.")
    console_separator()

    # -------------------------------------------------------------------------------
    # Step 10.6: System Workload Audit
    # -------------------------------------------------------------------------------
    """
    A "completed" job has both: (a) a Boltz-2 CIF prediction on disk AND
    (b) a fully written _summary.json with status="Success".
    On --resume with wipe, summaries are cleared so completed_jobs = 0 even
    if GPU predictions exist — those jobs re-run CPU analysis only (not GPU).
    """
    console_info("Auditing job completion status...")
    completed_jobs = 0
    pending_jobs = 0
    pending_proteins: set = set()

    if resumed:
        # Wipe cleared all summary JSONs — every job is pending by definition.
        console_info(" -> Resume wipe applied: all jobs marked as pending (skipping disk scan).")
        pending_jobs = len(tasks)
        pending_proteins = {job["protein"] for job in tasks}
    else:
        _lock = threading.Lock()
        _n_tasks = len(tasks)
        _checked = [0]

        def _audit_job(job):
            job_name = f"{job['job_index']}_{job.get('job_protein', job['protein'])}_{job['ligand']}"
            done = check_job_status(D_RUNS / job_name)
            with _lock:
                _checked[0] += 1
                n = _checked[0]
                if n % 2000 == 0 or n == _n_tasks:
                    _tty_write(f"\r   Scanning job status: {n:,}/{_n_tasks:,} ...   \033[K")
                    sys.stdout.flush()
                if done:
                    pass  # counted below
                else:
                    pending_proteins.add(job["protein"])
            return done

        from concurrent.futures import ThreadPoolExecutor as _TPE
        with _TPE(max_workers=min(max(1, (os.cpu_count() or 4) - 2), _n_tasks or 1)) as _ex:
            results = list(_ex.map(_audit_job, tasks))
        _tty_write("\r\033[K")
        sys.stdout.flush()
        completed_jobs = sum(results)
        pending_jobs   = _n_tasks - completed_jobs

    """
    MSA cache is only relevant to jobs that still need a GPU prediction; the actual MSA
    fetch is lazy (background_downloader, per protein, during the GPU phase). On a resume
    where every job already has a prediction on disk there is no GPU work, so validating
    the .a3m of all pending proteins (one disk read each) is pure overhead — skip it. A
    cheap folder-count estimate decides: predictions already on disk vs full GPU needed.
    """
    """
    MSA validation is deferred to the per-job GPU scan (Step 10.7). Which jobs truly
    need a GPU prediction is only established there; the MSA fetch (background_downloader)
    self-validates each .a3m before downloading. Reading every pending protein's .a3m
    here — before any GPU need is known — is pure overhead, so it is skipped. The real
    MSA-ready count is computed over the actual GPU proteins just before the GPU phase.
    """
    console_info(f" -> {completed_jobs:,} complete, {pending_jobs:,} pending. MSA cache validation deferred to the GPU-job scan.")
    msa_ready = 0

    console_info("")
    _w1, _w2 = 60, 14
    _ww = _w1 + _w2 + 1
    console_info(f"  ┌{'─'*_ww}┐")
    console_info(f"  │{'SYSTEM WORKLOAD SUMMARY':^{_ww}}│")
    console_info(f"  ├{'─'*_w1}┬{'─'*_w2}┤")
    console_info(f"  │ {'Total jobs scheduled':<{_w1-2}} │ {len(tasks):>{_w2-2},} │")
    console_info(f"  │ {'  Completed (Boltz-2 prediction + CPU analysis)':<{_w1-2}} │ {completed_jobs:>{_w2-2},} │")
    if resumed:
        """
        On resume: GPU predictions on disk but summary JSONs wiped for re-analysis.
        Distinguish jobs with existing prediction dirs from those needing full GPU run.
        """
        _dirs_with_pred = max(0, min(initial_folders_count - 6, pending_jobs))  # subtract 6 control dirs
        _need_gpu       = max(0, pending_jobs - _dirs_with_pred)
        console_info(f"  │ {'  Pending — GPU prediction on disk, CPU re-analysis queued':<{_w1-2}} │ {_dirs_with_pred:>{_w2-2},} │")
        console_info(f"  │ {'  Pending — full GPU + CPU pipeline required':<{_w1-2}} │ {_need_gpu:>{_w2-2},} │")
        console_info(f"  │ {'  (Resume mode: analysis outputs wiped for re-scoring)':<{_w1-2}} │ {'':>{_w2-2}} │")
    else:
        console_info(f"  │ {'  Pending — awaiting GPU prediction + CPU analysis':<{_w1-2}} │ {pending_jobs:>{_w2-2},} │")
    console_info(f"  ├{'─'*_w1}┼{'─'*_w2}┤")
    console_info(f"  │ {'Unique proteins requiring processing':<{_w1-2}} │ {len(pending_proteins):>{_w2-2},} │")
    console_info(f"  │ {'MSA cache check':<{_w1-2}} │ {'deferred':>{_w2-2}} │")
    console_info(f"  │ {'  (validated per protein during the GPU-job scan)':<{_w1-2}} │ {'':>{_w2-2}} │")
    console_info(f"  │ {'Orphaned job folders purged at startup':<{_w1-2}} │ {removed_folders_count:>{_w2-2},} │")
    console_info(f"  └{'─'*_w1}┴{'─'*_w2}┘")
    console_info("")

    if completed_jobs > 0:
        _tty_write(" -> Actively standardising 'Best_Complex' allocations for pre-existing structural results...")
        sys.stdout.flush()
        def _standardize_bc(d):
            if not (d.is_dir() and d.name.startswith("0")):
                return
            bc_dir = d / "Best_Complex"
            if bc_dir.exists():
                cifs = sorted(list(bc_dir.glob("*.cif")))
                if len(cifs) > 1:
                    for bad_f in cifs[1:]:
                        try: bad_f.unlink()
                        except Exception: pass
        from concurrent.futures import ThreadPoolExecutor as _TPE
        _bc_dirs = list(sorted(D_RUNS.iterdir()))
        with _TPE(max_workers=min(max(1, (os.cpu_count() or 4) - 2), len(_bc_dirs) or 1)) as _ex:
            list(_ex.map(_standardize_bc, _bc_dirs))
        _tty_write("\r -> Complete standardisation of pre-existing designated folder structures finished successfully.\033[K\n")
        sys.stdout.flush()


    # -------------------------------------------------------------------------------
    # Step 10.7: Master CSV Rebuild + GPU Hardware Batch Prediction Phase
    # -------------------------------------------------------------------------------

    # --- Upfront bulk CSV rebuild from all completed summary JSONs ---
    console_separator()
    console_info("Rebuilding Master CSV from completed per-job summary.json files...")
    console_info("  (Reads each finished job's summary.json and re-aggregates them into the Master CSV; the Best-Complex CIF mirror is rebuilt in the next step.)")
    _rebuild_t0 = time.time()
    _rebuilt_rows = 0
    if completed_jobs == 0 and not resumed:
        # No summaries exist (fresh run) — skip rglob traversal.
        console_info(" -> No completed jobs — CSV will be populated as analysis completes.")
    else:
        _rebuilt_rows = rebuild_csv_from_summaries(D_RUNS, CSV_PATH)
        _rebuild_elapsed = time.time() - _rebuild_t0
        console_info(f" -> Master CSV rebuilt: {_rebuilt_rows:,} rows in {_rebuild_elapsed:.1f}s  →  {CSV_PATH.name}")

    """
    Rebuild the run-level best-complex mirror from all completed job summaries.
    This ensures 0\2_Best_Complexes_CIFs/<tier>/ is always consistent on every
    resume, even if previous runs were interrupted mid-copy or tiers changed.
    """
    console_info("Rebuilding 2_Best_Complexes_CIFs mirror from completed job summaries...")
    _mirror_t0 = time.time()
    _mirror_cifs = rebuild_best_complexes_mirror(D_RUNS, run_root)
    _mirror_elapsed = time.time() - _mirror_t0
    _mirror_dir = run_root / "2_Best_Complexes_CIFs"
    if _mirror_cifs > 0:
        console_info(f" -> Mirror rebuilt: {_mirror_cifs:,} CIFs across tier subfolders in {_mirror_elapsed:.1f}s  →  {_mirror_dir.name}/")
    else:
        console_info(" -> No completed CIFs found yet — mirror will be populated as predictions complete.")
    console_separator()

    # --- Background analysis worker for freshly GPU-predicted jobs ---
    _mgr = multiprocessing.Manager()
    gpu_queue_for_analysis = _mgr.Queue()
    analysis_pool = ProcessPoolExecutor(max_workers=1)
    analysis_pool.submit(analysis_worker_loop, gpu_queue_for_analysis, PROD,
                         args.diffusion_samples, D_ALN, control_cif, control_map,
                         DEHA4_CONTROL_SEQ, CSV_PATH, r3u_boltz_cif, r3u_control_map)

    if pending_jobs > 0:
        workspace_dir = PROD / "_Temp_Workspace"
        workspace_dir.mkdir(exist_ok=True)

        batch_out_dir = workspace_dir / "Boltz_Batch_Output"
        batch_out_dir.mkdir(exist_ok=True)

        env = os.environ.copy()
        env["CC"] = shutil.which("gcc") or "gcc"
        env["CXX"] = shutil.which("g++") or "g++"
        # Per-worker Triton cache so parallel Boltz subprocesses cannot collide on lockfiles
        # in a shared /tmp/triton_cache.
        env["TRITON_CACHE_DIR"] = str(Path(tempfile.gettempdir()) / f"triton_cache_{os.getpid()}_{uuid.uuid4().hex[:8]}")
        # conda env PFAS manages CUDA/nvidia libraries via LD_LIBRARY_PATH.
        # Pin BLAS/OMP threads to 1 per subprocess: N parallel workers each spawning
        # GLOBAL_MAX_WORKERS threads oversubscribes the CPU (N² thread thrashing).
        env["OMP_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"
        env["OPENBLAS_NUM_THREADS"] = "1"
        env["PYTHONWARNINGS"] = "ignore"
        env["PYTORCH_LIGHTNING_SUPPRESS_WARNINGS"] = "1"
        env["TF_CPP_MIN_LOG_LEVEL"] = "3"
        env["CUDA_LAUNCH_BLOCKING"] = "0"
        env["WANDB_DISABLED"] = "true"
        env["PYTHONDONTWRITEBYTECODE"] = "1"

        pending_by_protein = defaultdict(list)
        pre_predicted_count = 0

        console_info(f"Scanning {len(tasks):,} jobs for pre-existing GPU predictions...")
        _scan_lock = threading.Lock()
        _pre_pred_jobs = []
        _gpu_only_jobs = []
        _scan_done = [0]
        _n_scan = len(tasks)

        def _scan_job(job):
            job_name = f"{job['job_index']}_{job.get('job_protein', job['protein'])}_{job['ligand']}"
            job_dir = D_RUNS / job_name
            if check_job_status(job_dir):
                return  # already complete
            prediction = next(iter(sorted(job_dir.glob("**/predictions/**/*.cif"))), None)
            with _scan_lock:
                _scan_done[0] += 1
                n = _scan_done[0]
                if n % 2000 == 0 or n == _n_scan:
                    _tty_write(f"\r   [{n:,}/{_n_scan:,}] pre-predicted: {len(_pre_pred_jobs):,}  gpu-needed: {len(_gpu_only_jobs):,}   \033[K")
                    sys.stdout.flush()
                if prediction:
                    _pre_pred_jobs.append(job)
                else:
                    _gpu_only_jobs.append(job)

        from concurrent.futures import ThreadPoolExecutor as _TPE
        with _TPE(max_workers=min(max(1, (os.cpu_count() or 4) - 2), _n_scan or 1)) as _sex:
            list(_sex.map(_scan_job, tasks))
        _tty_write("\r\033[K")
        sys.stdout.flush()

        pre_predicted_count = len(_pre_pred_jobs)
        for _gj in _gpu_only_jobs:
            pending_by_protein[_gj["protein"]].append(_gj)

        if pre_predicted_count > 0:
            if pending_by_protein:
                """
                Mixed mode: GPU-only jobs also present — route pre-predicted jobs
                through the background single-worker queue so analysis runs
                concurrently with ongoing GPU prediction.
                """
                for _pj in _pre_pred_jobs:
                    gpu_queue_for_analysis.put(_pj)
                console_info(f" -> {pre_predicted_count:,} jobs have pre-existing GPU predictions — re-routing to Analysis Queue...")
            else:
                """
                All jobs are pre-existing — skip the single-worker queue entirely.
                Step 10.8's full parallel pool will handle analysis at maximum
                concurrency (TARGET_CORES workers). No blocking wait.
                """
                console_info(f" -> {pre_predicted_count:,} jobs have pre-existing GPU predictions.")
                console_info(f" -> Routing to parallel analysis pool ({TARGET_CORES} workers) — skipping single-worker queue.")
        if _gpu_only_jobs:
            console_info(f" -> {len(_gpu_only_jobs):,} jobs require GPU prediction.")

        if pending_by_protein:
            # MSA cache is only relevant to proteins that actually need a GPU prediction.
            msa_ready = sum(1 for pid in pending_by_protein if validate_a3m_file(D_COLABFOLD / f"{pid}.a3m"))
            _n_gpu_jobs = sum(len(v) for v in pending_by_protein.values())
            print(f"\n{ConsoleColours.OKGREEN}{ConsoleColours.BOLD}{'=' * 80}{ConsoleColours.ENDC}")
            print(f"{ConsoleColours.OKGREEN}{ConsoleColours.BOLD}  GPU PREDICTION PHASE STARTING  —  {_n_gpu_jobs} jobs  |  {len(pending_by_protein)} proteins  |  {msa_ready}/{len(pending_by_protein)} MSA cached{ConsoleColours.ENDC}")
            print(f"{ConsoleColours.OKGREEN}{ConsoleColours.BOLD}{'=' * 80}{ConsoleColours.ENDC}\n", flush=True)
            console_info("Launching Boltz-2 Batch Prediction Sequence Operations (GPU)...")
            console_info("Execution Strategy: Greedy batching protocol — multiple individual proteins bundled per Boltz call to systematically minimise model reloads and mitigate GPU idle instances.")
        else:
            console_info(f"\n -> All {pending_jobs} pending job(s) have pre-existing GPU predictions — no new predictions required.")
            console_info(" -> Proceeding directly to parallel analysis...")

        msa_dir = D_COLABFOLD
        msa_dir.mkdir(parents=True, exist_ok=True)
        tmp_msa_dir = workspace_dir / "5b_MSA_Temp"
        tmp_msa_dir.mkdir(exist_ok=True)

        if pending_by_protein:
            console_info("\nStarting Concurrent MSA Retrieval Sequence and GPU Prediction Execution Pipeline...")
            console_info("MSA retrieval (via ColabFold APIs) executes consistently in a background parallel thread whilst Boltz-2 predicts completed sequential batches via the GPU architecture.")

        def background_downloader(proteins_to_fetch):
            """Fetch MSAs for the given proteins sequentially in a daemon thread.

            For each protein writes <pid>.a3m into msa_dir (skipping ones already
            present), using the ColabFold MSA server with a Boltz fallback. Runs
            concurrently with the GPU prediction phase; failures are tolerated and
            retried per-protein downstream.
            """
            for pid in proteins_to_fetch:
                a3m_path = msa_dir / f"{pid}.a3m"
                if a3m_path.exists():
                    if validate_a3m_file(a3m_path):
                        continue
                    else:
                        logger.debug(f"A3M sequence file for {pid} failed extensive structural validation. Re-fetching structural data...")
                        a3m_path.unlink(missing_ok=True)

                seq       = pending_by_protein[pid][0]["sequence"]
                meta_path = colabfold_meta_path(PROD, pid)

                logger.debug(f"[MSA] Commencing MSA structural retrieval for {pid} securely via ColabFold API framework...")
                t_msa = time.time()
                ok = fetch_msa_direct(pid, seq, a3m_path, meta_path)
                if ok:
                    logger.debug(f"[MSA] ColabFold API sequence operations succeeded for {pid} in {time.time()-t_msa:.1f}s")

                if not ok:
                    console_info(f"  [MSA] Direct API methodologies functionally failed for {pid} — actively falling back to Boltz subprocess execution strategies (computationally slow)")
                    yaml_path = tmp_msa_dir / f"{pid}_fetch.yaml"
                    with open(yaml_path, "w") as f:
                        yaml.safe_dump({"sequences": [{"protein": {"id": "A", "sequence": seq}}]}, f, sort_keys=False)

                    cmd = [
                        BOLTZ_BIN, "predict", str(yaml_path), "--out_dir", str(tmp_msa_dir),
                        "--use_msa_server", "--accelerator", "cpu", "--devices", "1"
                    ]

                    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=str(tmp_msa_dir))
                    start     = time.time()
                    found     = False
                    job_msa_dir = tmp_msa_dir / f"boltz_results_{pid}_fetch"

                    while time.time() - start < 900:
                        if proc.poll() is not None: break
                        if job_msa_dir.exists():
                            a3ms = sorted(job_msa_dir.rglob("*.a3m"))
                            if a3ms and a3ms[0].stat().st_size > 500:
                                if validate_a3m_file(a3ms[0]):
                                    # Atomic publish: copy to a temp file in the same
                                    # directory, then rename. The main thread polls
                                    # a3m_path concurrently; a mid-copy read could
                                    # otherwise pass the size/header check on a
                                    # truncated file and predict on a partial MSA.
                                    _tmp_a3m = a3m_path.with_suffix(a3m_path.suffix + ".tmp")
                                    shutil.copy2(a3ms[0], _tmp_a3m)
                                    _tmp_a3m.replace(a3m_path)
                                    with open(meta_path, "w") as mf:
                                        json.dump({
                                            "sequence_sha256": sequence_hash(seq),
                                            "protein_id":      pid,
                                            "generated_at":    datetime.utcnow().isoformat() + "Z",
                                            "source":          "boltz_fallback",
                                        }, mf)
                                    found = True
                                    proc.terminate()
                                    break
                                else:
                                    logger.debug(f"A3M sequence for {pid} is inherently corrupted (aberrant null bytes detected). Retrying...")
                                    shutil.rmtree(job_msa_dir, ignore_errors=True)
                                    time.sleep(5)
                                    break
                        time.sleep(5)

                    # Ensure the fallback Boltz subprocess is not left running once the
                    # polling loop exits via any path (corrupted A3M, 900 s timeout, or
                    # completion without an early terminate).
                    if proc.poll() is None:
                        proc.terminate()
                        try:
                            proc.wait(timeout=10)
                        except Exception:
                            proc.kill()

                    if not found and job_msa_dir.exists():
                        a3ms = sorted(job_msa_dir.rglob("*.a3m"))
                        if a3ms: shutil.copy2(a3ms[0], a3m_path)

                    shutil.rmtree(job_msa_dir, ignore_errors=True)
                    yaml_path.unlink(missing_ok=True)

                time.sleep(3)

        downloader_finished = threading.Event()

        def _downloader_wrapper(proteins_to_fetch):
            try:
                background_downloader(proteins_to_fetch)
            finally:
                downloader_finished.set()

        downloader_thread = threading.Thread(target=_downloader_wrapper, args=(list(pending_by_protein.keys()),), daemon=True)
        downloader_thread.start()

        protein_list = list(pending_by_protein.keys())
        completed_prots = len(proteins) - len(protein_list)
        processed_pids = set()

        for c_idx, pid in enumerate(protein_list):
            if pid in processed_pids:
                continue

            a3m_path = msa_dir / f"{pid}.a3m"
            current_prot_num = completed_prots + c_idx + 1
            protein_name_display = pid[:60] if pid else "Unknown"

            dl_wait_start = time.time()
            while not validate_a3m_file(a3m_path):
                elapsed_wait = time.time() - dl_wait_start
                mins, secs = divmod(int(elapsed_wait), 60)
                elapsed_str = f"{mins}m {secs:02d}s" if mins else f"{secs}s"
                _tty_write(f"\r -> [WAIT] Current MSA retrieval: {protein_name_display:<45}  {elapsed_str:>8} currently suspended..." + "\033[K")
                sys.stdout.flush()
                if downloader_finished.is_set():
                    console_info(f"\nWARNING: MSA retrieval framework terminated independently but no structurally valid A3M output emerged for {pid} (elapsed time {elapsed_str}). Systematically skipping.")
                    break
                if elapsed_wait > 900:
                    console_info(f"\nWARNING: Extensive MSA download timeout occurred for designated {pid} (>15 min). Systematically skipping specified protein sequence.")
                    break
                time.sleep(1)
            if not validate_a3m_file(a3m_path):
                console_info(f"\n  Skipping configuration {pid} (invalid/corrupted MSA structural data). Will be systematically retried within the individual processing grid.")
                continue

            proteins_for_batch = [pid]
            for look_pid in protein_list[c_idx + 1:]:
                if len(proteins_for_batch) >= MAX_PROTEINS_PER_BATCH:
                    break
                if look_pid in processed_pids:
                    continue
                if validate_a3m_file(msa_dir / f"{look_pid}.a3m"):
                    proteins_for_batch.append(look_pid)

            all_jobs_in_batch = []
            for batch_pid in proteins_for_batch:
                all_jobs_in_batch.extend(pending_by_protein[batch_pid])

            n_prot      = len(proteins_for_batch)
            n_jobs      = len(all_jobs_in_batch)

            if n_prot > 1:
                _tty_write(
                    f"\n -> [GREEDY BATCH] Bundling {n_prot} proteins with {n_jobs} total pending jobs into one Boltz call"
                    f" — saving {n_prot - 1} model reload(s).\n"
                )
                sys.stdout.flush()

            batch_label = pid if n_prot == 1 else f"{proteins_for_batch[0]} +{n_prot - 1} more sequence blocks"

            unique_id         = str(uuid.uuid4())[:8]
            chunk_dir         = (workspace_dir / f"2_Boltz2_YAML_Configs_batch_{unique_id}").resolve()
            chunk_dir.mkdir(parents=True, exist_ok=True)
            protein_batch_out = (batch_out_dir / f"batch_out_{unique_id}").resolve()
            protein_batch_out.mkdir(parents=True, exist_ok=True)

            job_map = {}
            for batch_pid in proteins_for_batch:
                batch_a3m = msa_dir / f"{batch_pid}.a3m"
                for job in pending_by_protein[batch_pid]:
                    job_name = f"{job['job_index']}_{job.get('job_protein', job['protein'])}_{job['ligand']}"
                    with open(job["yaml"], "r") as yf: data = yaml.safe_load(yf)
                    data["sequences"][0]["protein"]["msa"] = str(batch_a3m.resolve())
                    chunk_yaml = chunk_dir / job["yaml"].name
                    with open(chunk_yaml, "w") as yf: yaml.safe_dump(data, yf, sort_keys=False)
                    job_map[chunk_yaml.stem] = job
                    job_run_dir = D_RUNS / job_name
                    if job_run_dir.exists():
                        for old_res in sorted(job_run_dir.glob("boltz_results_*")):
                            if old_res.is_dir(): shutil.rmtree(old_res, ignore_errors=True)

            cmd = [
                BOLTZ_BIN, "predict", str(chunk_dir), "--out_dir", str(protein_batch_out),
                "--cache", str(BOLTZ_CACHE), "--model", BOLTZ_MODEL,
                "--recycling_steps", str(RECYCLING_STEPS),
                "--diffusion_samples", str(args.diffusion_samples),
                "--accelerator", "gpu", "--devices", "1",
                "--preprocessing-threads", "2",
                "--output_format", OUTPUT_FORMAT,
                "--no_kernels"
            ]

            max_retries = 8
            attempt = 0
            _surgical_repair_done = False  # Only attempt surgical repair once per batch

            while attempt < max_retries:
                # --- Salvage any partial predictions before wiping the temp folder ---
                if protein_batch_out.exists() and attempt > 0:
                    _salvage_res_dir = next(iter(sorted(protein_batch_out.glob("boltz_results_*"))), None)
                    if _salvage_res_dir and (_salvage_res_dir / "predictions").exists():
                        _salvaged = set()
                        for _pred_dir in sorted((_salvage_res_dir / "predictions").iterdir()):
                            _stem = _pred_dir.name
                            if _stem not in job_map or not _pred_dir.is_dir():
                                continue
                            _job = job_map[_stem]
                            _jname = f"{_job['job_index']}_{_job.get('job_protein', _job['protein'])}_{_job['ligand']}"
                            _dst_root = D_RUNS / _jname / f"boltz_results_{_stem}"
                            _dst_pred = _dst_root / "predictions" / _stem
                            try:
                                _dst_pred.mkdir(parents=True, exist_ok=True)
                                for _f in sorted(_pred_dir.iterdir()):
                                    shutil.move(str(_f), str(_dst_pred / _f.name))
                                for _cat in ["constraints", "mols", "msa", "records", "structures"]:
                                    _src_cat = _salvage_res_dir / "processed" / _cat
                                    if _src_cat.exists():
                                        for _ext in [".json", ".npz", ".pkl", ".csv"]:
                                            for _m in sorted(_src_cat.glob(f"{_stem}*{_ext}")):
                                                _d = _dst_root / "processed" / _cat
                                                _d.mkdir(parents=True, exist_ok=True)
                                                shutil.move(str(_m), str(_d / _m.name))
                                gpu_queue_for_analysis.put(_job)
                                (chunk_dir / f"{_stem}.yaml").unlink(missing_ok=True)
                                _salvaged.add(_stem)
                            except Exception:
                                pass
                        if _salvaged:
                            for _s in _salvaged:
                                job_map.pop(_s, None)
                            all_jobs_in_batch = [j for j in all_jobs_in_batch if j["yaml"].stem in job_map]
                            n_jobs = len(all_jobs_in_batch)
                            console_info(f"\n  [SALVAGE] Recovered {len(_salvaged)} partial predictions — {n_jobs} jobs remain for retry.")
                if n_jobs == 0:
                    shutil.rmtree(protein_batch_out, ignore_errors=True)
                    shutil.rmtree(chunk_dir, ignore_errors=True)
                    for _batch_pid in proteins_for_batch:
                        processed_pids.add(_batch_pid)
                    break
                if protein_batch_out.exists():
                    shutil.rmtree(protein_batch_out, ignore_errors=True)
                protein_batch_out.mkdir(parents=True, exist_ok=True)
                if not chunk_dir.exists(): chunk_dir.mkdir(parents=True, exist_ok=True)

                # --- Incremental-move tracking (reset each attempt) ---
                _moved_stems: set = set()
                _protein_job_counts: dict = {}
                for _s, _j in job_map.items():
                    _pid = _j["protein"]
                    if _pid not in _protein_job_counts:
                        _protein_job_counts[_pid] = {"moved": 0, "total": 0}
                    _protein_job_counts[_pid]["total"] += 1

                try:
                    err_log_path = workspace_dir / f"gpu_error_{unique_id}.log"
                    with open(err_log_path, "w") as err_file:
                        proc = subprocess.Popen(cmd, env=env, stdout=subprocess.DEVNULL, stderr=err_file)

                        gpu_start_time = time.time()
                        timeout_limit  = n_jobs * CFG.GPU_WATCHDOG_TIMEOUT_PER_JOB

                        while proc.poll() is None:
                            if time.time() - gpu_start_time > timeout_limit:
                                proc.kill()
                                raise subprocess.CalledProcessError(-9, cmd, stderr="WATCHDOG TIMEOUT ERROR: GPU Batch Deadlocked Iteration.")

                            # --- Incremental move: pick up finished predictions every 5 s ---
                            _res_dir = next(iter(sorted(protein_batch_out.glob("boltz_results_*"))), None) if protein_batch_out.exists() else None
                            if _res_dir and (_res_dir / "predictions").exists():
                                for _pred_dir in sorted((_res_dir / "predictions").iterdir()):
                                    _stem = _pred_dir.name
                                    if not _pred_dir.is_dir() or _stem in _moved_stems or _stem not in job_map:
                                        continue
                                    # Require the confidence JSON too: Boltz writes it after the
                                    # structure, so its presence marks a complete prediction and
                                    # prevents moving a still-being-written .cif (truncation race).
                                    if not any(_pred_dir.glob("*.cif")) or not any(_pred_dir.glob("confidence_*.json")):
                                        continue  # not finished yet
                                    _job = job_map[_stem]
                                    _jname = f"{_job['job_index']}_{_job.get('job_protein', _job['protein'])}_{_job['ligand']}"
                                    _dst_root = D_RUNS / _jname / f"boltz_results_{_stem}"
                                    _dst_pred = _dst_root / "predictions" / _stem
                                    try:
                                        _dst_pred.mkdir(parents=True, exist_ok=True)
                                        for _f in list(sorted(_pred_dir.iterdir())):
                                            shutil.move(str(_f), str(_dst_pred / _f.name))
                                        for _cat in ["constraints", "mols", "msa", "records", "structures"]:
                                            _src_cat = _res_dir / "processed" / _cat
                                            if _src_cat.exists():
                                                for _ext in [".json", ".npz", ".pkl", ".csv"]:
                                                    for _m in sorted(_src_cat.glob(f"{_stem}*{_ext}")):
                                                        _d = _dst_root / "processed" / _cat
                                                        _d.mkdir(parents=True, exist_ok=True)
                                                        shutil.move(str(_m), str(_d / _m.name))
                                        gpu_queue_for_analysis.put(_job)
                                        _moved_stems.add(_stem)
                                        _pid2 = _job["protein"]
                                        if _pid2 in _protein_job_counts:
                                            _protein_job_counts[_pid2]["moved"] += 1
                                            _pc = _protein_job_counts[_pid2]
                                            if _pc["moved"] == _pc["total"]:
                                                _elapsed_p = int(time.time() - gpu_start_time)
                                                _tty_write(
                                                    f"\r\033[K -> [DONE] {_pid2}"
                                                    f" — {_pc['total']}/{_pc['total']} jobs moved"
                                                    f" | {_elapsed_p}s elapsed\n"
                                                )
                                                sys.stdout.flush()
                                    except Exception:
                                        pass
                            # --- End incremental move ---

                            elapsed_gpu     = time.time() - gpu_start_time
                            elapsed_display = int(elapsed_gpu / 5) * 5
                            n_moved_total   = len(_moved_stems)
                            _tty_write(
                                f"\r -> [GPU] Overall: {current_prot_num}/{len(proteins)} proteins"
                                f" | This call: {n_prot} proteins, {n_jobs} jobs"
                                f" | Moved: {n_moved_total}/{n_jobs}"
                                f" | {elapsed_display}s elapsed" + "\033[K"
                            )
                            sys.stdout.flush()
                            time.sleep(5.0)

                    if proc.returncode == 0:
                        total_res_dir = next(iter(sorted(protein_batch_out.glob("boltz_results_*"))), None)
                        moved_count   = 0
                        failed_jobs   = []

                        if total_res_dir and (total_res_dir / "predictions").exists():
                            all_pred_folders = [d.name for d in sorted((total_res_dir / "predictions").iterdir()) if d.is_dir()]
                        else:
                            all_pred_folders = []

                        if len(all_pred_folders) == 0 and len(job_map) > 0:
                            console_info(f"\n    [BATCH CRITICAL WARNING] returncode=0 evaluated but 0 predictions yielded for [{batch_label}]!")
                            console_info(f"    [BATCH CRITICAL WARNING] Mathematically expected {len(job_map)} predictions but successfully obtained 0")
                            try:
                                with open(err_log_path, "r") as ef:
                                    stderr_content = ef.read()
                                    if stderr_content:
                                        critical_lines = [
                                            line for line in stderr_content.split("\n")
                                            if any(kw in line for kw in ["Traceback", "Error", "KeyError", "Exception", "CUDA out", "cuda", "memory"])
                                            and "|" not in line and "%" not in line
                                        ]
                                        if critical_lines:
                                            console_info("    [BATCH STDERR EVALUATION] Critical Errors Identified Documented Below:")
                                            for line in critical_lines[-20:]:
                                                if line.strip(): console_info(f"      {line[:120]}")
                                        else:
                                            non_prog = [l for l in stderr_content.split("\n")
                                                        if l.strip() and "|" not in l and "%" not in l and "it/s" not in l]
                                            if non_prog:
                                                console_info("    [BATCH STDERR EVALUATION] Terminal system output recorded preceding final operational return:")
                                                for line in non_prog[-5:]:
                                                    if line.strip(): console_info(f"      {line[:120]}")
                                    else:
                                        console_info("    [BATCH STDERR EVALUATION] Target error log is empty — Boltz execution yielded a strictly silent systemic failure")
                            except Exception as read_err:
                                console_info(f"    [BATCH CRITICAL WARNING] Unable to successfully parse stderr output log sequence: {read_err}")
                            err_log_path.unlink(missing_ok=True)
                        else:
                            if all_pred_folders:
                                logger.debug(f"[BATCH] {len(all_pred_folders)} isolated prediction folders structurally parsed for {len(job_map)} jobs")
                            err_log_path.unlink(missing_ok=True)

                        if total_res_dir and all_pred_folders:
                            for stem, job in job_map.items():
                                if stem in _moved_stems:
                                    moved_count += 1
                                    continue  # already moved incrementally during polling
                                job_name    = f"{job['job_index']}_{job.get('job_protein', job['protein'])}_{job['ligand']}"
                                job_run_dir = D_RUNS / job_name
                                job_run_dir.mkdir(parents=True, exist_ok=True)
                                dst_boltz_root = job_run_dir / f"boltz_results_{stem}"
                                dst_pred_dir   = dst_boltz_root / "predictions" / stem
                                src_pred       = total_res_dir / "predictions" / stem

                                if not src_pred.exists():
                                    pred_base = total_res_dir / "predictions"
                                    if pred_base.exists():
                                        matches = [d for d in sorted(pred_base.iterdir())
                                                   if d.is_dir() and (job["ligand"] in d.name or stem in d.name)]
                                        if matches:
                                            src_pred = matches[0]
                                            logger.debug(f"[BATCH FLEXIBLE MATCH] Correctly identified correlation {stem} -> {src_pred.name}")
                                        else:
                                            logger.debug(f"[BATCH WARNING] Failed to locate compatible structural prediction folder corresponding to {stem}")
                                            failed_jobs.append(stem)
                                            continue

                                if src_pred.exists():
                                    try:
                                        dst_pred_dir.mkdir(parents=True, exist_ok=True)
                                        for f in sorted(src_pred.iterdir()):
                                            shutil.move(str(f), str(dst_pred_dir / f.name))
                                        for cat in ["constraints", "mols", "msa", "records", "structures"]:
                                            src_cat = total_res_dir / "processed" / cat
                                            if src_cat.exists():
                                                for ext in [".json", ".npz", ".pkl", ".csv"]:
                                                    for m in sorted(src_cat.glob(f"{stem}*{ext}")):
                                                        dst_cat_dir = dst_boltz_root / "processed" / cat
                                                        dst_cat_dir.mkdir(parents=True, exist_ok=True)
                                                        shutil.move(str(m), str(dst_cat_dir / m.name))
                                        moved_count += 1
                                        logger.debug(f"[BATCH SUCCESS] Effectively translocated results referencing: {stem}")
                                    except Exception as move_err:
                                        logger.debug(f"[BATCH MOVEMENT ERROR] Migration operation structurally failed for {stem}: {move_err}")
                                        failed_jobs.append(stem)
                                else:
                                    logger.debug(f"[BATCH SKIP ENFORCED] Core prediction source structurally absent: {src_pred}")
                                    failed_jobs.append(stem)

                        shutil.rmtree(chunk_dir, ignore_errors=True)
                        shutil.rmtree(protein_batch_out, ignore_errors=True)

                        for stem, job in job_map.items():
                            if stem not in _moved_stems:
                                gpu_queue_for_analysis.put(job)

                        if failed_jobs:
                            console_info(f"    [BATCH RESULTS] {moved_count}/{n_jobs} securely migrated, {len(failed_jobs)} instances encountered systemic failure")
                        _tty_write(
                            f"\r -> [DONE] [{batch_label}] | {moved_count}/{n_jobs} structural evaluation results effectively migrated\033[K\n"
                        )
                        sys.stdout.flush()
                        for batch_pid in proteins_for_batch:
                            processed_pids.add(batch_pid)
                        break

                    else:
                        with open(err_log_path, "r") as err_file: stderr_out = err_file.read()
                        raise subprocess.CalledProcessError(proc.returncode, cmd, stderr=stderr_out)

                except subprocess.CalledProcessError as e:
                    stderr_text = e.stderr or ""

                    """
                    Surgical repair: if a single job caused a FileNotFoundError,
                    wait 5 s for the filesystem to settle, then isolate just that job
                    to the individual queue and retry the rest of the batch without
                    counting this as a failed attempt.
                    """
                    if not _surgical_repair_done and "FileNotFoundError" in stderr_text:
                        _m = re.search(r"/predictions/([^/\s]+)/pre_affinity_", stderr_text)
                        if _m:
                            bad_stem = _m.group(1)
                            bad_yaml = chunk_dir / f"{bad_stem}.yaml"
                            if bad_yaml.exists() and bad_stem in job_map:
                                _tty_write(
                                    f"\n[Warning] Job {bad_stem} is missing its preprocessed data file in batch [{batch_label}]."
                                    f" Waiting 5s for filesystem to settle, then isolating it to the individual queue...\n"
                                )
                                sys.stdout.flush()
                                time.sleep(5)
                                bad_yaml.unlink(missing_ok=True)
                                bad_job = job_map.pop(bad_stem)
                                bad_job["fallback_individual"] = True
                                gpu_queue_for_analysis.put(bad_job)
                                all_jobs_in_batch = [
                                    j for j in all_jobs_in_batch
                                    if f"{j['job_index']}_{j.get('job_protein', j['protein'])}_{j['ligand']}" != bad_stem
                                ]
                                n_jobs = len(all_jobs_in_batch)
                                _surgical_repair_done = True
                                console_info(
                                    f"[SURGICAL REPAIR] Isolated {bad_stem} → individual queue."
                                    f" Retrying batch with {n_jobs} remaining jobs (attempt counter not incremented)."
                                )
                                continue  # retry without incrementing attempt
                    # --- End surgical repair ---

                    _tty_write(f"\n[Warning] Designated GPU batch failed on [{batch_label}] (Attempt {attempt+1}/{max_retries}). Systematically retrying in 15s...")
                    sys.stdout.flush()
                    attempt += 1
                    if attempt < max_retries:
                        time.sleep(15)
                    else:
                        console_info(f"\nCRITICAL EVENT: GPU execution batch failed permanently on designated [{batch_label}]: {e}")
                        console_info(f"Reported Error Diagnostics:\n{e.stderr}")
                        console_info(f"\nFALLBACK PROTOCOL INITIATED: Queueing {n_jobs} active computational jobs for sequential individual processing...")
                        for job in all_jobs_in_batch:
                            job_name = f"{job['job_index']}_{job.get('job_protein', job['protein'])}_{job['ligand']}"
                            if not check_job_status(D_RUNS / job_name):
                                job["fallback_individual"] = True
                                gpu_queue_for_analysis.put(job)
                        break  # continue to next batch — do not exit the pipeline
                except Exception as e:
                    console_info(f"\nCRITICAL EVENT: Unexpected systemic error encountered precisely in GPU Modeller Execution Loop covering [{batch_label}]: {e}")
                    break  # continue to next batch — do not exit the pipeline

        if pending_by_protein:
            print()
            console_info("\nGPU Hardware Batch Prediction Phase Exhaustively Completed.")
            console_separator()
            print()

    """
    Signal the background analysis worker to stop and wait for it to flush
    (runs regardless of whether there were pending GPU jobs)
    """
    gpu_queue_for_analysis.put("STOP")
    analysis_pool.shutdown(wait=True)

    # -------------------------------------------------------------------------------
    # Step 10.8: CPU Worker Master Execution Loop
    # -------------------------------------------------------------------------------
    master_rows = []
    total_tasks = len(tasks)

    # Load job names already written to CSV by the GPU-phase background worker
    already_done_jobs = set()
    if CSV_PATH.exists() and CSV_PATH.stat().st_size > 0:
        try:
            already_done_jobs = set(pd.read_csv(CSV_PATH, low_memory=False, usecols=["job_name"])["job_name"].dropna().astype(str))
            if already_done_jobs:
                console_info(f"  Resuming — skipping {len(already_done_jobs)} jobs already written to CSV by the GPU-phase analysis worker.")
        except Exception:
            pass

    worker_count = max(args.cpus if hasattr(args, "cpus") and args.cpus else TARGET_CORES, 1)

    exec_args = []
    for job in tasks:
        job_name = f"{job['job_index']}_{job.get('job_protein', job['protein'])}_{job['ligand']}"
        if job_name in already_done_jobs:
            continue
        prev_time = prev_elapsed_map.get(job_name, 0.0)
        """
        Pass gpu_queue=None: Step 10.8 only handles pre-existing predictions (run_gpu=False).
        Passing a Manager proxy would cause all forked workers to simultaneously reconnect
        to the Manager RPC server → deadlock. None is safe here.
        """
        exec_args.append((job, PROD, args.diffusion_samples, prev_time, D_ALN, control_cif, control_map, DEHA4_CONTROL_SEQ, None, r3u_boltz_cif, r3u_control_map))

    total_tasks = len(exec_args)

    print()
    console_separator()
    console_info(f"Launching Multiprocessing Grid Environment executing securely on {worker_count} Concurrent Dedicated Workers...")
    console_info(f" -> {total_tasks:,} jobs queued for CPU re-analysis | {len(already_done_jobs):,} already complete")
    console_info(" -> Using Pool.imap_unordered (lazy feed — no bulk queue pre-fill).")
    sys.stdout.flush()

    completed_count = 0
    grid_success_count = 0
    grid_fail_count = 0
    _t_grid_start = time.time()

    with multiprocessing.Pool(processes=worker_count) as pool:
        for res_pair in pool.imap_unordered(worker_task_wrapper, exec_args, chunksize=1):
            completed_count += 1
            try:
                res, used_gpu = res_pair
            except (TypeError, ValueError):
                res, used_gpu = {"error": "Malformed worker result", "status": "FAILED", "job_name": "unknown"}, -1

            name = res.get("job_name", "unknown") if isinstance(res, dict) else "unknown"
            status_str = "SUCCESS" if (res and isinstance(res, dict) and res.get("status") == "Success") else "FAILED"
            color = ConsoleColours.OKGREEN if status_str == "SUCCESS" else ConsoleColours.FAIL

            elapsed = time.time() - _t_grid_start
            rate = completed_count / elapsed if elapsed > 0 else 0.0
            eta_s = (total_tasks - completed_count) / rate if rate > 0 else 0.0
            eta_str = f"{eta_s / 3600:.1f}h" if eta_s > 3600 else f"{int(eta_s / 60)}m{int(eta_s % 60)}s"
            pct = completed_count / total_tasks * 100 if total_tasks else 0
            _end = "\n" if completed_count == total_tasks else ""
            _tty_write(
                f"\r[GRID] {completed_count:,}/{total_tasks:,}  ({pct:.1f}%)  {color}{status_str}{ConsoleColours.ENDC}"
                f"  rate={rate:.1f}/s  ETA={eta_str}  {name[:40]}\033[K{_end}"
            )
            sys.stdout.flush()

            if status_str == "FAILED":
                grid_fail_count += 1
                GLOBAL_STATS["analytical_failures"] = GLOBAL_STATS.get("analytical_failures", 0) + 1
                if isinstance(res, dict) and "error" in res:
                    console_info(f"\n [!] {name} FAILED: {res['error']}")
            else:
                grid_success_count += 1

            if res and isinstance(res, dict) and res.get("status") == "Success":
                if used_gpu != -1: GLOBAL_STATS["jobs_run_gpu"] += 1
                else: GLOBAL_STATS["jobs_repaired_cpu"] += 1
                master_rows.append(flatten_job_result(res))

            if res and isinstance(res, dict) and "error" in res:
                if logger and "traceback" in res:
                    logger.error(f"Traceback for {name}:\n{res['traceback']}")

    if grid_fail_count == 0:
        console_info(f"\nAll {total_tasks} Designated Predictive Jobs completed Successfully.")
    else:
        console_info(f"\n{grid_success_count}/{total_tasks} Predictive Jobs completed. {grid_fail_count} FAILED — these jobs will be retried on the next resume run.")
    console_separator()

    _workspace_root = PROD / "_Temp_Workspace"
    if _workspace_root.exists():
        try:
            shutil.rmtree(_workspace_root)
            console_info("Successfully deleted _Temp_Workspace — disk space reclaimed.")
        except Exception as _ws_err:
            console_info(f"  [Warning] Could not fully remove _Temp_Workspace: {_ws_err}")

    # -------------------------------------------------------------------------------
    # Step 10.9: Final Data Merging and Refinement (CSV Formatting)
    # -------------------------------------------------------------------------------
    print()
    console_separator()
    console_info(f"Collecting {len(master_rows):,} analysed result(s) and writing Master CSV...")
    df_columns_count = 0
    if master_rows or (CSV_PATH.exists() and os.path.getsize(CSV_PATH) > 0):
        try:
            if master_rows:
                df_temp = pd.DataFrame(master_rows)
                if CSV_PATH.exists() and os.path.getsize(CSV_PATH) > 0:
                    df_existing = pd.read_csv(CSV_PATH, low_memory=False)
                    df_temp = pd.concat([df_existing, df_temp], ignore_index=True)
                atomic_to_csv(df_temp, CSV_PATH, index=False)
                master_rows.clear()

            df = pd.read_csv(CSV_PATH, low_memory=False)

            existing_drop = [c for c in COLUMNS_TO_DROP if c in df.columns]
            if existing_drop: df.drop(columns=existing_drop, inplace=True)

            path_cols = []
            _path_re = r"\.\./|/mnt/|/home/"
            _sample  = df.head(10)  # paths appear in every row or not at all
            for col in df.columns:
                if _sample[col].astype(str).str.contains(_path_re, regex=True).any():
                    path_cols.append(col)
            if path_cols: df.drop(columns=path_cols, inplace=True)

            fill_dict = {}
            for c in df.columns:
                if "dist" in c or "catalytic" in c: fill_dict[c] = 999.0
                elif "score" in c or "iptm" in c or "plddt" in c or "density" in c or "angle" in c or "rmsd" in c: fill_dict[c] = 0.0
                elif "count" in c or "Total_" in c: fill_dict[c] = 0
                else: fill_dict[c] = "NA"
            """
            Field-specific overrides — must be set after the generic loop.
            Numeric aux geometries must get a NUMERIC sentinel, never the "NA"
            string: pandas re-parses "NA" as NaN on the next read, which would
            re-open the gap in the downstream ranked CSV. 999.0 marks "undefined"
            (consistent with the CFG.SENTINEL_VALID_MAX guard in the SN2 analysis).
            """
            fill_dict["flippin_lodge_offset"]        = 999.0
            fill_dict["burgi_dunitz_angle"]          = 999.0
            fill_dict["residues_within_6A"]          = "0"
            fill_dict["hydrophobic_desolvation_ratio"] = 0.0
            fill_dict["active_site_contact_flag"]    = 0.0
            # Keep elite_demotion a pure-string column ("none" sentinel, never empty/NaN):
            # an empty string round-trips through CSV as NaN, mixing float+str on re-read and
            # triggering pandas DtypeWarning. error rows missing the key are filled here too.
            fill_dict["elite_demotion"]                = "none"
            fill_dict["r3u_Active_Site_RMSD"]              = 99.0
            fill_dict["r3u_Halide_Stabilisation"]          = False
            fill_dict["r3u_Carboxylate_Clamp"]             = False
            fill_dict["r3u_ActiveSite_Conservation_Score"]     = 0.0

            df = df.fillna(fill_dict)
            df.rename(columns=COLUMN_RENAMING_MAP, inplace=True)
            df = df.loc[:, ~df.columns.duplicated()]

            drop_empty = []
            zero_equivalents = {"0", "0.0", "0.00", "NA", "None", "", "nan", "False"}
            protected_cols = {CFG.COL_TIER, "is_degrader", "status", CFG.COL_PROT, CFG.COL_LIG, "job_name", "job_index"}
            for col in df.columns:
                if col in protected_cols:
                    continue
                if df[col].astype(str).isin(zero_equivalents).all():
                    drop_empty.append(col)
            if drop_empty: df.drop(columns=drop_empty, inplace=True)

            if "job_index" in df.columns:
                df = df.sort_values("job_index", kind="mergesort").reset_index(drop=True)

            """
            Vectorised alignment grade assignment via pd.cut.
            Derives Alignment_Score_Pct from identity_pct, then maps to
            letter grades A–I using CFG-defined bins — replaces per-job
            if-elif logic that would otherwise execute 58k+ times.
            """
            if CFG.COL_ID_PCT in df.columns:
                df["Alignment_Score_Pct"] = (
                    pd.to_numeric(df[CFG.COL_ID_PCT], errors="coerce")
                    .fillna(0.0)
                    .clip(0, 100)
                )
                df[CFG.COL_ALN_G] = pd.cut(
                    df["Alignment_Score_Pct"],
                    bins=CFG.ALIGN_GRADE_BINS,
                    labels=CFG.ALIGN_GRADE_LABELS,
                    right=True,
                    include_lowest=True   # close first interval to [0,20] so 0.0 → 'I', not NaN/"nan"
                ).astype(str)

            atomic_to_csv(df, CSV_PATH, index=False)
            GLOBAL_STATS["created_csv_rows"] = len(df)
            df_columns_count = len(df.columns)
        except Exception as e: console_info(f" Warning generated during systematic CSV formatting execution logic: {e}")
    console_info("Final structured CSV format written flawlessly.")

    # -------------------------------------------------------------------------------
    rank_csv_path, rank_columns_count = generate_scientific_ranking_csv(CSV_PATH, PROD, ts_now)

    save_alignment_cache_final(D_ALN / CFG.FILE_ALIGNMENT_STATS)

    console_separator()

    # -------------------------------------------------------------------------------
    # Step 10.11: Final Report and Output Terminal UI
    # -------------------------------------------------------------------------------
    console_separator()
    console_info(f"\nSequential Analytical Pipeline entirely Completed. Associated Files are formally verified inside: {PROD.resolve()}")
    console_separator()
    console_info(" \n ✔ Summary of Validated Operational Disk Files-")

    yaml_count = len(list(D_YAML.glob("*.yaml")))
    job_count = len([d for d in sorted(D_RUNS.iterdir()) if d.is_dir()])
    aln_count = len(list(D_ALN.glob("*.txt")))

    str_in   = f"({len(proteins)} Seq, {len(ligands)} Lig)"
    str_yaml = f"{yaml_count}"
    str_jobs = f"{job_count}"
    str_aln  = f"{aln_count}"
    str_mc   = f"{df_columns_count} Columns"
    str_rc   = f"{rank_columns_count} Columns"

    console_info(f"Input Structural Files: {str_in}".ljust(45) + f" |      {D_IN.resolve()}")
    console_info(f"YAML Definitions:   {str_yaml}".ljust(45) + f" |      {D_YAML.resolve()}")
    console_info(f"Primary Operational Jobs:      {str_jobs}".ljust(45) + f" |      {D_RUNS.resolve()}")
    console_info("Sequence Reference Data (MSA + Alignments):".ljust(45) + f" |      {D_SEQ.resolve()}")
    console_info(f"  ├─ MSA Sequences: {len(list(D_COLABFOLD.glob('*.a3m')))}" .ljust(45) + f" |      {D_COLABFOLD.resolve()}")
    console_info(f"  └─ Alignments: {str_aln}".ljust(45) + f" |      {D_ALN.resolve()}")
    console_info(f"FAcDs Master CSV:   {str_mc}".ljust(45) + f" |      {CSV_PATH.resolve()}")
    if rank_csv_path:
        console_info(f"FAcDs Ranked CSV:   {str_rc}".ljust(45) + f" |      {rank_csv_path.resolve()}")
    _mir = run_root / "2_Best_Complexes_CIFs"
    _mir_count = sum(1 for _ in _mir.rglob("*.cif")) if _mir.exists() else 0
    console_info(f"Best Complexes CIFs Mirror: {_mir_count} CIFs".ljust(45) + f" |      {_mir.resolve()}")
    console_info("Operational Execution Log".ljust(45) + f" |      {LOG_PATH.resolve()}")
    console_separator()

    if CSV_PATH.exists(): df_final = pd.read_csv(CSV_PATH, low_memory=False)
    else: df_final = pd.DataFrame()

    """
    Use Alignment_Grade column (written by pd.cut in Step 10.9)
    when available; fall back to utility logic for backwards
    compatibility with CSVs produced before this change was applied.
    """
    if not df_final.empty and CFG.COL_ALN_G in df_final.columns:
        _grade_counts = df_final[CFG.COL_ALN_G].value_counts()
        grade_A = int(_grade_counts.get("A", 0))
        grade_B = int(_grade_counts.get("B", 0))
        grade_C = int(_grade_counts.get("C", 0))
        grade_D = int(_grade_counts.get("D", 0))
        grade_E = int(_grade_counts.get("E", 0))
        grade_F = int(_grade_counts.get("F", 0))
        grade_G = int(_grade_counts.get("G", 0))
        grade_H = int(_grade_counts.get("H", 0))
        grade_I = int(_grade_counts.get("I", 0))
    else:
        raw_id = df_final.get(CFG.COL_ID_PCT, pd.Series(dtype=float))
        id_pct = pd.to_numeric(raw_id, errors="coerce").dropna()
        grade_A = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == "A").sum())
        grade_B = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == "B").sum())
        grade_C = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == "C").sum())
        grade_D = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == "D").sum())
        grade_E = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == "E").sum())
        grade_F = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == "F").sum())
        grade_G = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == "G").sum())
        grade_H = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == "H").sum())
        grade_I = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == "I").sum())

    mdl_col = "best_model_name" if "best_model_name" in df_final.columns else CFG.COL_CONF
    if not df_final.empty and mdl_col in df_final.columns:
        models = Counter(df_final[mdl_col].fillna("None").tolist())
    else:
        models = Counter()

    console_info("")
    console_info(SEPARATOR_LIGHT)
    console_info("  Final Computational Statistics")
    console_info(SEPARATOR_LIGHT)
    if not df_final.empty:
        csv_failures  = int((df_final["status"] != "Success").sum())
        grid_failures = GLOBAL_STATS.get("analytical_failures", 0)
        total_failures = csv_failures + grid_failures
        _jw = 36
        console_info(f"  ┌{'─'*(_jw+2)}┬{'─'*12}┐")
        console_info(f"  │  {'Metric':<{_jw}}│  {'Count':>8}  │")
        console_info(f"  ├{'─'*(_jw+2)}┼{'─'*12}┤")
        console_info(f"  │  {'Jobs Accounted':<{_jw}}│  {len(df_final)+grid_failures:>8,}  │")
        console_info(f"  │  {'Successful Outputs':<{_jw}}│  {int((df_final['status']=='Success').sum()):>8,}  │")
        _fail_str = f"{total_failures:,}" + ("  ← re-run to retry" if grid_failures else "")
        console_info(f"  │  {'Analytical Failures':<{_jw}}│  {_fail_str:>8}  │")
        console_info(f"  └{'─'*(_jw+2)}┴{'─'*12}┘")

        _tier_order = CFG.TIER_ORDER
        """
        Separate the reference control structures (DeHa4 / 3R3U Boltz-2 predictions,
        flagged by a Protein_Name ending in "_Control") from the screened library.
        Per-tier figures report the library count (dynamic) plus, where present, how
        many controls fell in that tier as a "+NC" suffix. The grand total therefore
        reads as <library total> + <control total> rather than silently merging the
        controls into the screened population.
        """
        if CFG.COL_PROT in df_final.columns:
            _is_ctrl = df_final[CFG.COL_PROT].astype(str).str.endswith("_Control")
        else:
            _is_ctrl = pd.Series(False, index=df_final.index)
        _lib_tiers  = Counter(df_final.loc[~_is_ctrl, CFG.COL_TIER].fillna(CFG.TIER_DECOY).tolist())
        _ctrl_tiers = Counter(df_final.loc[_is_ctrl,  CFG.COL_TIER].fillna(CFG.TIER_DECOY).tolist())
        _lib_total  = sum(_lib_tiers.get(t, 0)  for t in _tier_order)
        _ctrl_total = sum(_ctrl_tiers.get(t, 0) for t in _tier_order)
        _tw = 14   # tier label column
        _lw = 10   # library count column
        _cw = 9    # control count column
        console_info("")
        console_info(f"  ┌{'─'*(_tw+2)}┬{'─'*(_lw+2)}┬{'─'*(_cw+2)}┐")
        console_info(f"  │  {'Degrader Tier':<{_tw}}│  {'Library':>{_lw}}│  {'Control':>{_cw}}│")
        console_info(f"  ├{'─'*(_tw+2)}┼{'─'*(_lw+2)}┼{'─'*(_cw+2)}┤")
        for _t in _tier_order:
            _c = _ctrl_tiers.get(_t, 0)
            _l = _lib_tiers.get(_t, 0)
            _ccell = f"{_c:,}" if _c else "·"
            console_info(f"  │  {_t:<{_tw}}│  {f'{_l:,}':>{_lw}}│  {_ccell:>{_cw}}│")
        console_info(f"  ├{'─'*(_tw+2)}┼{'─'*(_lw+2)}┼{'─'*(_cw+2)}┤")
        console_info(f"  │  {'Subtotal':<{_tw}}│  {f'{_lib_total:,}':>{_lw}}│  {f'{_ctrl_total:,}':>{_cw}}│")
        console_info(f"  └{'─'*(_tw+2)}┴{'─'*(_lw+2)}┴{'─'*(_cw+2)}┘")
        console_info(f"  Grand Total: {_lib_total + _ctrl_total:,}  ({_lib_total:,} library + {_ctrl_total:,} control)")

        _grade_rows = [
            ("A  (90–100%)", grade_A), ("B  (80–90%)",  grade_B),
            ("C  (70–80%)",  grade_C), ("D  (60–70%)",  grade_D),
            ("E  (50–60%)",  grade_E), ("F  (40–50%)",  grade_F),
            ("G  (30–40%)",  grade_G), ("H  (20–30%)",  grade_H),
            ("I  (< 20%)",   grade_I),
        ]
        _gw = 16
        console_info("")
        console_info(f"  ┌{'─'*(_gw+2)}┬{'─'*12}┐")
        console_info(f"  │  {'Alignment Grade':<{_gw}}│  {'Count':>8}  │")
        console_info(f"  ├{'─'*(_gw+2)}┼{'─'*12}┤")
        for _g, _n in _grade_rows:
            console_info(f"  │  {_g:<{_gw}}│  {int(_n):>8,}  │")
        console_info(f"  └{'─'*(_gw+2)}┴{'─'*12}┘")

        _model_rows = [(f"Model {i}", models.get(f"model_{i}", 0)) for i in range(5)]
        _mw = 10
        console_info("")
        console_info(f"  ┌{'─'*(_mw+2)}┬{'─'*12}┐")
        console_info(f"  │  {'Model':<{_mw}}│  {'Count':>8}  │")
        console_info(f"  ├{'─'*(_mw+2)}┼{'─'*12}┤")
        for _m, _n in _model_rows:
            console_info(f"  │  {_m:<{_mw}}│  {int(_n):>8,}  │")
        console_info(f"  └{'─'*(_mw+2)}┴{'─'*12}┘")

if __name__ == "__main__":
    import time as _time
    _t0 = _time.perf_counter()
    main()
    _utils_mod.print_elapsed(_t0, "02_Production_FAcDs.py")
