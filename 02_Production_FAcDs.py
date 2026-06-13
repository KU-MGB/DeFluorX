#!/usr/bin/env python3

"""
===============================================================================
FAcDs Pipeline  |  Step 02  |  Boltz-2 Production, Analysis & FAcD Ranking
===============================================================================
Large-scale, resume-safe Boltz-2 protein-ligand predictions with deep
structural, geometric, and chemical scoring for FAcD SN2 degrader tiers.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 10 June 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 02_Production_FAcDs.py
  Role          : "Engine" — Boltz-2 prediction orchestrator and ranker.
  Imports from  : 00_02_Project_Config_FAcDs.py  (CFG — all geometric thresholds)
                  00_03_Project_Utils_FAcDs.py   (geometric utilities, console funcs)
  Reads         : 01_Merge_FAcDs.py output — merged *.fasta (protein sequences)
                  User-supplied *.smi (SMILES ligand file)
  Writes        : <Run>/1_Boltz2_Production/  (Boltz-2 CIF outputs)
                  <Run>/2_Best_Complexes_CIFs/ (top-model CIF selection)
                  <Run>/1_Boltz2_Production/7_Boltz2_FAcDs_Ranked_*.csv
                  <Run>/1_Boltz2_Production/6_Boltz2_FAcDs_Master_*.csv
  Upstream      : 01_Merge_FAcDs.py → writes the merged FASTA consumed here
  Downstream    : 03_Validation_Figures_FAcDs.py → reads ranked CSV
                  05_CIF-PDB_Preparation_FAcDs.py → reads Best_Complexes_CIFs
                  06_Top-N_Extraction_FAcDs.py → reads ranked CSV
───────────────────────────────────────────────────────────────────────────────

# ── The Critic's Corner: Known Limitations & Failure Points ──────────────────
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
       • Mechanistic Fingerprint >= 0.9 (Near-complete anchor set required — GATE)
       • (Active Site Conservation is reported downstream for ranking; it is NOT a tier gate.)
       • Nucleophile (Asp110): <= 3.0 A (Tight pre-reactive ground-state gate)
       • Base (His277):        <= 3.5 A (Optimal Proton Transfer Gate)
       • Acid (Asp134):        <= 4.5 A (Strict Internal Gate)
       • Attack Angle:         >= 175°  (Near-Ideal Linear Trajectory)
       • Stabilisation:        REQUIRED (Trp156/Tyr217 or Dynamic Polar Residue)

    2. Tier_1B (High Functional)
       • Mechanistic Fingerprint >= 0.7   
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

    2. Mechanistic Fingerprint (Secondary Sort Key):
       Within the same Tier, candidates are ranked by their physical integrity score.

    3. Likelihood Score (Tertiary Sort Key):
       Final tie-breaker based on the sequence/structure hybrid score.

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
    6b. FAcD SN2 Defluorination QM/MM Energetics (basis for QSite scan design):
       - Comprehensive understanding of FAcD-catalyzed degradation of fluorocarboxylic acids: a QM/MM approach.
       - Yue, Y. et al. (2021) Environ Sci Technol 55(14):9817–9825.
       - DOI: https://doi.org/10.1021/acs.est.0c08811
    7. Terminal CF3 Selectivity Against SN2 Attack:
       - Selective defluorination of trifluoromethyl substituents by conformationally induced remote substitution.
       - Jesani, M.H., Schwarz, M., Kim, S. et al. & Clayden, J. (2024) Angew Chem Int Ed 63:e202403477.
       - DOI: https://doi.org/10.1002/anie.202403477
    8. Directed Evolution of FAcD on Non-Natural Organofluorides:
       - Engineering fluoroacetate dehalogenase by growth-based selections on non-natural organofluorides.
       - Jansen, S.C., van Beers, P. & Mayer, C. (2026) Angew Chem Int Ed 65:e202524234.
       - DOI: https://doi.org/10.1002/anie.202524234

    9. Computational Tooling & Core Bioinformatics Stack:
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

===============================================================================
"""

# ===============================================================================
# SECTION 1: SYSTEM CONFIGURATION & IMPORTS
# ===============================================================================
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
import logging
import csv
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

from typing import Any, Dict, List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

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
    
    # Reduces the process scheduling priority so the pipeline yields CPU time to the operating system 
    # and interactive applications, preventing desktop freezes under sustained computational loads.
    try:
        p = psutil.Process(os.getpid())
        if sys.platform != 'win32':
            p.nice(15)  # POSIX niceness 15 designates a below-normal priority
    except Exception:
        pass

    # Enforces single-threaded BLAS/MKL execution per worker process.
    # When running numerous parallel workers, unconstrained BLAS thread pools cause severe
    # CPU over-subscription (where context-switching overhead dominates computational time).
    # Pinning each worker to one thread yields near-linear scaling with the worker count.
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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# -------------------------------------------------------------------------------
# Step 1.4: Bioinformatics & Chemistry
# -------------------------------------------------------------------------------
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
import gemmi                                           
from rdkit import Chem                                  
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

# -------------------------------------------------------------------------------
# Sub-Step 1.4.1: Torch Threading Configuration
# -------------------------------------------------------------------------------
try:
    import torch
    if torch is not None:
        # Restricts PyTorch intra-op parallelism to one thread per worker process,
        # ensuring consistency with the BLAS/MKL thread limits established above.
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
_cfg_mod   = _load_module("ProjectConfig", _REPO_DIR / "00_02_Project_Config_FAcDs.py")
_utils_mod = _load_module("ProjectUtils",  _REPO_DIR / "00_03_Project_Utils_FAcDs.py")
CFG        = _cfg_mod.CFG()

# Import centralised Ramachandran plotting helpers
compute_ramachandran_angles = _utils_mod.compute_ramachandran_angles
_rama_stats = _utils_mod._rama_stats
save_ramachandran_plot = _utils_mod.save_ramachandran_plot
save_ramachandran_comparison = _utils_mod.save_ramachandran_comparison
safe_name = _utils_mod.safe_name
# Auxiliary (non-gating) trajectory-geometry descriptors — see utils docstrings.
calculate_burgi_dunitz = _utils_mod.calculate_burgi_dunitz
calculate_flippin_lodge = _utils_mod.calculate_flippin_lodge
distance = _utils_mod.distance
calculate_angle = _utils_mod.calculate_angle


# ===============================================================================
# SECTION 2: PHYSICAL CONSTANTS & PARAMETERS (MAESTRO STANDARDS)
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 2.1: Path & Model Configuration
# -------------------------------------------------------------------------------
BOLTZ_BIN      = CFG.BOLTZ_EXECUTABLE
BOLTZ_MODEL    = CFG.BOLTZ_MODEL_VERSION
BOLTZ_CACHE    = Path.home() / ".cache" / "boltz"
OUTPUT_FORMAT  = CFG.BOLTZ_OUTPUT_FORMAT

# -------------------------------------------------------------------------------
# Step 2.2: Execution Settings
# -------------------------------------------------------------------------------
RECYCLING_STEPS           = CFG.BOLTZ_RECYCLING_STEPS
DIFFUSION_SAMPLES_DEFAULT = CFG.BOLTZ_DIFFUSION_SAMPLES
RETRY_ON_FAIL             = CFG.BOLTZ_RETRY_MAX
RETRY_SLEEP               = CFG.BOLTZ_RETRY_SLEEP
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
# Standard chemical bond definitions derived from the Schrödinger Maestro software suite.
# -------------------------------------------------------------------------------
# Sub-Step 2.4.1: Bond Distances
# -------------------------------------------------------------------------------
HBOND_MAX_DIST_HA        = CFG.THRESHOLD_HB_DIST_HA
HBOND_MAX_DIST_DA        = CFG.THRESHOLD_HB_DIST_MAX
HBOND_MIN_ANGLE          = CFG.THRESHOLD_HB_ANGLE_MIN
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

# Sets of amino acids grouped by their chemical properties for analytical purposes.
STANDARD_AA  = CFG.PREP_STANDARD_AA
POSITIVE_RES = {"ARG", "LYS", "HIS"}
NEGATIVE_RES = {"ASP", "GLU"}
METALS       = CFG.METALS

RES_PROPS = {
    "NonPolar": {"ALA", "VAL", "LEU", "ILE", "MET", "PRO", "GLY", "PHE", "TRP"},
    "Aromatic": {"PHE", "TYR", "TRP", "HIS"},
    "Polar":    {"SER", "THR", "ASN", "GLN", "CYS", "TYR", "HIS"},
    "Pos":      {"LYS", "ARG", "HIS"},
    "Neg":      {"ASP", "GLU"}
}
AROMATIC_RING_ATOMS = ["CG","CD1","CD2","CE1","CE2","CZ"]

# Functional sidechain atoms used to strictly validate dynamic stabilisation.
# This prevents backbone atoms from being incorrectly flagged as stabilising.
POLAR_SIDECHAIN_ATOMS = {
    "SER": {"OG"}, "THR": {"OG1"}, "ASN": {"OD1", "ND2"}, "GLN": {"OE1", "NE2"},
    "CYS": {"SG"}, "TYR": {"OH"}, "HIS": {"ND1", "NE2"}, 
    "TRP": {"NE1"}, "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"},
    "ARG": {"NH1", "NH2", "NE"}, "LYS": {"NZ"}
}

# Residue groups designed for Mechanistic Fingerprinting evaluations.
FINGERPRINT_GROUPS = {
    "AROMATIC": {"TRP", "TYR", "PHE", "HIS"},
    "POSITIVE": {"ARG", "LYS", "HIS"},
    "NEGATIVE": {"ASP", "GLU"},
    "POLAR":    {"SER", "THR", "ASN", "GLN", "TYR"}
}

# -------------------------------------------------------------------------------
# Step 2.6: The Geometric Reference Data (PDB 3R3U)
# -------------------------------------------------------------------------------
# The canonical sequence for DeHa4_[Delftia_acidovorans_D4B].
# This serves as a structural map. As target proteins vary in their numbering, 
# each is aligned to this canonical sequence to precisely map the positions of 
# key catalytic residues (such as Asp110, Asp134, His277).
# Source: https://doi.org/10.1021/acsomega.4c02517 (Position cross-validated, establishing HIS277 rather than HIS288 in DEHA4 sequences).
DEHA4_CONTROL_SEQ = CFG.DEHA4_CONTROL_SEQ
REF_SEQUENCE_STR = DEHA4_CONTROL_SEQ

# Fluoroacetate SMILES used for reference construction operations.
REF_LIGAND_SMILES = CFG.FLUOROACETATE_SMILES
DEHA4_CONTROL_SMI = REF_LIGAND_SMILES

# The three fluoroacetate control substrates used for DEHA4 and 3R3U calibration runs.
CTRL_LIGANDS = CFG.CTRL_LIGANDS

# Structural mapping based on PDB 3R3U (Rhodopseudomonas palustris FAcD,
# wild-type, 1.60 Å; carries Ni²⁺/Cl⁻ ions, no substrate bound),
# adapted from https://doi.org/10.1021/acsomega.4c02517 and matched against DEHA4 sequences.
REF_ACTIVE_SITE_MAP  = CFG.REF_ACTIVE_SITE_MAP
CATALYTIC_TRIAD_KEYS = CFG.CATALYTIC_TRIAD_KEYS

# -------------------------------------------------------------------------------
# Step 2.7: Output Formatting Configuration
# -------------------------------------------------------------------------------
# Mapping dictionary to convert internal variable names into user-friendly CSV column headers.
COLUMN_RENAMING_MAP = {
    "protein": "Protein_Name",
    "ligand": "Ligand_Name",
    "dist_Nuc": "Dist_Nucleophile_ASP110",
    "dist_Base": "Dist_Base_HIS277",
    "dist_Acid": "Dist_Acid_ASP134",
    "dist_Stab_W": "Dist_Stabiliser_TRP156",
    "dist_Stab_Y": "Dist_Stabiliser_TYR217",
    "dist_Carb1": "Dist_Clamp_ARG111",
    "dist_Carb2": "Dist_Clamp_ARG114",
    "binding_likelihood_computed": "Binding_Probability_Score",
    "confidence_score": "Boltz_Model_Confidence",
    "sn2_attack_angle": "SN2_Attack_Angle",
    "sn2_trajectory_dev": "SN2_Trajectory_Deviation_A",
    "Mapped_to_Control_All": "Full_Sequence_Alignment_Map",
    "Mapped_to_Control_Cat_Triad": "Active_Site_Triad_Map",
    "ActiveSite_Conservation_Score": "ActiveSite_Conservation_Score",
    "Mechanistic_Fingerprint_Score": "Mechanistic_Fingerprint_Score",
    "Active_Site_RMSD": "Active_Site_RMSD_to_Control",
    "Halide_Stabilisation": "Has_Halide_Stabilisation",
    "Carboxylate_Clamp": "Has_Carboxylate_Clamp"
}

# A list of internal columns to exclude from the final CSV output to maintain clarity.
COLUMNS_TO_DROP = [
    "protein_id", "catalytic_dist_A", "rc", "cross_interface_pae_mean", "completed_at",
    "active_site_mapping", "scientific_meaning", "constraint_check"
]

# -------------------------------------------------------------------------------
# Step 2.8: Default Metrics Template
# -------------------------------------------------------------------------------
DEFAULT_METRICS = {
    "status": "Unknown",
    "elapsed_seconds": 0.0,
    "degrader_tier": CFG.TIER_DECOY,
    "is_degrader": False,
    "ActiveSite_Conservation_Score": 0.0,
    "Mechanistic_Fingerprint_Score": 0.0,
    "Active_Site_RMSD": 999.0,
    "Identity_to_Control": 0.0,
    "Halide_Stabilisation": False,
    "Carboxylate_Clamp": False,
    "catalytic_dist_A": 999.0,
    "interaction_density": 0.0,
    "custom_affinity_score": 0.0,
    "binding_likelihood_computed": 0.0,
    "confidence_score": 0.0,
    "sn2_attack_angle": 0.0,
    "sn2_trajectory_dev": 999.0,
    "burgi_dunitz_angle": float("nan"),     # aux reference (carbonyl); non-gating
    "flippin_lodge_offset": float("nan"),   # aux reference (in-plane offset); non-gating
}
# Initialises all specific residue distances to an arbitrary maximum (999.0).
for k in REF_ACTIVE_SITE_MAP:
    DEFAULT_METRICS[f"dist_{k}"] = 999.0

# -------------------------------------------------------------------------------
# Step 2.9: Global State Variables
# -------------------------------------------------------------------------------

# ConsoleColours sourced from 00_03_Project_Utils (single canonical definition).
ConsoleColours  = _utils_mod.ConsoleColours
SEPARATOR_HEAVY = _utils_mod.SEPARATOR_HEAVY
SEPARATOR_LIGHT = _utils_mod.SEPARATOR_LIGHT

logger = None
_gpu_lock = threading.Lock() 
_gpu_free_mib: Dict[int, int] = {}
CACHED_ALIGNMENTS: Dict[str, Dict] = {} 

GLOBAL_STATS = {
    "renamed_yaml": 0, "renamed_dirs": 0, "renamed_files": 0,
    "created_yaml": 0, "created_aln": 0, "created_json": 0, "created_csv_rows": 0,
    "csv_columns": 0, "jobs_run_gpu": 0, "jobs_repaired_cpu": 0,
    "total_alignments_generated": 0,
    "complexes_copied": 0,
    "analytical_failures": 0,
}


# ===============================================================================
# SECTION 3: UTILITY FUNCTIONS (LOGGING & DISPLAY)
# ===============================================================================

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

_TTY_ANSI_RE = re.compile(r'\033\[[0-9;]*[mKABCDEFGHJSTfhilmnprsu]')

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
    
    for item in runs_dir.glob(f"*_{protein_id}_*"):
        if item.is_dir():
            target = trash_dir / f"{item.name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            try:
                shutil.move(str(item), str(target))
            except Exception:
                shutil.rmtree(item, ignore_errors=True)


def job_key_from_yaml(yaml_path: Path) -> Optional[Tuple[str, str]]:
    """Extracts the stable resume key from an existing YAML definition file."""
    try:
        with open(yaml_path, 'r') as yf:
            data = yaml.safe_load(yf)
        seq = data.get('sequences', [])[0].get('protein', {}).get('sequence', '')
        smi_raw = data.get('sequences', [])[1].get('ligand', {}).get('smiles', '')
        smi = smi_raw
        if smi_raw:
            mol = Chem.MolFromSmiles(smi_raw)
            if mol is not None:
                smi = Chem.MolToSmiles(mol)
        if seq and smi:
            return sequence_hash(seq), smi
    except Exception:
        pass
    return None


def reconcile_resume_job_names(prod_dir: Path, expected_job_stems: Dict[Tuple[str, str], str]):
    """Renames existing YAMLs and run folders to match the current input ordering during resume."""
    yaml_dir = prod_dir / "2_Boltz2_YAML_Configs"
    runs_dir = prod_dir / "4_Prediction_Jobs"
    if not yaml_dir.exists():
        return

    for yaml_path in sorted(yaml_dir.glob("*.yaml"), key=lambda p: p.name):
        key = job_key_from_yaml(yaml_path)
        if not key or key not in expected_job_stems:
            continue

        expected_stem = expected_job_stems[key]
        current_stem = yaml_path.stem
        if current_stem == expected_stem:
            continue

        target_yaml = yaml_dir / f"{expected_stem}.yaml"
        if target_yaml.exists():
            try:
                if yaml_path.read_text(encoding='utf-8') == target_yaml.read_text(encoding='utf-8'):
                    yaml_path.unlink(missing_ok=True)
                else:
                    alternate = yaml_dir / f"{expected_stem}_{int(time.time())}.yaml"
                    yaml_path.rename(alternate)
                    logger.debug(f"Preserved conflicting YAML {yaml_path.name} as {alternate.name}")
            except Exception:
                pass
        else:
            try:
                yaml_path.rename(target_yaml)
                GLOBAL_STATS["renamed_yaml"] += 1
            except Exception:
                pass

        old_run_dir = runs_dir / current_stem
        new_run_dir = runs_dir / expected_stem
        if old_run_dir.exists() and old_run_dir != new_run_dir:
            try:
                deep_rename_job_folder(old_run_dir, current_stem, expected_stem)
                if new_run_dir.exists():
                    shutil.rmtree(new_run_dir, ignore_errors=True)
                old_run_dir.rename(new_run_dir)
                GLOBAL_STATS["renamed_dirs"] += 1
            except Exception:
                pass


# ===============================================================================
# SECTION 4: REFERENCE & SELF-HEALING SYSTEM (DATA RECOVERY)
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 4.1: Geometric Truth Initialisation
# -------------------------------------------------------------------------------
def setup_reference_data(prod_dir: Path):
    """
    Initialises the reference directory. Checks for the presence of the reference PDB (3R3U)
    and associated sequence/SMILES files, subsequently downloading the PDB if absent.
    """
    ref_dir = prod_dir / CFG.REFERENCE_DIR_NAME
    ref_dir.mkdir(parents=True, exist_ok=True)

    pdb_path   = ref_dir / f"{CFG.REFERENCE_PDB_ID}.pdb"
    fasta_path = ref_dir / "DeHa4_Ref.fasta"
    smi_path   = ref_dir / "Control_Ligands_TFA_FA_DFA_Ref.smi"

    console_info(f"\n{SEPARATOR_LIGHT}")
    console_info(f"  Reference Data Initialisation")
    console_info(SEPARATOR_LIGHT)

    # 1. Download PDB 3R3U if missing
    if not pdb_path.exists():
        url = CFG.REFERENCE_PDB_URL
        try:
            console_info(f"  Downloading PDB 3R3U from {url}...")
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

    # 3. Create Reference SMILES (all 3 control ligands)
    # Legacy single-ligand file is superseded — always regenerate if stale or absent.
    _legacy_smi = ref_dir / "Fluoroacetate_Ref.smi"
    if _legacy_smi.exists():
        try: _legacy_smi.unlink()
        except Exception: pass
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
                        seq += code if code not in ('?', ' ') else 'X'
                if len(seq) >= 50:
                    return seq
    except Exception as e:
        if logger:
            logger.debug(f"Sequence extraction from {pdb_path.name} failed: {e}")
    return None


# -------------------------------------------------------------------------------
# Step 4.2: File Renaming & MSA Algorithms
# -------------------------------------------------------------------------------
def deep_rename_job_folder(job_dir: Path, old_stem: str, new_stem: str):
    """
    Recursively updates filenames within a job folder to reflect changes in identifier formats or protein names.
    Matches against the complete string stem to ensure thorough renaming across all file types (such as summary.json, .cif).
    """
    for root, dirs, files in os.walk(job_dir, topdown=False):
        for name in files:
            if old_stem in name:
                new_name = name.replace(old_stem, new_stem, 1)
                try: (Path(root)/name).rename(Path(root)/new_name)
                except Exception as e: logger.debug(f"Renaming operation deferred: {e}")
        for name in dirs:
            if old_stem in name:
                new_name = name.replace(old_stem, new_stem, 1)
                try: (Path(root)/name).rename(Path(root)/new_name)
                except Exception as e: logger.debug(f"Directory renaming operation deferred: {e}")

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
        with open(a3m_path, 'rb') as f:
            header = f.read(1024)
            # Checks for standard FASTA header initiation
            if not header.startswith(b'>'):
                return False
            # Scans header strings for anomalous null bytes
            if b'\x00' in header:
                return False
            # Iteratively scans the remainder of the file in chunks 
            while True:
                chunk = f.read(65536)
                if not chunk:
                    break
                if b'\x00' in chunk:
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
        # The ColabFold tarball contains uniref.a3m and bfd.mgnify30.metaeuk30.smag30.a3m files.
        # Merging both files yields maximum MSA coverage parameters.
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

        # Strips out problematic null bytes generated by underlying structural tools.
        a3m_content = b"\n".join(a3m_parts).replace(b"\x00", b"")

        # --- Step 5: Execute atomic write and metadata creation ---
        tmp_path = a3m_path.with_suffix(".tmp")
        tmp_path.write_bytes(a3m_content)

        if not validate_a3m_file(tmp_path):
            logger.debug(f"[MSA-direct] Downloaded A3M sequence failed post-processing validation for {pid}")
            tmp_path.unlink(missing_ok=True)
            return False

        tmp_path.rename(a3m_path)
        with open(meta_path, 'w') as mf:
            json.dump({
                'sequence_sha256': sequence_hash(seq),
                'protein_id':      pid,
                'generated_at':    datetime.utcnow().isoformat() + 'Z',
                'source':          'colabfold_api_direct',
                'a3m_bytes':       len(a3m_content),
            }, mf)

        logger.debug(f"[MSA-direct] MSA data saved for {pid} ({len(a3m_content)} bytes, {len(a3m_parts)} constituent files)")
        return True

    except Exception as e:
        logger.debug(f"[MSA-direct] Processing exception for {pid}: {type(e).__name__}: {e}")
        return False


# -------------------------------------------------------------------------------
# Step 4.4: Stale YAML and Run Folder Cleanup
# -------------------------------------------------------------------------------
def sanitize_dataset(prod_dir: Path):
    """Removes duplicate or incorrectly numbered YAML configuration files to prevent phantom tasks during resumed runs."""
    console_info("-------------------------------------------------------------------------------")
    console_info("YAML Duplicate Detection and Maintenance Cleanup")
    
    yaml_dir = prod_dir / "2_Boltz2_YAML_Configs"
    if yaml_dir.exists():
        all_yamls = sorted(list(yaml_dir.glob("*.yaml")))
        seen_stems = {} 
        
        for f in all_yamls:
            match = re.match(r"^(\d+)_(.+)$", f.name)
            if match:
                curr_id, rest = match.group(1), match.group(2)
                if rest in seen_stems:
                    old_f, old_id = seen_stems[rest]
                    if len(curr_id) > len(old_id):
                        if old_f.exists(): old_f.unlink()
                        seen_stems[rest] = (f, curr_id)
                    else:
                        if f.exists(): f.unlink()
                else:
                    seen_stems[rest] = (f, curr_id)
                    
                if len(curr_id) < 7:
                    new_name = f"{curr_id.zfill(7)}_{rest}"
                    target = yaml_dir / new_name
                    if not target.exists():
                        f.rename(target)
                        seen_stems[rest] = (target, curr_id.zfill(7))
    console_info("YAML Maintenance Completed. Standardised naming conventions enforced.")

def purge_orphans(prod_dir: Path, active_job_names: set) -> int:
    """Deletes YAML files and active run folders that are no longer referenced within the input dataset."""
    yaml_count, run_count = 0, 0
    
    # 1. Clean designated YAML structures
    yaml_dir = prod_dir / "2_Boltz2_YAML_Configs"
    if yaml_dir.exists():
        for f in yaml_dir.glob("*.yaml"):
            if f.stem not in active_job_names:
                f.unlink()
                yaml_count += 1
            
    # 2. Remove all active prediction run folders not mathematically associated with active jobs
    runs_dir = prod_dir / "4_Prediction_Jobs"
    if runs_dir.exists():
        for item in runs_dir.iterdir():
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


# ===============================================================================
# SECTION 5: BIOINFORMATICS (ALIGNMENT & SUPERIMPOSITION)
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 5.1: Caching Framework
# -------------------------------------------------------------------------------
def load_cached_alignments(csv_path: Path):
    """Loads previously calculated alignments into memory and ensures data types are strictly preserved."""
    global CACHED_ALIGNMENTS
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
    """Appends a newly computed alignment directly to the designated on-disk cache file."""
    if "protein" in data:
        CACHED_ALIGNMENTS[data["protein"]] = data
        
    keys = sorted(data.keys())
    file_exists = csv_path.exists()
    try:
        with open(csv_path, 'a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
            if not file_exists:
                writer.writeheader()
            writer.writerow(data)
    except Exception as e:
        if logger: logger.debug(f"Cache synchronisation deferred: {e}")

def save_alignment_cache_final(csv_path: Path):
    """Safely flushes the complete alignment cache to disk upon execution completion."""
    if not CACHED_ALIGNMENTS: return
    try:
        all_keys = set()
        for r in CACHED_ALIGNMENTS.values():
            all_keys.update(r.keys())
        keys = sorted(list(all_keys))
        
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
            writer.writeheader()
            for row in CACHED_ALIGNMENTS.values():
                writer.writerow(row)
        console_info(f"Canonical Alignment Cache persisted to {csv_path.name}")
    except Exception as e:
        console_info(f"Warning: Failed to save final alignment cache state: {e}")

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
        'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'Q':'GLN', 'E':'GLU',
        'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE',
        'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL', 'B':'ASX',
        'Z':'GLX', 'X':'UNK', '-': 'GAP'
    }

    # -------------------------------------------------------------------------------
    # Sub-Step 5.2.2: Configure Alignment Parameters
    # -------------------------------------------------------------------------------
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score   = CFG.ALIGN_OPEN_GAP_SCORE
    aligner.extend_gap_score = CFG.ALIGN_EXTEND_GAP_SCORE
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    
    if not target_seq or not target_seq.strip():
        return {}, {"error": "Target Sequence string is empty", "identity_pct": 0.0}, {}, "NA"

    target_seq = str(target_seq)

    # -------------------------------------------------------------------------------
    # Sub-Step 5.2.3: Execute Alignment Algorithm
    # -------------------------------------------------------------------------------
    alignments = aligner.align(REF_SEQUENCE_STR, target_seq)
    if not alignments: 
        return {}, {"error": "No viable alignment mapping identified", "identity_pct": 0.0}, {}, "NA"
        
    aln = alignments[0]
    
    if out_aln_path:
        out_aln_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_aln_path, "w") as f: f.write(format(aln))
        GLOBAL_STATS["total_alignments_generated"] += 1

    # -------------------------------------------------------------------------------
    # Sub-Step 5.2.4: Extract Statistics & Map Coordinate Frames
    # -------------------------------------------------------------------------------

    n_match, length = 0, 0
    seq1_str, seq2_str = str(aln[0]), str(aln[1])
    for c1, c2 in zip(seq1_str, seq2_str):
        length += 1
        if c1 == c2 and c1 != '-': n_match += 1
            
    pid = (n_match / length) * 100 if length > 0 else 0
    
    stats = {
        "align_score": aln.score, "identity_pct": round(pid, 2),
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
        if r_char != '-': ref_cursor += 1
        if t_char != '-': tgt_cursor += 1
        
        if r_char != '-':
            ref_aa_3 = AA_1_TO_3.get(r_char, 'UNK')
            if t_char != '-':
                tgt_aa_3 = AA_1_TO_3.get(t_char, 'UNK')
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
    
    for name, info in REF_ACTIVE_SITE_MAP.items():
        ref_idx = info["id"]
        tgt_idx = ref_to_tgt_idx.get(ref_idx, None)
        final_idx_map[name] = tgt_idx
        if tgt_idx is not None:
            aa_1 = tgt_seq_map.get(tgt_idx, 'X')
            aa_3 = AA_1_TO_3.get(aa_1, 'UNK')
            final_resname_map[name] = f"{aa_3}{tgt_idx}"
        else:
            final_resname_map[name] = "MISSING"
        
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


# ===============================================================================
# SECTION 6: CHEMINFORMATICS & PHYSICS ENGINE
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 6.1: Spatial Mathematical Utilities
# -------------------------------------------------------------------------------
# (Local distance and calculate_angle definitions removed; imported from central utils instead)
    
def get_plane_normal(atoms):
    """Calculates the optimal best-fit plane normal vector for Pi-stacking analysis using Singular Value Decomposition (SVD)."""
    if not atoms: return None, None
    coords = np.array([[a['x'], a['y'], a['z']] for a in atoms])
    centroid = np.mean(coords, axis=0)
    centred = coords - centroid
    u, s, vh = np.linalg.svd(centred)
    normal = vh[2, :] 
    return centroid, normal

def load_structure_safe(path: Path):
    """Safely loads standard PDB and mmCIF files using Gemmi's robust parsing capabilities."""
    try:
        if path.suffix == '.pdb': return gemmi.read_structure(str(path))
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
    params.randomSeed = 0xF00D
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
    U, S, Vt = np.linalg.svd(P_centred.T @ Q_centred)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = Qc - (R @ Pc)
    return R, t

def map_mmcif_to_rdkit(mmcif_atoms, rdkit_mol):
    """
    Intelligently maps Boltz-predicted atoms to RDKit atoms, effectively linking 
    physical structural positions to theoretical charge constraints.
    """
    if not mmcif_atoms or rdkit_mol is None: return {}
    mm_coords = np.array([[a["x"],a["y"],a["z"]] for a in mmcif_atoms])
    rd_coords = np.array(rdkit_coords_list(rdkit_mol))
    try:
        m = min(len(rd_coords), len(mmcif_atoms))
        if m >= 3:
            R, t = kabsch_transform(rd_coords[:m], mm_coords[:m])
            rd_rot = (R @ rd_coords.T).T + t
            cost = np.linalg.norm(mm_coords[:,None,:] - rd_rot[None,:,:], axis=2)
        else:
            cost = np.linalg.norm(mm_coords[:,None,:] - rd_coords[None,:,:], axis=2)
    except Exception:
        cost = np.linalg.norm(mm_coords[:,None,:] - rd_coords[None,:,:], axis=2)

    dynamic_penalty = (np.max(cost) * CFG.ALIGN_DYNAMIC_PENALTY_FACTOR
                       if cost.size > 0 else CFG.ALIGN_DYNAMIC_PENALTY_FALLBACK)
    for i, mm_atom in enumerate(mmcif_atoms):
        for j, rd_atom in enumerate(rdkit_mol.GetAtoms()):
            if mm_atom["element"] != rd_atom.GetSymbol():
                cost[i, j] += dynamic_penalty  

    row_ind, col_ind = linear_sum_assignment(cost)
    mapping = {}
    names = [a["atom_name"] for a in mmcif_atoms]
    for i, j in zip(row_ind, col_ind):
        if cost[i,j] <= COORD_MATCH_MAX_DIST: mapping[names[i]] = int(j) 
    return mapping

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
    formal charge alone would miss every carboxylate/sulfonate; treating ANY
    neutral N/O/S as charged (the previous heuristic) over-counts ethers and
    amides. This recognises the actual ionisable groups instead:
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
    strictly adhering to Maestro threshold limits.
    """
    types = []
    pel = p_atom["element"]
    lel = l_atom["element"]
    
    # -------------------------------------------------------------------------------
    # Sub-Step 6.3.1: Evaluator Block
    # -------------------------------------------------------------------------------
    # 1. Hydrogen Bonds (D-A distance only; D-H-A angle requires explicit H atoms not present in Boltz-2 CIF)
    if pel in {"N", "O"} and lel in {"N", "O"}:
        if d <= HBOND_MAX_DIST_DA: types.append("hydrogen_bond")
             
    # 2. Hydrophobic Bonds
    if pel=="C" and lel=="C" and d <= HYDROPHOBIC_MAX_DIST: types.append("hydrophobic")
        
    # 3. Salt Bridges (opposite-charge groups within SALT_BRIDGE_MAX_DIST).
    #    Ligand ionisability is resolved by functional group, not by the SMILES
    #    protonation state (see _ligand_ionisable), so neutral-form carboxylates
    #    still register while ethers/amides do not.
    if d <= SALT_BRIDGE_MAX_DIST:
        rd_atom = None
        idx = mm_to_rd.get(l_atom["atom_name"])
        if idx is not None and rdkit_mol:
            try: rd_atom = rdkit_mol.GetAtomWithIdx(idx)
            except Exception: rd_atom = None
        lig_ion = _ligand_ionisable(rd_atom)
        is_salt_bridge = False

        # Case A: cationic residue (Arg/Lys/His) vs anionic ligand group
        if p_atom["resname"] in POSITIVE_RES and lig_ion == "anion":
            is_salt_bridge = True

        # Case B: anionic residue (Asp/Glu) vs cationic ligand group
        if p_atom["resname"] in NEGATIVE_RES and lig_ion == "cation":
            is_salt_bridge = True

        if is_salt_bridge: types.append("salt_bridge")
        
    # 4. Halogen Contact Analysis & Specific Fluorous Metrics
    if lel.upper() == "F" and d <= HALOGEN_MAX_DIST:
        types.append("halogen_contact")
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
                    d = {'x':atom.pos.x, 'y':atom.pos.y, 'z':atom.pos.z, 
                         'element': (atom.element.name or atom.name[0]), 
                         'chain': chain.name, 'resname': res.name, 'resseq': res.seqid.num, 'atom_name': atom.name}
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
    for a in prot_atoms: res_map[(a['chain'], a['resseq'])].append(a)
    
    for (ch, seq), atoms in res_map.items():
        resname = atoms[0]['resname']
        if resname in RES_PROPS["Aromatic"]:
            ring_atoms = [a for a in atoms if a['atom_name'] in AROMATIC_RING_ATOMS]
            if len(ring_atoms) >= 5:
                cent, norm = get_plane_normal(ring_atoms)
                if cent is not None: prot_rings.append({'centroid': cent, 'normal': norm})

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
                    found = next((a for a in lig_atoms if a['atom_name'] == cif_name), None)
                    if found: cif_ring_atoms.append(found)
            
            if len(cif_ring_atoms) >= 5:
                cent, norm = get_plane_normal(cif_ring_atoms)
                if cent is not None: lig_rings.append({'centroid': cent, 'normal': norm})

    # -------------------------------------------------------------------------------
    # Sub-Step 6.5.3: Spatial Interference Verification Checks
    # -------------------------------------------------------------------------------
    for p_ring in prot_rings:
        p_c, p_n = p_ring['centroid'], p_ring['normal']
        
        # A. Pi-Stacking Logic Validation Framework
        for l_ring in lig_rings:
            l_c, l_n = l_ring['centroid'], l_ring['normal']
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
                l_pos = np.array([la['x'], la['y'], la['z']])
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
        else: lig_heavy_atoms = len([a for a in lig if a['element'].upper() != 'H'])
        
        # Structural safeguards implemented for single-atom ions (e.g. F-, Cl-)
        if lig_heavy_atoms == 0: lig_heavy_atoms = 1

        # -------------------------------------------------------------------------------
        # Sub-Step 6.5.5: Residue Demographics Census
        # -------------------------------------------------------------------------------
        census = {"Total_Residues": 0, "Total_Polar_Residues": 0, 
                  "Total_NonPolar_Residues": 0, "Total_Pos_Residues": 0, 
                  "Total_Neg_Residues": 0, "Total_Aromatic_Residues": 0}
        unique_residues = set()
        for p in prot: unique_residues.add((p['chain'], p['resseq'], p['resname']))
        census["Total_Residues"] = len(unique_residues)
        res_atom_map = defaultdict(list)
        for p in prot: res_atom_map[(p['chain'], p['resseq'], p['resname'])].append(p)

        for (chain, seq, resname), atoms in res_atom_map.items():
            if resname in RES_PROPS["Polar"]: census["Total_Polar_Residues"] += 1
            if resname in RES_PROPS["NonPolar"]: census["Total_NonPolar_Residues"] += 1
            if resname in RES_PROPS["Pos"]: census["Total_Pos_Residues"] += 1
            if resname in RES_PROPS["Neg"]: census["Total_Neg_Residues"] += 1
            if resname in RES_PROPS["Aromatic"]: census["Total_Aromatic_Residues"] += 1

        counts = defaultdict(int)
        dists = []
        rows = []
        mapped_count = 0
        total_f_atoms = len([a for a in lig if a['element'] == 'F'])
        interacting_f_set = set()
        
        # -------------------------------------------------------------------------------
        # Sub-Step 6.5.6: Optimised cKDTree Interaction Evaluation Loop
        # -------------------------------------------------------------------------------
        if prot and lig:
            prot_coords = [[p['x'], p['y'], p['z']] for p in prot]
            lig_coords  = [[l['x'], l['y'], l['z']] for l in lig]
            tree = cKDTree(prot_coords)
            indices_list = tree.query_ball_point(lig_coords, r=6.0)
            
            for i, p_indices in enumerate(indices_list):
                la = lig[i]
                l_pos = lig_coords[i]
                for p_idx in p_indices:
                    pa = prot[p_idx]
                    p_pos = prot_coords[p_idx]
                    d = distance(p_pos, l_pos)
                    labels = classify_pair(pa, la, d, rd_mol, mm_map)
                    if labels:
                        dists.append(d)
                        for l in labels: counts[l] += 1
                        if "fluorine_contact" in labels: interacting_f_set.add(la['atom_name']) 
                        rows.append({
                            "protein_chain": pa['chain'], "protein_resname": pa['resname'],
                            "protein_resseq": pa['resseq'], "protein_atom": pa['atom_name'],
                            "ligand_atom": la['atom_name'], "distance_A": round(d, 4),
                            "interaction_types": ";".join(labels)
                        })
                        if mm_map.get(la['atom_name']) is not None: mapped_count += 1
        
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


# ===============================================================================
# SECTION 7: GEOMETRIC ANALYSIS (CATALYSIS) & CALIBRATION
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 7.1: Reference Data Handling (Calibration Logic)
# -------------------------------------------------------------------------------
# --- Ramachandran helpers (relocated to 00_03_Project_Utils_FAcDs.py) ---------



def _log_ramachandran_table(stats_ref: Dict, stats_con: Dict,
                             label_ref: str = '3R3U (Crystal)',
                             label_con: str = 'DeHa4 Control (Modelled)') -> None:
    """Log a side-by-side Ramachandran statistics table to the console."""
    w = 28
    hdr = f"{'Category':<14}  {label_ref:>{w}}  {label_con:>{w}}"
    sep = '-' * len(hdr)
    console_info(f"\n  {'Ramachandran Statistics':^{len(hdr)}}")
    console_info(f"  {sep}")
    console_info(f"  {hdr}")
    console_info(f"  {sep}")
    for cat in ('Favored', 'Allowed', 'Outlier'):
        n_r = stats_ref['counts'][cat]
        p_r = stats_ref['pct'][cat]
        n_c = stats_con['counts'][cat]
        p_c = stats_con['pct'][cat]
        row = f"  {cat:<14}  {f'{p_r:5.1f}%  ({n_r:>4} res)':>{w}}  {f'{p_c:5.1f}%  ({n_c:>4} res)':>{w}}"
        console_info(row)
    console_info(f"  {sep}")
    t_r, t_c = stats_ref['total'], stats_con['total']
    console_info(f"  {'Total residues':<14}  {str(t_r):>{w}}  {str(t_c):>{w}}")
    console_info(f"  {sep}\n")

# --- End Ramachandran helpers ------------------------------------------------


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
        "rama_ref": None,
        "rama_con": None,
        "png_names": [],
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
            if hasattr(gemmi.SupSelect, 'Ca'):
                sup_select = gemmi.SupSelect.Ca
            elif hasattr(gemmi.SupSelect, 'CA'):
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

        if trust_score <= 2.0:
            verdict_colour = ConsoleColours.OKGREEN
            verdict_label  = "High Confidence"
            verdict_note   = "RMSD < 2.0 Å — Structures are reliable."
        elif trust_score <= 3.0:
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

        # --- RAMACHANDRAN ANALYSIS ---
        try:
            rama_ref  = compute_ramachandran_angles(st_ref)
            rama_con  = compute_ramachandran_angles(st_con)
            stats_ref = _rama_stats(rama_ref)
            stats_con = _rama_stats(rama_con)

            calib_data["rama_ref"] = stats_ref
            calib_data["rama_con"] = stats_con

            if not silent:
                _log_ramachandran_table(stats_ref, stats_con,
                                        label_con=f"{model_label} (Modelled)")

            ref_dir   = ref_pdb.parent
            safe_label = re.sub(r"[^A-Za-z0-9]", "_", model_label)
            path_ref_png = ref_dir / "Ramachandran_3R3U_Crystal.png"
            path_con_png = ref_dir / f"Ramachandran_{safe_label}.png"
            path_cmp_png = ref_dir / f"Ramachandran_{safe_label}_Comparison.png"

            save_ramachandran_plot(rama_ref, '3R3U Crystal Structure', path_ref_png, dpi=CFG.VIS_FIGURE_DPI)
            save_ramachandran_plot(rama_con, f'{model_label} (Boltz-2 Model)', path_con_png, dpi=CFG.VIS_FIGURE_DPI)
            save_ramachandran_comparison(rama_ref, rama_con,
                                         '3R3U (Crystal)', f'{model_label} (Boltz-2)',
                                         path_cmp_png, dpi=CFG.VIS_FIGURE_DPI)

            calib_data["png_names"] = [path_ref_png.name, path_con_png.name, path_cmp_png.name]
            _ci(f"  > Ramachandran plots saved for {path_ref_png.name}, {path_con_png.name} and {path_cmp_png.name}")
        except Exception as _rama_err:
            _ci(f"  [Warning] Ramachandran analysis failed: {_rama_err}")

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
                ref_pdb_id = site.get('pdb_id', site['id'])
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
                for _offset in (0, 1, -1, 2, -2):
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
                ref_pdb_id = site.get('pdb_id', site['id'])
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

                min_dist   = 6.0
                best_match = None
                best_match_resname = ""

                for res_inner in all_con_residues:
                    if res_inner.name[:3] != site['res'][:3]:
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
            _ci(f"    ! Critical Warning: Specified residues were not definitively mapped.")
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
    all_data  = deha4_data + r3u_data   # 6 items total

    W = 110
    sep_h  = "─" * W
    sep_eq = "=" * W

    # ── Header ──────────────────────────────────────────────────────────────────
    console_info(f"\n  {sep_eq}")
    console_info(f"  {'CONTROL CALIBRATION SUMMARY — DeHa4 × 3R3U  (3 Ligands Each)':^{W}}")
    console_info(f"  {sep_eq}\n")

    # ── Trust Scores ─────────────────────────────────────────────────────────────
    # Trust score = Cα RMSD (Å) between Boltz-2 prediction and 3R3U crystal structure.
    # Lower = better alignment. <2.0 Å = High Confidence, 2–3 Å = Moderate, >3 Å = Low.
    _ts_col = 16
    _ts_sep = "─" * (14 + 3 * (_ts_col + 2) + 5 + 3 * (_ts_col + 2))
    console_info(f"  Structural Trust Score  (Cα RMSD vs 3R3U crystal; <2.0 Å = high confidence)")
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
        col = ConsoleColours.OKGREEN if v <= 2.0 else (ConsoleColours.WARNING if v <= 3.0 else ConsoleColours.FAIL)
        ts_row += f"  {col}{f'{v:.3f} Å'.rjust(_ts_col)}{ConsoleColours.ENDC}"
    ts_row += "   │"
    for d in r3u_data:
        v   = d.get("trust_score", 99.0)
        col = ConsoleColours.OKGREEN if v <= 2.0 else (ConsoleColours.WARNING if v <= 3.0 else ConsoleColours.FAIL)
        ts_row += f"  {col}{f'{v:.3f} Å'.rjust(_ts_col)}{ConsoleColours.ENDC}"
    console_info(ts_row)
    # Confidence verdict row
    verd_row = f"  {'Confidence':<14}"
    for d in list(deha4_data) + ["│"] + list(r3u_data):
        if d == "│":
            verd_row += "   │"
            continue
        v = d.get("trust_score", 99.0)
        lbl = "High" if v <= 2.0 else ("Moderate" if v <= 3.0 else "Low")
        col = ConsoleColours.OKGREEN if v <= 2.0 else (ConsoleColours.WARNING if v <= 3.0 else ConsoleColours.FAIL)
        verd_row += f"  {col}{lbl:>{_ts_col}}{ConsoleColours.ENDC}"
    console_info(verd_row)
    console_info(f"  {_ts_sep}\n")

    # ── Ramachandran Table ───────────────────────────────────────────────────────
    col_labels = ["Crystal"] + [f"DeHa4_{s}" for s in lig_short] + [f"3R3U_{s}" for s in lig_short]
    col_data   = [None] + [d.get("rama_con") for d in deha4_data] + [d.get("rama_con") for d in r3u_data]
    # crystal stats: take from first deha4 item (which has rama_ref from crystal)
    crystal_stats = deha4_data[0].get("rama_ref") if deha4_data else None
    if crystal_stats:
        col_data[0] = crystal_stats

    if any(d is not None for d in col_data):
        cw = 13
        # Determine winner: highest Favored%, tie-broken by lowest Outlier%.
        # Crystal is index 0 and excluded from the competition (it is the reference).
        _ranked_cols = sorted(
            [(i, cd) for i, cd in enumerate(col_data) if i > 0 and cd is not None],
            key=lambda x: (x[1]['pct']['Favored'], -x[1]['pct']['Outlier']),
            reverse=True
        )
        _winner_idx = _ranked_cols[0][0] if _ranked_cols else -1

        console_info(f"  Ramachandran Backbone Quality — Crystal + All 6 Control Structures")
        console_info(f"  (Higher Favored% and lower Outlier% = better backbone geometry)")
        console_info(f"  {sep_h}")
        hdr = f"  {'Category':<10}" + "".join(f"  {lb:>{cw}}" for lb in col_labels)
        console_info(hdr)
        console_info(f"  {sep_h}")
        for cat in ("Favored", "Allowed", "Outlier"):
            row = f"  {cat:<10}"
            for i, cd in enumerate(col_data):
                if cd is None:
                    row += f"  {'—':>{cw}}"
                else:
                    pct = cd['pct'][cat]
                    n   = cd['counts'][cat]
                    cell = f'{pct:.1f}%({n})'
                    if i == _winner_idx and cat == "Favored":
                        row += f"  {ConsoleColours.OKGREEN}{cell:>{cw}}{ConsoleColours.ENDC}"
                    else:
                        row += f"  {cell:>{cw}}"
            console_info(row)
        console_info(f"  {sep_h}")
        row_tot = f"  {'Total res':<10}"
        for cd in col_data:
            row_tot += f"  {(str(cd['total']) if cd else '—'):>{cw}}"
        console_info(row_tot)
        console_info(f"  {sep_h}")
        # Winner row
        _winner_row = f"  {'Best model':<10}"
        for i, cd in enumerate(col_data):
            if i == 0:
                _winner_row += f"  {'(reference)':>{cw}}"
            elif i == _winner_idx:
                _winner_row += f"  {ConsoleColours.OKGREEN}{'← Winner':>{cw}}{ConsoleColours.ENDC}"
            else:
                _winner_row += f"  {'':>{cw}}"
        console_info(_winner_row)
        if _winner_idx >= 0:
            _wname = col_labels[_winner_idx]
            _wfav  = col_data[_winner_idx]['pct']['Favored']
            _wout  = col_data[_winner_idx]['pct']['Outlier']
            console_info(f"  {ConsoleColours.OKGREEN}  Best backbone geometry: {_wname}  "
                         f"(Favored {_wfav:.1f}%, Outlier {_wout:.1f}%){ConsoleColours.ENDC}")
        console_info(f"  {sep_h}\n")

    # ── Plots Saved ──────────────────────────────────────────────────────────────
    all_pngs_deha4 = []
    all_pngs_r3u   = []
    crystal_png    = None
    for d in deha4_data:
        pngs = d.get("png_names", [])
        if pngs:
            if not crystal_png:
                crystal_png = pngs[0]
            all_pngs_deha4.extend(pngs[1:])
    for d in r3u_data:
        pngs = d.get("png_names", [])
        all_pngs_r3u.extend(pngs[1:])

    console_info("  Ramachandran Plots Saved:")
    if crystal_png:
        console_info(f"    [Crystal]  {crystal_png}")
    for p in all_pngs_deha4:
        console_info(f"    [DeHa4 ]  {p}")
    for p in all_pngs_r3u:
        console_info(f"    [3R3U  ]  {p}")
    console_info("")

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


def analyse_candidate_structure(target_cif: Path, control_cif: Path, control_map: Dict[str, int]) -> Dict[str, Any]:
    """
    Superimposes the candidate structure onto the DeHa4 control model.
    Computes the likelihood score, RMSD, and mechanistic fingerprint by comparing targeted structural pockets.
    """
    result = {
        "Active_Site_RMSD": 99.0, 
        "Mechanistic_Fingerprint_Score": 0.0,
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
                if hasattr(gemmi.SupSelect, 'Ca'):
                    sup_select = gemmi.SupSelect.Ca
                elif hasattr(gemmi.SupSelect, 'CA'):
                    sup_select = gemmi.SupSelect.CA
                else:
                    sup_select = gemmi.SupSelect.All
            except Exception:
                sup_select = gemmi.SupSelect.All
            
            # Formally constructs a robust superposition mapping matrix
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
            # Sub-Step 7.1.2: Extract Structural Mechanistic Fingerprint Data
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
                
                for res in tar_chain:
                    tar_ca = res.find_atom("CA", "*")
                    if tar_ca:
                        d = con_pos.dist(tar_ca.pos)
                        if d < min_d: 
                            min_d = d
                            mapped_res = res
                
                # Accrue internal RMSD values representing core active site geometry
                if min_d < 5.0:
                    rmsd_sq_sum += (min_d * min_d)
                    count += 1
                else:
                    rmsd_sq_sum += 25.0
                    count += 1
                    
                # Validate exact spatial chemical identity requirements
                if mapped_res and min_d < 4.0:
                    res_name = mapped_res.name
                    total_checks += 1
                    ref_role = REF_ACTIVE_SITE_MAP[key]['role']
                    
                    if ref_role == "Fluoride_Cradle":
                        if res_name in FINGERPRINT_GROUPS["AROMATIC"]: hits += 1; has_halide_cradle = True
                    elif ref_role == "Carboxylate_Clamp":
                        if res_name in FINGERPRINT_GROUPS["POSITIVE"]: hits += 1; has_carb_clamp = True
                    elif ref_role == "Fluorine_Stabiliser":
                         if res_name in FINGERPRINT_GROUPS["AROMATIC"] or res_name in FINGERPRINT_GROUPS["POSITIVE"]: hits += 1
                    elif ref_role == "Nucleophile":
                         if res_name == "ASP": hits += 1
                    elif ref_role == "Acid_Catalyst":
                         if res_name == "ASP" or res_name == "GLU": hits += 1
                         
            result["Active_Site_RMSD"] = math.sqrt(rmsd_sq_sum / count) if count > 0 else 99.0
            result["Mechanistic_Fingerprint_Score"] = round(hits / max(total_checks, 1), 2)
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
            return f"Failure: Overall sequence alignment proved unreliable (ID {identity}% < 25%). Structure is likely invalid."
    return f"{meaning} ({constraint})"

def calculate_sn2_metrics(asp_atoms, lig_atoms, rd_mol=None, mm_map=None) -> Tuple[float, float, int, Any, Any, Any]:
    """
    Calculates the SN2 attack angle and the SN2 backside-attack trajectory deviation.
    An ideal SN2 backside attack strictly requires an angle of approximately 180 degrees.
    Deviation quantifies the perpendicular distance measured from the ideal C-X vector.
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
                if c.pos.dist(x.pos) < 2.2:
                    valid_cx_pairs.append((c, x))

    if not valid_cx_pairs: return 0.0, 999.0, 0, None, None, None

    # Pick the C-X pair whose carbon is closest to a nucleophile oxygen
    min_dist_O_C = 999.0
    best_O, best_C = None, None
    for o in asp_oxygens:
        for c, x in valid_cx_pairs:
            d = o.pos.dist(c.pos)
            if d < min_dist_O_C:
                min_dist_O_C = d
                best_O, best_C = o, c

    if not best_O or not best_C: return 0.0, 999.0, 0, None, None, None

    # Over the fluorines bonded to the chosen best_C, select the one
    # that maximises the Nu-C-F angle (backside attack directionality)
    best_X = None
    max_ang = -1.0
    for c, x in valid_cx_pairs:
        if c == best_C:
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
        cross_prod = np.cross(vec_cnu, vec_cx)
        deviation = np.linalg.norm(cross_prod) / (np.linalg.norm(vec_cx) + 1e-6)
    except Exception: deviation = 999.0
    
    # Teflon Shield Calculation: Count adjacent fluorines that sterically clash with Aspartate catalytic oxygen
    teflon_clashes = 0
    for f in lig_halogens:
        if f != best_X and f.element.name == "F":
            if best_O.pos.dist(f.pos) <= 2.5:
                teflon_clashes += 1

    # ── Auxiliary (NON-GATING) reference geometries — see utils docstrings ─────
    # Bürgi–Dunitz: computed on the substrate CARBONYL carbon (the –COO⁻ carbon),
    #   NOT the sp³ SN2 carbon — a nucleophile/clamp pre-organisation descriptor.
    # Flippin–Lodge: in-plane offset of the nucleophile at the attacked carbon,
    #   defined only when ≥2 heavy spectator substituents exist (else NaN).
    aux = {"burgi_dunitz_angle": float("nan"), "flippin_lodge_offset": float("nan")}
    if rd_mol and mm_map:
        try:
            inv_map = {v: k for k, v in mm_map.items()}
            _name2lig = {a.name: a for a in lig_atoms}
            # Bürgi–Dunitz on the carboxylate carbon
            for _rc in rd_mol.GetAtoms():
                if _rc.GetSymbol() != "C":
                    continue
                _o_nbrs = [n for n in _rc.GetNeighbors() if n.GetSymbol() == "O"]
                if len(_o_nbrs) < 2:
                    continue
                _cn = inv_map.get(_rc.GetIdx())
                _on = inv_map.get(_o_nbrs[0].GetIdx())
                _cC, _cO = _name2lig.get(_cn), _name2lig.get(_on)
                if _cC and _cO:
                    aux["burgi_dunitz_angle"] = round(
                        calculate_burgi_dunitz(best_O.pos, _cC.pos, _cO.pos), 1)
                break
            # Flippin–Lodge in-plane offset at the SN2 carbon (needs 2 heavy spectators)
            _c_rd_idx = mm_map.get(best_C.name)
            if _c_rd_idx is not None:
                _c_rd = rd_mol.GetAtomWithIdx(_c_rd_idx)
                _spect = [n for n in _c_rd.GetNeighbors()
                          if n.GetSymbol() not in ("H",)
                          and inv_map.get(n.GetIdx()) != best_X.name]
                _spect_pos = [_name2lig[inv_map[n.GetIdx()]].pos
                              for n in _spect
                              if inv_map.get(n.GetIdx()) in _name2lig]
                if len(_spect_pos) >= 2:
                    _fl = calculate_flippin_lodge(
                        best_O.pos, best_C.pos, best_X.pos, _spect_pos[0], _spect_pos[1])
                    if _fl < 990.0:
                        aux["flippin_lodge_offset"] = round(_fl, 1)
        except Exception:
            pass

    return angle, deviation, teflon_clashes, best_X.pos, best_C, aux

def sigmoid(x: float, k: float = 1.0, x0: float = 0.0) -> float:
    """Standard sigmoid function for soft-thresholding."""
    try:
        return 1 / (1 + math.exp(-k * (x - x0)))
    except OverflowError:
        return 0.0 if x - x0 < 0 else 1.0

def check_catalytic_geometry(cif_path: Path, mapped_sites: Dict[str, int], smiles_str: str = "") -> Dict[str, Any]:
    """
    Evaluates exact distance topologies within the active site against strictly defined tiers, 
    subsequently assigning a final categorical degrader tier.
    """
    results = {f"dist_{k}": 999.0 for k in REF_ACTIVE_SITE_MAP.keys()}
    results.update({"dist_ASP110": 999.0, "dist_HIS277": 999.0, "dist_ASP134": 999.0, 
                    "sn2_attack_angle": 0.0, "sn2_trajectory_dev": 999.0})

    try:
        # -------------------------------------------------------------------------------
        # Sub-Step 7.2.1: Formal Load of Positional Coordinate Data
        # -------------------------------------------------------------------------------
        mmcif_prot, mmcif_lig = load_atoms_from_structure(cif_path)
        rd_mol = rdkit_mol_from_smiles_with_3d(smiles_str) if smiles_str else None
        mm_map = map_mmcif_to_rdkit(mmcif_lig, rd_mol) if rd_mol else {}
        
        doc = gemmi.cif.read_file(str(cif_path))
        st = gemmi.make_structure_from_block(doc.sole_block())
        lig_atoms_obj, lig_coords = [], []
        site_atoms, site_atoms_obj = defaultdict(list), defaultdict(list)
        all_prot_atoms = []
        
        for model in st:
            for chain in model:
                for res in chain:
                    if res.name not in STANDARD_AA and res.name != "HOH":
                        for atom in res: 
                            lig_atoms_obj.append(atom)
                            lig_coords.append(atom.pos)
                    else:
                        for atom in res:
                            all_prot_atoms.append({
                                'pos': atom.pos, 'resname': res.name, 
                                'element': atom.element.name, 'atom': atom.name,
                                'atom_name': atom.name
                            })
                    if res.seqid.num in mapped_sites.values():
                        for name, idx in mapped_sites.items():
                            if idx == res.seqid.num:
                                for atom in res: 
                                    site_atoms[name].append(atom.pos)
                                    site_atoms_obj[name].append(atom)
        
        if not lig_coords: 
                    results.update({
                        "catalytic_dist_A": 999.0, "degrader_tier": CFG.TIER_DECOY, "is_degrader": False,
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

        results["dist_ASP110"] = results["dist_Nuc"]
        results["dist_HIS277"] = results["dist_Base"]
        results["dist_ASP134"] = results["dist_Acid"]
        
        dist_nuc_base = 999.0
        dist_base_acid = 999.0
        if "Nuc" in site_atoms and "Base" in site_atoms:
            dist_nuc_base = min([p1.dist(p2) for p1 in site_atoms["Nuc"] for p2 in site_atoms["Base"]])
        if "Base" in site_atoms and "Acid" in site_atoms:
            dist_base_acid = min([p1.dist(p2) for p1 in site_atoms["Base"] for p2 in site_atoms["Acid"]])
            
        results["dist_nuc_base_internal"] = round(dist_nuc_base, 2)
        results["dist_base_acid_internal"] = round(dist_base_acid, 2)


        sn2_angle, sn2_dev, steric_clashes, target_x_pos, best_c_atom, sn2_aux = calculate_sn2_metrics(site_atoms_obj.get("Nuc", []), lig_atoms_obj, rd_mol, mm_map)
        results["sn2_attack_angle"] = round(sn2_angle, 1)
        results["sn2_trajectory_dev"] = round(sn2_dev, 2)
        results["teflon_shield_clashes"] = steric_clashes
        # Auxiliary (non-gating) reference geometries
        _aux = sn2_aux or {}
        results["burgi_dunitz_angle"]   = _aux.get("burgi_dunitz_angle", float("nan"))
        results["flippin_lodge_offset"] = _aux.get("flippin_lodge_offset", float("nan"))

        d_nuc = results["dist_Nuc"]
        angle  = results["sn2_attack_angle"]
        d_trp, d_tyr = results.get("dist_Stab_W", 999.0), results.get("dist_Stab_Y", 999.0)
        
        # --- DYNAMIC ACTIVE SITE STABILISATION METRICS ---
        # Checks atom-specific physical proximity to accurately validate established site stabilisation.
        dynamic_stabilised = False
        target_stabilise_centre = target_x_pos
        if not target_stabilise_centre and lig_coords: target_stabilise_centre = lig_coords[0]

        if target_stabilise_centre:
            for p_at in all_prot_atoms:
                if p_at['resname'] in POLAR_SIDECHAIN_ATOMS:
                     valid_atoms = POLAR_SIDECHAIN_ATOMS[p_at['resname']]
                     if p_at['atom_name'] in valid_atoms:
                        d_stab = p_at['pos'].dist(target_stabilise_centre)
                        if d_stab <= CFG.MECH_STAB_RADIUS: dynamic_stabilised = True; break

        stabilised = (d_trp <= CFG.MECH_STAB_RADIUS or d_tyr <= CFG.MECH_STAB_RADIUS or dynamic_stabilised)
        results["halide_stabilisation_score"] = 1.0 if stabilised else 0.0

        clamp_ok = (results["dist_Carb1"] <= CFG.MECH_CLAMP_RADIUS or results["dist_Carb2"] <= CFG.MECH_CLAMP_RADIUS)
        results["carboxylate_clamp_integrity"] = 1.0 if clamp_ok else 0.0
        results["sn2_alignment_score"] = angle
        
        # --- SOFT SCORING ENGINE ---
        # Uses sigmoid functions to avoid binary threshold "cliffs"
        s_nuc = sigmoid(d_nuc, k=-4.0, x0=CFG.NAC_DIST_STRICT)   # High score for < NAC_DIST_STRICT (3.2 Å)
        s_ang = sigmoid(angle, k=0.15, x0=CFG.NAC_ANGLE_STRICT)  # High score for > NAC_ANGLE_STRICT (155°)
        s_int = sigmoid(dist_nuc_base, k=-2.0, x0=CFG.SOFT_NB_MIDPOINT) * sigmoid(dist_base_acid, k=-2.0, x0=CFG.SOFT_BA_MIDPOINT)
        soft_score = (s_nuc * 0.4) + (s_ang * 0.3) + (s_int * 0.3)
        results["soft_catalytic_score"] = round(soft_score, 3)

        # MECHANISTIC FINGERPRINT DATA SCORING
        # Scores specific physical requirements including Halide Stabilisation and Carboxylate Clamps.
        mech_score = 0.0
        if d_nuc <= CFG.NAC_DIST_STRICT: mech_score += 0.2
        if dist_nuc_base <= CFG.MECH_NB_GATE: mech_score += 0.1
        if dist_base_acid <= CFG.MECH_BA_GATE: mech_score += 0.1
        if clamp_ok: mech_score += 0.2 
        if stabilised: mech_score += 0.4
        if steric_clashes > 0: mech_score -= (1.0 * steric_clashes)

        # ELECTRONIC "DEAD-END" PENALTY EVALUATION
        terminal_f_count = 0
        if rd_mol and mm_map and best_c_atom:
            target_rd_idx = mm_map.get(best_c_atom.name)
            if target_rd_idx is not None:
                try:
                    target_rd_atom = rd_mol.GetAtomWithIdx(target_rd_idx)
                    for nbr in target_rd_atom.GetNeighbors():
                        if nbr.GetSymbol() == "F":
                            terminal_f_count += 1
                except Exception: pass
                
        # gem-difluoride (CF2 warhead) penalty: such carbons follow a hindered
        # SN2 path and engage a documented FAcD inhibition pathway.
        # Ref: Jansen, van Beers & Mayer (2026) Angew Chem Int Ed 65:e202524234,
        #      https://doi.org/10.1002/anie.202524234 (gem-difluoride inhibition).
        if terminal_f_count == 2:
            mech_score -= 0.5

        results["mechanistic_score"] = round(mech_score, 2)

        # DEEP-POCKET SIZE EXCLUSION (MAINCHAIN CLASHING METRICS)
        tail_clash_ratio = 0.0
        if rd_mol and mm_map and best_c_atom:
            target_rd_idx = mm_map.get(best_c_atom.name)
            if target_rd_idx is not None:
                distances = Chem.rdmolops.GetDistanceMatrix(rd_mol)
                tail_rd_indices = [i for i, d in enumerate(distances[target_rd_idx]) if d > 3]
                inv_map = {v: k for k, v in mm_map.items()}
                tail_cif_names = {inv_map.get(idx) for idx in tail_rd_indices if inv_map.get(idx)}
                tail_pos_list = [a.pos for a in lig_atoms_obj if a.name in tail_cif_names]
                
                if tail_pos_list:
                    clash_count = 0
                    bb_atoms = [pa for pa in all_prot_atoms if pa['atom_name'] in ["N", "CA", "C", "O"]]
                    for tpos in tail_pos_list:
                        is_clashing = False
                        for ba in bb_atoms:
                            if tpos.dist(ba['pos']) < 2.2:
                                is_clashing = True
                                break
                        if is_clashing: clash_count += 1
                    tail_clash_ratio = clash_count / len(lig_atoms_obj)
        
        results["mainchain_clash_ratio"] = round(tail_clash_ratio, 3)

        tier, is_degrader, meaning, constraint = CFG.TIER_DECOY, False, "No significant documented interactions.", "Fail"
        
        # -------------------------------------------------------------------------------
        # Sub-Step 7.2.3: ULTRA-STRICT FAcD HIERARCHY LOGIC
        # -------------------------------------------------------------------------------
        
        if tail_clash_ratio > 0.15:
            tier, meaning, is_degrader = CFG.TIER_DECOY, f"Failed: Severe Mainchain Structural Clashing (>{int(tail_clash_ratio*100)}% of tail atoms).", False
        
        elif terminal_f_count >= 3:
            # Allow short-chain substrates (≤2 C, e.g. TFA: CF3-COO⁻) — FAcD attacks alpha-CF3
            # Block only when CF3 is a terminus within a longer fluorocarbon chain
            _mol_c_count = sum(1 for _a in rd_mol.GetAtoms() if _a.GetSymbol() == "C") if rd_mol else 0
            if _mol_c_count > 2:
                tier, meaning, is_degrader = CFG.TIER_DECOY, "Failed: Electronic Dead-End Identified (Targeted Carbon is functionally a Terminal CF3 group).", False

        # Tier 1 (Elite High-Performance pool)
        # Pruned efficiently for near-ideal geometry to reduce secondary molecular dynamics workload.
        elif d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_TOP] and dist_nuc_base <= CFG.TIER_NB_MAX[CFG.TIER_TOP] and dist_base_acid <= CFG.TIER_BA_MAX[CFG.TIER_TOP] and angle >= CFG.TIER_ANGLE_MIN[CFG.TIER_TOP] and mech_score >= CFG.TIER_MECH_MIN[CFG.TIER_TOP]:
            tier, meaning, is_degrader = CFG.TIER_TOP, "Elite-Grade Analysis: Near-perfect SN2 Trajectory demonstrating absolute anchor integrity.", True
        
        # Tier 2 (High Functional grade)
        elif d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_ORDER[1]] and dist_nuc_base <= CFG.TIER_NB_MAX[CFG.TIER_ORDER[1]] and dist_base_acid <= CFG.TIER_BA_MAX[CFG.TIER_ORDER[1]] and angle >= CFG.TIER_ANGLE_MIN[CFG.TIER_ORDER[1]] and mech_score >= CFG.TIER_MECH_MIN[CFG.TIER_ORDER[1]]:
            tier, meaning, is_degrader = CFG.TIER_ORDER[1], "Crystal-Grade Analysis: Ideal ground-state contact sequence identified with strong anchoring profiles.", True
        
        # Tier 3 (Functional grade)
        elif d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_ORDER[2]] and dist_nuc_base <= CFG.TIER_NB_MAX[CFG.TIER_ORDER[2]] and dist_base_acid <= CFG.TIER_BA_MAX[CFG.TIER_ORDER[2]] and angle >= CFG.TIER_ANGLE_MIN[CFG.TIER_ORDER[2]] and mech_score >= CFG.TIER_MECH_MIN[CFG.TIER_ORDER[2]]:
            tier, meaning, is_degrader = CFG.TIER_ORDER[2], "Functional Analysis: Nucleophile located in tight contact, accompanied by acceptable target attack angles.", True
        
        # Tier 4 (Marginal grade)
        elif d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_ORDER[3]] and angle >= CFG.TIER_ANGLE_MIN[CFG.TIER_ORDER[3]]:
            tier, meaning, is_degrader = CFG.TIER_ORDER[3], "Marginal Analysis: Nucleophile indicates loose contact profiles paired with a marginal target attack angle.", True
        
        # Tier 5 (Loose/Non-functional grade)
        elif d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_ORDER[4]]:
            tier, meaning, is_degrader = CFG.TIER_ORDER[4], "Loose Alignment: Routine proximity searches identified the ligand, but orientation metrics strictly fail defined SN2 physical requirements.", False
        
        # Tier_4/Tier_5_Decoy
        else:
            tier = CFG.TIER_POOR if d_nuc <= CFG.TIER_NUC_DIST[CFG.TIER_POOR] else CFG.TIER_DECOY
            meaning = "Failed Analysis: Positional distance or established angle physically violates FAcD catalytic structural requirements."
            is_degrader = False
            
        # Derive constraint_check from the assigned tier — Tier_1/Tier_2 tiers pass,
        # Good is a partial pass (in pocket but wrong orientation), Poor/None fail.
        if is_degrader:
            constraint = "Pass"
        elif tier == CFG.TIER_ORDER[4]:
            constraint = "Partial"
        else:
            constraint = "Fail"

        close_residues = [k for k, d in dists.items() if d <= CATALYTIC_DIST_CUTOFF]
        results.update({
            "catalytic_dist_A": round(min(dists.values()), 2), "degrader_tier": tier,
            "is_degrader": is_degrader, "residues_within_6A": ";".join(close_residues) if close_residues else "None",
            "constraint_check": constraint, "scientific_meaning": meaning
        })
        return results
    except Exception as e:
        results.update({"catalytic_dist_A": 999.0, "error": str(e), "degrader_tier": "Error"})
        return results

# -------------------------------------------------------------------------------
# Step 7.3: Parallelised Model Analysis & Ranking Framework
# -------------------------------------------------------------------------------
def analyse_model_task(args: Tuple[Path, Dict[str, int], Path, str, str]) -> Tuple[str, Dict, float]:
    """Global worker routine functionally decoupled for robust serialisation capabilities across CPU threads."""
    cif_p, sites, json_p, m_name, smiles_str = args
    geom = check_catalytic_geometry(cif_p, sites, smiles_str)
    try:
        conf_data = json.loads(json_p.read_text())
        conf_score = float(conf_data.get("confidence_score", 0.0))
    except Exception:
        conf_score = 0.0
    return m_name, geom, conf_score

def select_best_degrader_model(br_dir: Path, mapped_sites: Dict[str, int], smiles_str: str = "") -> Tuple[str, Dict]:
    """Iterates iteratively through all predicted Boltz models, strictly selecting the candidate demonstrating optimal geometry."""
    best_model, best_rank_score, best_meta = "model_0", -100.0, {"degrader_tier": CFG.TIER_DECOY}
    TIER_SCORES = CFG.TIER_SCORE
    json_files = sorted(list(br_dir.rglob("confidence_*.json")))
    tasks = []
    for f in json_files:
        match = re.search(r"(model_\d+)", f.stem)
        if match:
            m_name = match.group(1)
            cif_path = f.parent / f"{f.parent.name}_{m_name}.cif"
            if not cif_path.exists(): cif_path = f.parent / f"{m_name}.cif"
            if cif_path.exists(): tasks.append((cif_path, mapped_sites, f, m_name, smiles_str))
    
    if tasks:
        results = [analyse_model_task(t) for t in tasks]
        for m_name, geom, conf in results:
            total_rank = (TIER_SCORES.get(geom.get("degrader_tier", CFG.TIER_DECOY), 0) * 10) + conf
            if total_rank > best_rank_score:
                best_rank_score = total_rank
                best_model = m_name
                best_meta = geom.copy()
                best_meta["confidence_score"] = conf
    return best_model, best_meta


# ===============================================================================
# SECTION 8: METRICS & GPU MANAGEMENT
# ===============================================================================

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
        conf_files = list(br_dir.rglob(f"confidence_*{model_name}.json"))
        if conf_files:
            d = json.loads(conf_files[0].read_text())
            out.update({k: d.get(k) for k in ["iptm", "confidence_score", "ptm", "ligand_iptm", "protein_iptm"]})
        plddt_files = list(br_dir.rglob(f"plddt_*{model_name}.npz"))
        if plddt_files:
            with np.load(str(plddt_files[0])) as data:
                if "plddt" in data: out["mean_plddt"] = float(np.mean(data["plddt"]))
        pae_files = list(br_dir.rglob(f"pae_*{model_name}.npz"))
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

def select_gpu():
    """Smart GPU routing logic: securely identifies the physical hardware block maintaining the most available VRAM."""
    try:
        info = query_gpu_memory() 
        if not info: return 0
        free = {i: total - used for i, (used, total) in enumerate(info)}
        return max(free.items(), key=lambda x: x[1])[0]
    except subprocess.TimeoutExpired:
        return 0 

def cpu_usage_summary():
    """Returns general CPU congestion levels."""
    try: return psutil.cpu_percent(interval=0.2), psutil.cpu_count(logical=True)
    except Exception: return 0.0, 1


# ===============================================================================
# SECTION 9: JOB PROCESSING LOGIC (CORE WORKER)
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 9.1: Status Validation Operations
# -------------------------------------------------------------------------------
def check_job_status(job_dir: Path) -> bool:
    """Verifies the finality of a job by checking for a complete and uncorrupted summary JSON file."""
    if not job_dir.exists(): return False
    summary = next(job_dir.glob("*_summary.json"), None)
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
            with open(yaml_path, 'r') as yf:
                data = yaml.safe_load(yf)
            if 'sequences' in data and len(data['sequences']) > 1:
                current_smi = data['sequences'][1].get('ligand', {}).get('smiles', '')
                if current_smi != new_smiles:
                    data['sequences'][1]['ligand']['smiles'] = new_smiles
                    with open(yaml_path, 'w') as yf:
                        yaml.safe_dump(data, yf, sort_keys=False, default_flow_style=False)
        except Exception:
            pass 

    if job_dir.exists():
        for cif in job_dir.rglob("*.cif"):
            try:
                with open(cif, 'r') as f:
                    lines = f.readlines()
                
                if not any("UPDATED_SMILES" in line for line in lines[:15]):
                    for i, line in enumerate(lines):
                        if line.startswith("data_model") or line.startswith("data_"):
                            lines.insert(i + 1, f"# UPDATED_SMILES {new_smiles}\n")
                            lines.insert(i + 2, f"# LIGAND_ID {lig_id}\n")
                            break
                    with open(cif, 'w') as f:
                        f.writelines(lines)
            except Exception:
                pass

# -------------------------------------------------------------------------------
# Step 9.3: Core Worker Engine Infrastructure
# -------------------------------------------------------------------------------

# Fields that cannot be recomputed from structure files — preserved across force-regen.
_PRESERVED_FIELDS = ["elapsed_seconds", "completed_at"]

def extract_preserved_fields(runs_dir: Path, backup_path: Path) -> dict:
    """Scans all completed *_summary.json files and extracts non-recoverable fields
    (elapsed_seconds, completed_at) keyed by job_name (folder name).
    Writes a backup JSON to backup_path before returning."""
    preserved = {}
    json_files = sorted(runs_dir.rglob("*_summary.json"))
    total = len(json_files)
    _pf_done = [0]

    def _load_pf(jf):
        try:
            with open(jf) as f:
                data = json.load(f)
            if data.get("status") != "Success":
                return None
            job_name = jf.parent.name
            return (job_name, {k: data[k] for k in _PRESERVED_FIELDS if k in data})
        except Exception:
            return None
        finally:
            _pf_done[0] += 1
            n = _pf_done[0]
            if n % 1000 == 0 or n == total:
                _tty_write(f"\r   Scanning summaries for preserved fields: {n}/{total} ...\033[K")
                sys.stdout.flush()

    from concurrent.futures import ThreadPoolExecutor as _TPE
    with _TPE(max_workers=min(max(1, (os.cpu_count() or 4) - 2), total or 1)) as ex:
        for item in ex.map(_load_pf, json_files):
            if item is not None:
                preserved[item[0]] = item[1]
    if total > 0:
        _tty_write("\r\033[K")
        sys.stdout.flush()
    try:
        backup_path.write_text(json.dumps(preserved, indent=2))
    except Exception as e:
        console_info(f"  Warning: could not write preserved-fields backup: {e}")
    return preserved


def process_single_job(job: Dict, prod_dir: Path, diffusion_samples: int, prev_elapsed: float, aln_dir: Path,
                       control_cif: Optional[Path], control_map: Dict, control_seq: str, gpu_queue=None,
                       preserved_fields: Optional[Dict] = None,
                       r3u_cif: Optional[Path] = None, r3u_map: Optional[Dict] = None) -> Tuple[Optional[Dict[str, Any]], int]:
    data = DEFAULT_METRICS.copy()
    jid = job["job_index"] 
    name = f"{jid}_{job.get('job_protein', job['protein'])}_{job['ligand']}"
    job_dir = prod_dir / "4_Prediction_Jobs" / name
    job_dir.mkdir(parents=True, exist_ok=True)
    summary_path = job_dir / f"{name}_summary.json"
    
    best_complex_dir = job_dir / "Best_Complex"
    if best_complex_dir.exists():
        shutil.rmtree(best_complex_dir, ignore_errors=True)
    # Retry mkdir — shutil.rmtree uses ignore_errors so the dir may briefly
    # still exist on a slow fs; exist_ok=True absorbs any EEXIST race.
    for _attempt in range(3):
        try:
            best_complex_dir.mkdir(parents=True, exist_ok=True)
            break
        except OSError:
            import time as _t; _t.sleep(0.05)
    
    aln_path = aln_dir / f"{name}_alignment.txt"
    stats_csv_path = aln_dir / "Alignment_Stats.csv"

    cached_aln = get_cached_alignment_for_protein(job["protein"])
    aln_inject = {
        "identity_pct": 0.0, "align_score": 0.0, 
        "Mapped_to_Control_All": "NA", "Mapped_to_Control_Cat_Triad": "NA"
    }
    if cached_aln:
        aln_inject.update({
            "identity_pct": cached_aln.get("identity_pct", 0.0),
            "align_score": cached_aln.get("align_score", 0.0),
            "Mapped_to_Control_All": cached_aln.get("Mapped_to_Control_All", "NA"),
            "Mapped_to_Control_Cat_Triad": cached_aln.get("Mapped_to_Control_Cat_Triad", "NA")
        })

    is_done = check_job_status(job_dir)
    
    if is_done and summary_path.exists() and preserved_fields is None:
        try:
            with open(summary_path, 'r') as f: data.update(json.load(f))
        except Exception: pass

    data["protein"] = job["protein"]
    data["ligand"]  = job["ligand"]
    data["job_name"] = name
    data["job_index"] = jid
    
    br_dir = next((p for p in job_dir.iterdir() if p.is_dir() and p.name.startswith("boltz_results_")), None)
    prediction_exists = False
    if br_dir and (br_dir / "predictions").exists():
        # A prediction is only considered complete when at least one confidence_*.json
        # exists alongside the CIF — Boltz writes these only after a model fully finishes.
        # A directory with CIFs but no confidence JSONs is a killed/partial run.
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
        br_dir_check = next((p for p in job_dir.iterdir() if p.is_dir() and p.name.startswith("boltz_results_")), None)
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
        env["TRITON_CACHE_DIR"] = str(Path(tempfile.gettempdir()) / "triton_cache")
        # conda env PFAS manages CUDA/nvidia libraries via LD_LIBRARY_PATH
        env["OMP_NUM_THREADS"] = "4"
        env["MKL_NUM_THREADS"] = "4"
        env["OPENBLAS_NUM_THREADS"] = "4"
        
        cmd = [
            BOLTZ_BIN, "predict", str(job["yaml"]), "--out_dir", str(job_dir),
            "--cache", str(BOLTZ_CACHE), "--model", BOLTZ_MODEL,
            "--recycling_steps", str(RECYCLING_STEPS),
            "--diffusion_samples", str(diffusion_samples),
            "--accelerator", "gpu", "--devices", "1",
            "--use_msa_server", "--output_format", OUTPUT_FORMAT,
            "--no_kernels"
        ]
        
        success = False
        
        for attempt in range(RETRY_ON_FAIL + 1):
            try:
                with open(job_dir / f"{jid}_stdout.log", "w") as so, open(job_dir / f"{jid}_stderr.log", "w") as se:
                    proc = subprocess.run(cmd, stdout=so, stderr=se, env=env)
                    if proc.returncode == 0:
                        success = True
                        break
                    else:
                        console_info(f"    ! Scheduled Job {name} subsequently failed on allocated GPU[{gpu_idx}] (Attempt {attempt+1}/{RETRY_ON_FAIL+1}). Initiating retry sequence...")
                        time.sleep(RETRY_SLEEP)
            except Exception as e:
                console_info(f"    ! Routine Exception Encountered: {e}")
                time.sleep(RETRY_SLEEP)
    
        if gpu_queue is not None:
            gpu_queue.put(gpu_idx) 
        
        if not success: return None, gpu_idx
        data["elapsed_seconds"] = time.time() - start_time 

    br_dir = next((p for p in job_dir.iterdir() if p.is_dir() and p.name.startswith("boltz_results_")), None)
    
    # -------------------------------------------------------------------------------
    # Sub-Step 9.3.2: Execution Physics & Structural Topology Engine (CPU)
    # -------------------------------------------------------------------------------
    data.update(aln_inject)
    if run_gpu or repair_cpu:
        mapped, aln_stats, resname_map, map_all_str = map_active_site_residues(
            job["protein"], job["sequence"], out_aln_path=aln_path, stats_csv=stats_csv_path
        )
        _, map_triad_str = format_control_mappings(resname_map, map_all_str)

        data.update(aln_stats)
        data["active_site_mapping"] = str(mapped) 
        data["Mapped_to_Control_All"] = map_all_str
        data["Mapped_to_Control_Cat_Triad"] = map_triad_str
        data["alignment_reliable"] = data.get("identity_pct", 0) >= 25.0

        if br_dir:
            best_model, best_meta = select_best_degrader_model(br_dir, mapped, job.get("smiles", ""))
            data["best_model_name"] = best_model
            data.update(best_meta)
        else:
            cifs = list(best_complex_dir.glob("*.cif"))
            if not cifs: cifs = list(job_dir.glob("*.cif"))
            if cifs:
                best_model = "model_0"
                if "model_" in cifs[0].name:
                    try: best_model = re.search(r"(model_\d+)", cifs[0].name).group(1)
                    except Exception: pass
                data["best_model_name"] = best_model
                geom = check_catalytic_geometry(cifs[0], mapped, job.get("smiles", ""))
                data.update(geom)

        cif_path = None
        if br_dir:
            cif_path = br_dir / "predictions" / name / f"{name}_{data.get('best_model_name','model_0')}.cif"
            if not cif_path.exists():
                cands = list(br_dir.rglob(f"*{data.get('best_model_name')}.cif"))
                if cands: cif_path = cands[0]
        
        if not cif_path:
             cifs = list(best_complex_dir.glob("*.cif"))
             if not cifs: cifs = list(job_dir.glob("*.cif"))
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
            data["interaction_density"] = round(interaction_density, 4)
            data["interaction_density_calc"] = f"{num_int} (Ints) / {lig_heavy_atoms} (HeavyAtoms) = {interaction_density:.2f}"

            iptm_v = data.get("iptm") or 0.0
            plddt_v = data.get("mean_plddt") or 0.0
            c_pae = data.get("cross_interface_pae_mean") or 0.0
            conf = data.get("confidence_score") or 0.0

            score_v = (W_IPTM * float(iptm_v)) + (W_PLDDT * (plddt_v/100.0)) + \
                      (W_INTERACTIONS * min(1.0, interaction_density / 2.0)) - \
                      (W_CROSS_PAE * min(1.0, c_pae/50.0)) + (W_CONF * conf)
            
            score_v = max(min(score_v, 50), -50)
            binding_prob = 1.0 / (1.0 + math.exp(-score_v))
            data["binding_likelihood_computed"] = binding_prob
            data["binding_likelihood_calc"] = (
                f"({W_IPTM} * {float(iptm_v):.2f}) + ({W_PLDDT} * {plddt_v/100.0:.2f}) + "
                f"({W_INTERACTIONS} * {min(1.0, interaction_density/2.0):.2f}) - "
                f"({W_CROSS_PAE} * {min(1.0, c_pae/50.0):.2f}) + ({W_CONF} * {conf:.2f}) = {score_v:.2f} -> Sigmoid = {binding_prob:.4f}"
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
            tier = data.get("degrader_tier", CFG.TIER_DECOY)
            
            base_justification = generate_rich_justification(
                tier=tier, meaning=scientific_meaning, constraint=constraint,
                aligned_ok=data.get("alignment_reliable", False), identity=data.get("identity_pct", 0.0)
            )
            data["Justification"] = f"{base_justification} | Tier: {tier} | ID: {data.get('identity_pct', 0)}%"

            # -------------------------------------------------------------------------------
            # SUPERIMPOSITION & LIKELIHOOD SCORING ARCHITECTURE
            # -------------------------------------------------------------------------------
            data["hydrophobic_desolvation_ratio"] = 0.0   # default; overwritten below if control_cif available
            data["product_inhibition_penalty"]    = 0.0   # default; overwritten below if inhibition condition met
            if control_cif and control_cif.exists() and control_map:
                mech_res = analyse_candidate_structure(cif_path, control_cif, control_map)
                data.update(mech_res)
                
                ident = data.get("identity_pct", 0.0)
                data["Identity_to_Control"] = round(ident, 2)
                
                rmsd = data.get("Active_Site_RMSD", 99.0)
                if rmsd >= 99.0: geo_fit = 0.0
                else: geo_fit = 100.0 / (1.0 + rmsd)
                
                likelihood = (0.25 * ident) + (0.75 * geo_fit)
                
                desolvation_ratio = fh / max((fh + fp), 1)
                data["hydrophobic_desolvation_ratio"] = round(desolvation_ratio, 3)
                if desolvation_ratio > CFG.SCORE_DESOLVATION_THRESHOLD:
                    likelihood += (CFG.SCORE_DESOLVATION_BONUS_FACTOR * desolvation_ratio)

                if interaction_density > CFG.SCORE_INHIBITION_DENSITY_THRESHOLD:
                    inhibition_penalty = CFG.SCORE_INHIBITION_PENALTY_FACTOR * (interaction_density - CFG.SCORE_INHIBITION_DENSITY_THRESHOLD)
                    likelihood -= inhibition_penalty
                    data["product_inhibition_penalty"] = round(inhibition_penalty, 2)
                else:
                    data["product_inhibition_penalty"] = 0.0

                data["ActiveSite_Conservation_Score"] = round(max(0, min(100, likelihood)), 2)

            if r3u_cif and r3u_cif.exists() and r3u_map:
                _r3u = analyse_candidate_structure(cif_path, r3u_cif, r3u_map)
                _r3u_rmsd  = _r3u.get("Active_Site_RMSD", 99.0)
                _r3u_ident = float(data.get("identity_pct", 0.0))
                _r3u_geo   = 0.0 if _r3u_rmsd >= 99.0 else 100.0 / (1.0 + _r3u_rmsd)
                data["r3u_Active_Site_RMSD"]              = _r3u_rmsd
                data["r3u_Mechanistic_Fingerprint_Score"] = _r3u.get("Mechanistic_Fingerprint_Score", 0.0)
                data["r3u_Halide_Stabilisation"]          = _r3u.get("Halide_Stabilisation", False)
                data["r3u_Carboxylate_Clamp"]             = _r3u.get("Carboxylate_Clamp", False)
                data["r3u_ActiveSite_Conservation_Score"]     = round(max(0.0, 0.25 * _r3u_ident + 0.75 * _r3u_geo), 2)
            else:
                data["r3u_Active_Site_RMSD"]              = 99.0
                data["r3u_Mechanistic_Fingerprint_Score"] = 0.0
                data["r3u_Halide_Stabilisation"]          = False
                data["r3u_Carboxylate_Clamp"]             = False
                data["r3u_ActiveSite_Conservation_Score"]     = 0.0

            if preserved_fields:
                for _pf_k, _pf_v in preserved_fields.items():
                    data[_pf_k] = _pf_v
            summary_path.write_text(json.dumps(data, indent=2))
            
            best_model = data.get("best_model_name", "model_0")
            cif_src = None
            if br_dir:
                cif_src = br_dir / "predictions" / name / f"{name}_{best_model}.cif"
                if not cif_src.exists():
                    cands = list(br_dir.rglob(f"*{best_model}.cif"))
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
    if "identity_pct" not in data or data.get("identity_pct", 0) == 0:
        mapped, aln_stats, resname_map, map_all_str = map_active_site_residues(
             job["protein"], job["sequence"], out_aln_path=aln_path, stats_csv=stats_csv_path
        )
        _, map_triad_str = format_control_mappings(resname_map, map_all_str)
        data.update(aln_stats)
        data["active_site_mapping"] = str(mapped) 
        data["Mapped_to_Control_All"] = map_all_str
        data["Mapped_to_Control_Cat_Triad"] = map_triad_str

    return data, gpu_idx

def worker_task_wrapper(args):
    """Secure pickle-safe wrapper routine intended strictly for Python multiprocessing environments."""
    try: return process_single_job(*args)
    except Exception as e:
        job = args[0] if args else {}
        job_name = f"{job.get('job_index', '?')}_{job.get('job_protein', job.get('protein', '?'))}_{job.get('ligand', '?')}"
        return {"error": str(e), "traceback": traceback.format_exc(), "job_name": job_name, "status": "FAILED"}, -1

def _reanalyse_job_r3u(args):
    """Worker: return 3R3U structural comparison metrics for a job.
    Reads from summary.json cache if available (populated during Step 10.8);
    falls back to re-reading the CIF and running analyse_candidate_structure."""
    job_name, runs_dir_str, r3u_cif_str, r3u_control_map = args
    try:
        job_dir = Path(runs_dir_str) / job_name
        summary_path = job_dir / f"{job_name}_summary.json"
        if summary_path.exists():
            try:
                with open(summary_path) as _f:
                    _s = json.load(_f)
                if "r3u_Active_Site_RMSD" in _s:
                    return {
                        "job_name": job_name,
                        "Active_Site_RMSD": _s["r3u_Active_Site_RMSD"],
                        "Mechanistic_Fingerprint_Score": _s.get("r3u_Mechanistic_Fingerprint_Score", 0.0),
                        "Halide_Stabilisation": _s.get("r3u_Halide_Stabilisation", False),
                        "Carboxylate_Clamp": _s.get("r3u_Carboxylate_Clamp", False),
                    }
            except Exception as _e:
                import logging as _logging
                _logging.warning(f"[_reanalyse_job_r3u] Failed to parse cached summary for {job_name}: {_e}")
        candidate_cif = None
        bc_dir = job_dir / "Best_Complex"
        if bc_dir.exists():
            cifs = sorted(bc_dir.glob("*.cif"))
            if cifs:
                candidate_cif = cifs[0]
        if not candidate_cif:
            candidate_cif = next(job_dir.glob("**/*_model_0.cif"), None)
        if not candidate_cif or not candidate_cif.exists():
            return None
        result = analyse_candidate_structure(candidate_cif, Path(r3u_cif_str), r3u_control_map)
        return {
            "job_name": job_name,
            "Active_Site_RMSD": result.get("Active_Site_RMSD", 99.0),
            "Mechanistic_Fingerprint_Score": result.get("Mechanistic_Fingerprint_Score", 0.0),
            "Halide_Stabilisation": result.get("Halide_Stabilisation", False),
            "Carboxylate_Clamp": result.get("Carboxylate_Clamp", False),
        }
    except Exception:
        return None

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
    "degrader_tier", "is_degrader",
    "ActiveSite_Conservation_Score", "Mechanistic_Fingerprint_Score",
    # --- Confidence scores ---
    "confidence_score", "iptm", "ptm", "ligand_iptm", "protein_iptm",
    "mean_plddt", "cross_interface_pae_mean",
    # --- Binding scores ---
    "binding_likelihood_computed", "custom_affinity_score", "interaction_density",
    # --- SN2 geometry (+ auxiliary non-gating reference angles) ---
    "sn2_attack_angle", "sn2_trajectory_dev",
    "burgi_dunitz_angle", "flippin_lodge_offset",
    "catalytic_dist_A", "Active_Site_RMSD", "Identity_to_Control",
    # --- Key active site distances ---
    "dist_Nuc", "dist_Carb1", "dist_Carb2", "dist_Acid",
    "dist_Stab_H", "dist_Stab_W", "dist_Stab_Y", "dist_Base",
    # --- Catalytic residue distances ---
    "dist_ASP110", "dist_HIS277", "dist_ASP134",
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
    "identity_pct", "align_score", "seq_length", "gap_count",
    "alignment_reliable", "target_sequence",
    # --- Active site mapping ---
    "active_site_mapping", "Mapped_to_Control_All", "Mapped_to_Control_Cat_Triad",
    "best_model_name",
    # --- Scoring / verdict ---
    "constraint_check", "scientific_meaning", "Justification",
    "sn2_alignment_score", "mechanistic_score", "residues_within_6A",
    # --- Verbose calculation strings ---
    "interaction_density_calc", "binding_likelihood_calc", "custom_affinity_calc",
    # --- Clash / penalty metrics ---
    "teflon_shield_clashes", "mainchain_clash_ratio",
    "hydrophobic_desolvation_ratio", "product_inhibition_penalty",
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

    job_dirs = [d for d in runs_dir.iterdir() if d.is_dir()] if runs_dir.exists() else []

    # Collect (src, dst) pairs in parallel — each dir check is independent I/O
    to_copy: list = []
    total_count = [0]
    _lock = threading.Lock()

    def _collect(job_dir):
        bc_dir = job_dir / "Best_Complex"
        if not bc_dir.exists():
            return
        local_copy = []
        local_total = 0
        for cif in bc_dir.glob("*.cif"):
            local_total += 1
            dst = dest_dir / cif.name
            if not (dst.exists() and dst.stat().st_size == cif.stat().st_size):
                local_copy.append((cif, dst))
        with _lock:
            total_count[0] += local_total
            to_copy.extend(local_copy)

    with ThreadPoolExecutor(max_workers=min(max(1, (os.cpu_count() or 4) - 2), len(job_dirs) or 1)) as ex:
        list(ex.map(_collect, job_dirs))
    total = total_count[0]

    def _do_copy(pair):
        src, dst = pair
        try:
            shutil.copy2(src, dst)
        except Exception:
            pass

    if to_copy:
        n_workers = min(max(1, (os.cpu_count() or 4) - 2), len(to_copy))
        with ThreadPoolExecutor(max_workers=n_workers) as ex:
            list(ex.map(_do_copy, to_copy))

    return total


def append_rows_to_csv(rows: list, csv_path: Path):
    """Appends a list of result rows to the master CSV, deduplicates, sorts by job_name,
    assigns a plain sequential job_index (0-based row counter), and reorders columns
    into the canonical human-readable CSV_COLUMN_ORDER."""
    try:
        df_new = pd.DataFrame(rows)
        if csv_path.exists() and csv_path.stat().st_size > 0:
            try:
                df_exist = pd.read_csv(csv_path, low_memory=False, on_bad_lines='warn')
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
        df_new.to_csv(csv_path, index=False)
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
                is_fallback = job.pop('fallback_individual', False)
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
                    local_rows.append(flatten_job_result(res))
                    if len(local_rows) >= 10:
                        append_rows_to_csv(local_rows, csv_path)
                        local_rows.clear()
        except Exception as e:
            time.sleep(5)
    if local_rows:
        append_rows_to_csv(local_rows, csv_path)


def rebuild_csv_from_summaries(runs_dir: Path, csv_path: Path,
                               extract_preserved: bool = False,
                               preserved_backup_path: Optional[Path] = None
                               ) -> tuple:
    """Single parallel scan over all *_summary.json files.

    Returns (n_rows_written, preserved_fields_map).
    When extract_preserved=False, preserved_fields_map is always {}.
    Pass extract_preserved=True (resume mode) to populate it in the same
    disk pass — avoids a second full rglob over 58k files.

    IMPORTANT — folder name is the canonical truth for job_index and job_name.
    """
    json_files = sorted(runs_dir.rglob("*_summary.json"))
    total = len(json_files)
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
            row = flatten_job_result(data)
            pf  = {k: data[k] for k in _PRESERVED_FIELDS if k in data} if extract_preserved else None
            return row, (folder_name, pf) if pf is not None else (folder_name, None)
        except Exception:
            return None, None
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

    rows = [row for row, _ in results if row is not None]
    preserved: dict = {}
    if extract_preserved:
        for _, pf_item in results:
            if pf_item is not None:
                job_name, pf_data = pf_item
                if pf_data:
                    preserved[job_name] = pf_data

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
        df.to_csv(csv_path, index=False)

    if extract_preserved and preserved_backup_path:
        try:
            preserved_backup_path.write_text(json.dumps(preserved, indent=2))
        except Exception as e:
            console_info(f"  Warning: could not write preserved-fields backup: {e}")

    return len(rows), preserved

# ===============================================================================
# SECTION 10: MAIN SYSTEM ENTRY POINT
# ===============================================================================

def main():
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

    # Step 1: migrate very-old flat folder names (legacy runs only)
    _folder_renames = {
        "1_Input_Data":     "1_Input_FASTA_and_SMILES",
        "2_YAML_Files":     "2_Boltz2_YAML_Configs",
        "3_Colabfold":      "3_MSA_Sequence_Data",   # kept so truly old runs cascade forward
        "4_Alignment_Data": "4_Active_Site_Alignments",
        "5_RUNS":           "4_Prediction_Jobs",
        "5_Prediction_Jobs": "4_Prediction_Jobs",
    }
    for _old, _new in _folder_renames.items():
        _old_path = PROD / _old
        _new_path = PROD / _new
        if _old_path.exists() and not _new_path.exists():
            _old_path.rename(_new_path)

    D_IN   = PROD / "1_Input_FASTA_and_SMILES"
    D_YAML = PROD / "2_Boltz2_YAML_Configs"
    D_RUNS = PROD / "4_Prediction_Jobs"

    # Step 2: migrate legacy folder names → current sequential naming scheme
    D_SEQ       = PROD / "3_Sequence_Reference_Data"
    D_COLABFOLD = D_SEQ / "MSA_Sequences"
    D_ALN       = D_SEQ / "Active_Site_Alignments"
    D_SEQ.mkdir(parents=True, exist_ok=True)

    _old_msa  = PROD / "3_MSA_Sequence_Data"
    _old_aln  = PROD / "4_Active_Site_Alignments"
    _old_runs = PROD / "5_Prediction_Jobs"
    if _old_msa.exists() and not D_COLABFOLD.exists():
        shutil.move(str(_old_msa), str(D_COLABFOLD))
        console_info("  [Migrate] 3_MSA_Sequence_Data  →  3_Sequence_Reference_Data/MSA_Sequences")
    if _old_aln.exists() and not D_ALN.exists():
        shutil.move(str(_old_aln), str(D_ALN))
        console_info("  [Migrate] 4_Active_Site_Alignments  →  3_Sequence_Reference_Data/Active_Site_Alignments")
    if _old_runs.exists() and not D_RUNS.exists():
        shutil.move(str(_old_runs), str(D_RUNS))
        console_info("  [Migrate] 5_Prediction_Jobs  →  4_Prediction_Jobs")

    for d in [D_IN, D_YAML, D_RUNS, D_COLABFOLD, D_ALN]: d.mkdir(parents=True, exist_ok=True)
    
    initial_folders_count = len([d for d in D_RUNS.iterdir() if d.is_dir()]) if D_RUNS.exists() else 0
    removed_folders_count = 0
    
    ts_now = datetime.now().strftime("%Y%m%d_%H%M%S")
    LOG_PATH = PROD / f"5_Boltz2_Log_{ts_now}.log"
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
            f_path = next(D_IN.glob("*.fasta"))
            s_path = next(D_IN.glob("*.smi"))
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
    
    files_on_disk = len(list(D_ALN.glob("*.txt")))

    proteins = []
    _MERGE_PREFIX = re.compile(r'^(\d+)_(.+?)_*$')  
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
    console_info(f"  SMILES Validation & Standardisation")
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
    console_info(f"  System Compute Resources")
    console_info(SEPARATOR_LIGHT)
    console_info(f"  ┌{'─'*(_ck+2)}┬{'─'*(_cv+4)}┐")
    for _k, _v in _cr_rows:
        console_info(f"  │  {_k:<{_ck}}│  {str(_v):>{_cv}}  │")
    console_info(f"  └{'─'*(_ck+2)}┴{'─'*(_cv+4)}┘")
    
    # -------------------------------------------------------------------------------
    # Step 10.3: System Reference Data Initialisation
    # -------------------------------------------------------------------------------
    ref_pdb_path = setup_reference_data(PROD)
    
    if resumed:
        sanitize_dataset(PROD)

    # Build an ordered expected mapping of sequence+ligand keys to current job stems.
    expected_job_stems = {}
    num_ligands = len(ligands)
    for protein_order, pid, seq, fasta_pos in proteins:
        protein_position = fasta_pos if fasta_pos > 0 else protein_order
        for lig_order, lid, smi in ligands:
            jid = str((protein_position - 1) * num_ligands + lig_order).zfill(7)
            job_pid = f"{protein_position}_{pid}"
            expected_job_stems[(sequence_hash(seq), smi)] = f"{jid}_{job_pid}_{lid}"

    if resumed:
        reconcile_resume_job_names(PROD, expected_job_stems)

    load_cached_alignments(D_ALN / "Alignment_Stats.csv")
    console_separator()

    # -------------------------------------------------------------------------------
    # Step 10.4: Construct Initial Task Matrix & Pre-Processing Evaluation
    # -------------------------------------------------------------------------------
    tasks = []
    task_count = 0
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
            cif_path = job_path / "boltz_results" / job_stem / f"{job_stem}_model_0.cif"
            is_completed = job_path.exists() and cif_path.exists()

            if y_name in existing_yamls:
                try:
                    with open(y_path, 'r') as yf:
                        y_data = yaml.safe_load(yf)
                    existing_seq = y_data.get('sequences', [])[0].get('protein', {}).get('sequence', '')
                    existing_smi = y_data.get('sequences', [])[1].get('ligand', {}).get('smiles', '')
                    
                    if existing_seq == seq and existing_smi == smi:
                        yaml_needs_write = False
                    else:
                        if is_completed:
                            console_info(f"[Resume] Structural metadata mismatch in {y_name} bypassed to preserve workflow continuity.")
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

                with open(y_path, 'w') as f: yaml.safe_dump(payload, f, sort_keys=False, default_flow_style=False)
                GLOBAL_STATS["created_yaml"] += 1

            tasks.append({
                "job_index": jid, "protein": pid, "job_protein": job_pid, "protein_idx": protein_order,
                "ligand": lid, "sequence": seq, "smiles": smi, "yaml": y_path
            })

    if resumed:
        for _, pid, seq, fasta_pos in proteins:
            a3m_file = colabfold_a3m_path(PROD, pid)
            meta_file = colabfold_meta_path(PROD, pid)
            if a3m_file.exists():
                if meta_file.exists():
                    try:
                        with open(meta_file, 'r') as f: meta = json.load(f)
                        if meta.get('sequence_sha256') != sequence_hash(seq):
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
        # Migrate very-old legacy DeHa4 folder names → _2_ intermediary (then caught below)
        _ctrl_legacy    = D_RUNS / "0000000_DeHa4_Control_26_Fluoroacetate"
        _ctrl_v0        = D_RUNS / "0000000_01_DeHa4_Control_26_Fluoroacetate"
        _ctrl_v2        = D_RUNS / "0000000_02_DeHa4_Control_26_Fluoroacetate"
        if _ctrl_legacy.exists() and not _ctrl_v2.exists():
            deep_rename_job_folder(_ctrl_legacy, "0000000_DeHa4_Control_26_Fluoroacetate",
                                   "0000000_02_DeHa4_Control_26_Fluoroacetate")
            _ctrl_legacy.rename(_ctrl_v2)
            console_info("  [Migrate] 0000000_DeHa4_Control_26_Fluoroacetate  →  0000000_02_DeHa4_Control_26_Fluoroacetate")
        elif _ctrl_v0.exists() and not _ctrl_v2.exists():
            deep_rename_job_folder(_ctrl_v0, "0000000_01_DeHa4_Control_26_Fluoroacetate",
                                   "0000000_02_DeHa4_Control_26_Fluoroacetate")
            _ctrl_v0.rename(_ctrl_v2)
            console_info("  [Migrate] 0000000_01_DeHa4_Control_26_Fluoroacetate  →  0000000_02_DeHa4_Control_26_Fluoroacetate")
        # Migrate very-old 3R3U crystal folder → _1_ intermediary (then caught below)
        _r3u_legacy     = D_RUNS / "0000001_0_3R3U_Crystal_26_Fluoroacetate"
        _r3u_v1         = D_RUNS / "0000000_02_3R3U_Control_26_Fluoroacetate"
        if _r3u_legacy.exists() and not _r3u_v1.exists():
            deep_rename_job_folder(_r3u_legacy, "0000001_0_3R3U_Crystal_26_Fluoroacetate",
                                   "0000000_02_3R3U_Control_26_Fluoroacetate")
            _r3u_legacy.rename(_r3u_v1)
            console_info("  [Migrate] 0000001_0_3R3U_Crystal_26_Fluoroacetate  →  0000000_02_3R3U_Control_26_Fluoroacetate")
        # Migrate old YAML files to _2_ / _1_ intermediary
        for _old_yaml, _new_yaml in [
            ("0000000_DeHa4_Control_26_Fluoroacetate.yaml",   "0000000_02_DeHa4_Control_26_Fluoroacetate.yaml"),
            ("0000000_01_DeHa4_Control_26_Fluoroacetate.yaml", "0000000_02_DeHa4_Control_26_Fluoroacetate.yaml"),
            ("0000001_0_3R3U_Crystal_26_Fluoroacetate.yaml",  "0000000_02_3R3U_Control_26_Fluoroacetate.yaml"),
        ]:
            _op = D_YAML / _old_yaml
            _np = D_YAML / _new_yaml
            if _op.exists() and not _np.exists():
                _op.rename(_np)
        # Migrate fixed-prefix naming (_2_DeHa4_ / _1_3R3U_) → sequential numbering (1–6)
        _n_ctrl_ligs = len(CTRL_LIGANDS)
        for _sq, (_mig_lig, _) in enumerate(CTRL_LIGANDS, start=1):
            _old_d = D_RUNS / f"0000000_02_DeHa4_Control_{_mig_lig}"
            _new_d = D_RUNS / f"0000000_{_sq}_DeHa4_Control_{_mig_lig}"
            if _old_d.exists() and not _new_d.exists():
                deep_rename_job_folder(_old_d, f"0000000_02_DeHa4_Control_{_mig_lig}", f"0000000_{_sq}_DeHa4_Control_{_mig_lig}")
                _old_d.rename(_new_d)
                console_info(f"  [Migrate] 0000000_02_DeHa4_Control_{_mig_lig}  →  0000000_{_sq}_DeHa4_Control_{_mig_lig}")
            _old_dy = D_YAML / f"0000000_02_DeHa4_Control_{_mig_lig}.yaml"
            _new_dy = D_YAML / f"0000000_{_sq}_DeHa4_Control_{_mig_lig}.yaml"
            if _old_dy.exists() and not _new_dy.exists():
                _old_dy.rename(_new_dy)
        for _sq, (_mig_lig, _) in enumerate(CTRL_LIGANDS, start=_n_ctrl_ligs + 1):
            _old_r = D_RUNS / f"0000000_02_3R3U_Control_{_mig_lig}"
            _new_r = D_RUNS / f"0000000_{_sq}_3R3U_Control_{_mig_lig}"
            if _old_r.exists() and not _new_r.exists():
                deep_rename_job_folder(_old_r, f"0000000_02_3R3U_Control_{_mig_lig}", f"0000000_{_sq}_3R3U_Control_{_mig_lig}")
                _old_r.rename(_new_r)
                console_info(f"  [Migrate] 0000000_02_3R3U_Control_{_mig_lig}  →  0000000_{_sq}_3R3U_Control_{_mig_lig}")
            _old_ry = D_YAML / f"0000000_02_3R3U_Control_{_mig_lig}.yaml"
            _new_ry = D_YAML / f"0000000_{_sq}_3R3U_Control_{_mig_lig}.yaml"
            if _old_ry.exists() and not _new_ry.exists():
                _old_ry.rename(_new_ry)
        for _sq, (_ctl_lig, _) in enumerate(CTRL_LIGANDS, start=1):
            active_job_stems.add(f"0000000_{_sq}_DeHa4_Control_{_ctl_lig}")
            active_job_stems.add(f"0000000_{_sq + _n_ctrl_ligs}_3R3U_Control_{_ctl_lig}")
        removed_folders_count = purge_orphans(PROD, active_job_stems)

        # Wipe per-job analysis outputs so re-analysis runs fresh on resume.
        # GPU predictions (boltz_results_*/) are preserved — only analysis
        # artifacts are removed so changed parameters take full effect.
        console_info("Wiping per-job analysis outputs for full re-analysis...")
        _wipe_dirs = [d for d in D_RUNS.iterdir() if d.is_dir() and d.name[0].isdigit()]
        _wiped = [0]

        def _wipe_one(d):
            for f in d.glob("*_summary.json"):
                try: f.unlink()
                except Exception: pass
            for f in d.glob("*_interactions.csv"):
                try: f.unlink()
                except Exception: pass
            bc = d / "Best_Complex"
            if bc.exists():
                try: shutil.rmtree(bc)
                except Exception: pass
            _wiped[0] += 1

        from concurrent.futures import ThreadPoolExecutor as _TPE
        with _TPE(max_workers=min(max(1, (os.cpu_count() or 4) - 2), len(_wipe_dirs) or 1)) as _wex:
            list(_wex.map(_wipe_one, _wipe_dirs))
        console_info(f" -> Wiped analysis outputs for {_wiped[0]:,} job directories.")

        # Delete stale master CSVs — will be rebuilt from fresh analysis.
        for _old_csv in PROD.glob("*_Master*.csv"):
            try: _old_csv.unlink()
            except Exception: pass

    console_separator()
    console_info("")

    # -------------------------------------------------------------------------------
    # Step 10.5: Control Job Execution & Calibration
    # DEHA4 (input enzyme) vs 3R3U (crystal reference) × 3 fluoroacetate control substrates.
    # Fluoroacetate map drives all 58k pipeline comparisons.
    # -------------------------------------------------------------------------------
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
    ctrl_result_data                       = None
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
            with open(_ctrl_yaml, 'w') as f:
                yaml.safe_dump(payload, f)
        _ctrl_idx[0] += 1
        _lig_short_name = _lig_name.split("_", 1)[-1]
        console_info(f"  [{_ctrl_idx[0]}/{_n_ctrl_total}]  DeHa4 × {_lig_short_name:<18}  →  {_ctrl_name}")
        _ctrl_res, _ = process_single_job(_ctrl_task, PROD, 5, 0.0, D_ALN, None, {}, "", None)
        ctrl_results[_lig_name] = _ctrl_res or {}
        _cifs = list((D_RUNS / _ctrl_name).glob("**/*.cif"))
        if _cifs:
            _ctrl_summary_p = D_RUNS / _ctrl_name / f"{_ctrl_name}_summary.json"
            _cached_ts = None
            _cached_cm = None
            if _ctrl_summary_p.exists():
                try:
                    _cs = json.load(open(_ctrl_summary_p))
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
                    "rama_ref": None, "rama_con": None, "png_names": [],
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
            ctrl_result_data = _ctrl_res
            if not _cifs:
                console_info("CRITICAL OPERATIONAL ERROR: Failed to generate DeHa4 Control (Fluoroacetate) CIF.")
                sys.exit(1)
            control_cif = _cifs[0]
            control_map = ctrl_maps.get("26_Fluoroacetate", {})

    # -------------------------------------------------------------------------------
    # Step 10.5a-post: Retroactive MFP computation for DeHa4 control jobs.
    #
    # DeHa4 is processed without a control_cif/control_map (they do not yet
    # exist when DeHa4 runs), so Mechanistic_Fingerprint_Score stays at 0.0.
    # Now that ref_pdb_path (3R3U crystal) is available, MFP is retroactively
    # computed by comparing each DeHa4 CIF against the crystal structure
    # reference positions. This is the same comparison applied to all
    # candidate proteins in the main pipeline (DeHa4 structure vs 3R3U crystal).
    #
    # crystal_map: residue numbers in the 3R3U crystal PDB (pdb_id from
    # REF_ACTIVE_SITE_MAP), used as anchor positions inside analyse_candidate_structure.
    # -------------------------------------------------------------------------------
    _crystal_map = {k: v["pdb_id"] for k, v in REF_ACTIVE_SITE_MAP.items()}
    for _lig_name, _lig_smiles in CTRL_LIGANDS:
        _ctrl_cif = ctrl_cifs.get(_lig_name)
        if _ctrl_cif and _ctrl_cif.exists() and ref_pdb_path and ref_pdb_path.exists():
            try:
                _deha4_mfp_res = analyse_candidate_structure(_ctrl_cif, ref_pdb_path, _crystal_map)
                if _lig_name in ctrl_results and ctrl_results[_lig_name]:
                    ctrl_results[_lig_name]["Mechanistic_Fingerprint_Score"] = _deha4_mfp_res.get(
                        "Mechanistic_Fingerprint_Score", 0.0
                    )
                    ctrl_results[_lig_name]["Halide_Stabilisation"] = _deha4_mfp_res.get(
                        "Halide_Stabilisation", False
                    )
                    ctrl_results[_lig_name]["Carboxylate_Clamp"] = _deha4_mfp_res.get(
                        "Carboxylate_Clamp", False
                    )
                    _rmsd = _deha4_mfp_res.get("Active_Site_RMSD", 99.0)
                    ctrl_results[_lig_name]["Active_Site_RMSD"] = _rmsd
                    _geo = 0.0 if _rmsd >= 99.0 else 100.0 / (1.0 + _rmsd)
                    ctrl_results[_lig_name]["ActiveSite_Conservation_Score"] = round(25.0 + 0.75 * _geo, 2)
            except Exception as _mfp_err:
                console_info(f"  [Warning] DeHa4 MFP retroactive computation failed for {_lig_name}: {_mfp_err}")

    # -------------------------------------------------------------------------------
    # Step 10.5b: 3R3U Crystal Reference Jobs (Experimental Positive Control)
    # Run 3R3U against all three fluoroacetate control substrates.
    # Per-ligand calibration + Ramachandran plots for each.
    # Fluoroacetate (26) provides r3u_boltz_cif / r3u_control_map for 58 k pipeline.
    # -------------------------------------------------------------------------------
    r3u_result_data:   Optional[dict]      = None
    r3u_csv_path                           = None
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
                with open(_r3u_yaml, 'w') as _f:
                    yaml.safe_dump(r3u_payload, _f, sort_keys=False, default_flow_style=False)
            _ctrl_idx[0] += 1
            _lig_short_name = _lig_name.split("_", 1)[-1]
            console_info(f"  [{_ctrl_idx[0]}/{_n_ctrl_total}]  3R3U  × {_lig_short_name:<18}  →  {_r3u_name}")
            _r3u_res, _ = process_single_job(
                _r3u_task, PROD, 5, 0.0, D_ALN, control_cif, control_map, DEHA4_CONTROL_SEQ, None, None
            )
            r3u_results[_lig_name] = _r3u_res or {}
            _r3u_cifs = list((D_RUNS / _r3u_name).glob("**/*.cif"))
            if _r3u_cifs:
                _r3u_summary_p = D_RUNS / _r3u_name / f"{_r3u_name}_summary.json"
                _r3u_cached_ts = None
                _r3u_cached_cm = None
                if _r3u_summary_p.exists():
                    try:
                        _rs = json.load(open(_r3u_summary_p))
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
                        "rama_ref": None, "rama_con": None, "png_names": [],
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

        # Print consolidated calibration summary (Ramachandran + mapping tables)
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
            _row("Degrader Tier",              "degrader_tier",                 fmt="{}",     unit="")
            _row("SN2 Attack Angle",           "sn2_attack_angle",              fmt="{:.1f}", unit="°")
            _row("Nucleophile Distance (ASP)", "dist_Nuc",                      fmt="{:.2f}", unit=" Å")
            _row("Base Distance (HIS)",        "dist_Base",                     fmt="{:.2f}", unit=" Å")
            _row("Acid Distance (ASP)",        "dist_Acid",                     fmt="{:.2f}", unit=" Å")
            _row("Mechanistic Fingerprint",    "Mechanistic_Fingerprint_Score", fmt="{:.2f}", unit="")
            _row("Active Site Conservation Score",  "ActiveSite_Conservation_Score",     fmt="{:.2f}", unit="%")
            _row("Boltz-2 Confidence",         "confidence_score",              fmt="{:.4f}", unit="")
            console_info(f"  {'─'*_tot}")

        # --- Ligand selectivity ranking ---
        # Rank by: (1) Degrader Tier, (2) Mechanistic Fingerprint Score, (3) SN2 angle proximity to 180°.
        # Tier hierarchy dynamically loaded.
        # SN2 reactions require a nucleophilic attack angle close to 180° (back-side attack).
        _TIER_RANK = CFG.TIER_RANK
        console_info(f"\n  {'─'*_tot}")
        console_info(f"  {'Substrate Selectivity Ranking':^{_tot}}")
        console_info(f"  {'Rank criteria: (1) Degrader Tier  (2) Mechanistic Fingerprint  (3) SN2 angle → 180°':^{_tot}}")
        console_info(f"  {'─'*_tot}")
        for _prot_s, _res_dict, _ts_dict in [
            ("3R3U  (Rhodopseudomonas palustris FAcD — crystal sequence)",  r3u_results,  r3u_trust_scores),
            ("DeHa4 (Input enzyme sequence — primary screening target)",    ctrl_results, ctrl_trust_scores),
        ]:
            ranked = sorted(
                [(ln, _res_dict.get(ln, {})) for ln, _ in CTRL_LIGANDS],
                key=lambda x: (
                    _TIER_RANK.get(x[1].get("degrader_tier", CFG.TIER_DECOY), 0),
                    float(x[1].get("Mechanistic_Fingerprint_Score") or 0.0),
                    -abs(float(x[1].get("sn2_attack_angle") or 0.0) - 180.0),
                ),
                reverse=True
            )
            console_info(f"\n  {_prot_s}:")
            for _rank, (ln, d) in enumerate(ranked, 1):
                mfp   = float(d.get("Mechanistic_Fingerprint_Score") or 0.0)
                sn2   = float(d.get("sn2_attack_angle") or 0.0)
                _tier = d.get("degrader_tier", CFG.TIER_DECOY)
                ts    = _ts_dict.get(ln)
                ts_s  = f"Trust: {ts:.3f} Å" if ts is not None else "Trust: —"
                _lig_short = ln.split("_", 1)[-1] if "_" in ln else ln
                if _rank == 1:
                    if _tier not in (CFG.TIER_DECOY, "Error", CFG.TIER_POOR, ""):
                        tag = f"  ← Preferred substrate"
                    else:
                        tag = f"  ← Top-ranked by MFP/SN2  [Tier: {_tier} — no degradation activity confirmed]"
                else:
                    tag = ""
                console_info(
                    f"    Rank {_rank}  │  {_lig_short:<20}  │  Tier: {_tier:<9}  │  "
                    f"MFP: {mfp:.2f}  │  SN2: {sn2:.1f}°  │  {ts_s}{tag}"
                )
        console_info(f"  {'─'*_tot}")
        console_info(f"  {'MFP=0.00: no mechanistic fingerprint match.  MFP=0.75: partial match (3/4 checks).  MFP=1.00: full match.':^{_tot}}")
        console_info(f"  {'─'*_tot}\n")

        r3u_csv_path = None  # reference data shown in tables above; no separate file written
    else:
        console_info("  [Warning] Could not extract sequence from 3R3U PDB — skipping reference jobs.")
    console_separator()

    # -------------------------------------------------------------------------------
    # Step 10.6: System Workload Audit
    # -------------------------------------------------------------------------------
    # A "completed" job has both: (a) a Boltz-2 CIF prediction on disk AND
    # (b) a fully written _summary.json with status="Success".
    # On --resume with wipe, summaries are cleared so completed_jobs = 0 even
    # if GPU predictions exist — those jobs re-run CPU analysis only (not GPU).
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

    # Count how many pending proteins already have a valid MSA cached
    console_info(f" -> {completed_jobs:,} complete, {pending_jobs:,} pending. Checking MSA cache...")
    msa_ready   = sum(1 for pid in pending_proteins if validate_a3m_file(D_COLABFOLD / f"{pid}.a3m"))
    msa_missing = len(pending_proteins) - msa_ready

    console_info("")
    _ww = 62
    console_info(f"  {'─'*_ww}")
    console_info(f"  {'SYSTEM WORKLOAD SUMMARY':^{_ww}}")
    console_info(f"  {'─'*_ww}")
    console_info(f"  {'Total jobs scheduled':<44}  {len(tasks):>8,}")
    console_info(f"  {'  Completed (Boltz-2 prediction + CPU analysis)':<44}  {completed_jobs:>8,}")
    if resumed:
        # On resume: GPU predictions on disk but summary JSONs wiped for re-analysis.
        # Distinguish jobs with existing prediction dirs from those needing full GPU run.
        _dirs_with_pred = max(0, min(initial_folders_count - 6, pending_jobs))  # subtract 6 control dirs
        _need_gpu       = max(0, pending_jobs - _dirs_with_pred)
        console_info(f"  {'  Pending — GPU prediction on disk, CPU re-analysis queued':<44}  {_dirs_with_pred:>8,}")
        console_info(f"  {'  Pending — full GPU + CPU pipeline required':<44}  {_need_gpu:>8,}")
        console_info(f"  {'  (Resume mode: analysis outputs wiped for re-scoring)':<44}")
    else:
        console_info(f"  {'  Pending — awaiting GPU prediction + CPU analysis':<44}  {pending_jobs:>8,}")
    console_info(f"  {'─'*_ww}")
    console_info(f"  {'Unique proteins requiring processing':<44}  {len(pending_proteins):>8,}")
    console_info(f"  {'MSA (ColabFold) ready':<44}  {msa_ready:>8,}")
    console_info(f"  {'MSA still required (will download during run)':<44}  {msa_missing:>8,}")
    console_info(f"  {'Orphaned job folders purged at startup':<44}  {removed_folders_count:>8,}")
    console_info(f"  {'─'*_ww}")
    console_info("")

    if completed_jobs > 0:
        _tty_write(" -> Actively standardizing 'Best_Complex' allocations for pre-existing structural results...")
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
        _bc_dirs = list(D_RUNS.iterdir())
        with _TPE(max_workers=min(max(1, (os.cpu_count() or 4) - 2), len(_bc_dirs) or 1)) as _ex:
            list(_ex.map(_standardize_bc, _bc_dirs))
        _tty_write("\r -> Complete standardization of pre-existing designated folder structures finished successfully.\033[K\n")
        sys.stdout.flush()


    # -------------------------------------------------------------------------------
    # Step 10.7: Master CSV Rebuild + GPU Hardware Batch Prediction Phase
    # -------------------------------------------------------------------------------

    # --- Upfront bulk CSV rebuild from all completed summary JSONs ---
    console_separator()
    console_info("Rebuilding Master CSV, summary JSONs and Best Complex CIFs from all completed jobs...")
    console_info("  (This regenerates the CSV from per-job summary.json files and mirrors Best_Complex CIFs.)")
    _rebuild_t0 = time.time()
    _pf_backup = PROD / "_preserved_fields_backup.json"
    preserved_fields_map: dict = {}
    _rebuilt_rows = 0
    if completed_jobs == 0 and not resumed:
        # No summaries exist (fresh run) — skip rglob traversal.
        console_info(" -> No completed jobs — CSV will be populated as analysis completes.")
    else:
        _rebuilt_rows, preserved_fields_map = rebuild_csv_from_summaries(
            D_RUNS, CSV_PATH,
            extract_preserved=bool(resumed),
            preserved_backup_path=_pf_backup if resumed else None,
        )
        _rebuild_elapsed = time.time() - _rebuild_t0
        console_info(f" -> Master CSV rebuilt: {_rebuilt_rows:,} rows in {_rebuild_elapsed:.1f}s  →  {CSV_PATH.name}")

    # Rebuild the run-level best-complex mirror from all completed job summaries.
    # This ensures 0\2_Best_Complexes_CIFs/<tier>/ is always consistent on every
    # resume, even if previous runs were interrupted mid-copy or tiers changed.
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

    if resumed and preserved_fields_map:
        console_info(f"  Preserved fields extracted for {len(preserved_fields_map):,} jobs  →  {_pf_backup.name}")
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
        for _legacy_ws in ["8_System_Workspace", "6_System_Workspace"]:
            _legacy_path = PROD / _legacy_ws
            if _legacy_path.exists() and not workspace_dir.exists():
                _legacy_path.rename(workspace_dir)
                break
        workspace_dir.mkdir(exist_ok=True)

        batch_out_dir = workspace_dir / "Boltz_Batch_Output"
        batch_out_dir.mkdir(exist_ok=True)

        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = "0"
        env["CC"] = shutil.which("gcc") or "gcc"
        env["CXX"] = shutil.which("g++") or "g++"
        env["TRITON_CACHE_DIR"] = str(Path(tempfile.gettempdir()) / "triton_cache")
        # conda env PFAS manages CUDA/nvidia libraries via LD_LIBRARY_PATH
        env["OMP_NUM_THREADS"] = "4"
        env["MKL_NUM_THREADS"] = "4"
        env["OPENBLAS_NUM_THREADS"] = "4"
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
            prediction = next(job_dir.glob("**/predictions/**/*.cif"), None)
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
                # Mixed mode: GPU-only jobs also present — route pre-predicted jobs
                # through the background single-worker queue so analysis runs
                # concurrently with ongoing GPU prediction.
                for _pj in _pre_pred_jobs:
                    gpu_queue_for_analysis.put(_pj)
                console_info(f" -> {pre_predicted_count:,} jobs have pre-existing GPU predictions — re-routing to Analysis Queue...")
            else:
                # All jobs are pre-existing — skip the single-worker queue entirely.
                # Step 10.8's full parallel pool will handle analysis at maximum
                # concurrency (TARGET_CORES={TARGET_CORES} workers). No blocking wait.
                console_info(f" -> {pre_predicted_count:,} jobs have pre-existing GPU predictions.")
                console_info(f" -> Routing to parallel analysis pool ({TARGET_CORES} workers) — skipping single-worker queue.")
        if _gpu_only_jobs:
            console_info(f" -> {len(_gpu_only_jobs):,} jobs require GPU prediction.")

        if pending_by_protein:
            _n_gpu_jobs = sum(len(v) for v in pending_by_protein.values())
            print(f"\n{ConsoleColours.OKGREEN}{ConsoleColours.BOLD}{'=' * 80}{ConsoleColours.ENDC}")
            print(f"{ConsoleColours.OKGREEN}{ConsoleColours.BOLD}  GPU PREDICTION PHASE STARTING  —  {_n_gpu_jobs} jobs  |  {len(pending_by_protein)} proteins  |  {msa_ready}/{len(pending_proteins)} MSA cached{ConsoleColours.ENDC}")
            print(f"{ConsoleColours.OKGREEN}{ConsoleColours.BOLD}{'=' * 80}{ConsoleColours.ENDC}\n", flush=True)
            console_info("Launching Boltz-2 Batch Prediction Sequence Operations (GPU)...")
            console_info("Execution Strategy: Greedy batching protocol — multiple individual proteins bundled per Boltz call to systematically minimise model reloads and mitigate GPU idle instances.")
        else:
            console_info(f"\n -> All {pending_jobs} pending job(s) have pre-existing GPU predictions — no new predictions required.")
            console_info(f" -> Proceeding directly to parallel analysis...")

        msa_dir = D_COLABFOLD
        msa_dir.mkdir(parents=True, exist_ok=True)
        tmp_msa_dir = workspace_dir / "5b_MSA_Temp"
        tmp_msa_dir.mkdir(exist_ok=True)

        if pending_by_protein:
            console_info("\nStarting Concurrent MSA Retrieval Sequence and GPU Prediction Execution Pipeline...")
            console_info("MSA retrieval (via ColabFold APIs) executes consistently in a background parallel thread whilst Boltz-2 predicts completed sequential batches via the GPU architecture.")

        def background_downloader(proteins_to_fetch):
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
                    with open(yaml_path, 'w') as f:
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
                            a3ms = list(job_msa_dir.rglob("*.a3m"))
                            if a3ms and a3ms[0].stat().st_size > 500:
                                if validate_a3m_file(a3ms[0]):
                                    shutil.copy2(a3ms[0], a3m_path)
                                    with open(meta_path, 'w') as mf:
                                        json.dump({
                                            'sequence_sha256': sequence_hash(seq),
                                            'protein_id':      pid,
                                            'generated_at':    datetime.utcnow().isoformat() + 'Z',
                                            'source':          'boltz_fallback',
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

                    if not found and job_msa_dir.exists():
                        a3ms = list(job_msa_dir.rglob("*.a3m"))
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
            n_ligs_each = len(pending_by_protein[pid]) 

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
                    with open(job["yaml"], 'r') as yf: data = yaml.safe_load(yf)
                    data["sequences"][0]["protein"]["msa"] = str(batch_a3m.resolve())
                    chunk_yaml = chunk_dir / job["yaml"].name
                    with open(chunk_yaml, 'w') as yf: yaml.safe_dump(data, yf, sort_keys=False)
                    job_map[chunk_yaml.stem] = job
                    job_run_dir = D_RUNS / job_name
                    if job_run_dir.exists():
                        for old_res in job_run_dir.glob("boltz_results_*"):
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
                    _salvage_res_dir = next(protein_batch_out.glob("boltz_results_*"), None)
                    if _salvage_res_dir and (_salvage_res_dir / "predictions").exists():
                        _salvaged = set()
                        for _pred_dir in (_salvage_res_dir / "predictions").iterdir():
                            _stem = _pred_dir.name
                            if _stem not in job_map or not _pred_dir.is_dir():
                                continue
                            _job = job_map[_stem]
                            _jname = f"{_job['job_index']}_{_job.get('job_protein', _job['protein'])}_{_job['ligand']}"
                            _dst_root = D_RUNS / _jname / f"boltz_results_{_stem}"
                            _dst_pred = _dst_root / "predictions" / _stem
                            try:
                                _dst_pred.mkdir(parents=True, exist_ok=True)
                                for _f in _pred_dir.iterdir():
                                    shutil.move(str(_f), str(_dst_pred / _f.name))
                                for _cat in ["constraints", "mols", "msa", "records", "structures"]:
                                    _src_cat = _salvage_res_dir / "processed" / _cat
                                    if _src_cat.exists():
                                        for _ext in [".json", ".npz", ".pkl", ".csv"]:
                                            for _m in _src_cat.glob(f"{_stem}*{_ext}"):
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
                    _pid = _j['protein']
                    if _pid not in _protein_job_counts:
                        _protein_job_counts[_pid] = {'moved': 0, 'total': 0}
                    _protein_job_counts[_pid]['total'] += 1

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
                            _res_dir = next(protein_batch_out.glob("boltz_results_*"), None) if protein_batch_out.exists() else None
                            if _res_dir and (_res_dir / "predictions").exists():
                                for _pred_dir in list((_res_dir / "predictions").iterdir()):
                                    _stem = _pred_dir.name
                                    if not _pred_dir.is_dir() or _stem in _moved_stems or _stem not in job_map:
                                        continue
                                    if not any(_pred_dir.glob("*.cif")):
                                        continue  # not finished yet
                                    _job = job_map[_stem]
                                    _jname = f"{_job['job_index']}_{_job.get('job_protein', _job['protein'])}_{_job['ligand']}"
                                    _dst_root = D_RUNS / _jname / f"boltz_results_{_stem}"
                                    _dst_pred = _dst_root / "predictions" / _stem
                                    try:
                                        _dst_pred.mkdir(parents=True, exist_ok=True)
                                        for _f in list(_pred_dir.iterdir()):
                                            shutil.move(str(_f), str(_dst_pred / _f.name))
                                        for _cat in ["constraints", "mols", "msa", "records", "structures"]:
                                            _src_cat = _res_dir / "processed" / _cat
                                            if _src_cat.exists():
                                                for _ext in [".json", ".npz", ".pkl", ".csv"]:
                                                    for _m in _src_cat.glob(f"{_stem}*{_ext}"):
                                                        _d = _dst_root / "processed" / _cat
                                                        _d.mkdir(parents=True, exist_ok=True)
                                                        shutil.move(str(_m), str(_d / _m.name))
                                        gpu_queue_for_analysis.put(_job)
                                        _moved_stems.add(_stem)
                                        _pid2 = _job['protein']
                                        if _pid2 in _protein_job_counts:
                                            _protein_job_counts[_pid2]['moved'] += 1
                                            _pc = _protein_job_counts[_pid2]
                                            if _pc['moved'] == _pc['total']:
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
                        total_res_dir = next(protein_batch_out.glob("boltz_results_*"), None)
                        moved_count   = 0
                        failed_jobs   = []

                        if total_res_dir and (total_res_dir / "predictions").exists():
                            all_pred_folders = [d.name for d in (total_res_dir / "predictions").iterdir() if d.is_dir()]
                        else:
                            all_pred_folders = []

                        if len(all_pred_folders) == 0 and len(job_map) > 0:
                            console_info(f"\n    [BATCH CRITICAL WARNING] returncode=0 evaluated but 0 predictions yielded for [{batch_label}]!")
                            console_info(f"    [BATCH CRITICAL WARNING] Mathematically expected {len(job_map)} predictions but successfully obtained 0")
                            try:
                                with open(err_log_path, 'r') as ef:
                                    stderr_content = ef.read()
                                    if stderr_content:
                                        critical_lines = [
                                            line for line in stderr_content.split('\n')
                                            if any(kw in line for kw in ['Traceback', 'Error', 'KeyError', 'Exception', 'CUDA out', 'cuda', 'memory'])
                                            and '|' not in line and '%' not in line
                                        ]
                                        if critical_lines:
                                            console_info("    [BATCH STDERR EVALUATION] Critical Errors Identified Documented Below:")
                                            for line in critical_lines[-20:]:
                                                if line.strip(): console_info(f"      {line[:120]}")
                                        else:
                                            non_prog = [l for l in stderr_content.split('\n')
                                                        if l.strip() and '|' not in l and '%' not in l and 'it/s' not in l]
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
                                        matches = [d for d in pred_base.iterdir()
                                                   if d.is_dir() and (job['ligand'] in d.name or stem in d.name)]
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
                                        for f in src_pred.iterdir():
                                            shutil.move(str(f), str(dst_pred_dir / f.name))
                                        for cat in ["constraints", "mols", "msa", "records", "structures"]:
                                            src_cat = total_res_dir / "processed" / cat
                                            if src_cat.exists():
                                                for ext in [".json", ".npz", ".pkl", ".csv"]:
                                                    for m in src_cat.glob(f"{stem}*{ext}"):
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

                    # --- Surgical repair: if a single job caused a FileNotFoundError,
                    #     wait 5 s for the filesystem to settle, then isolate just that
                    #     job to the individual queue and retry the rest of the batch
                    #     without counting this as a failed attempt. ---
                    if not _surgical_repair_done and "FileNotFoundError" in stderr_text:
                        _m = re.search(r'/predictions/([^/\s]+)/pre_affinity_', stderr_text)
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
                                bad_job['fallback_individual'] = True
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
                                job['fallback_individual'] = True
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

    # Signal the background analysis worker to stop and wait for it to flush
    # (runs regardless of whether there were pending GPU jobs)
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

    worker_count = max(args.cpus if hasattr(args, 'cpus') and args.cpus else TARGET_CORES, 1)

    exec_args = []
    for job in tasks:
        job_name = f"{job['job_index']}_{job.get('job_protein', job['protein'])}_{job['ligand']}"
        if job_name in already_done_jobs:
            continue
        prev_time = prev_elapsed_map.get(job_name, 0.0)
        pf = preserved_fields_map.get(job_name, None) if resumed else None
        # Pass gpu_queue=None: Step 10.8 only handles pre-existing predictions (run_gpu=False).
        # Passing a Manager proxy would cause all forked workers to simultaneously reconnect
        # to the Manager RPC server → deadlock. None is safe here.
        exec_args.append((job, PROD, args.diffusion_samples, prev_time, D_ALN, control_cif, control_map, DEHA4_CONTROL_SEQ, None, pf, r3u_boltz_cif, r3u_control_map))

    total_tasks = len(exec_args)

    print()
    console_separator()
    console_info(f"Launching Multiprocessing Grid Environment executing securely on {worker_count} Concurrent Dedicated Workers...")
    console_info(f" -> {total_tasks:,} jobs queued for CPU re-analysis | {len(already_done_jobs):,} already complete")
    console_info(f" -> Using Pool.imap_unordered (lazy feed — no bulk queue pre-fill).")
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

    _pf_backup_cleanup = PROD / "_preserved_fields_backup.json"
    if _pf_backup_cleanup.exists():
        try:
            _pf_backup_cleanup.unlink()
        except Exception:
            pass

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
                df_temp.to_csv(CSV_PATH, index=False)
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
            # Field-specific overrides — must be set after the generic loop.
            fill_dict["residues_within_6A"]          = "0"
            fill_dict["hydrophobic_desolvation_ratio"] = 0.0
            fill_dict["product_inhibition_penalty"]    = 0.0
            fill_dict["r3u_Active_Site_RMSD"]              = 99.0
            fill_dict["r3u_Mechanistic_Fingerprint_Score"] = 0.0
            fill_dict["r3u_Halide_Stabilisation"]          = False
            fill_dict["r3u_Carboxylate_Clamp"]             = False
            fill_dict["r3u_ActiveSite_Conservation_Score"]     = 0.0

            df = df.fillna(fill_dict)
            df.rename(columns=COLUMN_RENAMING_MAP, inplace=True)
            df = df.loc[:, ~df.columns.duplicated()]
            
            drop_empty = []
            zero_equivalents = {'0', '0.0', '0.00', 'NA', 'None', '', 'nan', 'False'}
            protected_cols = {'degrader_tier', 'is_degrader', 'status', 'Protein_Name', 'Ligand_Name', 'job_name', 'job_index'}
            for col in df.columns:
                if col in protected_cols:
                    continue
                if df[col].astype(str).isin(zero_equivalents).all():
                    drop_empty.append(col)
            if drop_empty: df.drop(columns=drop_empty, inplace=True)

            if "job_index" in df.columns:
                df = df.sort_values("job_index", kind="mergesort").reset_index(drop=True)

            # GEM-PERF-1: Vectorised alignment grade assignment via pd.cut.
            # Derives Alignment_Score_Pct from identity_pct, then maps to
            # letter grades A–I using CFG-defined bins — replaces per-job
            # if-elif logic that would otherwise execute 58k+ times.
            if "identity_pct" in df.columns:
                df["Alignment_Score_Pct"] = (
                    pd.to_numeric(df["identity_pct"], errors="coerce")
                    .fillna(0.0)
                    .clip(0, 100)
                )
                df["Alignment_Grade"] = pd.cut(
                    df["Alignment_Score_Pct"],
                    bins=CFG.ALIGN_GRADE_BINS,
                    labels=CFG.ALIGN_GRADE_LABELS,
                    right=True
                ).astype(str)

            df.to_csv(CSV_PATH, index=False)
            GLOBAL_STATS["created_csv_rows"] = len(df)
            df_columns_count = len(df.columns)
        except Exception as e: console_info(f" Warning generated during systematic CSV formatting execution logic: {e}")
    console_info("Final structured CSV format written flawlessly.")
    
    # -------------------------------------------------------------------------------
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
            
            bind_col = "Binding_Probability" if "Binding_Probability" in df_rank.columns else "Binding_Probability_Score"
            if bind_col not in df_rank.columns: df_rank[bind_col] = 0.0
            
            for c in ["ActiveSite_Conservation_Score", "Mechanistic_Fingerprint_Score"]:
                if c not in df_rank.columns: df_rank[c] = 0.0
            if 'degrader_tier' not in df_rank.columns: df_rank['degrader_tier'] = CFG.TIER_DECOY

            # SN2 angle column may be named 'sn2_attack_angle' (pre-rename) or
            # 'SN2_Attack_Angle' (post-rename) depending on pipeline path.
            _sn2_col = next((c for c in df_rank.columns
                             if c.lower() in ('sn2_attack_angle', 'sn2attackangle')), None)
            if _sn2_col:
                df_rank['sn2_score_norm'] = (df_rank[_sn2_col].fillna(0.0) / 180.0).clip(0.0, 1.0)
            else:
                df_rank['sn2_score_norm'] = 0.0

            tier_map = CFG.TIER_SORT_WEIGHT
            df_rank['tier_val'] = df_rank['degrader_tier'].map(tier_map).fillna(0)

            df_rank = df_rank.sort_values(
                by=["tier_val", "Mechanistic_Fingerprint_Score", "ActiveSite_Conservation_Score",
                    "sn2_score_norm", bind_col],
                ascending=[False, False, False, False, False]
            )
            df_rank.insert(0, 'Scientific_Rank', range(1, len(df_rank) + 1))
            df_rank['Ranking_Score_Calc'] = (
                "Tier:"  + df_rank['degrader_tier'].astype(str)
                + " | Mech:" + df_rank['Mechanistic_Fingerprint_Score'].map('{:.2f}'.format)
                + " | Like:" + df_rank['ActiveSite_Conservation_Score'].map('{:.2f}'.format)
            )

            cols_to_drop = ['tier_val', 'sn2_score_norm']
            df_rank = df_rank.drop(columns=[c for c in cols_to_drop if c in df_rank.columns])
            df_rank.to_csv(rank_csv_path, index=False)
            rank_columns_count = len(df_rank.columns)
            console_info(f"Scientific Ranking Output CSV generated flawlessly.")
        else: console_info(" Operation actively bypassed (No structured data detected).")
    except Exception as e: console_info(f" Mechanism ranking operation structurally failed ({e})")

    save_alignment_cache_final(D_ALN / "Alignment_Stats.csv")

    console_separator()

    # -------------------------------------------------------------------------------
    # Step 10.11: Final Report and Output Terminal UI
    # -------------------------------------------------------------------------------
    console_separator()
    console_info(f"\nSequential Analytical Pipeline entirely Completed. Associated Files are formally verified inside: {PROD.resolve()}")
    console_separator()
    console_info(" \n ✔ Summary of Validated Operational Disk Files-")
    
    yaml_count = len(list(D_YAML.glob("*.yaml")))
    job_count = len([d for d in D_RUNS.iterdir() if d.is_dir()])
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
    console_info(f"Sequence Reference Data (MSA + Alignments):".ljust(45) + f" |      {D_SEQ.resolve()}")
    console_info(f"  ├─ MSA Sequences: {len(list(D_COLABFOLD.glob('*.a3m')))}" .ljust(45) + f" |      {D_COLABFOLD.resolve()}")
    console_info(f"  └─ Alignments: {str_aln}".ljust(45) + f" |      {D_ALN.resolve()}")
    console_info(f"FAcDs Master CSV:   {str_mc}".ljust(45) + f" |      {CSV_PATH.resolve()}")
    if rank_csv_path:
        console_info(f"FAcDs Ranked CSV:   {str_rc}".ljust(45) + f" |      {rank_csv_path.resolve()}")
    _mir = run_root / "2_Best_Complexes_CIFs"
    _mir_count = sum(1 for _ in _mir.rglob("*.cif")) if _mir.exists() else 0
    console_info(f"Best Complexes CIFs Mirror: {_mir_count} CIFs".ljust(45) + f" |      {_mir.resolve()}")
    console_info(f"Operational Execution Log".ljust(45) + f" |      {LOG_PATH.resolve()}")
    console_separator()
    
    if CSV_PATH.exists(): df_final = pd.read_csv(CSV_PATH, low_memory=False)
    else: df_final = pd.DataFrame()
    
    if not df_final.empty and CFG.COL_TIER in df_final.columns:
        tiers = Counter(df_final[CFG.COL_TIER].fillna(CFG.TIER_DECOY).tolist())
    else:
        tiers = Counter()
    
    # Use Alignment_Grade column (written by pd.cut in Step 10.9)
    # when available; fall back to utility logic for backwards
    # compatibility with CSVs produced before this change was applied.
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
        grade_A = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == 'A').sum())
        grade_B = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == 'B').sum())
        grade_C = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == 'C').sum())
        grade_D = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == 'D').sum())
        grade_E = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == 'E').sum())
        grade_F = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == 'F').sum())
        grade_G = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == 'G').sum())
        grade_H = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == 'H').sum())
        grade_I = int(id_pct.apply(lambda x: _utils_mod.get_alignment_grade(x, CFG) == 'I').sum())

    mdl_col = "best_model_name" if "best_model_name" in df_final.columns else "Boltz_Model_Confidence"
    if not df_final.empty and mdl_col in df_final.columns:
        models = Counter(df_final[mdl_col].fillna('None').tolist())
    else:
        models = Counter()

    console_info("")
    console_info(SEPARATOR_LIGHT)
    console_info(f"  Final Computational Statistics")
    console_info(SEPARATOR_LIGHT)
    if not df_final.empty:
        csv_failures  = int((df_final['status'] != 'Success').sum())
        grid_failures = GLOBAL_STATS.get("analytical_failures", 0)
        total_failures = csv_failures + grid_failures
        _jw = 36
        console_info(f"  ┌{'─'*(_jw+2)}┬{'─'*12}┐")
        console_info(f"  │  {'Metric':<{_jw}}│  {'Count':>8}  │")
        console_info(f"  ├{'─'*(_jw+2)}┼{'─'*12}┤")
        console_info(f"  │  {'Jobs Accounted':<{_jw}}│  {len(df_final)+grid_failures:>8,}  │")
        console_info(f"  │  {'Successful Outputs':<{_jw}}│  {int((df_final['status']=='Success').sum()):>8,}  │")
        _fail_str = f"{total_failures}" + (f"  ← re-run to retry" if grid_failures else "")
        console_info(f"  │  {'Analytical Failures':<{_jw}}│  {total_failures:>8,}  │")
        console_info(f"  └{'─'*(_jw+2)}┴{'─'*12}┘")

        _tier_order = CFG.TIER_ORDER
        _tier_total = sum(tiers.get(t, 0) for t in _tier_order)
        _tw = 14
        console_info("")
        console_info(f"  ┌{'─'*(_tw+2)}┬{'─'*12}┐")
        console_info(f"  │  {'Degrader Tier':<{_tw}}│  {'Count':>8}  │")
        console_info(f"  ├{'─'*(_tw+2)}┼{'─'*12}┤")
        for _t in _tier_order:
            console_info(f"  │  {_t:<{_tw}}│  {tiers.get(_t,0):>8,}  │")
        console_info(f"  ├{'─'*(_tw+2)}┼{'─'*12}┤")
        console_info(f"  │  {'Total':<{_tw}}│  {_tier_total:>8,}  │")
        console_info(f"  └{'─'*(_tw+2)}┴{'─'*12}┘")

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
