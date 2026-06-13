"""
===============================================================================
FAcDs Pipeline  |  Module 00_02  |  Central Configuration (CFG)
===============================================================================
Single source of truth for every threshold, constant, weight, and parameter
used across the pipeline. Edit values here only — no other file should contain
hard-coded scientific values, configurable thresholds, or tunable settings.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 10 June 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Module        : 00_02_Project_Config_FAcDs.py
  Role          : "Blueprint" — pipeline-wide configuration repository.
  Imported by   : All pipeline scripts (00_02 through 07).
  Reads         : Nothing (pure Python dataclass; no file I/O).
  Writes        : Nothing.
  Upstream      : None (root module — must load before all others).
  Downstream    : 02_Production_FAcDs.py, 03_Validation_Figures_FAcDs.py,
                  04_Phylogeny_FAcDs.py, 05_CIF-PDB_Preparation_FAcDs.py,
                  06_Top-N_Extraction_FAcDs.py,
                  08_MD_Thermodynamics_QMMM_Engine_FAcDs.py
───────────────────────────────────────────────────────────────────────────────

# ── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Module Load Priority: This file MUST be imported before any other project
     modules to ensure CFG is available for dependent logic.
  2. Static Constancy: Changing thresholds (e.g., NAC_DIST) mid-pipeline will
     cause inconsistencies between stored CSV results and new analysis.
  3. Environment Sensitivity: Relies on `importlib.util` for path-safe imports;
     requires consistent relative directory structure.
───────────────────────────────────────────────────────────────────────────────

Usage:
    import importlib.util as _ilu
    _spec = _ilu.spec_from_file_location("ProjectConfig", path / "00_02_Project_Config_FAcDs.py")
    _mod  = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_mod)
    CFG   = _mod.CFG()   # instantiate once at module level

    CFG.NAC_DIST_STRICT   # access any attribute directly
    CFG.THRESHOLD_HB_DIST_MAX

Section map (prefix → section):
    BOLTZ_*        §2   Boltz-2 prediction engine & ColabFold MSA
    THRESHOLD_*    §4   Biochemical interaction geometry cutoffs
    NAC_*          §5   Near Attack Conformation geometry
    TIER_*         §9   Catalytic tier classification & display
    WATERMAP_*     §10  WaterMap hydration-site integration
    SCORE_*        §10/15  Likelihood scoring thresholds
    QSITE_*        §11  QM/MM extraction settings (QSite/Schrödinger)
    SMART_LOCK_*   §8   3D Smart-Lock residue-detection bias
    PROC_*         §13  Processing / file-I/O operational parameters
    VIS_*          §14  Visualisation rendering parameters
    ALIGN_*        §15  Sequence alignment penalty parameters
    GPU_*          §15  GPU batch watchdog timeout
    PREPWIZARD_*   §16  Schrödinger PrepWizard invocation parameters
    PREP_*         §16  PDB preparation chain & residue classification
───────────────────────────────────────────────────────────────────────────────

Scientific references
─────────────────────
  Boltz-2 structure prediction : Passaro, S. et al. (2025) bioRxiv 2025.06.14.659707. https://doi.org/10.1101/2025.06.14.659707
  ColabFold MSA server         : Mirdita, M. et al. (2022) Nature Methods 19:679–682. https://doi.org/10.1038/s41592-022-01488-1
  FAcD crystal structure (3R3U): Chan, P.W.Y. et al. (2011) JACS 133:7461–7468. https://doi.org/10.1021/ja200277d
  DEHA4 defluorination (D4B)   : Farajollahi, S. et al. (2024) ACS Omega 9(26):28546–28555. https://doi.org/10.1021/acsomega.4c02517
  H-bond geometry (D···A)      : Jeffrey, G.A. (1997) An Introduction to Hydrogen Bonding. Oxford University Press.
  H-bond geometry (H···A)      : Baker & Hubbard (1984) Prog Biophys Mol Biol 44:97–179. https://doi.org/10.1016/0079-6107(84)90007-5
  Salt bridge                  : Barlow & Thornton (1983) J Mol Biol 168:867–885.
  Hydrophobic contact          : Salentin et al. (2015) Nucleic Acids Res 43:W443–W447. https://doi.org/10.1093/nar/gkv315
  π–π stacking                 : McGaughey et al. (1998) J Biol Chem 273:15458–15463. https://doi.org/10.1074/jbc.273.25.15458
  π–cation                     : Gallivan & Dougherty (1999) PNAS 96:9459–9464. https://doi.org/10.1073/pnas.96.17.9459
  Halogen bond                 : Wilcken et al. (2013) J Med Chem 56:1363–1388. https://doi.org/10.1021/jm3012068
  Metal coordination           : Harding (2006) Acta Crystallogr D62:678–682. https://doi.org/10.1107/S0907444906014594
  NAC criteria                 : Lightstone & Bruice (1996) JACS 118:2595–2605. https://doi.org/10.1021/ja952589l
                                 Bruice (2002) Acc Chem Res 35:139–148. https://doi.org/10.1021/ar0001665
  Catalytic triad distances    : Holmquist (2000) Curr Protein Pept Sci 1:209–235. https://doi.org/10.2174/1389203003381405
  Bürgi–Dunitz angle (aux)     : Bürgi, Dunitz & Shefter (1973) JACS 95:5065–5067. https://doi.org/10.1021/ja00796a058
                                 Bürgi, Dunitz, Lehn & Wipff (1974) Tetrahedron 30:1563–1572. https://doi.org/10.1016/S0040-4020(01)90678-7
  WaterMap thermodynamics      : Abel, R. et al. (2008) JACS 130:2817–2831. https://doi.org/10.1021/ja0771033
  QSite DFT functional (M06-2X): Zhao, Y. & Truhlar, D.G. (2008) Theor Chem Acc 120:215–241. https://doi.org/10.1007/s00214-007-0310-x
  QSite QM/MM implementation   : Murphy, R.B. et al. (2000) J Comput Chem 21:1442–1457. https://doi.org/10.1002/1096-987X(200012)21:16<1442::AID-JCC3>3.0.CO;2-O
  QM/MM free-energy method     : Rosta, E. et al. (2006) J Phys Chem B 110:2934–2941. https://doi.org/10.1021/jp057109j  (general QM/MM benchmark, not QSite-specific)
  FAcD QM/MM defluorination    : Yue, Y. et al. (2021) Environ Sci Technol 55(14):9817–9825. https://doi.org/10.1021/acs.est.0c08811
  Sequence alignment           : Henikoff & Henikoff (1992) PNAS 89:10915–10919. https://doi.org/10.1073/pnas.89.22.10915
  PrepWizard protocol          : Madhavi Sastry et al. (2013) J Comput-Aided Mol Des 27:221–234. https://doi.org/10.1007/s10822-013-9644-8
===============================================================================
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class CFG:
    """
    Central Configuration Repository — FAcDs Pipeline.

    Attribute prefix → section
    ──────────────────────────
    BOLTZ_*        §2   Boltz-2 engine execution, sampling, ColabFold MSA
    THRESHOLD_*    §4   Biochemical interaction geometry cutoffs
    NAC_*          §5   Near Attack Conformation geometry
    TIER_*         §9   Catalytic tier classification & display
    WATERMAP_*     §10  WaterMap hydration-site analysis
    SCORE_*        §10/15  Likelihood scoring thresholds
    VIABILITY_*    §10  Catalytic viability display thresholds
    QSITE_*        §11  QM/MM extraction settings (QSite)
    SMART_LOCK_*   §8   3D Smart-Lock residue-detection bias
    PROC_*         §13  Processing / file-I/O operational parameters
    VIS_*          §14  Visualisation rendering parameters
    ALIGN_*        §15  Sequence alignment penalty parameters
    GPU_*          §15  GPU batch watchdog timeout
    PREPWIZARD_*   §16  Schrödinger PrepWizard invocation parameters
    PREP_*         §16  PDB preparation chain & residue classification
    """

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 1 ── PROJECT IDENTITY
    # ═════════════════════════════════════════════════════════════════════════════
    PROJECT_NAME: str   = "PFAS-27"
    PIPELINE_ROLE: str  = "Dehalogenase Evolution Tracking"

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 2 ── BOLTZ-2 PREDICTION ENGINE  (Step 02)
    # ═════════════════════════════════════════════════════════════════════════════

    # ── § 2.1  Executable & model ────────────────────────────────────────
    BOLTZ_EXECUTABLE: str          = "boltz"    # CLI binary name on PATH
    BOLTZ_MODEL_VERSION: str       = "boltz2"   # --model flag passed to boltz
    BOLTZ_OUTPUT_FORMAT: str       = "mmcif"    # --output_format flag

    # ── § 2.2  Sampling & batching ───────────────────────────────────────
    BOLTZ_RECYCLING_STEPS: int        = 1    # --recycling_steps per prediction
    BOLTZ_DIFFUSION_SAMPLES: int      = 5    # --diffusion_samples (structures per complex)
    BOLTZ_MAX_PROTEINS_PER_BATCH: int = 20   # proteins bundled in one boltz call
                                             # e.g. 20 proteins × 27 ligands = 540 jobs/batch

    # ── § 2.3  Retry & polling ───────────────────────────────────────────
    BOLTZ_RETRY_MAX: int     = 2   # max retries after GPU failure before abandoning job
    BOLTZ_RETRY_SLEEP: int   = 6   # seconds to wait between retry attempts
    BOLTZ_POLL_INTERVAL: int = 5   # seconds between job-status poll cycles

    # ── § 2.4  ColabFold MSA submission ──────────────────────────────────
    COLABFOLD_API_URL: str        = "https://api.colabfold.com"
    COLABFOLD_SUBMIT_RETRIES: int = 8    # max attempts to POST an MSA ticket
    COLABFOLD_MSA_SUBMIT_TIMEOUT: int       = 30   # seconds — HTTP POST timeout for ticket submission
    COLABFOLD_MSA_POLL_TIMEOUT: int         = 900  # seconds — max wait for MSA completion (15 min)
    COLABFOLD_MSA_POLL_REQUEST_TIMEOUT: int = 15   # seconds — HTTP GET timeout per polling cycle
    COLABFOLD_MSA_DOWNLOAD_TIMEOUT: int     = 120  # seconds — HTTP GET timeout for final tar.gz download

    # ── § 2.5  Confidence scoring weights ─────────────────────────────────
    # Weights for the composite Boltz confidence score used to rank models.
    # Higher weight = stronger contribution to final ranking.
    BOLTZ_W_IPTM: float         = 3.0   # interface predicted TM-score (primary signal)
    BOLTZ_W_CROSS_PAE: float    = 2.5   # cross-interface predicted aligned error
    BOLTZ_W_PLDDT: float        = 2.0   # per-residue local confidence
    BOLTZ_W_INTERACTIONS: float = 1.5   # density of physical interactions
    BOLTZ_W_CONF: float         = 1.0   # global structural confidence (ptm)

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 3 ── INPUT & REFERENCE DATA
    # ═════════════════════════════════════════════════════════════════════════════

    # ── § 3.1  Reference ligand ──────────────────────────────────────────
    # Fluoroacetate (FA) is the canonical FAcD substrate; benchmark for DEHA4 and 3R3U geometry calibration.
    FLUOROACETATE_SMILES: str = "C(C(=O)O)F"   # canonical SMILES for fluoroacetate
    INPUT_SMILES: str = "D_INP_PFAS-27_Ligands.smi"  # ligand SMILES panel

    # ── § 3.2  Reference protein sequence ────────────────────────────────
    # DEHA4 = DeHa4_[Delftia acidovorans D4B] — structural reference for all alignments.
    # Source: Farajollahi et al. (2024) ACS Omega https://doi.org/10.1021/acsomega.4c02517
    # Position cross-validated; establishes HIS277 (not HIS288) in DEHA4 sequences.
    DEHA4_CONTROL_SEQ: str = (
        "MHTDPWMPGLRQQRITVDDGVEINAWVGGQGPALLLVHGHPQTSAIWHRVAPRLAQQFTVVLADLRGYGDSSRPAGDPEH"
        "VNYSKRTMARDLLRLMARLGHEHFSVLAHDRGARVAHRLAMDYPASVQRLVLLDIAPTLAMYEQTGEAFARAYWHWFFLI"
        "QPAPLPERLIEADPAAYVREIMGRRSAGLAPFDPRALAEYQRCLALPGSAHGMCEDYRASAGIDLDHDREDRQLGRRLSM"
        "PLLVLWGEEGWHRCFDPLREWQLVADDVRGRPLACGHYIAEEAPDALLDAALPFLLQAG"
    )

    # ── § 3.3  Crystal structure reference — PDB 3R3U ────────────────────
    # Rhodopseudomonas palustris FAcD, wild-type, 1.60 Å resolution
    # (no substrate bound; structure carries Ni²⁺ and Cl⁻ ions).
    # Ref: Chan et al. (2011) JACS 133:7461–7468. DOI: 10.1021/ja200277d
    #      Jitsumori et al. (2009) J Bacteriol 191:2630–2637. DOI: 10.1128/JB.01654-08
    REFERENCE_PDB_ID: str    = "3R3U"
    REFERENCE_PDB_URL: str   = "https://files.rcsb.org/download/3R3U.pdb"
    REFERENCE_DIR_NAME: str  = "0_Reference_Crystal_3R3U"

    # ── § 3.4  Control ligand panel (fluoroacetate control substrates) ────
    # Three short-chain fluorinated substrates used for DEHA4 and 3R3U calibration runs.
    # Format: (job_suffix, SMILES)
    CTRL_LIGANDS: list = field(default_factory=lambda: [
        ("26_Fluoroacetate",   "C(C(=O)O)F"),         # fluoroacetate (canonical substrate)
        ("27_Difluoroacetate", "O=C(O)C(F)F"),        # difluoroacetate
        ("25_TFA",             "O=C(O)C(F)(F)F"),    # trifluoroacetate
    ])

    # ── § 3.5  Reference active-site mapping (3R3U / DEHA4 canonical) ────
    # Residue identities, sequential IDs, and PDB-numbered IDs for the
    # eight key catalytic / stabilising positions in the FAcD active site.
    REF_ACTIVE_SITE_MAP: dict = field(default_factory=lambda: {
        "Nuc":    {"res": "ASP", "id": 110, "pdb_id": 110, "role": "Nucleophile"},
        "Carb1":  {"res": "ARG", "id": 111, "pdb_id": 111, "role": "Carboxylate_Clamp"},
        "Carb2":  {"res": "ARG", "id": 114, "pdb_id": 114, "role": "Carboxylate_Clamp"},
        "Acid":   {"res": "ASP", "id": 134, "pdb_id": 134, "role": "Acid_Catalyst"},
        "Stab_H": {"res": "HIS", "id": 155, "pdb_id": 155, "role": "Fluorine_Stabiliser"},
        "Stab_W": {"res": "TRP", "id": 156, "pdb_id": 156, "role": "Fluoride_Cradle"},
        "Stab_Y": {"res": "TYR", "id": 217, "pdb_id": 219, "role": "Fluoride_Cradle"},
        "Base":   {"res": "HIS", "id": 277, "pdb_id": 280, "role": "Base_Catalyst"},
    })

    # ── § 3.6  Catalytic triad key set ────────────────────────────────────
    CATALYTIC_TRIAD_KEYS: list = field(default_factory=lambda: ["Nuc", "Base", "Acid"])

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 4 ── INTERACTION GEOMETRY THRESHOLDS  (Step 02)
    # All distances are heavy-atom unless labelled (HA = H-to-acceptor).
    # ═════════════════════════════════════════════════════════════════════════════

    # ── § 4.1  Hydrogen bonds ─────────────────────────────────────────────
    # Jeffrey (1997): D···A ≤ 3.5 Å; Baker & Hubbard (1984): H···A ≤ 2.8 Å
    THRESHOLD_HB_DIST_MAX: float   = 3.5    # Å  donor–acceptor (D···A) maximum
    THRESHOLD_HB_DIST_MIN: float   = 2.4    # Å  donor–acceptor (D···A) minimum
    THRESHOLD_HB_DIST_HA: float    = 2.8    # Å  hydrogen-to-acceptor (H···A) maximum
    THRESHOLD_HB_ANGLE_MIN: float  = 120.0  # °  D–H···A minimum angle

    # ── § 4.2  Salt bridges ───────────────────────────────────────────────
    # Barlow & Thornton (1983): opposite-charge heavy-atom ≤ 4.0 Å.
    THRESHOLD_SALT_BRIDGE: float   = 4.0    # Å  heavy-atom distance between charged groups

    # ── § 4.3  Hydrophobic contacts ───────────────────────────────────────
    THRESHOLD_HYDROPHOBIC_MAX: float = 4.0  # Å  carbon–carbon centroid distance

    # ── § 4.4  π–π stacking ───────────────────────────────────────────────
    # McGaughey et al. (1998): face–face ≤ 4.4 Å; edge–face ≤ 5.5 Å.
    THRESHOLD_PI_FACE: float       = 4.4    # Å  centroid–centroid, face–face (alias: PI_STACK_FACE_DIST_MAX)
    THRESHOLD_PI_EDGE: float       = 5.5    # Å  centroid–centroid, edge–face (alias: PI_STACK_EDGE_DIST_MAX)
    # Backward-compatibility aliases: reference the THRESHOLD_PI_* source above
    # (not independent literals) so the pair can never silently diverge.
    PI_STACK_FACE_DIST_MAX: float  = THRESHOLD_PI_FACE  # Å  alias of THRESHOLD_PI_FACE
    PI_STACK_FACE_ANGLE_MAX: float = 30.0   # °  maximum tilt angle, face–face
    PI_STACK_EDGE_DIST_MAX: float  = THRESHOLD_PI_EDGE  # Å  alias of THRESHOLD_PI_EDGE
    PI_STACK_EDGE_ANGLE_MIN: float = 60.0   # °  minimum tilt angle, edge–face

    # ── § 4.5  π–cation ───────────────────────────────────────────────────
    # Gallivan & Dougherty (1999): centroid-to-cation ≤ 6.0 Å; off-axis ≤ 30°.
    THRESHOLD_PI_CATION_MAX: float = 6.0    # Å  ring-centroid to cation distance
    PI_CATION_ANGLE_MAX: float     = 30.0   # °  maximum off-axis angle

    # ── § 4.6  Halogen bonds ──────────────────────────────────────────────
    # Wilcken et al. (2013): X···A ≤ 3.5 Å.
    HALOGEN_BOND_DIST_MAX: float   = 3.5    # Å  halogen···acceptor (X···A) maximum

    # ── § 4.7  Metal coordination ─────────────────────────────────────────
    # Harding (2006): metal–ligand bond ≤ 2.8 Å.
    METAL_COORD_DIST_MAX: float    = 2.8    # Å  metal–ligand coordination bond
    METALS: set = field(default_factory=lambda: {
        "MG", "ZN", "MN", "CA", "FE", "CO", "NI", "CU", "NA", "K",
    })

    # ── § 4.8  General catalytic site ────────────────────────────────────
    CATALYTIC_DIST_CUTOFF: float   = 6.0    # Å  residue included as "near active site"

    # ── § 4.9  Coordinate–structure match ────────────────────────────────
    COORD_MATCH_DIST_MAX: float    = 1.8    # Å  map residue to reference coordinate

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 5 ── NAC (Near Attack Conformation) GEOMETRY  (Steps 02, 08)
    # Substrate: fluoroacetate; electrophile: C–F carbon; nucleophile: Asp O.
    # Ideal SN2 backside attack: Nu–C–F collinear at 180°.
    # References: Lightstone & Bruice (1996); Bruice (2002).
    # ═════════════════════════════════════════════════════════════════════════════

    # ── § 5.1  Strict NAC — publication-grade catalytic viability ─────────
    NAC_DIST_STRICT: float   = 3.2    # Å  nucleophile O to electrophilic C
    NAC_ANGLE_STRICT: float  = 155.0  # °  O–C–F attack angle at C (180° = ideal backside)

    # ── § 5.2  Relaxed NAC — pre-reactive / entropy-inclusive sampling ────
    NAC_DIST_RELAXED: float  = 3.8    # Å  relaxed nucleophile–C distance
    NAC_ANGLE_RELAXED: float = 145.0  # °  relaxed attack angle

    # ── § 5.3  Pocket residency ───────────────────────────────────────────
    POCKET_RESIDENCY_DIST: float = 8.0   # Å  ligand considered "pocket-bound" when nuc–C ≤ this

    # ── § 5.4  Fluoride cradle detection (TRP / TYR aromatic basket) ──────
    # The fluoride cradle is the TRP/TYR aromatic shell that stabilises the
    # departing fluoride in FAcD-family dehalogenases.
    F_CRADLE_RADIUS: float   = 12.0   # Å  search radius from nucleophile to TRP/TYR heavy atoms

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 6 ── CATALYTIC TRIAD INTEGRITY  (Steps 02, 08)
    # Reference: fluoroacetate dehalogenase crystal structure (PDB 3R3U; Chan et al. 2011 JACS. DOI: 10.1021/ja200277d).
    # NB = nucleophile O to catalytic base N; BA = catalytic base N to acid O.
    #
    # Two threshold sets:
    #   Static (Step 02) — calibrated on crystal structures (energy-minimised).
    #   MD     (Step 08) — +2.0 Å / +2.0 Å buffer for 300 K thermal fluctuations in
    #                      solution; justified by Asp–His distance variance in FAcD MD
    #                      trajectories (σ ≈ 1–2 Å at 300 K).
    # ═════════════════════════════════════════════════════════════════════════════
    THRESHOLD_TRIAD_NB: float    = 4.5   # Å  crystal/static (Step 02)
    THRESHOLD_TRIAD_BA: float    = 7.0   # Å  crystal/static (Step 02)
    THRESHOLD_TRIAD_NB_MD: float = 6.5   # Å  MD-calibrated  (Step 08) = crystal + 2.0 Å
    THRESHOLD_TRIAD_BA_MD: float = 9.0   # Å  MD-calibrated  (Step 08) = crystal + 2.0 Å

    # ── § 6.1  Mechanistic-fingerprint & soft-score contact gates  (Step 02) ──
    # Used by 02_Production analyse_candidate_structure() mech_score / soft_score.
    # These are distinct from the §9 tier-cascade gates: they score the physical
    # anchor set (halide stabiliser, carboxylate clamp, triad proximity), whereas
    # §9 assigns the discrete tier. The Nuc–C gate reuses NAC_DIST_STRICT (§5.1)
    # and the angle sigmoid reuses NAC_ANGLE_STRICT so all three stay in lock-step.
    MECH_STAB_RADIUS: float  = 5.5   # Å  TRP/TYR (or dynamic polar) → F⁻ halide-stabilisation contact
    MECH_CLAMP_RADIUS: float = 5.0   # Å  ARG carboxylate clamp → ligand contact
    MECH_NB_GATE: float      = 5.0   # Å  Nuc–Base distance gate (+0.1 mech point)
    MECH_BA_GATE: float      = 5.5   # Å  Base–Acid distance gate (+0.1 mech point)
    # Soft-score (sigmoid) triad midpoints; nucleophile/angle sigmoids reuse the
    # strict NAC cutoffs directly (NAC_DIST_STRICT / NAC_ANGLE_STRICT).
    SOFT_NB_MIDPOINT: float  = 4.5   # Å  soft s_int Nuc–Base sigmoid midpoint (= THRESHOLD_TRIAD_NB)
    SOFT_BA_MIDPOINT: float  = 5.0   # Å  soft s_int Base–Acid sigmoid midpoint
                                     #    Intentionally < THRESHOLD_TRIAD_BA (7.0 Å): sigmoid midpoint sets
                                     #    the steepest scoring gradient in the 4–6 Å pre-reactive range;
                                     #    the hard 7.0 Å cutoff is the gate, not the sigmoid centre.

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 7 ── SN2 / WALDEN INVERSION GEOMETRY  (Step 08)
    # SN2 backside attack at sp³ C: ideal Walden-inversion trajectory = 180°.
    # Note: the Bürgi–Dunitz angle (107°) applies to nucleophilic addition at
    # sp² carbonyl carbons; it must NOT be conflated with the linear SN2 angle
    # here (see calculate_burgi_dunitz() in 00_03 — auxiliary metric only).
    # ═════════════════════════════════════════════════════════════════════════════
    WALDEN_IMPROPER_MAX: float = 15.0  # °  |improper dihedral| < this → TS-like (planar) geometry

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 8 ── SMART-LOCK RESIDUE DETECTION  (Step 08)
    # Geometry-biased scoring that steers 3D triad & fluoride-cradle detection
    # toward known residue positions from the 3R3U canonical mapping.
    # Negative values are distance bonuses (subtracted from the candidate distance;
    # lower effective score = better candidate).
    #
    # NB: the bias is a *soft preference*, not a hard lock. It is sized as a few Å —
    # comparable to inter-residue spacing — so a clearly closer geometric candidate
    # can still win over the sequence-aligned hint. (Previously −100/−50 Å, which
    # made the hint effectively deterministic and could force the wrong residue when
    # the alignment hint was wrong; SMART_LOCK_NUC_MAX_DIST and the Nuc–Base sanity
    # check remain the hard safety nets.)
    # ═════════════════════════════════════════════════════════════════════════════
    SMART_LOCK_BIAS_DIST: float    =   -5.0  # Å  soft bonus when residue matches mapped-hint position
    SMART_LOCK_CHAIN_BIAS: float   =   -2.5  # Å  additional soft bonus for same-chain residue
    SMART_LOCK_RESNUM_WINDOW: int  =    15   # residue-number window (±N) around each hint
    SMART_LOCK_NUC_MAX_DIST: float =   20.0  # Å  hard cutoff — no nucleophile candidate beyond this
    SMART_LOCK_OD_FALLBACK_DIST: float  =  6.0   # Å  Oδ–Cα search radius for nucleophile Oδ fallback
    SMART_LOCK_OD_FALLBACK_ANGLE: float = 120.0  # °  minimum backside approach angle (anti-F) in Oδ fallback
    SMART_LOCK_MAPPING_SANITY_DIST: float = 15.0  # Å  max plausible Nuc–Base dist; beyond = mapping error

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 9 ── CATALYTIC TIER CLASSIFICATION  (Steps 02, 03, 04, 06)
    #
    # The tier cascade is evaluated in strict descending order.
    # The top tier is index 0 in TIER_ORDER, the lowest is the last element.
    # A structure is assigned the highest tier whose ALL criteria are met.
    # Criteria: nucleophile distance + attack angle + triad distances + mech score.
    # ═════════════════════════════════════════════════════════════════════════════

    # --- Dynamic Tier Taxonomy ---
    # Alphanumeric sorting dictates the strict hierarchy.
    # Tier_1A automatically sorts to index 0 (Top Tier), Tier_5_Decoy to the bottom.
    TIER_NAMES: list[str] = field(default_factory=lambda: [
        "Tier_1A", "Tier_1B", "Tier_2A", "Tier_2B", 
        "Tier_3", "Tier_4", "Tier_5_Decoy"
    ])

    @property
    def TIER_ORDER(self) -> list[str]:
        # Sort alphanumerically to guarantee hierarchy without hardcoding names.
        return self.TIER_NAMES
    @property
    def TIER_TOP(self) -> str: return self.TIER_ORDER[0]
    @property
    def TIER_POOR(self) -> str: return self.TIER_ORDER[-2]
    @property
    def TIER_DECOY(self) -> str: return self.TIER_ORDER[-1]
    @property
    def TIER_HIGH_QUALITY(self) -> list: return self.TIER_ORDER[:4]


    # ── § 9.1  Nucleophile–C distance thresholds (Å, upper bound) ────────
    TIER_NUC_DIST: dict = field(default_factory=lambda: {
        "Tier_1A": 3.0,   # tight pre-reactive geometry; 2.7 Å unrealistically tight for Boltz-2 ground-state (true SN2 TS ~ 2.0–2.3 Å)
        "Tier_1B": 3.2,   # excellent pre-reactive geometry
        "Tier_2A":    3.2,   # = NAC_DIST_STRICT
        "Tier_2B":    3.8,   # = NAC_DIST_RELAXED
        "Tier_3":      4.2,   # productive but not pre-reactive
        "Tier_4":      8.0,   # pocket-bound, geometrically unproductive
    })

    # ── § 9.2  Attack angle minimum thresholds (°, lower bound) ──────────
    TIER_ANGLE_MIN: dict = field(default_factory=lambda: {
        "Tier_1A": 175.0,
        "Tier_1B": 165.0,
        "Tier_2A":    155.0,   # = NAC_ANGLE_STRICT
        "Tier_2B":    145.0,   # = NAC_ANGLE_RELAXED
    })

    # ── § 9.3  Catalytic triad distance thresholds (Å, upper bound) ──────
    TIER_NB_MAX: dict = field(default_factory=lambda: {
        "Tier_1A": 3.5,   # tightest triad — Nuc–Base ≤ 3.5 Å
        "Tier_1B": 4.0,
        "Tier_2A":    5.0,
    })
    TIER_BA_MAX: dict = field(default_factory=lambda: {
        "Tier_1A": 4.5,   # Base–Acid ≤ 4.5 Å
        "Tier_1B": 5.0,
        "Tier_2A":    6.0,
    })

    # ── § 9.4  Mechanistic fingerprint score thresholds (0–1, lower bound)
    TIER_MECH_MIN: dict = field(default_factory=lambda: {
        "Tier_1A": 0.9,   # 1.0 is unreachable when any single anchor (clamp/stabiliser) is absent
        "Tier_1B": 0.7,
        "Tier_2A":    0.5,
    })

    # ── § 9.4b  Substrate / inhibitor classification (Figs 19, 25) ─────────
    # Derived from the tier gates so every figure classifies identically.
    SUBSTRATE_ANGLE_MIN: float = 165.0   # = TIER_ANGLE_MIN['Tier_1B']: SN2 ≥ this → substrate geometry
    INHIBITOR_ANGLE_MAX: float = 145.0   # = TIER_ANGLE_MIN['Tier_2B']:   SN2 < this → potential inhibitor
    SUBSTRATE_CONF_MIN:  float = 0.75    # Boltz confidence ≥ this → AI-confident pose

    # ── § 9.4c  Confidence-vs-tier conflict thresholds (Fig 20 / Hidden-Gem rescue) ─
    # Used by 03_Validation_Figures analyse_conflicts() to split structures into
    # Consensus High / Hidden Gem / Consensus Low / Decoy. Kept separate from
    # SUBSTRATE_CONF_MIN: this axis is AI-vs-physics agreement, not substrate geometry.
    CONFLICT_CONF_HIGH: float = 0.70   # Boltz confidence ≥ this → AI-confident (Consensus High; high-tier below = no rescue)
    CONFLICT_CONF_LOW:  float = 0.60   # Boltz confidence < this → AI-doubtful (low-tier below = Consensus Low)

    # ── § 9.5  Tier scoring & ranking weights ─────────────────────────────
    # TIER_SCORE       — additive points for composite scoring
    # TIER_RANK        — integer rank for quality-sorted operations
    # TIER_SORT_WEIGHT — weight for multi-key DataFrame sorting
    TIER_SCORE: dict = field(default_factory=lambda: {
        "Tier_1A": 20, "Tier_1B": 15, "Tier_2A": 10,
        "Tier_2B":     8, "Tier_3":       4, "Tier_4":    1, "Tier_5_Decoy": 0,
    })
    TIER_RANK: dict = field(default_factory=lambda: {
        "Tier_1A": 6, "Tier_1B": 5, "Tier_2A": 4,
        "Tier_2B":    3, "Tier_3":      2, "Tier_4":   1, "Tier_5_Decoy": 0,
    })
    TIER_SORT_WEIGHT: dict = field(default_factory=lambda: {
        "Tier_1A": 50, "Tier_1B": 40, "Tier_2A": 30,
        "Tier_2B":    20, "Tier_3":      10, "Tier_4":    1, "Tier_5_Decoy": 0,
    })

    # ── § 9.6  Tier display colours (Okabe-Ito colourblind-safe palette) ──
    TIER_COLOUR: dict = field(default_factory=lambda: {
        "Tier_1A": "#009E73",   # green
        "Tier_1B": "#56B4E9",   # sky blue
        "Tier_2A":    "#0072B2",   # blue
        "Tier_2B":    "#CC79A7",   # pink
        "Tier_3":      "#E69F00",   # orange
        "Tier_4":      "#D55E00",   # vermillion
        "Tier_5_Decoy":     "#CBD5E1",   # light grey
        "Control":   "#333333",   # dark grey
    })

    # ── § 9.7  Conflict category colours ─────────────────────────────────
    CONFLICT_COLOUR: dict = field(default_factory=lambda: {
        "Consensus High": "#0072B2",
        "Hidden Gem":     "#CC79A7",
        "Decoy":          "#D55E00",
        "Consensus Low":  "#F0E442",
        "Ambiguous":      "#BBBBBB",
    })

    # ── § 9.8  Alignment grade colours (worst I = brown → best A = green) ─
    # Keys match single-letter ALIGN_GRADE_LABELS; tuple = (lo%, hi%, hex)
    GRADE_COLOUR: dict = field(default_factory=lambda: {
        'A': '#2E7D52', 'B': '#559B6A', 'C': '#88B88A',
        'D': '#C8C87A', 'E': '#C9A84C', 'F': '#C07C40',
        'G': '#B05535', 'H': '#8B3A2A', 'I': '#6B2A1F',
    })

    # Grade bands as list of (lo, hi, hex, letter) for Fig 02 histogram shading
    GRADE_BANDS: list = field(default_factory=lambda: [
        (90, 100, '#2E7D52', 'A'), (80, 90, '#559B6A', 'B'),
        (70, 80, '#88B88A', 'C'), (60, 70, '#C8C87A', 'D'),
        (50, 60, '#C9A84C', 'E'), (40, 50, '#C07C40', 'F'),
        (30, 40, '#B05535', 'G'), (20, 30, '#8B3A2A', 'H'),
        (0,  20, '#6B2A1F', 'I'),
    ])

    # Grade colours keyed by full label string (used in Fig 03 stacked bars)
    GRADE_COLOUR_FULL: dict = field(default_factory=lambda: {
        'I (<20%)':   '#6B2A1F', 'H (20–30%)': '#8B3A2A',
        'G (30–40%)': '#B05535', 'F (40–50%)': '#C07C40',
        'E (50–60%)': '#C9A84C', 'D (60–70%)': '#C8C87A',
        'C (70–80%)': '#88B88A', 'B (80–90%)': '#559B6A',
        'A (≥90%)':   '#2E7D52',
    })

    # ── § 9.9  Mechanistic outcome colours (Fig 25b) ─────────────────────
    # Deliberately a different hue family from the tier palette so the two
    # stacked bars in Fig 25b never read as the same encoding.
    OUTCOME_COLOUR: dict = field(default_factory=lambda: {
        'Substrate':           '#0E7C7B',   # teal
        'Borderline':          '#C9A227',   # mustard
        'Reactive (low conf)': '#8E44AD',   # purple
        'Non-reactive':        '#95A5A6',   # grey
        'Potential Inhibitor': '#6E2C00',   # brown-maroon
    })

    # ── § 9.9b  PFAS chain-length bin colours (Fig 25b master bars) ───────
    PFAS_SIZE_BIN_COLOUR: list = field(default_factory=lambda:
        ['#6A51A3', '#2171B5', '#238B45', '#D94801', '#A50F15'])

    # Best → worst 5-stop gradient (nucleophile-distance & SN2-angle bins).
    SANKEY_GRAD5: list = field(default_factory=lambda:
        ['#1B7837', '#5AAE61', '#D9EF8B', '#FDAE61', '#D73027'])
    SANKEY_CLAMP_COLOUR: dict = field(default_factory=lambda: {
        'Clamp intact': '#1B9E77', 'Clamp broken': '#D55E00'})
    SANKEY_STAB_COLOUR: dict = field(default_factory=lambda: {
        'Stabilised': '#1B9E77', 'Unstabilised': '#BDBDBD'})
    SANKEY_TRIAD_COLOUR: dict = field(default_factory=lambda: {
        'Triad tight': '#1B9E77', 'Triad moderate': '#E6A817', 'Triad loose': '#D55E00'})
    SANKEY_MECH_GRAD: list = field(default_factory=lambda:
        ['#D73027', '#FDAE61', '#A6D96A', '#1B7837'])   # worst→best, 4 bins

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 10 ── MD TRAJECTORY ANALYSIS  (Step 08)
    # ═════════════════════════════════════════════════════════════════════════════

    # ── § 10.1  Solvent / water residue names ────────────────────────────
    SOLVENT_RESTYPES: tuple = (
        "SPC", "TIP3P", "TIP4P", "SPCE", "OPC", "HOH", "WAT", "SOL",
        "T3P",   # Desmond internal name for TIP3P
    )

    # ── § 10.2  WaterMap integration ─────────────────────────────────────
    # Abel et al. (2008): site radius 5.0 Å captures active-site hydration shell.
    WATERMAP_SITE_RADIUS: float     = 5.0   # Å  radius to include WaterMap hydration sites
    WATERMAP_BLOCKADE_RADIUS: float = 3.0   # Å  cylinder radius for nucleophile-runway blockade check
    WATERMAP_MATCH_RADIUS: float    = 1.5   # Å  MD water ↔ WaterMap site assignment distance

    # ── § 10.3  Frame scoring weights (QM/MM frame selection only) ────────
    # These weights rank MD frames to select the single best input for QSite.
    # They do NOT affect NAC counts or catalytic viability percentages.
    SCORE_DIST_WEIGHT: float     = 100.0   # per Å below NAC_DIST_RELAXED
    SCORE_ANGLE_WEIGHT: float    =   5.0   # per degree above NAC_ANGLE_RELAXED
    SCORE_BLOCKADE_WEIGHT: float =  50.0   # per water blockade unit
    SCORE_WATERMAP_WEIGHT: float =  10.0   # per kcal mol⁻¹ WaterMap dG unit

    # ── § 10.4  Catalytic viability display thresholds (%) ───────────────
    # Used to colour-code console output and figures in Step 08.
    # Boltz-2 predicted structures (not crystal structures) typically show
    # lower NAC populations (0.01–2%) owing to the Boltz-2 starting geometry
    # not being pre-optimised for the reactive SN2 trajectory.
    VIABILITY_PASS_THRESHOLD: float = 0.05   # % — below this → "NAC FAIL" (red)
    VIABILITY_HIGH_THRESHOLD: float = 5.0    # % — above this → "NAC PASS – High" (green)

    # ── § 10.5  Force field — Desmond MD ────────────────────────────────────
    # System solvation, ion neutralisation, and NPT production MD all use the
    # OPLS4 force field (Roos et al. 2019; Lu et al. 2021, JCTC) as configured
    # in the Desmond .msj job file generated by the user via Maestro.
    # The Python script (08_MD_Thermodynamics_QMMM_Engine_FAcDs.py) reads the
    # completed trajectory — it does not control force-field selection.
    # OPLS4 reference: Lu et al. (2021) J Chem Theory Comput 17:4291–4300.
    #                  DOI: https://doi.org/10.1021/acs.jctc.1c00302

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 11 ── QM/MM EXTRACTION — QSite  (Step 08)
    # Level of theory: M06-2X / 6-31+G(d,p) — Zhao & Truhlar (2008) DFT
    # functional; Rosta et al. (2006) QM/MM free-energy methodology;
    # Murphy et al. (2000) QSite implementation.
    # ═════════════════════════════════════════════════════════════════════════════
    QSITE_FUNCTIONAL: str   = "M062X"        # DFT functional
    QSITE_BASIS_SET: str    = "6-31+G(d,p)"  # basis set
    QSITE_CHARGE: int       = -1             # default QM region charge (anionic carboxylate/sulfonate PFAS)
    # Neutral ligands (alcohols, non-ionised at pH 8) override the default charge.
    # Keys are substrings of the job-name ligand suffix (case-insensitive match).
    LIGAND_QM_CHARGES: dict = field(default_factory=lambda: {
        "FTOH": 0,   # fluorotelomer alcohols (6:2-FTOH, 8:2-FTOH) — neutral at pH 8
    })
    QSITE_MULT: int         = 1              # spin multiplicity (closed-shell singlet)
    # Implicit solvation for the coordinate scan. SGB (Surface Generalised Born) is a
    # SCREENING-grade choice: fast, but it can over-stabilise the charge-separated
    # SN2 C–F cleavage transition state, biasing barrier heights. For publication-grade
    # activation energies on top candidates, re-run with explicit-solvent QM/MM-FEP or
    # thermodynamic integration rather than trusting these implicit-solvent ΔG‡ values.
    QSITE_SOLVATION: str    = "sgb"          # implicit solvation model (screening-grade; see note)
    QSITE_SCAN_START: float = 3.5            # Å  scan start (pre-reaction approach; 3.5 → 1.3 Å = full SN2 coordinate)
    QSITE_SCAN_STEP: float  = -0.1           # Å  step per point (negative = bond compression)
    QSITE_SCAN_NSTEPS: int  = 23             # points total → covers 3.5 → 1.3 Å (last point: 3.5 + −0.1×22 = 1.3 Å)

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 12 ── CANONICAL RESIDUE MAPPING — 3R3U Reference  (Step 08)
    # Reference sequence positions in 3R3U / DEHA4_CONTROL_SEQ numbering.
    # Step 08 looks these up as keys in aln_dict (which is keyed by 3R3U
    # reference positions from the Full_Sequence_Alignment_Map column).
    # Must stay in sync with REF_ACTIVE_SITE_MAP (§3.5).
    # ═════════════════════════════════════════════════════════════════════════════
    DREAM_TEAM_REFS: dict = field(default_factory=lambda: {
        "Nuc":    110,   # Asp110 — nucleophile; forms covalent ester intermediate
        "Clamp1": 111,   # Arg111 — carboxylate clamp 1
        "Clamp2": 114,   # Arg114 — carboxylate clamp 2
        "Acid":   134,   # Asp134 — catalytic acid; protonates departing fluoride
        "Stab_H": 155,   # His155 — fluoride stabiliser
        "Stab_W": 156,   # Trp156 — fluoride cradle
        "Stab_Y": 217,   # Tyr217 — fluoride cradle (PDB 219)
        "Base":   277,   # His277 — base catalyst (PDB 280)
    })

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 13 ── PROCESSING PARAMETERS
    # ═════════════════════════════════════════════════════════════════════════════

    # ── § 13.1  File I/O & batching (Step 05) ────────────────────────────
    PROC_BATCH_SIZE: int         = 2000   # max futures submitted at once — caps peak memory
    PROC_HEADER_CHECK_LINES: int =   30   # lines to scan when reading a PDB header tag

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 14 ── VISUALISATION PARAMETERS
    # ═════════════════════════════════════════════════════════════════════════════

    # ── § 14.1  Global rendering (Step 08) ────────────────────────────────
    VIS_IMG_WIDTH: int   = 2400   # px  export width  (publication-quality figure)
    VIS_IMG_HEIGHT: int  = 2400   # px  export height
    VIS_RAY_TRACE: bool  = True   # enable PyMOL ray-tracing for publication quality
    VIS_FIGURE_DPI: int  = 300    # dots per inch for publication figures (minimum 300)

    # ── § 14.2  Per-structure rendering timeouts — seconds (Step 08) ─────
    VIS_TIMEOUT_PYMOL: int    = 600   # PyMOL render timeout (2× 2400-px ray traces ~5 min on CPU)
    VIS_TIMEOUT_CHIMERAX: int = 180   # ChimeraX render timeout
    VIS_TIMEOUT_MAESTRO: int  = 300   # Maestro render timeout
    VIS_TIMEOUT_PLIP: int     =  90   # PLIP interaction analysis timeout
    VIS_TIMEOUT_LIGPLOT: int     =  60   # LigPlot+ render timeout
    VIS_TIMEOUT_PYMOL_HEAVY: int = 1200  # s  PyMOL timeout for structures with ≥4 C–F bonds (polyfluorinated PFAS)

    # ── § 14.3  Geometric highlight radii (Step 08) ───────────────────────
    VIS_F_CONTACT_RADIUS: float = 4.0   # Å  fluorine contact highlight sphere
    VIS_POCKET_RADIUS: float    = 5.5   # Å  binding-pocket cartoon / surface shell

    # ── § 14.4  Tier summary figure layout (Step 03) ──────────────────────
    # Pixel dimensions and typography for the per-tier stacked-bar / star-plot
    # panels produced by 03_Validation_Figures_FAcDs.py.
    VIS_TT_STAR_SIZE: int    = 460    # px  star marker diameter
    VIS_TT_SHRINK_BORDER: int =  12   # px  border shrink for tight layout
    VIS_TT_IMG_PX: int       = 800    # px  panel image width
    VIS_TT_TEXT_PX: int      = 220    # px  text annotation column width
    VIS_TT_BORDER_PX: int    =  16    # px  outer border thickness
    VIS_TT_PAD: int          =  28    # px  inter-panel padding
    VIS_TT_FONT_SIZE: int    =  44    # pt  annotation font size
    VIS_TT_AX_WIDTH: float   =   0.12 # fraction of figure width for axis panel

    # ── § 14.5  Maximum display ranks (Step 03) ────────────────────────────
    VIS_MAX_RANKS_DISPLAY: int = 50   # top-N entries shown in ranked output plots
    VIS_TOP_N_AUTO_SELECT_TIMEOUT: int = 30  # s  auto-select countdown in Step 06 interactive tier selection
    VIS_FALLBACK_WEIGHTS: dict[str, float] = field(default_factory=lambda: {
        "Boltz_Model_Confidence": 0.45,
        "iptm":                   0.25,
        "Interaction_Density_Norm": 0.20,
        "mean_plddt":             0.10
    })
    VIS_PHYLO_COLUMN_MAP: dict[str, list[str]] = field(default_factory=lambda: {
        "Protein_Name": ["Protein_Name", "protein", "protein_id"],
        "Ligand_Name":  ["Ligand_Name",  "ligand",  "Ligand"],
        "Tier":         ["degrader_tier", "Degrader_Tier", "Tier"],
        "Score":        ["Binding_Probability", "ActiveSite_Conservation_Score", "binding_likelihood_computed"]
    })

    # ═════════════════════════════════════════════════════════════════════════════
    # SECTION 15 ── SEQUENCE ALIGNMENT & SCORING PARAMETERS  (Step 02)
    # ═════════════════════════════════════════════════════════════════════════════

    # ── § 15.1  Pairwise sequence alignment ──────────────────────────────
    # Biopython PairwiseAligner configured for global BLOSUM62 alignment.
    # Ref: Henikoff & Henikoff (1992) PNAS 89:10915–10919. DOI: 10.1073/pnas.89.22.10915 (BLOSUM62 substitution matrix)
    ALIGN_OPEN_GAP_SCORE: float   = -10.0  # open-gap penalty (full gap initiation cost)
    ALIGN_EXTEND_GAP_SCORE: float = -0.5   # gap-extension penalty (per-residue gap cost)

    # ── § 15.1b  Alignment grade bins (Step 02) ──────────────────────────
    # Used with pd.cut() to assign letter grades to Alignment_Score_Pct.
    # bins[i] < score ≤ bins[i+1] → label[i].  right=True (default for pd.cut).
    ALIGN_GRADE_BINS:   list = field(default_factory=lambda:
                              [0, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    ALIGN_GRADE_LABELS: list = field(default_factory=lambda:
                              ['I', 'H', 'G', 'F', 'E', 'D', 'C', 'B', 'A'])

    # ── § 15.2  Dynamic alignment cost-matrix scaling ────────────────────
    ALIGN_DYNAMIC_PENALTY_FACTOR: float   = 2.0  # scale factor applied to max(cost) for dynamic penalty
    ALIGN_DYNAMIC_PENALTY_FALLBACK: float = 5.0  # penalty value when cost matrix is empty

    # ── § 15.3  Likelihood scoring thresholds (Step 02) ──────────────────
    # Applied during Likelihood_Degrader_Score computation.
    SCORE_DESOLVATION_THRESHOLD: float        = 0.5   # hydrophobic desolvation ratio above which bonus applies
    SCORE_DESOLVATION_BONUS_FACTOR: float     = 10.0  # multiplier for desolvation bonus (× ratio)
    SCORE_INHIBITION_DENSITY_THRESHOLD: float = 1.5   # interaction density above which product-inhibition penalty fires
    SCORE_INHIBITION_PENALTY_FACTOR: float    = 15.0  # scaling factor for product-inhibition penalty

    # ── § 15.4  GPU batch watchdog ───────────────────────────────────────
    GPU_WATCHDOG_TIMEOUT_PER_JOB: int = 600  # seconds per job before GPU-batch watchdog kills the run

    # ── § 15.5  Affinity interaction contribution weights (Step 02) ──────────
    # custom_affinity_score = Σ w·interaction_count - w_dist·avg_dist
    SCORE_W_AFF_HB:   float = 1.5    # hydrogen bond
    SCORE_W_AFF_HP:   float = 1.0    # hydrophobic contact
    SCORE_W_AFF_SB:   float = 2.5    # salt bridge
    SCORE_W_AFF_HAL:  float = 2.0    # halogen bond
    SCORE_W_AFF_FL:   float = 1.0    # fluorine contact
    SCORE_W_AFF_FP:   float = 4.0    # fluorine polar (F···H-N/O)
    SCORE_W_AFF_FH:   float = 2.0    # fluorous–hydrophobic
    SCORE_W_AFF_FF:   float = 3.0    # fluorous–fluorous
    SCORE_W_AFF_MC:   float = 2.5    # metal coordination
    SCORE_W_AFF_DIST: float = 0.05   # average distance penalty

    # ═════════════════════════════════════════════════════════════════════════════

    # SECTION 16 ── PDB PREPARATION — Schrödinger PrepWizard  (Step 05)
    # ═════════════════════════════════════════════════════════════════════════════

    # ── § 16.1  PrepWizard protonation & minimisation ────────────────────
    # Ref: Madhavi Sastry et al. (2013) J Comput-Aided Mol Des 27:221–234. DOI: 10.1007/s10822-013-9644-8 (PrepWizard protocol)
    PREPWIZARD_PROPKA_PH: float      = 8.0  # protein protonation pH (FAcD physiological context)
    PREPWIZARD_EPIK_PH: float        = 8.0  # PFAS ligand protonation pH via Epik
    PREPWIZARD_RMSD_RESTRAIN: float  = 0.3  # Å — RMSD restraint for clash-resolving minimisation

    # ── § 16.2  Parallelism ───────────────────────────────────────────────
    PREP_CPU_RESERVE: int  = 2  # CPU cores to reserve for OS/desktop stability (not used by PrepWizard)

    # ── § 16.3  Chain assignment ──────────────────────────────────────────
    PREP_LIGAND_CHAIN: str = "L"  # chain identifier for non-protein (ligand) residues in converted PDB

    # ── § 16.4  Standard amino-acid residue set (used for ligand filtering) ─
    PREP_STANDARD_AA: set = field(default_factory=lambda: {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    })
    PREP_AMBIGUOUS_AA: set = field(default_factory=lambda: set("BXZJOU"))

    # ── § 16.5  Protein-associated residues (stay in protein chain) ───────
    # Modified amino acids and structural metals that must NOT be moved to Chain L.
    PREP_PROTEIN_ASSOCIATED: set = field(default_factory=lambda: {
        # Modified AAs
        "MSE", "SEP", "TPO", "PTR", "KCX", "CME", "CSO", "PCA", "LLP",
        # Structural / catalytic metals
        "ZN", "MG", "CA", "FE", "MN", "CO", "NI", "CU", "NA", "K",
    })

    # ===============================================================================
    # SECTION 17 ── DATA REGISTRY & VISUAL AESTHETICS
    # ===============================================================================

    # ── § 17.1  CSV Column Names ───────────────────────────────────────────
    # Registry for all master CSV columns to prevent hardcoding string keys.
    COL_TIER:     str = "degrader_tier"
    COL_PROT:     str = "Protein_Name"
    COL_LIG:      str = "Ligand_Name"
    COL_CONF:     str = "confidence_score"
    COL_PTM:      str = "ptm"
    COL_IPTM:     str = "iptm"
    COL_IDENS:    str = "interaction_density"
    COL_ID_PCT:   str = "identity_pct"
    COL_TRUST:    str = "Trust_Score"
    COL_SN2:      str = "sn2_attack_angle"
    COL_NUC_DIST: str = "NAC_dist_Nuc"
    COL_MECH_S:   str = "Mechanistic_Fingerprint_Score"
    COL_LIKE_S:   str = "ActiveSite_Conservation_Score"
    COL_ALN_G:    str = "Alignment_Grade"

    # ── § 17.2  Visual Plotting Properties ──────────────────────────────────
    # Canonical marker sizes (S) and opacities (A) for each tier.
    # Tier_1/Tier_2 tiers are drawn larger and more opaque.
    @property
    def VIS_TIER_SIZES(self) -> dict:
        return {
            self.TIER_ORDER[0]: 420, self.TIER_ORDER[1]: 45, self.TIER_ORDER[2]: 16,
            self.TIER_ORDER[3]: 9, self.TIER_ORDER[4]: 5, self.TIER_ORDER[5]: 3, self.TIER_ORDER[6]: 2
        }

    @property
    def VIS_TIER_ALPHAS(self) -> dict:
        return {
            self.TIER_ORDER[0]: 1.0, self.TIER_ORDER[1]: 0.88, self.TIER_ORDER[2]: 0.70,
            self.TIER_ORDER[3]: 0.38, self.TIER_ORDER[4]: 0.22, self.TIER_ORDER[5]: 0.18, self.TIER_ORDER[6]: 0.10
        }

    # ── § 17.3  Alignment Grade Bins ────────────────────────────────────────
    # Define identity % ranges for each grade letter.
    ALIGN_GRADE_DEFS: list = field(default_factory=lambda: [
        (90.0, 100.0, 'A'), (80.0, 90.0, 'B'), (70.0, 80.0, 'C'),
        (60.0, 70.0,  'D'), (50.0, 60.0, 'E'), (40.0, 50.0, 'F'),
        (30.0, 40.0,  'G'), (20.0, 30.0, 'H'), (0.0,  20.0, 'I')
    ])

    # ── § 9.8  Tier Meanings ──────────────────────────────────────────────────
    TIER_MEANING: dict = field(default_factory=lambda: {
        "Tier_1A": "Elite / Perfect_A (legacy)",
        "Tier_1B": "Excellent / Perfect_B (legacy)",
        "Tier_2A": "Strong / Best_A (legacy)",
        "Tier_2B": "Good / Best_B (legacy)",
        "Tier_3":  "Moderate / Good (legacy)",
        "Tier_4":  "Poor / Poor (legacy)",
        "Tier_5_Decoy": "Invalid / Decoy (legacy)",
        "Control": "Reference Control"
    })
