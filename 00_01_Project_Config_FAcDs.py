"""
===============================================================================
FAcDs Pipeline  |  Module 00_01  |  Central Configuration (CFG)
===============================================================================
Single source of truth for every threshold, constant, weight, and parameter
used across the pipeline. Edit values here only — no other file should contain
hard-coded scientific values, configurable thresholds, or tunable settings.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 10 July 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Module        : 00_01_Project_Config_FAcDs.py
  Role          : "Blueprint" — pipeline-wide configuration repository.
  Imported by   : All pipeline scripts (00_02 through 07).
  Reads         : Nothing (pure Python dataclass; no file I/O).
  Writes        : Nothing.
  Upstream      : None (root module — must load before all others).
  Downstream    : 02_Production_FAcDs.py, 03_Validation_Figures_FAcDs.py,
                  04_Dendrogram_FAcDs.py, 05_TopN_and_PDB_Preparation_FAcDs.py,
                  05_TopN_and_PDB_Preparation_FAcDs.py,
                  07_MD_QMMM_Defluorination_FAcDs.py
───────────────────────────────────────────────────────────────────────────────

── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Module Load Priority: This file MUST be imported before any other project
     modules to ensure CFG is available for dependent logic.
  2. Static Constancy: Changing thresholds (e.g., NAC_DIST) mid-pipeline will
     cause inconsistencies between stored CSV results and new analysis.
  3. Environment Sensitivity: Relies on `importlib.util` for path-safe imports;
     requires consistent relative directory structure.
───────────────────────────────────────────────────────────────────────────────

Usage:
    import importlib.util as _ilu
    _spec = _ilu.spec_from_file_location("ProjectConfig", path / "00_01_Project_Config_FAcDs.py")
    _mod  = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_mod)
    CFG   = _mod.CFG()   # instantiate once at module level

    CFG.NAC_DIST_STRICT   # access any attribute directly
    CFG.THRESHOLD_HB_DIST_MAX

Section map (prefix → section):
    BOLTZ_*        §1     Boltz-2 prediction engine & ColabFold MSA
    THRESHOLD_*    §3     Biochemical interaction geometry cutoffs
    NAC_*          §4     Near Attack Conformation geometry
    SMART_LOCK_*   §7     3D Smart-Lock residue-detection bias
    TIER_*         §8     Catalytic tier classification & display
    WATERMAP_*     §9     WaterMap hydration-site integration
    VIABILITY_*    §9     Catalytic viability display thresholds
    SCORE_*        §9/14  Likelihood scoring thresholds
    QSITE_*        §10    QM/MM extraction settings (QSite/Schrödinger)
    PROC_*         §12    Processing / file-I/O operational parameters
    VIS_*          §13    Visualisation rendering parameters
    ALIGN_*        §14    Sequence alignment penalty parameters
    GPU_*          §14    GPU batch watchdog timeout
    PREPWIZARD_*   §15    Schrödinger PrepWizard invocation parameters
    PREP_*         §15    PDB preparation chain & residue classification
    COL_*          §16    Master-CSV column-name registry
    MMGBSA_*       §17    Prime MM-GBSA binding free-energy parameters
    MD_*           §18    MD-ready cohort selection (gates Steps 05-07)
───────────────────────────────────────────────────────────────────────────────

-------------------------------------------------------------------------------
Scientific References:
    1. Structure prediction (Boltz-2):
       - Passaro, S. et al. (2025) bioRxiv 2025.06.14.659707.
       - DOI: https://doi.org/10.1101/2025.06.14.659707
    2. MSA generation (ColabFold server):
       - Mirdita, M. et al. (2022) Nature Methods 19:679–682.
       - DOI: https://doi.org/10.1038/s41592-022-01488-1
    3. FAcD crystal structure & catalytic mechanism (PDB 3R3U):
       - Chan, P.W.Y., Yakunin, A.F., Edwards, E.A. & Pai, E.F. (2011) JACS 133:7461–7468.
       - DOI: https://doi.org/10.1021/ja200277d
       - FAcD small-substrate scope + TFA recalcitrance (graded chemistry/containment penalties):
         Wackett, L.P. (2022) Microb Biotechnol 15(3):773–792. DOI: https://doi.org/10.1111/1751-7915.13928
    4. DEHA4 experimental defluorination (Delftia acidovorans D4B):
       - Farajollahi, S. et al. (2024) ACS Omega 9(26):28546–28555.
       - DOI: https://doi.org/10.1021/acsomega.4c02517
       - Jitsumori, K. et al. (2009) J Bacteriol 191:2630–2637.
       - DOI: https://doi.org/10.1128/JB.01654-08
    5. Interaction-geometry cutoffs — operational source (§3):
       Every interaction distance/angle cutoff in §3 reproduces the Schrödinger
       Maestro default interaction criteria (H-bonds, halogen bonds, salt bridges,
       aromatic H-bonds, π–π, π–cation, steric contacts), encoded as explicit
       constants so detection runs WITHOUT a Schrödinger installation (no
       $SCHRODINGER binary is called). Software default + primary literature:
       - Schrödinger Release 2026-1: Maestro, Schrödinger, LLC, New York, NY. https://www.schrodinger.com/maestro
       - H-bond geometry (H···A, angle): McDonald & Thornton (1994) J Mol Biol 238:777–793. DOI: https://doi.org/10.1006/jmbi.1994.1334
       - H-bond heavy-atom proxy (D···A ≈ H···A + ~1.0 Å, for H-free Boltz-2 CIF): Jeffrey, G.A. (1997) An Introduction to Hydrogen Bonding. Oxford University Press. ISBN 978-0-19-509549-4.
       - Salt bridge (≤5.0 Å): Barlow & Thornton (1983) J Mol Biol 168:867–885, DOI: https://doi.org/10.1016/S0022-2836(83)80079-5; Kumar & Nussinov (2002) ChemBioChem 3:604–617. DOI: https://doi.org/10.1002/1439-7633(20020703)3:7<604::AID-CBIC604>3.0.CO;2-X
       - Hydrophobic contact: Salentin, S. et al. (2015) Nucleic Acids Res 43:W443–W447. DOI: https://doi.org/10.1093/nar/gkv315
       - Aromatic (weak) H-bond: Levitt & Perutz (1988) J Mol Biol 201:751–754. DOI: https://doi.org/10.1016/0022-2836(88)90471-8
       - π–π stacking: McGaughey, G.B. et al. (1998) J Biol Chem 273:15458–15463. DOI: https://doi.org/10.1074/jbc.273.25.15458
       - π–cation (≤6.6 Å): Gallivan & Dougherty (1999) PNAS 96:9459–9464. DOI: https://doi.org/10.1073/pnas.96.17.9459
       - Halogen bond: Auffinger, P. et al. (2004) PNAS 101:16789–16794 (DOI: https://doi.org/10.1073/pnas.0407607101); Wilcken, R. et al. (2013) J Med Chem 56:1363–1388 (DOI: https://doi.org/10.1021/jm3012068).
       - Steric contact ratio (distance ÷ Σ vdW): vdW radii Bondi, A. (1964) J Phys Chem 68:441–451. DOI: https://doi.org/10.1021/j100785a001
       - Metal coordination (≤2.8 Å): Harding, M.M. (2006) Acta Crystallogr D62:678–682. DOI: https://doi.org/10.1107/S0907444906014594
    6. Near Attack Conformation (NAC) geometry:
       - Lightstone, F.C. & Bruice, T.C. (1996) JACS 118:2595–2605. DOI: https://doi.org/10.1021/ja952589l
       - Bruice, T.C. (2002) Acc Chem Res 35:139–148. DOI: https://doi.org/10.1021/ar0001665
       - Hur, S. & Bruice, T.C. (2003) PNAS 100:12015–12020. DOI: https://doi.org/10.1073/pnas.1534873100
    7. Catalytic triad distances:
       - Holmquist, M. (2000) Curr Protein Pept Sci 1:209–235. DOI: https://doi.org/10.2174/1389203003381405
    8. Bürgi–Dunitz trajectory (auxiliary, non-gating):
       - Bürgi, H.B., Dunitz, J.D. & Shefter, E. (1973) JACS 95:5065–5067. DOI: https://doi.org/10.1021/ja00796a058
       - Bürgi, H.B., Dunitz, J.D., Lehn, J.M. & Wipff, G. (1974) Tetrahedron 30:1563–1572. DOI: https://doi.org/10.1016/S0040-4020(01)90678-7
    9. SN2 dead-end consensus check (scissile C–F BDE + backside sterics):
       - O'Hagan, D. (2008) Chem Soc Rev 37:308–319. DOI: https://doi.org/10.1039/B711844A  (C–F bond strength)
       - Bento, A.P. & Bickelhaupt, F.M. (2008) J Org Chem 73:7290–7299. DOI: https://doi.org/10.1021/jo801215z  (backside SN2)
       - Bondi, A. (1964) J Phys Chem 68:441–451. DOI: https://doi.org/10.1021/j100785a001  (vdW radii)
   10. gem-Difluoroacetate is a genuine FAcD substrate (no CF2 penalty):
       - Khusnutdinova, A.N. et al. (2023) FEBS J 290:4966–4983. DOI: https://doi.org/10.1111/febs.16903
   11. WaterMap thermodynamic hydration-site analysis:
       - Abel, R., Young, T., Farid, R., Berne, B.J. & Friesner, R.A. (2008) JACS 130:2817–2831. DOI: https://doi.org/10.1021/ja0771033
   12. QSite QM/MM (B3LYP/6-31+G(d,p)) defluorination energetics:
       - Becke, A.D. (1993) J Chem Phys 98:5648–5652. DOI: https://doi.org/10.1063/1.464913
       - Lee, C., Yang, W. & Parr, R.G. (1988) Phys Rev B 37:785–789. DOI: https://doi.org/10.1103/PhysRevB.37.785
       - Murphy, R.B. et al. (2000) J Comput Chem 21:1442–1457. DOI: https://doi.org/10.1002/1096-987X(200012)21:16<1442::AID-JCC3>3.0.CO;2-O  (QSite implementation)
       - Rosta, E. et al. (2006) J Phys Chem B 110:2934–2941. DOI: https://doi.org/10.1021/jp057109j  (QM/MM free-energy benchmark)
       - Yue, Y. et al. (2021) Environ Sci Technol 55(14):9817–9825. DOI: https://doi.org/10.1021/acs.est.0c08811  (FAcD QM/MM)
   13. OPLS4 force field (MD parameterisation):
       - Lu, C. et al. (2021) J Chem Theory Comput 17:4291–4300. DOI: https://doi.org/10.1021/acs.jctc.1c00302
       - Roos, K. et al. (2019) J Chem Theory Comput 15:1863–1874. DOI: https://doi.org/10.1021/acs.jctc.8b01026
   14. Sequence alignment (Biopython PairwiseAligner, BLOSUM62 substitution matrix):
       - Cock, P.J.A. et al. (2009) Bioinformatics 25:1422–1423. DOI: https://doi.org/10.1093/bioinformatics/btp163
       - Henikoff, S. & Henikoff, J.G. (1992) PNAS 89:10915–10919. DOI: https://doi.org/10.1073/pnas.89.22.10915
   15. Protein preparation (Schrödinger PrepWizard protocol):
       - Sastry, G.M., Adzhigirey, M., Day, T., Annabhimoju, R. & Sherman, W. (2013) J Comput-Aided Mol Des 27:221–234. DOI: https://doi.org/10.1007/s10822-013-9644-8
===============================================================================
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class CFG:
    """
    Central Configuration Repository — FAcDs Pipeline.

    Attribute prefix → section
    ──────────────────────────
    BOLTZ_*        §1     Boltz-2 prediction engine & ColabFold MSA
    THRESHOLD_*    §3     Biochemical interaction geometry cutoffs
    NAC_*          §4     Near Attack Conformation geometry
    SMART_LOCK_*   §7     3D Smart-Lock residue-detection bias
    TIER_*         §8     Catalytic tier classification & display
    WATERMAP_*     §9     WaterMap hydration-site integration
    VIABILITY_*    §9     Catalytic viability display thresholds
    SCORE_*        §9/14  Likelihood scoring thresholds
    QSITE_*        §10    QM/MM extraction settings (QSite/Schrödinger)
    PROC_*         §12    Processing / file-I/O operational parameters
    VIS_*          §13    Visualisation rendering parameters
    ALIGN_*        §14    Sequence alignment penalty parameters
    GPU_*          §14    GPU batch watchdog timeout
    PREPWIZARD_*   §15    Schrödinger PrepWizard invocation parameters
    PREP_*         §15    PDB preparation chain & residue classification
    COL_*          §16    Master-CSV column-name registry
    MMGBSA_*       §17    Prime MM-GBSA binding free-energy parameters
    MD_*           §18    MD-ready cohort selection (gates Steps 05-07)
    """

    # ===============================================================================
    # SECTION 1: PROJECT IDENTITY & BOLTZ-2 PREDICTION ENGINE  (Step 02)
    # ===============================================================================
    PROJECT_NAME: str   = "PFAS-27"   # project identifier (used in CLI banners, e.g. Step 06)
    """
    Pipeline version — the single source of truth every step reports and every run stamps into its
    outputs. A result is only reproducible if the code that produced it can be named, and a bare
    timestamp cannot do that: two runs on the same day can straddle a change to the tier gates.
    Bump the MINOR component when a change moves the numbers (a gate, a weight, a metric definition)
    and the PATCH component for anything that cannot.
    """
    PIPELINE_VERSION: str = "2.1.0"   # 2.1.0 — protein-aware pocket containment, attacking-oxygen NAC, BDE elite ceiling, cradle-aware Šidák, active-site pLDDT gate

    # -------------------------------------------------------------------------------
    # Step 1.1: Executable & model
    # -------------------------------------------------------------------------------
    BOLTZ_EXECUTABLE: str          = "boltz"    # CLI binary name on PATH
    BOLTZ_MODEL_VERSION: str       = "boltz2"   # --model flag passed to boltz
    BOLTZ_OUTPUT_FORMAT: str       = "mmcif"    # --output_format flag

    # -------------------------------------------------------------------------------
    # Step 1.2: Sampling & batching
    # -------------------------------------------------------------------------------
    BOLTZ_RECYCLING_STEPS: int        = 1    # --recycling_steps per prediction
    BOLTZ_DIFFUSION_SAMPLES: int      = 5    # --diffusion_samples (structures per complex)
    BOLTZ_CACHE_DIR: str              = "~/.cache/boltz"   # Boltz-2 model cache (~ expanded at use); point to shared/scratch storage if needed
    BOLTZ_MAX_PROTEINS_PER_BATCH: int = 20   # proteins bundled in one boltz call
                                             # e.g. 20 proteins × 27 ligands = 540 jobs/batch

    # -------------------------------------------------------------------------------
    # Step 1.3: Retry & polling
    # -------------------------------------------------------------------------------
    BOLTZ_RETRY_MAX: int     = 2   # max retries after GPU failure before abandoning job
    BOLTZ_RETRY_SLEEP: int   = 6   # seconds to wait between retry attempts
    """
    Boltz-2 diffusion seed. The CLI default is None — NO seeding — so every prediction is stochastic
    and the same protein-ligand pair yields a different pose on every run. Measured on the DeHa4
    fluoroacetate control across two runs of identical code: the SN2 attack angle moved 145.5° → 132.4°
    and the nucleophile distance 3.77 → 3.03 Å. A tier assignment that changes when nothing changed
    is not a result, and it silently invalidates any comparison between runs — including the control
    read-out, which is the pipeline's own calibration check.

    Seeding costs nothing and makes a run reproducible from its config alone. It does NOT reduce
    conformational sampling: diffusion_samples still draws the full ensemble per complex; the seed
    only fixes where that ensemble starts, so the same command returns the same ensemble.
    """
    BOLTZ_SEED: int          = 42   # --seed passed to every Boltz prediction; fixes run-to-run drift
    CONTROL_RESIDUE_MATCH_RADIUS: float = 6.0   # Å — search radius when matching a control's catalytic residue onto its structural counterpart (Cα-anchored nearest same-type residue)
    BOLTZ_PREDICT_TIMEOUT_S: int = 3600  # per-prediction wall-clock ceiling; a single co-fold never approaches this, so a breach means a frozen GPU/CUDA driver → kill and retry rather than stall the pipeline
    BOLTZ_POLL_INTERVAL: int = 5   # seconds between job-status poll cycles

    # -------------------------------------------------------------------------------
    # Step 1.4: ColabFold MSA submission
    # -------------------------------------------------------------------------------
    COLABFOLD_API_URL: str        = "https://api.colabfold.com"
    COLABFOLD_SUBMIT_RETRIES: int = 8    # max attempts to POST an MSA ticket
    COLABFOLD_MSA_SUBMIT_TIMEOUT: int       = 30   # seconds — HTTP POST timeout for ticket submission
    COLABFOLD_MSA_POLL_TIMEOUT: int         = 900  # seconds — max wait for MSA completion (15 min)
    COLABFOLD_MSA_POLL_REQUEST_TIMEOUT: int = 15   # seconds — HTTP GET timeout per polling cycle
    COLABFOLD_MSA_DOWNLOAD_TIMEOUT: int     = 120  # seconds — HTTP GET timeout for final tar.gz download

    # -------------------------------------------------------------------------------
    # Step 1.5: Confidence scoring weights
    # -------------------------------------------------------------------------------
    """
    Weights for the composite Boltz confidence score used to rank models.
    Higher weight = stronger contribution to final ranking.
    """
    BOLTZ_W_IPTM: float         = 3.0   # interface predicted TM-score (primary signal)
    BOLTZ_W_CROSS_PAE: float    = 2.5   # cross-interface predicted aligned error
    BOLTZ_W_PLDDT: float        = 2.0   # per-residue local confidence
    BOLTZ_W_INTERACTIONS: float = 1.5   # density of physical interactions
    BOLTZ_W_CONF: float         = 1.0   # global structural confidence (ptm)

    # Binding-probability logit normalisers (Step 02 binding_likelihood_computed).
    BIND_INT_DENSITY_NORM: float = 2.0   # interaction-density saturation normaliser in the logit
    BIND_CROSS_PAE_NORM: float   = 50.0  # cross-interface PAE normaliser in the logit
    BIND_LOGIT_CLAMP: float      = 50.0  # ± clamp on the logit before the sigmoid

    # ===============================================================================
    # SECTION 2: INPUT & REFERENCE DATA
    # ===============================================================================

    # -------------------------------------------------------------------------------
    # Step 2.1: Reference ligand
    # -------------------------------------------------------------------------------
    # Fluoroacetate (FA) is the canonical FAcD substrate; benchmark for DEHA4 and 3R3U geometry calibration.
    FLUOROACETATE_SMILES: str = "C(C(=O)O)F"   # canonical SMILES for fluoroacetate
    INPUT_SMILES: str = "D_INP_PFAS-27_Ligands.smi"  # ligand SMILES panel

    # -------------------------------------------------------------------------------
    # Step 2.2: Reference protein sequence
    # -------------------------------------------------------------------------------
    """
    DEHA4 = DeHa4_[Delftia acidovorans D4B] — structural reference for all alignments.
    Source: Farajollahi et al. (2024) — see header Scientific References §4.
    Position cross-validated; establishes HIS277 (not HIS288) in DEHA4 sequences.
    """
    DEHA4_CONTROL_SEQ: str = (
        "MHTDPWMPGLRQQRITVDDGVEINAWVGGQGPALLLVHGHPQTSAIWHRVAPRLAQQFTVVLADLRGYGDSSRPAGDPEH"
        "VNYSKRTMARDLLRLMARLGHEHFSVLAHDRGARVAHRLAMDYPASVQRLVLLDIAPTLAMYEQTGEAFARAYWHWFFLI"
        "QPAPLPERLIEADPAAYVREIMGRRSAGLAPFDPRALAEYQRCLALPGSAHGMCEDYRASAGIDLDHDREDRQLGRRLSM"
        "PLLVLWGEEGWHRCFDPLREWQLVADDVRGRPLACGHYIAEEAPDALLDAALPFLLQAG"
    )
    DEHA4_REF_FASTA: str = "DeHa4_Ref.fasta"       # per-run reference FASTA filename (SSOT: 02 writes)

    # -------------------------------------------------------------------------------
    # Step 2.3: Crystal structure reference — PDB 3R3U
    # -------------------------------------------------------------------------------
    """
    Rhodopseudomonas palustris FAcD, wild-type, 1.60 Å resolution
    (no substrate bound; structure carries Ni²⁺ and Cl⁻ ions).
    Ref: Chan et al. (2011); Jitsumori et al. (2009) — see header Scientific References §3, §4.
    """
    REFERENCE_PDB_ID: str    = "3R3U"
    REFERENCE_PDB_URL: str   = "https://files.rcsb.org/download/3R3U.pdb"
    REFERENCE_PDB_FILE: str  = "3R3U.pdb"          # local copy filename (SSOT: 03 crystal read)
    REFERENCE_DIR_NAME: str  = "0_Reference_Crystal_3R3U"

    # -------------------------------------------------------------------------------
    # Step 2.4: Control ligand panel (fluoroacetate control substrates)
    # -------------------------------------------------------------------------------
    """
    Three short-chain fluorinated substrates used for DEHA4 and 3R3U calibration runs.
    Format: (job_suffix, SMILES)
    """
    CTRL_LIGANDS: list = field(default_factory=lambda: [
        ("26_Fluoroacetate",   "C(C(=O)O)F"),         # fluoroacetate (canonical substrate)
        ("27_Difluoroacetate", "O=C(O)C(F)F"),        # difluoroacetate
        ("25_TFA",             "O=C(O)C(F)(F)F"),    # trifluoroacetate
    ])

    # -------------------------------------------------------------------------------
    # Step 2.5: Reference active-site mapping (3R3U / DEHA4 canonical)
    # -------------------------------------------------------------------------------
    """
    Residue identities, sequential IDs, and PDB-numbered IDs for the eight key
    catalytic / stabilising positions in the FAcD (RPA1163) active site. Every
    identity and role below is taken from the primary literature — no inferred or
    unreferenced residues.

    Residue identities and roles are taken from the WILD-TYPE structure and the
    article text — never from the catalytically-dead mutant PDBs (3R3V/3R3W =
    Asp110Asn, 3R3Y/3R41 = His280Asn; those carry Asn, not the WT residue, and
    merely confirm the role by knockout). Positions are then verified against the
    DEHA4_CONTROL_SEQ (§2.2): every 'id' below holds the stated WT residue.

    Primary structural source — WT FAcD RPA1163, PDB 3R3U (3R3U PDB numbering)
    [Chan et al. (2011) — see header Scientific References §3]:
        • Asp110  nucleophile (forms the glycolyl-enzyme ester intermediate)
        • His280  general base (activates the hydrolytic water)
        • Arg111, Arg114             conserved substrate-carboxylate clamp
        • His155, Trp156, Tyr219     halide pocket — three H-bonds stabilising F⁻

    Mechanistic QM/MM source — catalytic-dyad acid + Tyr219 SN2 role
    [Yue et al. (2021) — see header Scientific References §12]:
        • Asp134  orients/polarises the His280 base (Asp–His catalytic dyad)
        • Tyr219  charge acceptor along the SN2 axis, easing electronic repulsion

    DEHA4 reference sequence (Delftia acidovorans D4B homolog; the alignment
    target for all 2150 sequences) [Farajollahi et al. (2024) — see header §4].

    Two numberings are kept because DEHA4 ≠ RPA1163 (different organism, slight
    offset past ~residue 215): 'id' = position in DEHA4_CONTROL_SEQ; 'pdb_id' =
    3R3U PDB residue number. Both verified to carry the canonical residue.
    role strings are consumed by ROLE_EXPECTED_RESIDUES (§2.7) — do not rename.
    """
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

    # -------------------------------------------------------------------------------
    # Step 2.6: Catalytic triad key set
    # -------------------------------------------------------------------------------
    CATALYTIC_TRIAD_KEYS: list = field(default_factory=lambda: ["Nuc", "Base", "Acid"])

    # -------------------------------------------------------------------------------
    # Step 2.7: Class-aware alignment search window
    # -------------------------------------------------------------------------------
    """
    When a catalytic residue does not map cleanly to the exact aligned target
    position (numbering shift, indel, or substitution), scan ±this many target
    sequence positions for a residue of the EXPECTED chemical class and map to
    the nearest match. Applies to every catalytic key (nucleophile, acid, base,
    clamp, cradle, stabiliser) so each is resolved the same dynamic way rather
    than relying on an exact-position hit.

    Expected classes are anchored to the canonical FAcD residues (§2.5, Chan
    et al. 2011, PDB 3R3U). Catalytic positions strictly conserved across the
    family — the Asp nucleophile (Asp110), Asp acid (Asp134), His base (His277)
    and Arg carboxylate clamp (Arg111/Arg114) — keep a narrow class (the
    canonical identity plus its protonation variants only). The halide-pocket
    positions tolerate the documented chemical alternatives: an H-bond-capable
    aromatic (Trp/Tyr/His) for the Trp156/Tyr217 cradle, and an H-bond donor /
    cationic group for the His155 fluoride stabiliser. PHE is excluded from both —
    its π-system offers no polar donor to stabilise the leaving F⁻.
    """
    RESIDUE_SEARCH_WINDOW: int = 5
    """
    Extra reach added to RESIDUE_SEARCH_WINDOW only when an ambiguous /
    non-standard residue (BXZJOU) sits inside the base window, so a genuine
    catalytic residue just beyond the ambiguous gap is still reachable.
    """
    RESIDUE_SEARCH_AMBIG_EXTENSION: int = 2
    ROLE_EXPECTED_RESIDUES: dict = field(default_factory=lambda: {
        "Nucleophile":       {"ASP", "ASH"},
        "Acid_Catalyst":     {"ASP", "ASH"},   # Asp134 — strictly conserved carboxylate dyad partner; no Glu drift (matches the Asp nucleophile)
        "Base_Catalyst":     {"HIS", "HID", "HIE", "HIP", "HSE", "HSD", "HSP"},
        "Carboxylate_Clamp": {"ARG"},   # Arg111/Arg114 — strictly conserved guanidinium clamp on the substrate carboxylate
        "Fluoride_Cradle":   {"TRP", "TYR", "TYM", "HIS", "HID", "HIE", "HIP", "HSE", "HSD", "HSP"},   # H-bond-capable aromatics only; PHE offers no polar donor to the leaving F⁻
        "Fluorine_Stabiliser": {"TRP", "TYR", "TYM", "HIS", "HID", "HIE", "HIP", "HSE", "HSD", "HSP", "ARG", "LYS", "LYN"},   # His155 stabilises leaving F⁻; PHE excluded — no polar H-bond donor
    })

    # -------------------------------------------------------------------------------
    # Step 2.8: Spatial validation cutoff for resolved residues
    # -------------------------------------------------------------------------------
    """
    After the class-aware ±window resolver maps a catalytic residue, the
    structure-level geometry check rejects any triad residue (Nuc/Base/Acid)
    whose nearest atom sits farther than this cutoff from the NEAREST LIGAND
    ATOM (the same metric the tier gate itself uses — ligand-size robust).
    A centroid-based test would falsely reject a genuine catalytic residue at
    the reactive head of a long-chain PFAS, whose centroid lies far down the
    chain. Prevents sequence-proximal but spatially remote residues (e.g. ASP
    on a distant loop) from inflating tier-gate inputs. Catalytic contact is
    typically ≤4–7 Å; 10 Å allows generous conformational flexibility while
    still excluding spatially absurd matches (e.g. ASP ~19 Å on another loop).
    """
    RESOLVER_SPATIAL_CUTOFF: float = 10.0   # Å (nearest residue-atom → nearest ligand-atom)

    # -------------------------------------------------------------------------------
    # Step 2.9: Protonation-variant → canonical residue name mapping
    # -------------------------------------------------------------------------------
    """
    Force-field engines (OPLS4, CHARMM, Amber) rename residues to encode
    protonation state. This map normalises them back to canonical 3-letter
    codes for class matching, structure filtering, and residue-class lookups.
    """
    PROTONATION_MAP: dict = field(default_factory=lambda: {
        "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
        "HSE": "HIS", "HSD": "HIS", "HSP": "HIS",
        "ASH": "ASP", "GLH": "GLU",
        "CYM": "CYS", "CYX": "CYS",
        "LYN": "LYS", "TYM": "TYR",
    })

    # -------------------------------------------------------------------------------
    # Step 2.10: Residue chemical-class sets (single source of truth)
    # -------------------------------------------------------------------------------
    """
    Amino-acid groupings by chemical property, consumed by Step 02 for residue-role
    classification and structure filtering. Defined here (not inline in the worker) so
    every script shares one definition. RESIDUE_CLASS_GROUPS sets include force-field
    protonation variants so HID/HIE/HIP etc. are recognised.
    """
    POSITIVE_RES: set = field(default_factory=lambda: {"ARG", "LYS", "HIS"})
    NEGATIVE_RES: set = field(default_factory=lambda: {"ASP", "GLU"})
    RES_PROPS: dict = field(default_factory=lambda: {
        "NonPolar": {"ALA", "VAL", "LEU", "ILE", "MET", "PRO", "GLY", "PHE", "TRP"},
        "Aromatic": {"PHE", "TYR", "TRP", "HIS"},
        "Polar":    {"SER", "THR", "ASN", "GLN", "CYS", "TYR", "HIS"},
        "Pos":      {"LYS", "ARG", "HIS"},
        "Neg":      {"ASP", "GLU"},
    })
    '''
    Single source of truth for π-stacking / π-cation ring-plane fitting. Covers ALL four
    aromatic residues so the ≥5-atom ring gate in Step 02 detects His (imidazole) and Trp
    (indole) — including the π interactions of the cradle residues Trp156/His155 — not only
    the Phe/Tyr six-membered carbons.
      PHE/TYR (6-ring): CG CD1 CD2 CE1 CE2 CZ
      HIS (imidazole):  CG ND1 CD2 CE1 NE2                 (5 atoms)
      TRP (indole):     CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 (9 atoms)
    '''
    AROMATIC_RING_ATOMS: list = field(default_factory=lambda:
        ["CG", "CD1", "CD2", "CE1", "CE2", "CZ",
         "ND1", "NE2", "NE1", "CE3", "CZ2", "CZ3", "CH2"])
    POLAR_SIDECHAIN_ATOMS: dict = field(default_factory=lambda: {
        "SER": {"OG"}, "THR": {"OG1"}, "ASN": {"OD1", "ND2"}, "GLN": {"OE1", "NE2"},
        "CYS": {"SG"}, "TYR": {"OH"}, "HIS": {"ND1", "NE2"},
        "TRP": {"NE1"}, "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"},
        "ARG": {"NH1", "NH2", "NE"}, "LYS": {"NZ"},
    })
    RESIDUE_CLASS_GROUPS: dict = field(default_factory=lambda: {
        "AROMATIC": {"TRP", "TYR", "TYM", "PHE", "HIS", "HID", "HIE", "HIP", "HSE", "HSD", "HSP"},
        "POSITIVE": {"ARG", "LYS", "LYN", "HIS", "HID", "HIE", "HIP", "HSE", "HSD", "HSP"},
        "NEGATIVE": {"ASP", "ASH", "GLU", "GLH"},
        "POLAR":    {"SER", "THR", "ASN", "GLN", "TYR", "TYM"},
    })

    # ===============================================================================
    # SECTION 3: INTERACTION GEOMETRY THRESHOLDS  (Step 02)
    # ===============================================================================
    """
    All distances are heavy-atom unless labelled (HA = H-to-acceptor).
    """

    # -------------------------------------------------------------------------------
    # Step 3.0: Spatial validation thresholds (Step 02/08)
    # -------------------------------------------------------------------------------
    CF_DIST_TOLERANCE: float       = 2.2    # Å  C–F bond-length ceiling for pairing a fluorine to its scissile carbon
    CLASH_DIST_TOLERANCE: float    = 2.5    # Å generic clash distance
    BACKBONE_CLASH_DIST: float     = 2.2    # Å ligand-tail atom vs protein backbone (N/CA/C/O) clash
    QC_MIN_DIST_5: float           = 5.0    # Å
    QC_MIN_DIST_4: float           = 4.0    # Å
    QC_RMSD_SQ_SUM: float          = 25.0   # RMSD squared
    TAIL_CLASH_RATIO: float        = 0.15   # size-fair clash-fraction VETO threshold (clashing tail atoms / ligand heavy atoms); paired with CLASH_MIN_FLOOR
    CLASH_MIN_FLOOR: int           = 3      # min absolute clash count required alongside the ratio veto — stops a 1–2 incidental clash vetoing a small ligand (size-fair both ways)
    MAINCHAIN_CLASH_PENALTY: float = 0.05   # competence_score penalty per clashing tail atom (graded, size-fair)
    CLASH_ABS_VETO: int            = 10     # reference clash-count ceiling reported as mainchain_clash_count; the pose veto uses the size-fair TAIL_CLASH_RATIO + CLASH_MIN_FLOOR (a flat count would bias against long PFAS)
    TRUST_SCORE_EXCELLENT: float   = 2.0    # Å Trust score limit 1
    TRUST_SCORE_ACCEPTABLE: float  = 3.0    # Å Trust score limit 2
    SOLVENT_SPHERE_RADIUS: float   = 20.0   # Å MD solvent sphere
    BOLTZ_THREADS: int             = 4      # Default thread limit
    PFAS_ORDER: list = field(default_factory=lambda: [
        "25_TFA", "26_Fluoroacetate", "27_Difluoroacetate", "7_PFBA", "8_PFPeA",
        "6_PFHxA", "13_PFHpA", "1_PFOA", "4_PFNA", "10_PFDA", "12_PFUnDA",
        "11_PFDoDA", "14_PFTrDA", "17_PFTeDA", "18_PFHxDA", "19_PFODA",
        "5_PFBS", "9_PFPeS", "3_PFHxS", "15_PFHpS", "2_PFOS", "16_PFDS",
        "20_GenX", "21_ADONA", "24_C6O4", "22_6-2-FTOH", "23_8-2-FTOH",
    ])

    # -------------------------------------------------------------------------------
    # Step 3.1: Hydrogen bonds (Schrödinger Maestro defaults)
    # -------------------------------------------------------------------------------
    """
    Maestro H-bond criteria: H···A ≤ 2.8 Å, donor (D–H···A) ≥ 120°, acceptor
    (H···A–X) ≥ 90°. Applied where explicit H are present (PrepWizard-prepared
    PDB). For H-FREE structures (Boltz-2 mmCIF) the heavy-atom donor–acceptor
    (D···A) proxy ≤ 3.5 Å is used instead (Jeffrey 1997; D···A ≈ H···A + ~1.0 Å).
    McDonald & Thornton (1994) for the H···A + angle geometry.
    """
    THRESHOLD_HB_DIST_HA: float    = 2.8    # Å  Maestro H···A maximum (explicit-H structures)
    THRESHOLD_HB_ANGLE_MIN: float  = 120.0  # °  Maestro donor minimum angle (D–H···A)
    THRESHOLD_HB_ANGLE_ACC: float  = 90.0   # °  Maestro acceptor minimum angle (H···A–X)
    THRESHOLD_HB_DIST_MAX: float   = 3.5    # Å  heavy-atom D···A proxy maximum (H-free Boltz-2 CIF)
    THRESHOLD_HB_DIST_MIN: float   = 2.4    # Å  donor–acceptor (D···A) minimum

    # -------------------------------------------------------------------------------
    # Step 3.2: Salt bridges (Schrödinger Maestro default)
    # -------------------------------------------------------------------------------
    """
    Maestro default: opposite-charge groups ≤ 5.0 Å. Underlying concept:
    Barlow & Thornton (1983) ≤ 4.0 Å; Kumar & Nussinov (2002) allow ~5–6 Å.
    """
    THRESHOLD_SALT_BRIDGE: float   = 5.0    # Å  Maestro salt-bridge maximum distance

    # -------------------------------------------------------------------------------
    # Step 3.3: Hydrophobic contacts
    # -------------------------------------------------------------------------------
    THRESHOLD_HYDROPHOBIC_MAX: float = 4.0  # Å  carbon–carbon centroid distance

    # -------------------------------------------------------------------------------
    # Step 3.4: π–π stacking (Schrödinger Maestro defaults)
    # -------------------------------------------------------------------------------
    """
    Maestro: face–face ≤ 4.4 Å, tilt ≤ 30°; edge–face ≤ 5.5 Å, tilt ≥ 60°.
    Coincides with McGaughey et al. (1998) aromatic-aromatic geometry.
    """
    THRESHOLD_PI_FACE: float       = 4.4    # Å  centroid–centroid, face–face (alias: PI_STACK_FACE_DIST_MAX)
    THRESHOLD_PI_EDGE: float       = 5.5    # Å  centroid–centroid, edge–face (alias: PI_STACK_EDGE_DIST_MAX)
    """
    Backward-compatibility aliases reference the THRESHOLD_PI_* source above
    (not independent literals) so the pair can never silently diverge.
    """
    PI_STACK_FACE_DIST_MAX: float  = THRESHOLD_PI_FACE  # Å  alias of THRESHOLD_PI_FACE
    PI_STACK_FACE_ANGLE_MAX: float = 30.0   # °  maximum tilt angle, face–face
    PI_STACK_EDGE_DIST_MAX: float  = THRESHOLD_PI_EDGE  # Å  alias of THRESHOLD_PI_EDGE
    PI_STACK_EDGE_ANGLE_MIN: float = 60.0   # °  minimum tilt angle, edge–face

    # -------------------------------------------------------------------------------
    # Step 3.5: π–cation (Schrödinger Maestro default)
    # -------------------------------------------------------------------------------
    """
    Maestro default: centroid-to-cation ≤ 6.6 Å; off-axis ≤ 30°.
    Underlying: Gallivan & Dougherty (1999).
    """
    THRESHOLD_PI_CATION_MAX: float = 6.6    # Å  Maestro ring-centroid to cation distance
    PI_CATION_ANGLE_MAX: float     = 30.0   # °  Maestro maximum off-axis angle

    # -------------------------------------------------------------------------------
    # Step 3.6: Halogen bonds (Schrödinger Maestro defaults)
    # -------------------------------------------------------------------------------
    """
    Maestro: X···A ≤ 3.5 Å. As donor: donor angle ≥ 140°, acceptor angle ≥ 90°.
    As acceptor: donor angle ≥ 120°, acceptor angle 90–170°. Auffinger (2004);
    Wilcken (2013). NB: aliphatic C–F (PFAS) is a weak σ-hole donor — F···N/O/S
    within range are recorded as fluorine contacts (§14.5 FP/FL/FF), NOT bona-
    fide halogen bonds; the angle gates apply to true C(sp²)–X donors only.
    """
    HALOGEN_BOND_DIST_MAX: float    = 3.5    # Å  halogen···acceptor (X···A) maximum
    HALOGEN_DON_ANGLE_MIN: float    = 140.0  # °  as-donor:    donor minimum angle (C–X···A)
    HALOGEN_DON_ACC_ANGLE_MIN: float = 90.0  # °  as-donor:    acceptor minimum angle
    HALOGEN_ACC_DON_ANGLE_MIN: float = 120.0 # °  as-acceptor: donor minimum angle
    HALOGEN_ACC_ANGLE_MIN: float    = 90.0   # °  as-acceptor: acceptor minimum angle
    HALOGEN_ACC_ANGLE_MAX: float    = 170.0  # °  as-acceptor: acceptor maximum angle

    # -------------------------------------------------------------------------------
    # Step 3.6b: Aromatic H-bonds (Schrödinger Maestro defaults)
    # -------------------------------------------------------------------------------
    """Weak H-bonds donated to aromatic π-acceptors. Levitt & Perutz (1988)."""
    AROM_HB_DIST_O_ACC: float       = 2.8    # Å  maximum distance, O acceptor
    AROM_HB_DIST_N_ACC: float       = 2.5    # Å  maximum distance, N= acceptor
    AROM_HB_DON_ANGLE_O: float      = 90.0   # °  donor minimum angle (O acceptor)
    AROM_HB_DON_ANGLE_N_MIN: float  = 108.0  # °  donor minimum angle (N= acceptor)
    AROM_HB_DON_ANGLE_N_MAX: float  = 130.0  # °  donor maximum angle (N= acceptor)
    AROM_HB_ACC_ANGLE_MIN: float    = 90.0   # °  acceptor minimum angle

    # -------------------------------------------------------------------------------
    # Step 3.6c: Water-mediated bridges
    # -------------------------------------------------------------------------------
    """
    Water-bridged H-bond: donor···Wat and Wat···acceptor both in range, with a
    permitted water-bridging angle (A···Wat···D). Maestro renders the same
    water-mediated H-bonds qualitatively, but the numeric ranges below are taken
    from Salentin et al. (2015) — Maestro does not expose explicit water-bridge
    cutoffs.
    """
    WATER_BRIDGE_DIST_MIN: float    = 2.5    # Å  min donor/acceptor···water distance
    WATER_BRIDGE_DIST_MAX: float    = 4.0    # Å  max donor/acceptor···water distance
    WATER_BRIDGE_OMEGA_MIN: float   = 75.0   # °  min water-bridging angle (A···Wat···D)
    WATER_BRIDGE_OMEGA_MAX: float   = 140.0  # °  max water-bridging angle
    WATER_BRIDGE_DON_ANGLE_MIN: float = 100.0  # °  donor minimum angle (D–H···Wat)

    # -------------------------------------------------------------------------------
    # Step 3.6d: Steric contacts (Schrödinger Maestro "Contacts")
    # -------------------------------------------------------------------------------
    """
    Contact quality = interatomic distance ÷ Σ van-der-Waals radii (Bondi 1964).
    ratio ≥ GOOD → acceptable; BAD/UGLY flag progressively severe clashes.
    H-bonds, salt bridges, and 1,4 interactions are excluded from the contact
    set (matching Maestro's "Exclude" checkboxes).
    """
    CONTACT_RATIO_GOOD: float       = 1.30   # ratio threshold — good contact
    CONTACT_RATIO_BAD: float        = 0.89   # ratio threshold — bad contact
    CONTACT_RATIO_UGLY: float       = 0.75   # ratio threshold — ugly (severe clash)
    CONTACT_EXCLUDE_HBOND: bool     = True    # exclude H-bonds from contacts
    CONTACT_EXCLUDE_SALT: bool      = True    # exclude salt bridges from contacts
    CONTACT_EXCLUDE_14: bool        = True    # exclude 1,4 interactions from contacts
    # van-der-Waals radii (Å) for the contact-ratio denominator. Bondi (1964).
    VDW_RADII: dict = field(default_factory=lambda: {
        "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
        "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98,
        "NA": 2.27, "K": 2.75, "MG": 1.73, "CA": 2.31, "ZN": 1.39,
    })
    VDW_RADIUS_DEFAULT: float       = 1.70    # Å  fallback vdW radius (carbon) for unlisted elements
    """
    Maximum donor–H bond length used to assign hydrogens to their heavy-atom
    donor when detecting true H···A geometry on protonated (PrepWizard) PDBs.
    """
    HB_DH_BOND_MAX: float           = 1.30    # Å  X–H covalent bond ceiling

    # -------------------------------------------------------------------------------
    # Step 3.7: Metal coordination
    # -------------------------------------------------------------------------------
    # Harding (2006): metal–ligand bond ≤ 2.8 Å.
    METAL_COORD_DIST_MAX: float    = 2.8    # Å  metal–ligand coordination bond
    METALS: set = field(default_factory=lambda: {
        "MG", "ZN", "MN", "CA", "FE", "CO", "NI", "CU", "NA", "K",
    })

    # -------------------------------------------------------------------------------
    # Step 3.8: General catalytic site
    # -------------------------------------------------------------------------------
    CATALYTIC_DIST_CUTOFF: float   = 6.0    # Å  residue included as "near active site"
    TAIL_MIN_BOND_DISTANCE: int    = 3      # topological bond distance (> this) defining ligand "tail" atoms for the mainchain-clash metric (Step 02)
    PLIP_CONTACT_FALLBACK_DIST: float = 5.0 # Å  fallback distance when a PLIP bs_residue lacks min_dist
    LIG_COVALENT_BOND_DIST: float  = 1.85   # Å  max inter-atom distance drawn as a covalent bond (2D interaction diagram)

    # -------------------------------------------------------------------------------
    # Step 3.9: Coordinate–structure match
    # -------------------------------------------------------------------------------
    COORD_MATCH_DIST_MAX: float    = 1.8    # Å  map residue to reference coordinate

    # ===============================================================================
    # SECTION 4: NAC (Near Attack Conformation) GEOMETRY  (Steps 02, 08)
    # ===============================================================================
    """
    Substrate: fluoroacetate; electrophile: C–F carbon; nucleophile: Asp O.
    Ideal SN2 backside attack: Nu–C–F collinear at 180°.
    References: Lightstone & Bruice (1996); Bruice (2002).
    """

    # -------------------------------------------------------------------------------
    # Step 4.1: Strict NAC — publication-grade catalytic viability
    # -------------------------------------------------------------------------------
    NAC_DIST_STRICT: float   = 3.2    # Å  nucleophile O to electrophilic C
    NAC_ANGLE_STRICT: float  = 155.0  # °  O–C–F attack angle at C (180° = ideal backside)

    # -------------------------------------------------------------------------------
    # Step 4.2: Relaxed NAC — pre-reactive / entropy-inclusive sampling
    # -------------------------------------------------------------------------------
    NAC_DIST_RELAXED: float  = 3.8    # Å  relaxed nucleophile–C distance
    NAC_ANGLE_RELAXED: float = 145.0  # °  relaxed attack angle

    # -------------------------------------------------------------------------------
    # Step 4.3: Pocket residency
    # -------------------------------------------------------------------------------
    POCKET_RESIDENCY_DIST: float = 8.0   # Å  ligand considered "pocket-bound" when nuc–C ≤ this

    # -------------------------------------------------------------------------------
    # Step 4.4: Fluoride cradle detection (TRP / TYR aromatic basket)
    # -------------------------------------------------------------------------------
    """
    The fluoride cradle is the TRP/TYR aromatic shell that stabilises the
    departing fluoride in FAcD-family dehalogenases.
    """
    F_CRADLE_RADIUS: float   = 12.0   # Å  search radius from nucleophile to TRP/TYR heavy atoms

    # ===============================================================================
    # SECTION 5: CATALYTIC TRIAD INTEGRITY  (Steps 02, 08)
    # ===============================================================================
    """
    Reference: FAcD crystal structure PDB 3R3U (Chan et al. 2011 — see header §3).
    NB = nucleophile O to catalytic base N; BA = catalytic base N to acid O.

    Two threshold sets:
      Static (Step 02) — calibrated on crystal structures (energy-minimised).
      MD     (Step 07) — +2.0 Å / +2.0 Å buffer for 300 K thermal fluctuations in
                         solution; justified by Asp–His distance variance in FAcD MD
                         trajectories (σ ≈ 1–2 Å at 300 K).
    """
    THRESHOLD_TRIAD_NB: float    = 4.5   # Å  crystal/static (Step 02)
    THRESHOLD_TRIAD_BA: float    = 7.0   # Å  crystal/static (Step 02)
    THRESHOLD_TRIAD_NB_MD: float = 6.5   # Å  MD-calibrated  (Step 07) = crystal + 2.0 Å
    THRESHOLD_TRIAD_BA_MD: float = 9.0   # Å  MD-calibrated  (Step 07) = crystal + 2.0 Å

    # -------------------------------------------------------------------------------
    # Step 5.1: Mechanistic-score & soft-score contact gates  (Step 02)
    # -------------------------------------------------------------------------------
    """
    Used by 02_Production analyse_candidate_structure() mech_score / soft_score.
    These are distinct from the §8 tier-cascade gates: they score the physical
    anchor set (halide stabiliser, carboxylate clamp, triad proximity) plus the
    graded SN2 attack angle, whereas §8 assigns the discrete tier. The Nuc–C gate
    reuses NAC_DIST_STRICT (§4.1) and the soft-score angle sigmoid reuses
    NAC_ANGLE_STRICT so the contact gates stay in lock-step. The mech_score angle
    term is graded separately via MECH_W_ANGLE (below).
    """
    MECH_STAB_RADIUS: float  = 5.5   # Å  TRP/TYR (or dynamic polar) → F⁻ halide-stabilisation contact
    MECH_CLAMP_RADIUS: float = 5.0   # Å  ARG carboxylate clamp → ligand contact
    """
    Mutation-tolerant fluoride-cradle detection. When the alignment-mapped
    canonical cradle residue (TRP/TYR/HIS) is absent (gap/substitution), an
    aromatic sidechain within this radius of the ligand's leaving halogen is
    accepted as a functional cradle — mirrors the dynamic nucleophile search.
    """
    MECH_CRADLE_RADIUS: float = 5.5  # Å  aromatic sidechain → leaving halogen (dynamic cradle)
    BOND_DIST_MAX: float     = 1.9   # Å  max heavy-atom separation treated as a covalent bond (structure-only geometry)
    MECH_NB_GATE: float      = 5.0   # Å  Nuc–Base distance gate (mech: MECH_W_NB)
    MECH_BA_GATE: float      = 5.5   # Å  Base–Acid distance gate (mech: MECH_W_BA)
    """
    Mechanistic-score component weights (sum = 1.00 → mech_score saturates at 1.0
    only for a complete anchor set AND an ideal 180° SN2 trajectory). The score is
    holistic: five binary anchor checks (weight 0.70 total) plus a graded SN2
    attack-angle term (weight 0.30), so a pose with intact catalytic machinery but a
    non-productive attack angle cannot read a perfect 1.00 — the angle is
    folded into mech, not scored separately. The angle term is Šidák multiplicity-
    corrected for the scissile C–F count: credit (1-p1)^n with p1=(1-cos δ)/2 and
    δ=180-angle, so a clean anti-attack on a mono-F carbon outscores the same angle
    reached as the best of two/three fluorines (CF2/CF3), which had more chances of
    one C–F landing near the 180° anti-axis. The discrete §8 tier angle-gates
    (174/165/155/145°) enforce SN2 linearity for tier entry. Each steric clash subtracts a small graded
    MECH_CLASH_PENALTY (capped via MECH_CLASH_PENALTY_MAX); the result is floored at 0.0.
    The penalty is deliberately small: a "teflon clash" is a ligand fluorine near the Asp
    oxygen, which is intrinsic to a fluorinated substrate sitting in the active site, so it
    nudges rather than nullifies the score.
    """
    MECH_W_NUC: float        = 0.15   # nucleophile NAC reach   (d_nuc ≤ NAC_DIST_STRICT)
    MECH_W_NB: float         = 0.10   # Nuc–Base relay          (≤ MECH_NB_GATE)
    MECH_W_BA: float         = 0.10   # Base–Acid relay         (≤ MECH_BA_GATE)
    MECH_W_CLAMP: float      = 0.10   # carboxylate clamp present
    MECH_W_STAB: float       = 0.25   # halide (F⁻) stabilisation present
    MECH_W_ANGLE: float      = 0.30   # graded SN2 attack angle, Šidák multiplicity-corrected (1-p1)^n for scissile C–F count
    MECH_CLASH_PENALTY: float = 0.03  # mech points subtracted per steric clash (small, graded)
    MECH_CLASH_PENALTY_MAX: float = 0.15  # cap on total clash penalty so clashes never dominate
    # Backside steric-occlusion (geometric SN2 determinant, Bento & Bickelhaupt 2008: the
    # nucleophile needs an open backside anti to the leaving F) is a single-source feasibility
    # penalty applied ONCE, in the graded angle-faded chemistry term (§5.2b, CHEM_PEN_W_OCCL),
    # which feeds mechanistic_score_effective. It is deliberately NOT also subtracted inside the
    # raw mechanistic_score below — that would double-penalise an occluded trajectory.

    # -------------------------------------------------------------------------------
    # Step 5.2: SN2 dead-end consensus check (scissile C–F energy + backside sterics)
    # -------------------------------------------------------------------------------
    """
    A pose is flagged a non-productive SN2 dead-end only when BOTH indicators agree
    (consensus → robust against single-signal false positives), and even then the final
    chemical verdict is deferred to Step-08 QM/MM + MD/WaterMap.

    A — Scissile C–F bond-dissociation energy (kcal/mol), keyed by the number of
        fluorines on the mapped attack carbon. α-Fluorination strengthens the C–F bond,
        so a CF3 carbon's C–F is too strong to cleave by SN2.
        Values: CH3F 109.9 · CH2F2 119.5 · CHF3 127.5 (fluoromethane series;
        O'Hagan 2008 — see header Scientific References §9).

    B — Backside steric occlusion: summed Bondi van-der-Waals radii of the halogen
        substituents on the attack carbon that crowd the nucleophile's backside approach
        (anti to the leaving F), measured from the docked pose. Above the cutoff the
        pentacoordinate SN2 transition state is sterically blocked
        (Bento & Bickelhaupt 2008; Bondi 1964 vdW radii — see header §9).
    """
    SCISSILE_CF_BDE: dict = field(default_factory=lambda: {1: 109.9, 2: 119.5, 3: 127.5})
    SCISSILE_CF_BDE_MAX: float   = 123.0   # kcal/mol; above → C–F too strong to cleave (3F=127.5 fails, 2F=119.5 passes)
    SN2_BACKSIDE_OCCL_MAX: float = 2.0     # Å (Σ vdW); above → backside SN2 approach sterically blocked (2×F=2.94 fails, 1×F=1.47 passes)

    # Graded chemical-feasibility factor (feasibility_factor) — continuous, never a veto.
    FEASIBILITY_FLOOR: float     = 0.50    # lowest the feasibility multiplier can reach (a recalcitrant but real substrate such as TFA is down-ranked, not buried — QM/MM is the true arbiter)
    FEAS_BDE_LO: float           = 120.0   # kcal/mol; scissile C–F BDE ≤ this → no BDE penalty (FA 109.9, DFA 119.5 pass)
    FEAS_BDE_HI: float           = 132.0   # kcal/mol; BDE ≥ this → full BDE penalty. TFA (127.5) lands graded (~0.38), not the floor, so it still ranks as a lead — but the Tier_1A bond-strength ceiling (TIER_ELITE_BDE_MAX, §8.5) caps it at Tier_1B regardless of its pose
    FEAS_BETA_PER_F: float       = 0.35    # per-β-fluorine penalty: f_beta = 1/(1 + this·β_F) (FA/DFA β=0; PFAS β≥2)

    # -------------------------------------------------------------------------------
    # Step 5.2b: Graded feasibility + pocket-fit penalties folded into the tier-gate mech
    # -------------------------------------------------------------------------------
    """
    Continuous penalties subtracted from the geometric mechanistic_score to form the
    feasibility-weighted score that gates the degrader tier (competence/ranking keep the
    raw geometry). Both engage only past a chemistry/steric threshold, so genuine
    substrates are unpenalised and the demotion is graded, never a hard class veto.

      Chemistry  : penalty = CHEM_PEN_W_BDE·max(0, BDE − SCISSILE_CF_BDE_MAX)          [no fade]
                           + CHEM_PEN_W_OCCL·max(0, occlusion − SN2_BACKSIDE_OCCL_MAX) [angle-faded]
                   The two terms answer different questions and therefore behave differently.
                   Backside occlusion is a STERIC obstruction of the attack trajectory: an
                   enzyme that organises the substrate into a near-linear near-attack geometry
                   has, by construction, cleared that trajectory, so the occlusion term fades to
                   zero as the angle approaches 180° (CHEM_PEN_ANGLE_* below).
                   The C–F bond dissociation energy is an INTRINSIC property of the bond being
                   broken. No approach geometry lowers it: a perfect 180° trajectory onto a
                   127 kcal/mol C–F still has to break a 127 kcal/mol C–F. The BDE term is
                   therefore flat in angle, and a substrate whose scissile C–F exceeds
                   TIER_ELITE_BDE_MAX cannot hold the elite tier at any angle (§8.5) — it may
                   still surface as a lower-tier discovery lead, where Step-07 QM/MM is the
                   arbiter of whether the barrier is in fact surmountable.
      Containment: penalty = CONTAIN_PEN_W·max(0, CONTAIN_PEN_TARGET − pocket_containment_cavity)
                   FAcD is a small-substrate (haloacetate) hydrolase (Wackett 2022; Chan 2011),
                   but the penalty must express the POCKET, not the ligand: the cavity metric is
                   measured against the protein the ligand is actually docked into (§5.2c), so a
                   genuine wide-pocket homolog that really does enclose a longer chain is spared
                   and the same ligand can score differently in different enzymes. Long-PFAS
                   hydrolytic-SN2 hits remain EXPLORATORY, not degraders.
    """
    """
    CHEM_PEN_W_BDE is the flat cost per kcal/mol of scissile C–F strength above SCISSILE_CF_BDE_MAX.
    It answers one question — how much harder is this bond to break — and nothing else.

    It is deliberately modest (an α-CF3 pays ~0.07 of mech score) because the discrimination that
    keeps a poly-fluorinated substrate honest is now made where it belongs: in the multiplicity-
    corrected attack angle (§5.2d), which strips the best-of-N geometric advantage a CF3 carbon
    enjoys before the tier ladder ever sees it. A heavier BDE penalty would double-charge the same
    substrate — once for its bond strength and again, implicitly, for a pose advantage that has
    already been removed — and would bar an α-CF3 from the elite tier by arithmetic rather than by
    evidence. The bond-strength ceiling (TIER_ELITE_BDE_MAX, §8.5) remains the hard limit, and
    Step-07 QM/MM remains the arbiter of whether the barrier is actually surmountable.
    """
    CHEM_PEN_W_BDE: float    = 0.015   # penalty per kcal/mol of scissile C–F BDE above SCISSILE_CF_BDE_MAX
    CHEM_PEN_W_OCCL: float   = 0.13    # penalty per Å of backside occlusion above SN2_BACKSIDE_OCCL_MAX
    CHEM_PEN_W_BETA: float   = 0.08    # penalty per β-fluorine on the attack-carbon chain: β-fluorination inductively withdraws electron density from the α-C–F, raising its cleavage barrier beyond the raw α-F-count BDE. Continuous and pose-independent; a substrate with no β-fluorine (β_F=0) receives no β term and is governed by the α-BDE/occlusion penalty above
    # The angle fade applies to the backside-occlusion term ONLY: full at/below
    # CHEM_PEN_ANGLE_FULL, zero at/above CHEM_PEN_ANGLE_NONE, linear between. The BDE,
    # β-fluorination and containment terms do not fade — none of them is a trajectory
    # obstruction that a good angle can relieve. Step-07 QM/MM remains the final arbiter.
    CHEM_PEN_ANGLE_FULL: float = 175.0   # ° SN2 angle at/below which the occlusion penalty applies in full
    CHEM_PEN_ANGLE_NONE: float = 180.0   # ° SN2 angle at/above which the occlusion penalty is fully waived

    # -------------------------------------------------------------------------------
    # Step 5.2c: Pocket containment — TWO protein-aware measurements
    # -------------------------------------------------------------------------------
    """
    Containment asks whether the ENZYME holds the ligand, so both measurements are made
    against protein coordinates. They answer two different questions and are reported
    side by side for every pose.

    (A) pocket_containment_cavity — does the protein cavity enclose the ligand?
        Each ligand heavy atom casts BURIAL_RAYS rays over the sphere; a ray is blocked when
        a protein heavy atom obstructs it within BURIAL_PROBE_A. The atom's buriedness is the
        blocked fraction, and the atom counts as contained at/above BURIAL_MIN. The metric is
        the mean over ligand heavy atoms. This is the term the tier gate penalises: it is the
        physical 'can this enzyme hold this molecule' question, and it varies with the protein,
        so a wide-pocket homolog is not punished for the ligand's intrinsic length.

    (B) pocket_containment_site8 — how much of the ligand sits inside the catalytic
        constellation? Fraction of ligand heavy atoms within SITE8_SHELL_A of any heavy atom
        of the eight mapped active-site residues (§2.5). This measures catalytic ENGAGEMENT,
        not cavity fit: a long tail leaving the shell says the tail is outside the reactive
        machinery, which is a mechanistic statement rather than a size penalty. Reported and
        plotted; it does not gate.
    """
    BURIAL_RAYS: int          = 42     # rays per ligand heavy atom (icosphere-like Fibonacci sphere)
    BURIAL_PROBE_A: float     = 8.0    # Å; ray length searched for a blocking protein atom
    BURIAL_RAY_CLEARANCE: float = 1.8  # Å; a protein heavy atom within this of the ray axis blocks it
    """
    BURIAL_MIN is calibrated on the corpus, not assumed. Measured mean contained-atom fraction over
    5 poses per ligand at four candidate thresholds:

        buriedness ≥      0.50    0.70    0.80    0.90
        fluoroacetate     1.00    0.92    0.76    0.64
        difluoroacetate   1.00    1.00    0.97    0.67
        TFA               1.00    1.00    1.00    0.71
        PFBA  (C4)        1.00    0.99    0.63    0.26
        PFHxA (C6)        1.00    0.99    0.85    0.21
        PFOA  (C8)        1.00    0.91    0.66    0.08
        PFTeDA(C14)       0.76    0.36    0.15    0.02
        PFODA (C18)       0.62    0.35    0.12    0.03

    At 0.50 the metric saturates — every pose reads 1.00 and it discriminates nothing. At 0.80 and
    above it goes non-monotonic in size (fluoroacetate 0.76 falls BELOW PFHxA 0.85): those thresholds
    are measuring how deep a small ligand sits, not whether the pocket holds it. 0.70 keeps the
    native substrates high, holds the mid-chain PFCAs near 0.9, and separates the genuinely oversized
    C14+ chains that spill out of any pocket.
    """
    BURIAL_MIN: float         = 0.70   # buriedness at/above which a ligand heavy atom counts as cavity-contained
    SITE8_SHELL_A: float      = 5.0    # Å; ligand heavy atom within this of an active-site residue atom is engaged
    CONTAIN_PEN_TARGET: float = 0.85   # cavity containment at/above this → no penalty
    CONTAIN_PEN_W: float      = 1.00   # penalty per unit of cavity-containment shortfall below the target

    # -------------------------------------------------------------------------------
    # Step 5.2d: Angle multiplicity — the Šidák exponent and the effective attack angle
    # -------------------------------------------------------------------------------
    """
    A poly-fluorinated attack carbon gets more than one chance at a near-linear backside angle, and
    that advantage is geometric, not catalytic. An α-CF3 has three equivalent C–F bonds arranged
    about the Cα–COO⁻ axis: rotate the head group and SOME fluorine always lands roughly opposite the
    nucleophile. Fluoroacetate has a single C–F and must be oriented exactly.

    Measured on the DeHa4 control across all five Boltz diffusion samples, trifluoroacetate
    out-angles the native substrate in EVERY sample (151–159° vs 95–145°). Left uncorrected, that
    best-of-3 buys TFA a higher tier than fluoroacetate on the very enzyme that is known not to turn
    TFA over (Wackett 2022) — the control inverts, and the ladder is measuring fluorine count rather
    than catalytic competence.

    sn2_effective_angle() removes the inflation: it converts the observed angle into the angle a
    SINGLE-C–F substrate would have to show to be equally improbable. p1 = (1-cos δ)/2 is the
    single-bond chance of landing within δ of linear; the pose had n such chances, so its Šidák
    survival is (1-p1)^n, and the effective angle is the one whose single-bond probability equals it.
    A mono-fluoro substrate is unchanged by construction; a CF3 pose is deflated in proportion to how
    mediocre it is, and a genuinely near-ideal CF3 pose (≈178°) barely moves — which is the intent.
    The tier ladder gates on this effective angle; the raw angle stays reported.
    """
    def sn2_effective_angle(self, angle: float, scissile_f_count: int = 1) -> float:
        import math
        _n = max(1, int(scissile_f_count))
        if _n == 1:
            return float(angle)
        _delta = max(0.0, 180.0 - float(angle))
        _p1 = (1.0 - math.cos(math.radians(_delta))) / 2.0
        _q = max(0.0, min(1.0, (1.0 - _p1) ** _n))          # Šidák survival across the n bonds
        _p_eff = 1.0 - _q                                    # equivalent single-bond probability
        _cos = max(-1.0, min(1.0, 1.0 - 2.0 * _p_eff))
        return float(180.0 - math.degrees(math.acos(_cos)))

    # -------------------------------------------------------------------------------
    # Step 5.3: Reactive-centre gating — α-carbon attack + bidentate carboxylate clamp
    # -------------------------------------------------------------------------------
    """
    Reactive-centre gating ties every ligand-side gate to the same mechanistically
    correct atom — the α-carbon (the carbon bonded to the substrate carboxylate, the
    position FAcD defluorinates) — and requires the substrate carboxylate to be held by
    BOTH clamp arginines. A large PFAS therefore cannot satisfy the nucleophile-reach or
    clamp gates with stray fluorines, nor score an SN2 angle on an incidental C–F.

      • SCISSILE_REQUIRE_ALPHA — the SN2 attack carbon used for the distance/angle gates
        must be the α-carbon adjacent to the ligand carboxylate; a mid-chain CF2 near the
        nucleophile does not qualify. FAcD attacks Cα of a 2-haloalkanoate (Chan et al.
        2011; Kurihara & Esaki 2008 — see header §3/§6).
      • CLAMP_REQUIRE_BIDENTATE — the elite (Tier_1A) gate requires the ligand carboxylate
        oxygens to salt-bridge BOTH distinct clamp arginines (Arg111/Arg114 equivalents).
        Non-carboxylate heads (sulfonate –SO3⁻, ether) cannot satisfy this and so cannot
        reach the elite tier via the FAcD carboxylate-anchoring mechanism.
      • CLAMP_SALT_BRIDGE_DIST — max carboxylate-O ↔ Arg-guanidinium-N separation for a
        clamp arm to count as engaged (Maestro salt-bridge geometry; Donald et al. 2011).
    """
    SCISSILE_REQUIRE_ALPHA: bool  = True
    CLAMP_REQUIRE_BIDENTATE: bool = True
    CLAMP_SALT_BRIDGE_DIST: float = 4.0   # Å  carboxylate O ↔ Arg guanidinium N (salt-bridge contact)

    # -------------------------------------------------------------------------------
    # Step 5.4: Nucleophile alignment-rescue audit (±window resolver QC)
    # -------------------------------------------------------------------------------
    """
    The ±RESIDUE_SEARCH_WINDOW resolver (§2.7) maps the canonical Asp110 column through the
    sequence alignment, then — if that column is a gap/substitution — scans nearby target
    positions for an Asp. Because Asp is common, a wide rescue can latch onto a
    NON-catalytic Asp. The resolution method and offset are recorded per residue
    (master-CSV column 'nuc_resolution') for QC, and the elite tier requires the
    nucleophile to be a direct alignment hit or a TIGHT rescue (offset ≤
    NUC_RESCUE_MAX_OFFSET_ELITE) so a far-fetched rescued Asp cannot seed a Tier_1A call.
    The intra-protein Nuc–Base distance gate rejects a spatially remote rescue; this is an
    explicit, auditable second guard.
    """
    NUC_RESCUE_MAX_OFFSET_ELITE: int = 2   # residues; max |rescue − aligned column| for elite eligibility

    # -------------------------------------------------------------------------------
    # Step 5.5: Model consensus (robust pose across Boltz diffusion samples)
    # -------------------------------------------------------------------------------
    """
    model_degrader_consensus is the fraction of a candidate's Boltz diffusion samples that
    independently reach a degrader tier. The representative pose is the best-tier model
    (see select_best_degrader_model). The consensus fraction is REPORTED per candidate
    (master-CSV column 'model_degrader_consensus' and the Ranking_Score_Calc string) and is
    used as the FINAL Scientific_Rank tiebreaker — applied only after tier, competence and
    active-site conservation, so it breaks ties between geometrically equal poses (a
    reproducible pose above a single-frame fluke) without ever crossing a tier or competence
    boundary. It is NEITHER a hard degrader gate NOR a filter: a conformationally flexible
    true substrate keeps its tier and is never dropped for its sampling spread.
    """

    def mechanistic_score(self, d_nuc: float, dist_nuc_base: float,
                          dist_base_acid: float, clamp_ok: bool, stabilised: bool,
                          angle: float, steric_clashes: int = 0,
                          scissile_cf_bde: float = 0.0,
                          backside_occlusion: float = 0.0,
                          beta_f_count: int = 0,
                          angle_multiplicity: int = 1) -> float:
        """
        Single source of truth for the holistic mechanistic score (0–1) — the raw geometry
        term behind the tier-gate key (TIER_MECH_MIN, via mechanistic_score_effective). Five
        binary anchor checks (weights total 0.70) plus a graded SN2 attack-angle term
        (MECH_W_ANGLE), less MECH_CLASH_PENALTY per steric clash. GEOMETRY + active-site-
        machinery only: NO chemical-feasibility (C–F BDE, backside occlusion, β-fluorination)
        — those are the graded penalties in mechanistic_score_effective (§5.2b), applied once
        each. Tiering stays geometric here, no ligand excluded a priori; chemistry rides in
        the effective score, competence/diagnostics and Step-08 QMMM.

        The angle term is Šidák multiplicity-corrected by angle_multiplicity — the number of
        independent chances the pose had at a near-linear angle (§5.2d), which the caller computes:
        the two aspartate oxygens always, times the fluorines on the scissile carbon only when the
        cradle did NOT fix the leaving F. A CF2/CF3 attack carbon presents more equivalent C–F bonds, so
        more chances of one landing near the 180° anti-axis. The credit (1-p1)^n, with
        p1=(1-cos δ)/2 and δ=180-angle, removes that best-of-N inflation from tiering — the
        same statistic used in competence_score, applied here so the tier gate is not gamed
        by fluorine multiplicity. Backside occlusion is NOT applied here — it is a single-
        source feasibility penalty in mechanistic_score_effective (§5.2b). Floored at 0,
        rounded to 2 dp.
        """
        import math
        s = 0.0
        if d_nuc <= self.NAC_DIST_STRICT:        s += self.MECH_W_NUC
        if dist_nuc_base <= self.MECH_NB_GATE:   s += self.MECH_W_NB
        if dist_base_acid <= self.MECH_BA_GATE:  s += self.MECH_W_BA
        if clamp_ok:                             s += self.MECH_W_CLAMP
        if stabilised:                           s += self.MECH_W_STAB
        """
        Šidák multiplicity-corrected attack-angle credit. The exponent is the number of equivalent
        C–F bonds on the scissile carbon (§5.2d) — the number of chances the POSE had at presenting
        some fluorine anti-periplanar to the nucleophile.

        The multiplicity is a property of the substrate's geometry, not of how the leaving fluorine
        is later identified. An α-CF3 carbon has three-fold symmetry about the Cα–COO⁻ axis: rotate
        the head group and SOME fluorine always ends up reasonably backside. Fluoroacetate has one
        C–F and two hydrogens, so its single fluorine must be oriented exactly. Resolving the leaving
        F by the fluoride cradle tells us WHICH bond breaks; it does not undo the fact that the CF3
        had three ways to look good. Measured on the DeHa4 control across all five diffusion samples,
        TFA out-angles fluoroacetate in every one (151–159° vs 95–145°) — a best-of-3 advantage that
        is not a catalytic one, and precisely the inflation this correction exists to remove.
        """
        _delta = max(0.0, 180.0 - float(angle))
        _p1    = (1.0 - math.cos(math.radians(_delta))) / 2.0
        _n     = max(1, int(angle_multiplicity))
        _q     = max(0.0, min(1.0, (1.0 - _p1) ** _n))
        s += self.MECH_W_ANGLE * _q
        if steric_clashes > 0:
            s -= min(self.MECH_CLASH_PENALTY * steric_clashes, self.MECH_CLASH_PENALTY_MAX)
        # Backside occlusion is deliberately NOT subtracted here: it is applied once, as the
        # angle-faded chemistry penalty in mechanistic_score_effective (§5.2b). The unused
        # backside_occlusion parameter is retained for call-site signature stability.
        return round(max(0.0, s), 2)

    def feasibility_factor(self, scissile_cf_bde: float = 0.0, beta_f_count: int = 0) -> float:
        """
        Graded chemical-feasibility multiplier ∈ [FEASIBILITY_FLOOR, 1.0] for enzymatic
        defluorination. Two continuous, physically-grounded penalties:
          • α-fluorine C–F bond strength — ramps 1.0 (≤ FEAS_BDE_LO) down to the floor
            (≥ FEAS_BDE_HI). O'Hagan 2008 (C–F is the strongest single bond; geminal F
            raises it): FA 109.9 → 1.0, DFA 119.5 → ~1.0, TFA 127.5 → floor.
          • β-fluorination — chain/ether perfluorination withdrawing density from the
            reactive centre and marking a recalcitrant perfluoroalkyl substrate. FA/DFA
            (β = 0) → 1.0; PFAS/GenX (β ≥ 2) → penalised by FEAS_BETA_PER_F per β-F.
        Product of the two, floored. Substrate-class agnostic otherwise.
        """
        _bde_rng = self.FEAS_BDE_HI - self.FEAS_BDE_LO
        f_bde = 1.0 if _bde_rng <= 0 else (self.FEAS_BDE_HI - scissile_cf_bde) / _bde_rng
        f_bde = max(0.0, min(1.0, f_bde))
        f_beta = 1.0 / (1.0 + self.FEAS_BETA_PER_F * max(0, int(beta_f_count)))
        return max(self.FEASIBILITY_FLOOR, f_bde * f_beta)

    # -------------------------------------------------------------------------------
    # Step 5.6: Gated continuous competence score (Step 02 — Scientific-ranking key)
    # -------------------------------------------------------------------------------
    """
    The Scientific-ranking key. A single continuous 0–1 score in which every independent
    catalytic axis enters exactly ONCE, under hard chemistry gates. Unlike mechanistic_score
    (binary anchors → saturates near 1.0 within a tier) and soft_catalytic_score (geometry
    only → high even for dead-ends), competence_score discriminates across the whole range
    and collapses to 0 for a non-productive pose:

      hard gate (→ 0.0): the SN2 attack carbon must be the α-carbon adjacent to a ligand
        carboxylate (a non-carboxylate head or mid-chain attack scores 0). The A+B SN2
        dead-end is NOT gated here — it rides in the continuous feasibility_factor scaling
        below (which discounts a strong-C–F / β-fluorinated centre in the rank) and is decided
        finally by Step-08 QM/MM, matching the diagnostic-only dead-end treatment in
        02_Production §7.2.3.
      graded terms (each 0–1, weights sum to 1.0):
        • angle  — Šidák multiplicity-corrected (1-p1)^n for scissile C–F count, the Walden backside trajectory
        • dist   — Asp-Oδ → α-carbon, closer within [COMP_DIST_MIN, COMP_DIST_MAX] = higher
        • clamp  — carboxylate_clamp_integrity (0 / 0.5 / 1.0)
        • traj   — SN2 trajectory deviation, smaller = higher
        • triad  — Nuc–Base and Base–Acid relay closeness (mean)
        • halide — fluoride-cradle stabilisation present
    """
    COMP_DIST_MIN: float   = 2.4    # Å  ideal Asp-Oδ → α-carbon (full distance credit at/below)
    COMP_DIST_MAX: float   = 3.8    # Å  relaxed reach (zero distance credit at/above)
    COMP_TRAJ_MAX: float   = 2.5    # Å  trajectory deviation giving zero traj credit
    COMP_RELAY_MIN: float  = 2.5    # Å  relay distance giving full triad credit at/below
    COMP_W_ANGLE: float    = 0.30
    COMP_W_DIST: float     = 0.25
    COMP_W_CLAMP: float    = 0.20
    COMP_W_TRAJ: float     = 0.10
    COMP_W_TRIAD: float    = 0.10
    COMP_W_HALIDE: float   = 0.05

    def competence_score(self, scissile_is_alpha: bool, head_is_carboxylate: bool,
                         angle: float, d_nuc: float,
                         clamp_integrity: float, sn2_trajectory_dev: float,
                         stabilised: bool, dist_nuc_base: float,
                         dist_base_acid: float,
                         scissile_cf_bde: float = 0.0, beta_f_count: int = 0,
                         scissile_f_count: int = 1) -> float:
        """
        Single source of truth for the gated continuous competence score (0–1), the
        within-tier ranking key. Returns 0.0 for any non-productive pose (not α-carbon
        attack or non-carboxylate head). Six graded catalytic-geometry axes (angle,
        nucleophile distance, clamp, trajectory, triad relay, halide) are summed.

        Two physically-grounded layers, applied to the WITHIN-TIER ranking only (tiers
        stay pure geometry — feasibility never gates tier entry):
          1. Geometry: six graded catalytic axes, with the attack-angle term Šidák
             multiplicity-corrected for the scissile C–F count (scissile_f_count). A carbon
             bearing more fluorines (CF2/CF3) has more independent chances of one C–F landing
             anti-periplanar by chance; credit (1-p1)^n, p1=(1-cos δ)/2, rewards an alignment
             that is unlikely given n bonds. FA (1 F) uncorrected; DFA/PFCA-α (2 F), TFA (3 F)
             corrected on geometric multiplicity.
          2. Reactivity: the geometric score is scaled by feasibility_factor (scissile C–F
             BDE + β-fluorination). This is the measured physical reactivity barrier — the
             reason FAcD turns over FA/DFA but not TFA/long-PFCAs (O'Hagan 2008; β-F electron
             withdrawal). At the α-CF2 reactive centre a long PFCA is geometrically identical
             to DFA, so geometry ALONE ranks them equal; β-fluorination is the only
             discriminator, and it belongs in the rank, not hidden. NOT a substrate-class
             label — a continuous, physically-derived factor; QM/MM (Step 07) is the final
             arbiter. Applied to competence (rank), NOT to mech_score (tier gate), so it
             orders FA/DFA above recalcitrant poses WITHIN a tier without excluding any ligand
             from a tier a priori. Any script calls CFG.competence_score(...).
        """
        if (not scissile_is_alpha) or (not head_is_carboxylate):
            return 0.0
        import math
        def _closer(x, lo, hi):
            if hi <= lo: return 0.0
            return max(0.0, min(1.0, (hi - x) / (hi - lo)))
        '''
        Attack-angle credit, Šidák-corrected for scissile C–F multiplicity.
        δ  = deviation from the ideal 180° anti-periplanar (backside) attack
        p1 = fraction of random directions within δ of the anti-axis = (1-cos δ)/2
             (the chance a single C–F lands this close to ideal by chance)
        p_any = 1-(1-p1)^n = chance ANY of the n equivalent C–F bonds aligns this well
        credit = 1 - p_any = (1-p1)^n → high only when the alignment is UNLIKELY by chance
             given n bonds: a clean anti-attack on a 1-F carbon outscores the same angle
             reached as the best of 2/3 fluorines (CF2/CF3), which had more chances.
        '''
        _delta = max(0.0, 180.0 - float(angle))
        _p1    = (1.0 - math.cos(math.radians(_delta))) / 2.0
        _n     = max(1, int(scissile_f_count))
        angle_term  = max(0.0, min(1.0, (1.0 - _p1) ** _n))
        dist_term   = _closer(d_nuc, self.COMP_DIST_MIN, self.COMP_DIST_MAX)
        clamp_term  = max(0.0, min(1.0, clamp_integrity))
        traj_term   = max(0.0, min(1.0, (self.COMP_TRAJ_MAX - sn2_trajectory_dev) / self.COMP_TRAJ_MAX)) if self.COMP_TRAJ_MAX > 0 else 0.0
        triad_term  = 0.5 * (_closer(dist_nuc_base, self.COMP_RELAY_MIN, self.MECH_NB_GATE)
                             + _closer(dist_base_acid, self.COMP_RELAY_MIN, self.MECH_BA_GATE))
        halide_term = 1.0 if stabilised else 0.0
        s = (self.COMP_W_ANGLE * angle_term + self.COMP_W_DIST * dist_term
             + self.COMP_W_CLAMP * clamp_term + self.COMP_W_TRAJ * traj_term
             + self.COMP_W_TRIAD * triad_term + self.COMP_W_HALIDE * halide_term)
        # Scale geometric competence by the physical reactivity barrier (within-tier rank
        # only; never a tier gate). Separates a long PFCA from DFA — geometrically identical
        # at the α-CF2 centre, distinguished solely by β-fluorination / C–F BDE.
        s *= self.feasibility_factor(scissile_cf_bde, beta_f_count)
        return round(max(0.0, min(1.0, s)), 3)
    """
    Soft-score (sigmoid) triad midpoints; nucleophile/angle sigmoids reuse the
    strict NAC cutoffs directly (NAC_DIST_STRICT / NAC_ANGLE_STRICT).
    """
    SOFT_NB_MIDPOINT: float  = 4.5   # Å  soft s_int Nuc–Base sigmoid midpoint (= THRESHOLD_TRIAD_NB)
    SOFT_BA_MIDPOINT: float  = 5.0   # Å  soft s_int Base–Acid sigmoid midpoint
    # soft_catalytic_score component weights (must sum to 1.0): nucleophile reach, SN2 angle, triad integrity.
    SOFT_W_NUC: float        = 0.4   # weight on s_nuc (nucleophile NAC reach)
    SOFT_W_ANG: float        = 0.3   # weight on s_ang (SN2 attack angle)
    SOFT_W_INT: float        = 0.3   # weight on s_int (triad-integrity product)
    # Sigmoid steepness (k) for the soft-score terms; sign sets direction (negative = higher score below x0).
    SOFT_K_NUC: float        = -4.0  # s_nuc nucleophile-distance sigmoid steepness
    SOFT_K_ANG: float        = 0.15  # s_ang SN2-angle sigmoid steepness
    SOFT_K_TRIAD: float      = -2.0  # s_int triad-relay sigmoid steepness (Nuc–Base and Base–Acid)
    """
    Intentionally < THRESHOLD_TRIAD_BA (7.0 Å): sigmoid midpoint sets
    the steepest scoring gradient in the 4–6 Å pre-reactive range;
    the hard 7.0 Å cutoff is the gate, not the sigmoid centre.
    """

    # ===============================================================================
    # SECTION 6: SN2 / WALDEN INVERSION GEOMETRY  (Step 07)
    # ===============================================================================
    """
    SN2 backside attack at sp³ C: ideal Walden-inversion trajectory = 180°.
    Note: the Bürgi–Dunitz angle (107°) applies to nucleophilic addition at
    sp² carbonyl carbons; it must NOT be conflated with the linear SN2 angle
    here (see calculate_burgi_dunitz() in 00_03 — auxiliary metric only).
    """
    WALDEN_IMPROPER_MAX: float = 15.0  # °  |improper dihedral| < this → TS-like (planar) geometry
    WALDEN_TS_FRAME_BONUS: float = 1.1  # MD frame-score multiplier for a TS-flat (Walden) frame (Step 07 frame selection)
    SN2_ANGLE_MARGINAL_MIN: float = 120.0  # °  lower bound of the marginal SN2-angle band for figure colour-coding (Step 07); ≥ NAC_ANGLE_RELAXED is favourable
    MD_EAF_SMOOTH_WINDOW: int = 50      # frames — rolling-average window for EAF-MSA trajectory smoothing (Step 07; frame count, independent of stride)

    # ===============================================================================
    # SECTION 7: SMART-LOCK RESIDUE DETECTION  (Step 07)
    # ===============================================================================
    """
    Geometry-biased scoring that steers 3D triad & fluoride-cradle detection
    toward known residue positions from the 3R3U canonical mapping.
    Negative values are distance bonuses (subtracted from the candidate distance;
    lower effective score = better candidate).

    NB: the bias is a *soft preference*, not a hard lock. It is sized as a few Å —
    comparable to inter-residue spacing — so a clearly closer geometric candidate
    can still win over the sequence-aligned hint. SMART_LOCK_NUC_MAX_DIST and the
    Nuc–Base sanity check remain the hard safety nets.
    """
    SMART_LOCK_BIAS_DIST: float    =   -5.0  # Å  soft bonus when residue matches mapped-hint position
    SMART_LOCK_CHAIN_BIAS: float   =   -2.5  # Å  additional soft bonus for same-chain residue
    SMART_LOCK_RESNUM_WINDOW: int  =    15   # residue-number window (±N) around each hint
    SMART_LOCK_NUC_MAX_DIST: float =   20.0  # Å  hard cutoff — no nucleophile candidate beyond this
    SMART_LOCK_OD_FALLBACK_DIST: float  =  6.0   # Å  Oδ–Cα search radius for nucleophile Oδ fallback
    SMART_LOCK_OD_FALLBACK_ANGLE: float = 120.0  # °  minimum backside approach angle (anti-F) in Oδ fallback
    SMART_LOCK_MAPPING_SANITY_DIST: float = 15.0  # Å  max plausible Nuc–Base dist; beyond = mapping error

    # ===============================================================================
    # SECTION 8: CATALYTIC TIER CLASSIFICATION  (Steps 02, 03, 04, 06)
    # ===============================================================================
    """
    The tier cascade is evaluated in strict descending order.
    The top tier is index 0 in TIER_ORDER, the lowest is the last element.
    A structure is assigned the highest tier whose ALL criteria are met.
    Criteria: nucleophile distance + attack angle + triad distances + mech score.
    """

    # --- Dynamic Tier Taxonomy ---
    """
    Alphanumeric sorting dictates the strict hierarchy.
    Tier_1A automatically sorts to index 0 (Top Tier), Tier_5_Decoy to the bottom.
    """
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


    # -------------------------------------------------------------------------------
    # Step 8.1: Nucleophile–C distance thresholds (Å, upper bound)
    # -------------------------------------------------------------------------------
    TIER_NUC_DIST: dict = field(default_factory=lambda: {
        "Tier_1A": 3.0,   # tight pre-reactive geometry; 2.7 Å unrealistically tight for Boltz-2 ground-state (true SN2 TS ~ 2.0–2.3 Å)
        "Tier_1B": 3.2,   # excellent pre-reactive geometry
        "Tier_2A":    3.2,   # = NAC_DIST_STRICT
        "Tier_2B":    3.8,   # = NAC_DIST_RELAXED
        "Tier_3":      4.2,   # productive but not pre-reactive
        "Tier_4":      8.0,   # pocket-bound, geometrically unproductive
    })

    # -------------------------------------------------------------------------------
    # Step 8.2: Attack angle minimum thresholds (°, lower bound)
    # -------------------------------------------------------------------------------
    TIER_ANGLE_MIN: dict = field(default_factory=lambda: {
        "Tier_1A": 170.0,   # near-ideal linear SN2 trajectory (within ~10° of the 180° Walden-inversion TS). Surfaces a chemotype-clean set of native-substrate (FA/DFA) elite candidates on diverse enzymes. Perfluoroalkyl (β-fluorinated) chains are kept out of Tier_1A by the β-feasibility + containment penalties, not by this angle gate. An α-CF3 substrate (TFA, scissile C–F 127.5) is barred from Tier_1A by the bond-strength ceiling TIER_ELITE_BDE_MAX (§8.5) at ANY angle — geometry cannot repeal thermochemistry — and is capped at Tier_1B as a discovery lead for Step-07 QM/MM.
        "Tier_1B": 165.0,
        "Tier_2A":    155.0,   # = NAC_ANGLE_STRICT
        "Tier_2B":    145.0,   # = NAC_ANGLE_RELAXED
    })

    # -------------------------------------------------------------------------------
    # Step 8.3: Catalytic triad distance thresholds (Å, upper bound)
    # -------------------------------------------------------------------------------
    TIER_NB_MAX: dict = field(default_factory=lambda: {
        "Tier_1A": 3.5,   # tightest triad — Nuc–Base ≤ 3.5 Å
        "Tier_1B": 4.0,
        "Tier_2A":    5.0,
        "Tier_2B":    6.0,   # loosest still-connected proton relay; beyond this the base is dissociated → not a degrader
    })
    TIER_BA_MAX: dict = field(default_factory=lambda: {
        "Tier_1A": 4.5,   # Base–Acid ≤ 4.5 Å
        "Tier_1B": 5.0,
        "Tier_2A":    6.0,
        "Tier_2B":    7.0,   # loosest still-connected Base–Acid pair; beyond this the acid is dissociated → not a degrader
    })

    # -------------------------------------------------------------------------------
    # Step 8.4: Mechanistic score thresholds (0–1, lower bound)
    # -------------------------------------------------------------------------------
    """
    Calibrated against the six control jobs (DeHa4 / 3R3U × FA / DFA / TFA) under the
    holistic mech_score (anchors 0.70 + graded Šidák-corrected SN2 angle 0.30, §5.1).
    A full-machinery mono-F pose scores 0.70 + 0.30·(1-p1), p1=(1-cos(180-angle))/2, so
    at each tier's angle floor it reaches: Tier_1A(170°)≈0.99, Tier_1B(165°)≈0.975,
    Tier_2A(155°)≈0.958. The
    minima below sit under those so the gate stays meaningful (binding only when an
    anchor is missing) while the control tiers are reproduced.
    """
    TIER_MECH_MIN: dict = field(default_factory=lambda: {
        "Tier_1A": 0.85,   # ladder floor; the elite tier is further refined by the coupled machinery gate below
        "Tier_1B": 0.85,   # complete anchors + strong (≥165°) trajectory
        "Tier_2A":    0.70,   # tolerates one missing anchor with a good (≥155°) trajectory
    })

    """
    Coupled elite-machinery gate (Tier_1A refinement, Step 02 — downgrade-only).
    A pose holds Tier_1A only if it has either (a) complete catalytic machinery and an
    open SN2 backside (mech ≥ MECH_ELITE_HI), OR (b) near-complete machinery
    (mech ≥ MECH_ELITE_LO) redeemed by a crystal-exact catalytic constellation
    (catalytic_constellation_score ≥ MECH_ELITE_CONSTELLATION, RMSD ≲ 0.25 Å). Clause (b)
    is the principled route by which a slightly backside-occluded native substrate such
    as α-CF3 trifluoroacetate can still register as elite — but only its single most
    crystal-perfect pose, not every mediocre one. This is geometry/machinery only; no
    substrate-class label is used.
    """
    MECH_ELITE_HI: float            = 0.90   # mech at/above which the backside is open enough for elite on its own
    MECH_ELITE_LO: float            = 0.85   # mech floor for the constellation-compensated elite route
    MECH_ELITE_CONSTELLATION: float = 0.74   # constellation that compensates a mech in [LO, HI) for Tier_1A (RMSD ≲ 0.35 Å, crystal-grade)

    """
    Elite bond-strength ceiling (Tier_1A, downgrade-only). A scissile C–F above TIER_ELITE_BDE_MAX
    cannot hold the elite tier at any attack angle: the strength of the bond being broken is a
    property of the bond, not of the approach geometry, so no near-linear trajectory lowers it. Such
    a pose is capped at Tier_1B and remains a discovery lead for Step-07 QM/MM to adjudicate.

    The ceiling sits above the α-CF3 class (127.5 kcal/mol), so an α-CF3 substrate such as
    trifluoroacetate is eligible for Tier_1A when its geometry earns it. TFA is the hardest substrate
    the hydrolytic SN2 mechanism can plausibly reach, and the pipeline surfaces it as a lead rather
    than ruling it out a priori. The bond strength is not waived: it is carried as a flat, non-fading
    penalty on mechanistic_score_effective (CHEM_PEN_W_BDE × the excess over SCISSILE_CF_BDE_MAX),
    which costs an α-CF3 pose ~0.11 of mech score. Only the most crystal-perfect α-CF3 pose therefore
    survives the coupled elite gate — one lead, not a chemotype.

    The ceiling is chemically specific: only an α-CF3 carbon reaches 127.5. The long perfluoro
    carboxylates present an α-CF2 (119.5) and are held below the elite tier by the β-fluorine and
    containment penalties, not by this bond-strength ceiling.
    """
    TIER_ELITE_BDE_MAX: float       = 128.0  # kcal/mol; a scissile C–F above this cannot hold Tier_1A at any attack angle

    """
    Top-tier confidence guard. The catalytic machinery for Tier_1A is enforced by the tier
    ladder itself — all eight catalytic residues mapped to the correct type at catalytic
    distances with the fluoride cradle engaging the leaving F (Chan 2011). Global sequence
    identity is therefore NOT a tier gate: a distant homolog with a valid active site keeps
    its tier. The residual risk is an unconfident predicted fold.

    The confidence that matters is LOCAL. Every geometric quantity the elite tier rests on — the
    SN2 angle, the nucleophile distance, the triad relay, the fluoride cradle — is measured on the
    eight catalytic residues and on nothing else. A GLOBAL confidence score averages those eight
    residues with hundreds of loop and surface residues the tier decision never touches: a protein
    with a crisply resolved active site and disordered termini is punished for the termini, while a
    globally confident fold with a smeared active site sails through. The gate therefore keys on
    active_site_plddt — the mean Boltz pLDDT over exactly the eight mapped residues.

    TIER_ELITE_AS_PLDDT_MIN = 90 is the conventional 'very high confidence' pLDDT band (the
    AlphaFold/Boltz interpretation scale), not a corpus-fitted number. Measured on this corpus, the
    active-site pLDDT of the current Tier_1A set runs 92.1–98.9 (median 97.9) against a corpus 1st
    percentile of 90.1, so the floor disqualifies a genuinely smeared active site without pruning
    the well-resolved elite set. A Tier_1A hit below it is demoted ONE notch to Tier_1B (still elite
    geometry, still reviewed); the geometric tier is retained in its own column.

    TIER_ELITE_CONF_MIN stays as the FALLBACK for a pose that carries no per-residue confidence, so
    such a pose is still checked rather than waved through.
    """
    TIER_ELITE_AS_PLDDT_MIN: float = 90.0   # active-site pLDDT (0–100) floor for the Tier_1A top label
    TIER_ELITE_CONF_MIN: float = 0.85       # global Boltz confidence floor — fallback when active_site_plddt is absent

    # -------------------------------------------------------------------------------
    # Step 8.4a: Catalytic-constellation (Criterion B) per-tier floors (0–1, lower bound)
    # -------------------------------------------------------------------------------
    """
    TWO active-site criteria. Criterion A — active_site_integrity — asks only whether the eight
    catalytic residues are PRESENT and correctly typed (sequence/identity). Criterion B —
    catalytic_constellation_score — asks whether they are GEOMETRICALLY ASSEMBLED like the 3R3U
    crystal, derived from the all-eight-residue Cα-superposition RMSD (B = 1/(1+Active_Site_RMSD),
    0–1). A alone is geometry-blind (eight correctly-typed residues scattered across the fold
    still score 1.0); B closes that gap.

    B CAPS the tier — downgrade only, never a promoter: the SN2 attack trajectory remains the
    positive catalytic signal, so a good constellation cannot manufacture a top tier on weak
    geometry. A pose whose B is below its assigned tier's floor is moved to the highest tier it
    qualifies for; B below the Tier_2B floor → TIER_CONSTELLATION_FLOOR_TIER (B>0) or
    Tier_5_Decoy (B==0, constellation unmeasurable), with is_degrader cleared. geometric_tier
    preserves the pre-demotion call. Calibrated on the control run: FA/DFA and the near-ideal
    Tier_1A poses (B 0.62–0.81) are retained; only genuinely mis-assembled sites are downgraded.
    """
    TIER_CONSTELLATION_MIN: dict = field(default_factory=lambda: {
        "Tier_1A": 0.55,   # RMSD ≲ 0.8 Å — eight-residue constellation floor (the coupled elite-machinery gate + pocket-fit cap do the elite separation; this floor admits all three controls FA/DFA/TFA whose best poses sit at B≈0.665–0.83)
        "Tier_1B": 0.45,   # RMSD ≲ 1.2 Å
        "Tier_2A": 0.35,   # RMSD ≲ 1.9 Å
        "Tier_2B": 0.25,   # RMSD ≲ 3.0 Å
    })
    TIER_CONSTELLATION_FLOOR_TIER: str = "Tier_4"   # B>0 but below every degrader floor → here (decoy if B==0)

    '''
    Pocket-fit elite ceiling (Step 02). FAcD's catalytic pocket evolved for a
    2-carbon haloacetate; a ligand whose longest interatomic extent exceeds this
    cannot occupy the active site productively even with a near-attack angle, so
    it is barred from Tier_1A (downgraded one notch, never below — discovery-open
    in the lower degrader tiers). Controls FA/DFA/TFA ≈ 3.5 Å; short PFCAs ≥ 5.7 Å.
    '''
    TIER_1A_MAX_LIGAND_EXTENT: float = 5.0   # Å — threshold defining the controls-only validation subset (ligands ≤ this, ≈2-carbon haloacetate controls ~3.5 Å; long-chain PFAS exceed it). The tier ladder itself is size-agnostic.

    # -------------------------------------------------------------------------------
    # Step 8.4b: Substrate / inhibitor classification (Figs 19, 25)
    # -------------------------------------------------------------------------------
    # Derived from the tier gates so every figure classifies identically.
    SUBSTRATE_ANGLE_MIN: float = 165.0   # = TIER_ANGLE_MIN['Tier_1B']: SN2 ≥ this → substrate geometry
    INHIBITOR_ANGLE_MAX: float = 145.0   # = TIER_ANGLE_MIN['Tier_2B']:   SN2 < this → potential inhibitor
    SUBSTRATE_CONF_MIN:  float = 0.75    # Boltz confidence ≥ this → AI-confident pose

    # -------------------------------------------------------------------------------
    # Step 8.4c: Confidence-vs-tier conflict thresholds (Fig 20 / Hidden-Gem rescue)
    # -------------------------------------------------------------------------------
    """
    Used by 03_Validation_Figures analyse_conflicts() to split structures into
    Consensus High / Hidden Gem / Consensus Low / Decoy. Kept separate from
    SUBSTRATE_CONF_MIN: this axis is AI-vs-physics agreement, not substrate geometry.
    """
    CONFLICT_CONF_HIGH: float = 0.70   # Boltz confidence ≥ this → AI-confident (Consensus High; high-tier below = no rescue)
    CONFLICT_CONF_LOW:  float = 0.60   # Boltz confidence < this → AI-doubtful (low-tier below = Consensus Low)

    # -------------------------------------------------------------------------------
    # Step 8.5: Tier scoring & ranking weights
    # -------------------------------------------------------------------------------
    """
    TIER_SCORE       — additive points for composite scoring
    TIER_RANK        — integer rank for quality-sorted operations
    TIER_SORT_WEIGHT — weight for multi-key DataFrame sorting
    """
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

    # -------------------------------------------------------------------------------
    # Step 8.6: Tier display colours (Okabe-Ito colourblind-safe palette)
    # -------------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------------
    # Step 8.7: Conflict category colours
    # -------------------------------------------------------------------------------
    CONFLICT_COLOUR: dict = field(default_factory=lambda: {
        "Consensus High": "#0072B2",
        "Hidden Gem":     "#CC79A7",
        "Decoy":          "#D55E00",
        "Consensus Low":  "#F0E442",
        "Ambiguous":      "#BBBBBB",
    })

    # -------------------------------------------------------------------------------
    # Step 8.8: Alignment grade colours (worst I = brown → best A = green)
    # -------------------------------------------------------------------------------
    # Keys match single-letter ALIGN_GRADE_LABELS; tuple = (lo%, hi%, hex)
    GRADE_COLOUR: dict = field(default_factory=lambda: {
        "A": "#2E7D52", "B": "#559B6A", "C": "#88B88A",
        "D": "#C8C87A", "E": "#C9A84C", "F": "#C07C40",
        "G": "#B05535", "H": "#8B3A2A", "I": "#6B2A1F",
    })

    # Grade bands as list of (lo, hi, hex, letter) for Fig 02 histogram shading
    GRADE_BANDS: list = field(default_factory=lambda: [
        (90, 100, "#2E7D52", "A"), (80, 90, "#559B6A", "B"),
        (70, 80, "#88B88A", "C"), (60, 70, "#C8C87A", "D"),
        (50, 60, "#C9A84C", "E"), (40, 50, "#C07C40", "F"),
        (30, 40, "#B05535", "G"), (20, 30, "#8B3A2A", "H"),
        (0,  20, "#6B2A1F", "I"),
    ])

    # Grade colours keyed by full label string (used in Fig 03 stacked bars)
    GRADE_COLOUR_FULL: dict = field(default_factory=lambda: {
        "I (<20%)":   "#6B2A1F", "H (20–30%)": "#8B3A2A",
        "G (30–40%)": "#B05535", "F (40–50%)": "#C07C40",
        "E (50–60%)": "#C9A84C", "D (60–70%)": "#C8C87A",
        "C (70–80%)": "#88B88A", "B (80–90%)": "#559B6A",
        "A (≥90%)":   "#2E7D52",
    })

    # -------------------------------------------------------------------------------
    # Step 8.8b: Active-site residue display order + role-group colours (Fig 01)
    # -------------------------------------------------------------------------------
    # Mechanistic role sequence for the active-site mapping-coverage figure:
    # nucleophile → acid/base catalysis → carboxylate clamp → fluoride pocket.
    ACTIVE_SITE_ROLE_ORDER: list = field(default_factory=lambda: [
        "Nuc", "Acid", "Base", "Carb1", "Carb2", "Stab_H", "Stab_W", "Stab_Y",
    ])
    # Residue key → functional group label (residues of the same role share a colour).
    ACTIVE_SITE_ROLE_GROUP: dict = field(default_factory=lambda: {
        "Nuc":    "Nucleophile",
        "Acid":   "Acid/base catalysis", "Base": "Acid/base catalysis",
        "Carb1":  "Carboxylate clamp",   "Carb2": "Carboxylate clamp",
        "Stab_H": "Fluoride pocket", "Stab_W": "Fluoride pocket", "Stab_Y": "Fluoride pocket",
    })
    # Role-group colours (Okabe–Ito colourblind-safe palette, matching TIER_COLOUR).
    ACTIVE_SITE_ROLE_GROUP_COLOUR: dict = field(default_factory=lambda: {
        "Nucleophile":          "#0072B2",
        "Acid/base catalysis":  "#D55E00",
        "Carboxylate clamp":    "#009E73",
        "Fluoride pocket":      "#CC79A7",
    })

    # -------------------------------------------------------------------------------
    # Step 8.9: Mechanistic outcome colours (Fig 25b)
    # -------------------------------------------------------------------------------
    """
    Deliberately a different hue family from the tier palette so the two
    stacked bars in Fig 25b never read as the same encoding.
    """
    OUTCOME_COLOUR: dict = field(default_factory=lambda: {
        "Substrate":           "#0E7C7B",   # teal
        "Borderline":          "#C9A227",   # mustard
        "Reactive (low conf)": "#8E44AD",   # purple
        "Non-reactive":        "#95A5A6",   # grey
        "Potential Inhibitor": "#6E2C00",   # brown-maroon
    })

    # -------------------------------------------------------------------------------
    # Step 8.9b: PFAS chain-length bin colours (Fig 25b master bars)
    # -------------------------------------------------------------------------------
    PFAS_SIZE_BIN_COLOUR: list = field(default_factory=lambda:
        ["#6A51A3", "#2171B5", "#238B45", "#D94801", "#A50F15"])

    # Best → worst 5-stop gradient (nucleophile-distance & SN2-angle bins).
    SANKEY_GRAD5: list = field(default_factory=lambda:
        ["#1B7837", "#5AAE61", "#D9EF8B", "#FDAE61", "#D73027"])
    SANKEY_CLAMP_COLOUR: dict = field(default_factory=lambda: {
        "Clamp intact": "#1B9E77", "Clamp broken": "#D55E00"})
    SANKEY_STAB_COLOUR: dict = field(default_factory=lambda: {
        "Stabilised": "#1B9E77", "Unstabilised": "#BDBDBD"})
    SANKEY_TRIAD_COLOUR: dict = field(default_factory=lambda: {
        "Triad tight": "#1B9E77", "Triad moderate": "#E6A817", "Triad loose": "#D55E00"})
    SANKEY_MECH_GRAD: list = field(default_factory=lambda:
        ["#D73027", "#FDAE61", "#A6D96A", "#1B7837"])   # worst→best, 4 bins

    # ===============================================================================
    # SECTION 9: MD TRAJECTORY ANALYSIS  (Step 07)
    # ===============================================================================

    # -------------------------------------------------------------------------------
    # Step 9.1: Solvent / water residue names
    # -------------------------------------------------------------------------------
    SOLVENT_RESTYPES: tuple = (
        "SPC", "TIP3P", "TIP4P", "SPCE", "OPC", "HOH", "WAT", "SOL",
        "T3P",   # Desmond internal name for TIP3P
    )

    # -------------------------------------------------------------------------------
    # Step 9.2: WaterMap integration
    # -------------------------------------------------------------------------------
    # Abel et al. (2008): site radius 5.0 Å captures active-site hydration shell.
    WATERMAP_SITE_RADIUS: float     = 5.0   # Å  radius to include WaterMap hydration sites
    WATERMAP_BLOCKADE_RADIUS: float = 3.0   # Å  cylinder radius for nucleophile-runway blockade check
    WATERMAP_MATCH_RADIUS: float    = 1.5   # Å  MD water ↔ WaterMap site assignment distance
    """
    Fold-integrity guard on the frame→reference superposition. The WaterMap sites are computed in
    one reference frame and carried into each trajectory frame by a rigid Kabsch transform. A rigid
    transform exists between ANY two point sets, so a frame whose fold has drifted, unfolded, or
    whose atom correspondence has broken still returns a rotation — and the sites it carries land
    at arbitrary positions. Above this Cα RMSD the frame's WaterMap term is withheld rather than
    trusted. 3.0 Å is a loose fold-identity bound: normal thermal breathing of a folded protein
    stays well below it, so only genuinely broken frames are rejected.
    """
    MD_FOLD_RMSD_MAX: float         = 3.0   # Å  Cα RMSD above which a frame's WaterMap mapping is discarded
    """
    Counter-ion capping of the reactive centre. A PFAS carboxylate pairs strongly with Na⁺, and the
    FAcD active site is an anion trap (Asp nucleophile, Asp acid, substrate carboxylate). System
    Builder's ion-exclusion region keeps counter-ions out at BUILD time, but one can diffuse in during
    the run, coordinate the nucleophile Oδ or the ligand head, and screen the charge the SN2 depends
    on — while the NAC distance and angle still look perfectly productive.

    A frame is capped when a cation sits within this distance of either. 3.0 Å is the inner-sphere
    bound: direct Na⁺–carboxylate-O coordination is ~2.4 Å, whereas a solvent-separated ion pair sits
    beyond ~4.5 Å and does not screen the site. Capped frames are flagged in the per-frame CSV,
    counted in the stats (Frames_Cation_Capped / Cation_Capped_Pct), and barred from QM/MM frame
    selection — a QM region containing a Na⁺ on the nucleophile computes that ion pair's barrier,
    not the enzyme's.
    """
    CATION_CAP_DIST: float          = 3.0   # Å  cation ↔ nucleophile Oδ / ligand carboxylate O: inner-sphere coordination

    """
    NPT equilibration verification (Step 07, from the Desmond <job>.ene stream).

    Every quantity Step 07 reports — NAC occupancy, strict-NAC dwell, MM-GBSA — is an equilibrium
    average. An average taken over a system that is still relaxing is not an average of anything, and
    nothing else in the pipeline would notice: a drifting box produces perfectly well-formed NAC
    statistics. The box VOLUME is the slow coordinate (the barostat is still working on it long after
    the thermostat has settled), so equilibration is declared from V and T is a separate thermostat
    sanity check. Frames before the equilibration time are excluded from the sampled statistics.
    """
    MD_EQUIL_SKIP_FRAC: float          = 0.20   # leading fraction of the run treated as relaxation when forming the production-window reference
    MD_EQUIL_BLOCKS: int               = 200    # volume is block-averaged into this many blocks before the settled test — an equilibrated NPT box still spikes instantaneously (measured: 0.02 % of points exceed 1 % of the mean), so a point-wise test would reject a settled trajectory
    MD_EQUIL_V_TOL_PCT: float          = 1.0    # %  block-mean box volume within this of the production mean counts as settled
    MD_EQUIL_V_DRIFT_MAX_PCT_NS: float = 0.05   # %/ns  residual volume drift above this = not equilibrated (a box still shrinking or swelling)
    MD_EQUIL_TARGET_T: float           = 300.0  # K  the thermostat set point
    MD_EQUIL_T_TOL_K: float            = 3.0    # K  mean temperature may deviate from the set point by at most this

    # -------------------------------------------------------------------------------
    # Step 9.3: Frame scoring weights (QM/MM frame selection only)
    # -------------------------------------------------------------------------------
    """
    These weights rank MD frames to select the single best input for QSite.
    They do NOT affect NAC counts or catalytic viability percentages.
    """
    SCORE_DIST_WEIGHT: float     = 100.0   # per Å below NAC_DIST_RELAXED
    SCORE_ANGLE_WEIGHT: float    =   5.0   # per degree above NAC_ANGLE_RELAXED
    SCORE_BLOCKADE_WEIGHT: float =  50.0   # per water blockade unit
    SCORE_WATERMAP_WEIGHT: float =  10.0   # per kcal mol⁻¹ WaterMap dG unit

    # -------------------------------------------------------------------------------
    # Step 9.4: Catalytic viability display thresholds (%)
    # -------------------------------------------------------------------------------
    """
    Used to colour-code console output and figures in Step 07.
    Boltz-2 predicted structures (not crystal structures) typically show
    lower NAC populations (0.01–2%) owing to the Boltz-2 starting geometry
    not being pre-optimised for the reactive SN2 trajectory.
    """
    VIABILITY_PASS_THRESHOLD: float = 0.05   # % — below this → "NAC FAIL" (red)
    VIABILITY_HIGH_THRESHOLD: float = 5.0    # % — above this → "NAC PASS – High" (green)

    # -------------------------------------------------------------------------------
    # Step 9.5: Force field — Desmond MD
    # -------------------------------------------------------------------------------
    """
    System solvation, ion neutralisation, and NPT production MD all use the
    OPLS4 force field (Roos et al. 2019; Lu et al. 2021, JCTC) as configured
    in the Desmond .msj job file generated by the user via Maestro.
    The Python script (07_MD_QMMM_Defluorination_FAcDs.py) reads the
    completed trajectory — it does not control force-field selection.
    OPLS4 force field: Lu et al. (2021); Roos et al. (2019) — see header §13.
    """

    # ===============================================================================
    # SECTION 10: QM/MM EXTRACTION — QSite  (Step 07)
    # ===============================================================================
    """
    Level of theory: B3LYP / 6-31+G(d,p) — Becke (1993) + Lee, Yang & Parr (1988)
    hybrid functional; Rosta et al. (2006) QM/MM free-energy methodology;
    Murphy et al. (2000) QSite implementation.
    NB: QSite frozen-orbital QM/MM cuts (used to place the catalytic sidechains in
    the QM region) only support a limited set of plain functionals — B3LYP/HF.
    Meta-GGA hybrids (e.g. M06-2X) and dispersion-corrected variants (B3LYP-D3)
    are rejected by QSite with frozen cuts, hence B3LYP for the residue-selective
    QM/MM coordinate scan.
    """
    """
    DFT functional. Plain B3LYP carries no dispersion term, and that IS a real limitation here: the
    substrate is polyfluorinated and the C–F···π contacts holding it against the Trp/Tyr cradle are
    dispersion-bound, so the barrier is computed without them.

    It cannot be fixed by switching functional. QSite's frozen-orbital cuts — required for a
    residue-selective QM region — reject dispersion-corrected and meta-GGA functionals outright;
    Jaguar aborts with "ERROR 5029: Disallowed QM Method for QSite with Frozen Orbital Cuts:
    DFT(b3lyp-d3)" (verified against a real QM/MM input, 13 July 2026). Dispersion would require
    abandoning the frozen-cut QM region, which is a larger change than it buys. The limitation is
    declared in the paper rather than hidden.
    """
    QSITE_FUNCTIONAL: str   = "b3lyp"        # DFT functional (the only family QSite frozen cuts accept)
    QSITE_BASIS_SET: str    = "6-31+G(d,p)"  # basis for single-point energies (the relaxed scan
                                             # runs QSITE_SCAN_BASIS, which is non-diffuse for SCF
                                             # stability — see §10.2)
    QSITE_CHARGE: int       = -1             # default QM region charge (anionic carboxylate/sulfonate PFAS)
    """
    Neutral ligands (alcohols, non-ionised at pH 8) override the default charge.
    Keys are substrings of the job-name ligand suffix (case-insensitive match).
    """
    LIGAND_QM_CHARGES: dict = field(default_factory=lambda: {
        "FTOH": 0,   # fluorotelomer alcohols (6:2-FTOH, 8:2-FTOH) — neutral at pH 8
    })
    QSITE_MULT: int         = 1              # spin multiplicity (closed-shell singlet)
    """
    Force field for the MM half of the QM/MM Hamiltonian, emitted into &mmkey.

    The keyword is qsite_ff and it takes a STRING. Confirmed against Schrodinger's own QSite driver
    (mmshare .../common/qsite_binding_energies.py), which offers exactly two force fields —

        parser.add_argument('-ffield', choices=['OPLS_2005', 'OPLS3e'], default='OPLS_2005')
        if cmd_args.ffield == 'OPLS3e':
            qs_mmkey_dict['qsite_ff'] = 'opls3e'

    — so QSite's classical region CANNOT run OPLS4. The best available is OPLS3e, and the DEFAULT is
    OPLS_2005: left unset, the MM region would score an OPLS4 trajectory on a force field two
    generations older, and the barrier would carry an energetic discontinuity that is not chemistry.
    PrepWizard was found doing exactly this (defaulting to OPLS_2005 while everything downstream ran
    OPLS4), so the assumption is worth stating rather than trusting.

    This is a DECLARED LIMITATION, not a fix: OPLS3e is closer to the trajectory's OPLS4 than
    OPLS_2005 is, but it is not the same force field. The QM region — where the bond actually breaks —
    is unaffected; the mismatch sits in the classical environment around it.
    """
    QSITE_MM_FF: str        = "qsite_ff=opls3e"   # &mmkey MM force field; OPLS4 is not offered by QSite
    '''
    Solvation: the extracted frame carries its explicit TIP3P water box in the MM
    region, so no implicit-solvation keyword is emitted (an implicit model would
    double-count). Publication-grade activation energies on top candidates use
    explicit-solvent QM/MM-FEP or thermodynamic integration downstream.
    '''
    QSITE_QM_INCLUDE_CRADLE: bool = True     # include the full fluoride cradle (Trp156 + Tyr217 alongside His155) in the QM region — the departing F⁻ is a hard base whose charge-transfer/polarisation with the aromatic cradle is poorly captured by MM point charges; only the handful of Tier_1A candidates reach QSite, so the added DFT cost is bounded. Set False for the cheaper His155-only QM region.
    QSITE_SCAN_START: float = 3.5            # Å  scan start (pre-reaction approach; 3.5 → 1.3 Å = full SN2 coordinate)
    QSITE_SCAN_STEP: float  = -0.1           # Å  step per point (negative = bond compression)
    QSITE_SCAN_NSTEPS: int  = 23             # points total → covers 3.5 → 1.3 Å (last point: 3.5 + −0.1×22 = 1.3 Å)
    """
    Execution of the generated QSite jobs from Step 07. When True, Step 07
    launches `$SCHRODINGER/qsite` on each freshly extracted frame, writing all
    output inside a per-job folder; jobs whose folder already exists are skipped
    (idempotent, mirroring the PDB-preparation cache). QM/MM relaxed scans are
    expensive — disable with --no-run-qsite to only write the .in/.mae inputs.
    """
    QSITE_IMPVERSION: str   = "huge"         # Jaguar &gen impversion (memory/architecture tier; tune per cluster). igeopt=1 (relaxed scan) and mmqm=1 (QM/MM) are required mode flags for this calculation and stay fixed in the writer.
    QSITE_RUN: bool         = True
    QSITE_PROCS: int        = 10             # CPUs per QSite job (qsite -PARALLEL). Jaguar's SCF
                                            # scales usefully to ~8-16 cores; 3 concurrent jobs at 10
                                            # leaves headroom on a 32-core box
    QSITE_PROGRESS_INTERVAL_SEC: int = 20    # heartbeat cadence while a QSite job runs (live progress, prevents "frozen" look)
    # --- Step 10.1: Post-scan QM/MM plotting / rendering (Step 07 figures) ---
    HARTREE_TO_KCAL: float  = 627.509474     # Eh → kcal mol⁻¹ (relative scan energies)
    QSITE_MAESTRO_RENDER: bool = True        # render reactant/TS/product QM-region images via headless Maestro (best-effort; degrades gracefully if no $SCHRODINGER/display)
    QSITE_RENDER_TIMEOUT_SEC: int = 600      # hard timeout for each headless Maestro render call

    # --- Step 10.2: Multi-frame QM/MM barrier (defensible ensemble, not a single-frame lower bound) ---
    QSITE_N_FRAMES: int = 3                  # number of top pre-organised NAC frames to run the QM/MM SN2 scan on; the reported ΔE‡ is min/mean/σ over them. 1 scans only the single best frame, which reports a lower bound rather than an ensemble
    QSITE_MAX_QM_RESIDUES: int = 8           # cap on catalytic residues in the QM region (nucleophile/base/acid/stab first, then nearest cradle). A very large QM region (e.g. 17 residues) inflates the electron count and makes molchg/electron-parity errors likely → Jaguar 'incorrect molecular charge' and every scan point skipped. 0 = no cap.
    """
    The QM/MM system around the reaction centre. The droplet is the MM shell the QSite job keeps
    around the QM region; the QM waters are the few molecules close enough to the reactive centre
    that leaving them classical would misdescribe the fluoride's first solvation shell.
    """
    """
    The MM solvation droplet kept around the QM region. The periodic water box is trimmed to this
    radius around the ligand, which leaves the droplet with a FREE SURFACE: QSite has no periodic
    boundary and exposes no frozen-shell or boundary-potential keyword (its `&mmkey` section takes
    MMIM keywords such as `use_nb_cutoff`; nothing in the QSite API, the Jaguar binaries or the
    shipped examples provides an atom freeze). The surface waters are therefore unrestrained during
    the relaxed scan.

    The radius is the control that matters. At 8 Å the free surface sits directly on the reaction
    centre's first and second solvation shells, where any splaying distorts the electrostatics the
    QM region feels. At 15 Å the surface is two solvation shells away from the ligand — far enough
    that its relaxation cannot reach the scissile bond over the short scan — and the cost is MM-only
    (a few thousand extra classical waters), which is negligible beside the DFT.
    """
    QSITE_DROPLET_RADIUS: float = 15.0       # Å  MM droplet retained around the QM region
    """
    The droplet's boundary, in three zones around the ligand. QSite reads them from a per-atom
    property on the structure (`i_i_constraint`), and Jaguar reports back exactly how many atoms it
    took as frozen and as constrained — verified against a live QM/MM job (13 July 2026), where the
    property values map as 0 = free, 1 = FROZEN, 2 = CONSTRAINED.

    Without it the droplet has a free surface: the periodic box is cut to a finite ball of water,
    nothing holds the outermost molecules, and they relax into vacuum during the scan — distorting
    the electrostatics the QM region sits in. Freezing the outer shell removes that surface, while
    the buffer stays restrained (not rigid) so the solvation shell around the reaction centre can
    still respond to the reaction.
    """
    QSITE_FREE_RADIUS: float   = 5.0         # Å  ≤ this from the ligand: fully mobile
    QSITE_BUFFER_RADIUS: float = 8.0         # Å  ≤ this: restrained; beyond: frozen (the surface)
    QSITE_QM_WATER_RADIUS: float = 3.5       # Å  a water within this of the reactive centre goes QM
    QSITE_QM_WATER_MAX: int = 3              # cap on QM waters (each adds electrons to the SCF)
    """
    SCF accuracy grid for the relaxed scan. 1 is the fast/robust grid; 2-3 densify the integration
    grid, which matters for an anionic leaving group (F⁻ has a diffuse, slowly-decaying density) at
    a cost the scan can ill afford — each point is already ~3 h. Raise it for a final single-point
    re-evaluation of the barrier rather than for the scan itself.
    """
    QSITE_SCF_IACC: int = 1
    QSITE_SCAN_BASIS: str = "6-31G**"        # basis for the relaxed scan points: the scan is many
                                             # SCF cycles, so it runs a cheaper basis than the
                                             # single-point QSITE_BASIS_SET
    SOLVENT_SPHERE_SAMPLE_FRAMES: int = 12   # frames sampled when measuring the solvent sphere

    # --- Step 10.3: Defluorination verdict — the concrete "does it defluorinate?" gate (Step 07) ---
    """
    Binding (MM-GBSA) proves a Michaelis complex, not turnover. A candidate is called
    defluorination-competent only when it (i) persists in a STRICT near-attack
    conformation for a real dwell, (ii) surmounts a QM/MM SN2 barrier at body
    temperature, and (iii) the SN2 product (F⁻ displaced onto the Asp nucleophile) is
    not uphill. Defluorination_Propensity fuses persistence and barrier into a single
    kcat-like rate proxy:  P(strict-NAC) · exp(−ΔE‡ / RT)  (RT from GAS_CONSTANT_KCAL ×
    MMGBSA_TEMPERATURE_K). Rank candidates by that proxy, not by ΔG_bind.
    """
    DEFLUOR_STRICT_VIABILITY_MIN_PCT: float = 1.0    # Xs — min % of pocket-bound frames in the STRICT NAC
    DEFLUOR_DWELL_MIN_NS: float             = 1.0    # Y  — min longest CONTINUOUS strict-NAC residence (ns)
    DEFLUOR_BARRIER_MAX_KCAL: float         = 22.0   # Z  — max surmountable QM/MM SN2 barrier ΔE‡ (kcal/mol)
    DEFLUOR_DERXN_MAX_KCAL: float           = 0.0    # SN2 reaction energy ceiling — product must be ≤ reactant (ΔE_rxn ≤ this)

    # --- Step 10.4: Defluorination figure parameters (Step 07 reaction-profile,
    #     MM-GBSA-decomposition, and landscape plots) — all colours + thresholds
    #     here so no Step-07 figure hard-codes them (SSOT). ---
    DEFLUOR_FLUORIDE_CHARGE_MIN: float  = -1.2   # Mulliken-charge window low bound when parsing the departing-F charge
    DEFLUOR_FLUORIDE_CHARGE_MAX: float  = -0.4   # window high bound (rejects O / still-bonded F so only near-fluoride is tracked)
    DEFLUOR_FIG_COLOUR: dict = field(default_factory=lambda: {
        "pes":        "#1D4ED8",   # QM/MM potential-energy-surface line
        "ts":         "#DC2626",   # transition state (barrier peak)
        "product":    "#16A34A",   # product well
        "f_charge":   "#B45309",   # departing-fluoride Mulliken-charge curve
        "gate":       "#22C55E",   # competence-quadrant shading (landscape)
        "gate_line":  "#16A34A",   # gate threshold lines
        "gate_text":  "#15803D",   # gate annotation text
        "scatter":    "#3B82F6",   # landscape default marker (no propensity colour)
        "edge":       "#334155",   # marker / bar edge
        "zone_relaxed": "#74C476", # NAC dashboard: relaxed SN2 zone (distance × angle)
        "zone_strict":  "#006D2C", # NAC dashboard: strict SN2 zone
        "warhead":      "#D55E00", # NAC dashboard: nucleophile → warhead-C trace and its criterion line
        "tail":         "#0072B2", # NAC dashboard: cradle → ligand-F trace and its criterion line
        "kde":          "#111111", # NAC dashboard: density contours over the scatter
        "legend_edge":  "#E2E8F0", # legend frame
    })
    # Comparative residue-engagement heatmap (11_): colour ramp and the colour a missing residue
    # takes — a homolog that has no such residue must read as absent, never as a distance.
    ENGAGE_HEATMAP_CMAP: str = "RdYlGn_r"
    ENGAGE_HEATMAP_NAN: str  = "#E5E7EB"

    # ===============================================================================
    # SECTION 11: CANONICAL RESIDUE MAPPING — 3R3U Reference  (Step 07)
    # ===============================================================================
    """
    Reference sequence positions in DEHA4_CONTROL_SEQ numbering (the alignment
    reference, REF_SEQUENCE_STR = DEHA4_CONTROL_SEQ). Step 07 looks these up as
    keys in aln_dict, which is keyed by DEHA4 reference positions from the
    Full_Sequence_Alignment_Map column. (DEHA4 numbering Asp110/Asp134/His277;
    the 3R3U crystal carries the +3-shifted His280 etc.)
    Must stay in sync with REF_ACTIVE_SITE_MAP (§2.5).

    Catalytic triad = Asp110–His277–Asp134 (DEHA4 numbering; His280 in the 3R3U
    crystal of Chan et al. (2011) JACS 133:7461). Mechanism: Asp110 SN2 attack on
    Cα displaces F⁻ (→ glycolyl-enzyme ester); His277 activates water to hydrolyse
    the ester; Asp134 is the third triad residue that orients/polarises His277 (it
    does NOT protonate fluoride — F⁻ leaves stabilised as the anion by the pocket).
    Halide pocket (3 H-bonds to F⁻): His155, Trp156, Tyr217. Oxyanion hole:
    backbone amides of Phe40 + Arg111. Carboxylate clamp: Arg111, Arg114.
    """
    DREAM_TEAM_REFS: dict = field(default_factory=lambda: {
        "Nuc":    110,   # Asp110 — nucleophile; SN2 attack, forms covalent glycolyl-ester intermediate
        "Clamp1": 111,   # Arg111 — carboxylate clamp 1 (also oxyanion-hole backbone amide)
        "Clamp2": 114,   # Arg114 — carboxylate clamp 2
        "Acid":   134,   # Asp134 — third catalytic-triad residue; orients/polarises His277 base
        "Stab_H": 155,   # His155 — halide-pocket fluoride stabiliser (H-bond to F⁻)
        "Stab_W": 156,   # Trp156 — halide-pocket fluoride stabiliser (cradle)
        "Stab_Y": 217,   # Tyr217 — halide-pocket fluoride stabiliser / charge acceptor on SN2 axis (DEHA4 seq 217; 3R3U PDB Tyr219)
        "Base":   277,   # His277 — general base; activates hydrolytic water (DEHA4 seq 277; 3R3U PDB His280)
    })

    # ===============================================================================
    # SECTION 12: PROCESSING PARAMETERS
    # ===============================================================================

    @property
    def GLOBAL_MAX_WORKERS(self) -> int:
        import multiprocessing
        return max(1, multiprocessing.cpu_count() - self.PREP_CPU_RESERVE)
    PROC_BATCH_SIZE: int         = 2000   # max futures submitted at once — caps peak memory
    PROC_HEADER_CHECK_LINES: int =   30   # lines to scan when reading a PDB header tag
    # Step 01 merge — length window applied to the secondary (BLAST/UniProt/NCBI) FASTA only;
    # spans the FAcD single-domain α/β-hydrolase fold (~300 aa) with headroom for partial hits.
    MERGE_SECONDARY_LEN_MIN: int =  250   # aa — discard secondary sequences shorter than this
    MERGE_SECONDARY_LEN_MAX: int =  360   # aa — discard secondary sequences longer than this

    # ===============================================================================
    # SECTION 13: VISUALISATION PARAMETERS
    # ===============================================================================

    # -------------------------------------------------------------------------------
    # Step 13.1: Global rendering (Step 07)
    # -------------------------------------------------------------------------------
    VIS_IMG_WIDTH: int   = 2400   # px  export width  (publication-quality figure)
    VIS_IMG_HEIGHT: int  = 2400   # px  export height
    VIS_RAY_TRACE: bool  = True   # enable PyMOL ray-tracing for publication quality
    VIS_FIGURE_DPI: int  = 300    # dots per inch for publication figures (minimum 300)
    VIS_NAC_DIST_WARN_MAX: float = 5.0   # Å — dashboard caution band: NAC_DIST_RELAXED ≤ dist < this → orange, ≥ this → red
    """
    Typography and canvas — ONE definition for every figure the pipeline draws, applied through
    utils.apply_figure_style(). The point sizes carry the hierarchy on their own: axis labels are
    set in plain weight, since bolding every label emphasises nothing. Weight is spent only where
    it must be read against a filled bar (the value written inside it) or where a tick label
    doubles as a legend (a colour-coded role or tier name).
    """
    """
    Thresholds the DIAGNOSTIC figures judge by. They live here so a figure cannot disagree with the
    gate it is drawn to illustrate. The reactive-engagement panel calls a candidate 'ready' by the
    same nucleophile distance and attack angle used everywhere else, and a distance at or beyond the
    sentinel means the measurement is absent, not that the contact is long.
    """
    VIS_DIAG_READY_DIST_A: float   = 3.5    # Å  nucleophile-C distance at/below which a pose is 'ready'
    VIS_DIAG_ANGLE_MIN_DEG: float  = 150.0  # °  backside attack angle required alongside it
    VIS_DIAG_PROD_LO_A: float      = 2.5    # Å  productive window, lower edge (shading only)
    VIS_DIAG_PROD_HI_A: float      = 3.5    # Å  productive window, upper edge
    VIS_DIAG_DIST_SENTINEL_A: float = 20.0  # Å  at/above this a distance means 'not measured'
    VIS_DIAG_MIN_N_PER_BIN: int    = 30     # a bin below this many complexes is not plotted
    VIS_DIAG_MIN_N_FOR_BINNING: int = 50    # below this many rows the binned panel is not attempted
    SCORE_PCA_NEUTRAL: float       = 50.0   # the neutral PCA score assigned when PCA cannot run
                                            # (mid-scale of 0-100 — never a silent zero)
    VIS_FONT_FAMILY: tuple = ("Arial", "Helvetica", "DejaVu Sans")
    VIS_FONT_AXIS_LABEL: float  = 11.0   # x/y axis labels — plain weight
    VIS_FONT_TICK: float        = 9.0    # tick labels
    VIS_FONT_LEGEND: float      = 8.5    # legend entries
    VIS_FONT_ANNOT: float       = 7.5    # in-figure annotations (values on/inside bars)
    VIS_GRID_COLOUR: str        = "#EBEBEB"
    VIS_GRID_LINEWIDTH: float   = 0.6
    VIS_GRID_ALPHA: float       = 0.25   # the grid is a reading aid, never a mark competing with the data
    VIS_LEGEND_FRAME_ALPHA: float = 0.92
    """
    ── Step 03 ink: every mark colour in the validation figures that is not a tier, grade, conflict or
    role colour (those have their own dicts above). The candidate/tier colours say WHAT a mark is; the
    ink here says how it is drawn — medians, reference lines, annotation text, box edges, fills.

    Kept in one place for the reason the rest of this file exists: a restyle must be one edit, not a
    hunt through ten thousand lines of plotting code. Keys are by ROLE, not by hue, so changing the
    'reference line' colour changes every reference line and nothing else.
    """
    VIS_INK: dict = field(default_factory=lambda: {
        # neutrals — text, strokes, edges, fills
        "black":      "#000000",
        "near_black": "#111111",
        "outline":    "#222222",
        "dark":       "#333333",
        "soft":       "#444444",
        "muted":      "#555555",
        "mid":        "#666666",
        "grey":       "#777777",
        "ghost":      "#888888",
        "faint":      "#999999",
        "pale":       "#AAAAAA",
        "paler":      "#BBBBBB",
        "palest":     "#CCCCCC",
        "hairline":   "#DDDDDD",
        "tick":       "#DCDCDC",
        "wash":       "#E0E0E0",
        "mist":       "#ECECEC",
        "smoke":      "#EEEEEE",
        "grid":       "#EBEBEB",
        "panel":      "#FAFAFA",
        "canvas":     "#FBFBFB",
        "white":      "#FFFFFF",
        "slate":      "#9E9E9E",
        "silver":     "#BDBDBD",
        # slate-blue neutrals used for structural annotation frames
        "steel":      "#C9D2D9",
        "stone":      "#AAB7B8",
        "shadow":     "#7F8C8D",
        "charcoal":   "#566573",
        "ink_navy":   "#2E4053",
        "ink_deep":   "#2C3E50",
        "ink_pure":   "#1A1A1A",
    })
    """
    Qualitative accents. The base six are Okabe-Ito (colour-blind safe) and carry the same meaning
    wherever they appear; the rest are role-specific accents the validation figures need — the MD
    star, the twin-axis pair, the pass/fail marks.
    """
    VIS_ACCENT: dict = field(default_factory=lambda: {
        # Okabe-Ito qualitative base
        "blue":       "#0072B2",
        "vermillion": "#D55E00",
        "green":      "#009E73",
        "amber":      "#E69F00",
        "magenta":    "#CC79A7",
        "sky":        "#56B4E9",
        "yellow":     "#F0E442",
        # the MD-selected star (fill + stroke) — it must read instantly at a glance
        "star":       "#FFD400",
        "star_edge":  "#B8860B",
        # twin-axis pair: each axis carries a different quantity, so its ticks, label and gridlines
        # take the axis colour and cannot be misread as one another
        "axis_left":  "#2C6FAC",
        "axis_right": "#A06000",
        # pass / warn / fail
        "good":       "#1B7837",
        "warn":       "#B22222",
        "bad":        "#C0392B",
        "alert":      "#CC0000",
        "error":      "#FF6B6B",
    })
    """
    Sequential ramps and the pale tints used behind annotations. A ramp is ordered — light to dark —
    and is indexed, never picked from by name, so a figure cannot silently reorder its own scale.
    """
    VIS_RAMP: dict = field(default_factory=lambda: {
        "green":  ("#D4EFDF", "#A1D99B", "#238B45", "#1B7837", "#0B5345"),
        "orange": ("#FCE4D0", "#E08A3C", "#D94801", "#B84000", "#7E3E00"),
        "blue":   ("#D6EAF8", "#9DB8D2", "#4C8BC9", "#2C6FAC", "#1B3A5E"),
        "purple": ("#F9E4F2", "#B15FBF", "#8E44AD", "#6A51A3", "#3A1E4A"),
    })
    VIS_TINT: dict = field(default_factory=lambda: {
        "red":    "#FDECEA",   # behind a warning annotation
        "green":  "#D6F5EB",
        "amber":  "#FFF8E7",
        "blue":   "#CFE0EA",
        "cream":  "#FFF3CC",
    })
    """
    The three qualitative bands — strong / moderate / weak — used by every figure that shades a zone
    or labels a threshold region (confidence bands, engagement bands, mechanistic zones, tertiles).

    One green, one gold, one red, for ALL of them. The same band label must not be one green in the
    confidence figure and a slightly different green in the engagement figure: a reader who sees two
    greens is entitled to assume they mean two different things.
    """
    VIS_BAND: dict = field(default_factory=lambda: {
        "high":      "#007A50",
        "moderate":  "#8A6000",
        "low":       "#CC2222",
        "high_fill": "#D4EFDF",
        "mod_fill":  "#FFF3CC",
        "low_fill":  "#FADBD8",
    })
    """
    Interaction (bond) types, one colour each, wherever an interaction profile is stacked or split.
    Sourced here so the H-bond in one figure is the H-bond in every other.
    """
    BOND_TYPE_COLOUR: dict = field(default_factory=lambda: {
        "H-Bond":        "#4C72B0",
        "Salt Bridge":   "#DD8452",
        "Halogen":       "#55A868",
        "F-Polar":       "#C44E52",
        "F-Hydrophobic": "#8172B2",
        "Hydrophobic":   "#937860",
    })
    """
    Ordered series palettes. Indexed, never name-picked, so a figure cannot silently reorder its scale.
    RADAR carries one colour per plotted ligand; TREND one per metric on a multi-metric trend panel.
    """
    VIS_RADAR_SERIES: tuple = ("#057759", "#0BF1E2", "#E69F00", "#CC79A7", "#0072B2",
                               "#56B4E9", "#F0E442", "#009E73", "#D55E00", "#CC79A7")
    VIS_TREND_SERIES: tuple = ("#0072B2", "#E69F00", "#9467BD", "#009E73", "#D55E00")
    """
    Deep/secondary shades of the accents, for the marks that must sit ON a filled band of the same
    hue and still be legible (a dark-green label on the pale-green 'strong' fill, and so on).
    """
    VIS_ACCENT_DEEP: dict = field(default_factory=lambda: {
        "green":      "#005840",
        "green_alt":  "#1B4D2E",
        "teal":       "#0E7C7B",
        "blue":       "#00408B",
        "blue_alt":   "#1A4080",
        "blue_mid":   "#2c7fb8",
        "blue_light": "#2E86C1",
        "orange":     "#A03000",
        "orange_alt": "#803000",
        "orange_mid": "#B35400",
        "orange_hot": "#C04000",
        "gold":       "#D4A000",
        "gold_deep":  "#9A6B00",
        "gold_dark":  "#7A5C00",
        "gold_hot":   "#FFC300",
        "red":        "#A93226",
        "red_deep":   "#8B0000",
        "red_dark":   "#990000",
        "red_bright": "#D62728",
        "purple":     "#7B1FA2",
        "purple_deep":"#7A0177",
        "purple_dark":"#2E1B5E",
        "emerald":    "#2ECC71",
        "emerald_alt":"#27AE60",
        "leaf":       "#4CAF50",
        "moss":       "#2C4A1E",
        "pine":       "#1E4A3A",
        "jade":       "#17A589",
        "sea":        "#1B9E77",
        "cyan":       "#1B9E9E",
        "clover":     "#2ca02c",
        "sun":        "#F1C40F",
        "tangerine":  "#E67E22",
        "brick":      "#E74C3C",
        "rose_fill":  "#F8D7DA",
        "rose_edge":  "#FFF5F5",
        "mint_fill":  "#D4EDDA",
        "cream_fill": "#FFF3CD",
        "peach_fill": "#FBEEE6",
        "orange_deepest": "#6E2C00",   # the no-fit / potential-inhibitor mark
        "gold_muted":     "#C9A227",   # edge of an amber annotation box
    })
    VIS_PFAS_FCOUNT_BINS: list = field(default_factory=lambda: [0, 8, 13, 18, 24, float("inf")])  # total-fluorine-count bin edges for the PFAS chain-length size figure (Step 03)
    # Bin labels paired with VIS_PFAS_FCOUNT_BINS (must have len(bins)-1 entries; kept beside the
    # edges so they never drift apart). Multi-line for legend, short for axis ticks.
    VIS_PFAS_SIZE_LABELS: list = field(default_factory=lambda: [
        "Very Short\n(≤C4, ≤8F)", "Short-chain\n(C5–6, 9–13F)", "PFOA/PFOS\n(C7–9, 14–18F)",
        "Long-chain\n(C10–12, 19–24F)", "Ultra-long\n(≥C13, ≥25F)"])
    VIS_PFAS_SIZE_LABELS_SHORT: list = field(default_factory=lambda: [
        "≤C4 (≤8F)", "C5–6 (9–13F)", "C7–9 (14–18F)", "C10–12 (19–24F)", "≥C13 (≥25F)"])
    # PyMOL render palette (Step 05) — centralised so the presentation theme is tunable.
    VIS_PYMOL_BG_COLOR: str            = "gray60"
    VIS_PYMOL_SURFACE_TRANSPARENCY: float = 0.50
    VIS_PYMOL_CARTOON_TRANSPARENCY: float = 0.35
    VIS_PYMOL_LABEL_COLOR: str         = "white"
    VIS_PYMOL_SURFACE_RADIUS: float    = 18.0  # Å  render surface only within this radius of the ligand; the far protein is off-frame after zoom/clip, so this cuts ray time with no visible change

    # -------------------------------------------------------------------------------
    # Step 13.1b: Boltz confidence quality bands (figure shading — Step 03)
    # -------------------------------------------------------------------------------
    """
    Confidence-score band edges and Okabe–Ito colour-blind-safe band colours used by
    Step 03 figures to shade high / acceptable / below-threshold confidence zones.
    Single source so band edges and colours stay consistent across every figure.
    """
    CONF_BAND_HIGH: float       = 0.90   # confidence ≥ this → high-quality zone
    CONF_BAND_ACCEPTABLE: float = 0.80   # confidence ≥ this → acceptable zone
    CONF_BAND_COLOURS: dict = field(default_factory=lambda: {
        "high":       "#009E73",   # green  — high confidence
        "acceptable": "#E69F00",   # amber  — acceptable
        "below":      "#D55E00",   # vermillion — below threshold
    })
    """
    Active-site RMSD-to-crystal quality bands (Å) for Step 03 figure shading.
    ≤EXCELLENT green · ≤ACCEPTABLE amber · ≤DIVERGED vermillion · above → severe.
    """
    """
    Geometry helpers used while READING structures, not while judging them. The C–O cut-off simply
    asks 'are these two atoms bonded' (a C–O bond is 1.21-1.43 Å; 1.6 Å separates bonded from
    non-bonded with room to spare), and the degeneracy tolerance decides when two principal axes of
    a pocket are too close in length to be told apart.
    """
    BOND_CO_MAX_A: float          = 1.6    # Å  above this, a C and an O are not bonded
    PCA_DEGENERACY_TOL: float     = 0.15   # relative gap below which two axes are degenerate
    RMSD_BAND_EXCELLENT: float  = 1.0   # Å  ≤ this → excellent
    RMSD_BAND_ACCEPTABLE: float = 2.0   # Å  ≤ this → acceptable
    RMSD_BAND_DIVERGED: float   = 3.0   # Å  ≤ this → diverged (above → severe)
    """
    Mechanistic-score quality bands (0–1) — Step 03 figure shading. Rescaled for the
    holistic mech_score (full machinery alone = 0.70; the linear SN2 angle term lifts
    it toward 1.0): strong = complete anchors + a near-linear productive trajectory.
    """
    MECH_FP_BAND_STRONG: float   = 0.90   # ≥ this → strong
    MECH_FP_BAND_MODERATE: float = 0.75   # ≥ this → moderate (below → weak)
    """
    Active-site engagement fraction quality bands (0–1) — Step 03 figure shading.
    """
    ENGAGEMENT_BAND_HIGH: float     = 0.75   # ≥ this → high engagement
    ENGAGEMENT_BAND_MODERATE: float = 0.50   # ≥ this → moderate (below → low)

    # -------------------------------------------------------------------------------
    # Step 13.2: Per-structure rendering timeouts — seconds (Step 07)
    # -------------------------------------------------------------------------------
    VIS_TIMEOUT_PYMOL: int    = 600   # PyMOL render timeout (2× 2400-px ray traces, pocket-local surface — a few min on CPU)
    VIS_TIMEOUT_CHIMERAX: int = 180   # ChimeraX render timeout
    VIS_TIMEOUT_MAESTRO: int  = 300   # Maestro render timeout
    VIS_TIMEOUT_PLIP: int     =  90   # PLIP interaction analysis timeout
    VIS_TIMEOUT_LIGPLOT: int     =  60   # LigPlot+ render timeout
    VIS_TIMEOUT_PYMOL_HEAVY: int = 1200  # s  PyMOL timeout for structures with ≥4 C–F bonds (polyfluorinated PFAS)

    # -------------------------------------------------------------------------------
    # Step 13.3: Geometric highlight radii (Step 07)
    # -------------------------------------------------------------------------------
    VIS_F_CONTACT_RADIUS: float = 4.0   # Å  fluorine contact highlight sphere
    VIS_POCKET_RADIUS: float    = 5.5   # Å  binding-pocket cartoon / surface shell

    # -------------------------------------------------------------------------------
    # Step 13.4: Tier summary figure layout (Step 03)
    # -------------------------------------------------------------------------------
    """
    Pixel dimensions and typography for the per-tier stacked-bar / star-plot
    panels produced by 03_Validation_Figures_FAcDs.py.
    """
    VIS_TT_STAR_SIZE: int    = 460    # px  star marker diameter
    VIS_TT_SHRINK_BORDER: int =  12   # px  border shrink for tight layout
    VIS_TT_IMG_PX: int       = 800    # px  panel image width
    VIS_TT_TEXT_PX: int      = 220    # px  text annotation column width
    VIS_TT_BORDER_PX: int    =  16    # px  outer border thickness
    VIS_TT_PAD: int          =  28    # px  inter-panel padding
    VIS_TT_FONT_SIZE: int    =  44    # pt  annotation font size
    VIS_TT_AX_WIDTH: float   =   0.12 # fraction of figure width for axis panel
    VIS_TT_ZOOM_BUFFER: float =  1.0  # Å  PyMOL zoom padding around the pocket (lower = more zoomed-in)

    # -------------------------------------------------------------------------------
    # Step 13.5: Maximum display ranks (Step 03)
    # -------------------------------------------------------------------------------
    VIS_MAX_RANKS_DISPLAY: int = 50   # top-N entries shown in ranked output plots
    VIS_MAX_THUMBNAILS: int = 10   # max structure thumbnails embedded in a quality-space panel (top-N by Scientific_Rank)
    VIS_RADAR_MAX_HITS: int = 10   # max series plotted on a radar / spider chart (top-N by ranking key)
    # Slice colours for the model-selection pie inset (Figure 02); one per Boltz model index.
    VIS_PIE_MODEL_COLOURS: list = field(default_factory=lambda: [
        "#8ECFC9", "#A5C8E1", "#FABEBE", "#FFF4B8", "#C9E8C4",
        "#F8C8A8", "#D5C8E8", "#B8D4E8",
    ])
    '''
    Bar value-label segment colours (Figure 01): count, "|" separator and percentage
    are drawn in three distinct hues so each datum reads separately. Two variants keyed
    on the bar-fill luminance keep every segment legible against any role-group colour —
    the light triple sits on dark fills, the dark triple on light fills (see _text_color).
    '''
    VIS_BAR_LABEL_COLOURS_ON_DARK: list = field(default_factory=lambda: [
        "#FFFFFF",   # count       — white
        "#FFE066",   # separator   — light amber
        "#9AD0FF",   # percentage  — light azure
    ])
    VIS_BAR_LABEL_COLOURS_ON_LIGHT: list = field(default_factory=lambda: [
        "#1A1A1A",   # count       — near-black
        "#B05A00",   # separator   — burnt amber
        "#0050A0",   # percentage  — deep azure
    ])
    VIS_TOP_N_AUTO_SELECT_TIMEOUT: int = 30  # s  auto-select countdown in Step 05 interactive tier selection
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
        "Score":        ["Binding_Probability_Score", "Binding_Probability", "ActiveSite_Conservation_Score", "binding_likelihood_computed"]
    })

    # ===============================================================================
    # SECTION 14: SEQUENCE ALIGNMENT & SCORING PARAMETERS  (Step 02)
    # ===============================================================================

    # -------------------------------------------------------------------------------
    # Step 14.1: Pairwise sequence alignment
    # -------------------------------------------------------------------------------
    """
    Biopython PairwiseAligner configured for global BLOSUM62 alignment.
    Ref: Cock et al. (2009); Henikoff & Henikoff (1992) — see header §14.
    """
    ALIGN_OPEN_GAP_SCORE: float   = -10.0  # open-gap penalty (full gap initiation cost)
    ALIGN_EXTEND_GAP_SCORE: float = -0.5   # gap-extension penalty (per-residue gap cost)
    ALIGN_MIN_SEQ_IDENTITY: float = 25.0   # %  identity_pct ≥ this → active-site mapping flagged reliable
    # ActiveSite_Conservation_Score component weights (must sum to 1.0): machinery-led, not identity-led.
    CONSERV_W_INTEGRITY: float = 0.60  # Criterion-A active-site integrity (8 catalytic residues correctly mapped)
    CONSERV_W_GEO: float       = 0.30  # active-site geometric fit to control (RMSD-derived)
    CONSERV_W_IDENT: float     = 0.10  # global sequence identity (light corroborating signal)
    # Reference-structure conservation score (3R3U crystal + DeHa4 control): a 2-term
    # identity/geometry blend, distinct from the 3-term candidate score above. Control
    # rows pass identity = 100.0, so the identity term contributes CONSERV_REF_W_IDENT*100.
    CONSERV_REF_W_IDENT: float = 0.25  # reference-structure conservation: sequence-identity weight
    CONSERV_REF_W_GEO: float   = 0.75  # reference-structure conservation: active-site geometric-fit weight

    # -------------------------------------------------------------------------------
    # Step 14.1b: Alignment grade bins (Step 02)
    # -------------------------------------------------------------------------------
    """
    Used with pd.cut() to assign letter grades to Alignment_Score_Pct.
    bins[i] < score ≤ bins[i+1] → label[i].  right=True (default for pd.cut).
    """
    ALIGN_GRADE_BINS:   list = field(default_factory=lambda:
                              [0, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    ALIGN_GRADE_LABELS: list = field(default_factory=lambda:
                              ["I", "H", "G", "F", "E", "D", "C", "B", "A"])

    # -------------------------------------------------------------------------------
    # Step 14.2: Dynamic alignment cost-matrix scaling
    # -------------------------------------------------------------------------------
    ALIGN_DYNAMIC_PENALTY_FACTOR: float   = 2.0  # scale factor applied to max(cost) for dynamic penalty
    ALIGN_DYNAMIC_PENALTY_FALLBACK: float = 5.0  # penalty value when cost matrix is empty

    # -------------------------------------------------------------------------------
    # Step 14.2b: Alignment-free k-mer dendrogram (Step 04)
    # -------------------------------------------------------------------------------
    DENDRO_KMER_SIZE: int          = 3          # k-mer length for sequence frequency profiles
    DENDRO_DISTANCE_METRIC: str    = "cosine"   # pairwise distance metric (scipy pdist)
    DENDRO_LINKAGE_METHOD: str     = "average"  # hierarchical linkage (average = UPGMA)

    # -------------------------------------------------------------------------------
    # Step 14.3: Likelihood scoring thresholds (Step 02)
    # -------------------------------------------------------------------------------
    # Applied during Likelihood_Degrader_Score computation.
    SCORE_DESOLVATION_THRESHOLD: float        = 0.5   # hydrophobic desolvation ratio above which bonus applies
    SCORE_DESOLVATION_BONUS_FACTOR: float     = 10.0  # multiplier for desolvation bonus (× ratio)
    SCORE_INHIBITION_DENSITY_THRESHOLD: float = 1.5   # interaction-density above which the active-site contact flag fires
    SCORE_INHIBITION_PENALTY_FACTOR: float    = 15.0  # scaling factor for the active-site contact flag

    # -------------------------------------------------------------------------------
    # Step 14.4: GPU batch watchdog
    # -------------------------------------------------------------------------------
    GPU_WATCHDOG_TIMEOUT_PER_JOB: int = 600  # seconds per job before GPU-batch watchdog kills the run

    # -------------------------------------------------------------------------------
    # Step 14.5: Affinity interaction contribution weights (Step 02)
    # -------------------------------------------------------------------------------
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

    # ===============================================================================
    # SECTION 15: PDB PREPARATION — Schrödinger PrepWizard  (Step 05)
    # ===============================================================================

    # -------------------------------------------------------------------------------
    # Step 15.1: PrepWizard protonation & minimisation
    # -------------------------------------------------------------------------------
    # PrepWizard protocol — Sastry et al. (2013) — see header §15.
    PREPWIZARD_PROPKA_PH: float      = 8.0  # protein protonation pH (FAcD physiological context)
    PREPWIZARD_EPIK_PH: float        = 8.0  # PFAS ligand protonation pH via Epik
    """
    Catalytic protonation, enforced by ROLE after PrepWizard.

    PropKa assigns protonation per structure from a pKa prediction, and it does so without knowing
    which residue is the nucleophile. That is not a tuning parameter — it decides whether the
    chemistry can happen at all, and a pKa heuristic gets it wrong: across the three MD systems the
    SAME catalytic aspartate was given different charges, and every catalytic histidine came out as
    HIP (+1), a 'base' with no lone pair to accept a proton.

    The consequence is measurable, not theoretical. A controlled experiment (two systems identical
    but for ONE hydrogen, same build, same 5-stage Desmond relaxation) gave:
        HIP (+1)  : attack angle 150° → 98°   — the reactive pose is destroyed during relaxation
        HID (0)   : attack angle 150° → 170°  at 3.77 Å — a pose that passes relaxed NAC
    So the state each catalytic residue takes is fixed here by the MECHANISM it must perform, and
    PREPWIZARD_ENFORCE_PROTONATION applies it after PrepWizard, before anything is simulated.

    Each entry: the hydrogens to REMOVE so the residue takes its required state. Deleting the proton
    is what changes the state — Desmond's force-field templates read the charged (ASP) or neutral
    (ASH) form from the hydrogens present.
    """
    PREPWIZARD_ENFORCE_PROTONATION: bool = True
    CATALYTIC_PROTONATION_POLICY: dict = field(default_factory=lambda: {
        "Nuc":  {"state": "deprotonated", "charge": -1, "strip_H": ("HD2", "HE2"),
                 "why": "it attacks the warhead carbon; a neutral COOH cannot"},
        "Acid": {"state": "deprotonated", "charge": -1, "strip_H": ("HD2", "HE2"),
                 "why": "the dyad aspartate polarises the histidine; it accepts charge"},
        "Base": {"state": "neutral HID",  "charge": 0,  "strip_H": ("HE2",),
                 "why": "Nd1-H points at the dyad; the Ne2 lone pair takes the proton"},
    })
    PREPWIZARD_RMSD_RESTRAIN: float  = 0.3  # Å — RMSD restraint for clash-resolving minimisation
    """
    PrepWizard's restrained-minimisation force field. Its own default is OPLS_2005 — a different
    force field from the one every downstream stage uses (the Desmond system build, the MD production
    run, WaterMap and Prime MM-GBSA all run the OPLS4 family, exposed by prepwizard as 'S-OPLS').
    Left at the default, the structure is minimised on one potential and then simulated on another,
    so the prepared geometry the MD starts from is not a minimum of the MD's own potential and
    relaxes the moment the simulation begins. The only two values prepwizard accepts are S-OPLS and
    OPLS_2005.
    """
    PREPWIZARD_FORCEFIELD: str       = "S-OPLS"   # OPLS4-family; matches Desmond/WaterMap/Prime

    """
    QM (Jaguar ESP) partial charges for the MD-ready ligands — Step 05b.

    Desmond takes its ligand charges from OPLS4. A fixed-charge force field represents fluorine's low
    polarisability poorly, and the worst case is precisely this chemistry: a perfluoroalkyl carboxylate,
    a hard-charged head on a long, weakly-polarisable-modelled tail. The error lands on the alpha-carbon's
    electrophilicity — the one quantity the SN2 turns on — and nothing downstream compensates for it.

    icfit=1 fits the charges to the electrostatic potential. Verified by running it, not assumed:
    Jaguar prints 'Atomic charges from electrostatic potential' and the fitted charges sum to the
    formal charge (fluoroacetate: -1.000). The fit is done on the PREPARED geometry — the structure
    the MD actually starts from — because a charge set derived from a different conformer is a charge
    set for a different molecule.
    """
    ESP_CHARGE_BASIS: str            = "6-31G**"   # basis for the ESP single point
    ESP_CHARGE_DFT: str              = "b3lyp"     # functional for the ESP single point

    # -------------------------------------------------------------------------------
    # Step 15.2: Parallelism
    # -------------------------------------------------------------------------------
    PREP_CPU_RESERVE: int  = 2  # CPU cores to reserve for OS/desktop stability (not used by PrepWizard)

    # -------------------------------------------------------------------------------
    # Step 15.3: Chain assignment
    # -------------------------------------------------------------------------------
    PREP_LIGAND_CHAIN: str = "L"  # chain identifier for non-protein (ligand) residues in converted PDB

    # -------------------------------------------------------------------------------
    # Step 15.4: Standard amino-acid residue set (used for ligand filtering)
    # -------------------------------------------------------------------------------
    PREP_STANDARD_AA: set = field(default_factory=lambda: {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        # Protonation variants (force-field naming: OPLS4 / CHARMM / Amber)
        "HID", "HIE", "HIP", "HSE", "HSD", "HSP",
        "ASH", "GLH", "CYM", "CYX", "LYN", "TYM",
    })
    PREP_AMBIGUOUS_AA: set = field(default_factory=lambda: set("BXZJOU"))
    PREP_MAX_AMBIGUOUS_PCT: float = 5.0   # %  sequence rejected if ambiguous-residue fraction exceeds this

    # -------------------------------------------------------------------------------
    # Step 15.5: Protein-associated residues (stay in protein chain)
    # -------------------------------------------------------------------------------
    # Modified amino acids and structural metals that must NOT be moved to Chain L.
    PREP_PROTEIN_ASSOCIATED: set = field(default_factory=lambda: {
        # Modified AAs
        "MSE", "SEP", "TPO", "PTR", "KCX", "CME", "CSO", "PCA", "LLP",
        # Structural / catalytic metals
        "ZN", "MG", "CA", "FE", "MN", "CO", "NI", "CU", "NA", "K",
    })

    # ===============================================================================
    # SECTION 16: DATA REGISTRY & VISUAL AESTHETICS
    # ===============================================================================

    # -------------------------------------------------------------------------------
    # Step 16.1: CSV Column Names
    # -------------------------------------------------------------------------------
    # Registry for all master CSV columns to prevent hardcoding string keys.
    # Controls are emitted by 02 with a reserved zero job index, so a job name beginning with this
    # prefix IS a control. Kept here because several steps test for it and a literal "0000000" in
    # three scripts is three places to get it wrong.
    CONTROL_JOB_PREFIX: str = "0000000"
    COL_TIER:     str = "degrader_tier"
    COL_PROT:     str = "Protein_Name"
    COL_LIG:      str = "Ligand_Name"
    COL_CONF:     str = "Boltz_Model_Confidence"
    COL_PTM:      str = "ptm"
    COL_IPTM:     str = "iptm"
    COL_IDENS:    str = "interaction_density"
    COL_ID_PCT:   str = "identity_pct"
    COL_TRUST:    str = "Trust_Score"
    COL_SN2:      str = "SN2_Attack_Angle"
    COL_NUC_DIST: str = "Dist_Nucleophile"
    COL_MECH_S:   str = "mechanistic_score"
    COL_LIKE_S:   str = "ActiveSite_Conservation_Score"
    COL_ALN_G:    str = "Alignment_Grade"

    # -------------------------------------------------------------------------------
    # Step 16.2: Visual Plotting Properties
    # -------------------------------------------------------------------------------
    """
    Canonical marker sizes (S) and opacities (A) for each tier.
    Tier_1/Tier_2 tiers are drawn larger and more opaque.
    """
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

    # -------------------------------------------------------------------------------
    # Step 16.3: Alignment Grade Bins
    # -------------------------------------------------------------------------------
    # Define identity % ranges for each grade letter.
    ALIGN_GRADE_DEFS: list = field(default_factory=lambda: [
        (90.0, 100.0, "A"), (80.0, 90.0, "B"), (70.0, 80.0, "C"),
        (60.0, 70.0,  "D"), (50.0, 60.0, "E"), (40.0, 50.0, "F"),
        (30.0, 40.0,  "G"), (20.0, 30.0, "H"), (0.0,  20.0, "I")
    ])

    # -------------------------------------------------------------------------------
    # Step 16.4: Tier meanings
    # -------------------------------------------------------------------------------
    TIER_MEANING: dict = field(default_factory=lambda: {
        "Tier_1A": "Elite",
        "Tier_1B": "Excellent",
        "Tier_2A": "Strong",
        "Tier_2B": "Good",
        "Tier_3":  "Moderate",
        "Tier_4":  "Poor",
        "Tier_5_Decoy": "Invalid / Decoy",
        "Control": "Reference Control"
    })

    # -------------------------------------------------------------------------------
    # Step 16.5: Canonical pipeline artefact filenames (cross-script SSOT)
    # -------------------------------------------------------------------------------
    # Written by one stage, read by another — centralised so the producer and
    # consumer filename cannot drift apart.
    """
    File-name masks. Every one of these is a CONTRACT between two steps — 05 globs what 03 wrote,
    07 globs what 06 wrote — and a mask spelled out at the call site is a contract only one side
    can see. A rename then breaks the consumer silently (an empty glob reads as 'nothing to do',
    not as an error).
    """
    GLOB_RANKED_CSV:    str = "7_Boltz2_FAcDs_Ranked_*.csv"   # 02 writes → 05/07 read
    SUFFIX_RAW_PDB:     str = "_RAW.pdb"                      # 05 internal (pre-prep structure)
    SUFFIX_PLIP_LOG:    str = "_PLIP.log"                     # 05 internal (interaction run log)
    SUFFIX_SID_EAF:     str = "_SID-out.eaf"                  # Desmond SID → 06 reads
    SUFFIX_MMGBSA_CSV:  str = "-out-prime-mmgbsa.csv"         # Prime → 06 reads
    SUFFIX_CMS_OUT:     str = "-out.cms"                      # Desmond → 06/07 read
    SUFFIX_NAC_DATA:    str = "_NAC_Data.csv"                 # 07 writes → 07 figures read
    FILE_ALIGNMENT_STATS:  str = "Alignment_Stats.csv"           # 02 writes → 03 reads
    FILE_MMGBSA_SUMMARY:   str = "00_MMGBSA_Summary.csv"         # 06 writes → 07 reads
    FILE_VALIDATED_MASTER: str = "03_Final_Validated_Master.csv" # 03 writes → 04 reads

    # ===============================================================================
    # SECTION 17: PRIME MM-GBSA  (Step 06 — end-state binding free energy)
    # ===============================================================================
    """
    Runs $SCHRODINGER/run thermal_mmgbsa.py <job>-out.cms per completed MD job.
    MM-GBSA gives the ensemble ligand BINDING free energy (non-covalent,
    end-state) — complementary to, NOT a substitute for, the QSite QM/MM
    reaction barrier. GB implicit solvent overstabilises anionic PFAS, so treat
    ΔG_bind as a RELATIVE ranking, never an absolute value.
    """
    MMGBSA_RUN: bool        = True           # run thermal_mmgbsa.py in Step 06
    """
    Frame subsampling for Prime. 0 or 1 = every frame; N>1 appends `-step_size N`, scoring
    every Nth frame. Prime cost is linear in the number of structures, so this is the only
    knob that changes MM-GBSA wall-clock by an order of magnitude.

    Why a stride costs almost nothing statistically: a long MD run samples frames far more
    finely than the ~0.1–1 ns decorrelation time of a bound pose, so neighbouring frames are
    near-duplicates rather than independent samples. The precision of ⟨ΔG_bind⟩ is set by the
    number of INDEPENDENT samples (N_eff ≈ T_total / 2τ) — a function of the trajectory's
    LENGTH, not of how finely it is sliced. Any stride whose spacing stays below τ therefore
    discards redundancy, not information, and leaves ⟨ΔG_bind⟩, its spread and the ranking
    statistically indistinguishable. Step 06 derives and prints the resulting frame spacing
    from the trajectory itself, so this value can be changed freely without touching code.
    """
    MMGBSA_STEP_SIZE: int   = 10
    MMGBSA_LIGAND_ASL: str  = "res.ptype LIG"   # ASL passed to thermal_mmgbsa via its -lig_asl flag so Prime scores the correct molecule (matches Step 07's --lig LIG convention). A heavily fluorinated PFAS can be misassigned as solvent by auto-detection; empty string "" reverts to auto-detect.
    MMGBSA_DG_COLUMN: str   = "r_psp_MMGBSA_dG_Bind"   # primary per-frame dG_bind column in the thermal_mmgbsa CSV
    MMGBSA_TIMEOUT_SEC: int = 0               # 0 = no timeout (Prime can run for hours); >0 caps each job
    MMGBSA_OUTPUT_SUBDIR: str = "Prime_MMGBSA"  # figures folder under <run>/6_Physics_Validation/MolecularDynamics/ (path derived, not hardcoded)
    SCHRODINGER_SCRATCH_SUBDIR: str = "_Schrodinger_Scratch"  # job scratch dir under the run's MolecularDynamics working folder (large disk); keeps Prime's hundreds-of-GB per-subjob staging off /tmp on the OS disk
    SCHRODINGER_SCRATCH_COOLDOWN_SEC: int = 60  # settle window after the last job finishes: the job server is still flushing/copying outputs back to the working folders, so wait before deleting the scratch dir
    """
    End-of-run housekeeping for the two Schrödinger directories.

    _Schrodinger_Scratch (SCHRODINGER_SCRATCH_SUBDIR) is the CLIENT's TMPDIR — pure transient
    working space, deleted after the settle window above.

    _Schrodinger_JobServer (SCHRODINGER_JOBSERVER_SUBDIR) is NOT scratch: it is the job-server
    daemon's home (jobdb.sqlite, filestore/, logs/, bin/). Deleting it under a live daemon is
    what produces 'Error locating localhost job server config' and kills every submission. What
    actually GROWS inside it is the per-job filestore and logs, so the safe cleanup is to delete
    COMPLETED JOBS from the server (`jsc delete`), leaving the daemon and its database intact.
    Removing the whole directory is only ever safe with the server stopped and idle; it is
    off by default (the next run recreates it, but it costs a server restart).
    """
    SCHRODINGER_JOBSERVER_PURGE_COMPLETED: bool = True   # `jsc delete` finished jobs → frees filestore/logs, daemon stays up
    SCHRODINGER_JOBSERVER_REMOVE_WHEN_IDLE: bool = False # stop the idle server and delete its whole home dir (only when nothing is queued)
    SCHRODINGER_JOBSERVER_SUBDIR: str = "_Schrodinger_JobServer"  # working-disk directory (at the project root) that becomes the Schrödinger job server's scratch `tmpdir`. A Prime MM-GBSA / QSite run stages hundreds of GB of per-subjob scratch under the server's tmpdir (the localhost entry of $SCHRODINGER/schrodinger.hosts); its /tmp default (OS disk) is what overflows and kills the job (rc=1, copy_file_range: no space left on device). ensure_jobserver_on_working_disk (00_03) points that tmpdir here and applies it LIVE via `jsc admin reload-hosts` — no server stop, no killed job. Env vars (SCHRODINGER_TMPDIR/TMPDIR) alone do NOT move an already-running server's scratch.
    MMGBSA_PROGRESS_INTERVAL_SEC: int = 30    # heartbeat cadence for the in-place (\r) MM-GBSA progress ticker
    """
    Failed Prime minimisations. A small PFAS ligand cannot bind at −1000 kcal/mol; frames that
    far outside the ensemble are minimisation artefacts (a blown-up structure), not physics. They
    are a fraction of a percent but they wreck an axis and drag the arithmetic mean, so the per-job
    figure scales to the ROBUST core (percentile clip) and flags them explicitly rather than
    silently deleting them. Frames beyond K × IQR from the quartiles are counted as failures.
    """
    MMGBSA_DG_OUTLIER_IQR_K: float = 3.0   # Tukey fence multiplier used to FLAG failed minimisations
    MMGBSA_PLOT_CLIP_PCT: float = 0.5      # per-job figure y-axis spans this to (100 − this) percentile
    """
    Statistics for the cross-rank comparison. MD frames are NOT independent: neighbouring frames
    are the same configuration re-measured, so a t-test or Mann-Whitney over ~10⁴ frames returns
    p ≈ 0 for any difference whatsoever and is meaningless. The defensible test is a MOVING-BLOCK
    BOOTSTRAP — resample contiguous blocks longer than the correlation time, so each block is
    effectively one independent draw — giving a confidence interval on the median. Effect size is
    reported as Cliff's delta, which needs no distributional assumption and is unmoved by outliers.
    """
    MMGBSA_BOOTSTRAP_N: int = 2000          # bootstrap resamples for the median confidence interval
    MMGBSA_BOOTSTRAP_BLOCK: int = 50        # FLOOR on frames per block; the actual block is derived per series
    MMGBSA_BOOTSTRAP_BLOCK_TAU_MULT: float = 2.0   # block = this × the series' measured autocorrelation time τ
                                                   # (a fixed block shorter than τ leaves the blocks correlated and
                                                   #  yields a falsely narrow CI — the interval must reflect N_eff)
    MMGBSA_BOOTSTRAP_CI: float = 95.0       # confidence level (%)
    MMGBSA_LEGEND_MAX_ROWS: int = 3         # legend column count is DERIVED from this: ncol = ceil(entries / rows),
                                            # so the columns stay balanced as entries are added or candidates change
    MMGBSA_TIME_WINDOW_NS: float = 100.0    # width of the time windows the ΔG distribution is resolved into
                                            # (panel B): shows drift, stability and any change of binding mode
                                            # across the trajectory, which a single pooled violin hides
    MMGBSA_RANK_PALETTE: tuple = ("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                                  "#E69F00", "#56B4E9")   # colour-blind-safe, one per compared rank
    # Grid colours for the twin-axis panel: the two axes carry different quantities, so their
    # gridlines are coloured to match the axis they belong to and cannot be misread as one another.
    MMGBSA_GRID_LEFT: str  = "#4C72B0"   # bottom/left axis  (time course)
    MMGBSA_GRID_RIGHT: str = "#937860"   # top/right axis    (cumulative fraction)
    """
    Ink palette — every non-candidate mark colour in Step 06's figures. The candidate colours come
    from MMGBSA_RANK_PALETTE; these are the neutrals and accents used for medians, means, outlier
    flags, strokes and legend proxies. Kept here so a restyle is one edit, never a hunt through the
    plotting code.
    """
    MMGBSA_INK: dict = field(default_factory=lambda: {
        "dark":    "#111111",   # medians, in-figure text, box edges
        "light":   "#FFFFFF",   # box fills and the stroke behind outlined text
        "outline": "#222222",   # bar edges, zero line
        "muted":   "#555555",   # legend proxy marks, whiskers in proxies
        "soft":    "#333333",   # secondary annotation text
        "faint":   "#BBBBBB",   # band proxies, note borders
        "ghost":   "#999999",   # the mean line, where it is a secondary reading
        "mean":    "#D55E00",   # mean marker/line (it is an accent, not a candidate colour)
        "warn":    "#CC3311",   # failed minimisations, outlier counts
        "series":  "#0072B2",   # the single-series per-rank profile (no candidate to colour by)
        "hist":    "#1B9E77",   # the per-rank profile's histogram
    })
    """
    ── Step 07: the reactive-pose figures (MM-GBSA decomposition · machinery engagement) ──
    Both figures compare the whole trajectory against the frames in which the pose is reactive, so
    every constant that decides what is DRAWN lives here rather than in the plotting code.

    The occupancy criterion is not a new number: a residue 'holds contact' when it is within
    THRESHOLD_SALT_BRIDGE of the warhead carbon — the outermost of the four criterion bands the
    figure already draws — so the percentage written in each bar is judged by the same cut-off the
    band shows.
    """
    DEFLUOR_COMPONENT_MIN_KCAL: float = 0.5     # an MM-GBSA term is plotted only if it reaches this
                                                # magnitude in at least one candidate; below it the
                                                # term is numerically dead and the box is empty
    DEFLUOR_COMPONENT_TICK_KCAL: float = 5.0    # ΔG axis tick step: the components span 60 down to
                                                # 1 kcal/mol, so the axis must resolve the small ones
    DEFLUOR_ENGAGE_Y_MIN_TOP: float = 10.0      # floor on the engagement axis ceiling (the ceiling
                                                # itself is taken from the data), so the criterion
                                                # bands always have room
    """
    The four criterion bands, in the order the engagement figure stacks them: reactive contact,
    H-bond range, relaxed NAC, electrostatic range. The ramp runs green → amber → red, so a bar's
    height alone says whether the residue sits where the mechanism works or is merely in
    electrostatic reach.
    """
    DEFLUOR_BAND_COLOURS: tuple = ("#2E9E5B", "#7FBF3F", "#E8A33D", "#D1495B")
    DEFLUOR_BAND_ALPHA: float = 0.13            # band fill: under the bars, still four distinct bands
    """
    Ligand abbreviations for every figure in the pipeline. The full names ('Difluoroacetate') are
    long enough to push a legend across a panel, and the abbreviations are the names actually used
    in the chemistry. A ligand not listed here keeps its full name rather than being mangled.
    """
    VIS_LIGAND_SHORT: dict = field(default_factory=lambda: {
        "fluoroacetate":    "FA",
        "difluoroacetate":  "DFA",
        "trifluoroacetate": "TFA",
        "tfa":              "TFA",
    })
    """
    Prime subjob parallelism. The target is GLOBAL_MAX_WORKERS (cpu_count − PREP_CPU_RESERVE),
    exactly like every other step; MMGBSA_MAX_NJOBS = 0 means 'no extra ceiling'. Step 06 only
    drops below that target when the disk holding the Schrödinger scratch cannot hold one
    ~MMGBSA_SCRATCH_GB_PER_SUBJOB staging copy per subjob within MMGBSA_SCRATCH_HEADROOM_FRAC
    of its free space — the condition that fills the disk and kills the job.
    """
    """
    Tokens the SID-out.eaf Result vector carries per trajectory frame. Step 06's completion gate is
    EAF_TOKENS_PER_FRAME · traj_frames: 1 is the standard per-frame scalar series. If a SID analysis
    writes k > 1 tokens per frame, a partial EAF can still reach traj_frames tokens and be misjudged
    complete under a bare per-frame test; scaling the threshold by the true k closes that gap.
    """
    EAF_TOKENS_PER_FRAME: int = 1
    MMGBSA_MAX_NJOBS: int = 0                      # 0 = no ceiling beyond GLOBAL_MAX_WORKERS
    MMGBSA_SCRATCH_GB_PER_SUBJOB: float = 22.0     # measured staging footprint of one Prime subjob (full complexes copy)
    MMGBSA_SCRATCH_HEADROOM_FRAC: float = 0.75     # fraction of the scratch disk's free space the run may occupy
    """
    Sharded MM-GBSA. A bare thermal_mmgbsa.py run is serial where it hurts: ONE process
    reads every trajectory frame (one core, hours, and its RSS grows ~0.36 GB per 1000
    frames — ~40 GB by frame 100k) and only then hands the whole ensemble to Prime, so the
    other cores idle throughout. Step 06 instead splits the trajectory into contiguous frame
    shards (-start_frame/-end_frame) and runs MMGBSA_SHARD_CONCURRENCY of them at once, each
    with its own Prime subjobs. Reading is then parallel, one shard minimises while the next
    still reads, no core is left idle, and per-process memory is bounded by the shard size.
    Shards are disjoint and cover the full trajectory, so every frame is still scored — this
    is a scheduling change, not a subsampling one. Shard CSVs are concatenated into the
    single per-job CSV the rest of the step already reads.
    Set MMGBSA_SHARD_FRAMES = 0 to fall back to one plain serial thermal_mmgbsa run.
    """
    MMGBSA_SHARD_FRAMES: int = 2000        # trajectory frames per shard (0 = no sharding)
    MMGBSA_SHARD_PRIME_NJOBS: int = 6      # Prime subjobs inside one shard
    MMGBSA_SHARD_CONCURRENCY: int = 0      # shards in flight (0 = auto: GLOBAL_MAX_WORKERS // MMGBSA_SHARD_PRIME_NJOBS, then clamped to free RAM)
    """
    Memory budget for the shard plan: a shard costs one reader plus its MMGBSA_SHARD_PRIME_NJOBS
    Prime subjobs, and Prime is the dominant term (measured ~1.8 GB × 30 subjobs ≈ 54 GB — enough
    to exhaust a 60 GB box on its own and push it into swap). Step 06 shrinks the shard
    concurrency until reader + Prime fit in free RAM. Raise MMGBSA_RAM_HEADROOM_FRAC (or lower
    MMGBSA_PRIME_RAM_GB) to run more subjobs and accept some swapping; lower it to be safer.
    A reader's RSS is a base plus growth with the frames it reads, so a stride shrinks it.
    """
    MMGBSA_READER_RAM_BASE_GB: float = 3.6          # reader RSS before any frames are read
    MMGBSA_READER_RAM_PER_1K_FRAMES_GB: float = 0.36  # reader RSS growth per 1000 frames read
    MMGBSA_PRIME_RAM_GB: float = 1.8                # measured peak RSS of ONE Prime subjob
    MMGBSA_RAM_HEADROOM_FRAC: float = 0.85          # fraction of free RAM the run may occupy
    MMGBSA_SHARD_SUBDIR: str = "_MMGBSA_Shards"   # per-shard logs and Prime outputs live here, out of the job folder's glob path
    """
    Frame-ensemble averaging estimator for the headline per-job ΔG_bind.
    "mean" = arithmetic ⟨ΔGᵢ⟩ (default; the standard thermal MM-GBSA estimate — MD
    frames are already Boltzmann-sampled, so an unweighted mean is the ensemble
    average); "median" = robust to per-frame outliers/failed Prime frames;
    "boltzmann" = −RT ln⟨exp(−ΔGᵢ/RT)⟩ — reserve for non-Boltzmann-sampled ensembles
    only, as it re-weights an already-canonical ensemble and collapses toward the
    single most negative frame. The summary CSV always reports all three columns;
    this knob only selects which one drives the combined bar plot and the headline value.
    """
    MMGBSA_AVERAGING: str   = "mean"          # "mean" | "median" | "boltzmann"
    GAS_CONSTANT_KCAL: float = 1.9872036e-3   # kcal/mol/K — R for the boltzmann estimator
    MMGBSA_TEMPERATURE_K: float = 298.15      # K — ensemble temperature for the boltzmann estimator

    # ===============================================================================
    # SECTION 18: MD-READY SELECTION  (gates heavy downstream compute — Steps 05→07)
    # ===============================================================================
    """
    Single source of truth for WHICH complexes receive the expensive downstream
    pipeline (CIF→PDB, PrepWizard, MM-GBSA, MD). Step 02 writes two columns into the
    ranked CSV — MD_Selected (bool) and MD_Rank (1..N over the selected set, ordered
    by Scientific_Rank) — and Steps 05/06/07/08 prepare/simulate ONLY those rows.
    Analysis and figures still run over the full population (all 58k); only heavy
    compute is gated, so the tier/decoy statistics are unaffected.

    Modes:
      • "tier"       — every complex whose degrader_tier is in MD_TIERS.
      • "topN"       — the top MD_TOP_N complexes by Scientific_Rank.
      • "per_ligand" — the single best (lowest Scientific_Rank) complex per ligand,
                       restricted to tier ≥ MD_PER_LIGAND_TIER. When MD_PER_LIGAND_AUTO
                       is True the ligand roster is DATA-DRIVEN: every unique ligand that
                       reached the tier is represented (no hardcoded selection). When
                       False the explicit MD_PER_LIGAND panel is used instead.
    """
    MD_SELECTION_MODE: str  = "per_ligand"        # "tier" | "topN" | "per_ligand"
    MD_TIERS: list          = field(default_factory=lambda: ["Tier_1A"])
    MD_TOP_N: int           = 10                   # used when MD_SELECTION_MODE == "topN"
    MD_PER_LIGAND_TIER: list = field(default_factory=lambda: ["Tier_1A"])   # per_ligand MD reps: one best complex per unique ligand, drawn from Tier_1A only (accepts a single tier string too). Keeps the MD cohort strictly elite — a small, honest set of the highest-confidence degraders — rather than diluting it with lower-tier leads
    # Data-driven roster (default): one best complex per UNIQUE ligand that reached
    # MD_PER_LIGAND_TIER — no hardcoded ligand list. Set False to use the explicit
    # MD_PER_LIGAND panel below (curated chemotype-stratified subset).
    MD_PER_LIGAND_AUTO: bool = True
    MD_PER_LIGAND: list     = field(default_factory=lambda: [
        "Fluoroacetate", "Difluoroacetate", "TFA",        # α-fluorination ladder (mono/di/tri)
        "PFBA", "PFHxA", "PFOA", "PFDA", "PFTrDA",         # PFCA chain-length gradient (C4→C13)
        "GenX",                                            # branched perfluoroether
    ])
    MD_SELECTED_COL: str = "MD_Selected"          # bool column written to the ranked CSV
    MD_RANK_COL: str     = "MD_Rank"              # int (1..N over selected), NaN otherwise

    # ===============================================================================
    # SECTION 19: SSOT FILE-NAMING EXTENSIONS
    # ===============================================================================
    """
    Standard suffixes appended to file basenames for downstream tasks.
    Centralized here to ensure uniformity across the Preparation and MD stages.
    """
    EXT_QMMM_READY: str  = "_QMMM_Ready.pdb"
    EXT_LIGAND_SDF: str  = "_Ligand.sdf"
    EXT_COMPLEX_PDB: str = "_Complex.pdb"

    def __post_init__(self):
        """
        Internal-consistency guards (read-only; frozen-dataclass safe). Several constants
        are intentionally equal to another CFG value or must sum to 1.0. A field default
        edit that breaks one of these couplings fails loudly at import rather than silently
        drifting apart. Assertions only — no mutation.
        """
        import math as _math
        _isclose = lambda a, b: _math.isclose(float(a), float(b), abs_tol=1e-9)
        # Re-typed literals that must track their documented source value.
        assert _isclose(self.SOFT_NB_MIDPOINT, self.THRESHOLD_TRIAD_NB), "SOFT_NB_MIDPOINT must equal THRESHOLD_TRIAD_NB"
        assert _isclose(self.TIER_NUC_DIST["Tier_2A"], self.NAC_DIST_STRICT), "TIER_NUC_DIST['Tier_2A'] must equal NAC_DIST_STRICT"
        assert _isclose(self.TIER_NUC_DIST["Tier_2B"], self.NAC_DIST_RELAXED), "TIER_NUC_DIST['Tier_2B'] must equal NAC_DIST_RELAXED"
        assert _isclose(self.SUBSTRATE_ANGLE_MIN, self.TIER_ANGLE_MIN["Tier_1B"]), "SUBSTRATE_ANGLE_MIN must equal TIER_ANGLE_MIN['Tier_1B']"
        assert _isclose(self.INHIBITOR_ANGLE_MAX, self.NAC_ANGLE_RELAXED), "INHIBITOR_ANGLE_MAX must equal NAC_ANGLE_RELAXED"
        # Component-weight sets that must sum to 1.0.
        assert _isclose(self.MECH_W_NUC + self.MECH_W_NB + self.MECH_W_BA + self.MECH_W_CLAMP + self.MECH_W_STAB + self.MECH_W_ANGLE, 1.0), "mechanistic_score weights must sum to 1.0"
        assert _isclose(self.COMP_W_ANGLE + self.COMP_W_DIST + self.COMP_W_CLAMP + self.COMP_W_TRAJ + self.COMP_W_TRIAD + self.COMP_W_HALIDE, 1.0), "competence_score weights must sum to 1.0"
        assert _isclose(self.SOFT_W_NUC + self.SOFT_W_ANG + self.SOFT_W_INT, 1.0), "soft_catalytic_score weights must sum to 1.0"
        assert _isclose(self.CONSERV_W_INTEGRITY + self.CONSERV_W_GEO + self.CONSERV_W_IDENT, 1.0), "candidate conservation weights must sum to 1.0"
        assert _isclose(self.CONSERV_REF_W_IDENT + self.CONSERV_REF_W_GEO, 1.0), "reference conservation weights must sum to 1.0"
