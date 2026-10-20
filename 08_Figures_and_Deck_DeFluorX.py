#!/usr/bin/env python3
"""
===============================================================================
DeFluorX Pipeline  |  Step 08  |  Figure & Deck Assembly
===============================================================================

Final, non-computational step: assembles the publication-grade panels that the
earlier steps (02-07) wrote across the run into a small set of MERGED, Nature-style
manuscript figures + an editable PowerPoint deck - and builds the one figure that
is composed here from scratch: the cross-enzyme DEFLUORINATION CASCADE scheme.

It does three things:
  1. Cascade scheme - selects the best per-class degrader (TFA, DFA, MFA) from this
     run's 01_Master_Ranking.csv and draws the hydrolytic double-displacement for
     each (Asp110 nucleophile -> glycolyl-enzyme ester -> His280/water hydrolysis),
     directly from the real QM .mae geometry, with a per-band QM/MM reaction-
     coordinate plot (scissile C-F length + departing-F charge) and a joined
     QM/MM . geometry . WaterMap . MD data strip. Built as panel a of the
     mechanism / QM-MM-validation figure.
  2. Merged composites - reads the step panels IN PLACE (never copies the sources)
     and tiles related panels into portrait, tight, panel-lettered figures.
  3. Editable deck - rebuilds the .pptx from the merged SVGs (python-pptx + rsvg).

Panels missing (a step still running) are skipped and fold in on a re-run; a
coverage audit at the end lists any source panel no figure placed. Output format
follows the FIG_FORMAT decider (default "svg", mirroring CFG.VIS_FIGURE_FORMAT).

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 20 October 2026 <───────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 08_Figures_and_Deck_DeFluorX.py
  Role          : Manuscript figure + deck assembly; final non-computational step.
  Imports from  : (self-contained; standard scientific stack only - Pillow, lxml,
                  numpy, matplotlib, python-pptx. No pipeline module is imported.)
  Reads         : <Run>/3_Validation_Figures/ · 4_Dendrogram/ · 5_TopN_and_Preparation/
                  <Run>/6_Physics_Validation/ · 7_QMMM_Defluorination/ (step panels, in place)
                  <Run>/7_QMMM_Defluorination/01_Master_Ranking.csv          (cascade selection)
                  <Run>/7_QMMM_Defluorination/Rank_N/10_QSite_Scan_Data.csv   (scan plot)
                  <Run>/7_QMMM_Defluorination/Rank_N/01_MD_Stats.json         (WaterMap/MD strip)
                  <Run>/7_QMMM_Defluorination/Rank_N/QSite_SN2/*/*.01.mae     (substrate geometry)
                  01_Journal_Files/Figures_HQ/ (CURATED: Figure-1 workflow HTML, _assets)
  Writes        : <Run>/8_Figures_and_Deck_HQ/
                    - _cascade/Reaction_Cascade.{svg,png}          (cascade source)
                    - Figure_NN_<Name>.svg                         (merged main figures; the
                                                                   cascade is panel a of the
                                                                   mechanism / QM-MM figure)
                    - Supplementary_Figure_NN_<Name>.svg
                    - DeFluorX_Figures.pptx                         (editable deck)
  Env           : base python3 (python-pptx + rsvg-convert for the deck); the figure
                  build itself needs only Pillow + lxml + numpy + matplotlib.

All thresholds, colours and the output format decider mirror CFG (SSOT).
"""
from __future__ import annotations
from pathlib import Path
from io import BytesIO
import os, subprocess, tempfile, shutil, time, base64, re

from PIL import Image, ImageDraw, ImageFont, ImageChops
Image.MAX_IMAGE_PIXELS = None
from lxml import etree   # true-vector SVG compositing (nest panel SVGs, namespace their ids)

import sys as _sys, re as _cre

'''
Console styling - a small, self-contained colour layer so this step's log matches the rest of the
pipeline (coloured markers, headers, banners and a summary table) without importing a pipeline module.
_say() colourises a leading status marker (ok / !! / .. / -- / ::) and otherwise behaves like print;
it writes through sys.stdout directly so a global print -> _say rename cannot recurse. Colour is
disabled when NO_COLOR is set.
'''
class _C:
    OKGREEN = "\033[92m"; WARNING = "\033[93m"; FAIL = "\033[91m"; OKBLUE = "\033[94m"
    MAGENTA = "\033[95m"; CYAN = "\033[96m"; BOLD = "\033[1m"; DIM = "\033[2m"; ENDC = "\033[0m"
_USE_COLOUR = os.environ.get("NO_COLOR") is None
_MARK = {"ok": (_C.OKGREEN, "✔"), "!!": (_C.FAIL, "✖"), "..": (_C.OKBLUE, "·"),
         "--": (_C.DIM, "·"), "::": (_C.CYAN, "•")}
def _col(text: str, colour: str) -> str:
    return f"{colour}{text}{_C.ENDC}" if _USE_COLOUR else text
def _colourise(s: str) -> str:
    m = _cre.match(r'^(\s*)(ok|!!|\.\.|--|::)(\s.*)$', s, _cre.S)
    if m and _USE_COLOUR:
        ind, mk, rest = m.groups(); colour, sym = _MARK[mk]
        return f"{ind}{colour}{sym}{_C.ENDC}{rest}"
    return s
def _say(*args, **kw):
    end = kw.get("end", "\n"); sep = kw.get("sep", " ")
    out = _colourise(args[0]) if (len(args) == 1 and isinstance(args[0], str)) else sep.join(str(a) for a in args)
    _sys.stdout.write(out + end); _sys.stdout.flush()
_SEP = "═" * 80
_ARROW = "▶"
_PIPE = "│"
def _head(text: str) -> None:
    _say(_col(text, _C.CYAN + _C.BOLD))
def _banner(title: str, subtitle: str = "") -> None:
    _title = _col(_ARROW + "  " + title, _C.MAGENTA + _C.BOLD)
    _now = time.strftime("%Y-%m-%d %H:%M:%S")
    _say("\n" + _SEP)
    _say("  " + _title + "  " + _PIPE + "  DeFluorX Pipeline")
    if subtitle:
        _say("  " + subtitle)
    _say("  Started : " + _now)
    _say(_SEP + "\n")

# ----------------------------------------------------------------------------------------------------
# Output format - single decider (mirrors the pipeline SSOT CFG.VIS_FIGURE_FORMAT). "svg" => true-vector
# composites: each vector panel is nested as a live <svg> (crisp at any zoom); panels whose SOURCE is
# irreducibly raster (PyMOL renders, Desmond SID, dendrogram matrix-exports) are embedded as <image>.
# Set to "png" for the Pillow raster-montage path instead of the vector SVG montage.
# ----------------------------------------------------------------------------------------------------
FIG_FORMAT = "svg"
SVG_NS   = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
_SVG_IDREF_ATTRS = ("clip-path", "mask", "fill", "stroke", "filter",
                    "marker-start", "marker-mid", "marker-end")

SRC_ROOT = Path("/mnt/wdpassport/02_DeFluorX/02_Data")
# Curated manuscript inputs (the hand-authored Figure-1 workflow HTML, cropped substitute panels and
# the _assets folder) live beside the journal manuscript, OUTSIDE the pipeline run; read-only here.
CURATED = Path("/mnt/wdpassport/02_DeFluorX/01_Journal_Files/Figures_HQ")


def _latest_run() -> Path:
    # The committed controls-only demo (an EXAMPLE run, or the _DeFluorX-suffixed tracked copy that
    # ships in the repo) is never a figure source - only the full local run supplies data panels.
    runs = sorted((p for p in SRC_ROOT.glob("Boltz-2_Run_*")
                   if "EXAMPLE" not in p.name and not p.name.endswith("_DeFluorX")),
                  key=lambda p: p.stat().st_mtime)
    return runs[-1] if runs else SRC_ROOT / "Boltz-2_Run_20260309T085406Z"


def _resolve_run() -> Path:
    """Run to assemble: an explicit run id or path on the command line (as the other steps take
    `<Run>`), else the newest full local run. The deck/figures mode tokens are not run ids."""
    import sys
    for a in sys.argv[1:]:
        if a in ("pptx", "figures") or a.startswith("-"):
            continue
        p = Path(a)
        if p.is_dir() and p.name.startswith("Boltz-2_Run"):
            return p
        cand = SRC_ROOT / a
        if cand.is_dir():
            return cand
    return _latest_run()


RUN = _resolve_run()
# Step-08 output: the HQ merged journal figures + deck land inside the run as a numbered step folder.
OUT_DIR = RUN / "8_Figures_and_Deck_HQ"
# Figure 01 is the mechanism-guided workflow / graphical overview, authored as a self-consistent HTML
# (curated input) and rendered here to PNG so it sits in the same numbered series as the composites.
WORKFLOW_HTML = CURATED / "Figure_1" / "DeFluorX_Figure_1.html"
_V = RUN / "3_Validation_Figures"
_D = RUN / "4_Dendrogram"            # per-tier PFAS x protein matrix-export PNGs (trees themselves are interactive HTML)
_T = RUN / "5_TopN_and_Preparation"
_P = RUN / "6_Physics_Validation"
_M = RUN / "7_QMMM_Defluorination"

# Sankey funnel promoted to Figure 1 (beside the workflow schematic, main() builds it): reserved here
# so the PFAS-scope supplementary glob does not reclaim it once it leaves the PFAS main figure.
SANKEY_PANEL = _V / "07_PFAS_Scope_and_Synthesis" / "10_Sankey_Workflow.png"

# -----------------------------------------------------------------------------
# MAIN figures - curated composites (Figure 02 onward; Figure 01 = workflow + Sankey, built in main()).
# Panel order below IS the montage order, so caption letters follow the list positions directly.
# Layout maps one figure per Abstract claim (Nature-Comms style): the confidence-vs-competence hook is
# its own Figure 3, the fluorine-interaction panels fold into Supplementary 3, and preparation controls
# move to Supplementary 6. The used-set (main() reserves these before the supplementary globs run) makes
# every move a genuine swap, so no panel is duplicated.
# -----------------------------------------------------------------------------
_ROWBREAK = object()   # sentinel in a panel list: force a new montage row here (manual art-direction)

MAIN_FIGS = [
    # Fig 2: the screen and its inputs - a homogeneous, phylogenetically diverse library with the
    # catalytic triad mapped across every variant, a broad competence-tier spread, uniformly high AI
    # confidence. Residue-mapping coverage + phylogeny promoted from Supplementary 1 (foundational to
    # the screen, and they fill the figure).
    ("Screening_and_Model_Quality", [
        SRC_ROOT / "C_INP_Merged_for_Boltz-2.png",
        _V / "03_Dataset_and_Alignment_Overview" / "01_Active_Site_Residue_Mapping_Coverage.png",
        _V / "03_Dataset_and_Alignment_Overview" / "05_Evolutionary_Phylogeny.png",
        _V / "03_Dataset_and_Alignment_Overview" / "02_Tier_Distribution.png",
        _V / "03_Dataset_and_Alignment_Overview" / "04_Tier_Grade_Distribution.png",
        _V / "04_AI_Confidence_Quality" / "01_AI_Quality_Assessment.png",
        _V / "04_AI_Confidence_Quality" / "03_pTM_vs_ipTM_by_Tier.png",
    ]),
    # Fig 3 (NEW): AI co-folding confidence does not predict catalytic competence - the methods hook.
    # Conflict composition (Hidden Gem / Decoy) + confidence-vs-SN2 landscape + multi-model consensus/agreement.
    ("Confidence_Does_Not_Predict_Competence", [
        _V / "07_PFAS_Scope_and_Synthesis" / "05_Conflict_Composition.png",
        _V / "07_PFAS_Scope_and_Synthesis" / "04_Conf_SN2_Landscape.png",
        _V / "08_Diagnostic_and_MultiModel_Trends" / "03_MultiModel_Consensus_by_Tier.png",
        _V / "08_Diagnostic_and_MultiModel_Trends" / "08_Model_Agreement.png",
        _V / "08_Diagnostic_and_MultiModel_Trends" / "05_Quality_and_Competence_Diagnostics.png",
        # robustness of the geometry call across the Boltz-2 model ensemble (two subset views)
        _V / "08_Diagnostic_and_MultiModel_Trends" / "13_MultiModel_Geometry_Variance_by_Model.png",
        _V / "08_Diagnostic_and_MultiModel_Trends" / "14_MultiModel_Geometry_Variance_conf90.png",
    ]),
    # Fig 4: mechanism-first geometry separates what confidence cannot. Score/angle/geometry + the gate
    # UpSet (graded score adds information beyond the binary gates) + fluorine engagement by tier.
    ("Catalytic_Mechanism_and_SN2_Geometry", [
        _V / "05_Catalytic_Geometry_and_Mechanism" / "05_Mechanistic_Score_by_Tier.png",
        _V / "05_Catalytic_Geometry_and_Mechanism" / "07_SN2_Angle_by_Tier.png",
        _V / "05_Catalytic_Geometry_and_Mechanism" / "08_Mechanism_Geometry_Scatter.png",
        _V / "05_Catalytic_Geometry_and_Mechanism" / "01_ActiveSite_RMSD_by_Tier.png",
        _V / "06_Ligand_Interactions_and_Chemical_Space" / "03_Fluorine_Engagement_by_Tier.png",
        _V / "05_Catalytic_Geometry_and_Mechanism" / "14_Mechanistic_Breakdown_by_Tier.png",
        # 06_Mechanistic_Fingerprint (radar) dropped from main - echoes the Figure 5 radars; it stays in
        # Supplementary Figure 2 via the glob.
    ]),
    # Fig 5: the screen recovers FAcD chemistry unprompted, and names the leads. Success + chain-length +
    # size landscape/composition + carbon-vs-confidence (emergent size preference) + Tier-1A lead heatmap.
    # Fig 5: PFAS scope and lead candidates - which enzymes degrade which PFAS, and the elite pairs.
    # The molecular-SIZE / chain-length preference panels (11-14) are intentionally NOT here - they are
    # extracted for a separate manuscript (see _DROP_PANELS) - so this figure is pure scope + leads.
    ("PFAS_Scope_and_Lead_Candidates", [
        _V / "07_PFAS_Scope_and_Synthesis" / "03_Tier_Success_Rates.png",
        _V / "07_PFAS_Scope_and_Synthesis" / "15_Tier1A_Cross_Ligand_Heatmap.png",
        # 16_Swimmer_Top_Per_Tier dropped from main (reads oddly) - stays in Supplementary Figure 1.
        _V / "07_PFAS_Scope_and_Synthesis" / "01_Radar_TopHits.png",
        _V / "07_PFAS_Scope_and_Synthesis" / "02_Radar_TierReps.png",
        _V / "07_PFAS_Scope_and_Synthesis" / "08_Top25_Multitarget_Proteins.png",
        _V / "07_PFAS_Scope_and_Synthesis" / "18_Pocket_Residue_Reaction_Readiness.png",
        _V / "07_PFAS_Scope_and_Synthesis" / "17_Circos_Ligand_Tier_Assignment.png",
    ]),
    # Fig 6: reactive pose stable, near-attack populated, hydration + binding characterised. Solvation QC,
    # MD QC, WaterMap landscape (the hydration panel §3.8 needs), MM-GBSA and defluorination summaries.
    # (Physics folder was renumbered when WaterMap became 02; MD QC is 03, MM-GBSA 04, defluorination 05.)
    ("Physics_MD_and_Hydration", [
        _P / "06_Analysis" / "01_Physics_Build_Solvation_QC.png",
        _P / "06_Analysis" / "03_MD_Trajectory_QC.png",
        _P / "06_Analysis" / "02_WaterMap_Landscapes_AllRanks.png",
        _P / "06_Analysis" / "04_MMGBSA_Combined_AllRanks.png",
        _P / "06_Analysis" / "05_Defluorination_Combined_AllRanks.png",
        # active-site ligand-protein 2D interaction maps (Desmond, contact % per residue) for the four
        # simulated systems - MD results, so they live here rather than in a standalone figure.
        _P / "05_MD_Simulations" / "Analysis_MDS_R_1" / "images" / "LP-Contacts_2d-Summary.png",
        _P / "05_MD_Simulations" / "Analysis_MDS_R_2" / "images" / "LP-Contacts_2d-Summary.png",
        _P / "05_MD_Simulations" / "Analysis_MDS_R_8" / "images" / "LP-Contacts_2d-Summary.png",
        _P / "05_MD_Simulations" / "Analysis_MDS_R_69" / "images" / "LP-Contacts_2d-Summary.png",
    ]),
    # QM/MM validation of C-F cleavage (Figure 07): the reaction-mechanism scene, the all-jobs PES
    # overlay and the four per-rank relaxed-scan profiles, all emitted by the pipeline (Step 07) as SVG.
    # Mechanism figure: (a) the cross-enzyme defluorination-cascade scheme (built in main() before the
    # montage), (b) the all-jobs QM/MM PES overlay. The four per-rank reaction profiles move to the
    # supplementary ("QSite_Reaction_Profiles") so this stays a two-panel main figure.
    ("Mechanism_and_QSite_Validation", [
        OUT_DIR / "_cascade" / "Reaction_Cascade.svg", _ROWBREAK,   # a: cascade, full-width top row
        _M / "08_QSite_Profiles_AllJobs.png",                       # b: all-jobs PES, full-width below
    ]),
    # Cross-rank defluorination outcome (promoted from the old Supplementary 7 - these are key results,
    # not diagnostics): which systems are viable, the defluorination landscape, the comparative analysis
    # and residue engagement, plus the cross-rank MM-GBSA decomposition and machinery engagement.
    # All six panels are wide strips, so they are stacked full-width (A4 portrait); the residue-engagement
    # heatmap (d), being the widest, then occupies the shortest band - i.e. the smallest panel.
    ("Cross_Rank_Defluorination_Outcome", [
        # Art-directed: C (Comparative) and F (Machinery) each get a full-width row (the key results,
        # most space); D (Residue engagement, least important) shares a row with E so it stays small.
        _M / "03_Viability_Summary.png",
        _M / "05_Defluorination_Landscape.png", _ROWBREAK,
        _M / "02_Comparative_Analysis.png", _ROWBREAK,                 # c - big, full width
        _M / "04_Comparative_Residue_Engagement.png",                  # d - small (shares row with e)
        _M / "06_MMGBSA_Decomposition_AllRanks.png", _ROWBREAK,
        _M / "07_Machinery_Engagement_AllRanks.png",                   # f - big, full width
    ]),
]

# -----------------------------------------------------------------------------
# SUPPLEMENTARY glob groups - one merged figure per group; each can draw from several source dirs so no
# figure is ever a lone panel. (figure_name, [(base_dir, glob), ...]).
# The control Ramachandran (3_Validation/02_Ramachandran) is deliberately NOT collected - the
# MD-selected comparative Ramachandran (5_TopN/3_Comparative) covers it; both would duplicate.
# From 5_TopN only the prepared-pose interaction maps, comparative analysis and ESP charges are kept;
# the 1_Converted raw renders and the PLIP/PyMOL renders are left out.
# Per-rank 06/07 panels are NOT globbed here - they are assembled into per-rank figures below, so the
# only 6_Physics/7_MD panels a glob touches are the cross-rank ones under 7_MD's top level.
# -----------------------------------------------------------------------------
# Panels EXCLUDED from every manuscript figure (neither a main panel nor swept into a supplementary
# glob). Two reasons: (1) the molecular-SIZE / chain-length preference story is a SEPARATE manuscript
# (extracted to Future_Paper_Size_Preference/); (2) a few low-value panels dropped on review.
_DROP_PANELS = (
    # size-preference paper
    "11_PFAS_Size_Hexbin_Landscape", "12_PFAS_Size_Composition_Merged", "13_PFAS_Carbon_Confidence",
    "14_Chain_Length_by_Tier", "01_Pocket_vs_Ligand_Volume", "02_Pocket_Occupancy_by_Carbon_Number",
    "06_Size_Preference_Containment",
    # dropped as low-value (Manhattan duplicates the Volcano significance panel)
    "11_Manhattan_Candidate_Significance",
    # prepared-pose Ramachandran comparison - not needed in the manuscript
    "Rama_Comparison",
)

SUPP_GROUPS = [
    # Old thin "Dataset composition" (Supp 1) and "PFAS-scope diagnostics" (Supp 4) merged into one
    # screening/scope diagnostics figure so neither is a near-empty standalone supplementary.
    ("Screening_and_Scope_Diagnostics", [
        (_V / "03_Dataset_and_Alignment_Overview", "*.png"),
        (_V / "04_AI_Confidence_Quality", "*.png"),
        (_D, "**/*.png"),
        (_V / "07_PFAS_Scope_and_Synthesis", "*.png"),
    ]),
    ("Catalytic_Mechanism_Diagnostics", [
        (_V / "05_Catalytic_Geometry_and_Mechanism", "*.png"),
    ]),
    ("Ligand_Interactions_and_Chemical_Space", [
        (_V / "06_Ligand_Interactions_and_Chemical_Space", "*.png"),
    ]),
    ("MultiModel_and_Diagnostic_Trends", [
        (_V / "08_Diagnostic_and_MultiModel_Trends", "*.png"),
    ]),
    ("Preparation_ESP_and_Interaction_Maps", [
        (_T / "3_Comparative_Analysis", "**/*.png"),
        (_T / "2_Prepared_PDBs" / "Figures" / "InteractionMap", "*.png"),
        (_T / "4_Ligand_ESP_Charges", "*.png"),
    ]),
    # (The cross-rank QSite / defluorination summaries are a MAIN figure - see MAIN_FIGS
    # 'Cross_Rank_Defluorination_Outcome' - because they are key results, not diagnostics.)
]

# Per-supplementary ROW layout: name -> list of rows, each a list of filename-stem fragments. The
# collected panels are reordered into these rows (a _ROWBREAK between rows); any unmatched panel is
# appended as a final row. Lets a globbed supplementary be art-directed without losing the glob.
_SUPP_LAYOUT = {
    # ESP (01_ESP..) is very wide and was landing alone at the bottom; pair it with the other wide
    # panel (Pose drift). Interaction maps fill the lower rows.
    "Preparation_ESP_and_Interaction_Maps": [
        ["02_Pose_Drift", "01_ESP_Alpha_Carbon_Charge"],
        ["03_Machinery_Engagement_Distribution", "3R3U_Control"],
        ["A0A2U3PT06", "METTU", "GJE15807"],
    ],
}

# Panel source substitutions: a run-generated panel whose stem matches is replaced by a locally
# cropped/fixed copy under Figures_HQ/_assets (the run source is never modified).
_PANEL_SUBSTITUTE = {
    # Browser-exported dendrogram SVG had a too-tall viewBox (content top, big empty bottom);
    # use the viewBox-cropped copy so it fills its panel instead of floating with white space.
    "Tier_1A_PFAS_Matrix_Export": OUT_DIR / "_assets" / "Tier_1A_PFAS_Matrix_Export.svg",
}

def _substitute_panel(p: "Path") -> "Path":
    for _stem, _sub in _PANEL_SUBSTITUTE.items():
        if _stem in p.name and _sub.exists():
            return _sub
    return p

def _apply_supp_layout(name: str, pngs: list) -> list:
    rows = _SUPP_LAYOUT.get(name)
    if not rows:
        return pngs
    remaining = list(pngs); out = []
    for ri, row in enumerate(rows):
        placed_any = False
        for frag in row:
            hit = next((p for p in remaining if frag in p.name), None)
            if hit is not None:
                out.append(hit); remaining.remove(hit); placed_any = True
        if placed_any and ri < len(rows) - 1:
            out.append(_ROWBREAK)
    if remaining:                       # any panel the layout did not name -> its own final row
        if out and out[-1] is not _ROWBREAK:
            out.append(_ROWBREAK)
        out.extend(remaining)
    return out

# -----------------------------------------------------------------------------
# PER-RANK deep-dive figures - one merged figure per MD rank, in panel reading order:
# structure/pose -> dynamics -> energetics -> QSite reaction. The last two (QSite) are placeholders
# that do not exist until Step 07 finishes; they are reserved here so they slot in on a re-run.
# -----------------------------------------------------------------------------
def _rank_panels(rank_dir: Path, n: str) -> list:
    # Note the differing per-rank dir/file conventions across steps:
    #   06 defluorination : Defluorination_R{n}          07 MD : Rank_{n} (= rank_dir)
    #   06 Desmond SID     : Analysis_MDS_R_{n}/images    06 MM-GBSA : Prime-MMGBSA/MMGBSA_Profile_R{n}
    defl = _P / "06_Analysis" / "Defluorination" / f"Defluorination_R{n}"
    # The following per-rank panels are deliberately OMITTED because a cross-rank version already tells
    # the same story once, so per-rank copies would duplicate it:
    #   RMSD / RMSF (PL-RMSD, P-RMSF, L-RMSF)       -> Figure 6 panel b (03_MD_Trajectory_QC)
    #   MM-GBSA profile (MMGBSA_Profile_R{n})        -> Figure 6 panel d (04_MMGBSA_Combined_AllRanks)
    #   MM-GBSA decomposition (06_..._and_Trace)     -> Supplementary Figure 7 (06_MMGBSA_Decomposition_AllRanks)
    #   machinery engagement (07_Machinery...)       -> Supplementary Figure 7 (07_Machinery_Engagement_AllRanks)
    #   QSite reaction profile (08_QSite_...)        -> Figure 7 (per-system profiles)
    # The Desmond SID MD-stability panels (SSE timeline/histogram, contacts timeline/histogram, ligand
    # 2D schematic, ligand properties) are dropped: RMSD/RMSF are in Figure 6b and the ligand-protein 2D
    # interaction map in Figure 6f-i, so those per-rank panels add nothing. Only the reactive-pose
    # trajectory is kept from the old Part 1, and it folds into this single per-rank figure.
    # The per-system near-attack / competence story, trimmed to the five panels that each add a distinct
    # view. Dropped as redundant/low-signal: 04_Active_Site_Dynamics (interaction distances over time -
    # overlaps the reactive-pose panel and the cross-rank engagement) and 10_MD_QSite_Timeline (attack
    # angle over time - already carried by the reactive-pose series and the NAC-dashboard density).
    # Order puts the two crucial panels (reactive-pose + NAC) first so they share a 2-panel row and stay
    # large; the free-energy landscapes go LAST so they cap to a smaller lone-last slot rather than hog
    # space.
    # This rank's prepared-pose comparative Ramachandran (variable job-id in the filename -> glob).
    _rama_dir = _T / "3_Comparative_Analysis" / "04_MD_Selected_4hits_Comparative_Ramachandran_Plots"
    _rama = sorted(_rama_dir.glob(f"Rank_{n}_*_Rama_Comparison.svg")) if _rama_dir.exists() else []
    # Row-controlled 2-per-row layout: a+b crucial pair (row 1); c+d (row 2); e WaterMap + f free-energy
    # (row 3); g Ramachandran + h QM/MM reaction profile (row 4). The QM/MM profile (h) is wide, so in a
    # shared two-panel row it still takes the larger share of the width - prominent, without the excessive
    # height a lone full-width band gave it.
    qsite = _M / f"Rank_{n}" / "08_QSite_Reaction_Profile.png"            # h: QM/MM PES for this rank
    panels = [
        defl / "01_Reactive_Pose_Trajectory.png",             # a: reactive-pose geometry over the run
        rank_dir / "03_NAC_Dashboard.png", _ROWBREAK,         # b: near-attack basin + anchoring  (row 1 = a,b)
        defl / "02_Reactive_Summary.png",                     # c: per-residue occupancy + NAC-criteria summary
        defl / "03_Binding_vs_Reactivity.png", _ROWBREAK,     # d: binding (MM-GBSA) vs reactivity  (row 2 = c,d)
        _P / "03_WaterMaps" / f"watermap_R_{n}.png",          # e: WaterMap hydration render
        rank_dir / "05_Free_Energy_Landscapes.png", _ROWBREAK,  # f: free-energy landscapes  (row 3 = e,f)
    ]
    if _rama:
        panels.append(_rama[0])                               # g: prepared-pose Ramachandran (row 4, beside h)
    panels.append(qsite)                                      # h: QM/MM reaction profile (wide -> dominant in row 4)
    return panels


_EXCLUDE_PARTS = {"01_Analysis_Data"}      # data CSV folder - never a figure
# Deliberate coverage-audit skips: raw pre-prep renders and the per-engine render dumps.
# These per-rank panels are deliberately unplaced because a cross-rank version tells the story once:
#   PL-RMSD / P-RMSF / L-RMSF           -> Figure 6b (03_MD_Trajectory_QC)
#   MMGBSA_Profile_R{n}                 -> Figure 6d (04_MMGBSA_Combined_AllRanks)
#   06_MMGBSA_Decomposition_and_Trace   -> Supplementary Figure 7 (06_MMGBSA_Decomposition_AllRanks)
#   07_Machinery_Engagement.           -> Supplementary Figure 7 (07_Machinery_Engagement_AllRanks)
# The per-rank MD-stability panels below are dropped entirely (RMSD/RMSF -> Fig 6b, ligand-protein 2D
# interaction -> Fig 6f-i), so SSE / contacts timeline+histogram / ligand-2D / ligand-properties are
# deliberately unplaced too.
_AUDIT_SKIP = ("1_Converted_Raw_PDB", "02_Ramachandran", "/PLIP/", "/PyMOL/", "/Figures/PyMOL",
               "PL-RMSD", "P-RMSF", "L-RMSF", "MMGBSA_Profile_R",
               "06_MMGBSA_Decomposition_and_Trace", "07_Machinery_Engagement.",
               "P-SSE_Timeline", "P-SSE_Histogram", "PL-Contacts_Timeline", "PL-Contacts_Histogram",
               "L_2d_main", "L-Properties", "04_Active_Site_Dynamics", "10_MD_QSite_Timeline")

# ---- montage geometry (portrait, Nature-style; justified rows to minimise white space) ----
TARGET_W = 2400        # figure width (px); rows are scaled to exactly this -> no ragged right edge
NOMINAL_H = 900        # target row height (px) - larger => fewer, bigger panels per row => taller/portrait
PAD = 14               # gap between panels (px)
A4_RATIO = 297.0 / 210.0    # PORTRAIT A4 target height/width (~1.414): every merged figure is laid out to
                            # fill one A4 portrait sheet with minimal white space.
FONT_BOOST = 1.00           # merger-level enlargement of every nested-panel <text> (labels are editable,
                            # so this is applied in place). Set to 1.0 to disable.
# Panels that are already finalised and must NOT have their fonts touched (by source stem).
_FONT_NO_BOOST = {"C_INP_Merged_for_Boltz-2"}
# Per-figure manual nudges of the panel LETTER (a, b, ...) when it sits on panel content. Keyed by output
# figure stem -> {panel index (0=a): dx in px (negative = left, positive = right)}.
# Each entry is {panel index (0=a): dx | (dx, dy)} in px (negative = left/up). Usually empty: the panel
# content is inset down-right so the letter sits in a clear top-left margin; add a nudge only for a panel
# whose own art still reaches into that corner.
_LETTER_OFFSETS = {}
LABEL_RGB = (0, 114, 178)   # panel-letter colour: Okabe-Ito blue #0072B2 (Nature colourblind-safe palette)


def _partition_rows(ars: list, r: int) -> list:
    """Split the ORDERED panels into r contiguous rows so the row aspect-sums are as balanced as
    possible (minimise the largest row's total aspect ratio). Order is preserved, so the caption
    letters a, b, c ... stay in reading order. Classic 'split array into k parts, minimise max sum' DP."""
    n = len(ars)
    pre = [0.0]
    for a in ars:
        pre.append(pre[-1] + a)
    INF = float("inf")
    dp = [[INF] * (r + 1) for _ in range(n + 1)]
    cut = [[0] * (r + 1) for _ in range(n + 1)]
    dp[0][0] = 0.0
    for i in range(1, n + 1):
        for k in range(1, min(r, i) + 1):
            for j in range(k - 1, i):
                val = max(dp[j][k - 1], pre[i] - pre[j])
                if val < dp[i][k]:
                    dp[i][k] = val; cut[i][k] = j
    rows, i, k = [], n, r
    while k > 0:
        j = cut[i][k]
        rows.append(list(range(j, i)))
        i, k = j, k - 1
    rows.reverse()
    return rows


def _pack_rows(ars: list, forced: list = None):
    """Lay the panels out as justified rows that FILL the page (like a packed collage), fitting one A4
    portrait sheet. The row count is chosen so the figure aspect is close to A4 AND so no row is a lone
    panel (a single panel in a row is what leaves the big side white-space / tiny-bottom look). Every row
    is then justified to the FULL width - no row is shrunk or centred - so the sheet is fully filled; if
    the whole thing is taller than A4 it is scaled down uniformly. Returns (layout, H).
    `forced` (a list of index-rows, from explicit _ROWBREAK markers) overrides the automatic row search
    so a figure can be art-directed - e.g. crucial panels alone in a short row so they stay large."""
    n = len(ars)
    if n == 0:
        return [], float(PAD)

    if forced is not None:
        rows = forced
        hs = [(TARGET_W - PAD * (len(row) + 1)) / sum(ars[i] for i in row) for row in rows]
        best = (0.0, rows, hs)
    else:
        best = None
    _H_MAX = TARGET_W * A4_RATIO                           # one A4 portrait page (full width)
    for r in ([] if forced is not None else range(1, n + 1)):
        rows = _partition_rows(ars, r)
        hs = [(TARGET_W - PAD * (len(row) + 1)) / sum(ars[i] for i in row) for row in rows]
        H = sum(hs) + PAD * (len(rows) + 1)
        # Every row is justified to the FULL width, so width is always filled (no side white-space). The
        # goal is a page that FITS one A4 portrait sheet (H <= H_MAX) and fills as much of it as possible:
        #  - overflow past A4 height is penalised hard (it would force a shrink -> side margins);
        #  - among layouts that fit, prefer the one whose height is closest to a full A4 page;
        #  - a lone panel in a row (full-width band) is nudged against.
        lone = sum(1 for row in rows if len(row) == 1)
        over = max(0.0, H - _H_MAX) / _H_MAX
        fill_gap = max(0.0, _H_MAX - H) / _H_MAX
        score = 6.0 * over + fill_gap + 0.18 * lone
        if best is None or score < best[0]:
            best = (score, rows, hs)
    _, rows, hs = best

    raw = [[row, h] for row, h in zip(rows, hs)]      # EVERY row justifies to the full width (incl. a
                                                      # lone last panel) so there is no side white-space
    sum_h = sum(h for _, h in raw)
    pad_total = PAD * (len(raw) + 1)
    H = sum_h + pad_total
    H_MAX = TARGET_W * A4_RATIO
    # Auto layouts are fitted to one A4 portrait page. For FORCED (art-directed) rows the author chose the
    # stacking deliberately, so we keep every row justified to the FULL width (no side white-space) and let
    # the page run taller than A4 rather than shrinking the panels inward.
    if forced is None and H > H_MAX and sum_h > 0:    # too tall -> scale down uniformly to one A4 page
        s = (H_MAX - pad_total) / sum_h
        for row in raw:
            row[1] *= s
        H = H_MAX
    return [(row, h) for row, h in raw], H


def _font(size: int, bold: bool = False):
    try:
        from matplotlib.font_manager import FontProperties, findfont
        path = findfont(FontProperties(family="DejaVu Sans", weight="bold" if bold else "normal"))
        return ImageFont.truetype(path, size)
    except Exception:
        return ImageFont.load_default()


def _montage_png(paths: list, out_path: Path) -> None:
    """Justified-row (masonry) raster montage (legacy FIG_FORMAT='png' fallback): read PNGs in place,
    keep each panel's aspect ratio, and scale every row to the full figure width so panels pack tight
    with minimal white space. Portrait-biased. Panel letters (a, b, ...) on each panel's top-left."""
    # Read panels defensively: a source PNG still being written (or corrupt) must never crash the
    # whole montage - skip it with a warning and tile the rest.
    imgs, kept = [], []
    for p in paths:
        try:
            imgs.append(Image.open(p).convert("RGB"))
            kept.append(p)
        except Exception as e:
            _say(f"  !! skipped unreadable panel: {p.name} ({type(e).__name__})")
    paths = kept
    n = len(imgs)
    if n == 0:
        return
    ars = [im.width / im.height for im in imgs]

    # A4-targeted justified layout (row count chosen so the figure fits one A4 portrait page).
    layout, Hf = _pack_rows(ars)
    rows = [row for row, _ in layout]

    H = int(round(Hf))
    canvas = Image.new("RGB", (TARGET_W, H), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    overlay = Image.new("RGBA", (TARGET_W, H), (0, 0, 0, 0))   # panel labels drawn here, then alpha-composited
    odraw = ImageDraw.Draw(overlay)
    idx = 0; y = PAD
    for row, h in layout:
        hh = int(round(h))
        ws = [max(1, int(round(ars[i] * h))) for i in row]
        row_w = sum(ws) + PAD * (len(row) + 1)
        x = PAD + max(0, (TARGET_W - row_w) // 2)   # centre a narrower row so any gap is symmetric, not bottom-right
        for i, w in zip(row, ws):
            canvas.paste(imgs[i].resize((w, hh), Image.LANCZOS), (x, y))
            # Bold panel letter in a white badge with a dark border: reads instantly on any panel
            # background (light or dark), unlike a thin coloured glyph.
            lab = chr(97 + idx) if idx < 26 else str(idx + 1)
            fsz = max(30, int(h * 0.05))
            lfont = _font(fsz, bold=True)
            _tb = draw.textbbox((0, 0), lab, font=lfont)
            _tw, _th = _tb[2] - _tb[0], _tb[3] - _tb[1]
            _pad = max(6, int(fsz * 0.30))
            _bx, _by = x + 8, y + 8
            odraw.rounded_rectangle([_bx, _by, _bx + _tw + 2 * _pad, _by + _th + 2 * _pad],
                                    radius=int(fsz * 0.28), fill=(255, 255, 255, 115),
                                    outline=(206, 34, 34, 240), width=2)
            odraw.text((_bx + _pad - _tb[0], _by + _pad - _tb[1]), lab, fill=(15, 15, 15, 255), font=lfont)
            idx += 1; x += w + PAD
        y += hh + PAD
    canvas = Image.alpha_composite(canvas.convert("RGBA"), overlay).convert("RGB")
    canvas.save(out_path, dpi=(300, 300))
    _say(f"  ok  {out_path.name}  ({n} panels, {len(rows)} rows)")


# ----------------------------------------------------------------------------------------------------
# True-vector SVG compositing
# ----------------------------------------------------------------------------------------------------
_SVG_MAX_BYTES = 4 * 1024 * 1024    # a panel SVG larger than this is dense (per-point vectors); embed it as raster rather than millions of vector elements
_RASTER_DIR = Path(tempfile.mkdtemp(prefix="dfx_raster_"))
def _find_rsvg() -> str:
    """Locate rsvg-convert even when the active conda env (e.g. PFAS) does not carry it on PATH; the deck
    rasteriser and the oversized-panel fallback both need it, and it ships with the base miniconda env."""
    cand = shutil.which("rsvg-convert")
    if cand:
        return cand
    for p in (Path.home() / "miniconda3" / "bin" / "rsvg-convert",
              Path.home() / "anaconda3" / "bin" / "rsvg-convert",
              Path("/usr/local/bin/rsvg-convert"), Path("/usr/bin/rsvg-convert")):
        if p.exists():
            return str(p)
    for p in sorted(Path("/opt").glob("*conda*/bin/rsvg-convert")):
        return str(p)
    return "rsvg-convert"      # last resort; the callers handle a missing renderer gracefully
_RSVG_BIN = _find_rsvg()

_RASTER_CACHE: dict = {}

def _rasterise_svg(svg: Path, width: int = 1800):
    """A compact PNG of a dense panel SVG (rsvg-convert), so the montage embeds a small image instead of
    the panel's millions of vector elements. Cached per source path; returns None if the renderer fails."""
    key = str(svg)
    if key not in _RASTER_CACHE:
        out = _RASTER_DIR / (svg.stem + "_r.png")
        try:
            subprocess.run([_RSVG_BIN, "-b", "white", "-w", str(width), str(svg), "-o", str(out)],
                           check=True, timeout=240)
            _RASTER_CACHE[key] = out
        except Exception as e:
            _say(f"  !! could not rasterise oversized panel {svg.name} ({type(e).__name__})")
            _RASTER_CACHE[key] = None
    return _RASTER_CACHE[key]

def _resolved(p: Path) -> Path:
    """Prefer a panel's .svg sibling over its .png (the pipeline emits SVG via the CFG SSOT), UNLESS that
    .svg is large: dense per-point panels (Desmond SID timelines, 58k-point scatters) balloon a merged
    figure past submission limits, so an oversized panel falls back to its small PNG sibling or, failing
    that, an on-the-fly raster. PyMOL/dendrogram panels have no .svg and stay raster."""
    if FIG_FORMAT == "svg":
        svg = p.with_suffix(".svg")
        if svg.exists():
            if svg.stat().st_size <= _SVG_MAX_BYTES:
                return svg
            png = p.with_suffix(".png")
            if png.exists() and png.stat().st_size <= _SVG_MAX_BYTES:
                return png
            r = _rasterise_svg(svg)
            return r if r is not None else svg
    return p


def _svg_viewbox_wh(root) -> tuple[float, float]:
    """Intrinsic (viewBox) width/height of an SVG root, so a nested panel keeps its true aspect ratio."""
    vb = root.get("viewBox")
    if vb:
        parts = [float(x) for x in vb.replace(",", " ").split()]
        if len(parts) == 4 and parts[2] > 0 and parts[3] > 0:
            return parts[2], parts[3]
    def _u(v):
        return float(re.sub(r"[a-z%]+$", "", v)) if v else 0.0
    w, h = _u(root.get("width")), _u(root.get("height"))
    return (w or 1.0), (h or 1.0)


def _svg_namespace_ids(root, prefix: str) -> None:
    """Prefix every id and every reference to it so panels never collide on shared matplotlib ids
    (clip-paths, gradient/pattern defs, glyph <use> refs) once nested into one parent SVG."""
    ids = {el.get("id") for el in root.iter() if el.get("id")}
    if not ids:
        return
    def _sub_url(text: str) -> str:
        for _id in ids:
            text = text.replace(f"url(#{_id})", f"url(#{prefix}{_id})")
        return text
    xlink_href = f"{{{XLINK_NS}}}href"
    for el in root.iter():
        _id = el.get("id")
        if _id in ids:
            el.set("id", prefix + _id)
        for attr in _SVG_IDREF_ATTRS:
            v = el.get(attr)
            if v and "url(#" in v:
                el.set(attr, _sub_url(v))
        st = el.get("style")
        if st and "url(#" in st:
            el.set("style", _sub_url(st))
        for href_attr in (xlink_href, "href"):
            hv = el.get(href_attr)
            if hv and hv.startswith("#") and hv[1:] in ids:
                el.set(href_attr, "#" + prefix + hv[1:])


def _svg_aspect(path: Path) -> float:
    """Panel aspect ratio (w/h): from the SVG viewBox for vector panels, from pixel size for raster."""
    if path.suffix.lower() == ".svg":
        w, h = _svg_viewbox_wh(etree.parse(str(path)).getroot())
        return w / h
    with Image.open(path) as im:
        return im.width / im.height


def _boost_fonts(root, f: float) -> None:
    """Enlarge every <text>/<tspan> font-size by factor f in place (labels are editable <text>)."""
    if f == 1.0:
        return
    for el in root.iter():
        tag = el.tag if isinstance(el.tag, str) else ""
        if tag.endswith("}text") or tag.endswith("}tspan"):
            fs = el.get("font-size")
            if fs:
                try:
                    el.set("font-size", f"{float(re.sub('[a-z%]+$', '', fs)) * f:.2f}")
                except ValueError:
                    pass
            st = el.get("style")
            if st and "font-size" in st:
                el.set("style", re.sub(r"font-size:\s*([0-9.]+)([a-z%]*)",
                                       lambda m: f"font-size:{float(m.group(1)) * f:.2f}{m.group(2)}", st))


def _svg_place(path: Path, x: float, y: float, w: float, h: float, idx: int) -> etree._Element:
    """One placed panel: a nested live <svg> for a vector source, an <image> for a raster source."""
    if path.suffix.lower() == ".svg":
        panel = etree.parse(str(path)).getroot()
        vbw, vbh = _svg_viewbox_wh(panel)
        # Preserve the source viewBox ORIGIN (min-x/min-y). Matplotlib panels start at 0,0, but a
        # browser-exported SVG (e.g. the dendrogram) can have a negative origin; dropping it shifted the
        # content (top clipped, white padding below). Keep the real origin so content fills the panel.
        _vb = panel.get("viewBox")
        _minx, _miny = 0.0, 0.0
        if _vb:
            _pp = [float(v) for v in _vb.replace(",", " ").split()]
            if len(_pp) == 4:
                _minx, _miny = _pp[0], _pp[1]
        _svg_namespace_ids(panel, f"pnl{idx}_")
        _boost_fonts(panel, 1.0 if path.stem in _FONT_NO_BOOST else FONT_BOOST)
        nested = etree.Element(f"{{{SVG_NS}}}svg")
        nested.set("x", f"{x:.2f}"); nested.set("y", f"{y:.2f}")
        nested.set("width", f"{w:.2f}"); nested.set("height", f"{h:.2f}")
        nested.set("viewBox", f"{_minx:.3f} {_miny:.3f} {vbw:.3f} {vbh:.3f}")
        nested.set("preserveAspectRatio", "xMidYMid meet")
        for child in panel:
            nested.append(child)
        return nested
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    mime = "image/jpeg" if path.suffix.lower() in (".jpg", ".jpeg") else "image/png"
    img = etree.Element(f"{{{SVG_NS}}}image")
    img.set("x", f"{x:.2f}"); img.set("y", f"{y:.2f}")
    img.set("width", f"{w:.2f}"); img.set("height", f"{h:.2f}")
    img.set("preserveAspectRatio", "xMidYMid meet")
    img.set(f"{{{XLINK_NS}}}href", f"data:{mime};base64,{data}")
    return img


def _svg_label(parent, lab: str, x: float, y: float, fsz: float) -> None:
    """Vector panel letter: a bold near-black glyph with a thin white halo for legibility on any panel
    background, placed top-left. No pill or border, so nothing extends below the letter.

    The halo and the glyph are emitted as two stacked <text> nodes (halo first, glyph second) rather
    than one node carrying both a fill and a stroke under paint-order="stroke". Document order puts the
    glyph on top in every renderer, including the MS Word SVG engine, which ignores paint-order and
    would otherwise paint the white stroke over the glyph and erase the letter."""
    common = {
        "x": f"{x + fsz * 0.35:.2f}", "y": f"{y + fsz * 1.00:.2f}", "text-anchor": "start",
        "font-family": "DejaVu Sans, Arial, sans-serif", "font-weight": "bold",
        "font-size": f"{fsz:.1f}",
    }
    halo = etree.SubElement(parent, f"{{{SVG_NS}}}text")
    for k, v in common.items():
        halo.set(k, v)
    halo.set("fill", "none"); halo.set("stroke", "#ffffff"); halo.set("stroke-width", f"{fsz * 0.14:.2f}")
    halo.text = lab
    glyph = etree.SubElement(parent, f"{{{SVG_NS}}}text")
    for k, v in common.items():
        glyph.set(k, v)
    glyph.set("fill", "#0f0f0f")
    glyph.text = lab


def _montage_svg(paths: list, out_path: Path) -> None:
    """Justified-row (masonry) TRUE-VECTOR montage: identical packing geometry to the raster montage,
    but each panel is placed as a nested <svg> (vector) or <image> (raster source). Panel letters are
    vector. Output is one self-contained .svg that stays crisp at any zoom."""
    # explicit row breaks (manual art-direction) -> remember the group sizes, then drop the sentinels
    sizes = None
    if any(p is _ROWBREAK for p in paths):
        sizes, c = [], 0
        for p in paths:
            if p is _ROWBREAK:
                if c:
                    sizes.append(c); c = 0
            else:
                c += 1
        if c:
            sizes.append(c)
        paths = [p for p in paths if p is not _ROWBREAK]
    kept, ars = [], []
    for p in paths:
        try:
            ars.append(_svg_aspect(p)); kept.append(p)
        except Exception as e:
            _say(f"  !! skipped unreadable panel: {p.name} ({type(e).__name__})")
    paths = kept
    n = len(paths)
    if n == 0:
        return
    forced = None
    if sizes and sum(sizes) == n:                      # honour the explicit rows only if nothing dropped
        forced, i = [], 0
        for s in sizes:
            forced.append(list(range(i, i + s))); i += s
    # A4-targeted justified layout (row count chosen so the figure fits one A4 portrait page) - same
    # geometry as the raster path - unless explicit rows were given.
    layout, H = _pack_rows(ars, forced)
    rows = [row for row, _ in layout]

    root = etree.Element(f"{{{SVG_NS}}}svg", nsmap={None: SVG_NS, "xlink": XLINK_NS})
    root.set("width", f"{TARGET_W}"); root.set("height", f"{H:.2f}")
    root.set("viewBox", f"0 0 {TARGET_W} {H:.2f}")
    bg = etree.SubElement(root, f"{{{SVG_NS}}}rect")
    bg.set("width", f"{TARGET_W}"); bg.set("height", f"{H:.2f}"); bg.set("fill", "#ffffff")
    fsz = TARGET_W * 0.017
    # Each panel box has the SAME aspect as its content (box width = ar * h), so a 'meet' fit leaves no
    # gutter and a corner letter would sit on the y-axis / title. Reserve a small top-left margin by
    # insetting the content down-and-right by one letter strip; the letter then sits in that clear corner.
    lins = fsz * 1.05
    idx = 0; y = float(PAD)
    for row, h in layout:
        ws = [max(1.0, ars[i] * h) for i in row]
        row_w = sum(ws) + PAD * (len(row) + 1)
        x = PAD + max(0.0, (TARGET_W - row_w) / 2)
        for i, w in zip(row, ws):
            root.append(_svg_place(paths[i], x + lins, y + lins, w - lins, h - lins, idx))
            lab = chr(97 + idx) if idx < 26 else str(idx + 1)
            off = _LETTER_OFFSETS.get(out_path.stem, {}).get(idx, 0.0)   # manual per-panel letter nudge
            ldx, ldy = off if isinstance(off, (tuple, list)) else (off, 0.0)
            _svg_label(root, lab, x + ldx, y + ldy, fsz)
            idx += 1; x += w + PAD
        y += h + PAD
    out_path.write_bytes(etree.tostring(root, xml_declaration=True, encoding="utf-8"))
    n_vec = sum(1 for p in paths if p.suffix.lower() == ".svg")
    _say(f"  ok  {out_path.name}  ({n} panels, {len(rows)} rows; {n_vec} vector, {n - n_vec} raster)")


_FIG_MAX_BYTES = 20 * 1024 * 1024   # per-figure submission size ceiling

def montage(paths: list, out_path: Path) -> None:
    """Dispatch to the vector or raster montage by the output extension (driven by FIG_FORMAT). If the
    merged SVG still exceeds the submission size ceiling, rasterise every panel and rebuild so the figure
    meets the limit."""
    if out_path.suffix.lower() == ".svg":
        _montage_svg([p if p is _ROWBREAK else _resolved(p) for p in paths], out_path)
        if out_path.stat().st_size > _FIG_MAX_BYTES:
            big = [p if p is _ROWBREAK else (_rasterise_svg(p.with_suffix(".svg")) or _resolved(p)) for p in paths]
            _montage_svg(big, out_path)
    else:
        _montage_png([p for p in paths if p is not _ROWBREAK], out_path)


def _render_workflow(html_path: Path, out_path: Path, target_w: int = 2126) -> bool:
    """Place the workflow / graphical-overview (Figure 1) into the numbered series. In SVG mode the
    HTML's own vector export (Figure_01_Workflow_Overview.svg) is used verbatim - true vector, crisp at
    any zoom, no rasterising. In PNG mode it falls back to a pre-rendered PNG (the HTML's 'Export
    high-res PNG' button) or, failing that, a headless-firefox screenshot upscaled to 180 mm / 300 dpi.
    A missing export is skipped and folds in on a re-run, like every other panel."""
    # SVG mode: take the curated vector export directly; never rasterise the graphical overview.
    if out_path.suffix.lower() == ".svg":
        # curated vector export (Figure 01), used verbatim; accept it at OUT_DIR (the destination, where
        # the export is copied) or beside the HTML. Never rasterise or regenerate the graphical overview.
        srcs = [html_path.parent / "Figure_01_Workflow_Overview.svg"]
        if out_path not in srcs:
            srcs.insert(0, out_path)
        for _svg in srcs:
            if not _svg.exists():
                continue
            # Cap the workflow to the A4 box (2400 x 2400*A4_RATIO) so it never crosses the A4 boundary;
            # it is vector, so this only sets the width/height attributes (aspect preserved).
            try:
                root = etree.parse(str(_svg)).getroot()
                vbw, vbh = _svg_viewbox_wh(root)
                s = min(TARGET_W / vbw, (TARGET_W * A4_RATIO) / vbh)
                root.set("width", f"{vbw * s:.1f}"); root.set("height", f"{vbh * s:.1f}")
                out_path.write_bytes(etree.tostring(root, xml_declaration=True, encoding="utf-8"))
            except Exception:
                if _svg != out_path:
                    shutil.copyfile(_svg, out_path)
            _say(f"  ok  {out_path.name}  (vector export, A4-fitted)")
            return True
        return False      # no user-supplied workflow SVG present; main() logs the reserved-slot note
    # Prefer a pre-rendered PNG next to the HTML (produced by the HTML's own 'Export high-res PNG'
    # button, or a prior render) - reliable and firefox-free. Headless firefox is only the fallback.
    _pre = html_path.parent / "Figure_01_Workflow_Overview.png"
    if _pre.exists():
        im = Image.open(_pre).convert("RGB")
        if im.width != target_w:
            im = im.resize((target_w, round(im.height * target_w / im.width)), Image.LANCZOS)
        im.save(out_path, dpi=(300, 300))
        _say(f"  ok  {out_path.name}  (from pre-rendered {_pre.name}, {im.width}x{im.height})")
        return True
    if not html_path.exists():
        _say(f"  --  workflow HTML pending: {html_path.name}")
        return False
    tmp = Path(tempfile.mkdtemp())
    # Force software GL: headless firefox on this host otherwise trips a RenderCompositorSWGL
    # framebuffer failure and hangs. llvmpipe + MOZ_DISABLE_GPU renders at 1x device pixels, but the
    # failure is intermittent, so each attempt gets a fresh isolated profile (never collides with a
    # running firefox's default profile) and a stuck headless firefox is reaped before the next try.
    _env = dict(os.environ, LIBGL_ALWAYS_SOFTWARE="1", GALLIUM_DRIVER="llvmpipe",
                MOZ_ENABLE_WAYLAND="0", MOZ_DISABLE_GPU="1")
    try:
        shot = None
        for _attempt in range(1, 5):
            _cand = tmp / f"shot{_attempt}.png"
            # Clear any leftover headless firefox first so this render owns the (default) profile;
            # -f matches the full command line, which safely targets firefox and never this process.
            subprocess.run(["pkill", "-9", "-f", "firefox --headless"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            time.sleep(3)   # let the killed instance release the GL/framebuffer before relaunching
            try:
                subprocess.run(["firefox", "--headless", "--window-size=1712,1780",
                                "--screenshot", str(_cand), html_path.as_uri()],
                               timeout=110, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=_env)
            except subprocess.TimeoutExpired:
                pass
            if _cand.exists() and _cand.stat().st_size > 2000:
                shot = _cand; break
            _say(f"  ..  workflow render attempt {_attempt} produced nothing; retrying")
        if shot is None:
            _say("  !! workflow render failed after 4 attempts (firefox headless SWGL)"); return False
        im = Image.open(shot).convert("RGB")
        # trim the page background (#f2f3f5) down to the figure's own bordered frame
        bg = Image.new("RGB", im.size, (242, 243, 245))
        bbox = ImageChops.difference(im, bg).convert("L").point(lambda v: 255 if v > 6 else 0).getbbox()
        if bbox:
            im = im.crop(bbox)
        if im.width != target_w:
            im = im.resize((target_w, round(im.height * target_w / im.width)), Image.LANCZOS)
        im.save(out_path, dpi=(300, 300))
        _say(f"  ok  {out_path.name}  (workflow HTML render, {im.width}x{im.height})")
        return True
    except Exception as e:
        _say(f"  !! workflow render failed: {type(e).__name__}: {e}")
        return False
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _reset_outputs() -> None:
    """Replace ONLY this step's own figure families so a re-run never leaves stale outputs - and never
    touches anything else in the folder (the user may keep their own Figure 01 panel/folder, notes, etc.
    here). Figure 01 is always the user-owned workflow slot, so any 'Figure_01*' is skipped outright."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for stem in ("Figure_*", "Supplementary_Figure_*", "00_*_Collage"):
        for ext in (".png", ".svg"):
            for p in OUT_DIR.glob(stem + ext):
                if p.name.startswith("Figure_01"):
                    continue   # Figure 01 = user-curated workflow slot; never auto-deleted
                p.unlink()


def _md_ranks() -> list:
    """Discovered MD ranks (7_MD/Rank_N), sorted numerically by N."""
    def _n(p: Path) -> int:
        tail = p.name.split("_", 1)[1]
        return int(tail) if tail.isdigit() else 1_000_000
    return sorted((p for p in _M.glob("Rank_*") if p.is_dir()), key=_n)


def _key(p: Path) -> Path:
    """Suffix-stripped resolved path - a panel's .png and .svg share one key, so a panel is reserved
    and placed exactly once regardless of which format is on disk."""
    return p.resolve().with_suffix("")


# -----------------------------------------------------------------------------
# Figure legends (written to Figure_Legends.md alongside the figures). Keyed by the montage NAME so a
# figure's number can shift (e.g. when the workflow panel is present) without breaking the mapping.
# -----------------------------------------------------------------------------
MAIN_LEGENDS = {
    "Workflow_Overview": "Mechanism-guided discovery workflow. Graphical overview of the DeFluorX pipeline: Boltz-2 co-folding of every bacterial FAcD variant against the 27-PFAS panel, catalytic-geometry tiering, MD-based near-attack sampling, and QM/MM defluorination verdicts, shown beside the screen funnel (Sankey).",
    "Screening_and_Model_Quality": "Screening library and model quality. (a) Composition of the merged FAcD-variant and 27-PFAS input library prepared for Boltz-2 co-folding. (b) Active-site catalytic-residue mapping coverage across the library. (c) Evolutionary phylogeny of the variants. (d) Competence-tier distribution of the screened complexes. (e) Tier-by-grade distribution. (f) AI co-folding quality (pTM, ipTM and confidence). (g) pTM against ipTM by tier.",
    "Confidence_Does_Not_Predict_Competence": "Confidence does not predict competence. AI co-folding confidence is decoupled from mechanistic competence, so a high-confidence pose is not necessarily catalytically competent and ranking is geometry-first rather than confidence-first. (a) Conflict composition (hidden-gem and decoy complexes where confidence and competence disagree). (b) Confidence against SN2-geometry landscape. (c) Multi-model consensus by tier. (d) Boltz-2 model agreement. (e) Joint quality-and-competence diagnostics. (f) Multi-model geometry variance by model. (g) Multi-model geometry variance for the ipTM-90 subset.",
    "Catalytic_Mechanism_and_SN2_Geometry": "Catalytic mechanism and SN2 geometry. The mechanism-first score and its geometric terms separate the elite (Tier_1A) degraders. (a) Mechanistic score by tier. (b) SN2 backside-attack angle by tier. (c) Mechanistic score against geometry. (d) Active-site RMSD by tier. (e) Fluorine engagement by tier. (f) Mechanistic-term breakdown by tier.",
    "PFAS_Scope_and_Lead_Candidates": "PFAS scope and lead candidates. Substrate scope across the 27-PFAS panel and the leads advanced to physics validation. (a) Tier success rates. (b) Tier_1A enzyme-by-PFAS cross-ligand heatmap. (c) Radar of the top hits. (d) Radar of the tier representatives. (e) Top multi-target proteins. (f) Per-PFAS pocket reaction-readiness. (g) Circos map of ligand-tier assignment.",
    "Physics_MD_and_Hydration": "Physics and MD overview for the MD-selected complexes. (a) Build and solvation quality control. (b) Cross-rank trajectory quality control (Ca-RMSD, ligand RMSD, RMSF and temperature). (c) WaterMap hydration landscapes across ranks. (d) Prime MM-GBSA binding free energies across ranks. (e) Cross-rank defluorination summary. (f-i) Ligand-protein 2D interaction maps, in the order R1 (fluoroacetate lead), R2 (difluoroacetate lead), R8 (trifluoroacetate lead) and the 3R3U/R69 control.",
    "Mechanism_and_QSite_Validation": "Enzymatic defluorination mechanism and QM/MM validation. (a) The hydrolytic double-displacement mechanism (Asp110 backside SN2 -> glycolyl-enzyme ester -> His280/water hydrolysis), drawn from the real QM geometry for the three independent lead FAcD variants, each on its own substrate (trifluoroacetate, Methylocella tundrae; difluoroacetate, Methylobacterium marchantiae; monofluoroacetate, Bradyrhizobium sp.), with per-enzyme QM/MM reaction-coordinate plots (scissile C-F length and departing-fluoride charge) and a joined QM/MM-geometry-WaterMap-MD read-out; the enzyme replaces one C-F with C-OH, and for the di- and mono-fluoro alpha-hydroxy products a non-enzymatic HF elimination follows. These are three independent enzymes, one per substrate, not a single-substrate cascade. (b) All-jobs QM/MM potential-energy-surface overlay over the scissile C-F bond length, with the geometric cleavage verdict per system. Per-system relaxed-scan reaction profiles, including the RPA1163/3R3U reference, are panel (h) of the per-system near-attack supplementary figures (Supplementary Figures 6-9).",
    "Cross_Rank_Defluorination_Outcome": "Cross-rank defluorination outcome. Whole-story view linking viability, near-attack persistence, QM/MM energetics and MM-GBSA binding across ranks, with the geometric C-F cleavage verdict per system. (a) Viability summary. (b) Defluorination landscape. (c) Comparative analysis across ranks. (d) Comparative residue engagement. (e) Cross-rank MM-GBSA decomposition. (f) Cross-rank machinery engagement.",
}

SUPP_LEGENDS = {
    "Screening_and_Scope_Diagnostics": "Screening and scope diagnostics. Further characterisation of the input library and the PFAS scope (the catalytic-residue mapping coverage and the evolutionary phylogeny are promoted to Figure 2b,c). (a) Distribution of active-site alignment grades across the library. (b) Treemap of the tier-by-ligand composition of the screened complexes. (c) AI-quality feature space (pTM, ipTM and confidence). (d) Tier_1A enzyme-by-PFAS dendrogram and cross-ligand tier matrix for fluoroacetate, trifluoroacetate and difluoroacetate. (e) Euler diagram of substrate-category overlap. (f) Top-tier protein-by-PFAS breakdown. (g) Swimmer plot of the top hits per tier.",
    "Catalytic_Mechanism_Diagnostics": "Catalytic-mechanism diagnostics. Analysis of the mechanism-first score that expands Figure 4. (a) Pairwise correlations among the mechanistic features. (b) Tier-quality dot plot. (c) Cross-tabulation of discrete mechanistic states. (d) Mechanistic fingerprint of the tiers. (e) Mechanistic-quality feature space. (f) Criterion A gating criterion B. (g) Empirical cumulative distribution of criterion B by tier. (h) SN2 dead-end gate. (i) Joint geometry-and-uncertainty view. (j) Metric co-variation network. (k) Mechanistic-gate intersection (UpSet) analysis.",
    "Ligand_Interactions_and_Chemical_Space": "Ligand interactions and chemical space. Expanded interaction and chemical-space analysis supporting Figures 4 and 5. (a) Molecular interaction profile. (b) Interaction-quality feature space. (c) Catalytic quality against predicted inhibition, separating productive from non-productive binding. (d) Active-site contact density by tier. (e) Binding energetics. (f) Chemical-space map. (g) Binding-affinity metric distributions across the tiers.",
    "MultiModel_and_Diagnostic_Trends": "Multi-model and diagnostic trends. Diagnostic trends across the screen and the robustness of the ranking across the Boltz-2 model ensemble. (a) Confidence against consensus. (b) Reactive engagement. (c) Pillar divergence by tier. (d) Volcano plot of metric significance. (e) Feature-importance tier streams. (f) Multi-model geometry variance for the ipTM-90 subset.",
    "Preparation_ESP_and_Interaction_Maps": "Preparation, ESP and interaction maps. Preparation quality control for the systems advanced to physics validation. (a) Pose drift from the co-fold CIF to the prepared structure. (b) Electrostatic-potential-fitted alpha-carbon charge of the reactive ligands. (c) Distribution of machinery engagement. (d-g) Prepared-complex protein-ligand interaction maps, in the order the 3R3U control, R1 (the fluoroacetate lead), R8 (the trifluoroacetate lead) and R2 (the difluoroacetate lead).",
    "MD_Rank1_Dynamics_and_Reaction": "Near-attack dynamics and C-F cleavage, R1. Per-system near-attack dynamics and the C-F cleavage reaction for R1, the fluoroacetate lead (Bradyrhizobium sp., UniProt A0A2U3PT06), over the one-microsecond trajectory. (a) Reactive-pose trajectory (attack angle, halide-cradle and carboxylate-clamp distances over time). (b) Near-attack-conformation dashboard. (c) Per-residue occupancy and NAC-criteria summary. (d) Binding (MM-GBSA) against reactivity, per frame. (e) WaterMap hydration sites. (f) Free-energy landscapes of the reactive coordinate. (g) Prepared-pose comparative Ramachandran plot. (h) QM/MM relaxed-scan reaction profile over the scissile C-F bond (per-frame and ensemble barrier, reaction energy, departing-fluoride charge and the geometric cleavage verdict).",
    "MD_Rank2_Dynamics_and_Reaction": "Near-attack dynamics and C-F cleavage, R2. Per-system near-attack dynamics and the C-F cleavage reaction for R2, the difluoroacetate lead (Methylobacterium marchantiae, GenBank GJE15807), over the one-microsecond trajectory. (a) Reactive-pose trajectory. (b) Near-attack-conformation dashboard. (c) Per-residue occupancy and NAC-criteria summary. (d) Binding (MM-GBSA) against reactivity, per frame. (e) WaterMap hydration sites. (f) Free-energy landscapes of the reactive coordinate. (g) Prepared-pose comparative Ramachandran plot. (h) QM/MM relaxed-scan reaction profile over the scissile C-F bond (per-frame and ensemble barrier, reaction energy, departing-fluoride charge and the geometric cleavage verdict).",
    "MD_Rank8_Dynamics_and_Reaction": "Near-attack dynamics and C-F cleavage, R8. Per-system near-attack dynamics and the C-F cleavage reaction for R8, the trifluoroacetate lead (Methylocella tundrae, UniProt A0A4V6IMR0), over the one-microsecond trajectory. (a) Reactive-pose trajectory. (b) Near-attack-conformation dashboard. (c) Per-residue occupancy and NAC-criteria summary. (d) Binding (MM-GBSA) against reactivity, per frame. (e) WaterMap hydration sites. (f) Free-energy landscapes of the reactive coordinate. (g) Prepared-pose comparative Ramachandran plot. (h) QM/MM relaxed-scan reaction profile over the scissile C-F bond (per-frame and ensemble barrier, reaction energy, departing-fluoride charge and the geometric cleavage verdict).",
    "MD_Rank69_Dynamics_and_Reaction": "Near-attack dynamics and C-F cleavage, R69. Per-system near-attack dynamics and the C-F cleavage reaction for R69, the RPA1163/3R3U reference on its native substrate fluoroacetate, over the one-microsecond trajectory. (a) Reactive-pose trajectory. (b) Near-attack-conformation dashboard. (c) Per-residue occupancy and NAC-criteria summary. (d) Binding (MM-GBSA) against reactivity, per frame. (e) WaterMap hydration sites. (f) Free-energy landscapes of the reactive coordinate. (g) Prepared-pose comparative Ramachandran plot. (h) QM/MM relaxed-scan reaction profile over the scissile C-F bond (per-frame and ensemble barrier, reaction energy, departing-fluoride charge and the geometric cleavage verdict). Residues are numbered as in the prepared model, which runs two positions higher than the deposited 3R3U numbering used elsewhere (Asp112 here corresponds to the canonical Asp110).",
}


def write_legends(main_figs, supp_figs):
    """Write Figure_Legends.md (main + supplementary), matching the figures just emitted."""
    out = ["# DeFluorX - figure legends", "",
           "_Auto-generated by Step 08; one entry per emitted figure._", ""]
    out.append("## Main figures\n")
    for num, name in main_figs:
        leg = MAIN_LEGENDS.get(name, name.replace("_", " ") + ".")
        out.append(f"**Figure {num}** | {leg}\n")
    out.append("\n## Supplementary figures\n")
    for num, name in supp_figs:
        leg = SUPP_LEGENDS.get(name, name.replace("_", " ") + ".")
        out.append(f"**Supplementary Figure {num}** | {leg}\n")
    (OUT_DIR / "Figure_Legends.md").write_text("\n".join(out), encoding="utf-8")
    _say(f"  ok  Figure_Legends.md ({len(main_figs)} main + {len(supp_figs)} supplementary)")


def main() -> int:
    _banner("08_Figures_and_Deck_DeFluorX.py", "Figure & deck assembly (merged composites + defluorination cascade)")
    _reset_outputs()
    used: set[Path] = set()
    main_figs: list = []        # (number, name) for Figure_Legends.md
    supp_figs: list = []

    # Build the cross-enzyme defluorination cascade (best per-class TFA/DFA/MFA degrader, merged from the
    # former Step-07 reaction-mechanism scene). It is drawn to its own SVG first, then placed as panel a
    # of the mechanism / QM-MM-validation main figure (panel b = the all-jobs PES overlay).
    _casc = None
    try:
        _casc = build_cascade(_M, OUT_DIR / "_cascade")
        _say(f"  ok  cascade source built: {_casc.name}" if _casc else "  --  cascade: no cleaving per-class enzymes found")
    except Exception as _ce:
        _say(f"  !!  cascade figure skipped ({_ce})")

    _say(f"Run: {RUN.name}  (format: {FIG_FORMAT})")
    _head("MAIN manuscript figures:")
    # Figure 01 is always reserved for the mechanism-guided workflow / graphical overview (the user keeps
    # and places that panel themselves); the data composites below start at Figure 02 whether or not the
    # curated HTML renders here (headless browser may be absent on this machine).
    fi = 1
    if _render_workflow(WORKFLOW_HTML, OUT_DIR / f"Figure_01_Workflow_Overview.{FIG_FORMAT}"):
        main_figs.append((1, "Workflow_Overview"))
        _say("  ok  Figure 01 = workflow overview (user-supplied panel kept and counted as Figure 1)")
    else:
        _say("  ::  Figure 01 reserved for the user-supplied workflow overview "
              f"(drop Figure_01_Workflow_Overview.{FIG_FORMAT} into the output folder); "
              "data figures start at Figure 02 by design")
    # The screen-funnel Sankey belongs in Figure 1, beside the workflow. Figure 1 is a curated HTML
    # export (its own sub-panels; a foreignObject that will not rasterise), so the Sankey is composited
    # in by hand rather than montaged here. Reserve it so no supplementary glob reclaims it.
    used.add(_key(SANKEY_PANEL))
    if _resolved(SANKEY_PANEL).exists():
        _say(f"  ..  Sankey reserved for Figure 1 (place beside the workflow by hand): {_resolved(SANKEY_PANEL).name}")
    else:
        _say("  --  Sankey panel pending (10_Sankey_Workflow)")
    for name, panels in MAIN_FIGS:
        present = [p for p in panels if p is _ROWBREAK or _resolved(p).exists()]
        for p in panels:
            if p is not _ROWBREAK:
                used.add(_key(p))                       # reserve pending panels too, so supp never grabs them
        missing = [p.name for p in panels if p is not _ROWBREAK and not _resolved(p).exists()]
        if not present:
            _say(f"  --  (all panels pending) {name}"); continue
        fi += 1
        montage(present, OUT_DIR / f"Figure_{fi:02d}_{name}.{FIG_FORMAT}")
        main_figs.append((fi, name))
        if missing:
            _say(f"        (pending, fold in on re-run: {', '.join(missing)})")

    _head("SUPPLEMENTARY figures:")
    si = 0
    for name, sources in SUPP_GROUPS:
        pngs = []
        for base, pat in sources:
            if not base.exists():
                continue
            # gather both formats so an SVG-only panel (PNG deleted after the SSOT switch) is still found
            cand = sorted(base.glob(pat)) + sorted(base.glob(pat.replace(".png", ".svg")))
            for p in cand:
                k = _key(p)
                if k in used or (_EXCLUDE_PARTS & set(p.parts)):
                    continue
                if any(d in p.name for d in _DROP_PANELS):     # extracted for the size-preference paper
                    continue
                used.add(k); pngs.append(_substitute_panel(_resolved(p)))   # vector sibling, cropped copy if any
        if not pngs:
            continue
        pngs = _apply_supp_layout(name, pngs)   # optional per-figure row art-direction
        si += 1
        montage(pngs, OUT_DIR / f"Supplementary_Figure_{si:02d}_{name}.{FIG_FORMAT}")
        supp_figs.append((si, name))

    # Per-rank deep-dive: ONE A4 figure per MD rank (reactive pose + near-attack dynamics + reaction).
    # The old Part 1 MD-stability figures are gone - RMSD/RMSF are Figure 6b, the ligand-protein 2D
    # interaction map is Figure 6f-i, so only the reactive-pose panel survived and now folds in here.
    for rp in _md_ranks():
        n = rp.name.split("_", 1)[1]
        panels = _rank_panels(rp, n)
        for p in panels:
            if p is not _ROWBREAK:
                used.add(_key(p))
        present = [p for p in panels if p is _ROWBREAK or _resolved(p).exists()]
        if not any(p is not _ROWBREAK for p in present):
            continue
        si += 1
        montage(present, OUT_DIR / f"Supplementary_Figure_{si:02d}_MD_Rank{n}_Dynamics_and_Reaction.{FIG_FORMAT}")
        supp_figs.append((si, f"MD_Rank{n}_Dynamics_and_Reaction"))
        missing = [p.name for p in panels if p is not _ROWBREAK and not _resolved(p).exists()]
        if missing:
            _say(f"        (pending QSite, fold in on re-run: {', '.join(missing)})")

    # Coverage audit: list any source panel that no figure placed (excluding deliberate skips).
    orphans, seen_orphan = [], set()
    for base in (_V, _T, _P, _M, _D):
        if not base.exists():
            continue
        for p in list(base.rglob("*.png")) + list(base.rglob("*.svg")):
            k = _key(p)
            if k in used or k in seen_orphan or (_EXCLUDE_PARTS & set(p.parts)):
                continue
            if any(s in str(p) for s in _AUDIT_SKIP) or any(d in p.name for d in _DROP_PANELS):
                continue
            seen_orphan.add(k); orphans.append(p)
    if orphans:
        _say(f"\nCoverage audit: {len(orphans)} UNPLACED panel(s):")
        for p in sorted(orphans):
            _say(f"    {p.relative_to(RUN)}")
    else:
        _say("\nCoverage audit: every source panel is placed (or a deliberate skip).")

    write_legends(main_figs, supp_figs)
    # The cascade source folder is only an intermediate for the mechanism figure's panel a; drop it.
    shutil.rmtree(OUT_DIR / "_cascade", ignore_errors=True)

    _rows = [("Main figures", str(fi)), ("Supplementary figures", str(si)),
             ("Legends", "Figure_Legends.md"), ("Output folder", OUT_DIR.name)]
    _kw = max(len(k) for k, _ in _rows); _vw = max(len(v) for _, v in _rows)
    _h = "─"; _bar = lambda: _col("│", _C.CYAN)
    _say("")
    _say(_col("  ┌" + _h * (_kw + 2) + "┬" + _h * (_vw + 2) + "┐", _C.CYAN))
    for _k, _v in _rows:
        _say("  " + _bar() + f" {_k.ljust(_kw)} " + _bar() + f" {_v.ljust(_vw)} " + _bar())
    _say(_col("  └" + _h * (_kw + 2) + "┴" + _h * (_vw + 2) + "┘", _C.CYAN))
    return 0


# ============================================================================================
# EDITABLE POWERPOINT DECK  (merged in from the former svg_to_pptx_DeFluorX.py)
# Reads the merged Figure_*.svg / Supplementary_Figure_*.svg this script writes and rebuilds each as
# ONE slide of fully editable objects: every sub-figure is its own movable/selectable picture, every
# panel letter a real text box, and every vector panel is embedded as native SVG (with a PNG fallback)
# so PowerPoint can 'Convert to Shape' it. The .svg files are only READ, never modified.
# Needs python-pptx + rsvg-convert on PATH; run under base python3 (the PFAS env lacks rsvg-convert).
# Invoke:  python DeFluorX_Figures.py           (build the figures, THEN rebuild the deck fresh - default)
#          python DeFluorX_Figures.py figures   (figures only, skip the deck)
#          python DeFluorX_Figures.py pptx      (deck only, from the SVGs already on disk)
# ============================================================================================
ASVG_NS = "http://schemas.microsoft.com/office/drawing/2016/SVG/main"
OUT_PPTX = OUT_DIR / "DeFluorX_Figures.pptx"           # -> beside the figures, in Figures_HQ/

EMU_IN     = 914400
SLIDE_W_IN = 13.333
SLIDE_H_IN = 20.0
MARGIN_IN  = 0.5
USABLE_W   = SLIDE_W_IN - 2 * MARGIN_IN
USABLE_H   = SLIDE_H_IN - 2 * MARGIN_IN
LABEL_HEX     = "0F0F0F"          # near-black panel letter, matching the SVG figures
SVG_EMBED_MAX = 1_500_000         # embed panels as vector below this; larger -> high-res PNG only

try:
    from pptx import Presentation
    from pptx.util import Emu, Pt
    from pptx.dml.color import RGBColor
    from pptx.enum.text import MSO_ANCHOR
    from pptx.opc.package import Part
    from pptx.opc.constants import RELATIONSHIP_TYPE as RT
    from pptx.oxml.ns import qn
    _HAS_PPTX = True
except ImportError:
    _HAS_PPTX = False

_PPTX_TMP = None


def _ptmp() -> Path:
    global _PPTX_TMP
    if _PPTX_TMP is None:
        _PPTX_TMP = Path(tempfile.mkdtemp(prefix="dfx_pptx_"))
    return _PPTX_TMP


def _emu(inches: float) -> "Emu":
    return Emu(int(round(inches * EMU_IN)))


def _vb_wh(root) -> tuple[float, float]:
    vb = root.get("viewBox")
    if vb:
        p = [float(x) for x in vb.replace(",", " ").split()]
        if len(p) == 4:
            return p[2], p[3]
    def _u(v):
        return float(re.sub(r"[a-z%]+$", "", v)) if v else 1.0
    return _u(root.get("width")), _u(root.get("height"))


def _rasterise(svg_path: Path, width: int) -> "bytes | None":
    out = _ptmp() / (svg_path.stem + f"_{width}.png")
    try:
        subprocess.run([_RSVG_BIN, "-b", "white", "-w", str(width),
                        str(svg_path), "-o", str(out)], check=True, timeout=300)
        return out.read_bytes()
    except Exception as e:
        _say(f"     !! rasterise failed {svg_path.name}: {type(e).__name__}")
        return None


def _add_svg_picture(slide, svg_bytes: bytes, png_bytes: bytes, x, y, cx, cy):
    """Picture whose primary blip is the PNG fallback and whose svgBlip extension is the true vector.
    PowerPoint renders the vector and offers 'Convert to Shape'; anything older shows the PNG."""
    pic = slide.shapes.add_picture(BytesIO(png_bytes), x, y, cx, cy)
    pkg = slide.part.package
    part = Part(pkg.next_partname("/ppt/media/image%d.svg"), "image/svg+xml", pkg, svg_bytes)
    rId = slide.part.relate_to(part, RT.IMAGE)
    blip = pic._element.xpath(".//a:blip")[0]
    extLst = blip.makeelement(qn("a:extLst"), {})
    ext = extLst.makeelement(qn("a:ext"), {"uri": "{96DAC541-7B7A-43D3-8B79-37D633B846F1}"})
    ext.append(ext.makeelement(f"{{{ASVG_NS}}}svgBlip", {qn("r:embed"): rId}))
    extLst.append(ext)
    blip.append(extLst)
    return pic


def _add_label(slide, letter: str, x, y, pt: float):
    tb = slide.shapes.add_textbox(x, y, _emu(0.5), _emu(0.5))
    tf = tb.text_frame
    tf.word_wrap = False
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = MSO_ANCHOR.TOP
    r = tf.paragraphs[0].add_run()
    r.text = letter
    r.font.bold = True
    r.font.size = Pt(max(11.0, pt))
    r.font.name = "Arial"
    r.font.color.rgb = RGBColor.from_string(LABEL_HEX)


def _standalone_panel(el) -> bytes:
    """A nested panel <svg> serialised as its own document (position attrs stripped)."""
    clone = etree.fromstring(etree.tostring(el))
    for a in ("x", "y"):
        clone.attrib.pop(a, None)
    return etree.tostring(clone, xml_declaration=True, encoding="utf-8")


def _decode_data_uri(href: str) -> "bytes | None":
    m = re.match(r"data:[^;]+;base64,(.*)$", href, re.S)
    return base64.b64decode(m.group(1)) if m else None


def _pptx_convert(svg_path: Path, prs) -> int:
    root = etree.parse(str(svg_path)).getroot()
    W, H = _vb_wh(root)
    scale = min(USABLE_W / W, USABLE_H / H)          # inches per SVG unit
    off_x = MARGIN_IN + (USABLE_W - W * scale) / 2
    off_y = MARGIN_IN + (USABLE_H - H * scale) / 2

    def X(px): return _emu(off_x + px * scale)
    def Y(px): return _emu(off_y + px * scale)
    def L(px): return _emu(px * scale)

    slide = prs.slides.add_slide(prs.slide_layouts[6])
    panels = 0

    children = [c for c in root if isinstance(c.tag, str)]
    vec_or_img = [c for c in children if etree.QName(c).localname in ("svg", "image")]

    if not vec_or_img:
        # Figure 01 (HTML foreignObject workflow): no separable panels -> whole graphic as one raster.
        png = _rasterise(svg_path, 2600)
        if png:
            slide.shapes.add_picture(BytesIO(png), X(0), Y(0), L(W), L(H))
            panels = 1
        else:
            _say(f"     -- {svg_path.name}: no rasterisable content, empty slide")
        return panels

    for c in children:
        ln = etree.QName(c).localname
        if ln == "rect":
            continue
        if ln == "text":
            # Panel letters are written as two stacked nodes (a white halo with fill="none" behind the
            # dark glyph) so they survive the MS Word SVG engine; the deck needs only one editable text
            # box per letter, so skip the halo node and keep the glyph.
            if (c.get("fill") or "").strip().lower() == "none":
                continue
            x = float(c.get("x", 0)); y = float(c.get("y", 0))
            fs = float(c.get("font-size", 40))
            _add_label(slide, (c.text or "").strip(),
                       X(x - fs * 0.35), Y(y - fs), fs * scale * 72)
            continue
        x = float(c.get("x", 0)); y = float(c.get("y", 0))
        w = float(c.get("width", 0)); h = float(c.get("height", 0))
        if ln == "image":
            data = _decode_data_uri(c.get(f"{{{XLINK_NS}}}href") or c.get("href") or "")
            if data:
                slide.shapes.add_picture(BytesIO(data), X(x), Y(y), L(w), L(h))
                panels += 1
            continue
        # nested vector panel
        svg_bytes = _standalone_panel(c)
        tmp = _ptmp() / f"panel_{svg_path.stem}_{panels}.svg"
        tmp.write_bytes(svg_bytes)
        raster_w = int(min(2400, max(800, w * 1.5)))
        png = _rasterise(tmp, raster_w)
        if png is None:
            continue
        if len(svg_bytes) <= SVG_EMBED_MAX:
            _add_svg_picture(slide, svg_bytes, png, X(x), Y(y), L(w), L(h))
        else:
            hi = _rasterise(tmp, int(min(3000, max(1400, w * 2)))) or png
            slide.shapes.add_picture(BytesIO(hi), X(x), Y(y), L(w), L(h))
        panels += 1
    return panels


def build_pptx() -> int:
    """Assemble the editable PowerPoint deck from the Figure_*.svg already written next to this script."""
    if not _HAS_PPTX:
        _say("python-pptx not available - skipping deck (pip install python-pptx, run under base python3).")
        return 1
    figs = sorted(OUT_DIR.glob("Figure_*.svg")) + sorted(OUT_DIR.glob("Supplementary_Figure_*.svg"))
    if not figs:
        _say("No Figure_*.svg found - run the figure build first."); return 1
    prs = Presentation()
    prs.slide_width  = _emu(SLIDE_W_IN)
    prs.slide_height = _emu(SLIDE_H_IN)
    total = 0
    for f in figs:
        n = _pptx_convert(f, prs)
        total += n
        _say(f"  ok  {f.name}  ({n} objects)")
    prs.save(str(OUT_PPTX))
    size_mb = OUT_PPTX.stat().st_size / 1e6
    _say(f"\nDone: {len(figs)} slides, {total} panel objects -> {OUT_PPTX.name}  ({size_mb:.1f} MB)")
    Presentation(str(OUT_PPTX))         # sanity re-open
    _say("Re-opened cleanly.")
    return 0


# =============================================================================
#  DEFLUORINATION-CASCADE FIGURE  (merged from the former Step-07 reaction-mechanism scene)
#  One tall figure: the three best per-class FAcD variants (TFA -> DFA -> MFA -> glycolate),
#  each a condensed hydrolytic double-displacement defluorination drawn from the real QM .mae
#  geometry, with a per-band QM/MM reaction-coordinate plot + QM/MM.geometry.WaterMap.MD strip.
#  Self-contained: nothing from the pipeline is imported or executed; reads only this run's
#  QSite/MD outputs (01_Master_Ranking.csv, Rank_*/{10_QSite_Scan_Data.csv,01_MD_Stats.json,
#  QSite_SN2/*/*.01.mae}). Residue numbering = CFG SSOT (3R3U).
# =============================================================================
import math as _c_math
import csv as _c_csv, json as _c_json
import numpy as _c_np
import matplotlib as _c_mpl
_c_mpl.use("Agg")
import matplotlib.pyplot as _c_plt
from matplotlib.patches import FancyArrowPatch as _FAP, Circle as _Circle, FancyBboxPatch as _FBP

# short local aliases so the ported body reads unchanged
np = _c_np; math = _c_math; plt = _c_plt
FancyArrowPatch = _FAP; Circle = _Circle; FancyBboxPatch = _FBP

C_C, C_O, C_F, C_H = "#3A3A3A", "#C0392B", "#1592A8", "#9AA0A6"
PUR, ORG, GRN, PINK = "#6A2CBF", "#D55E00", "#009E73", "#CC79A7"
BLUE, INK, MUT = "#0072B2", "#1A1A1A", "#5A5A5A"
OXY = "#8C6D4F"
GEOM = "#5B4B8A"
GREY = "#8C8C8C"
_3TO1 = {"ASP":"D","GLU":"E","HIS":"H","ARG":"R","LYS":"K","TRP":"W","TYR":"Y","PHE":"F",
         "SER":"S","THR":"T","CYS":"C","ASN":"N","GLN":"Q","GLY":"G","ALA":"A","VAL":"V",
         "LEU":"L","ILE":"I","MET":"M","PRO":"P"}
def _short(code):
    code = str(code).strip()
    for k, v in _3TO1.items():
        if code.upper().startswith(k):
            return v + code[len(k):]
    return code
def _title(code):
    code = str(code).strip()
    for k in _3TO1:
        if code.upper().startswith(k):
            return k.title() + code[len(k):]
    return code
_RES_DEFAULT = dict(nuc="D110", nuc3="Asp110", base="H280", acid="D134", c1="R111", c2="R114",
                    sH="H155", sW="W156", sY="Y219")
_RES = dict(_RES_DEFAULT)
def load_residues(acc):
    """Per-enzyme catalytic residues (from this run's Ranked CSV Mapped_* columns), short-coded."""
    import csv as _c
    base = ROOT.parent / "1_Boltz2_Production"
    files = sorted(base.glob("*Ranked*.csv"), key=lambda p: p.stat().st_mtime) if base.exists() else []
    if not files:
        return dict(_RES_DEFAULT)
    for row in _c.DictReader(open(files[-1])):
        jn = row.get("job_name", row.get("Job_Name", ""))
        if acc in jn:
            g = lambda k: row.get(k, "")
            return dict(nuc=_short(g("Mapped_Nucleophile")), nuc3=_title(g("Mapped_Nucleophile")),
                        base=_short(g("Mapped_Base")), acid=_short(g("Mapped_Acid")),
                        c1=_short(g("Mapped_Clamp1")), c2=_short(g("Mapped_Clamp2")),
                        sH=_short(g("Mapped_Stabiliser_H")), sW=_short(g("Mapped_Stabiliser_W")),
                        sY=_short(g("Mapped_Stabiliser_Y")))
    return dict(_RES_DEFAULT)
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 16: "S", 17: "Cl"}
COV = {1: 0.31, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57, 16: 1.05, 17: 1.02}

# run root (the 7_QMMM_Defluorination dir) and output dir, set by build_cascade()
ROOT = None
_CASC_OUT = None
CASCADE = []

# chemistry fixed per substrate class: (cascade substrate, cascade product, direct hydrolytic product)
_CLASS_CHEM = {
    "TFA": ("trifluoroacetate",  "difluoroacetate",   "difluoroglycolate"),
    "DFA": ("difluoroacetate",   "monofluoroacetate", "fluoroglycolate"),
    "MFA": ("monofluoroacetate", "glycolate",         "glycolate"),
}
# known accession -> scientific species (preprint systems); fall back to the Job_Name organism mnemonic
SPECIES_MAP = {
    "A0A4V6IMR0": "Methylocella tundrae",
    "GJE15807":   "Methylobacterium marchantiae",
    "A0A2U3PT06": "Bradyrhizobium vignae",
}


# ----------------- geometry parsing (from the real QM .mae) -----------------
def _parse(mp):
    txt = Path(mp).read_text(errors="ignore").splitlines(); cts = []; i = 0; n = len(txt)
    while i < n:
        if txt[i].strip().startswith("m_atom["):
            j = i + 1; sch = []
            while j < n and ":::" not in txt[j]:
                s = txt[j].strip()
                if s and not s.startswith("#"): sch.append(s)
                j += 1
            j += 1; atoms = []
            while j < n and ":::" not in txt[j]:
                if txt[j].strip():
                    tok = re.findall(r'"[^"]*"|\S+', txt[j].strip()); dd = {"idx": tok[0]}
                    for k, nm in enumerate(sch): dd[nm] = tok[k + 1] if k + 1 < len(tok) else None
                    atoms.append(dd)
                j += 1
            cts.append(atoms); i = j
        else: i += 1
    return cts


def _av(a):
    q = lambda v: v.strip().strip('"').strip() if v else v
    return dict(idx=int(a["idx"]), xyz=np.array([float(a["r_m_x_coord"]), float(a["r_m_y_coord"]), float(a["r_m_z_coord"])]),
                resname=q(a.get("s_m_pdb_residue_name")), Z=int(a["i_m_atomic_number"]),
                mull=(None if a.get("r_qs_Mulliken_Charges") in (None, "<>") else float(a["r_qs_Mulliken_Charges"])))


def bonded(a, b):
    return float(np.linalg.norm(a["xyz"] - b["xyz"])) < 1.3 * (COV.get(a["Z"], 0.7) + COV.get(b["Z"], 0.7))


def retained_substituents(r):
    mae = next((ROOT / f"Rank_{r}/QSite_SN2").glob("*/*.01.mae"))
    cts = [[_av(a) for a in c] for c in _parse(mae)]
    cts = [ct for ct in cts if any(a["resname"] == "LIG" and a["Z"] == 9 and a["mull"] is not None for a in ct)] or cts
    ct0, scene = cts[0], cts[-1]
    lig0 = [a for a in ct0 if a["resname"] == "LIG"]; Cs = [a for a in lig0 if a["Z"] == 6]; Fs = [a for a in lig0 if a["Z"] == 9]
    fnd = lambda ct, ref: next(a for a in ct if a["idx"] == ref["idx"])
    carbox = lambda c: sum(1 for a in lig0 if a["Z"] == 8 and bonded(a, c)) >= 2
    Ca = next((c for c in Cs if any(bonded(f, c) for f in Fs) and any(bonded(c, cc) for cc in Cs if cc["idx"] != c["idx"] and carbox(cc))),
              max(Cs, key=lambda c: sum(bonded(f, c) for f in Fs)))
    caL = fnd(scene, Ca)
    Fsc = max(Fs, key=lambda f: float(np.linalg.norm(fnd(scene, f)["xyz"] - caL["xyz"])) - float(np.linalg.norm(f["xyz"] - Ca["xyz"])))
    carbC = next((cc for cc in Cs if cc["idx"] != Ca["idx"] and carbox(cc) and bonded(cc, Ca)), None)
    ret = [a for a in lig0 if a["idx"] not in (Ca["idx"], Fsc["idx"]) and (carbC is None or a["idx"] != carbC["idx"]) and bonded(a, Ca)]
    labels = []
    for a in sorted(ret, key=lambda x: -x["Z"]):
        labels.append("R" if a["Z"] == 6 else SYM.get(a["Z"], "?"))
    return mae, labels


# ----------------- drawing primitives -----------------
def atom(ax, p, txt, col, fs=12, rad=0.12, fc="white", lw=1.6, z=6, tcol=None):
    ax.add_patch(Circle(p, rad, fc=fc, ec=col, lw=lw, zorder=z))
    ax.text(p[0], p[1], txt, color=(tcol or col), fontsize=fs, ha="center", va="center", fontweight="bold", zorder=z + 1)


def bond(ax, a, b, col=C_C, lw=2.4, ls="-", z=3, shrink=0.14):
    d = (b - a); L = np.linalg.norm(d); u = d / (L + 1e-9)
    ax.plot([a[0] + u[0] * shrink, b[0] - u[0] * shrink], [a[1] + u[1] * shrink, b[1] - u[1] * shrink],
            color=col, lw=lw, ls=ls, solid_capstyle="round", zorder=z)


def dbl_bond(ax, a, b, col, lw=2.1, z=3):
    d = (b - a); u = d / (np.linalg.norm(d) + 1e-9); n = np.array([-u[1], u[0]]) * 0.05
    for s in (n, -n):
        ax.plot([a[0] + u[0] * 0.14 + s[0], b[0] - u[0] * 0.14 + s[0]],
                [a[1] + u[1] * 0.14 + s[1], b[1] - u[1] * 0.14 + s[1]], color=col, lw=lw, zorder=z)


def carboxylate(ax, cpos, direc, col=C_O, labelminus=True):
    cc = np.array(cpos)
    o1 = cc + np.array([math.cos(math.radians(direc - 30)), math.sin(math.radians(direc - 30))]) * 0.55
    o2 = cc + np.array([math.cos(math.radians(direc + 30)), math.sin(math.radians(direc + 30))]) * 0.55
    atom(ax, cc, "C", C_C, fs=9, rad=0.11)
    dbl_bond(ax, cc, o1, col); bond(ax, cc, o2, col, lw=2.0)
    atom(ax, o1, "O", col, fs=10); atom(ax, o2, "O", col, fs=10)
    if labelminus:
        ax.text(o2[0] + 0.12, o2[1] + 0.12, "$-$", color=col, fontsize=13, ha="center", va="center", zorder=8)
    return cc, o1, o2


def substrate_core(ax, c, ret, right_group, right_col, right_bond="-", left_O=None, left_bond=None):
    c = np.array(c); L = 0.62
    up = c + np.array([math.cos(math.radians(78)), math.sin(math.radians(78))]) * L
    dl = c + np.array([math.cos(math.radians(212)), math.sin(math.radians(212))]) * L
    cooC = c + np.array([math.cos(math.radians(312)), math.sin(math.radians(312))]) * L
    rightp = c + np.array([L, 0.0]); leftp = c + np.array([-L, 0.0])
    bond(ax, c, cooC, C_C, lw=2.4); _cc, oA, oB = carboxylate(ax, cooC, 312, C_O)
    cols = {"F": C_F, "H": C_H, "R": "#6E6E6E"}
    for sym, pos in zip((ret + ["", ""])[:2], (up, dl)):
        if not sym: continue
        col = cols.get(sym, C_H)
        bond(ax, c, pos, col, lw=(1.7 if sym == "H" else 2.3))
        atom(ax, pos, ("R" if sym == "R" else sym), col, fs=(9 if sym == "H" else 11), rad=(0.11 if sym == "H" else 0.13))
    if right_group is not None:
        bond(ax, c, rightp, right_col, lw=2.3, ls=right_bond)
        atom(ax, rightp, right_group, right_col, fs=12, rad=0.14)
    if left_O is not None:
        bond(ax, c, leftp, left_O[1], lw=2.3, ls=(left_bond or "-"))
        fcl = left_O[2] if len(left_O) > 2 else "white"
        atom(ax, leftp, left_O[0], left_O[1], fs=11, rad=0.15, fc=fcl, tcol=("white" if fcl != "white" else left_O[1]))
    atom(ax, c, r"C$_\alpha$", C_C, fs=10, rad=0.16, lw=1.9)
    return dict(c=c, up=up, dl=dl, cooC=cooC, right=rightp, left=leftp, oA=oA, oB=oB)


# ----------------- data loaders (this run only) -----------------
def geom_attack(r):
    mae = next((ROOT / f"Rank_{r}/QSite_SN2").glob("*/*.01.mae"))
    cts = [[_av(a) for a in c] for c in _parse(mae)]
    cts = [ct for ct in cts if any(a["resname"] == "LIG" and a["Z"] == 9 and a["mull"] is not None for a in ct)] or cts
    ct0, scene = cts[0], cts[-1]
    lig0 = [a for a in ct0 if a["resname"] == "LIG"]; Cs = [a for a in lig0 if a["Z"] == 6]; Fs = [a for a in lig0 if a["Z"] == 9]
    fnd = lambda ct, ref: next(a for a in ct if a["idx"] == ref["idx"])
    carbox = lambda c: sum(1 for a in lig0 if a["Z"] == 8 and bonded(a, c)) >= 2
    Ca = next((c for c in Cs if any(bonded(f, c) for f in Fs) and any(bonded(c, cc) for cc in Cs if cc["idx"] != c["idx"] and carbox(cc))),
              max(Cs, key=lambda c: sum(bonded(f, c) for f in Fs)))
    caL = fnd(scene, Ca)
    Fsc = max(Fs, key=lambda f: float(np.linalg.norm(fnd(scene, f)["xyz"] - caL["xyz"])) - float(np.linalg.norm(f["xyz"] - Ca["xyz"])))
    protO = [a for a in ct0 if a["resname"] != "LIG" and a["Z"] == 8]
    nuc = min(protO, key=lambda o: float(np.linalg.norm(o["xyz"] - Ca["xyz"])))
    dOC = float(np.linalg.norm(nuc["xyz"] - Ca["xyz"])); dCF = float(np.linalg.norm(Fsc["xyz"] - Ca["xyz"]))
    v1 = nuc["xyz"] - Ca["xyz"]; v2 = Fsc["xyz"] - Ca["xyz"]
    ang = float(np.degrees(np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))))
    return dOC, dCF, ang


def _master_rows():
    return list(_c_csv.DictReader(open(ROOT / "01_Master_Ranking.csv")))


def _master_row(acc):
    for row in _master_rows():
        if acc in row.get("Job_Name", ""):
            return row
    return {}


def load_scan(rank):
    rows = list(_c_csv.DictReader(open(ROOT / f"Rank_{rank}/10_QSite_Scan_Data.csv")))
    best = [x for x in rows if x.get("Is_Best_Frame", "").strip().lower() in ("true", "1", "yes")] or rows
    ff = best[0]["Frame_Folder"]
    pts = sorted([x for x in best if x["Frame_Folder"] == ff], key=lambda x: float(x["Scan_Point"]))
    sp = [float(x["Scan_Point"]) for x in pts]
    cf = [float(x["CF_Dist_A"]) for x in pts]
    fq = [float(x["F_Charge_Mulliken"]) if x["F_Charge_Mulliken"] not in ("", "<>") else None for x in pts]
    return sp, cf, fq


def load_qm(rank, acc):
    m = _master_row(acc)
    st = _c_json.load(open(ROOT / f"Rank_{rank}/01_MD_Stats.json"))
    g = lambda d, k: d.get(k) if d.get(k) not in (None, "", "<>") else None
    f = lambda v, fmt="{:.2f}": (fmt.format(float(v)) if v not in (None, "", "<>") else "n/a")
    ncl = g(m, "QSite_NFrames_Cleaved"); ntot = g(m, "QSite_NFrames_Scored") or ncl
    return {
        "cleaved": f"{ncl}/{ntot}" if ncl is not None else "n/a",
        "angle": f(g(m, "SN2_Attack_Angle"), "{:.0f}") + "°",
        "barrier": f(g(m, "QSite_Barrier_Mean_kcal"), "{:.1f}"),
        "derxn": (lambda v: ("+" if v and float(v) >= 0 else "") + f(v, "{:.1f}"))(g(m, "QSite_dErxn_kcal")),
        "gbsa": f(g(m, "MMGBSA_dG_Boltzmann_kcal"), "{:.1f}"),
        "wm_sites": g(st, "WM_N_Sites"), "wm_disp": g(st, "WM_N_Unstable"),
        "wm_dg": f(g(st, "WM_Mean_dG"), "{:.1f}"),
        "dwell": f(g(st, "NAC_Total_ns"), "{:.1f}"),
    }


def _job_class(jn: str):
    """Substrate fluorine class from a Job_Name (None for a control or non-fluoroacetate)."""
    s = jn.lower()
    if "control" in s or "3r3u" in s or "deha4" in s:
        return None
    if "trifluoro" in s or s.endswith("_tfa") or "_tfa" in s:
        return "TFA"
    if "difluoro" in s or "_dfa" in s:
        return "DFA"
    if "fluoroacetate" in s or "_mfa" in s or "monofluoro" in s:
        return "MFA"
    return None


def _acc_of(jn: str):
    parts = jn.split("_")
    return parts[2] if len(parts) > 2 else jn


def _species_of(acc: str, jn: str):
    if acc in SPECIES_MAP:
        return SPECIES_MAP[acc]
    parts = jn.split("_")
    org = parts[3] if len(parts) > 3 else ""
    return org if (org and org[0].isalpha()) else acc


def select_cascade(root):
    """Best per-class degrader (TFA, DFA, MFA) among cleaving enzymes, ordered by descending fluorine
    count so the figure reads TFA -> DFA -> MFA -> glycolate. Iterates the actual Rank_N QSite folders
    (the folder number is the authoritative directory, NOT any ranking column), reads each folder's
    Job_Name from its 01_MD_Stats.json, and keeps the best per class by Defluor_Rank (lower = better)."""
    global ROOT
    ROOT = root
    master = {_acc_of(r.get("Job_Name", "")): r for r in _master_rows()}
    best = {}   # class -> (defluor_rank, dirnum, acc, job_name)
    for d in sorted(root.glob("Rank_*")):
        if not (d / "QSite_SN2").exists():
            continue
        try:
            dirnum = int(d.name.split("_")[1])
        except Exception:
            continue
        stj = d / "01_MD_Stats.json"
        jn = _c_json.load(open(stj)).get("Job_Name", "") if stj.exists() else ""
        cls = _job_class(jn)
        if cls is None:
            continue
        acc = _acc_of(jn); mr = master.get(acc, {})
        if str(mr.get("Is_Defluorinating", "")).strip() not in ("1", "1.0", "True", "true", "yes"):
            continue
        try:
            dr = int(float(mr.get("Defluor_Rank", "")))
        except Exception:
            dr = dirnum
        if cls not in best or dr < best[cls][0]:
            best[cls] = (dr, dirnum, acc, jn)
    out = []
    for cls in ("TFA", "DFA", "MFA"):
        if cls not in best:
            continue
        _dr, dirnum, acc, jn = best[cls]
        sub_cas, prod_cas, prod_direct = _CLASS_CHEM[cls]
        out.append((dirnum, sub_cas, prod_cas, prod_direct, acc, _species_of(acc, jn)))
    return out


# ----------------- scene helpers -----------------
def _c_disc(ax, p, col, label, rad=0.33, fs=10.2):
    ax.add_patch(Circle(p, rad, fc=col, ec="white", lw=1.7, zorder=6))
    ax.text(p[0], p[1], label, color="white", fontsize=fs, ha="center", va="center", fontweight="bold", zorder=7)


def _c_hbond(ax, a, b, col, lw=1.4):
    d = (b - a); L = np.linalg.norm(d); u = d / (L + 1e-9)
    ax.plot([a[0] + u[0] * 0.30, b[0] - u[0] * 0.28], [a[1] + u[1] * 0.30, b[1] - u[1] * 0.28],
            ls=(0, (1.5, 1.5)), lw=lw, color=col, zorder=2)


def _c_numbox(ax, p, n, rad=0.25, fs=12.0):
    ax.add_patch(Circle(p, rad, fc=INK, ec="white", lw=1.9, zorder=11))
    ax.text(p[0], p[1], str(n), color="white", fontsize=fs, ha="center", va="center", fontweight="bold", zorder=12)


def _c_earrow(ax, a, b, col=INK, rad=-0.32, lw=1.9, ms=14):
    ax.add_patch(FancyArrowPatch(a, b, connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
                                 mutation_scale=ms, lw=lw, color=col, zorder=9))


def _c_clamp(ax, gcore):
    o1, o2 = gcore["oA"], gcore["oB"]
    r1 = o2 + np.array([0.30, -0.70]); r2 = o1 + np.array([-0.48, -0.70])
    for rp, op, lab in [(r1, o2, _RES["c1"]), (r2, o1, _RES["c2"])]:
        _c_hbond(ax, rp, op, GRN); _c_disc(ax, rp, GRN, lab)


def _c_dyad(ax, cx, ytop, active=False, water=None, target=None):
    h = np.array([cx + 0.04, ytop]); d = np.array([cx - 0.88, ytop - 0.05])
    _c_hbond(ax, h, d, ORG); _c_disc(ax, d, ORG, _RES["acid"]); _c_disc(ax, h, ORG, _RES["base"])
    if active and water is not None and target is not None:
        atom(ax, water, "OH$_2$", BLUE, fs=8.5, rad=0.19, fc="#EAF3FB")
        _c_earrow(ax, h + np.array([0.16, -0.06]), water + np.array([-0.2, 0.07]), col=ORG, rad=0.25, lw=1.5, ms=11)
        ax.add_patch(FancyArrowPatch(water + np.array([0.0, -0.23]), target + np.array([0.0, 0.26]),
                                     arrowstyle="-|>", mutation_scale=12, lw=1.6, color=BLUE, ls=(0, (3, 2)), zorder=4))
    return h


def _c_nuc_node(ax, gcore):
    O = gcore["left"] + np.array([-0.52, 0.0])
    atom(ax, O, "O", PUR, fs=12, rad=0.16, fc=PUR, tcol="white")
    ax.text(O[0] - 0.02, O[1] + 0.36, "$-$", color=PUR, fontsize=16, ha="center", zorder=10)
    ax.text(O[0] - 0.05, O[1] - 0.46, _RES["nuc3"], color=PUR, fontsize=11.0, ha="center", va="top", fontweight="bold", zorder=10)
    return O


def _c_rxn_arrow(ax, x0, x1, y, top=None, bot=None):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>", mutation_scale=18, lw=2.4, color=INK, zorder=1))
    if top:
        ax.text((x0 + x1) / 2, y + 0.36, top, color=INK, fontsize=12.5, ha="center", va="bottom", fontweight="bold", zorder=10)
    if bot:
        ax.text((x0 + x1) / 2, y - 0.24, bot, color=INK, fontsize=10.5, ha="center", va="top", fontstyle="italic", zorder=10)


def _c_scan_inset(fig, ax, rank, dx0, dx1, dy0, dy1):
    X0, X1 = ax.get_xlim(); Y0, Y1 = ax.get_ylim()
    pos = ax.get_position()
    fx = lambda x: pos.x0 + (x - X0) / (X1 - X0) * pos.width
    fy = lambda y: pos.y0 + (y - Y0) / (Y1 - Y0) * pos.height
    rect = [fx(dx0), fy(dy0), fx(dx1) - fx(dx0), fy(dy1) - fy(dy0)]
    sp, cf, fq = load_scan(rank)
    ia = fig.add_axes(rect, zorder=9)
    ia.plot(sp, cf, "-o", color=BLUE, lw=1.8, ms=3.0)
    ia.axhline(1.8, color=GRN, ls=(0, (4, 3)), lw=1.1)
    ia.set_ylabel("C$\\alpha$-F (Å)", color=BLUE, fontsize=9.6); ia.tick_params(labelsize=8.4)
    ia.tick_params(axis="y", colors=BLUE)
    ia.set_xticks([])
    ib = ia.twinx()
    ib.plot(sp, [q if q is not None else np.nan for q in fq], "-s", color="#8A6D1F", lw=1.6, ms=2.8)
    ib.set_ylabel("F charge (e)", color="#8A6D1F", fontsize=9.6); ib.tick_params(labelsize=8.4)
    ib.tick_params(axis="y", colors="#8A6D1F")
    ia.text(0.5, -0.09, "QM scan: reactant → product", transform=ia.transAxes, fontsize=9.0, color=INK, ha="center", va="top")
    for s in ia.spines.values(): s.set_edgecolor("#BBBBBB")
    ia.set_facecolor("white")
    return ia


# --------------------------- one cascade band ---------------------------
def _c_band(ax, fig, yc, step, rank, ret, sub_cas, prod_cas, prod_direct, acc, species, geom):
    global _RES
    _RES = load_residues(acc)
    dOC, dCF, angA = geom
    cap_t, cap_s = yc - 2.42, yc - 2.70
    ax.add_patch(Circle((XSTEP, yc + 0.05), 0.44, fc=INK, ec="none", zorder=5))
    ax.text(XSTEP, yc + 0.05, str(step), color="white", fontsize=20, ha="center", va="center", fontweight="bold", zorder=6)
    ax.text(XSPEC + 0.30, yc + 0.05, acc, color=INK, fontsize=13.0, rotation=90, ha="center", va="center", fontweight="bold", zorder=10)
    ax.text(XSPEC - 0.22, yc + 0.05, species, color=GEOM, fontsize=11.0, rotation=90, ha="center", va="center", fontstyle="italic", fontweight="bold", zorder=10)

    g1 = substrate_core(ax, (XA, yc), ret, right_group="F", right_col=C_F)
    O1 = _c_nuc_node(ax, g1); _c_clamp(ax, g1); _c_dyad(ax, XA, yc + DYH)
    _c_earrow(ax, O1 + np.array([0.22, 0.08]), g1["c"] + np.array([-0.30, 0.02]), rad=-0.30)
    mid1 = (g1["c"] + g1["right"]) / 2
    _c_earrow(ax, mid1 + np.array([0.04, 0.12]), g1["right"] + np.array([0.36, 0.36]), rad=-0.34)
    _c_numbox(ax, (O1[0] - 0.02, yc + 0.48), 1)
    ax.add_patch(FancyBboxPatch((XA - 1.80, yc + 0.82), 1.82, 1.02, boxstyle="round,pad=0.02,rounding_size=0.10",
                                fc="white", ec=PUR, lw=1.8, zorder=8))
    ax.text(XA - 0.89, yc + 1.58, f"O$\\delta$2$\\cdots$C$_\\alpha$  {dOC:.2f} Å", color=PUR, fontsize=11.0, ha="center", va="center", fontweight="bold", zorder=10)
    ax.text(XA - 0.89, yc + 1.28, f"$\\angle$ O–C$_\\alpha$–F  {angA:.0f}°", color=PUR, fontsize=11.0, ha="center", va="center", fontweight="bold", zorder=10)
    ax.text(XA - 0.89, yc + 0.98, f"C$_\\alpha\\cdots$F  {dCF:.2f} Å", color=C_F, fontsize=11.0, ha="center", va="center", fontweight="bold", zorder=10)
    ax.text(XA, cap_t, "Michaelis complex", color=INK, fontsize=13.0, ha="center", fontweight="bold", zorder=10)
    ax.text(XA, cap_s, f"{sub_cas} bound", color=PUR, fontsize=11.0, ha="center", zorder=10)

    g2 = substrate_core(ax, (XB, yc), ret, right_group=None, right_col=C_F, left_O=("O", PUR, PUR))
    _c_clamp(ax, g2)
    EstC = g2["left"] + np.array([-0.62, 0.0]); Od1 = EstC + np.array([-0.42, 0.28])
    atom(ax, EstC, "C", C_C, fs=9.5, rad=0.12); dbl_bond(ax, EstC, Od1, PUR); atom(ax, Od1, "O", PUR, fs=9.5)
    ax.text(EstC[0] - 0.12, yc - 0.52, f'{_RES["nuc3"]} ester', color=PUR, fontsize=11.0, ha="center", va="top", fontweight="bold", zorder=10)
    _c_numbox(ax, ((g2["left"][0] + g2["c"][0]) / 2, yc + 0.34), 3)
    Wat = EstC + np.array([0.15, 1.02])
    h4 = _c_dyad(ax, XB - 0.6, yc + DYH, active=True, water=Wat, target=EstC)
    ax.text(Wat[0] + 0.32, Wat[1], "catalytic\nwater", color=BLUE, fontsize=10.0, ha="left", va="center", zorder=10)
    _c_numbox(ax, ((h4[0] + Wat[0]) / 2 + 0.05, (h4[1] + Wat[1]) / 2), 4)
    oxp = Od1 + np.array([-0.20, 0.86]); _c_disc(ax, oxp, OXY, "Phe", rad=0.31, fs=9.5); _c_hbond(ax, oxp, Od1, OXY)
    ax.text(oxp[0] - 0.62, oxp[1], "oxyanion\nhole", color=OXY, fontsize=9.5, ha="center", va="center", fontweight="bold", zorder=10)
    Fm = g2["c"] + np.array([1.28, 1.24])
    atom(ax, Fm, "F$^-$", C_F, fs=13.5, rad=0.29, fc="#E7F5F8", lw=1.8)
    for aa, lab in [(150, _RES["sH"]), (90, _RES["sW"]), (33, _RES["sY"])]:
        p = Fm + np.array([math.cos(math.radians(aa)), math.sin(math.radians(aa))]) * 1.00
        _c_hbond(ax, p, Fm, PINK); _c_disc(ax, p, PINK, lab)
    _c_numbox(ax, (Fm[0] - 0.52, Fm[1] - 0.50), 2)
    ax.text(XB + 0.2, cap_t, "glycolyl-enzyme ester + F$^-$", color=INK, fontsize=13.0, ha="center", fontweight="bold", zorder=10)
    ax.text(XB + 0.2, cap_s, "water engages ester carbonyl", color=ORG, fontsize=11.0, ha="center", zorder=10)

    g3 = substrate_core(ax, (XC, yc), ret, right_group=None, right_col=C_F, left_O=("OH", C_O, "white"))
    _c_clamp(ax, g3); _c_dyad(ax, XC, yc + DYH)
    _c_numbox(ax, (g3["left"][0] - 0.02, yc + 0.60), 5)
    ax.text(XC, cap_t, f"{prod_direct}", color=INK, fontsize=13.0, ha="center", fontweight="bold", zorder=10)
    ax.text(XC, cap_s, f'{_RES["nuc3"]} regenerated', color=PUR, fontsize=11.0, ha="center", zorder=10)
    _hf = {"difluoroglycolate": "then spontaneous −HF → oxalate",
           "fluoroglycolate": "then spontaneous −HF → glyoxylate"}.get(prod_direct)
    if _hf:
        ax.text(XROLE - 0.60, yc + 0.05, _hf, color=MUT, fontsize=9.0, rotation=270,
                ha="center", va="center", fontstyle="italic", zorder=10)

    _c_rxn_arrow(ax, XA + 1.15, XB - 1.35, yc, top="$-$ F$^-$")
    _c_rxn_arrow(ax, XB + 1.80, XC - 1.3, yc, top="$+$ H$_2$O", bot="hydrolysis")

    _c_scan_inset(fig, ax, rank, 12.10, 15.65, yc + 0.55, yc + 2.60)
    qm = load_qm(rank, acc)
    sx, sy, dyl = 12.15, yc - 0.04, 0.325
    strip = [("C–F cleaved", qm["cleaved"] + " frames", GRN),
             ("S$_N$2 attack (QM)", qm["angle"], INK),
             ("ΔE$^\\ddagger$ (1-D scan, mean)", qm["barrier"] + " kcal/mol", INK),
             ("ΔE$_{rxn}$ (scan)", qm["derxn"] + " kcal/mol", INK),
             ("MM-GBSA ΔG$_{bind}$", qm["gbsa"] + " kcal/mol", PUR),
             ("WaterMap sites", f"{qm['wm_disp']}/{qm['wm_sites']} displ.", PINK),
             ("mean site ΔG", qm["wm_dg"] + " kcal/mol", PINK),
             ("reactive NAC dwell", qm["dwell"] + " ns", BLUE)]
    for i, (lab, val, col) in enumerate(strip):
        yy = sy - i * dyl
        ax.text(sx, yy, lab, color=INK, fontsize=10.2, ha="left", va="center", zorder=10)
        ax.text(sx + 2.55, yy, val, color=col, fontsize=10.6, ha="left", va="center", fontweight="bold", zorder=10)

    ax.text(XROLE, yc + 0.05, f"{sub_cas}  →  {prod_direct}", color=GEOM, fontsize=13.0, rotation=270, ha="center", va="center", fontstyle="italic", fontweight="bold", zorder=10)


def _c_legend_arrow(ax, x, y, style, col, label):
    ls = {"solid": "-", "dash": (0, (3, 2)), "curl": "-"}[style]
    if style == "curl":
        ax.add_patch(FancyArrowPatch((x, y - 0.05), (x + 0.42, y + 0.05), connectionstyle="arc3,rad=-0.5",
                                     arrowstyle="-|>", mutation_scale=10, lw=1.7, color=col, zorder=3))
    else:
        ax.add_patch(FancyArrowPatch((x, y), (x + 0.46, y), arrowstyle="-|>", mutation_scale=11,
                                     lw=(2.1 if style == "solid" else 1.5), color=col, ls=ls, zorder=3))
    ax.text(x + 0.60, y, label, color=INK, fontsize=11.0, ha="left", va="center", zorder=10)


# -------- layout --------
XA, XB, XC = 2.6, 6.4, 10.2
DYH = 2.4
XSTEP = 0.35
XSPEC = -0.85
XROLE = 16.75


def _cascade_build():
    fig = plt.figure(figsize=(16.4, 15.9)); fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.0, 0.006, 1.0, 0.988]); ax.axis("off"); ax.set_aspect("equal")
    ax.set_xlim(-1.4, 17.3); ax.set_ylim(-3.72, 14.2)
    n = len(CASCADE)
    YC = [11.5, 5.9, 0.3] if n == 3 else list(np.linspace(11.5, 0.3, max(n, 1)))
    CXC = (-1.4 + 17.3) / 2

    for i, (rank, sub_cas, prod_cas, prod_direct, acc, species) in enumerate(CASCADE):
        _, ret = retained_substituents(rank)
        _c_band(ax, fig, YC[i], i + 1, rank, ret, sub_cas, prod_cas, prod_direct, acc, species, geom_attack(rank))


    # Flow-wrap legend: all items (role chips, arrow key, step key) sit tightly next to each other and
    # wrap to as many rows as fit the width - with this content that settles to two closely-spaced rows.
    FS = 11.0
    GAP = 0.34
    _tlen = lambda s: len(s.replace("$", "").replace("^", "").replace("_", ""))   # rendered length (drop math)
    chips = [(PUR, "nucleophile (Asp)"), (ORG, "general base (His / Asp)"),
             (PINK, "F$^-$ cradle (His / Trp / Tyr)"), (GRN, "COO$^-$ clamp (Arg / Arg)"),
             (OXY, "oxyanion hole (Phe / Arg N-H)")]
    arrs = [("solid", INK, "reaction step"), ("curl", INK, "electron shift"), ("dash", ORG, "proton / water")]
    steps = ["steps:", "1 attack", "2 C-F cleaved", "3 ester", "4 water activated", "5 hydrolysis"]
    items = ([("chip", c, l) for c, l in chips] + [("arr", s, c, l) for s, c, l in arrs]
             + [("txt", t) for t in steps])
    def _iw(it):                                        # item width incl. its trailing gap
        if it[0] == "chip": return 0.34 + _tlen(it[2]) * 0.116 + GAP
        if it[0] == "arr":  return 0.62 + _tlen(it[3]) * 0.116 + GAP
        return _tlen(it[1]) * 0.106 + GAP
    Wmax = 17.4                                         # wrap width (fits inside the x-range with margins)
    rows, cw = [[]], 0.0
    for it in items:
        w = _iw(it)
        if cw + w > Wmax and rows[-1]:
            rows.append([]); cw = 0.0
        rows[-1].append((it, w)); cw += w
    dyr, yTop = 0.48, -2.82
    row_ws = [sum(w for _, w in r) - GAP for r in rows]
    halfw = max(row_ws) / 2 + 0.40
    yBot = yTop - (len(rows) - 1) * dyr
    ax.add_patch(FancyBboxPatch((CXC - halfw, yBot - 0.23), 2 * halfw, (yTop + 0.23) - (yBot - 0.23),
                                boxstyle="round,pad=0.02,rounding_size=0.10", fc="none", ec="#B4ABCB", lw=1.1, zorder=0))
    for ri, row in enumerate(rows):
        y = yTop - ri * dyr
        x = CXC - row_ws[ri] / 2
        for it, w in row:
            if it[0] == "chip":
                _, col, lab = it
                ax.add_patch(Circle((x + 0.14, y), 0.15, fc=col, ec="white", lw=1.0, zorder=4))
                ax.text(x + 0.36, y, lab, color=INK, fontsize=FS, ha="left", va="center", zorder=10)
            elif it[0] == "arr":
                _, style, col, lab = it
                _c_legend_arrow(ax, x, y, style, col, lab)
            else:
                ax.text(x, y, it[1], color=INK, fontsize=FS, ha="left", va="center", zorder=10)
            x += w

    _CASC_OUT.mkdir(parents=True, exist_ok=True)
    out = _CASC_OUT / "Reaction_Cascade"
    fig.savefig(f"{out}.png", dpi=175, facecolor="white")
    fig.savefig(f"{out}.svg", facecolor="white")
    plt.close(fig)
    return Path(f"{out}.svg")


def build_cascade(root: "Path", out_dir: "Path"):
    """Entry point: select best-per-class cascade from this run and render the figure (svg+png)."""
    global ROOT, _CASC_OUT, CASCADE
    ROOT = Path(root)
    _CASC_OUT = Path(out_dir)
    CASCADE = select_cascade(ROOT)
    if not CASCADE:
        return None
    return _cascade_build()


if __name__ == "__main__":
    import sys
    argv = sys.argv[1:]                   # may carry a <Run> id/path and/or a mode token (pptx|figures)
    if "pptx" in argv:                   # deck only, from the SVGs already on disk
        raise SystemExit(build_pptx())
    _main_t0 = time.perf_counter()
    rc = main()                          # build the merged figures
    '''
    Default run ALSO rebuilds the editable deck fresh (pass "figures" to skip it). Needs python-pptx
    and rsvg-convert; if unavailable, build_pptx warns and the figure build still succeeds.
    '''
    if rc == 0 and "figures" not in argv:
        try:
            build_pptx()
        except Exception as e:
            _say(f"  !! pptx build skipped: {type(e).__name__}: {e}")
    if rc == 0:
        _done = _col("✔  08_Figures_and_Deck_DeFluorX.py  -  Figures & deck complete", _C.OKGREEN + _C.BOLD)
        _say("\n" + _SEP)
        _say("  " + _done + "  " + _PIPE + f"  Total Elapsed: {time.perf_counter() - _main_t0:.0f}s")
        _say(_SEP)
    raise SystemExit(rc)
