#!/usr/bin/env python3
"""
===============================================================================
FAcDs Pipeline  |  Step 06  |  Desmond SID + Prime MM-GBSA Post-Processing
===============================================================================
Runs Schrödinger Event Analysis (event_analysis.py) and Simulation Interaction
Diagram analysis (analyze_simulation.py) on completed Desmond molecular
dynamics trajectories, producing the *_SID-out.eaf files that are consumed
downstream by 07_MD_QMMM_Defluorination_FAcDs.py.

It then runs Prime MM-GBSA (thermal_mmgbsa.py <job>-out.cms) on every completed
MD job — the end-state ligand binding free energy over the MD ensemble — and
plots per-job + combined ΔG_bind. MM-GBSA is complementary to the QSite QM/MM
reaction barrier (Step 07): it scores BINDING, not C–F bond cleavage.

Uses the central CFG / ProjectUtils modules for logging, console styling, and
conventions shared across the pipeline. Runs under the project 'PFAS' conda
environment and shells out to $SCHRODINGER/run for the Schrödinger interpreter
— it does NOT need to be launched with $SCHRODINGER/run.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 15 July 2026
===============================================================================
Usage (standalone):
  python 06_SID_Prime-MMGBSA_FAcDs.py [Boltz-2_Run_Directory]

Usage (from pipeline runner — oomd already managed by 00_00_run_pipeline):
  python 06_SID_Prime-MMGBSA_FAcDs.py [Boltz-2_Run_Directory] --pipeline-mode

-------------------------------------------------------------------------------
Dependency Map
-------------------------------------------------------------------------------
  Script        : 06_SID_Prime-MMGBSA_FAcDs.py
  Role          : Step 06 — Desmond post-simulation post-processing.
                  Produces EAF interaction files for the Step 07 MD engine.
  Imports from  : 00_01_Project_Config_FAcDs.py  (CFG — project metadata)
                  00_02_Project_Utils_FAcDs.py   (ConsoleColours, logging, banners)
  Reads         : Boltz-2_Run_X/6_Physics_Validation/MolecularDynamics/desmond_md_job_R_N/*-out.cms
                  Boltz-2_Run_X/6_Physics_Validation/MolecularDynamics/desmond_md_job_R_N/*_trj
  Writes        : .../desmond_md_job_R_N/*_SID-in.eaf, *_SID-out.eaf, *.log
                  .../desmond_md_job_R_N/*_mmgbsa-prime-out.csv  (per-frame ΔG_bind,
                      carrying a `Frame` column: the trajectory frame each scored
                      structure came from — Step 07 joins on it)
                  .../desmond_md_job_R_N/_MMGBSA_Shards/  (per-shard logs + CSVs; kept
                      so an interrupted MM-GBSA resumes instead of restarting)
                  .../MolecularDynamics/Prime_MMGBSA/00_MMGBSA_Summary.csv
                  .../MolecularDynamics/Prime_MMGBSA/01_MMGBSA_Combined_AllRanks.png
                  .../MolecularDynamics/Prime_MMGBSA/Rank_NN_MMGBSA_Profile_*.png
  Upstream      : Desmond molecular dynamics simulations (manual Maestro step).
  Downstream    : 07_MD_QMMM_Defluorination_FAcDs.py (Step 07; consumes EAF output).
-------------------------------------------------------------------------------
The Critic's Corner: Known Limitations & Failure Points
-------------------------------------------------------------------------------
  1. Sudo (OPTIONAL): Masking systemd-oomd needs root. Standalone, the script
     offers to prime sudo; if sudo is unavailable or skipped it continues in
     NORMAL mode (oomd not masked). Under the pipeline runner (--pipeline-mode)
     oomd masking is owned by the runner across Steps 07 and 08 (also optional).
  2. Trajectory Volume: Processing 100,000 frames sequentially is time-intensive
     (typically 10–15 hours). Do not run multiple instances concurrently.
  3. Local I/O Capping: Relies on -LOCAL to prevent huge tmp partition writes;
     requires sufficient disk space in the destination filesystem.
  4. No Desmond Jobs: If no completed Desmond jobs are found, the script exits
     with status 0 (non-fatal). The pipeline runner continues to Step 07,
     which will then report that no EAF files are present.
  5. MM-GBSA Cost: Prime minimises every scored structure, and the cost is linear
     in their number — CFG.MMGBSA_STEP_SIZE is the only knob that changes the
     wall-clock by an order of magnitude. Step 06 runs the trajectory as concurrent
     frame-range shards, so the read is parallel and overlaps Prime; concurrency is
     capped by free RAM (Prime dominates it), not by core count alone.
  6. MM-GBSA Validity: GB implicit solvent overstabilises anionic PFAS, so ΔG_bind is
     a RELATIVE ranking only; it scores binding, not the QSite reaction barrier. A
     fraction of a percent of frames are failed minimisations (ΔG of hundreds of
     kcal/mol); they are flagged (N_Failed_Minimisations) and excluded from the
     figure's scale, and the MEDIAN is reported because they drag the mean.
===============================================================================
-------------------------------------------------------------------------------
Scientific References:
    1. Protein refinement engine (Schrödinger Prime):
       - Jacobson, M.P. et al. (2004) A hierarchical approach to all-atom protein
         loop prediction. Proteins 55:351–367. DOI: https://doi.org/10.1002/prot.10613
    2. MM-GBSA end-state binding free energy (VSGB 2.0 implicit solvation):
       - Li, J., Abel, R., Zhu, K., Cao, Y., Zhao, S. & Friesner, R.A. (2011)
         The VSGB 2.0 model. Proteins 79:2794–2812. DOI: https://doi.org/10.1002/prot.23106
    3. Molecular dynamics engine + SID/Event Analysis (Schrödinger Desmond):
       - Bowers, K.J. et al. (2006) SC'06: Proc ACM/IEEE Conf Supercomputing.
       - DOI: https://doi.org/10.1109/SC.2006.54
    4. MD force field (OPLS4):
       - Lu, C. et al. (2021) J Chem Theory Comput 17:4291–4300.
       - DOI: https://doi.org/10.1021/acs.jctc.1c00302
    5. Plotting:
       - Hunter, J.D. (2007) Matplotlib. Comput Sci Eng 9:90–95. DOI: https://doi.org/10.1109/MCSE.2007.55
===============================================================================
"""

import argparse
import concurrent.futures as cf
import importlib.util as _ilu
import math
import os
import shutil
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402

"""
CPU usage cap (total cores - 2; mirrors CFG.PREP_CPU_RESERVE). Reserve 2 cores
for OS/desktop stability by limiting the thread-pool maths libraries (BLAS /
MKL / OpenMP / NumExpr) used by downstream tools and any numerical imports.
setdefault() preserves any value exported by the caller or pipeline runner.
"""
_CPU_CAP = str(max(1, (os.cpu_count() or 4) - 2))
for _tv in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_tv, _CPU_CAP)

# -------------------------------------------------------------------------------
# SECTION 1: PIPELINE MODULES (00_01 CFG, 00_02 utils) via importlib
# -------------------------------------------------------------------------------
def _load_module(name: str, path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Required module not found: {path}")
    spec = _ilu.spec_from_file_location(name, str(path))
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_SCRIPT_DIR = Path(__file__).resolve().parent
_utils_mod = _load_module("ProjectUtils", _SCRIPT_DIR / "00_02_Project_Utils_FAcDs.py")
_cfg_mod = _load_module("ProjectConfig", _SCRIPT_DIR / "00_01_Project_Config_FAcDs.py")

CFG = _cfg_mod.CFG()
print_script_banner = _utils_mod.print_script_banner
print_elapsed = _utils_mod.print_elapsed
apply_figure_style = _utils_mod.apply_figure_style
auto_label_colour = _utils_mod.auto_label_colour
apply_figure_style(CFG)   # one typography definition for every figure the pipeline draws

# Auto-set SCHRODINGER if the env var is absent (mirrors the 08 engine default).
os.environ.setdefault("SCHRODINGER", "/opt/schrodinger")
SCHRODINGER = os.environ["SCHRODINGER"]
SCHROD_RUN = os.path.join(SCHRODINGER, "run")

_SEP = "============================================================================="
_RULE = "─────────────────────────────────────────────────────────────────────────"
_DEFAULT_FRAME_TOTAL = 100_000   # heartbeat fallback when the trajectory length is unreadable
EXIT_WARN = 3   # step completed but a complementary part (MM-GBSA) was deferred/failed; the pipeline runner renders WARN and continues (0=PASS, 1=hard error, 3=warn)

# Tokens the SID-out.eaf Result vector carries per trajectory frame. Governs the
# completion gate (is_eaf_complete): threshold = EAF_TOKENS_PER_FRAME · traj_frames.
# 1 is the standard per-frame scalar series; set to the measured ratio once read
# off a known-good, fully-analysed EAF so a partial high-k EAF cannot pass.
EAF_TOKENS_PER_FRAME = int(CFG.EAF_TOKENS_PER_FRAME)


# =============================================================================
# SECTION 2: SMALL HELPERS
# =============================================================================
# Main-log file handle for this step, opened in main() once the run directory is known. Every _echo
# line is mirrored to it (ANSI stripped), so Step 06 has the same 00_<StepName>.log every other step
# writes. The per-\r progress bars (written straight to sys.stdout) are deliberately NOT mirrored —
# a log does not want carriage-return redraws.
_LOG_FH = None
_ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def _open_step_log(physics_dir: Path) -> None:
    """Open 6_Physics_Validation/00_SID_MMGBSA.log for this run (fresh each run)."""
    global _LOG_FH
    try:
        physics_dir.mkdir(parents=True, exist_ok=True)
        _LOG_FH = open(physics_dir / "00_SID_MMGBSA.log", "w", encoding="utf-8")
        import atexit
        atexit.register(lambda: _LOG_FH and not _LOG_FH.closed and _LOG_FH.close())
    except Exception:
        _LOG_FH = None


def _echo(msg: str = "") -> None:
    """Print to terminal immediately (flush) — keeps live progress visible — and mirror to the
    step log file with ANSI colour codes stripped."""
    print(msg, flush=True)
    if _LOG_FH is not None:
        try:
            _LOG_FH.write(_ANSI_RE.sub("", str(msg)) + "\n")
            _LOG_FH.flush()
        except Exception:
            pass


def out_eaf_frames(eaf_path: Path) -> int:
    """Count frames written into a SID-out.eaf (the analysed result vector).

    Parses the ``Result = [ ... ]`` block and counts whitespace-separated
    tokens. Returns 0 if the file is missing or the block is absent.
    """
    try:
        content = Path(eaf_path).read_text(errors="ignore")
    except Exception:
        return 0
    m = re.search(r"Result\s*=\s*\[([\d\.\-\s+eE]+)\]", content)
    return len(m.group(1).split()) if m else 0


def traj_frame_count(trj_dir: Path) -> int:
    """Ground-truth frame count of a Desmond trajectory via the Schrödinger
    traj API. This is the denominator for completion — it adapts automatically
    whether the trajectory holds 1,000 or 100,000 frames. Returns 0 if
    unreadable (e.g. the directory is absent or Schrödinger cannot read it).
    """
    trj_dir = Path(trj_dir)
    if not trj_dir.is_dir():
        return 0
    pycode = (
        "import sys\n"
        "try:\n"
        "    from schrodinger.application.desmond.packages import traj\n"
        "    print(len(traj.read_traj(sys.argv[1])))\n"
        "except Exception:\n"
        "    print(0)\n"
    )
    try:
        res = subprocess.run(
            [SCHROD_RUN, "python3", "-c", pycode, str(trj_dir)],
            capture_output=True, text=True, timeout=600,
        )
        return int((res.stdout or "0").strip().split()[0])
    except Exception:
        return 0


def traj_span_ns(trj_dir: Path) -> float:
    """Wall-clock length of a Desmond trajectory in ns, read from the frames' own
    timestamps (chemical time is stored in ps). Nothing about the simulation length is
    assumed — a 100 ns and a 1000 ns run both report themselves correctly. Returns 0.0
    when unreadable, in which case callers simply omit the timing commentary.
    """
    trj_dir = Path(trj_dir)
    if not trj_dir.is_dir():
        return 0.0
    pycode = (
        "import sys\n"
        "try:\n"
        "    from schrodinger.application.desmond.packages import traj\n"
        "    fr = traj.read_traj(sys.argv[1])\n"
        "    print((fr[-1].time - fr[0].time) / 1000.0)\n"
        "except Exception:\n"
        "    print(0.0)\n"
    )
    try:
        res = subprocess.run(
            [SCHROD_RUN, "python3", "-c", pycode, str(trj_dir)],
            capture_output=True, text=True, timeout=600,
        )
        return float((res.stdout or "0").strip().split()[0])
    except Exception:
        return 0.0


def is_eaf_complete(out_frames: int, traj_frames: int,
                    tokens_per_frame: int = EAF_TOKENS_PER_FRAME) -> bool:
    """A SID-out.eaf is complete when its ``Result=[…]`` vector holds one full
    token per trajectory frame, scaled by ``tokens_per_frame``.

    The completion threshold is ``tokens_per_frame · traj_frames``. With the
    default ratio of 1 this reduces to the plain per-frame test. If a SID
    analysis writes k>1 tokens per frame, a partial EAF at <100 % progress can
    still reach ``traj_frames`` tokens and be misjudged complete under a bare
    ``>=`` test; scaling the denominator by the true k closes that gap.

    Set ``EAF_TOKENS_PER_FRAME`` (module constant) to the measured ratio once it
    has been read off a known-good, fully-analysed EAF; until then it is 1.

    If the trajectory count cannot be read (traj_frames == 0), an existing
    non-empty EAF is treated as complete — never destructively re-run on doubt.
    """
    if traj_frames > 0:
        return out_frames >= max(1, tokens_per_frame) * traj_frames
    return out_frames > 0


def _natural_rank(job_dir: Path) -> int:
    """Sort key for desmond_md_job_R_N directories (numeric, -V style).

    Matches the '_R_N' rank token; the '_Rank_N' spelling is also accepted.
    """
    m = re.search(r"_R(?:ank)?_(\d+)", job_dir.name)
    return int(m.group(1)) if m else 0


def _rank_of(job_name: str) -> str:
    """Extract the rank token after the '_md_job_R_' (or '_md_job_Rank_') prefix."""
    m = re.search(r"_md_job_R(?:ank)?_(.+)$", job_name)
    return m.group(1) if m else job_name.split("_md_job_")[-1]


def _proc_alive(pattern: str) -> bool:
    """True if a process whose command line matches `pattern` is running
    (equivalent to `pgrep -f`). Excludes this very process.
    """
    try:
        res = subprocess.run(["pgrep", "-f", pattern], capture_output=True, text=True)
    except FileNotFoundError:
        return False
    if res.returncode != 0:
        return False
    me = str(os.getpid())
    return any(pid.strip() and pid.strip() != me for pid in res.stdout.split())


# =============================================================================
# SECTION 3: SYSTEMD-OOMD MANAGEMENT (standalone mode only)
# =============================================================================
class OomdGuard:
    """Optionally mask systemd-oomd to prevent Out-Of-Memory kills during long SID runs.

    Standalone mode only, and OPTIONAL: sudo is offered but not required. If sudo
    is unavailable or the user skips it, the guard self-disables (`active=False`)
    and the run proceeds in NORMAL mode (oomd not masked). When called via the
    pipeline runner the runner owns oomd masking across Step 06 (this script) and
    Step 07, so this guard is a no-op. A background thread keeps the credential warm.
    """

    def __init__(self, active: bool):
        self.active = active
        self._restored = False
        self._keepalive_stop = threading.Event()
        self._keepalive_thread: threading.Thread | None = None

    def __enter__(self):
        if not self.active:
            _echo("  [PIPELINE-MODE] oomd management delegated to pipeline runner.")
            return self
        _echo("OPTIONAL — protect this run from the Linux out-of-memory killer.")
        _echo("  SID and Prime MM-GBSA both hold large trajectories in memory for hours, and "
              "systemd-oomd can kill them mid-run. Masking it needs root.")
        _echo("  Enter your sudo password to mask systemd-oomd, or press Enter / Ctrl-D to skip "
              "and run unprotected.")
        # Prime the sudo credential interactively against the controlling tty (optional).
        _sudo_ok = False
        try:
            with open("/dev/tty") as _tty:
                _sudo_ok = subprocess.run(["sudo", "-v"], stdin=_tty).returncode == 0
        except Exception:
            try:
                _sudo_ok = subprocess.run(["sudo", "-v"]).returncode == 0
            except Exception:
                _sudo_ok = False
        if not _sudo_ok:
            _echo("  [NORMAL MODE] sudo unavailable/skipped — systemd-oomd NOT masked (run unprotected from the OOM-killer).")
            self.active = False
            return self

        def _keepalive():
            while not self._keepalive_stop.wait(60):
                subprocess.run(["sudo", "-vn"], capture_output=True)

        self._keepalive_thread = threading.Thread(target=_keepalive, daemon=True)
        self._keepalive_thread.start()

        _echo("  Masking systemd-oomd — it will be restored automatically when Step 06 exits.")
        subprocess.run(["sudo", "systemctl", "stop", "systemd-oomd"], capture_output=True)
        subprocess.run(["sudo", "systemctl", "mask", "systemd-oomd.socket"], capture_output=True)
        return self

    def restore(self):
        if not self.active or self._restored:
            return
        _echo("Restoring systemd-oomd services...")
        subprocess.run(["sudo", "systemctl", "unmask", "systemd-oomd.socket"], capture_output=True)
        subprocess.run(["sudo", "systemctl", "start", "systemd-oomd"], capture_output=True)
        self._restored = True

    def __exit__(self, exc_type, exc, tb):
        self._keepalive_stop.set()
        self.restore()
        return False


# =============================================================================
# SECTION 4: HEARTBEATS (in-place progress for the long-running steps)
#   Heartbeat      — one process, one log  (SID; serial MM-GBSA)
#   ShardHeartbeat — many shard logs at once (sharded MM-GBSA)
# =============================================================================
class Heartbeat:
    """Periodically report a long Schrödinger step's progress by tailing its log.

    Progress stays on ONE continuously-updating line: each "[PROGRESS] …" tick is
    rewritten in place with a carriage return (\\r), never a new line, so the frame
    count / percentage / elapsed minutes advance on a single refreshing line. SID
    and MM-GBSA run sequentially (one Schrödinger step at a time), so nothing
    interleaves. The ONLY newline breaks are: a one-time colour banner announcing
    the MM-GBSA PHASE 1→2 transition (trajectory read done → Prime minimisation),
    and a closing newline when the step ENDS (completion or failure, via __exit__).
    A live terminal (06's stdout is teed, so \\r reaches it) shows the banner then a
    single advancing line for the whole Prime phase.
    Each step prints its own per-frame marker; pass the matching regex(es):
      • analyze_simulation.py (SID)  → "analysing frame# N"   (default)
      • thermal_mmgbsa.py (MM-GBSA)  → "Reading frame N" / "Structure N"
    The latest captured N is shown as a percentage of the known trajectory length.
    """

    # ANSI styling for the PHASE 2 banner (passes through the pipeline tee to a
    # live terminal; harmless byte-noise in a pure-file redirect).
    _C_PHASE2 = "\033[1;36m"   # bold cyan
    _C_DIM    = "\033[2m"      # dim
    _C_RST    = "\033[0m"      # reset

    def __init__(self, log_file: Path, label: str, total: int,
                 patterns: list[str] | None = None, interval: int = 120,
                 read_phase: str = ""):
        self.log_file = Path(log_file)
        self.label = label
        # Wording for the frame-reading phase. SID reads frames for its whole run, so it
        # leaves this empty; MM-GBSA labels it "PHASE 1/2" to distinguish it from the Prime
        # minimisation that follows.
        self.read_phase = f"{read_phase} " if read_phase else ""
        self.total = max(int(total), 1)
        self.patterns = [re.compile(p) for p in (patterns or [r"analyzing frame# (\d+)"])]
        self.interval = max(5, int(interval))
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._start = 0.0
        self._line_open = False         # an in-place \r progress line is currently open
        self._prime_announced = False   # PHASE 2 banner printed once

    def _progress(self, msg: str) -> None:
        """Rewrite the single progress line in place with a carriage return."""
        sys.stdout.write(f"\r    [PROGRESS] {msg}\033[K")
        sys.stdout.flush()
        self._line_open = True

    def _end_line(self) -> None:
        """Close the open in-place \r line with a newline (phase change / exit)."""
        if self._line_open:
            sys.stdout.write("\n")
            sys.stdout.flush()
            self._line_open = False

    def _latest_frame(self, text: str) -> int:
        """Largest last-match across all configured patterns (phase-tolerant)."""
        best = 0
        for pat in self.patterns:
            m = pat.findall(text)
            if m:
                best = max(best, int(m[-1]))
        return best

    # thermal_mmgbsa hands the read frames to a separate Prime MM-GBSA job whose
    # progress is written elsewhere; the main log falls silent. Detect that hand-off
    # so the heartbeat reports the Prime scoring phase instead of a frozen frame count.
    _PRIME_HANDOFF = re.compile(r"Passing\s+\d+\s+structures\s+to\s+Prime|Running Prime MMGBSA", re.IGNORECASE)

    def _prime_progress(self) -> "tuple[int, int] | None":
        """Best-effort (done, total) for the Prime phase from a ``*-prime*.log`` in
        the job directory, if the job server wrote one there; otherwise ``None``."""
        try:
            for _plog in sorted(self.log_file.parent.glob("*-prime*.log")):
                _t = _plog.read_text(errors="ignore")
                m = re.findall(r"(\d+)\s*(?:/|of)\s*(\d+)\s*(?:sub)?jobs?", _t, re.IGNORECASE)
                if m:
                    _d, _tot = m[-1]
                    return int(_d), int(_tot)
        except Exception:
            pass
        return None

    # Prime per-structure marker, so PHASE 2 can report WHICH structure/frame is
    # being minimised (like PHASE 1's frame count) rather than an opaque "working".
    _STRUCT_RE = re.compile(r"[Ss]tructure\s+(?:#\s*)?(\d+)(?:\s+of\s+(\d+))?")

    def _latest_structure(self, text: str) -> "tuple[int, int] | None":
        """Highest 'Structure N (of M)' seen in the main log or any ``*-prime*.log``
        the job server writes. Returns (n, total) — total falls back to the frame
        count — or ``None`` if Prime has not logged a per-structure marker yet."""
        best_n, best_tot = 0, 0
        sources = [text]
        try:
            for _p in sorted(self.log_file.parent.glob("*-prime*.log")):
                sources.append(_p.read_text(errors="ignore"))
        except Exception:
            pass
        for _s in sources:
            for _m in self._STRUCT_RE.finditer(_s):
                _n = int(_m.group(1))
                if _n > best_n:
                    best_n = _n
                    if _m.group(2):
                        best_tot = int(_m.group(2))
        if best_n > 0:
            return best_n, (best_tot or self.total)
        return None

    def _run(self):
        while not self._stop.wait(self.interval):
            elapsed_min = int((time.time() - self._start) // 60)
            try:
                text = self.log_file.read_text(errors="ignore")
            except Exception:
                text = ""
            # Prime MM-GBSA scoring phase: trajectory read is complete; Prime now
            # minimises every structure. Announce the transition once (colour banner
            # on its own line), then keep the single \r line refreshing — reporting
            # WHICH structure is being minimised (structure N/total, like PHASE 1's
            # frame count), falling back to subjob count, then elapsed time.
            if self._PRIME_HANDOFF.search(text):
                if not self._prime_announced:
                    self._prime_announced = True
                    self._end_line()   # close the PHASE 1 \r line with a newline
                    _echo(f"{self._C_PHASE2}    ══════ [PHASE 2/2] {self.label}: every trajectory "
                          f"frame has been read. Prime is now minimising and scoring "
                          f"{self.total:,} structures ══════{self._C_RST}")
                    _echo(f"{self._C_DIM}    This is the long phase — it can run for hours. The "
                          f"line below refreshes in place with the structure Prime is on."
                          f"{self._C_RST}")
                _struct = self._latest_structure(text)
                if _struct is not None:
                    _n, _tot = _struct
                    _pct = min(100, _n * 100 // max(_tot, 1))
                    self._progress(f"{self.label}: PHASE 2/2 Prime minimised structure {_n:,} of "
                                   f"{_tot:,} ({_pct}%, {elapsed_min}m elapsed)")
                else:
                    prog = self._prime_progress()
                    if prog is not None:
                        _done, _tot = prog
                        _pct = min(100, _done * 100 // max(_tot, 1))
                        self._progress(f"{self.label}: PHASE 2/2 Prime finished {_done} of {_tot} "
                                       f"subjobs ({_pct}%, {elapsed_min}m elapsed)")
                    else:
                        # Prime spinning up — no per-structure marker logged yet.
                        self._progress(f"{self.label}: PHASE 2/2 Prime is starting up on "
                                       f"{self.total:,} structures ({elapsed_min}m elapsed)")
                continue
            # Frame-read progress (SID's whole run, and MM-GBSA PHASE 1). Kept
            # phase-neutral so it reads correctly for SID, which has no PHASE 2.
            frame = self._latest_frame(text)
            pct = min(100, frame * 100 // self.total)
            self._progress(f"{self.label}: {self.read_phase}read frame {frame:,} of "
                           f"{self.total:,} ({pct}%, {elapsed_min}m elapsed)")

    def __enter__(self):
        self._start = time.time()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1)
        self._end_line()   # close the single \r progress line (completion or failure)
        return False


class ShardHeartbeat:
    """One in-place progress line for a sharded MM-GBSA run.

    Aggregates across the shard logs rather than following a single process: frames read
    (summed over every shard), how many shards have already been scored by Prime, and how
    many are minimising right now. Mirrors Heartbeat's single-\\r-line contract.
    """

    def __init__(self, shard_dir: Path, label: str, n_shards: int, total_frames: int,
                 done_fn, interval: int = 30):
        self.shard_dir = Path(shard_dir)
        self.label = label
        self.n_shards = max(1, int(n_shards))
        self.total_frames = max(1, int(total_frames))
        self.done_fn = done_fn
        self.interval = max(5, int(interval))
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._start = 0.0
        self._line_open = False

    def _counts(self) -> "tuple[int, int]":
        """(frames read across all shards, shards currently in Prime minimisation)."""
        read = priming = 0
        for log in self.shard_dir.glob("*_shard*.log"):
            try:
                text = log.read_text(errors="ignore")
            except Exception:
                continue
            read += text.count("Reading frame")
            if "Running Prime MMGBSA job" in text:
                priming += 1
        return read, priming

    def _run(self) -> None:
        while not self._stop.wait(self.interval):
            mins = int((time.time() - self._start) // 60)
            read, priming = self._counts()
            scored = self.done_fn()
            pct = min(100, read * 100 // self.total_frames)
            priming = max(0, priming - scored)
            sys.stdout.write(
                f"\r    [PROGRESS] {self.label}: read {read:,} of {self.total_frames:,} frames "
                f"({pct}%) · {scored}/{self.n_shards} shards scored · {priming} minimising in "
                f"Prime · {mins}m elapsed\033[K")
            sys.stdout.flush()
            self._line_open = True

    def __enter__(self):
        self._start = time.time()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1)
        if self._line_open:
            sys.stdout.write("\n")
            sys.stdout.flush()
        return False


# =============================================================================
# SECTION 5: RUN-DIRECTORY RESOLUTION
# =============================================================================
def resolve_run_dir(run_arg: str | None) -> str:
    """Return the Boltz-2_Run_* directory name, auto-detecting the latest if
    none was supplied. Mirrors the shell auto-detect (sort | tail -1).
    """
    if run_arg:
        return run_arg
    candidates = sorted(
        d.name for d in _SCRIPT_DIR.iterdir()
        if d.is_dir() and d.name.startswith("Boltz-2_Run_")
    )
    if candidates:
        chosen = candidates[-1]
        _echo(f"Auto-detected run directory: {chosen}")
        return chosen
    _echo("ERROR: Run directory not specified and no Boltz-2_Run_* found.")
    _echo("Usage: python 06_SID_Prime-MMGBSA_FAcDs.py [Boltz-2_Run_Directory]")
    sys.exit(1)


# =============================================================================
# SECTION 6: SCAN PHASE
# =============================================================================
def scan_jobs(job_dirs: list[Path], md_dir: Path) -> list[Path]:
    """Classify every desmond_md_job_R_* directory and return the subset
    that still needs SID analysis (READY or INCOMPLETE).
    """
    _echo("Scanning MD job directories to see which ones still need SID analysis...")
    to_run: list[Path] = []
    completed = running = pending = 0

    for d in job_dirs:
        job_name = d.name
        rank = _rank_of(job_name)
        out_eaf = d / f"{job_name}{CFG.SUFFIX_SID_EAF}"
        cms_file = d / f"{job_name}{CFG.SUFFIX_CMS_OUT}"
        trj_dir = d / f"{job_name}_trj"

        """
        1. No -out.cms => MD simulation not finished. A live process here means
        the SIMULATION is running (job-name match avoids SID processes).
        """
        if not cms_file.is_file():
            if _proc_alive(job_name):
                running += 1
                _echo(f"  - Rank {rank}: MD simulation is still running — SID must wait for it.")
            else:
                pending += 1
                _echo(f"  - Rank {rank}: no MD output (-out.cms) and nothing running — SID cannot start.")
            continue

        """
        2. -out.cms present => MD done. Classify by SID-out.eaf completeness,
        measured against the actual trajectory length.
        """
        if out_eaf.is_file():
            of = out_eaf_frames(out_eaf)
            tf = traj_frame_count(trj_dir)
            if is_eaf_complete(of, tf):
                completed += 1
                _echo(f"  - Rank {rank}: SID already complete — all {tf:,} trajectory frames analysed.")
            else:
                to_run.append(d)
                _echo(f"  - Rank {rank}: SID stopped early — only {of:,} of {tf:,} frames analysed, "
                      f"so it will be re-run from scratch.")
        else:
            to_run.append(d)
            _echo(f"  - Rank {rank}: MD finished but SID has never been run — queued.")

    _echo("")
    _echo("Scan summary:")
    _echo(f"  MD jobs found                  : {len(job_dirs)}")
    _echo(f"  SID already complete           : {completed}")
    _echo(f"  MD simulation still running    : {running}")
    _echo(f"  MD simulation not started yet  : {pending}")
    _echo(f"  SID analysis to run now        : {len(to_run)}")
    _echo(_SEP)
    return to_run


# =============================================================================
# SECTION 7: SID ANALYSIS PER JOB
# =============================================================================
def run_event_analysis(job_dir: Path, job_name: str, cms_file: Path, in_eaf: Path) -> None:
    """Step 1: event_analysis.py — generates the SID-in.eaf descriptor."""
    if in_eaf.is_file():
        _echo(f"  [SKIP] SID-in.eaf already exists: {in_eaf}")
        return
    _echo("  Running event_analysis.py to generate SID-in.eaf...")
    log = job_dir / f"{job_name}_event_analysis.log"
    with open(log, "w") as fh:
        subprocess.run(
            [SCHROD_RUN, "event_analysis.py", "analyze", str(cms_file),
             "-out", f"{job_name}_SID"],
            cwd=str(job_dir), stdout=fh, stderr=subprocess.STDOUT, check=True,
        )


def run_analyze_simulation(job_dir: Path, job_name: str, cms_file: Path,
                           trj_dir: Path, out_eaf: Path, label: str) -> int:
    """Step 2: analyze_simulation.py — generates the SID-out.eaf result vector.

    Uses -LOCAL to avoid remote-server overhead and cap memory. Returns the
    trajectory frame total used as the completion denominator.
    """
    total_frames = traj_frame_count(trj_dir) or _DEFAULT_FRAME_TOTAL
    _echo(f"  Running analyze_simulation.py (long step, ~10-15h; {total_frames} frames)...")
    log = job_dir / f"{job_name}_analyze_simulation.log"
    with Heartbeat(log, label, total_frames):
        with open(log, "w") as fh:
            subprocess.run(
                [SCHROD_RUN, "analyze_simulation.py", "-NOJOBID", "-LOCAL",
                 str(cms_file), str(trj_dir),
                 f"{job_name}_SID-out.eaf", f"{job_name}_SID-in.eaf"],
                cwd=str(job_dir), stdout=fh, stderr=subprocess.STDOUT, check=True,
            )
    return total_frames


def process_jobs(to_run: list[Path]) -> None:
    """Iterate the work list, generating SID-in then SID-out for each job.

    Each rank is isolated: a Schrödinger failure (event_analysis/analyze_simulation)
    or an incomplete EAF is logged and skipped so the remaining ranks still run.
    Step 07 consumes whichever *_SID-out.eaf files completed.
    """
    total = len(to_run)
    failures: list[str] = []
    for index, d in enumerate(to_run, start=1):
        job_name = d.name
        rank = _rank_of(job_name)
        out_eaf = d / f"{job_name}_SID-out.eaf"
        in_eaf = d / f"{job_name}_SID-in.eaf"
        cms_file = d / f"{job_name}-out.cms"
        trj_dir = d / f"{job_name}_trj"

        _echo("")
        _echo(_SEP)
        # Re-run (incomplete) vs. fresh processing — preserve the partial EAF.
        if out_eaf.is_file():
            of = out_eaf_frames(out_eaf)
            tf = traj_frame_count(trj_dir)
            _echo(f"[Job {index}/{total}] Re-processing (Incomplete): {job_name} (Rank {rank})")
            _echo(_SEP)
            _echo(f"  [RE-RUN] Output EAF is incomplete ({of}/{tf} frames). Re-running analysis...")
            # A misjudged "incomplete" must never destroy a good result.
            stamp = time.strftime("%Y%m%d%H%M%S")
            bak = out_eaf.with_name(f"{out_eaf.name}.incomplete.{stamp}.bak")
            out_eaf.replace(bak)
            _echo(f"  [BACKUP] Previous EAF moved to: {bak}")
        else:
            _echo(f"[Job {index}/{total}] Processing: {job_name} (Rank {rank})")
            _echo(_SEP)

        # Verify the trajectory folder is present (dir or packed .xtc).
        if not trj_dir.is_dir() and not (trj_dir.parent / f"{trj_dir.name}.xtc").is_file():
            _echo(f"  [WARNING] Trajectory folder not found: {trj_dir}. Skipping Rank {rank}.")
            continue

        label = f"Job {index}/{total} Rank {rank}"
        try:
            run_event_analysis(d, job_name, cms_file, in_eaf)
            total_frames = run_analyze_simulation(d, job_name, cms_file, trj_dir, out_eaf, label)
        except (subprocess.CalledProcessError, OSError) as exc:
            _echo(f"  [ERROR] SID analysis failed for Rank {rank} ({job_name}): {exc}. Skipping rank.")
            failures.append(f"Rank {rank} ({job_name})")
            continue

        if out_eaf.is_file() and is_eaf_complete(out_eaf_frames(out_eaf), total_frames):
            of = out_eaf_frames(out_eaf)
            _echo(f"  [SUCCESS] Created output: {out_eaf} ({of}/{total_frames} frames)")
            _echo(f"  [LOGS] Logs written to {job_name}_event_analysis.log "
                  f"and {job_name}_analyze_simulation.log")
        else:
            _echo(f"  [ERROR] Incomplete output for Rank {rank}: {out_eaf}. Skipping rank "
                  f"(re-run to retry — the partial EAF is preserved).")
            failures.append(f"Rank {rank} ({job_name})")
            continue

    if failures:
        _echo("")
        _echo(_SEP)
        _echo(f"  [SUMMARY] {len(failures)}/{total} rank(s) did not produce a complete SID-out.eaf:")
        for fitem in failures:
            _echo(f"            - {fitem}")
        _echo("  Step 07 will process only the ranks that completed.")


# =============================================================================
# SECTION 8: PRIME MM-GBSA (end-state binding free energy over the MD ensemble)
#   8.1  result discovery and failure diagnosis
#   8.2  execution — sharded (default) and serial
#   8.3  statistics — ΔG estimators, failed-minimisation flagging
#   8.4  figures — per job and combined
#   8.5  phase driver
# =============================================================================
"""
Mirrors the manual workflow:
    cd <desmond_md_job_R_N> ; $SCHRODINGER/run thermal_mmgbsa.py <job>-out.cms
run folder-wise and idempotent (skip when a valid results CSV already exists),
then plots per-job + combined ΔG_bind. MM-GBSA is complementary to the QSite
QM/MM reaction barrier (Step 07), not a replacement: it scores binding, not
bond cleavage. GB implicit solvent overstabilises anionic PFAS → relative only.
"""

# ── 8.1  Result discovery and failure diagnosis ──────────────────────────────
def _mmgbsa_csv(job_dir: Path, job_name: str):
    """Locate a thermal_mmgbsa results CSV in `job_dir` (name varies by version).
    Returns the newest match or None."""
    hits = []
    for pat in (f"{job_name}{CFG.SUFFIX_MMGBSA_CSV}", f"{job_name}*prime*mmgbsa*.csv",
                "*prime*mmgbsa*.csv", "*mmgbsa*.csv"):
        hits += list(job_dir.glob(pat))
    hits = [h for h in dict.fromkeys(hits) if h.is_file() and h.stat().st_size > 0]
    return sorted(hits, key=lambda p: (-p.stat().st_mtime, p.name))[0] if hits else None


def _diagnose_mmgbsa_failure(job_dir: Path, job_name: str) -> "str | None":
    """Best-effort human-readable cause when MM-GBSA exits non-zero, read from the
    thermal_mmgbsa log and any Prime subjob logs. Turns a bare 'rc=1' into an
    actionable line (the disk-full case is silent in the summary otherwise)."""
    blob = []
    for cand in (job_dir / f"{job_name}_mmgbsa.log", *sorted(job_dir.glob("*-prime*.log"))):
        try:
            blob.append(cand.read_text(errors="ignore"))
        except Exception:
            pass
    text = "\n".join(blob)
    if re.search(r"no space left on device|copy_file_range.*no space", text, re.I):
        return ("DISK FULL — the Schrödinger job server ran out of space staging "
                "per-subjob scratch (each Prime subjob copies the multi-GB complexes "
                "file). This is the job-SERVER directory filling up, not SCHRODINGER_TMPDIR; "
                "relocate it onto the working disk with "
                "`jsc local-server-dir --set <working-disk>` (server stopped), then retry.")
    if re.search(r"licen[sc]e", text, re.I) and re.search(r"error|fail|not available|checkout", text, re.I):
        return "LICENSE — a Prime/PLOP (PSP_PLOP) license was unavailable; check FlexLM."
    m = re.search(r"^ERROR:.*$", text, re.M)
    return m.group(0).strip() if m else None


# ── 8.2  Execution: shard planning, frame stamping, the sharded/serial runners ─
def _avail_ram_gb() -> float:
    """Free RAM the kernel expects to hand out without swapping (MemAvailable)."""
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 2**20
    except Exception:
        pass
    return 0.0


def _shard_plan(total: int, ncpu: int, step: int = 1) -> "tuple[list[tuple[int, int]], int, int]":
    """Split [0, total) into contiguous frame shards and decide how many run at once.

    Returns (ranges, concurrency, prime_njobs_per_shard). Concurrency × prime_njobs is
    kept at ncpu so no core sits idle, then clamped so the concurrent shard readers fit
    in free RAM (a reader's RSS grows with the frames it has read).
    """
    size = max(1, int(getattr(CFG, "MMGBSA_SHARD_FRAMES", 2000)))
    ranges = [(a, min(a + size, total)) for a in range(0, total, size)]

    njobs = max(1, int(getattr(CFG, "MMGBSA_SHARD_PRIME_NJOBS", 6)))
    conc = int(getattr(CFG, "MMGBSA_SHARD_CONCURRENCY", 0) or 0) or max(1, ncpu // njobs)
    conc = min(conc, len(ranges))

    """
    RAM budget. A shard costs its reader PLUS the njobs Prime subjobs it spawns, and Prime is
    the expensive half (~1.8 GB per subjob × 30 = ~54 GB — enough to exhaust a 60 GB box on its
    own). Budgeting only the readers is how a 30-subjob run ends up with 2 GB of headroom and
    starts swapping. Shrink concurrency until reader + Prime fit inside the free RAM, keeping
    cores busy only to the extent memory allows: an OOM-killed subjob costs a whole shard.
    """
    _ram = _avail_ram_gb() * float(getattr(CFG, "MMGBSA_RAM_HEADROOM_FRAC", 0.85))
    """
    A reader's RSS is a fixed base plus growth with the frames it actually reads — and a stride
    means it reads only size/step of them. Budgeting the every-frame figure would reserve memory
    that Prime could otherwise use.
    """
    _read_n = max(1, size // max(1, step))
    _reader = (float(getattr(CFG, "MMGBSA_READER_RAM_BASE_GB", 3.6))
               + float(getattr(CFG, "MMGBSA_READER_RAM_PER_1K_FRAMES_GB", 0.36)) * _read_n / 1000.0)
    _prime = max(0.2, float(getattr(CFG, "MMGBSA_PRIME_RAM_GB", 1.8)))
    if _ram > 0:
        while conc > 1 and conc * (_reader + njobs * _prime) > _ram:
            conc -= 1
        # A single shard that still cannot fit trims its own Prime subjobs instead.
        while njobs > 1 and conc * (_reader + njobs * _prime) > _ram:
            njobs -= 1
    else:
        njobs = max(njobs, ncpu // max(1, conc))   # no RAM reading: fall back to filling cores
    return ranges, conc, njobs


# Every mark colour in this step's figures comes from CFG (candidate colours from
# MMGBSA_RANK_PALETTE, everything else from MMGBSA_INK) — none is written here.
_INK = CFG.MMGBSA_INK

MMGBSA_FRAME_COL = "Frame"   # trajectory frame index each scored structure came from


def _stamp_frames(df: "pd.DataFrame", start: int, end: int, step: int,
                  label: str = "") -> "pd.DataFrame":
    """Record which trajectory frame each Prime row actually came from.

    thermal_mmgbsa scores frames range(start, end, step) in order, and Prime preserves that
    order, so row i is frame start + i·step. Writing that out explicitly is the whole point:
    downstream code (07's NAC-conditioned MM-GBSA) must never infer the frame from the row
    POSITION, because with a stride row 5,000 is frame 50,000 — a silent misattribution.
    The column is only stamped when the row count matches the frame range exactly; a
    mismatch means the assumption is broken and a wrong index is worse than none.
    """
    frames = list(range(int(start), int(end), max(1, int(step))))
    if len(df) == len(frames):
        df.insert(0, MMGBSA_FRAME_COL, frames)
    else:
        _echo(f"  [WARN] {label}: Prime returned {len(df)} rows for {len(frames)} requested "
              f"frames — frame indices NOT stamped. 07's NAC-conditioned MM-GBSA will fall "
              f"back to positional alignment, which is only valid at step_size ≤ 1.")
    return df


def _retrofit_frame_stamps(csv: Path, job_dir: Path, job_name: str, rank: str) -> None:
    """Add the `Frame` column to an MM-GBSA CSV that was written before stamping existed.

    Reconstructs the frame list exactly as the sharded run generated it — contiguous
    CFG.MMGBSA_SHARD_FRAMES blocks, every CFG.MMGBSA_STEP_SIZE-th frame within each, concatenated
    in shard order — and writes it in place, keeping the unstamped file as a .bak. Idempotent: a
    CSV that already carries the column is left alone. Refuses to act if the row count does not
    match the reconstructed list, since a wrong frame index is worse than a missing one.
    """
    try:
        df = pd.read_csv(csv)
    except Exception:
        return
    if MMGBSA_FRAME_COL in df.columns:
        return
    total = traj_frame_count(job_dir / f"{job_name}_trj")
    step  = max(1, int(getattr(CFG, "MMGBSA_STEP_SIZE", 0) or 1))
    size  = max(1, int(getattr(CFG, "MMGBSA_SHARD_FRAMES", 2000)))
    if not total:
        return
    frames: list[int] = []
    for a in range(0, total, size):
        frames += list(range(a, min(a + size, total), step))
    if len(df) != len(frames):
        _echo(f"    [WARN] Rank {rank}: {csv.name} has {len(df):,} rows but the current shard/stride "
              f"plan implies {len(frames):,} frames — NOT stamping (a wrong frame index is worse "
              f"than a missing one). Step 07 will skip the NAC-conditioned ΔG for this job.")
        return
    shutil.copy2(csv, csv.with_suffix(".csv.unstamped.bak"))
    df.insert(0, MMGBSA_FRAME_COL, frames)
    _tmp = csv.with_suffix(".csv.tmp")
    df.to_csv(_tmp, index=False)
    _tmp.replace(csv)
    _echo(f"    ✔ Frame-stamped: {csv.name} ({len(df):,} rows) — it predates the `Frame` column; "
          f"Step 07 needs it to align NAC frames with their binding energies.")


def _diagnose_shard_failure(slog: Path) -> "str | None":
    """Turn a shard's bare rc=1 into the actual cause, read from its log."""
    try:
        text = slog.read_text(errors="ignore")
    except Exception:
        return None
    if re.search(r"requires a locally running job server|Unable to submit job", text, re.I):
        return ("the Schrödinger local job server is DOWN, so Prime could not be submitted. "
                "Start it with `$SCHRODINGER/jsc local-server-start` and re-run — the frames "
                "were read fine, only the hand-off failed.")
    if re.search(r"no space left on device", text, re.I):
        return "the disk holding the job scratch filled up."
    if re.search(r"licen[sc]e", text, re.I) and re.search(r"error|fail|not available|checkout", text, re.I):
        return "a Prime/PSP license could not be checked out."
    m = re.search(r"^(?:MMGBSA )?Error:.*$", text, re.M)
    return m.group(0).strip() if m else None


def run_mmgbsa_sharded(job_dir: Path, job_name: str, rank: str, cms_file: Path,
                       lig_asl: str, step: int, total: int, ncpu: int,
                       interval: int) -> Path | None:
    """Score the whole trajectory as concurrent contiguous frame shards.

    Each shard is a full thermal_mmgbsa run over -start_frame..-end_frame, so it reads its
    own slice and then runs its own Prime minimisation. Several shards are in flight at
    once: reading runs in parallel across shards, one shard minimises while the next is still
    reading, and a shard's memory is bounded by its slice. The shards are disjoint
    and cover [0, total), so every frame is scored — identical coverage to one serial run.
    Shard CSVs are concatenated into the per-job CSV the rest of Step 06 reads.
    """
    shard_dir = job_dir / getattr(CFG, "MMGBSA_SHARD_SUBDIR", "_MMGBSA_Shards")
    shard_dir.mkdir(exist_ok=True)
    """
    The .cms names its trajectory relatively (s_chorus_trajectory_file = '<job>_trj'), so a
    shard run from shard_dir needs both the .cms and the _trj folder visible beside it.
    Symlinks keep the shard outputs out of the job folder without copying 33 GB.
    """
    for src in (cms_file, job_dir / f"{job_name}_trj"):
        link = shard_dir / src.name
        if not link.exists():
            link.symlink_to(os.path.relpath(src, shard_dir))

    ranges, conc, njobs = _shard_plan(total, ncpu, step)
    _scored = total if step <= 1 else -(-total // step)     # structures Prime will minimise
    """
    Report the plan as aligned fields rather than prose: the numbers (frames, stride, shards,
    cores, memory) are what a reader actually checks, and they are all derived — the trajectory
    span comes from the frames' own timestamps, so nothing here assumes a run length.
    """
    _span = traj_span_ns(job_dir / f"{job_name}_trj")
    _echo(f"    Trajectory     : {total:,} frames" + (f" · {_span:,.0f} ns" if _span > 0 else ""))
    if step > 1:
        _ps = (_span * 1000.0 * step / max(total, 1)) if _span > 0 else 0.0
        _echo(f"    Sampling       : every {step}th frame → {_scored:,} structures to Prime"
              + (f", one per {_ps:,.0f} ps" if _ps else ""))
        _echo(f"                     (below the ~0.1–1 ns decorrelation time of a bound pose, so the "
              f"independent-sample count — and ⟨ΔG_bind⟩ — is unchanged)")
    else:
        _echo(f"    Sampling       : every frame → {_scored:,} structures to Prime")
    _echo(f"    Shards         : {len(ranges)} × {getattr(CFG, 'MMGBSA_SHARD_FRAMES', 2000):,} frames, "
          f"contiguous and disjoint (the whole trajectory is covered)")
    _echo(f"    Parallelism    : {conc} shards at once × {njobs} Prime subjobs = {conc * njobs} cores; "
          f"reading and minimisation overlap")
    _echo(f"    Memory         : {_avail_ram_gb():.0f} GB free · budgeted "
          f"{conc * (float(getattr(CFG, 'MMGBSA_READER_RAM_BASE_GB', 3.6)) + njobs * float(getattr(CFG, 'MMGBSA_PRIME_RAM_GB', 1.8))):.0f} GB "
          f"(reader + Prime subjobs per shard)")

    done: dict[int, Path] = {}
    failed: list[int] = []
    _lock = threading.Lock()
    """
    Fail fast. A shard that dies for an environmental reason — job server down, no license,
    disk full — will kill every other shard the same way, each only AFTER re-reading its
    frames. Abort the rank once the first full wave has failed without a single success,
    so the cause is reported in a minute instead of an hour.
    """
    _abort = threading.Event()

    def _shard_csv(i: int) -> Path:
        return shard_dir / f"{job_name}_mmgbsa_shard{i:03d}-prime-out.csv"

    def _run_shard(i: int, a: int, b: int) -> None:
        csv = _shard_csv(i)
        if csv.is_file() and csv.stat().st_size > 0:      # resume: shard already scored
            with _lock:
                done[i] = csv
            return
        if _abort.is_set():
            return
        sname = f"{job_name}_mmgbsa_shard{i:03d}"
        cmd = [SCHROD_RUN, "thermal_mmgbsa.py", cms_file.name,
               "-j", sname, "-HOST", f"localhost:{njobs}",
               "-start_frame", str(a), "-end_frame", str(b)]
        if lig_asl:
            cmd += ["-lig_asl", lig_asl]
        if step > 0:
            cmd += ["-step_size", str(step)]
        slog = shard_dir / f"{sname}.log"
        try:
            with open(slog, "w") as fh:
                rc = subprocess.Popen(cmd, cwd=str(shard_dir), stdout=fh,
                                      stderr=subprocess.STDOUT).wait()
        except Exception as e:
            _echo(f"\n  [Rank {rank}] shard {i:03d} (frames {a:,}–{b:,}) failed to launch: {e}")
            with _lock:
                failed.append(i)
            return
        if rc != 0 or not (csv.is_file() and csv.stat().st_size > 0):
            _echo(f"\n  [Rank {rank}] shard {i:03d} (frames {a:,}–{b:,}) failed rc={rc} — see {slog.name}.")
            _why = _diagnose_shard_failure(slog)
            if _why:
                _echo(f"  [Rank {rank}] ↳ cause: {_why}")
            with _lock:
                failed.append(i)
                if not done and len(failed) >= conc:
                    _abort.set()
                    _echo(f"  [Rank {rank}] Aborting: the first {len(failed)} shards all failed and "
                          f"none succeeded — this is an environment problem, not a bad frame range. "
                          f"Remaining shards skipped.")
            return
        with _lock:
            done[i] = csv

    # Denominator is the number of frames Prime will actually see (stride applied), not the
    # raw trajectory length — otherwise the read percentage caps at 100/step and never reaches 100.
    hb = ShardHeartbeat(shard_dir, f"MM-GBSA Rank {rank}", len(ranges), _scored,
                        lambda: len(done), interval=interval)
    with hb, cf.ThreadPoolExecutor(max_workers=conc) as pool:
        list(pool.map(lambda r: _run_shard(r[0], r[1][0], r[1][1]), list(enumerate(ranges))))

    if failed or len(done) != len(ranges):
        _echo(f"    ✘ Incomplete   : {len(done)}/{len(ranges)} shards scored, "
              f"{len(failed)} failed. Re-running Step 06 resumes from the finished shards "
              f"(their CSVs are kept in {shard_dir.name}/).")
        return None

    frames: list[pd.DataFrame] = []
    for i in sorted(done):
        try:
            _df = pd.read_csv(done[i])
        except Exception as e:
            _echo(f"  [Rank {rank}] shard {i:03d} CSV unreadable ({e}).")
            return None
        _a, _b = ranges[i]
        frames.append(_stamp_frames(_df, _a, _b, max(1, step), f"Rank {rank} shard {i:03d}"))
    merged = job_dir / f"{job_name}_mmgbsa-prime-out.csv"
    _all = pd.concat(frames, ignore_index=True)
    _tmp = merged.with_suffix(".csv.tmp")
    _all.to_csv(_tmp, index=False)
    _tmp.replace(merged)                  # atomic: a reader never sees a half-written CSV
    _stamped = "frame-stamped" if MMGBSA_FRAME_COL in _all.columns else "NOT frame-stamped"
    _echo(f"    ✔ Scored       : all {len(ranges)} shards → {merged.name} "
          f"({len(_all):,} rows, {_stamped})")
    return merged


def run_mmgbsa(job_dir: Path, job_name: str, rank: str) -> Path | None:
    """Run thermal_mmgbsa.py on <job>-out.cms inside its MD folder. Idempotent:
    returns the existing CSV when already computed. Returns the results CSV path
    or None on failure."""
    cms_file = job_dir / f"{job_name}-out.cms"
    if not cms_file.is_file():
        _echo(f"    ✘ Skipped      : no {cms_file.name} — the MD simulation has not finished.")
        return None
    existing = _mmgbsa_csv(job_dir, job_name)
    if existing is not None:
        """
        A CSV written before frame stamping existed carries no `Frame` column, and Step 07 must not
        infer the frame from the row position (row i is frame i·step under a stride). Stamp it in
        place — the frame list is reconstructible from the shard plan — so an already-scored run is
        brought up to the current format without re-scoring anything.
        """
        _retrofit_frame_stamps(existing, job_dir, job_name, rank)
        _echo(f"    ✔ Already done : reusing {existing.name} (delete it to force a re-score).")
        return existing

    """
    Parallelise frame subjobs across cores: total cores − reserve (same cap as
    the rest of the pipeline, CFG.GLOBAL_MAX_WORKERS = cpu_count − PREP_CPU_RESERVE).
    """
    _ncpu = max(1, int(getattr(CFG, "GLOBAL_MAX_WORKERS", max(1, (os.cpu_count() or 4) - 2))))
    """
    Disk-aware safety net on the Prime-subjob count. Each Prime subjob stages its own
    copy of the complexes file (~22 GB here), so N concurrent subjobs need ~N × 22 GB of
    scratch. The cores-minus-reserve figure above is the target; it is only reduced when
    the disk that actually holds the scratch cannot hold that many copies. Probe THAT disk
    (SCHRODINGER_TMPDIR, relocated to the run's working disk at start-up) — not /tmp, which
    is irrelevant once the job-server has been relocated and would otherwise throttle the
    run for no reason.
    """
    try:
        import shutil as _sh
        _probe = os.environ.get("SCHRODINGER_TMPDIR") or "/tmp"
        if not os.path.isdir(_probe):
            _probe = "/tmp"
        _free_gb = _sh.disk_usage(_probe).free / 2**30
        _gb_per_subjob = float(getattr(CFG, "MMGBSA_SCRATCH_GB_PER_SUBJOB", 22))
        _headroom = float(getattr(CFG, "MMGBSA_SCRATCH_HEADROOM_FRAC", 0.75))
        _cap = int(getattr(CFG, "MMGBSA_MAX_NJOBS", 0)) or _ncpu   # 0 = no ceiling beyond GLOBAL_MAX_WORKERS
        _disk_cap = max(2, min(_cap, int(_headroom * _free_gb / _gb_per_subjob)))
        if _disk_cap < _ncpu:
            _echo(f"    Scratch disk   : ~{_free_gb:,.0f} GB free — only enough for {_disk_cap} subjobs "
                  f"(~{_gb_per_subjob:.0f} GB each), so Prime is reduced {_ncpu} → {_disk_cap}. "
                  f"Running more would fill the disk and kill the job.")
            _ncpu = _disk_cap
        else:
            _echo(f"    Scratch disk   : ~{_free_gb:,.0f} GB free — room for {_ncpu} Prime subjobs "
                  f"(~{_gb_per_subjob:.0f} GB each)")
    except Exception:
        pass

    _lig_asl_cfg = str(getattr(CFG, "MMGBSA_LIGAND_ASL", "") or "").strip()
    _step_cfg = int(getattr(CFG, "MMGBSA_STEP_SIZE", 0) or 0)
    _interval = max(5, int(getattr(CFG, "MMGBSA_PROGRESS_INTERVAL_SEC", 30)))
    # Ground-truth frame count (same source the SID phase uses).
    _total = traj_frame_count(job_dir / f"{job_name}_trj") or _DEFAULT_FRAME_TOTAL

    """
    Sharded path (default): concurrent frame-range shards, so the trajectory read is
    parallel and overlaps Prime instead of blocking it on one core. Same frames, same
    scores — see run_mmgbsa_sharded. CFG.MMGBSA_SHARD_FRAMES = 0 keeps the plain single
    serial thermal_mmgbsa run below.
    """
    _shard_frames = int(getattr(CFG, "MMGBSA_SHARD_FRAMES", 0) or 0)
    if _shard_frames > 0 and _total > _shard_frames:
        return run_mmgbsa_sharded(job_dir, job_name, rank, cms_file, _lig_asl_cfg,
                                  _step_cfg, _total, _ncpu, _interval)

    cmd = [SCHROD_RUN, "thermal_mmgbsa.py", cms_file.name,
           "-j", f"{job_name}_mmgbsa", "-HOST", f"localhost:{_ncpu}"]
    # Pin the ligand explicitly (ASL from CFG, via thermal_mmgbsa's -lig_asl flag) so
    # Prime scores the PFAS molecule; a small/heavily-fluorinated ligand can otherwise
    # be misassigned as solvent by auto-detection. Empty CFG value → auto-detect.
    _lig_asl = str(getattr(CFG, "MMGBSA_LIGAND_ASL", "") or "").strip()
    if _lig_asl:
        cmd += ["-lig_asl", _lig_asl]
    if getattr(CFG, "MMGBSA_STEP_SIZE", 0) and CFG.MMGBSA_STEP_SIZE > 0:
        cmd += ["-step_size", str(CFG.MMGBSA_STEP_SIZE)]
    _echo(f"  [Rank {rank}] Launching MM-GBSA across {_ncpu} parallel Prime subjobs.")
    _echo(f"  [Rank {rank}] Command: {' '.join(cmd[1:])}")
    _echo(f"  [Rank {rank}] [PHASE 1/2] thermal_mmgbsa is reading the MD trajectory frame by "
          f"frame (single-threaded, so this part is slow but steady). When every frame has "
          f"been read it hands the structures to Prime, and PHASE 2/2 begins.")
    log = job_dir / f"{job_name}_mmgbsa.log"
    _timeout = getattr(CFG, "MMGBSA_TIMEOUT_SEC", 0) or None
    _interval = max(5, int(getattr(CFG, "MMGBSA_PROGRESS_INTERVAL_SEC", 30)))
    # Ground-truth denominator (same source the SID phase uses).
    _total = traj_frame_count(job_dir / f"{job_name}_trj") or _DEFAULT_FRAME_TOTAL
    """
    thermal_mmgbsa logs "Reading frame N..." during the (serial) trajectory
    extraction phase and "Structure N" during the parallel Prime phase. The
    Heartbeat tails the log and renders ONE "[PROGRESS] …" line, rewritten in
    place (\\r) from start to finish across both phases, with a single closing
    newline when the step ends (completion or failure).
    """
    _hb_patterns = [r"Reading frame (\d+)", r"[Ss]tructure[:\s]+(\d+)",
                    r"[Ff]rame\s+(\d+)", r"Processing\s+(\d+)"]
    try:
        with open(log, "w") as fh:
            proc = subprocess.Popen(cmd, cwd=str(job_dir), stdout=fh,
                                    stderr=subprocess.STDOUT)
            with Heartbeat(log, f"MM-GBSA Rank {rank}", _total,
                           patterns=_hb_patterns, interval=_interval,
                           read_phase="PHASE 1/2"):
                if _timeout:
                    _t0 = time.time()
                    while True:
                        try:
                            proc.wait(timeout=_interval)
                            break
                        except subprocess.TimeoutExpired:
                            if (time.time() - _t0) > _timeout:
                                proc.kill()
                                try:
                                    proc.wait(timeout=10)   # reap the killed child (no zombie)
                                except Exception:
                                    pass
                                _echo(f"  [Rank {rank}] MM-GBSA timed out — skipped (see {log.name}).")
                                return None
                else:
                    proc.wait()
        if proc.returncode != 0:
            _echo(f"  [Rank {rank}] MM-GBSA exited rc={proc.returncode} — see {log.name}.")
            _diag = _diagnose_mmgbsa_failure(job_dir, job_name)
            if _diag:
                _echo(f"  [Rank {rank}] ↳ cause: {_diag}")
            return None
    except Exception as e:
        _echo(f"  [Rank {rank}] MM-GBSA failed ({e}) — see {log.name}.")
        return None
    """
    Stamp the frame index on the serial path too, so both paths hand downstream code a CSV
    that says which frame each row came from instead of leaving it to be inferred.
    """
    _csv = _mmgbsa_csv(job_dir, job_name)
    if _csv is not None:
        try:
            _df = pd.read_csv(_csv)
            if MMGBSA_FRAME_COL not in _df.columns:
                _stamp_frames(_df, 0, _total, max(1, _step_cfg), f"Rank {rank}").to_csv(_csv, index=False)
        except Exception as e:
            _echo(f"  [Rank {rank}] could not stamp frame indices on {_csv.name} ({e}).")
    return _csv


# ── 8.3  Statistics: ΔG estimators and failed-minimisation flagging ──────────
def _mmgbsa_dg_series(csv_path: Path) -> "pd.Series":
    """Extract the per-frame ΔG_bind series, tolerant of column-name variants."""
    try:
        df = pd.read_csv(csv_path)
    except Exception:
        return pd.Series(dtype=float)
    col = getattr(CFG, "MMGBSA_DG_COLUMN", "r_psp_MMGBSA_dG_Bind")
    if col not in df.columns:
        cands = [c for c in df.columns
                 if re.search(r"dg.?bind", c, re.I) or re.search(r"mmgbsa.*bind", c, re.I)]
        if not cands:
            return pd.Series(dtype=float)
        col = cands[0]
    _s = pd.to_numeric(df[col], errors="coerce")
    """
    Index the series by the TRAJECTORY frame each structure came from, not by row number.
    Under a stride, row 150 is frame 1,500 — plotting against the row index would compress a
    1000 ns run onto a meaningless 0–10,000 axis and hide where in the simulation an event
    happened. Falls back to the row index only for an unstamped (every-frame) CSV.
    """
    if MMGBSA_FRAME_COL in df.columns:
        _s.index = pd.to_numeric(df[MMGBSA_FRAME_COL], errors="coerce")
        _s.index.name = MMGBSA_FRAME_COL
    return _s.dropna()


def _boltzmann_mean_dg(dg, T: float = CFG.MMGBSA_TEMPERATURE_K) -> float:
    """Log-sum-exp ("Boltzmann") ensemble mean binding free energy over frames.

        ⟨ΔG⟩ = −RT ln( (1/N) Σ exp(−ΔGᵢ / RT) ),  via a max-shifted log-sum-exp.

    Reserve this for a NON-Boltzmann-sampled ensemble. MD frames are already drawn
    from the canonical (Boltzmann) distribution, so re-weighting them by exp(−ΔGᵢ/RT)
    double-counts the Boltzmann factor and collapses the estimate toward the single
    most negative frame; for a standard thermal MM-GBSA ensemble the arithmetic mean
    (CFG.MMGBSA_AVERAGING="mean") is the correct ensemble average. R and T come from
    CFG (GAS_CONSTANT_KCAL, MMGBSA_TEMPERATURE_K).
    """
    vals = [float(x) for x in dg if x == x]   # drop NaN
    if not vals:
        return float("nan")
    RT  = CFG.GAS_CONSTANT_KCAL * T            # kcal mol⁻¹
    xs  = [-v / RT for v in vals]              # Boltzmann exponents
    m   = max(xs)
    lse = m + math.log(sum(math.exp(x - m) for x in xs))
    return -RT * (lse - math.log(len(vals)))


def _avg_dg(dg) -> float:
    """Headline per-job ΔG_bind estimator selected by CFG.MMGBSA_AVERAGING.

    "mean" (default) is the arithmetic ensemble average — the standard thermal
    MM-GBSA estimate for an already-Boltzmann-sampled MD trajectory; "median" is
    robust to per-frame outliers (e.g. failed Prime frames); "boltzmann" is the
    log-sum-exp mean, reserved for non-canonical ensembles. Unknown value → mean.
    """
    mode = str(getattr(CFG, "MMGBSA_AVERAGING", "mean")).lower()
    if mode == "median":
        return float(pd.to_numeric(pd.Series(list(dg)), errors="coerce").median())
    if mode == "boltzmann":
        return _boltzmann_mean_dg(dg)
    return float(pd.to_numeric(pd.Series(list(dg)), errors="coerce").mean())


def _lookup_tiers(run_root: Path) -> dict:
    """Map Scientific_Rank → degrader_tier from the Step 02 ranked CSV, for tier
    colouring. Returns {} if not found (figures fall back to a neutral colour)."""
    prod = run_root / "1_Boltz2_Production"
    hits = sorted(prod.glob("*Ranked*.csv")) if prod.is_dir() else []
    if not hits:
        return {}
    try:
        df = pd.read_csv(hits[-1], usecols=["Scientific_Rank", "degrader_tier"])
        return {int(r): str(t) for r, t in zip(df["Scientific_Rank"], df["degrader_tier"])}
    except Exception:
        return {}


def _lookup_ligands(run_root: Path) -> dict:
    """Map Scientific_Rank → Ligand_Name from the Step 02 ranked CSV, so the figures name the
    PFAS species rather than only its rank. Returns {} if unavailable (labels fall back to
    the bare rank)."""
    prod = run_root / "1_Boltz2_Production"
    hits = sorted(prod.glob("*Ranked*.csv")) if prod.is_dir() else []
    if not hits:
        return {}
    try:
        df = pd.read_csv(hits[-1], usecols=["Scientific_Rank", "Ligand_Name"])
        return {int(r): str(l) for r, l in zip(df["Scientific_Rank"], df["Ligand_Name"])}
    except Exception:
        return {}


def _dg_failures(dg: "pd.Series") -> "pd.Series":
    """Frames whose ΔG_bind is a failed Prime minimisation rather than physics.

    A small PFAS ligand cannot bind at hundreds of kcal/mol; such values are blown-up
    structures. Flagged by a Tukey fence (K × IQR) so the threshold adapts to each job's own
    spread instead of being a arbitrary constant.
    """
    if dg.empty:
        return dg.iloc[0:0]
    q1, q3 = dg.quantile(0.25), dg.quantile(0.75)
    k = float(getattr(CFG, "MMGBSA_DG_OUTLIER_IQR_K", 3.0))
    lo, hi = q1 - k * (q3 - q1), q3 + k * (q3 - q1)
    return dg[(dg < lo) | (dg > hi)]


# ── 8.4  Figures ─────────────────────────────────────────────────────────────
def _block_bootstrap_median_ci(dg: "pd.Series") -> "tuple[float, float]":
    """Confidence interval on the median by MOVING-BLOCK bootstrap.

    Frames are autocorrelated, so resampling them individually would pretend there are ~10⁴
    independent samples and return an absurdly tight interval. Resampling contiguous BLOCKS
    longer than the pose correlation time keeps the correlation inside the block, so each block
    counts as roughly one independent draw and the interval is honest.
    """
    v = dg.to_numpy(dtype=float)
    n = len(v)
    """
    Block length must EXCEED the series' own correlation time, or the resampled blocks are still
    correlated with one another, the bootstrap behaves as if there were far more independent
    samples than there are, and the interval comes out too narrow. Derive it from the measured
    autocorrelation (τ ≈ n / 2N_eff) and take a multiple of that; CFG.MMGBSA_BOOTSTRAP_BLOCK is
    only a floor. A slowly-drifting ΔG series (N_eff of a few dozen) therefore gets a long block
    and an honestly wide interval, instead of a falsely precise one.
    """
    _neff = max(1.0, _effective_n(dg))
    _tau  = max(1.0, n / (2.0 * _neff))
    blk = int(max(int(getattr(CFG, "MMGBSA_BOOTSTRAP_BLOCK", 50)),
                  math.ceil(getattr(CFG, "MMGBSA_BOOTSTRAP_BLOCK_TAU_MULT", 2.0) * _tau)))
    blk = min(blk, max(2, n // 8))        # keep at least ~8 blocks, else the bootstrap degenerates
    if n < 2 * blk:
        return float("nan"), float("nan")
    reps = max(200, int(getattr(CFG, "MMGBSA_BOOTSTRAP_N", 2000)))
    ci = float(getattr(CFG, "MMGBSA_BOOTSTRAP_CI", 95.0))
    nblocks = max(1, n // blk)
    starts = np.random.default_rng(0).integers(0, n - blk, size=(reps, nblocks))
    idx = (starts[:, :, None] + np.arange(blk)[None, None, :]).reshape(reps, -1)
    meds = np.median(v[idx], axis=1)
    return float(np.percentile(meds, (100 - ci) / 2)), float(np.percentile(meds, 100 - (100 - ci) / 2))


def _cliffs_delta(a: "pd.Series", b: "pd.Series", sample: int = 4000) -> float:
    """Cliff's delta: P(a > b) − P(a < b). A non-parametric effect size that needs no
    distributional assumption and is unmoved by the outlier frames. |δ| ≥ 0.474 is
    conventionally 'large', 0.33 medium, 0.15 small."""
    rng = np.random.default_rng(0)
    x = a.to_numpy(dtype=float); y = b.to_numpy(dtype=float)
    if len(x) > sample:
        x = rng.choice(x, sample, replace=False)
    if len(y) > sample:
        y = rng.choice(y, sample, replace=False)
    gt = (x[:, None] > y[None, :]).sum()
    lt = (x[:, None] < y[None, :]).sum()
    return float((gt - lt) / (len(x) * len(y)))


def _delta_word(d: float) -> str:
    """Plain-language size of a Cliff's delta, by the conventional thresholds."""
    a = abs(d)
    return "negligible" if a < 0.15 else "small" if a < 0.33 else "medium" if a < 0.474 else "large"


def _failure_windows(bad: "pd.Series", ns_per_frame: float = 0.0, gap_ns: float = 5.0) -> str:
    """Group failed frames into the time windows they occupy.

    Whether the failures cluster or scatter is diagnostic: a cluster means the trajectory itself
    went bad over that stretch (worth inspecting), whereas an even scatter is ordinary Prime
    convergence noise. Frames closer than `gap_ns` are treated as one window.
    """
    if bad.empty or ns_per_frame <= 0:
        return ""
    t = sorted(float(f) * ns_per_frame for f in bad.index)
    wins, start, prev = [], t[0], t[0]
    for x in t[1:]:
        if x - prev > gap_ns:
            wins.append((start, prev))
            start = x
        prev = x
    wins.append((start, prev))
    return ", ".join(f"{a:.0f} ns" if b - a < 1 else f"{a:.0f}–{b:.0f} ns" for a, b in wins)


def _short_ligand(name) -> str:
    """The ligand's short name (FA / DFA / TFA) from CFG — the same abbreviations every figure uses.

    A ligand absent from the map keeps its full name rather than being silently mangled.
    """
    _n = re.sub(r"^\d+_", "", str(name or "")).replace("_", " ").strip()
    return CFG.VIS_LIGAND_SHORT.get(_n.lower(), _n)


def plot_mmgbsa_individual(out_dir: Path, job_name: str, rank: str, dg: "pd.Series",
                           ns_per_frame: float = 0.0) -> None:
    """Per-job MM-GBSA: ΔG_bind against simulation time, plus its distribution.

    Two things this figure must get right:
      • The x-axis is the SIMULATION, not the row number. The series is indexed by trajectory
        frame, so a strided run still spans the full trajectory (with gaps), and 'when' an
        event happened is readable.
      • A handful of failed minimisations (ΔG ≈ −1000 kcal/mol) must not set the y-scale, or
        the 99%+ of frames that carry the actual signal collapse onto a flat line. The axis is
        therefore scaled to the robust core and the failures are reported explicitly at the
        panel edge — flagged, never silently dropped.
    """
    if dg.empty:
        return
    out_dir.mkdir(parents=True, exist_ok=True)

    _bad = _dg_failures(dg)
    _core = dg.drop(_bad.index)
    _x = dg.index.to_numpy(dtype=float)
    if ns_per_frame > 0:
        _x = _x * ns_per_frame
        _xlabel = "Simulation time (ns)"
    else:
        _xlabel = "Trajectory frame"

    """
    Y-window from the CORE ensemble (failures already excluded), not from the full series —
    otherwise a handful of near-failures at −100 still stretch the axis and squash the band
    where 99% of the data lives. The failures are not hidden: they are drawn at the floor and
    counted in the legend.
    """
    _clip = float(getattr(CFG, "MMGBSA_PLOT_CLIP_PCT", 0.5))
    _ref = _core if len(_core) > 10 else dg
    _lo, _hi = _ref.quantile(_clip / 100.0), _ref.quantile(1 - _clip / 100.0)
    _pad = 0.15 * max(_hi - _lo, 1.0)
    _ylo, _yhi = _lo - _pad, _hi + _pad

    _med = float(dg.median())
    _central = _avg_dg(dg)
    _est_label = str(getattr(CFG, "MMGBSA_AVERAGING", "mean")).lower()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.2),
                                   gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05},
                                   sharey=True)

    # ── Panel 1: ΔG_bind over the simulation ──────────────────────────────────────────────
    # Plot the core points only; the failures appear once, as markers at the floor, so a blown-up
    # frame is never shown twice (as a stray dot AND as a flagged failure).
    _cx = _core.index.to_numpy(dtype=float) * (ns_per_frame if ns_per_frame > 0 else 1.0)
    ax1.scatter(_cx, _core.values, s=3, color=_INK["series"], alpha=0.28, linewidths=0, zorder=2)
    """
    Rolling trend on the CORE ensemble. Computing it on the full series lets a single blown-up
    frame drag the window's lower quartile hundreds of kcal/mol, flaring the band downward
    wherever a failure sits — an artefact of the failure, not a change in binding.
    """
    _w = max(5, len(_core) // 50)
    _roll = _core.rolling(_w, min_periods=max(1, _w // 4), center=True)
    ax1.fill_between(_cx, _roll.quantile(0.25).values, _roll.quantile(0.75).values,
                     color=_INK["series"], alpha=0.22, linewidth=0, zorder=3,
                     label="rolling interquartile range")
    ax1.plot(_cx, _roll.median().values, color=_INK["dark"], linewidth=2.0, zorder=4,
             label="rolling median")
    ax1.axhline(_med, color=_INK["mean"], linestyle="--", linewidth=1.6, zorder=5,
                label=f"median = {_med:.1f} kcal/mol")
    if abs(_central - _med) > 0.05:
        ax1.axhline(_central, color=_INK["ghost"], linestyle=":", linewidth=1.4, zorder=5,
                    label=f"{_est_label} = {_central:.1f} (outlier-pulled)")
    if len(_bad):
        _bx = _bad.index.to_numpy(dtype=float) * (ns_per_frame if ns_per_frame > 0 else 1.0)
        ax1.scatter(_bx, np.full(len(_bad), _ylo), marker="v", s=34, color=_INK["warn"],
                    clip_on=False, zorder=6,
                    label=f"{len(_bad)} outlier frames (to {_bad.min():.0f})")
    ax1.set_xlabel(_xlabel)
    ax1.set_ylabel("ΔG$_{bind}$ (kcal/mol)  ·  Prime MM-GBSA")
    ax1.set_ylim(_ylo, _yhi)
    """
    Span the WHOLE simulation, not merely the last sampled frame: with a stride the final frame
    falls short of the run's end (999.8 of 1000 ns), which would clip the axis just before its
    round final tick. Round the upper limit up so the axis reads 0 → the true trajectory length.
    """
    _xmax = float(np.ceil(_x.max() / 100.0) * 100) if ns_per_frame > 0 else float(_x.max())
    ax1.set_xlim(0, _xmax)
    # Compact 3-column legend inside the panel: tight column/handle spacing keeps it narrow
    # enough to sit in the sparse lower-left corner without covering the data.
    ax1.legend(frameon=True, fontsize=8, ncol=3, loc="lower left", framealpha=0.92,
               columnspacing=0.9, handlelength=1.6, handletextpad=0.5, borderpad=0.4)
    ax1.grid(alpha=0.25, linewidth=0.5)

    # ── Panel 2: distribution of the core ensemble (shares the y-axis) ────────────────────
    ax2.hist(_core.values, bins=45, range=(_ylo, _yhi), color=_INK["hist"], alpha=0.85,
             orientation="horizontal", zorder=2)
    ax2.axhline(_med, color=_INK["mean"], linestyle="--", linewidth=1.6, zorder=3)
    ax2.axhspan(dg.quantile(0.25), dg.quantile(0.75), color=_INK["mean"], alpha=0.10, zorder=1)
    ax2.set_xlabel(f"Frames scored\n(n = {len(dg):,}; IQR shaded)")
    # Same y gridlines and tick marks as the time panel (labels suppressed, since the axis is
    # shared) — so a ΔG value can be read straight across from one panel to the other.
    ax2.tick_params(labelleft=False, left=True)
    ax2.grid(alpha=0.25, linewidth=0.5)

    _rk = f"{int(rank):02d}" if str(rank).isdigit() else str(rank)
    out_path = out_dir / f"Rank_{_rk}_MMGBSA_Profile_{job_name}.png"
    plt.savefig(out_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    _echo(f"    ✔ Figure       : {out_path.name}")


def plot_mmgbsa_combined(out_dir: Path, per_job: list, tiers: dict,
                         ligands: "dict | None" = None, nspf: "dict | None" = None) -> None:
    """Compare ΔG_bind across the ranks, three ways.

      left   — the full distribution each summary is drawn from (violin + IQR box), with the
               mean tracked alongside the median so the gap between them exposes which ensembles
               are skewed by failed minimisations.
      middle — the cumulative distributions, which show HOW the ensembles differ (a shift in
               location versus a difference in spread) rather than only THAT they differ.
      right  — the ranking, as the median with a moving-block bootstrap confidence interval,
               plus Cliff's delta against the tightest binder.

    Statistics note: MD frames are autocorrelated, so a t-test or Mann-Whitney over ~10⁴ frames
    would report p ≈ 0 for any difference at all and mean nothing. The bootstrap resamples
    contiguous blocks (CFG.MMGBSA_BOOTSTRAP_BLOCK) so each block counts as one independent draw.
    Colours come from CFG.MMGBSA_RANK_PALETTE.
    """
    per_job = [(r, dg) for r, dg in per_job if not dg.empty]
    if not per_job:
        _echo("  [!] MM-GBSA combined skipped — no parsed ΔG_bind series.")
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    per_job.sort(key=lambda x: int(x[0]))
    ligands, nspf = ligands or {}, nspf or {}

    _pal    = list(CFG.MMGBSA_RANK_PALETTE)
    _cols   = [_pal[i % len(_pal)] for i in range(len(per_job))]
    _cores  = [dg.drop(_dg_failures(dg).index) for _, dg in per_job]
    _nbad   = [len(_dg_failures(dg)) for _, dg in per_job]
    _meds   = [float(dg.median()) for _, dg in per_job]
    _means  = [float(dg.mean())   for _, dg in per_job]
    _q1     = [float(dg.quantile(.25)) for _, dg in per_job]
    _q3     = [float(dg.quantile(.75)) for _, dg in per_job]
    _cis    = [_block_bootstrap_median_ci(dg) for _, dg in per_job]
    _x      = np.arange(len(per_job))

    def _label(r):
        lig = _short_ligand(ligands.get(int(r), ""))
        return f"#{r}\n{lig}" if lig else f"#{r}"
    _labels = [_label(r) for r, _ in per_job]

    _clip = float(getattr(CFG, "MMGBSA_PLOT_CLIP_PCT", 0.5))
    _pool = pd.concat(_cores)
    _lo, _hi = _pool.quantile(_clip / 100.0), _pool.quantile(1 - _clip / 100.0)
    _pad = 0.14 * max(_hi - _lo, 1.0)

    fig = plt.figure(figsize=(17, 13))
    # The top panels carry their definitions in their own legends, so the row gap only needs to
    # clear the tick labels.
    gs  = fig.add_gridspec(2, 2, height_ratios=[1, 1.1], hspace=0.22, wspace=0.22)
    ax1 = fig.add_subplot(gs[0, 0])     # distribution
    ax3 = fig.add_subplot(gs[0, 1])     # ranking
    ax4 = fig.add_subplot(gs[1, :])     # time + cumulative, spanning the row

    # ── Panel 1: distribution, with median vs mean ───────────────────────────────────────
    parts = ax1.violinplot([c.values for c in _cores], positions=_x,
                           showmeans=False, showextrema=False, widths=0.85)
    for i, b in enumerate(parts["bodies"]):
        b.set_facecolor(_cols[i]); b.set_alpha(0.45)
        b.set_edgecolor(_cols[i]); b.set_linewidth(1.2)
    # The box is the INTERQUARTILE RANGE: its edges are the 25th/75th percentiles (half of all
    # frames lie inside it), the bar is the median and the whiskers reach 1.5 × IQR. White fill
    # so it reads against the violin.
    ax1.boxplot([c.values for c in _cores], positions=_x, widths=0.18, showfliers=False,
                     patch_artist=True,
                     medianprops=dict(color=_INK["dark"], linewidth=2.2),
                     boxprops=dict(facecolor=_INK["light"], edgecolor=_INK["dark"], linewidth=1.0),
                     whiskerprops=dict(color=_INK["dark"], linewidth=1.0),
                     capprops=dict(color=_INK["dark"], linewidth=1.0))
    ax1.plot(_x, _means, marker="D", markersize=7, color=_INK["warn"], linewidth=1.6,
             linestyle="--", zorder=6)
    ax1.plot(_x, _meds, marker="o", markersize=7, color=_INK["dark"], linewidth=1.6,
             zorder=6)
    """
    The median's value belongs ON the median, not floating above the panel where the reader has to
    match it back by position. It is written beside the marker in the candidate's colour, with a
    white stroke around the glyphs so it stays legible wherever the violin is dense. The outlier
    count sits below the distribution, out of the data's way.
    """
    _stroke = [pe.withStroke(linewidth=2.6, foreground=_INK["light"])]
    for i, (m, n) in enumerate(zip(_meds, _nbad)):
        ax1.annotate(f"{m:.1f}", xy=(_x[i], m), xytext=(_x[i] + 0.30, m),
                     textcoords="data", va="center", ha="left", fontsize=10.5,
                     fontweight="bold", color=_cols[i], zorder=9, path_effects=_stroke)
        if n:
            ax1.text(_x[i], _lo - _pad * 0.55, f"{n} outlier{'s' if n > 1 else ''}",
                     ha="center", fontsize=8, color=_cols[i], fontweight="bold",
                     zorder=9, path_effects=_stroke)
    ax1.set_xticks(_x); ax1.set_xticklabels(_labels)
    ax1.set_ylabel("ΔG$_{bind}$ (kcal/mol)  ·  per-frame Prime MM-GBSA")
    ax1.set_ylim(_lo - _pad, _hi + _pad)
    """
    One legend, inside, lower right — no caption under the panel. The violin needs no entry: it IS
    the distribution and reads as one on sight. What does need saying is what the white box and its
    whiskers mean, so the box appears as a SWATCH rather than as a sentence.
    """
    _hdlA = [
        Patch(facecolor=_INK["light"], edgecolor=_INK["dark"], linewidth=1.2,
              label="box = IQR · whiskers = 1.5 × IQR"),
        Line2D([], [], color=_INK["dark"], marker="o", linewidth=1.6, markersize=7,
               label="median (robust)"),
        Line2D([], [], color=_INK["warn"], marker="D", linestyle="--", linewidth=1.6, markersize=7,
               label="mean (outlier-dragged)"),
    ]
    ax1.legend(handles=_hdlA, loc="lower right", frameon=True, fontsize=8.2, framealpha=0.93,
               ncol=1, borderpad=0.6, labelspacing=0.5, handlelength=1.8, handletextpad=0.6)
    ax1.grid(alpha=0.25, linewidth=0.5, axis="y")

    # ── Panel 3: the ranking, with an honest interval and an effect size ─────────────────
    _err = np.array([[m - a for m, a in zip(_meds, _q1)], [b - m for m, b in zip(_meds, _q3)]])
    ax3.bar(_x, _meds, color=_cols, alpha=0.85, edgecolor=_INK["outline"], linewidth=0.7,
                    yerr=_err, capsize=6, error_kw=dict(ecolor=_INK["muted"], lw=1.0))
    # Bootstrap CI drawn as a thick inner bar: far narrower than the IQR, because it is the
    # uncertainty ON THE MEDIAN, not the spread of the frames.
    for i, (lo_ci, hi_ci) in enumerate(_cis):
        if lo_ci == lo_ci:
            ax3.plot([_x[i], _x[i]], [lo_ci, hi_ci], color=_INK["dark"], linewidth=4.5,
                     solid_capstyle="butt", zorder=5)
    # Label inside each bar, positioned in DATA coordinates: a padding offset is flipped by the
    # inverted y-axis and lands the text outside the bar.
    for i, m in enumerate(_meds):
        ax3.text(_x[i], m * 0.5, f"{m:.1f}", ha="center", va="center",
                 fontsize=11.5, fontweight="bold", color=_INK["dark"], zorder=9,
                 path_effects=[pe.withStroke(linewidth=3.0, foreground=_INK["light"])])
    # Effect size vs the tightest binder — the number that says whether a gap MATTERS.
    """
    Effect size, reported PAIRWISE. There is no control or reference candidate here — these are
    three substrates compared against each other — so singling one out as a baseline would
    invent a hierarchy the design does not have. Cliff's delta is computed for every pair
    instead: δ(X vs Y) > 0 means X's frames sit at HIGHER (weaker-binding) ΔG than Y's.
    """
    _pairs = []
    for i in range(len(per_job)):
        for j in range(i + 1, len(per_job)):
            d = _cliffs_delta(per_job[i][1], per_job[j][1])
            _pairs.append(f"#{per_job[i][0]} vs #{per_job[j][0]}: δ = {d:+.2f} ({_delta_word(d)})")
    # The pairwise deltas become the legend's title — inside the panel, never a caption below it.
    # One pair per line keeps the legend narrow enough to sit clear of the tallest bar.
    _delta_title = ("Cliff's δ — δ > 0: the first binds more weakly\n"
                    + "\n".join(_pairs)) if _pairs else None
    ax3.axhline(0.0, color=_INK["outline"], linewidth=0.8)
    ax3.set_xticks(_x); ax3.set_xticklabels(_labels)
    ax3.set_ylabel("Median ΔG$_{bind}$ (kcal/mol)")
    """
    The marks are labelled ON the bars rather than in a legend: the median is already printed inside
    each bar, and the two error bars are named once, in place, on the first bar — a legend for two
    marks that are visible at a glance is wasted panel space. Only the pairwise effect sizes, which
    have nowhere natural to sit, remain as a boxed note.
    """
    _lo_ci0, _hi_ci0 = _cis[0]
    # Named in place, INSIDE the bar, beside the mark each label refers to — no leader lines to
    # cross a neighbouring bar, and no legend entry for something already visible.
    """
    Each label sits beside the mark it names. They are BLACK on a white stroke, so they read the same
    whether they land on the coloured bar or on the panel behind it — a single colour that depends on
    the background would vanish on one of the two.
    """
    _stroke_b = [pe.withStroke(linewidth=2.6, foreground=_INK["light"])]
    ax3.text(_x[0] + 0.16, _q1[0], "IQR", fontsize=8.2, color=_INK["dark"], fontweight="bold",
             va="center", ha="left", zorder=9, path_effects=_stroke_b)
    if _lo_ci0 == _lo_ci0:
        ax3.text(_x[0] + 0.16, (_lo_ci0 + _hi_ci0) / 2,
                 f"{getattr(CFG, 'MMGBSA_BOOTSTRAP_CI', 95):.0f}% CI", fontsize=8.2,
                 color=_INK["dark"], fontweight="bold", va="center", ha="left", zorder=9,
                 path_effects=_stroke_b)
    if _delta_title:
        ax3.text(0.03, 0.97, _delta_title, transform=ax3.transAxes, va="top", ha="left",
                 fontsize=7.8, color=_INK["soft"],
                 bbox=dict(boxstyle="round,pad=0.45", facecolor=_INK["light"],
                           edgecolor=_INK["faint"], linewidth=0.8, alpha=0.93))
    ax3.grid(alpha=0.25, linewidth=0.5, axis="y")
    ax3.invert_yaxis()

    # ── Panel 4: time course + cumulative distribution, on one shared ΔG axis ───────────
    _draw_time_cumulative(ax4, per_job, _cols, ligands, nspf)

    """
    Publication finish: label the panels A/B/C (a reader cites 'panel B', not 'the top-right one'),
    drop the top/right spines that carry no data, and give every axes the same tick geometry. The
    time+cumulative panel keeps its top spine — a second x-axis lives there.
    """
    # Tick labels in each candidate's own colour: the label and its violin/bar then read as one
    # object, and the eye does not have to match them through a legend.
    for _ax in (ax1, ax3):
        for _lbl, _c in zip(_ax.get_xticklabels(), _cols):
            _lbl.set_color(_c)
            _lbl.set_fontweight("bold")

    # Uniform spines and tick geometry across the panels (the top+right spines carry no data; the
    # time+cumulative panel keeps its top spine, where its second x-axis lives).
    for _ax in (ax1, ax3, ax4):
        _ax.tick_params(direction="out", length=4.5, width=1.0)
        _ax.spines["right"].set_visible(False)
        if _ax is not ax4:
            _ax.spines["top"].set_visible(False)
        for _sp in ("left", "bottom"):
            _ax.spines[_sp].set_linewidth(1.0)

    out_path = out_dir / "01_MMGBSA_Combined_AllRanks.png"
    plt.savefig(out_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    _echo(f"  MM-GBSA combined figure saved: {out_path.resolve()}")


def _effective_n(dg: "pd.Series") -> float:
    """Effective number of INDEPENDENT samples in an autocorrelated series.

    N_eff = N / (1 + 2·Σρ_k), summing the autocorrelation until it first goes negative
    (the standard initial-positive-sequence cut-off). This is the number that governs how
    precisely ⟨ΔG_bind⟩ is known — not the raw frame count, which merely counts how often the
    same configuration was re-measured.
    """
    v = dg.to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    n = len(v)
    if n < 10:
        return float(n)
    v = v - v.mean()
    denom = float((v * v).sum())
    if denom <= 0:
        return float(n)
    ac = np.correlate(v, v, mode="full")[n - 1:] / denom
    tau = 0.0
    for k in range(1, min(n // 2, 2000)):
        if ac[k] <= 0:
            break
        tau += ac[k]
    return float(n / (1.0 + 2.0 * tau))


def _draw_time_cumulative(axT, per_job: list, cols: list, ligands: dict, nspf: dict) -> None:
    """Draw the time course and the cumulative distribution into ONE axes, on a shared ΔG axis.

      • BOTTOM x-axis — simulation time, in CFG.MMGBSA_TIME_WINDOW_NS windows. Solid line with
        markers = window median; band = window IQR; pale band = the candidate's overall IQR;
        ▼ = a window whose minimisations blew up.
      • TOP x-axis — cumulative fraction of frames (0→1). Dashed line = the candidate's cumulative
        curve, against the SAME ΔG axis. Where it is steep, that ΔG level is heavily populated —
        and it should be the level the solid time course keeps returning to.

    Reading the two together answers what neither can alone: whether a candidate's typical binding
    energy is a stable property (steep curve + flat time course) or an average over states it keeps
    moving between (stretched curve + wandering time course). One ΔG axis, labelled once; the two
    x-axes carry different quantities and are named in a single label, the top one in brackets.
    """
    _gl = str(CFG.MMGBSA_GRID_LEFT)
    _gr = str(CFG.MMGBSA_GRID_RIGHT)
    _win = float(getattr(CFG, "MMGBSA_TIME_WINDOW_NS", 100.0))

    cores = [dg.drop(_dg_failures(dg).index) for _, dg in per_job]
    meds  = [float(dg.median()) for _, dg in per_job]
    q1s_o = [float(dg.quantile(.25)) for _, dg in per_job]
    q3s_o = [float(dg.quantile(.75)) for _, dg in per_job]
    cis   = [_block_bootstrap_median_ci(dg) for _, dg in per_job]

    def _label(r):
        return f"#{r} {_short_ligand(ligands.get(int(r), ''))}".strip()
    flat = [_label(r) for r, _ in per_job]

    _clip = float(getattr(CFG, "MMGBSA_PLOT_CLIP_PCT", 0.5))
    pool = pd.concat(cores)
    lo, hi = float(pool.quantile(_clip / 100.0)), float(pool.quantile(1 - _clip / 100.0))
    pad = 0.14 * max(hi - lo, 1.0)

    axC = axT.twiny()          # same plot area, same ΔG; its own x = cumulative fraction
    summary, tmax = [], 0.0
    for i, (r, dg) in enumerate(per_job):
        ns, core = float(nspf.get(int(r), 0.0)), cores[i]
        if ns > 0:
            tc = core.index.to_numpy(dtype=float) * ns
            tmax = max(tmax, float(tc.max()))
            edges = np.arange(0, tc.max() + _win, _win)
            axT.axhspan(q1s_o[i], q3s_o[i], color=cols[i], alpha=0.06, linewidth=0, zorder=1)
            cents, wm, wq1, wq3 = [], [], [], []
            for k in range(len(edges) - 1):
                sel = core.to_numpy(dtype=float)[(tc >= edges[k]) & (tc < edges[k + 1])]
                if sel.size < 5:
                    continue
                cents.append((edges[k] + edges[k + 1]) / 2.0)
                wm.append(float(np.median(sel)))
                wq1.append(float(np.percentile(sel, 25)))
                wq3.append(float(np.percentile(sel, 75)))
            if cents:
                axT.fill_between(cents, wq1, wq3, color=cols[i], alpha=0.20, linewidth=0, zorder=2)
                axT.plot(cents, wm, color=cols[i], linewidth=2.4, marker="o", markersize=5.5,
                         markerfacecolor=_INK["light"], markeredgewidth=1.6, zorder=5)
            bad = _dg_failures(dg)
            if len(bad):
                tb = bad.index.to_numpy(dtype=float) * ns
                for k in range(len(edges) - 1):
                    if ((tb >= edges[k]) & (tb < edges[k + 1])).any():
                        axT.plot((edges[k] + edges[k + 1]) / 2.0, hi + pad * 0.16, marker="v",
                                 markersize=8, color=cols[i], zorder=6, clip_on=False)
        v = np.sort(core.to_numpy(dtype=float))
        frac = np.arange(1, len(v) + 1) / len(v)
        axC.plot(frac, v, color=cols[i], linewidth=2.0, linestyle="--", alpha=0.75, zorder=4)
        axT.axhline(meds[i], color=cols[i], linewidth=0.9, linestyle=":", alpha=0.7, zorder=3)
        summary.append(f"{flat[i]}  ·  median {meds[i]:.1f} [{cis[i][0]:.1f}, {cis[i][1]:.1f}]"
                       f"  ·  N$_{{eff}}$ ≈ {_effective_n(dg):,.0f}")

    tend = float(np.ceil(tmax / _win) * _win) if tmax > 0 else 1.0
    axT.set_ylim(lo - pad, hi + pad)
    axT.set_xlim(0, tend)
    axT.set_xticks(np.arange(0, tend + 1, max(_win, tend / 10)))
    axT.set_xlabel(f"Simulation time (ns), in {_win:.0f} ns windows — solid line, markers"
                   f"        [top axis: cumulative fraction of frames — dashed line]")
    axT.set_ylabel("ΔG$_{bind}$ (kcal/mol)  ·  per-frame Prime MM-GBSA")
    axT.grid(color=_gl, alpha=0.25, linewidth=0.6)
    axT.set_axisbelow(True)
    axC.set_xlim(0, 1)
    axC.tick_params(axis="x", colors=_gr)
    axC.grid(axis="x", color=_gr, alpha=0.30, linewidth=0.8, linestyle="--")
    axC.set_axisbelow(True)

    hdl = [Patch(facecolor=cols[i], alpha=0.75, edgecolor=cols[i], label=summary[i])
           for i in range(len(per_job))]
    hdl += [
        Line2D([], [], color=_INK["muted"], linewidth=2.4, marker="o", markerfacecolor=_INK["light"],
               markersize=7,
               label="window median   ·   band = window IQR   ·   pale band = overall IQR"),
        Line2D([], [], color=_INK["muted"], linewidth=2.0, linestyle="--",
               label="cumulative curve (top axis)   ·   ⋯ = overall median"),
        Line2D([], [], marker="v", color=_INK["muted"], linestyle="none", markersize=8,
               label="window containing failed minimisations"),
    ]
    # Column count derived from the entry count, so the legend stays balanced as candidates change.
    _rows = max(1, int(getattr(CFG, "MMGBSA_LEGEND_MAX_ROWS", 3)))
    axT.legend(handles=hdl, loc="lower right", frameon=True, fontsize=8.5, borderpad=0.7,
               ncol=max(1, math.ceil(len(hdl) / _rows)), columnspacing=1.1, handlelength=1.8,
               handletextpad=0.5, labelspacing=0.5, framealpha=0.93,
               title="Relative ranking only — GB overstabilises anionic PFAS: compare candidates, "
                     "never absolute values",
               title_fontsize=8.5)


# ── 8.5  Phase driver ────────────────────────────────────────────────────────
def run_mmgbsa_phase(md_dir: Path, run_root: Path, scratch_ok: bool = True) -> str:
    """Run + plot MM-GBSA for every completed MD job (idempotent).

    Returns a status the caller maps to the pipeline step result:
      • "ok"       — every eligible job produced a ΔG_bind CSV (or none eligible);
      • "deferred" — the phase was skipped because scratch could not be guaranteed
                     off /tmp (a live job blocked relocation) — nothing computed;
      • "warn"     — the phase ran but ≥1 job's MM-GBSA failed (e.g. Prime rc=1).
    "deferred"/"warn" let the pipeline report WARN (not a false PASS) without
    hard-aborting: MM-GBSA is complementary to the QSite barrier, so Step 07 still
    runs. "disabled" (CFG.MMGBSA_RUN=False) is an intentional no-op → "ok"."""
    if not getattr(CFG, "MMGBSA_RUN", False):
        _echo("  MM-GBSA disabled (CFG.MMGBSA_RUN = False) — skipped.")
        return "ok"
    if not scratch_ok:
        # The job-server scratch dir is still on the OS disk (relocation deferred —
        # typically because another job is running). A 100k-frame Prime run would
        # copy hundreds of GB there and die after hours; skip it rather than burn
        # the time. SID above has already run. Re-run this step once the disk is
        # relocated (no live jobs) to compute MM-GBSA.
        _echo("  MM-GBSA SKIPPED — Schrödinger job-server scratch is not on the "
              "working disk (see the [job-server] note above). Re-run 06 once the "
              "server dir is relocated and no other job is running.")
        return "deferred"
    job_dirs = sorted(
        (d for d in md_dir.iterdir()
         if d.is_dir() and re.match(r"desmond_md_job_R(?:ank)?_\d", d.name)
         and (d / f"{d.name}-out.cms").is_file()),
        key=_natural_rank,
    )
    if not job_dirs:
        _echo("  No completed MD jobs (-out.cms) for MM-GBSA. Skipping.")
        return "ok"
    _echo(_SEP)
    _echo(f"Prime MM-GBSA — rescoring binding free energy for {len(job_dirs)} completed MD job(s)")
    _echo(_SEP)
    # Figures folder derived from the resolved MD dir (not a hardcoded path).
    out_dir = md_dir / getattr(CFG, "MMGBSA_OUTPUT_SUBDIR", "Prime_MMGBSA")
    tiers = _lookup_tiers(run_root)
    ligands = _lookup_ligands(run_root)
    nspf: dict = {}   # rank → ns per trajectory frame, for the time-resolved panel
    per_job = []
    summary_rows = []
    _failed = 0
    for _i, d in enumerate(job_dirs, start=1):
        job_name = d.name
        rank = _rank_of(job_name)
        _echo("")
        _echo(f"  {_RULE}")
        _echo(f"  Rank {rank}  ·  {job_name}   ({_i} of {len(job_dirs)})")
        _echo(f"  {_RULE}")
        csv = run_mmgbsa(d, job_name, rank)
        if csv is None:
            _failed += 1
            _echo(f"    ✘ MM-GBSA did not complete for Rank {rank} — see the messages above.")
            continue
        dg = _mmgbsa_dg_series(csv)
        _bad = _dg_failures(dg)
        # Frame index → ns, read from the trajectory itself (no assumed run length).
        _span = traj_span_ns(d / f"{job_name}_trj")
        _nfr = traj_frame_count(d / f"{job_name}_trj")
        _nspf = (_span / (_nfr - 1)) if (_span > 0 and _nfr > 1) else 0.0
        if not dg.empty:
            _echo(f"    ΔG_bind        : median {float(dg.median()):.2f} kcal/mol "
                  f"(IQR {float(dg.quantile(.25)):.1f} to {float(dg.quantile(.75)):.1f}) "
                  f"over {dg.size:,} scored frames")
            if len(_bad):
                """
                Only the strongly-negative tail is a blown-up minimisation. Outliers on the WEAK
                side (a frame where the ligand had drifted out of the pocket) are ordinary
                unbinding events, not numerical failures — calling them 'blown-up' and quoting
                their minimum would misdescribe them.
                """
                _q1, _q3 = dg.quantile(.25), dg.quantile(.75)
                _lo_tail = _bad[_bad < _q1]
                _hi_tail = _bad[_bad > _q3]
                _parts = []
                if len(_lo_tail):
                    _parts.append(f"{len(_lo_tail)} unphysically strong (to {float(_lo_tail.min()):.0f} "
                                  f"kcal/mol — blown-up minimisations)")
                if len(_hi_tail):
                    _parts.append(f"{len(_hi_tail)} weakly bound (to {float(_hi_tail.max()):.0f} "
                                  f"kcal/mol — ligand loosely held, not a numerical failure)")
                _echo(f"    Outlier frames : {len(_bad)} of {dg.size:,} "
                      f"({100 * len(_bad) / dg.size:.2f}%) — {'; '.join(_parts)}. "
                      f"They shift the arithmetic mean ({float(dg.mean()):.2f}), so the median "
                      f"above is the honest estimate.")
                _wins = _failure_windows(_bad, _nspf)
                if _wins:
                    _echo(f"    ↳ when         : {_wins}. Failures that CLUSTER in time point at "
                          f"unstable stretches of the trajectory (inspect those frames); failures "
                          f"scattered evenly are ordinary Prime convergence noise.")
            _echo(f"    Interpretation : GB implicit solvent overstabilises anionic PFAS — compare "
                  f"ΔG_bind BETWEEN ranks, never as an absolute affinity.")
        try:
            plot_mmgbsa_individual(out_dir, job_name, rank, dg, ns_per_frame=_nspf)
        except Exception as e:
            _echo(f"    ✘ per-job plot failed ({e}) — skipped.")
        per_job.append((rank, dg))
        if str(rank).isdigit():
            nspf[int(rank)] = _nspf
        if not dg.empty:
            summary_rows.append({
                "Scientific_Rank": rank, "Job_Name": job_name,
                "degrader_tier": tiers.get(int(rank), "") if str(rank).isdigit() else "",
                "N_Frames": int(dg.size),
                "N_Failed_Minimisations": int(len(_bad)),
                "MMGBSA_dG_Boltzmann_kcal": round(_boltzmann_mean_dg(dg), 2),
                "MMGBSA_dG_ArithMean_kcal": round(float(dg.mean()), 2),
                "MMGBSA_dG_Median_kcal": round(float(dg.median()), 2),
                "MMGBSA_dG_Min_kcal":  round(float(dg.min()), 2),
                "MMGBSA_dG_Max_kcal":  round(float(dg.max()), 2),
                "MMGBSA_dG_Std_kcal":  round(float(dg.std()), 2),
                "Source_CSV": csv.name,
            })
    # Write the consolidated data table alongside the figures.
    if summary_rows:
        out_dir.mkdir(parents=True, exist_ok=True)
        _sdf = pd.DataFrame(summary_rows)
        if "Scientific_Rank" in _sdf.columns:
            _sdf = _sdf.sort_values("Scientific_Rank",
                                    key=lambda s: pd.to_numeric(s, errors="coerce"))
        _csv_out = out_dir / CFG.FILE_MMGBSA_SUMMARY
        _sdf.to_csv(_csv_out, index=False)
        _echo(f"  MM-GBSA summary table saved : {_csv_out.resolve()}")
    try:
        plot_mmgbsa_combined(out_dir, per_job, tiers, ligands, nspf)
    except Exception as e:
        _echo(f"  [!] MM-GBSA combined plot failed ({e}) — skipped.")
    if _failed:
        _echo(f"  [!] MM-GBSA: {_failed} of {len(job_dirs)} job(s) failed "
              f"(see per-rank rc/cause above) — step will report WARN, not PASS.")
        return "warn"
    return "ok"


# =============================================================================
# SECTION 9: END-OF-RUN HOUSEKEEPING (Schrödinger scratch + job-server home)
# =============================================================================
def _jsc_jobs() -> "tuple[list[str], int]":
    """(completed job IDs, number still alive). Alive = Running/Waiting/anything not finished.

    Reads `jsc list`; an unreachable server yields (no ids, 0 alive) so callers treat it as
    'nothing of ours is running' rather than crashing.
    """
    jsc = Path(SCHRODINGER) / "jsc"
    if not jsc.is_file():
        return [], 0
    try:
        r = subprocess.run([str(jsc), "list"], capture_output=True, text=True, timeout=60)
    except Exception:
        return [], 0
    """
    `jsc list` is single-space aligned and its Status may itself contain a space ("1% done"),
    so neither a two-space split nor a plain field split is safe. Anchor on the Created date
    (Mon-DD) and take everything between the job name and that date as the status:

        ee9d84d8 desmond_md_job_1IVO_extend Running  Jul-13 01:46 3h24m shark
        1c987b28 ..._mmgbsa_shard002-prime  1% done  Jul-12 11:07 7m41s shark
    """
    done, alive = [], 0
    for line in (r.stdout or "").splitlines():
        m = re.match(r"^([0-9a-f]{8})\s+(\S+)\s+(.*?)\s+[A-Z][a-z]{2}-\d{2}\s", line)
        if not m:
            continue
        jid, status = m.group(1), m.group(3).strip().lower()
        if re.search(r"running|waiting|launched|submitted|% done", status):
            alive += 1
        elif re.search(r"completed|finished|died|killed|stopped|incorporated", status):
            done.append(jid)
    return done, alive


def cleanup_schrodinger_dirs(scratch: Path) -> None:
    """End-of-run housekeeping for the two Schrödinger directories.

    Only ever runs when NOTHING is left alive on the job server — a live job still writes into
    both of these, and removing them under it is exactly what produces 'Error locating localhost
    job server config' (and kills every subsequent submission). A settle window first lets the
    server flush and copy outputs back to the working folders before anything is deleted.
    """
    _cool = int(getattr(CFG, "SCHRODINGER_SCRATCH_COOLDOWN_SEC", 60))
    _echo("")
    _echo(f"Housekeeping — waiting {_cool}s for the job server to flush its outputs before cleaning up.")
    time.sleep(_cool)

    _done, _alive = _jsc_jobs()
    if _alive:
        _echo(f"  ⚠ {_alive} job(s) still on the server — nothing removed. Both Schrödinger "
              f"directories are left intact (deleting them under a live job breaks it).")
        return

    # 1. Client scratch: pure transient working space. Name-guarded so only the dedicated
    #    folder can ever be removed, never an arbitrary path.
    if (scratch.name == getattr(CFG, "SCHRODINGER_SCRATCH_SUBDIR", "_Schrodinger_Scratch")
            and scratch.exists()):
        shutil.rmtree(scratch, ignore_errors=True)
        _echo(f"  ✔ Removed the client scratch directory ({scratch.name}) — transient working space.")

    # 2. Job-server home: NOT scratch. It holds the daemon's jobdb, filestore, logs and binaries,
    #    and it must survive. What grows is the per-job filestore/logs, so delete the COMPLETED
    #    jobs from the server and leave the daemon running.
    _js = _SCRIPT_DIR / getattr(CFG, "SCHRODINGER_JOBSERVER_SUBDIR", "_Schrodinger_JobServer")
    if getattr(CFG, "SCHRODINGER_JOBSERVER_PURGE_COMPLETED", True) and _done:
        try:
            subprocess.run([str(Path(SCHRODINGER) / "jsc"), "delete", *_done],
                           capture_output=True, text=True, timeout=300)
            _echo(f"  ✔ Purged {len(_done)} completed job(s) from the job server — frees its "
                  f"filestore and logs, daemon stays up.")
        except Exception as e:
            _echo(f"  ⚠ Could not purge completed jobs ({e}) — {_js.name} keeps their files.")

    """
    Removing the job server's home entirely is only safe with the server STOPPED and idle; the
    next run recreates it (at the cost of a server restart). Off by default — the directory is
    the daemon's install, not junk.
    """
    if getattr(CFG, "SCHRODINGER_JOBSERVER_REMOVE_WHEN_IDLE", False) and _js.exists():
        try:
            subprocess.run([str(Path(SCHRODINGER) / "jsc"), "local-server-stop"],
                           capture_output=True, text=True, timeout=120)
            time.sleep(3)
            shutil.rmtree(_js, ignore_errors=True)
            _echo(f"  ✔ Stopped the idle job server and removed {_js.name}; the next run recreates it.")
        except Exception as e:
            _echo(f"  ⚠ Could not remove {_js.name} ({e}) — left in place.")
    elif _js.exists():
        try:
            _sz = sum(f.stat().st_size for f in _js.rglob("*") if f.is_file()) / 2**30
            _echo(f"  ℹ Kept {_js.name} ({_sz:.1f} GB) — this is the job server's home (jobdb, "
                  f"filestore, logs, binaries), not scratch. Deleting it under a live daemon "
                  f"breaks every submission. Set CFG.SCHRODINGER_JOBSERVER_REMOVE_WHEN_IDLE = True "
                  f"to have it stopped and removed when idle.")
        except Exception:
            pass


# =============================================================================
# SECTION 10: MAIN
# =============================================================================
def main() -> int:
    # Collapse stacked separator rules: a caller prints a rule, a helper prints its own,
    # and the log grows triple bars with nothing between them.
    _utils_mod.install_console_rule_filter()
    parser = argparse.ArgumentParser(
        description=f"{CFG.PROJECT_NAME} Step 06 — Desmond SID + Prime MM-GBSA Post-Processing")
    parser.add_argument("run_dir", nargs="?", default=None,
                        help="Boltz-2_Run_* directory (auto-detect latest if omitted)")
    parser.add_argument("--pipeline-mode", action="store_true",
                        help="Set when called from 00_00_run_pipeline_FAcDs.sh; "
                             "delegates systemd-oomd masking to the pipeline runner.")
    args = parser.parse_args()

    t0 = time.perf_counter()
    print_script_banner("06_SID_Prime-MMGBSA_FAcDs.py",
                        "Step 06 — Desmond SID + Prime MM-GBSA Post-Processing")

    if not os.path.isfile(SCHROD_RUN):
        _echo(f"ERROR: Schrödinger 'run' not found at {SCHROD_RUN}. "
              f"Set the SCHRODINGER environment variable.")
        return 1

    run_dir = resolve_run_dir(args.run_dir)
    md_dir = _SCRIPT_DIR / run_dir / "6_Physics_Validation" / "MolecularDynamics"
    _open_step_log(md_dir.parent)     # 6_Physics_Validation/00_SID_MMGBSA.log — mirrors this run's output

    if not md_dir.is_dir():
        _echo(f"ERROR: MolecularDynamics directory not found: {md_dir}")
        return 1

    """
    Keep ALL Schrödinger job scratch on the run's own (large) working disk, never
    /tmp on the OS disk. Prime MM-GBSA replicates the multi-GB complexes file into
    every subjob's scratch dir — hundreds of GB on a 100k-frame trajectory — which
    exhausts a small /tmp mid-job (the cause of a silent Prime rc=1 after hours).
    TWO mechanisms:
      1. The authoritative one — set the job server's scratch `tmpdir` (the
         localhost entry of $SCHRODINGER/schrodinger.hosts) to the working disk and
         apply it LIVE with `jsc admin reload-hosts`. This works on a RUNNING server
         (no stop, no killed job), so MM-GBSA still runs when another job happens to be
         active. This is what actually prevents the /tmp overflow.
      2. Belt-and-braces — point SCHRODINGER_TMPDIR/TMPDIR at the working disk for
         any tool that still honours them (driver-side temp, non-server steps).
    """
    _scratch = md_dir / getattr(CFG, "SCHRODINGER_SCRATCH_SUBDIR", "_Schrodinger_Scratch")
    _scratch.mkdir(parents=True, exist_ok=True)
    # Export the working-disk scratch BEFORE ensure_jobserver: if it (re)starts the job-server
    # daemon (only ever when idle), the daemon inherits these and runs every Prime subjob under
    # the USB scratch instead of /tmp. The daemon's own TMPDIR — not local-server-dir and not a
    # client-set var on a live daemon — is what governs where subjob working dirs are created.
    os.environ["SCHRODINGER_TMPDIR"] = str(_scratch)
    os.environ["TMPDIR"] = str(_scratch)
    _scratch_ok = _utils_mod.ensure_jobserver_on_working_disk(
        SCHRODINGER, _SCRIPT_DIR,
        getattr(CFG, "SCHRODINGER_JOBSERVER_SUBDIR", "_Schrodinger_JobServer"),
        echo=lambda m="": _echo(f"  {m}"))

    _echo(_SEP)
    _echo("Step 06 — Desmond SID analysis, then Prime MM-GBSA rescoring")
    _echo(f"  Run directory        : {run_dir}")
    _echo(f"  MD trajectories      : {md_dir}")
    _echo(f"  Schrödinger scratch  : {_scratch}")
    _echo("                         (on the run's working disk — never /tmp, so a large "
          "Prime job cannot fill the OS disk)")
    _echo(_SEP)

    job_dirs = sorted(
        (d for d in md_dir.iterdir()
         if d.is_dir() and re.match(r"desmond_md_job_R(?:ank)?_\d", d.name)),
        key=_natural_rank,
    )
    if not job_dirs:
        _echo("SKIP — no desmond_md_job_R_* directories in MolecularDynamics yet; "
              "MD not run for this cohort. Nothing to post-process (SID + MM-GBSA skipped).")
        return 0

    to_run = scan_jobs(job_dirs, md_dir)

    """
    Mask oomd (standalone only); the context manager restores it on exit. Both
    SID and MM-GBSA run inside the guard. MM-GBSA runs even when no SID work is
    pending, since it only needs completed MD jobs (-out.cms).
    """
    run_root = _SCRIPT_DIR / run_dir
    with OomdGuard(active=not args.pipeline_mode):
        if to_run:
            process_jobs(to_run)
        else:
            _echo("SID analysis: nothing to do — every MD job is either already analysed or "
                  "not finished simulating. Moving on to Prime MM-GBSA.")
        _mmgbsa_status = run_mmgbsa_phase(md_dir, run_root, scratch_ok=_scratch_ok)

    _echo("")
    cleanup_schrodinger_dirs(_scratch)
    _echo(_SEP)
    # Report the true outcome. SID always ran; MM-GBSA is complementary, so a
    # deferred/failed MM-GBSA is a WARNING (distinct exit EXIT_WARN), not a hard
    # error — the pipeline runner renders WARN and continues to Step 07 rather than
    # printing a false PASS or aborting.
    if _mmgbsa_status == "deferred":
        _echo("Desmond SID post-processing complete; MM-GBSA DEFERRED (scratch not on "
              "the working disk — re-run once no job is running and the server dir is "
              "relocated).")
        _echo(_SEP)
        print_elapsed(t0, "06_SID_Prime-MMGBSA_FAcDs.py")
        return EXIT_WARN
    if _mmgbsa_status == "warn":
        _echo("Desmond SID post-processing complete; MM-GBSA completed with FAILURES "
              "(see the per-rank cause above).")
        _echo(_SEP)
        print_elapsed(t0, "06_SID_Prime-MMGBSA_FAcDs.py")
        return EXIT_WARN
    _echo("All available Desmond SID + MM-GBSA post-processing jobs completed successfully.")
    _echo(_SEP)
    print_elapsed(t0, "06_SID_Prime-MMGBSA_FAcDs.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
