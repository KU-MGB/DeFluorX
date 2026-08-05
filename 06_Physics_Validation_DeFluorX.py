#!/usr/bin/env python3
"""
===============================================================================
DeFluorX Pipeline  |  Step 06  |  ESP Physics: WaterMap → System Builder → MD → SID → MM-GBSA → Defluorination
===============================================================================
Builds and runs the full explicit-solvent physics for every MD-selected complex,
using the Jaguar ESP partial charges on the ligand so the reactive α-carbon
carries its true electrophilicity (the property SN2 defluorination depends on),
then post-processes each trajectory - all in one pass, one merged log.

Per run, in phases (all complexes at each phase before the next):
  1. import   - copy the prepared complex (05's R{N}_<stem>.pdb handover) as-is and
                write the ESP charges onto it (→ 01_Prepared_Proteins, 02_ESP_Charged_Complexes).
  2. WaterMap - hydration-site thermodynamics around the ligand, holo (→ 03_WaterMaps).
                Each WaterMap is tried up to 3× (GCMC is stochastic); after 3 it is
                SKIPPED (red) and the run continues - a missing WaterMap never fails its complex.
  3. build    - Desmond System Builder: minimise-volume, orthorhombic TIP3P box (10 Å
                buffer, OPLS4), auto-neutralise + 0.15 M NaCl, then write the ESP charges
                into the built .cms force field and HARD-VERIFY (→ 04_System_Builder).
  4. MD       - Desmond MD (relax + NPT production) → 05_MD_Simulations, pipelined GPU→CPU per rank.
                The ligand + backbone pose restraint spans the FINAL relaxation stage and all of
                production (not production alone): the stock relaxation's last stage runs unrestrained,
                which would let the substrate relax out of the near-attack docked pose before
                production's restraint engages, so it is injected there too (see _md_msj).
                Each rank's MD runs on the GPU; the instant it lands and its files settle, that rank's
                Extraction (unpack the production _trj/.ene) → SID (event_analysis + analyze_simulation
                → *_SID-out.eaf) → Prime MM-GBSA (thermal_mmgbsa → per-frame ΔG_bind) → defluorination
                geometry (SN2 attack pose + NAC + fluoride cradle + carboxylate clamp + MM-GBSA drivers,
                read natively from the cms + _trj) is queued to a single CPU worker, and the GPU starts
                the next rank's MD immediately. The worker drains one rank at a time, so exactly one
                Prime batch touches the scratch disk at once (no disk contention) while the GPU never
                idles. A rank's OWN figures (MM-GBSA profile + defluorination geometry) are drawn the
                moment that rank's post-processing lands, so each rank is readable while the later
                ranks are still on the GPU; only the cross-rank merged figures wait for the final
                pass, which needs every rank present. All figures land in 06_Analysis.

MM-GBSA (end-state binding ΔG over the ensemble) is complementary to the QSite QM/MM
reaction barrier (Step 07): it scores BINDING, not C–F cleavage. The defluorination step adds
the reactive geometry (does the substrate reach the in-line attack pose) alongside binding -
both are necessary for turnover; Step 07 delivers the QM/MM verdict.

Uses the central CFG / ProjectUtils modules. The ESP/build/WaterMap/MD stage bodies
import `schrodinger` in-process, so the script runs under the Schrödinger Python; it
may be launched either as `$SCHRODINGER/run 06_...py` or as a plain `python 06_...py`
(project conda env) - in the latter case it transparently re-execs under $SCHRODINGER/run.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 05 August 2026 <─────────────────────────────────────────────────────────
===============================================================================
Usage:
  python 06_Physics_Validation_DeFluorX.py [Boltz-2_Run_Directory] [options]

  --test            quick run: WaterMap 2 ns · MD 5 ns/500 frames (~30 min end-to-end)
  --md-ns   N       MD production length (ns)          [default 1000]
  --md-frames N     MD trajectory frames               [default 100000]
  --wm-ns   N       WaterMap production length (ns)     [default 5]
  --lig-dist N      WaterMap active-site radius (Å)     [default 10]
  --stages  a,b,c   subset of {merge,watermap,build,md} [default all]  (Extraction + SID + MM-GBSA + Defluorination run on CPU after each MD)
  --out     DIR     output root                         [default <run>/6_Physics_Validation]
  --pipeline-mode   called from 00_00_run_pipeline_DeFluorX.sh (delegates oomd masking to the runner)

-------------------------------------------------------------------------------
Dependency Map
-------------------------------------------------------------------------------
  Script        : 06_Physics_Validation_DeFluorX.py
  Role          : Step 06 - build + run the ESP-charged explicit-solvent physics
                  (WaterMap, System Builder, MD) and post-process it (SID + MM-GBSA + Defluorination).
  Imports from  : 00_01_Project_Config_DeFluorX.py  (CFG), 00_02_Project_Utils_DeFluorX.py (utils)
  Reads         : <Run>/5_TopN_and_Preparation/3_Comparative_Analysis/
                       06_<tier>_<count>hits_Molecular_Handover_Files/R{N}_<stem>.pdb  (N=Scientific_Rank, SSOT)
                  <Run>/5_TopN_and_Preparation/4_Ligand_ESP_Charges/<stem>_ESP.mae
                  <Run>/1_Boltz2_Production/*Ranked*.csv  (job_name → Mapped_Base)
  Writes        : <out>/01_Prepared_Proteins/R{N}_<stem>.pdb
                  <out>/02_ESP_Charged_Complexes/R_N_<stem>_ESP_Complex.mae
                  <out>/03_WaterMaps/watermap_R_N/*_wm.maegz + *-in.maegz (aligning frame) + watermap_R_N.csv
                  <out>/04_System_Builder/desmond_setup_R_N/desmond_setup_R_N-out.cms
                  <out>/05_MD_Simulations/desmond_md_job_R_N/{-out.cms, _trj/, .ene, *_SID-out.eaf,
                       *_mmgbsa-prime-out.csv (per-frame ΔG_bind + Frame column)}
                  <out>/06_Analysis/{00_MMGBSA_Summary.csv, 01_Physics_Build_Solvation_QC.svg,
                       02_WaterMap_Landscapes_AllRanks.svg, 03_MD_Trajectory_QC.svg,
                       04_MMGBSA_Combined_AllRanks.svg, 05_Defluorination_Combined_AllRanks.svg,
                       Prime-MMGBSA/MMGBSA_Profile_R{N}.svg,
                       Defluorination/Defluorination_R{N}/01_Reactive_Pose_Trajectory.svg … 06_Figure_Descriptions.txt}
                  <out>/00_Physics_Validation.log  (single merged, colour-preserving log; `tail -f` it)
  Upstream      : 05_TopN_and_PDB_Preparation_DeFluorX.py (prepared PDBs + ESP charges).
  Downstream    : 07_MD_QMMM_Defluorination_DeFluorX.py (reads 05_MD_Simulations + 03_WaterMaps).

── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Pipelined GPU→CPU per rank: MD runs on the GPU; the instant it lands (and its files settle)
     that rank's Extraction → SID → MM-GBSA → Defluorination is queued to a single CPU worker while
     the GPU starts the next rank's MD. The worker processes one rank at a time, so two Prime batches
     never share the scratch disk and a running SID/MM-GBSA is never pre-empted, yet the GPU does not
     idle through the hours of CPU post-processing. A later MD that finishes first waits in the queue.
  2. WaterMap ligand is NOT ESP-charged: WaterMap rebuilds ligand charges with its own
     S-OPLS/TIP4P (its GCMC μ_excess is calibrated only for TIP4P). Correct for water
     thermodynamics; ESP lives in the MD/QSite branch where C–F electrophilicity matters.
  3. WaterMap launch: run via `bash -lc` with a detached login-shell env + a /tmp scratch
     (copied back); a job inheriting the $SCHRODINGER/run session, or given an absolute
     input path, cannot stage its GCMC ligand companion and dies at stage 8.
  4. MM-GBSA cost: Prime minimises every scored structure; cost is linear in their number
     (CFG.MMGBSA_STEP_SIZE is the order-of-magnitude knob). Concurrent frame-shards are
     capped by free /tmp (each stages ~22 GB) - the disk-aware subjob cap, not USB relocation.
  5. MM-GBSA validity: GB implicit solvent overstabilises anionic PFAS → ΔG_bind is a
     RELATIVE ranking only, complementary to the QSite barrier. A fraction of a percent of
     frames are failed minimisations (flagged, excluded from the figure scale; MEDIAN reported).
  6. Sudo (OPTIONAL): masking systemd-oomd needs root; standalone offers to prime sudo,
     else continues in NORMAL mode. Under --pipeline-mode the runner owns oomd masking.
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
import csv
import glob
import importlib.util as _ilu
import math
import os
import shutil
import shlex
import re
import stat
import subprocess
import sys
import tempfile
import threading
import queue
import time
from contextlib import contextmanager
from pathlib import Path

# Re-exec under $SCHRODINGER/run when Schrödinger's Python is not the interpreter. The ESP/build/
# WaterMap/MD stage bodies import `schrodinger` in-process, so the script must run under
# $SCHRODINGER/run; this lets it also be launched as a plain `python 06_...py` (project conda env).
try:
    import schrodinger  # noqa: F401
except ModuleNotFoundError:
    _schro = os.environ.get("SCHRODINGER", "/opt/schrodinger")
    os.execv(f"{_schro}/run", [f"{_schro}/run", "python3", os.path.abspath(__file__), *sys.argv[1:]])

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
plt.rcParams["axes.labelpad"] = 8.0   # gap between axis labels and tick values (per-call labelpad still overrides)
import matplotlib.colors as _mcolors  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
# --- consolidated top-level imports (optional/heavy + Schrodinger stay function-local) ---
import atexit
import math as _math
import shutil as _sh
import signal

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
_utils_mod = _load_module("ProjectUtils", _SCRIPT_DIR / "00_02_Project_Utils_DeFluorX.py")
_cfg_mod = _load_module("ProjectConfig", _SCRIPT_DIR / "00_01_Project_Config_DeFluorX.py")

CFG = _cfg_mod.CFG()
print_script_banner = _utils_mod.print_script_banner
print_elapsed = _utils_mod.print_elapsed
apply_figure_style = _utils_mod.apply_figure_style
auto_label_colour = _utils_mod.auto_label_colour
apply_figure_style(CFG)   # one typography definition for every figure the pipeline draws

# Auto-set SCHRODINGER if the env var is absent (Step 07 QM/MM uses the same default).
os.environ.setdefault("SCHRODINGER", "/opt/schrodinger")
SCHRODINGER = os.environ["SCHRODINGER"]
SCHROD_RUN = os.path.join(SCHRODINGER, "run")

_SEP  = _utils_mod.SEPARATOR_HEAVY   # ═×80 - house major boundary
_RULE = _utils_mod.SEPARATOR_LIGHT   # ─×80 - house step / subsection rule
_DEFAULT_FRAME_TOTAL = 100_000   # heartbeat fallback when the trajectory length is unreadable
EXIT_WARN = 3   # step completed but a complementary part (MM-GBSA) was deferred/failed; the pipeline runner renders WARN and continues (0=PASS, 1=hard error, 3=warn)

# Schrödinger jobs run under jobserverd, independent of this process - so on Ctrl-C / kill they would
# outlive the script. Every WaterMap/build/MD job name (and each MM-GBSA shard's Prime subjob) is registered
# here BEFORE submit and removed only on clean completion; a signal/atexit handler cancels whatever is still
# registered (see _install_job_cleanup).
_LAUNCHED_JOBS: "set[str]" = set()
# Local Schrödinger subprocesses that do NOT go through the job server - SID (analyze_simulation/event_analysis,
# run -LOCAL) and the MM-GBSA thermal_mmgbsa drivers. These reparent to init on a kill and keep running, so the
# same handler terminates each one's process group. Launched via _run_tracked, which registers and removes them.
_LAUNCHED_PROCS: "set[subprocess.Popen]" = set()
_CLEANUP_DONE = False

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
# writes. The per-\r progress bars (written straight to sys.stdout) are deliberately NOT mirrored -
# a log does not want carriage-return redraws.
_LOG_FH = None
_ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def _open_step_log(physics_dir: Path) -> None:
    """Open 6_Physics_Validation/00_Physics_Validation.log for this run (fresh each run)."""
    global _LOG_FH
    try:
        physics_dir.mkdir(parents=True, exist_ok=True)
        _LOG_FH = open(physics_dir / "00_Physics_Validation.log", "w", encoding="utf-8")
        atexit.register(lambda: _LOG_FH and not _LOG_FH.closed and _LOG_FH.close())
    except Exception:
        _LOG_FH = None


# One live \r progress line at a time across all threads. The MD heartbeat (main thread) and the CPU
# post-processing bars (worker thread - extraction, SID, MM-GBSA) share the terminal, so every \r write
# and every full-line _echo takes this lock: a full line first clears the open bar, and no two writers
# interleave bytes on the same line. Progress bars are never mirrored to the log file - only _echo lines
# are - so the on-disk log stays clean of carriage-return redraws.
_CONSOLE_LOCK = threading.RLock()
_BAR_OPEN = False


def _progress_line(msg: str) -> None:
    """Refresh the single shared \r progress line (terminal only, never the log)."""
    global _BAR_OPEN
    with _CONSOLE_LOCK:
        sys.stdout.write(f"\r    [PROGRESS] {msg}\033[K")
        sys.stdout.flush()
        _BAR_OPEN = True


def _close_bar() -> None:
    """Close the open \r progress line with a newline (phase change / completion)."""
    global _BAR_OPEN
    with _CONSOLE_LOCK:
        if _BAR_OPEN:
            sys.stdout.write("\n")
            sys.stdout.flush()
            _BAR_OPEN = False


def _worker_progress(msg: str) -> None:
    """Refresh the shared single \r progress line from a CPU worker (extraction / SID / MM-GBSA).

    Uses the same one-line \r bar as the MD heartbeat (_progress_line), so a long SID or MM-GBSA read
    updates IN PLACE instead of scrolling a new line every tick. When a CPU worker and the GPU MD
    heartbeat run concurrently they share that one line (serialised by _CONSOLE_LOCK): the line shows
    whichever ticked last. Terminal only - never the log, so the log keeps only the phase start/finish
    lines, not thousands of progress ticks."""
    _progress_line(msg)


def _echo(msg: str = "") -> None:
    """Print to terminal immediately (flush) - keeps live progress visible - and mirror to the
    step log file KEEPING ANSI colour, so `tail -f` of the merged log shows the same green/red.
    Clears any open \r progress line first so a full line never appends to a live bar."""
    global _BAR_OPEN
    with _CONSOLE_LOCK:
        if _BAR_OPEN:
            sys.stdout.write("\r\033[K")
            _BAR_OPEN = False
        print(msg, flush=True)
        if _LOG_FH is not None:
            try:
                _LOG_FH.write(str(msg) + "\n")
                _LOG_FH.flush()
            except Exception:
                pass


# Colour helpers used by the ESP/build/WaterMap/MD phase functions (single merged log).
SCHRO = SCHRODINGER
_C = _utils_mod.ConsoleColours


def _log(m: str = "") -> None:
    _echo(f"  {m}")


def _ok(m: str) -> None:
    _echo(f"  {_C.OKGREEN}{m}{_C.ENDC}")


def _fail(m: str) -> None:
    _echo(f"  {_C.FAIL}{m}{_C.ENDC}")


def _warn(m: str) -> None:
    _echo(f"  {_C.WARNING}[!] {m}{_C.ENDC}")


def _section(title: str) -> None:
    """House-style section header - bold title + light rule (matches ReportManager.section in 00_02)."""
    _echo(f"\n{_C.BOLD}{title}{_C.ENDC}")
    _echo(_RULE)


# ── Per-phase / per-job wall-clock timing ─────────────────────────────────────────────────────────
"""
Every heavy sub-job (one WaterMap, one System Builder, one MD, one SID, one MM-GBSA) records its
own wall-clock here via _timed(). The Summary then prints the individual job times, the per-phase
totals (all WaterMaps, all builds, …) and a grand total, and writes them to 00_Phase_Timings.csv.
"""
_TIMINGS: list = []
"""
Settle gate after an MD job lands (and after each rank's CPU post-processing finishes): the job
directory is re-scanned until no file's size or mtime changes for _SETTLE_STABLE_SEC, so a fresh run
never collides with an in-flight write. Watching the files (rather than sleeping a fixed duration)
covers a large multisim that takes a variable time to flush, archive its production tgz and unpack
the trajectory. _SETTLE_TIMEOUT_SEC caps the wait so a stuck archive never blocks forever.
"""
_SETTLE_STABLE_SEC = 45
_SETTLE_POLL_SEC = 15
_SETTLE_TIMEOUT_SEC = 3600


def _wait_settled(watch: "Path", stable_sec: int = _SETTLE_STABLE_SEC,
                  poll: int = _SETTLE_POLL_SEC, timeout: int = _SETTLE_TIMEOUT_SEC) -> bool:
    """Poll a job directory until its file set stops changing for stable_sec, then return True, so an
    MD job still flushing/archiving (multisim tgz, trajectory unpack) is never read or overwritten
    mid-write. Returns False on timeout (caller proceeds anyway)."""
    _t0 = time.time()
    _prev = None
    _stable = 0
    while time.time() - _t0 < timeout:
        _snap = []
        try:
            for _p in watch.rglob("*"):
                try:
                    _st = _p.stat()
                except OSError:
                    continue
                if stat.S_ISREG(_st.st_mode):
                    _snap.append((str(_p), _st.st_size, int(_st.st_mtime)))
        except OSError:
            _snap = None
        if _snap is not None:
            _snap.sort()
            if _snap == _prev:
                _stable += poll
                if _stable >= stable_sec:
                    return True
            else:
                _prev = _snap
                _stable = 0
        time.sleep(poll)
    return False


def _fmt_dur(sec: float) -> str:
    """Seconds → compact H/M/S (2h07m03s · 8m12s · 41s)."""
    sec = int(round(sec))
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m{s:02d}s"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def _eta_str(done: float, total: float, elapsed_sec: float) -> str:
    """A ' · ETA <dur>' suffix from a linear extrapolation of the current rate - empty until there is
    enough progress to extrapolate (done in (0, total)), so a just-started or finished step shows none."""
    if done <= 0 or done >= total or elapsed_sec <= 0:
        return ""
    return f" · ⏱ ETA {_fmt_dur(elapsed_sec * (total - done) / done)}"


@contextmanager
def _timed(phase: str, rank):
    """Time one sub-job, log its wall-clock, and register it for the phase totals + timings CSV."""
    _t = time.perf_counter()
    try:
        yield
    finally:
        dt = time.perf_counter() - _t
        _TIMINGS.append({"phase": phase, "rank": rank, "seconds": round(dt, 1)})
        if dt >= 1.0:   # skip sub-second echoes (all-cached reruns); still recorded for the phase totals
            _log(f"       ⏱ {phase} R_{rank} took {_fmt_dur(dt)}")


def _emit_timings(out_root: Path) -> None:
    """Print per-phase job times + totals and write 00_Phase_Timings.csv (atomic). No-op if nothing ran."""
    if not _TIMINGS:
        return
    order = ["watermap", "build", "md", "sid", "mmgbsa"]
    label = {"watermap": "WaterMap", "build": "System Builder", "md": "MD Simulation",
             "sid": "SID Analysis", "mmgbsa": "MM-GBSA"}
    by: dict = {}
    for t in _TIMINGS:
        by.setdefault(t["phase"], []).append(t)

    _section("Timing - per job (individual) + per-phase totals")
    grand = 0.0
    for ph in order:
        rows = by.get(ph)
        if not rows:
            continue
        tot = sum(r["seconds"] for r in rows)
        grand += tot
        per = " · ".join(f"R_{r['rank']} {_fmt_dur(r['seconds'])}" for r in rows)
        _echo(f"  {label[ph]:<15} {_fmt_dur(tot):>10}   ({len(rows)} job{'s' if len(rows) != 1 else ''}: {per})")
    _echo(_RULE)
    _echo(f"  {'TOTAL (jobs)':<15} {_fmt_dur(grand):>10}")

    # CSV: one row per job, then a per-phase TOTAL row (rank=ALL) and a grand-TOTAL row.
    csv_path = out_root / "00_Phase_Timings.csv"
    try:
        tmp = csv_path.with_suffix(".csv.tmp")
        with tmp.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["phase", "rank", "seconds", "hms"])
            for ph in order:
                for r in by.get(ph, []):
                    w.writerow([label[ph], f"R_{r['rank']}", f"{r['seconds']:.1f}", _fmt_dur(r["seconds"])])
            for ph in order:
                rows = by.get(ph)
                if rows:
                    tot = sum(r["seconds"] for r in rows)
                    w.writerow([label[ph], "ALL", f"{tot:.1f}", _fmt_dur(tot)])
            w.writerow(["TOTAL", "ALL", f"{grand:.1f}", _fmt_dur(grand)])
        os.replace(tmp, csv_path)                          # atomic - a kill never leaves a half-written CSV
        _echo(f"  {_C.OKGREEN}✔{_C.ENDC} timings CSV  → {csv_path.name}")
    except Exception as exc:
        _warn(f"[timing] could not write {csv_path.name}: {str(exc).splitlines()[0]}")


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
    traj API. This is the denominator for completion - it adapts automatically
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
    assumed - a 100 ns and a 1000 ns run both report themselves correctly. Returns 0.0
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
    non-empty EAF is treated as complete - never destructively re-run on doubt.
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
        self._entered = False
        self._keepalive_stop = threading.Event()
        self._keepalive_thread: threading.Thread | None = None

    def __enter__(self):
        if self._entered:          # idempotent: primed once at start, reused by the MD-phase `with`
            return self
        self._entered = True
        if not self.active:
            _echo("  [PIPELINE-MODE] oomd management delegated to pipeline runner.")
            return self
        _echo(f"{_C.BOLD}OPTIONAL - protect this run from the Linux out-of-memory (OOM) killer.{_C.ENDC}")
        _echo("  SID and Prime MM-GBSA hold large trajectories in memory for hours; systemd-oomd can "
              "kill them mid-run. Masking systemd-oomd needs root.")
        _echo(f"  {_C.FAIL}{_C.BOLD}Enter your sudo password to mask systemd-oomd, or press Enter / "
              f"Ctrl-D to skip and run unprotected:{_C.ENDC}")
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
            _echo("  [NORMAL MODE] sudo unavailable/skipped - systemd-oomd NOT masked (run unprotected from the OOM-killer).")
            self.active = False
            return self

        def _keepalive():
            while not self._keepalive_stop.wait(60):
                subprocess.run(["sudo", "-vn"], capture_output=True)

        self._keepalive_thread = threading.Thread(target=_keepalive, daemon=True)
        self._keepalive_thread.start()

        _echo("  Masking systemd-oomd - it will be restored automatically when Step 06 exits.")
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
#   Heartbeat      - one process, one log  (SID; serial MM-GBSA)
#   ShardHeartbeat - many shard logs at once (sharded MM-GBSA)
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

    # PHASE 2 banner styling from the shared ConsoleColours palette (passes through the pipeline
    # tee to a live terminal; harmless byte-noise in a pure-file redirect).
    _C_PHASE2 = _C.BOLD + _C.CYAN   # bold cyan
    _C_DIM    = _C.DIM              # dim
    _C_RST    = _C.ENDC             # reset

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
        self._prime_announced = False   # PHASE 2 banner printed once

    def _progress(self, msg: str) -> None:
        """Report progress as a full logged line (runs on the CPU worker, concurrent with the MD bar)."""
        _worker_progress(msg)

    def _end_line(self) -> None:
        """Close the shared \r progress line at phase change / completion (newline)."""
        _close_bar()

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
        the job server writes. Returns (n, total) - total falls back to the frame
        count - or ``None`` if Prime has not logged a per-structure marker yet."""
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
            elapsed_sec = time.time() - self._start
            elapsed_min = int(elapsed_sec // 60)
            try:
                text = self.log_file.read_text(errors="ignore")
            except Exception:
                text = ""
            # Prime MM-GBSA scoring phase: trajectory read is complete; Prime now
            # minimises every structure. Announce the transition once (colour banner
            # on its own line), then keep the single \r line refreshing - reporting
            # WHICH structure is being minimised (structure N/total, like PHASE 1's
            # frame count), falling back to subjob count, then elapsed time.
            if self._PRIME_HANDOFF.search(text):
                if not self._prime_announced:
                    self._prime_announced = True
                    self._end_line()   # close the PHASE 1 \r line with a newline
                    _echo(f"{self._C_PHASE2}    ══════ [PHASE 2/2] {self.label}: every trajectory "
                          f"frame has been read. Prime is now minimising and scoring "
                          f"{self.total:,} structures ══════{self._C_RST}")
                    _echo(f"{self._C_DIM}    This is the long phase - it can run for hours. The "
                          f"line below refreshes in place with the structure Prime is on."
                          f"{self._C_RST}")
                _struct = self._latest_structure(text)
                if _struct is not None:
                    _n, _tot = _struct
                    _pct = min(100, _n * 100 // max(_tot, 1))
                    self._progress(f"{self.label}: PHASE 2/2 Prime minimised structure {_n:,} of "
                                   f"{_tot:,} ({_pct}%, {elapsed_min}m elapsed"
                                   f"{_eta_str(_n, _tot, elapsed_sec)})")
                else:
                    prog = self._prime_progress()
                    if prog is not None:
                        _done, _tot = prog
                        _pct = min(100, _done * 100 // max(_tot, 1))
                        self._progress(f"{self.label}: PHASE 2/2 Prime finished {_done} of {_tot} "
                                       f"subjobs ({_pct}%, {elapsed_min}m elapsed)")
                    else:
                        # Prime spinning up - no per-structure marker logged yet.
                        self._progress(f"{self.label}: PHASE 2/2 Prime is starting up on "
                                       f"{self.total:,} structures ({elapsed_min}m elapsed)")
                continue
            # Frame-read progress (SID's whole run, and MM-GBSA PHASE 1). Kept
            # phase-neutral so it reads correctly for SID, which has no PHASE 2.
            frame = self._latest_frame(text)
            pct = min(100, frame * 100 // self.total)
            self._progress(f"{self.label}: {self.read_phase}read frame {frame:,} of "
                           f"{self.total:,} ({pct}%, {elapsed_min}m elapsed"
                           f"{_eta_str(frame, self.total, elapsed_sec)})")

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
        self._scored0 = 0        # shards already scored when this run began; the ETA rate ignores them
                                 # so a resume does not credit the previous run's shards to this elapsed

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
            elapsed_sec = time.time() - self._start
            mins = int(elapsed_sec // 60)
            read, priming = self._counts()
            scored = self.done_fn()
            pct = min(100, read * 100 // self.total_frames)
            priming = max(0, priming - scored)
            # ETA off shards completed since this run began, so a resume's pre-scored shards do not
            # inflate the rate; the remaining count is still measured to the full shard total.
            eta = _eta_str(scored - self._scored0, self.n_shards - self._scored0, elapsed_sec)
            _worker_progress(
                f"{self.label}: read {read:,} of {self.total_frames:,} frames "
                f"({pct}%) · {scored}/{self.n_shards} shards scored · {priming} minimising in "
                f"Prime · {mins}m elapsed{eta}")

    def __enter__(self):
        self._start = time.time()
        # Baseline = shards already scored on disk from a previous run (resumed). done_fn() is 0 at this
        # instant - the worker pool has not registered the resume-skips yet - so reading it here would
        # credit those instant resumes to this run's clock and make the ETA far too short. Count the CSVs.
        try:
            self._scored0 = sum(1 for _ in self.shard_dir.glob("*_shard*-prime-out.csv"))
        except Exception:
            self._scored0 = 0
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1)
        return False


class MDHeartbeat:
    """One in-place progress line for a running Desmond MD, driven by the live ``.ene``.

    Desmond appends one row per recorded step to ``<stage>.ene`` whose first column is the chemical
    (simulated) time in ps; the actively-growing .ene (newest mtime) is the stage running right now.
    Its last time, against the known production length, gives current ns / total ns and a percentage;
    successive samples give the throughput (ns/day) and an ETA. Before production starts (the short
    relaxation stages), it reports which relaxation stage multisim is on instead of a percentage.
    Mirrors Heartbeat's single-``\\r``-line contract - nothing else runs concurrently to interleave.
    """

    _HDR_PROD = re.compile(r"stage\s+(\d+)\s*-\s*.*Production", re.IGNORECASE)
    _HDR_ANY  = re.compile(r"stage\s+(\d+)\s*-", re.IGNORECASE)
    _STAGE_RUN = re.compile(r"Stage\s+(\d+)\s*-\s*(?:simulate|task)", re.IGNORECASE)
    _STAGE_DONE = re.compile(r"Stage\s+(\d+)\s+completed", re.IGNORECASE)

    def __init__(self, wd: Path, multisim_log: Path, label: str, total_ns: float, interval: int = 30):
        self.wd = Path(wd)
        self.multisim_log = Path(multisim_log)
        self.label = label
        self.total_ns = max(float(total_ns), 1e-6)
        self.total_ps = self.total_ns * 1000.0
        self.interval = max(5, int(interval))
        self._stop = threading.Event()
        self._thread: "threading.Thread | None" = None
        self._start = 0.0
        self._line_open = False
        self._last: "tuple[float, float] | None" = None   # (wall_s, t_ps) of the previous sample, for the rate
        # multisim runs the job under jsc, so the LIVE multisim log and the production .ene do not sit in
        # wd - they live in the job-server scratch, /tmp/<user>/jobs/<jobid>/ (and the production stage is
        # its own subjob, in a sibling <jobid>/ dir). Search those too, keyed by jobname, newest mtime wins.
        self._jobname = self.wd.name
        _user = os.environ.get("USER") or Path.home().name
        _roots = {os.environ.get("SCHRODINGER_TMPDIR") or "", tempfile.gettempdir(), "/tmp"}
        self._scratch_bases = [Path(r) / _user / "jobs" for r in _roots if r]

    def _log_text(self) -> str:
        """Text of the live multisim log: wd if present, else the newest {jobname}_multisim.log in the
        job-server scratch (where a jsc-controlled run keeps it until the job finishes)."""
        cands = [self.multisim_log] if self.multisim_log.exists() else []
        for base in self._scratch_bases:
            try:
                cands.extend(p for p in base.glob(f"*/{self._jobname}_multisim.log") if p.is_file())
            except Exception:
                pass
        if not cands:
            return ""
        try:
            return max(cands, key=lambda p: p.stat().st_mtime).read_text(errors="ignore")
        except Exception:
            return ""

    def _newest_ene_time(self) -> "float | None":
        """Last chemical time (ps) in the most-recently-written .ene for this job - searched in wd AND the
        job-server scratch (the production stage writes its .ene in a subjob scratch dir, never in wd)."""
        enes = [p for p in self.wd.rglob(f"{self._jobname}*.ene") if p.is_file()]
        for base in self._scratch_bases:
            try:
                enes.extend(p for p in base.glob(f"*/{self._jobname}*.ene") if p.is_file())          # subjob-root .ene
                enes.extend(p for p in base.glob(f"*/{self._jobname}_*/{self._jobname}*.ene") if p.is_file())  # stage subdir
            except Exception:
                pass
        if not enes:
            return None
        newest = max(enes, key=lambda p: p.stat().st_mtime)
        try:
            tail = newest.read_bytes()[-4096:].decode(errors="ignore").splitlines()
        except Exception:
            return None
        for ln in reversed(tail):
            ln = ln.strip()
            if ln and not ln.startswith("#"):
                try:
                    return float(ln.split()[0])
                except (ValueError, IndexError):
                    continue
        return None

    def _in_production(self) -> bool:
        """True once multisim has entered the Production stage (and not yet completed it)."""
        text = self._log_text()
        if not text:
            return False
        prod = self._HDR_PROD.search(text)
        stages = [int(m) for m in self._HDR_ANY.findall(text)]
        prod_stage = int(prod.group(1)) if prod else (max(stages) if stages else None)
        if prod_stage is None:
            return False
        started = prod_stage in {int(m) for m in self._STAGE_RUN.findall(text)}
        done = prod_stage in {int(m) for m in self._STAGE_DONE.findall(text)}
        return started and not done

    def _relax_stage(self) -> "tuple[int, int] | None":
        """(current relaxation stage, production stage number) for the pre-production phase."""
        text = self._log_text()
        if not text:
            return None
        stages = [int(m) for m in self._HDR_ANY.findall(text)]
        prod = self._HDR_PROD.search(text)
        prod_stage = int(prod.group(1)) if prod else (max(stages) if stages else 7)
        run = [int(m) for m in self._STAGE_RUN.findall(text)]
        return (max(run) if run else 1), prod_stage

    def _progress(self, msg: str) -> None:
        _progress_line(msg)
        self._line_open = True

    def _run(self) -> None:
        while not self._stop.wait(self.interval):
            mins = int((time.time() - self._start) // 60)
            if not self._in_production():
                rs = self._relax_stage()
                if rs is not None:
                    st, prod = rs
                    self._progress(f"{self.label}: equilibrating - relaxation stage {st}/{prod} "
                                   f"({mins}m elapsed)")
                else:
                    self._progress(f"{self.label}: starting up ({mins}m elapsed)")
                continue
            t_ps = self._newest_ene_time()
            if t_ps is None:
                self._progress(f"{self.label}: production starting ({mins}m elapsed)")
                continue
            now = time.time()
            ns_done = min(t_ps, self.total_ps) / 1000.0
            pct = min(100.0, 100.0 * t_ps / self.total_ps)
            rate = ""      # ns/day between the last two samples
            eta = ""
            if self._last is not None:
                dw, dpt = now - self._last[0], t_ps - self._last[1]
                if dw > 0 and dpt > 0:
                    ns_per_day = (dpt / 1000.0) / (dw / 86400.0)
                    rate = f" · {ns_per_day:,.0f} ns/day"
                    remain_ns = max(0.0, self.total_ns - ns_done)
                    eta_s = remain_ns / max(ns_per_day, 1e-9) * 86400.0
                    eta = f" · ⏱ ETA {_fmt_dur(eta_s)}"
            self._last = (now, t_ps)
            self._progress(f"{self.label}: {ns_done:.2f}/{self.total_ns:.2f} ns ({pct:.1f}%)"
                           f"{rate}{eta} · {mins}m elapsed")

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
            _close_bar()
            self._line_open = False
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
    _echo("Usage: python 06_Physics_Validation_DeFluorX.py <Boltz-2_Run_Directory>/ [--test]")
    _echo("  quick test (WaterMap 2 ns · MD 5 ns): python 06_Physics_Validation_DeFluorX.py Boltz-2_Run_20260309T085406Z/ --test")
    sys.exit(1)


# =============================================================================
# SECTION 6: SCAN PHASE
# =============================================================================
def scan_jobs(job_dirs: list[Path]) -> list[Path]:
    """Classify every desmond_md_job_R_* directory and return the subset
    that still needs SID analysis (READY or INCOMPLETE).
    """
    _multi = len(job_dirs) > 1            # the full scan header/summary is noise for a single-rank scan
    if _multi:
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
                _echo(f"  - Rank {rank}: MD simulation is still running - SID must wait for it.")
            else:
                pending += 1
                _echo(f"  - Rank {rank}: no MD output (-out.cms) and nothing running - SID cannot start.")
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
                _echo(f"  - Rank {rank}: SID already complete - all {tf:,} trajectory frames analysed.")
            else:
                to_run.append(d)
                _echo(f"  - Rank {rank}: SID stopped early - only {of:,} of {tf:,} frames analysed, "
                      f"so it will be re-run from scratch.")
        else:
            to_run.append(d)
            _echo(f"  - Rank {rank}: MD finished but SID has never been run - queued.")

    if _multi:
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
    """Step 1: event_analysis.py - generates the SID-in.eaf descriptor."""
    if in_eaf.is_file():
        _echo(f"  [SKIP] SID-in.eaf already exists: {in_eaf}")
        return
    _echo("  Running event_analysis.py to generate SID-in.eaf...")
    log = job_dir / f"{job_name}_event_analysis.log"
    with open(log, "w") as fh:
        _run_tracked(
            [SCHROD_RUN, "event_analysis.py", "analyze", str(cms_file),
             "-out", f"{job_name}_SID"],
            cwd=str(job_dir), stdout=fh, stderr=subprocess.STDOUT, check=True,
        )


def run_analyze_simulation(job_dir: Path, job_name: str, cms_file: Path,
                           trj_dir: Path, out_eaf: Path, label: str) -> int:
    """Step 2: analyze_simulation.py - generates the SID-out.eaf result vector.

    Uses -LOCAL to avoid remote-server overhead and cap memory. Returns the
    trajectory frame total used as the completion denominator.
    """
    total_frames = traj_frame_count(trj_dir) or _DEFAULT_FRAME_TOTAL
    _echo(f"  Running analyze_simulation.py (SID over {total_frames:,} frames; live progress below)…")
    log = job_dir / f"{job_name}_analyze_simulation.log"
    with Heartbeat(log, label, total_frames):
        with open(log, "w") as fh:
            _run_tracked(
                [SCHROD_RUN, "analyze_simulation.py", "-NOJOBID", "-LOCAL",
                 str(cms_file), str(trj_dir),
                 f"{job_name}_SID-out.eaf", f"{job_name}_SID-in.eaf"],
                cwd=str(job_dir), stdout=fh, stderr=subprocess.STDOUT, check=True,
            )
    return total_frames


def process_jobs(to_run: list[Path], job_index: int | None = None,
                 job_total: int | None = None) -> None:
    """Iterate the work list, generating SID-in then SID-out for each job.

    Each rank is isolated: a Schrödinger failure (event_analysis/analyze_simulation)
    or an incomplete EAF is logged and skipped so the remaining ranks still run.
    Step 07 consumes whichever *_SID-out.eaf files completed.

    job_index/job_total, when supplied by the sequential driver, label the banner
    with the rank's position across the whole run (e.g. Job 2/3). Without them the
    counter falls back to this call's local work list.
    """
    total = len(to_run)
    failures: list[str] = []
    for index, d in enumerate(to_run, start=1):
        job_name = d.name
        rank = _rank_of(job_name)
        disp_index = job_index if job_index is not None else index
        disp_total = job_total if job_total is not None else total
        out_eaf = d / f"{job_name}_SID-out.eaf"
        in_eaf = d / f"{job_name}_SID-in.eaf"
        cms_file = d / f"{job_name}-out.cms"
        trj_dir = d / f"{job_name}_trj"

        _echo("")
        _echo(_SEP)
        # Re-run (incomplete) vs. fresh processing - preserve the partial EAF.
        if out_eaf.is_file():
            of = out_eaf_frames(out_eaf)
            tf = traj_frame_count(trj_dir)
            _echo(f"[Job {disp_index}/{disp_total}] Re-processing (Incomplete): {job_name} (Rank {rank})")
            _echo(_SEP)
            _echo(f"  [RE-RUN] Output EAF is incomplete ({of}/{tf} frames). Re-running analysis...")
            # A misjudged "incomplete" must never destroy a good result.
            stamp = time.strftime("%Y%m%d%H%M%S")
            bak = out_eaf.with_name(f"{out_eaf.name}.incomplete.{stamp}.bak")
            out_eaf.replace(bak)
            _echo(f"  [BACKUP] Previous EAF moved to: {bak}")
        else:
            _echo(f"[Job {disp_index}/{disp_total}] Processing: {job_name} (Rank {rank})")
            _echo(_SEP)

        # Verify the trajectory folder is present (dir or packed .xtc).
        if not trj_dir.is_dir() and not (trj_dir.parent / f"{trj_dir.name}.xtc").is_file():
            _echo(f"  [WARNING] Trajectory folder not found: {trj_dir}. Skipping Rank {rank}.")
            continue

        label = f"Job {disp_index}/{disp_total} Rank {rank}"
        try:
            run_event_analysis(d, job_name, cms_file, in_eaf)
            total_frames = run_analyze_simulation(d, job_name, cms_file, trj_dir, out_eaf, label)
        except (subprocess.CalledProcessError, OSError) as exc:
            _echo(f"  [ERROR] SID analysis failed for Rank {rank} ({job_name}): {exc}. Skipping rank.")
            failures.append(f"Rank {rank} ({job_name})")
            continue

        if out_eaf.is_file() and is_eaf_complete(out_eaf_frames(out_eaf), total_frames):
            of = out_eaf_frames(out_eaf)
            _echo(f"  {_C.OKGREEN}[SID DONE] Rank {rank}: {out_eaf.name} complete ({of:,}/{total_frames:,} frames).{_C.ENDC}")
            _echo(f"  [LOGS] {job_name}_event_analysis.log · {job_name}_analyze_simulation.log")
            _echo(f"  [VIEW] Open {out_eaf.name} in Maestro's Simulation Interaction Diagram "
                  f"(Tasks → Analyze → Simulation Interaction Diagram) to inspect the protein Cα-RMSD, "
                  f"per-residue RMSF, ligand RMSD and the protein-ligand interaction timeline.")
            if getattr(CFG, "MD_RESTRAIN_LIGAND", False):
                _echo(f"  {_C.WARNING}[NOTE] Production ran under a positional restraint (ligand heavy "
                      f"atoms k={CFG.MD_RESTRAIN_LIG_FORCE_K}, backbone k={CFG.MD_RESTRAIN_BB_FORCE_K} "
                      f"kcal/mol/Å²): the ligand-RMSD / pocket retention shown in SID is RESTRAINT-ENFORCED, "
                      f"not spontaneous. The unrestrained reactivity verdict is the Step-07 QM/MM barrier.{_C.ENDC}")
        else:
            _echo(f"  [ERROR] Incomplete output for Rank {rank}: {out_eaf}. Skipping rank "
                  f"(re-run to retry - the partial EAF is preserved).")
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
#   8.2  execution - sharded (default) and serial
#   8.3  statistics - ΔG estimators, failed-minimisation flagging
#   8.4  figures - per job and combined
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


def _mmgbsa_expected_rows(job_dir: Path, job_name: str) -> int:
    """Structures Prime should have scored = trajectory frames folded by the stride."""
    total = traj_frame_count(job_dir / f"{job_name}_trj") or _DEFAULT_FRAME_TOTAL
    step = max(1, int(getattr(CFG, "MMGBSA_STEP_SIZE", 0) or 1))
    return total if step <= 1 else -(-total // step)


def _mmgbsa_complete(csv: Path, job_dir: Path, job_name: str) -> bool:
    """A merged MM-GBSA CSV is 100% complete when its row count covers every
    strided structure of the trajectory. A short CSV (older/partial) is not."""
    try:
        return len(pd.read_csv(csv)) >= _mmgbsa_expected_rows(job_dir, job_name)
    except Exception:
        return False


def _cleanup_mmgbsa_shards(job_dir: Path, job_name: str, rank: str) -> None:
    """Remove the per-job _MMGBSA_Shards scratch once the merged CSV is settled.

    The merged <job>_mmgbsa-prime-out.csv (beside the job folder) holds every ΔG and is
    the only MM-GBSA product Steps 06/07 read. The shard directory (per-shard maegz
    complexes, logs, symlinks) is intermediate - several GB per rank - and nothing
    downstream consumes it. Called only after a verified-complete merge, at the very end
    of a rank's MM-GBSA. Symlinks inside are unlinked, not followed, so the real cms/_trj
    are untouched."""
    shard_dir = job_dir / getattr(CFG, "MMGBSA_SHARD_SUBDIR", "_MMGBSA_Shards")
    if not shard_dir.is_dir():
        return
    freed = 0
    for p in shard_dir.rglob("*"):
        try:
            if p.is_file() and not p.is_symlink():
                freed += p.stat().st_size
        except OSError:
            pass
    shutil.rmtree(shard_dir, ignore_errors=True)
    _echo(f"    ✔ Cleanup      : removed {shard_dir.name}/ "
          f"({freed / 2**30:.1f} GB of shard scratch) - merged CSV kept.")


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
        return ("DISK FULL - the Schrödinger job server ran out of space staging "
                "per-subjob scratch (each Prime subjob copies the multi-GB complexes "
                "file). Free space on the job server's scratch disk (its tmpdir, /tmp by "
                "default) or lower the concurrent subjob count, then retry.")
    if re.search(r"licen[sc]e", text, re.I) and re.search(r"error|fail|not available|checkout", text, re.I):
        return "LICENSE - a Prime/PLOP (PSP_PLOP) license was unavailable; check FlexLM."
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


def _free_swap_gb() -> float:
    """Free swap (SwapFree). Only added to the MM-GBSA budget when MMGBSA_RAM_SWAP_FRAC > 0 -
    Prime whose working set lands on a swapfile runs at disk speed and can trip the OOM killer."""
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("SwapFree:"):
                return int(line.split()[1]) / 2**20
    except Exception:
        pass
    return 0.0


def _shard_plan(total: int, ncpu: int, step: int = 1) -> "tuple[list[tuple[int, int]], int, int]":
    """Split [0, total) into contiguous frame shards and decide how many run at once.

    Returns (ranges, concurrency, prime_njobs_per_shard). The pair is chosen to keep the most
    cores busy on this machine's CPU and free RAM, then a final pass guarantees reader + Prime
    fit the RAM budget. With MMGBSA_SHARD_PRIME_NJOBS = 0 the njobs is auto-picked too; a positive
    value pins it and only concurrency is derived.
    """
    size = max(1, int(getattr(CFG, "MMGBSA_SHARD_FRAMES", 2000)))
    ranges = [(a, min(a + size, total)) for a in range(0, total, size)]

    """
    Per-shard cost feeds both the RAM fit and the core count. A shard is one frame reader (a fixed
    RSS base plus growth with the strided frames it reads) plus its Prime subjobs - the reader is
    one core, the Prime workers the rest. The budget is physical RAM the kernel hands out without
    swapping (MemAvailable), optionally plus a slice of free swap: Prime whose working set spills to
    a swapfile runs at disk speed and can trip the OOM killer, so swap is opt-in (MMGBSA_RAM_SWAP_FRAC,
    default 0). thermal_mmgbsa does not saturate a large -NJOBS - it runs about
    MMGBSA_PRIME_EFFECTIVE_CORES live workers per shard - so a shard's core benefit is scored against
    that effective figure, not the request, and njobs above it would only reserve RAM for idle cores.
    """
    _read_n = max(1, size // max(1, step))
    _reader = (float(getattr(CFG, "MMGBSA_READER_RAM_BASE_GB", 3.6))
               + float(getattr(CFG, "MMGBSA_READER_RAM_PER_1K_FRAMES_GB", 0.36)) * _read_n / 1000.0)
    _prime = max(0.2, float(getattr(CFG, "MMGBSA_PRIME_RAM_GB", 1.8)))
    _eff = max(1, int(getattr(CFG, "MMGBSA_PRIME_EFFECTIVE_CORES", 3)))
    _ram = (_avail_ram_gb()
            + max(0.0, float(getattr(CFG, "MMGBSA_RAM_SWAP_FRAC", 0.0))) * _free_swap_gb()
            ) * float(getattr(CFG, "MMGBSA_RAM_HEADROOM_FRAC", 0.85))

    def _fit_conc(nj: int) -> int:
        """Concurrent shards this njobs allows - the tighter of the CPU and RAM limits."""
        c_cpu = max(1, ncpu // (1 + nj))                        # one reader core + nj Prime cores
        c_ram = int(_ram // (_reader + nj * _prime)) if _ram > 0 else c_cpu
        return max(1, min(c_cpu, max(1, c_ram), len(ranges)))

    njobs_cfg = int(getattr(CFG, "MMGBSA_SHARD_PRIME_NJOBS", 0) or 0)
    conc_cfg = int(getattr(CFG, "MMGBSA_SHARD_CONCURRENCY", 0) or 0)
    if njobs_cfg > 0:
        njobs = njobs_cfg
        conc = conc_cfg or _fit_conc(njobs)
    else:
        """
        Auto. Search njobs and take the (njobs, concurrency) pair that keeps the most cores busy -
        more Prime subjobs per shard means fewer shards fit, so the two trade off. Busy cores are
        scored against the effective Prime figure, and ties break toward more concurrent shards,
        since overlapping readers are what hide the serial per-shard read phase.
        """
        best = None
        for nj in range(1, min(ncpu, 8) + 1):
            c = _fit_conc(nj)
            key = (c * (1 + min(nj, _eff)), c)
            if best is None or key > best[0]:
                best = (key, c, nj)
        _, conc, njobs = best
    conc = max(1, min(conc, len(ranges)))

    # Final RAM safety: trim concurrency, then njobs, until reader + Prime fit the budget.
    if _ram > 0:
        while conc > 1 and conc * (_reader + njobs * _prime) > _ram:
            conc -= 1
        while njobs > 1 and conc * (_reader + njobs * _prime) > _ram:
            njobs -= 1
    return ranges, conc, njobs


# Every mark colour in this step's figures comes from CFG (candidate colours from
# MMGBSA_RANK_PALETTE, everything else from MMGBSA_INK) - none is written here.
_INK = CFG.MMGBSA_INK

MMGBSA_FRAME_COL = "Frame"   # trajectory frame index each scored structure came from


def _stamp_frames(df: "pd.DataFrame", start: int, end: int, step: int,
                  label: str = "") -> "pd.DataFrame":
    """Record which trajectory frame each Prime row actually came from.

    thermal_mmgbsa scores frames range(start, end, step) in order, and Prime preserves that
    order, so row i is frame start + i·step. Writing that out explicitly is the whole point:
    downstream code (07's NAC-conditioned MM-GBSA) must never infer the frame from the row
    POSITION, because with a stride row 5,000 is frame 50,000 - a silent misattribution.
    The column is only stamped when the row count matches the frame range exactly; a
    mismatch means the assumption is broken and a wrong index is worse than none.
    """
    frames = list(range(int(start), int(end), max(1, int(step))))
    if len(df) == len(frames):
        df.insert(0, MMGBSA_FRAME_COL, frames)
    else:
        _echo(f"  [WARN] {label}: Prime returned {len(df)} rows for {len(frames)} requested "
              f"frames - frame indices NOT stamped. 07's NAC-conditioned MM-GBSA will fall "
              f"back to positional alignment, which is only valid at step_size ≤ 1.")
    return df


def _retrofit_frame_stamps(csv: Path, job_dir: Path, job_name: str, rank: str) -> None:
    """Add the `Frame` column to an MM-GBSA CSV that was written before stamping existed.

    Reconstructs the frame list exactly as the sharded run generated it - contiguous
    CFG.MMGBSA_SHARD_FRAMES blocks, every CFG.MMGBSA_STEP_SIZE-th frame within each, concatenated
    in shard order - and writes it in place, keeping the unstamped file as a .bak. Idempotent: a
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
              f"plan implies {len(frames):,} frames - NOT stamping (a wrong frame index is worse "
              f"than a missing one). Step 07 will skip the NAC-conditioned ΔG for this job.")
        return
    shutil.copy2(csv, csv.with_suffix(".csv.unstamped.bak"))
    df.insert(0, MMGBSA_FRAME_COL, frames)
    _tmp = csv.with_suffix(".csv.tmp")
    df.to_csv(_tmp, index=False)
    _tmp.replace(csv)
    _echo(f"    ✔ Frame-stamped: {csv.name} ({len(df):,} rows) - it predates the `Frame` column; "
          f"Step 07 needs it to align NAC frames with their binding energies.")


def _csv_nonempty(p: Path) -> bool:
    """True if p exists with content. A file that vanishes between the exists-check and the stat is
    treated as absent rather than crashing the pool with FileNotFoundError (TOCTOU-safe)."""
    try:
        return p.is_file() and p.stat().st_size > 0
    except OSError:
        return False


def _diagnose_shard_failure(slog: Path) -> "str | None":
    """Turn a shard's bare rc=1 into the actual cause, read from its log."""
    try:
        text = slog.read_text(errors="ignore")
    except Exception:
        return None
    if re.search(r"requires a locally running job server|Unable to submit job", text, re.I):
        return ("the Schrödinger local job server is DOWN, so Prime could not be submitted. "
                "Start it with `$SCHRODINGER/jsc local-server-start` and re-run - the frames "
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
    and cover [0, total), so every frame is scored - identical coverage to one serial run.
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
    cores, memory) are what a reader actually checks, and they are all derived - the trajectory
    span comes from the frames' own timestamps, so nothing here assumes a run length.
    """
    _span = traj_span_ns(job_dir / f"{job_name}_trj")
    _echo(f"    Trajectory     : {total:,} frames" + (f" · {_span:,.0f} ns" if _span > 0 else ""))
    if step > 1:
        _ps = (_span * 1000.0 * step / max(total, 1)) if _span > 0 else 0.0
        _echo(f"    Sampling       : every {step}th frame → {_scored:,} structures to Prime"
              + (f", one per {_ps:,.0f} ps" if _ps else ""))
        _echo("                     (below the ~0.1–1 ns decorrelation time of a bound pose, so the "
              "independent-sample count - and ⟨ΔG_bind⟩ - is unchanged)")
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
    Fail fast. A shard that dies for an environmental reason - job server down, no license,
    disk full - will kill every other shard the same way, each only AFTER re-reading its
    frames. Abort the rank once the first full wave has failed without a single success,
    so the cause is reported in a minute instead of an hour.
    """
    _abort = threading.Event()

    def _shard_csv(i: int) -> Path:
        return shard_dir / f"{job_name}_mmgbsa_shard{i:03d}-prime-out.csv"

    def _shard_complete(csv: Path, a: int, b: int) -> bool:
        """Accept a resumed shard only if it parses AND has EXACTLY the row count _stamp_frames expects
        for its frame range. thermal_mmgbsa writes non-atomically, so a kill mid-write leaves a truncated
        shard that is size>0 but short; concatenating it stamps NaN frames (invisible in plots, dropped
        wrongly by dg.drop on a NaN label). A short shard is re-run instead."""
        try:
            return len(pd.read_csv(csv)) == len(range(int(a), int(b), max(1, int(step))))
        except Exception:
            return False

    def _run_shard(i: int, a: int, b: int) -> None:
        csv = _shard_csv(i)
        if _csv_nonempty(csv) and _shard_complete(csv, a, b):   # resume: already scored AND complete
            with _lock:
                done[i] = csv
            return
        if _abort.is_set():
            return
        sname = f"{job_name}_mmgbsa_shard{i:03d}"
        cmd = [SCHROD_RUN, "thermal_mmgbsa.py", cms_file.name,
               "-j", sname, "-NJOBS", str(njobs),
               "-start_frame", str(a), "-end_frame", str(b)]
        if lig_asl:
            cmd += ["-lig_asl", lig_asl]
        if step > 0:
            cmd += ["-step_size", str(step)]
        slog = shard_dir / f"{sname}.log"
        # thermal_mmgbsa drives a Prime subjob named "{sname}-prime" on the job server; register it so a
        # kill cancels it (like MD), and run the driver itself through _run_tracked so its process group
        # is terminated too. Both are removed in the finally, so a completed shard leaves nothing registered.
        _prime_job = f"{sname}-prime"
        _LAUNCHED_JOBS.add(_prime_job)
        try:
            with open(slog, "w") as fh:
                rc = _run_tracked(cmd, cwd=str(shard_dir), stdout=fh, stderr=subprocess.STDOUT)
        except Exception as e:
            _echo(f"\n  [Rank {rank}] shard {i:03d} (frames {a:,}–{b:,}) failed to launch: {e}")
            with _lock:
                failed.append(i)
            return
        finally:
            _LAUNCHED_JOBS.discard(_prime_job)
        if rc != 0 or not _csv_nonempty(csv):
            _echo(f"\n  [Rank {rank}] shard {i:03d} (frames {a:,}–{b:,}) failed rc={rc} - see {slog.name}.")
            _why = _diagnose_shard_failure(slog)
            if _why:
                _echo(f"  [Rank {rank}] ↳ cause: {_why}")
            with _lock:
                failed.append(i)
                if not done and len(failed) >= conc:
                    _abort.set()
                    _echo(f"  [Rank {rank}] Aborting: the first {len(failed)} shards all failed and "
                          f"none succeeded - this is an environment problem, not a bad frame range. "
                          f"Remaining shards skipped.")
            return
        with _lock:
            done[i] = csv

    # Denominator is the number of frames Prime will actually see (stride applied), not the
    # raw trajectory length - otherwise the read percentage caps at 100/step and never reaches 100.
    def _done_count():
        with _lock:
            return len(done)
    hb = ShardHeartbeat(shard_dir, f"MM-GBSA Rank {rank}", len(ranges), _scored,
                        _done_count, interval=interval)
    with hb, cf.ThreadPoolExecutor(max_workers=conc) as pool:
        list(pool.map(lambda r: _run_shard(r[0], r[1][0], r[1][1]), list(enumerate(ranges))))

    # Each shard's thermal_mmgbsa -NJOBS should block, but if a build submits its Prime batch to
    # jobserverd and returns early, a shard's Prime subjobs can outlive the pool. Drain any job
    # whose name carries this rank's prefix before returning, so the next rank's MD never shares the
    # scratch disk with a still-running Prime batch (the serial path drains at run_mmgbsa's tail).
    _await_mmgbsa_jobserver(f"{job_name}_mmgbsa")

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
    merged = job_dir / f"{job_name}{CFG.SUFFIX_MMGBSA_CSV}"
    _all = pd.concat(frames, ignore_index=True)
    _tmp = merged.with_suffix(".csv.tmp")
    _all.to_csv(_tmp, index=False)
    _tmp.replace(merged)                  # atomic: a reader never sees a half-written CSV
    _stamped = "frame-stamped" if MMGBSA_FRAME_COL in _all.columns else "NOT frame-stamped"
    _echo(f"    ✔ Scored       : all {len(ranges)} shards → {merged.name} "
          f"({len(_all):,} rows, {_stamped})")
    _cleanup_mmgbsa_shards(job_dir, job_name, rank)   # merge settled → drop the shard scratch
    return merged


def _await_mmgbsa_jobserver(job_prefix: str, poll: int = 30, max_wait: int = 172800) -> None:
    """Block until no active job-server job whose name contains `job_prefix` remains (master + subjobs).

    Belt-and-braces after the thermal_mmgbsa driver returns: it should already have waited, but if a
    Schrödinger build submits the Prime batch to jobserverd and returns early, this keeps the caller
    blocked so the strictly-sequential MD → SID → MM-GBSA order holds and two Prime batches never share
    the scratch disk. Silent no-op if jsc is unavailable; bounded by max_wait so it can never hang forever.
    """
    _t0 = time.time()
    while time.time() - _t0 < max_wait:
        try:
            out = subprocess.run([f"{SCHRO}/jsc", "list", "-j"], capture_output=True,
                                 text=True, timeout=30).stdout
        except Exception:
            return
        if not any(job_prefix in ln for ln in out.splitlines()):
            return
        time.sleep(poll)


def _resolve_mmgbsa_lig_asl(cms_file: Path) -> str:
    """Resolve the -lig_asl handed to thermal_mmgbsa, size-independently.

    Preference order:
      1. CFG.MMGBSA_LIGAND_ASL (default "res.ptype LIG") when it selects ≥1 atom in the
         built complex - the reliable path: a 5-atom fluoroacetate is pinned exactly like a
         large PFAS, with no dependence on molecule size.
      2. Fallback when that name is absent (a rebuild under a different resname): the smallest
         non-protein / non-solvent molecule carrying at least CFG.LIGAND_MIN_ATOMS atoms,
         pinned by an explicit molecule ASL. This deliberately avoids thermal_mmgbsa's own
         AslLigandSearcher, whose 5-atom floor would drop a bare fluoroacetate.
      3. Empty string → let thermal_mmgbsa auto-detect (last resort)."""
    _cfg_asl = str(getattr(CFG, "MMGBSA_LIGAND_ASL", "") or "").strip()
    try:
        from schrodinger.structure import StructureReader
        from schrodinger.structutils import analyze
        fs = max(StructureReader(str(cms_file)), key=lambda c: c.atom_total)
    except Exception:
        return _cfg_asl   # cannot introspect the complex → trust the configured ASL
    if _cfg_asl:
        try:
            if len(analyze.evaluate_asl(fs, _cfg_asl)) > 0:
                return _cfg_asl
        except Exception:
            pass
    _min_atoms = int(getattr(CFG, "LIGAND_MIN_ATOMS", 5))
    _std = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS HID HIE HIP ILE LEU LYS "
               "MET PHE PRO SER THR TRP TYR VAL".split())
    _solvent = {"T3P", "SPC", "HOH", "WAT", "NA", "CL", "K", "POT", "SOD", "CLA", "NA+", "CL-"}
    _cands = []
    for mol in fs.molecule:
        _resns = {a.pdbres.strip() for a in mol.atom}
        if (_resns & _std) or (_resns <= _solvent):
            continue
        if mol.atom_total >= _min_atoms:
            _cands.append(mol)
    if _cands:
        _lig = min(_cands, key=lambda m: m.atom_total)
        _rn = next(iter({a.pdbres.strip() for a in _lig.atom}))
        _echo(f"    [!] MM-GBSA: '{_cfg_asl}' matched no atoms - pinning smallest non-solvent "
              f"ligand '{_rn}' ({_lig.atom_total} atoms) via 'mol.num {_lig.number}'.")
        return f"mol.num {_lig.number}"
    _echo(f"    [!] MM-GBSA: '{_cfg_asl}' matched no atoms and no fallback ligand "
          f"≥{_min_atoms} atoms found - reverting to thermal_mmgbsa auto-detect.")
    return ""


def run_mmgbsa(job_dir: Path, job_name: str, rank: str) -> Path | None:
    """Run thermal_mmgbsa.py on <job>-out.cms inside its MD folder. Idempotent:
    returns the existing CSV when already computed. Returns the results CSV path
    or None on failure."""
    # The disable gate lives HERE, not only in run_mmgbsa_phase: the sequential MD loop calls
    # run_mmgbsa() directly, so gating only in the phase wrapper would run the whole Prime batch and
    # then print "disabled - skipped".
    if not getattr(CFG, "MMGBSA_RUN", True):
        _echo("    ✘ Skipped      : MM-GBSA disabled (CFG.MMGBSA_RUN = False).")
        return None
    cms_file = job_dir / f"{job_name}-out.cms"
    if not cms_file.is_file():
        _echo(f"    ✘ Skipped      : no {cms_file.name} - the MD simulation has not finished.")
        return None
    existing = _mmgbsa_csv(job_dir, job_name)
    if existing is not None and _mmgbsa_complete(existing, job_dir, job_name):
        """
        A CSV written before frame stamping existed carries no `Frame` column, and Step 07 must not
        infer the frame from the row position (row i is frame i·step under a stride). Stamp it in
        place - the frame list is reconstructible from the shard plan - so an already-scored run is
        brought up to the current format without re-scoring anything. A 100%-complete CSV here means
        MM-GBSA never re-runs; any leftover shard scratch is tidied away.
        """
        _retrofit_frame_stamps(existing, job_dir, job_name, rank)
        _echo(f"    ✔ Already done : reusing {existing.name} (delete it to force a re-score).")
        _cleanup_mmgbsa_shards(job_dir, job_name, rank)
        return existing
    if existing is not None:
        _echo(f"    ⚠ Partial CSV  : {existing.name} has fewer rows than the trajectory expects "
              f"- re-scoring to complete it.")

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
    (SCHRODINGER_TMPDIR if the environment sets it, else /tmp - the job server's default
    scratch location).
    """
    try:
        _probe = os.environ.get("SCHRODINGER_TMPDIR") or "/tmp"
        if not os.path.isdir(_probe):
            _probe = "/tmp"
        _free_gb = _sh.disk_usage(_probe).free / 2**30
        _gb_per_subjob = float(getattr(CFG, "MMGBSA_SCRATCH_GB_PER_SUBJOB", 22))
        _headroom = float(getattr(CFG, "MMGBSA_SCRATCH_HEADROOM_FRAC", 0.75))
        _cap = int(getattr(CFG, "MMGBSA_MAX_NJOBS", 0)) or _ncpu   # 0 = no ceiling beyond GLOBAL_MAX_WORKERS
        _disk_cap = max(2, min(_cap, int(_headroom * _free_gb / _gb_per_subjob)))
        if _disk_cap < _ncpu:
            _echo(f"    Scratch disk   : ~{_free_gb:,.0f} GB free - only enough for {_disk_cap} subjobs "
                  f"(~{_gb_per_subjob:.0f} GB each), so Prime is reduced {_ncpu} → {_disk_cap}. "
                  f"Running more would fill the disk and kill the job.")
            _ncpu = _disk_cap
        else:
            _echo(f"    Scratch disk   : ~{_free_gb:,.0f} GB free - room for {_ncpu} Prime subjobs "
                  f"(~{_gb_per_subjob:.0f} GB each)")
    except Exception:
        pass

    _lig_asl_cfg = _resolve_mmgbsa_lig_asl(cms_file)   # size-independent LIG resolve + ≥5-atom fallback
    _step_cfg = int(getattr(CFG, "MMGBSA_STEP_SIZE", 0) or 0)
    _interval = max(5, int(getattr(CFG, "MMGBSA_PROGRESS_INTERVAL_SEC", 30)))
    # Ground-truth frame count (same source the SID phase uses).
    _total = traj_frame_count(job_dir / f"{job_name}_trj") or _DEFAULT_FRAME_TOTAL

    """
    Sharded path (default): concurrent frame-range shards, so the trajectory read is
    parallel and overlaps Prime instead of blocking it on one core. Same frames, same
    scores - see run_mmgbsa_sharded. CFG.MMGBSA_SHARD_FRAMES = 0 keeps the plain single
    serial thermal_mmgbsa run below.
    """
    _shard_frames = int(getattr(CFG, "MMGBSA_SHARD_FRAMES", 0) or 0)
    if _shard_frames > 0 and _total > _shard_frames:
        return run_mmgbsa_sharded(job_dir, job_name, rank, cms_file, _lig_asl_cfg,
                                  _step_cfg, _total, _ncpu, _interval)

    cmd = [SCHROD_RUN, "thermal_mmgbsa.py", cms_file.name,
           "-j", f"{job_name}_mmgbsa", "-NJOBS", str(_ncpu)]
    # Pin the ligand explicitly via thermal_mmgbsa's -lig_asl flag so Prime scores the substrate;
    # a small or heavily-fluorinated ligand is otherwise misassigned by auto-detection. The ASL was
    # resolved once above (configured name if it matches, else the ≥5-atom fallback, else empty).
    _lig_asl = _lig_asl_cfg
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
                                _echo(f"  [Rank {rank}] MM-GBSA timed out - skipped (see {log.name}).")
                                return None
                else:
                    proc.wait()
        # thermal_mmgbsa -NJOBS should block until its Prime subjobs finish; if a build submits them
        # async and returns early, wait for the job server to drain this rank's batch before returning,
        # so the next rank's MD never shares the scratch disk with a still-running Prime batch.
        _await_mmgbsa_jobserver(f"{job_name}_mmgbsa")
        if proc.returncode != 0:
            _echo(f"  [Rank {rank}] MM-GBSA exited rc={proc.returncode} - see {log.name}.")
            _diag = _diagnose_mmgbsa_failure(job_dir, job_name)
            if _diag:
                _echo(f"  [Rank {rank}] ↳ cause: {_diag}")
            return None
    except Exception as e:
        _echo(f"  [Rank {rank}] MM-GBSA failed ({e}) - see {log.name}.")
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
                _stamped = _stamp_frames(_df, 0, _total, max(1, _step_cfg), f"Rank {rank}")
                _tmp = _csv.with_suffix(".csv.tmp")
                _stamped.to_csv(_tmp, index=False)
                _tmp.replace(_csv)                 # atomic: a kill mid-write cannot corrupt the CSV
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
    Under a stride, row 150 is frame 1,500 - plotting against the row index would compress a
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

    "mean" (default) is the arithmetic ensemble average - the standard thermal
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
    hit = (_utils_mod.latest_by_mtime(prod.glob(CFG.GLOB_RANKED_CSV))
           or _utils_mod.latest_by_mtime(prod.glob("*Ranked*.csv"))) if prod.is_dir() else None
    if hit is None:
        return {}
    try:
        df = pd.read_csv(hit, usecols=["Scientific_Rank", "degrader_tier"])
        return {int(r): str(t) for r, t in zip(df["Scientific_Rank"], df["degrader_tier"])}
    except Exception:
        return {}


def _lookup_ligands(run_root: Path) -> dict:
    """Map Scientific_Rank → Ligand_Name from the Step 02 ranked CSV, so the figures name the
    PFAS species rather than only its rank. Returns {} if unavailable (labels fall back to
    the bare rank)."""
    prod = run_root / "1_Boltz2_Production"
    hit = (_utils_mod.latest_by_mtime(prod.glob(CFG.GLOB_RANKED_CSV))
           or _utils_mod.latest_by_mtime(prod.glob("*Ranked*.csv"))) if prod.is_dir() else None
    if hit is None:
        return {}
    try:
        df = pd.read_csv(hit, usecols=["Scientific_Rank", "Ligand_Name"])
        return {int(r): str(l) for r, l in zip(df["Scientific_Rank"], df["Ligand_Name"])}
    except Exception:
        return {}


def _short_ligand(name) -> str:
    """The ligand's short name (FA / DFA / TFA) from CFG - the same abbreviations every figure uses.

    A ligand absent from the map keeps its full name rather than being silently mangled.
    """
    _n = re.sub(r"^\d+_", "", str(name or "")).replace("_", " ").strip()
    return CFG.VIS_LIGAND_SHORT.get(_n.lower(), _n)


_CTRL_LABEL = f"{CFG.REFERENCE_PDB_ID}-{_short_ligand(CFG.CONTROL_MD_LIGANDS[0])}"   # distinct label for the reference × control-ligand positive control (e.g. 3R3U-FA), derived from CFG
_CTRL_COLOUR = CFG.VIS_ACCENT["control"]    # one distinct colour for the control across every 06/07 figure
_QC_SEP_KW = dict(color=CFG.VIS_INK["mid"], linestyle=":", alpha=0.5, linewidth=1.0, zorder=0)   # per-column separator, matching the Step 03 figures


def _qc_group_seps(axobj, n: int) -> None:
    """Vertical dotted separators between the n categorical ligand groups (lines at 0.5 … n-1.5),
    the same per-column look the Step 03 figures use, so group boundaries read consistently."""
    for _i in range(int(n) - 1):
        axobj.axvline(_i + 0.5, **_QC_SEP_KW)


def _lookup_controls(run_root: Path) -> set:
    """Scientific_Ranks that are the 3R3U control (is_control, or job_name reserved index
    CONTROL_JOB_PREFIX). Lets every physics figure paint the control distinctly and label it
    '3R3U-FA' instead of a bare 'FA' shared with the candidate fluoroacetate. Empty if unavailable."""
    prod = run_root / "1_Boltz2_Production"
    hit = (_utils_mod.latest_by_mtime(prod.glob(CFG.GLOB_RANKED_CSV))
           or _utils_mod.latest_by_mtime(prod.glob("*Ranked*.csv"))) if prod.is_dir() else None
    if hit is None:
        return set()
    try:
        df = pd.read_csv(hit, usecols=lambda c: c in {"Scientific_Rank", "is_control", "job_name"})
        _pref = str(getattr(CFG, "CONTROL_JOB_PREFIX", "0000000"))
        _is_ctrl = pd.Series(False, index=df.index)
        if "is_control" in df.columns:
            _is_ctrl |= df["is_control"].astype(str).str.strip().str.lower().isin(["true", "1", "1.0", "yes"])
        if "job_name" in df.columns:
            _is_ctrl |= df["job_name"].astype(str).str.startswith(_pref)
        return {int(r) for r in df.loc[_is_ctrl, "Scientific_Rank"]}
    except Exception:
        return set()


def _ctrl_label(rk, ligands: dict, controls: set) -> str:
    """'3R3U-FA' for a control rank, else the short PFAS name (fallback R_<rk>)."""
    if rk in controls:
        return _CTRL_LABEL
    return _short_ligand(ligands.get(rk, "")) or f"R_{rk}"


def _ctrl_palette(ranks, controls: set, base_palette) -> list:
    """Per-rank colours: the control gets the one distinct control colour, candidates cycle the
    base palette (skipping the control so its colour is never reused for a candidate)."""
    out, _bi = [], 0
    for rk in ranks:
        if rk in controls:
            out.append(_CTRL_COLOUR)
        else:
            out.append(base_palette[_bi % len(base_palette)]); _bi += 1
    return out


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


def plot_mmgbsa_individual(out_dir: Path, job_name: str, rank: str, dg: "pd.Series",
                           ns_per_frame: float = 0.0, quiet: bool = False) -> None:
    """Per-job MM-GBSA: ΔG_bind against simulation time, plus its distribution.

    Two things this figure must get right:
      • The x-axis is the SIMULATION, not the row number. The series is indexed by trajectory
        frame, so a strided run still spans the full trajectory (with gaps), and 'when' an
        event happened is readable.
      • A handful of failed minimisations (ΔG ≈ −1000 kcal/mol) must not set the y-scale, or
        the 99%+ of frames that carry the actual signal collapse onto a flat line. The axis is
        therefore scaled to the robust core and the failures are reported explicitly at the
        panel edge - flagged, never silently dropped.
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
    Y-window from the CORE ensemble (failures already excluded), not from the full series -
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
    wherever a failure sits - an artefact of the failure, not a change in binding.
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
    # Single-row legend inside the panel: one column per entry so it never wraps to a second row;
    # tight column/handle spacing keeps it narrow enough to sit in the sparse lower-left corner.
    _lh, _ll = ax1.get_legend_handles_labels()
    ax1.legend(_lh, _ll, frameon=True, ncol=max(1, len(_lh)), loc="lower left",
               columnspacing=1.0, handletextpad=0.5)
    ax1.grid(alpha=0.25, linewidth=0.5)

    # ── Panel 2: distribution of the core ensemble (shares the y-axis) ────────────────────
    ax2.hist(_core.values, bins=45, range=(_ylo, _yhi), color=_INK["hist"], alpha=0.85,
             orientation="horizontal", zorder=2)
    ax2.axhline(_med, color=_INK["mean"], linestyle="--", linewidth=1.6, zorder=3)
    ax2.axhspan(dg.quantile(0.25), dg.quantile(0.75), color=_INK["mean"], alpha=0.10, zorder=1)
    ax2.set_xlabel(f"Frames scored\n(n = {len(dg):,}; IQR shaded)")
    # Same y gridlines and tick marks as the time panel (labels suppressed, since the axis is
    # shared) - so a ΔG value can be read straight across from one panel to the other.
    ax2.tick_params(labelleft=False, left=True)
    ax2.grid(alpha=0.25, linewidth=0.5)
    out_path = out_dir / f"MMGBSA_Profile_R{rank}.svg"
    plt.savefig(out_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    if not quiet:
        _echo(f"    ✔ Figure       : {_utils_mod.deflx_fig_name(out_path.name)}")


def plot_mmgbsa_combined(out_dir: Path, per_job: list,
                         ligands: "dict | None" = None, nspf: "dict | None" = None,
                         controls: "set | None" = None) -> None:
    """Compare ΔG_bind across the ranks, three ways.

      left   - the full distribution each summary is drawn from (violin + IQR box), with the
               mean tracked alongside the median so the gap between them exposes which ensembles
               are skewed by failed minimisations.
      middle - the cumulative distributions, which show HOW the ensembles differ (a shift in
               location versus a difference in spread) rather than only THAT they differ.
      right  - the ranking, as the median with a moving-block bootstrap confidence interval,
               plus Cliff's delta against the tightest binder.

    Statistics note: MD frames are autocorrelated, so a t-test or Mann-Whitney over ~10⁴ frames
    would report p ≈ 0 for any difference at all and mean nothing. The bootstrap resamples
    contiguous blocks (CFG.MMGBSA_BOOTSTRAP_BLOCK) so each block counts as one independent draw.
    Colours come from CFG.MMGBSA_RANK_PALETTE.
    """
    per_job = [(r, dg) for r, dg in per_job if not dg.empty]
    if not per_job:
        _echo("  [!] MM-GBSA combined skipped - no parsed ΔG_bind series.")
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    # Consistent order across every combined figure: non-controls by rank, the 3R3U control LAST.
    per_job.sort(key=lambda x: (int(x[0]) in controls, int(x[0])))
    ligands, nspf, controls = ligands or {}, nspf or {}, controls or set()

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
        lig = _ctrl_label(int(r), ligands, controls)
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
    ax1 = fig.add_subplot(gs[0, 0])     # distribution (top-left)
    ax3 = fig.add_subplot(gs[0, 1])     # ranking (top-right)
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
                     textcoords="data", va="center", ha="left", fontsize=CFG.VIS_FONT_AXIS_LABEL,
                     fontweight="bold", color=_cols[i], zorder=9, path_effects=_stroke)
        if n:
            ax1.text(_x[i], _lo - _pad * 0.55, f"{n} outlier{'s' if n > 1 else ''}",
                     ha="center", fontsize=CFG.VIS_FONT_LEGEND, color=_cols[i], fontweight="bold",
                     zorder=9, path_effects=_stroke)
    ax1.set_xticks(_x); ax1.set_xticklabels(_labels)
    ax1.set_ylabel("ΔG$_{bind}$ (kcal/mol)  ·  per-frame Prime MM-GBSA")
    ax1.set_ylim(_lo - _pad, _hi + _pad)
    """
    One legend, inside, lower right - no caption under the panel. The violin needs no entry: it IS
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
    # Top-left, single column (one row per entry).
    ax1.legend(handles=_hdlA, loc="upper left", frameon=True, ncol=1)
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
                 fontsize=CFG.VIS_FONT_AXIS_LABEL, fontweight="bold", color=_INK["dark"], zorder=9,
                 path_effects=[pe.withStroke(linewidth=3.0, foreground=_INK["light"])])
    # Effect size vs the tightest binder - the number that says whether a gap MATTERS.
    """
    Effect size, reported PAIRWISE. There is no control or reference candidate here - these are
    three substrates compared against each other - so singling one out as a baseline would
    invent a hierarchy the design does not have. Cliff's delta is computed for every pair
    instead: δ(X vs Y) > 0 means X's frames sit at HIGHER (weaker-binding) ΔG than Y's.
    """
    _pairs = []
    for i in range(len(per_job)):
        for j in range(i + 1, len(per_job)):
            d = _cliffs_delta(per_job[i][1], per_job[j][1])
            _pairs.append(f"#{per_job[i][0]} vs #{per_job[j][0]}: δ = {d:+.2f} ({_delta_word(d)})")
    # The pairwise deltas are the legend's entries (one pair per line) under the "Cliff's δ" title -
    # inside the panel, never a caption below it, and narrow enough to sit clear of the tallest bar.
    ax3.axhline(0.0, color=_INK["outline"], linewidth=0.8)
    ax3.set_xticks(_x); ax3.set_xticklabels(_labels)
    ax3.set_ylabel("Median ΔG$_{bind}$ (kcal/mol)")
    """
    The marks are labelled ON the bars rather than in a legend: the median is already printed inside
    each bar, and the two error bars are named once, in place, on the first bar - a legend for two
    marks that are visible at a glance is wasted panel space. Only the pairwise effect sizes, which
    have nowhere natural to sit, remain as a boxed note.
    """
    _lo_ci0, _hi_ci0 = _cis[0]
    # Named in place, INSIDE the bar, beside the mark each label refers to - no leader lines to
    # cross a neighbouring bar, and no legend entry for something already visible.
    """
    Each label sits beside the mark it names. They are BLACK on a white stroke, so they read the same
    whether they land on the coloured bar or on the panel behind it - a single colour that depends on
    the background would vanish on one of the two.
    """
    _stroke_b = [pe.withStroke(linewidth=2.6, foreground=_INK["light"])]
    ax3.text(_x[0] + 0.16, _q1[0], "IQR", fontsize=CFG.VIS_FONT_LEGEND, color=_INK["dark"], fontweight="bold",
             va="center", ha="left", zorder=9, path_effects=_stroke_b)
    if _lo_ci0 == _lo_ci0:
        ax3.text(_x[0] + 0.16, (_lo_ci0 + _hi_ci0) / 2,
                 f"{getattr(CFG, 'MMGBSA_BOOTSTRAP_CI', 95):.0f}% CI", fontsize=CFG.VIS_FONT_LEGEND,
                 color=_INK["dark"], fontweight="bold", va="center", ha="left", zorder=9,
                 path_effects=_stroke_b)
    if _pairs:
        # Two-column, text-only legend (blank handles): the pairwise deltas read as a compact grid
        # rather than a tall single column. Header is the legend title.
        _dh3 = [Line2D([], [], color="none") for _ in _pairs]
        _leg3 = ax3.legend(_dh3, _pairs, loc="upper left", ncol=2,
                           title="Cliff's δ - δ > 0: the first binds more weakly",
                           fontsize=CFG.VIS_FONT_ANNOT, title_fontsize=CFG.VIS_FONT_ANNOT,
                           frameon=True, handlelength=0, handletextpad=0,
                           columnspacing=1.2, labelspacing=0.3, borderpad=0.5)
        _leg3.get_frame().set_facecolor(_INK["light"]); _leg3.get_frame().set_edgecolor(_INK["faint"])
        _leg3.get_frame().set_alpha(0.93); _leg3.get_frame().set_linewidth(0.8)
    ax3.grid(alpha=0.25, linewidth=0.5, axis="y")
    ax3.invert_yaxis()

    # ── Panel 4: time course + cumulative distribution, on one shared ΔG axis ───────────
    _draw_time_cumulative(ax4, per_job, _cols, ligands, nspf, controls)

    """
    Publication finish: label the panels A/B/C (a reader cites 'panel B', not 'the top-right one'),
    drop the top/right spines that carry no data, and give every axes the same tick geometry. The
    time+cumulative panel keeps its top spine - a second x-axis lives there.
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

    out_path = out_dir / "04_MMGBSA_Combined_AllRanks.svg"
    plt.savefig(out_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    _echo(f"  MM-GBSA combined figure saved: {_utils_mod.deflx_fig_name(out_path.resolve())}")


def _effective_n(dg: "pd.Series") -> float:
    """Effective number of INDEPENDENT samples in an autocorrelated series.

    N_eff = N / (1 + 2·Σρ_k), summing the autocorrelation until it first goes negative
    (the standard initial-positive-sequence cut-off). This is the number that governs how
    precisely ⟨ΔG_bind⟩ is known - not the raw frame count, which merely counts how often the
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


def _draw_time_cumulative(axT, per_job: list, cols: list, ligands: dict, nspf: dict,
                          controls: "set | None" = None) -> None:
    """Draw the time course and the cumulative distribution into ONE axes, on a shared ΔG axis.

      • BOTTOM x-axis - simulation time, in CFG.MMGBSA_TIME_WINDOW_NS windows. Solid line with
        markers = window median; band = window IQR; pale band = the candidate's overall IQR;
        ▼ = a window whose minimisations blew up.
      • TOP x-axis - cumulative fraction of frames (0→1). Dashed line = the candidate's cumulative
        curve, against the SAME ΔG axis. Where it is steep, that ΔG level is heavily populated -
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
        return f"#{r} {_ctrl_label(int(r), ligands, controls)}".strip()
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
    axT.set_xlabel(f"Simulation time (ns), in {_win:.0f} ns windows - solid line, markers"
                   f"        [top axis: cumulative fraction of frames - dashed line]")
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
    axT.legend(handles=hdl, loc="lower right", frameon=True,
               ncol=max(1, math.ceil(len(hdl) / _rows)),

               title="Relative ranking only - GB overstabilises anionic PFAS: compare candidates, "
                     "never absolute values",
               )


# ── 8.5  Phase driver ────────────────────────────────────────────────────────
def _plot_rank_mmgbsa(md_dir: Path, job_dir: Path, job_name: str, rank: str, csv: Path, quiet: bool = False) -> None:
    """Draw ONE rank's MM-GBSA profile from an existing CSV - never re-scores.

    A rank is finished the moment its own MD → SID → MM-GBSA has landed, so its profile is readable
    then; making it wait for the cross-rank Finalise pass would hold a completed rank's figure
    hostage to every other rank's MD. Only the COMBINED plots, which genuinely need all ranks,
    belong in Finalise. Both callers route through here so the output path is defined once.
    """
    dg = _mmgbsa_dg_series(csv)
    if dg.empty:
        return
    _trj = job_dir / f"{job_name}_trj"
    _span, _nfr = traj_span_ns(_trj), traj_frame_count(_trj)
    _nspf = (_span / (_nfr - 1)) if (_span > 0 and _nfr > 1) else 0.0
    _out = _analysis_dir(md_dir.parent) / getattr(CFG, "MMGBSA_OUTPUT_SUBDIR", "Prime-MMGBSA")
    _out.mkdir(parents=True, exist_ok=True)
    plot_mmgbsa_individual(_out, job_name, rank, dg, ns_per_frame=_nspf, quiet=quiet)


def run_mmgbsa_phase(md_dir: Path, run_root: Path, plots_only: bool = False) -> str:
    """Run + plot MM-GBSA for every completed MD job (idempotent).

    Returns a status the caller maps to the pipeline step result:
      • "ok"   - every eligible job produced a ΔG_bind CSV (or none eligible);
      • "warn" - the phase ran but ≥1 job's MM-GBSA failed (e.g. Prime rc=1).
    "warn" lets the pipeline report WARN (not a false PASS) without hard-aborting:
    MM-GBSA is complementary to the QSite barrier, so Step 07 still runs.
    "disabled" (CFG.MMGBSA_RUN=False) is an intentional no-op → "ok"."""
    if not getattr(CFG, "MMGBSA_RUN", True):
        _echo("  MM-GBSA disabled (CFG.MMGBSA_RUN = False) - skipped.")
        return "ok"
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
    _echo(f"Prime MM-GBSA - rescoring binding free energy for {len(job_dirs)} completed MD job(s)")
    _echo(_SEP)
    _echo(f"  {_C.DIM}Note: GB implicit solvent overstabilises anionic PFAS - compare ΔG_bind BETWEEN "
          f"ranks, never as an absolute affinity.{_C.ENDC}")
    # MM-GBSA data (the summary CSV) stays in its own compute subdir; the FIGURES go to the single
    # 06_Analysis folder alongside every other Step-06 figure. md_dir.parent is 6_Physics_Validation.
    fig_dir = _analysis_dir(md_dir.parent)                                  # 06_Analysis (all Step-06 figures)
    out_dir = fig_dir / getattr(CFG, "MMGBSA_OUTPUT_SUBDIR", "Prime-MMGBSA")  # one flat profile per rank lands here
    tiers = _lookup_tiers(run_root)
    ligands = _lookup_ligands(run_root)
    controls = _lookup_controls(run_root)   # ranks labelled with the CFG control tag (e.g. 3R3U-FA)
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
        # plots_only (the Finalise pass) must NEVER re-run Prime: a rank whose MM-GBSA failed or timed
        # out has no CSV, and run_mmgbsa would restart the multi-hour batch. Read the existing CSV only.
        csv = _mmgbsa_csv(d, job_name) if plots_only else run_mmgbsa(d, job_name, rank)
        if csv is None:
            _failed += 1
            _echo(f"    ✘ MM-GBSA did not complete for Rank {rank} - see the messages above.")
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
                unbinding events, not numerical failures - calling them 'blown-up' and quoting
                their minimum would misdescribe them.
                """
                _q1, _q3 = dg.quantile(.25), dg.quantile(.75)
                _lo_tail = _bad[_bad < _q1]
                _hi_tail = _bad[_bad > _q3]
                _parts = []
                if len(_lo_tail):
                    _parts.append(f"{len(_lo_tail)} unphysically strong (to {float(_lo_tail.min()):.0f} "
                                  f"kcal/mol - blown-up minimisations)")
                if len(_hi_tail):
                    _parts.append(f"{len(_hi_tail)} weakly bound (to {float(_hi_tail.max()):.0f} "
                                  f"kcal/mol - ligand loosely held, not a numerical failure)")
                _echo(f"    Outlier frames : {len(_bad)} of {dg.size:,} "
                      f"({100 * len(_bad) / dg.size:.2f}%) - {'; '.join(_parts)}. "
                      f"They shift the arithmetic mean ({float(dg.mean()):.2f}), so the median "
                      f"above is the honest estimate.")
                _wins = _failure_windows(_bad, _nspf)
                if _wins:
                    _echo(f"    ↳ when         : {_wins}. Failures that CLUSTER in time point at "
                          f"unstable stretches of the trajectory (inspect those frames); failures "
                          f"scattered evenly are ordinary Prime convergence noise.")
        try:
            _plot_rank_mmgbsa(md_dir, d, job_name, rank, csv)
        except Exception as e:
            _echo(f"    ✘ per-job plot failed ({e}) - skipped.")
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
        _csv_out = fig_dir / CFG.FILE_MMGBSA_SUMMARY
        _utils_mod.atomic_write_csv(_sdf, _csv_out)
        _echo(f"  MM-GBSA summary table saved : {_csv_out.resolve()}")
    try:
        plot_mmgbsa_combined(fig_dir, per_job, ligands, nspf, controls)   # combined at 06_Analysis root
    except Exception as e:
        _echo(f"  [!] MM-GBSA combined plot failed ({e}) - skipped.")
    if _failed:
        _echo(f"  [!] MM-GBSA: {_failed} of {len(job_dirs)} job(s) failed "
              f"(see per-rank rc/cause above) - step will report WARN, not PASS.")
        return "warn"
    return "ok"


# =============================================================================
# SECTION 8b: DEFLUORINATION GEOMETRY (native per-rank trajectory analysis)
# =============================================================================
"""SN2-defluorination geometry per rank, read straight from the Desmond cms + _trj
with Schrodinger's own traj API (no format conversion, no external MD toolkit). The
eight catalytic residues come from the ranked CSV per rank (their numbering shifts
between ranks); the ligand's alpha-carbon, leaving fluoride and carboxylate oxygens
are auto-detected from its bond graph, so any carboxylate PFAS works. Writes three
figures + one merged CSV + a log + descriptions to Defluorination/Defluorination_R{N}/.
"""
# All colours from CFG (single source of truth); order preserved for _DEFL_OKABE[i] usage.
_DEFL_OKABE = [CFG.VIS_ACCENT["blue"], CFG.VIS_ACCENT["amber"], CFG.VIS_ACCENT["green"],
               CFG.VIS_ACCENT["magenta"], CFG.VIS_ACCENT["vermillion"], CFG.VIS_ACCENT["sky"],
               CFG.VIS_ACCENT["yellow"], CFG.VIS_INK["near_black"]]
_DEFL_GREEN = CFG.VIS_ACCENT["green"]        # NAC shading / reactive markers
_DEFL_GREY = CFG.VIS_INK["faint"]            # cutoff lines, non-reactive points, "other" bars
_DEFL_MUTE = CFG.VIS_INK["mid"]              # zero-reference line
_DEFL_GRID_X = CFG.VIS_INK["faint"]          # neutral grey gridlines (never a data-line colour)
_DEFL_GRID_Y = CFG.VIS_INK["faint"]          # neutral grey gridlines (never a data-line colour)
_DEFL_AXTXT = CFG.VIS_INK["near_black"]       # axis labels + tick values (black)
_DEFL_TEXT_MAX_LUM = float(getattr(CFG, "VIS_TEXT_MAX_LUM", 0.42))   # bar colour used as TEXT darkened to at most this luminance (CFG SSOT)


def _readable(colour, max_lum: float = _DEFL_TEXT_MAX_LUM) -> str:
    """A bar's own colour, darkened toward its hue until it reads on white - so the same colour can name
    the bar's y-label AND its value without a light hue (yellow, sky, amber) vanishing. Colours already
    dark enough (blue, vermillion, green, near-black) are returned unchanged. Relative luminance is
    linear in a uniform RGB scale, so one factor both lowers it to the target and preserves the hue."""
    r, g, b = _mcolors.to_rgb(colour)
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    if lum <= max_lum or lum == 0.0:
        return _mcolors.to_hex((r, g, b))
    k = max_lum / lum
    return _mcolors.to_hex((r * k, g * k, b * k))
_DEFL_NAC_DIST = float(CFG.DEFLUOR_NAC_DIST_A)      # Od...C(alpha) near-attack distance
_DEFL_NAC_ANGLE = float(CFG.DEFLUOR_NAC_ANGLE_DEG)  # Od-C(alpha)-F in-line attack angle
_DEFL_ENGAGE = float(getattr(CFG, "DEFLUOR_ENGAGE_A", 4.0))          # residue engaged with ligand within this
_DEFL_POCKET = float(getattr(CFG, "DEFLUOR_POCKET_RADIUS_A", 8.0))   # frame-0 pocket radius for the COM reference


def _defl_resid(mapped) -> int:
    return int("".join(c for c in str(mapped) if c.isdigit()) or 0)


def _defl_min_image(dvec, box):
    """Orthorhombic minimum-image on a displacement array (n,3)."""
    L = np.array([box[0][0], box[1][1], box[2][2]], float)
    L[L <= 0] = 1e9
    return dvec - L * np.round(dvec / L)


def _defl_ligand_atoms(fs):
    """(alpha_C aid, [alpha-C F aids], [carboxylate O aids]) from the LIG bond graph.

    Every fluorine on the scissile alpha-carbon is returned, not just one: on a poly-fluorinated
    carbon (CF2 in difluoroacetate, CF3 in trifluoroacetate) the fluorine aligned for backside SN2
    attack rotates frame-to-frame, so the leaving fluoride must be chosen per frame from the whole
    set rather than fixed here."""
    lig = [a for a in fs.atom if a.pdbres.strip() == "LIG" and a.element.strip() != "H"]
    ele = {int(a): a.element.strip() for a in lig}
    nb = {int(a): [] for a in lig}
    for b in fs.bond:
        i, j = int(b.atom1), int(b.atom2)
        if i in nb and j in nb:
            nb[i].append(j); nb[j].append(i)
    cC = next((i for i in ele if ele[i] == "C" and sum(ele[k] == "O" for k in nb[i]) >= 2), None)
    if cC is None:
        raise RuntimeError("no carboxylate carbon in ligand")
    cox = [k for k in nb[cC] if ele[k] == "O"]
    cn = [k for k in nb[cC] if ele[k] == "C"]
    aC = next((c for c in cn if any(ele[k] == "F" for k in nb[c])), cn[0] if cn else None)
    if aC is None:
        raise RuntimeError("no alpha-carbon in ligand")
    Fs = [k for k in nb[aC] if ele[k] == "F"]
    if not Fs:
        raise RuntimeError("no fluorine on the alpha-carbon (not a scissile C-F ligand)")
    return aC, Fs, cox


def _defl_mapped_residues(run_root: Path, rank: str) -> dict:
    """Per-rank catalytic residue numbers from the ranked CSV (Mapped_* columns)."""
    try:
        rk = pd.read_csv(newest_ranked_csv(run_root), low_memory=False)
        m = rk[rk["Scientific_Rank"].astype(str) == str(rank)]
        if not len(m):
            return {}
        r = m.iloc[0]
        return {k: _defl_resid(r.get(c)) for k, c in (
            ("nuc", "Mapped_Nucleophile"), ("base", "Mapped_Base"), ("acid", "Mapped_Acid"),
            ("c1", "Mapped_Clamp1"), ("c2", "Mapped_Clamp2"), ("sh", "Mapped_Stabiliser_H"),
            ("sw", "Mapped_Stabiliser_W"), ("sy", "Mapped_Stabiliser_Y"))}
    except Exception:
        return {}


def run_defluorination(job_dir: Path, job_name: str, rank: str, md_dir: Path,
                       run_root: Path, md_ns: float) -> "Path | None":
    """SN2-defluorination geometry for one rank, native from cms + _trj."""
    if not getattr(CFG, "DEFLUOR_RUN", True):
        _echo("    ✘ Skipped      : defluorination disabled (CFG.DEFLUOR_RUN = False).")
        return None
    cms_file = job_dir / f"{job_name}-out.cms"
    trj = job_dir / f"{job_name}_trj"
    if not cms_file.is_file() or not trj.is_dir():
        _echo("    ✘ Skipped      : no cms/_trj for defluorination.")
        return None
    R = _defl_mapped_residues(run_root, rank)
    if not R.get("nuc"):
        _echo("    ✘ Skipped      : no mapped catalytic residues for this rank.")
        return None
    out = _analysis_dir(md_dir.parent) / getattr(CFG, "DEFLUOR_OUTPUT_SUBDIR", "Defluorination") / f"Defluorination_R{rank}"
    out.mkdir(parents=True, exist_ok=True)
    _echo("")
    _echo(f"  Defluorination geometry - Rank {rank}  ·  SN2 attack pose · NAC · fluoride cradle · "
          f"carboxylate clamp · MM-GBSA drivers  (native from cms + _trj)")

    # ── atom selections as trajectory gids ────────────────────────────────────────────
    # Catalytic residues come from the ranked CSV (per rank); the ligand's alpha-carbon,
    # leaving fluoride and carboxylate oxygens are auto-detected from its bond graph. Every
    # selection is mapped to trajectory gids so positions can be read straight from a frame.
    from schrodinger.application.desmond.packages import topo, traj
    msys, cms = topo.read_cms(str(cms_file))
    fs = cms.fsys_ct
    n_real = fs.atom_total
    aid2gid = dict(zip(range(1, n_real + 1), topo.aids2gids(cms, list(range(1, n_real + 1)))))

    def G(resnums, names):
        rs, ns = set(resnums), set(names)
        return [aid2gid[int(a)] for a in fs.atom
                if a.resnum in rs and a.pdbname.strip() in ns and int(a) in aid2gid]

    aC_aid, F_aids, cox_aids = _defl_ligand_atoms(fs)
    gC = aid2gid[aC_aid]; gFs = [aid2gid[i] for i in F_aids]
    gCox = [aid2gid[i] for i in cox_aids]
    gOd = G([R["nuc"]], {"OD1", "OD2", "OE1", "OE2"})
    _STD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS HID HIE HIP ILE LEU LYS MET PHE "
               "PRO SER THR TRP TYR VAL".split())
    gCA = [aid2gid[int(a)] for a in fs.atom
           if a.pdbname.strip() == "CA" and a.pdbres.strip() in _STD and int(a) in aid2gid]
    gLig = [aid2gid[int(a)] for a in fs.atom
            if a.pdbres.strip() == "LIG" and a.element.strip() != "H" and int(a) in aid2gid]
    if len(gOd) == 0 or len(gLig) == 0:
        _echo("    ✘ Skipped      : nucleophile has no OD/OE atoms, or no ligand heavy atoms, in the prepared cms.")
        return None
    cradle = {f"HIS{R['sh']}": G([R["sh"]], {"ND1", "NE2"}),
              f"TRP{R['sw']}": G([R["sw"]], {"NE1"}),
              f"TYR{R['sy']}": G([R["sy"]], {"OH"}),
              f"ARG{R['c1']}": G([R["c1"]], {"NH1", "NH2", "NE"}),
              f"ARG{R['c2']}": G([R["c2"]], {"NH1", "NH2", "NE"})}
    clamp = {f"ARG{R['c1']}": cradle[f"ARG{R['c1']}"], f"ARG{R['c2']}": cradle[f"ARG{R['c2']}"]}
    res8 = {f"ASP{R['nuc']}·Nu": [R["nuc"]], f"HIS{R['base']}·base": [R["base"]],
            f"ASP{R['acid']}·acid": [R["acid"]], f"ARG{R['c1']}·clamp": [R["c1"]],
            f"ARG{R['c2']}·clamp": [R["c2"]], f"HIS{R['sh']}·stab": [R["sh"]],
            f"TRP{R['sw']}·stab": [R["sw"]], f"TYR{R['sy']}·stab": [R["sy"]]}
    res8_g = {k: [aid2gid[int(a)] for a in fs.atom
                  if a.resnum in set(v) and a.element.strip() != "H" and int(a) in aid2gid]
              for k, v in res8.items()}

    # ── per-frame geometry over the strided trajectory ────────────────────────────────
    # attack distance + SN2 angle, fluoride cradle, carboxylate clamp, per-residue
    # engagement, protein Rg and ligand-COM displacement - all minimum-image (PBC-aware).
    tr = traj.read_traj(str(trj))
    stride = max(1, int(getattr(CFG, "DEFLUOR_STRIDE", 100)))
    idx = list(range(0, len(tr), stride))
    nfr = len(idx)
    if nfr == 0:
        _echo("    ✘ Skipped      : trajectory recorded zero frames (MD-kill / empty _trj).")
        return None
    total_ns = md_ns if md_ns and md_ns > 0 else 1000.0
    t = np.linspace(0, total_ns, nfr)

    def _mind(pa, pb, box):
        d = _defl_min_image(pa[:, None, :] - pb[None, :, :], box).reshape(-1, 3)
        return np.linalg.norm(d, axis=1).min()

    attack = np.empty(nfr); angle = np.empty(nfr); rg = np.empty(nfr); com = np.empty(nfr)
    cradle_d = {k: np.empty(nfr) for k in cradle}
    clamp_d = {k: np.empty(nfr) for k in clamp}
    eng = np.empty((len(res8), nfr))
    # pocket COM reference: protein Cα atoms within DEFLUOR_POCKET_RADIUS_A of the ligand at frame 0
    p0 = tr[idx[0]].pos()
    lig0, ca0 = p0[gLig], p0[gCA]
    dmask = np.linalg.norm(ca0[:, None, :] - lig0[None, :, :], axis=2).min(axis=1) < _DEFL_POCKET
    poc_ref = ca0[dmask].mean(axis=0) if dmask.any() else ca0.mean(axis=0)

    for m, fi in enumerate(idx):
        fr = tr[fi]; pos = fr.pos(); box = fr.box
        od = pos[gOd]; c = pos[gC][None, :]
        dd = np.linalg.norm(_defl_min_image(od - c, box), axis=1)
        attack[m] = dd.min()
        o = od[np.argmin(dd)]
        v1 = _defl_min_image((o - pos[gC])[None, :], box)[0]
        n1 = np.linalg.norm(v1)
        # Leaving fluoride, chosen per frame: SN2 breaks the C-F bond anti-periplanar to the
        # incoming nucleophile, so among the alpha-carbon fluorines the reactive one is whichever
        # is most backside-aligned (largest Od-Ca-F angle, nearest the collinear 180 deg attack
        # trajectory). On a rotating CF2/CF3 this identity changes frame-to-frame; a fixed choice
        # would track a spectator F and read a spuriously low angle. The chosen F also sets the
        # fluoride-cradle distance for that frame.
        best = -1.0; f = pos[gFs[0]]
        for gFi in gFs:
            vf = _defl_min_image((pos[gFi] - pos[gC])[None, :], box)[0]
            a = np.degrees(np.arccos(np.clip(v1 @ vf / (n1 * np.linalg.norm(vf)), -1, 1)))
            if a > best:
                best = a; f = pos[gFi]
        angle[m] = best
        for k, g in cradle.items():
            cradle_d[k][m] = _mind(f[None, :], pos[g], box) if g else np.nan
        for k, g in clamp.items():
            clamp_d[k][m] = _mind(pos[g], pos[gCox], box) if g else np.nan
        for j, g in enumerate(res8_g.values()):
            eng[j, m] = _mind(pos[g], pos[gLig], box) if g else np.nan
        ca = pos[gCA]
        rg[m] = np.sqrt(((ca - ca.mean(0)) ** 2).sum(1).mean())
        com[m] = np.linalg.norm(pos[gLig].mean(0) - poc_ref)

    nac = (attack < _DEFL_NAC_DIST) & (angle > _DEFL_NAC_ANGLE)
    nac_pct = 100 * nac.mean()
    occ = {name: 100 * np.nanmean(eng[j] < _DEFL_ENGAGE) for j, name in enumerate(res8)}
    dist_pct = 100 * (attack < _DEFL_NAC_DIST).mean()
    angle_pct = 100 * (angle > _DEFL_NAC_ANGLE).mean()

    # ── MM-GBSA binding joined per frame ──────────────────────────────────────────────
    # Frame-stamped CSV; the defluorination (strided) frames are a subset of the scored
    # frames, so they join directly. Reactivity (NAC) and binding (ΔG_bind) together gate
    # turnover - a ligand must reach the attack pose AND stay bound. GB overstabilises
    # anionic PFAS, so ΔG_bind (total + per-component terms) is relative-only.
    _COMPS = [("Coulomb", "_Coulomb"), ("vdW", "_vdW"), ("Hbond", "_Hbond"), ("Lipo", "_Lipo"),
              ("Packing", "_Packing"), ("SelfCont", "_SelfCont"), ("Solv_GB", "_Solv_GB"),
              ("Solv_SA", "_Solv_SA"), ("Covalent", "_Covalent")]
    dG = np.full(nfr, np.nan)
    comp = {}                          # component label → per-frame ΔG contribution (joined)
    mmgbsa_median = np.nan; mmgbsa_cond = np.nan; _cond_label = "n/a"
    _mmcsv = _mmgbsa_csv(job_dir, job_name)
    if _mmcsv is not None:
        try:
            _mm = pd.read_csv(_mmcsv)
            if MMGBSA_FRAME_COL in _mm.columns and CFG.MMGBSA_DG_COLUMN in _mm.columns:
                _mm.index = pd.to_numeric(_mm[MMGBSA_FRAME_COL], errors="coerce")
                _sel = _mm.reindex(idx)          # rows for our (strided) frames; NaN where unscored
                dG = pd.to_numeric(_sel[CFG.MMGBSA_DG_COLUMN], errors="coerce").to_numpy(float)
                for _lab, _suf in _COMPS:
                    if CFG.MMGBSA_DG_COLUMN + _suf in _sel.columns:
                        _arr = pd.to_numeric(_sel[CFG.MMGBSA_DG_COLUMN + _suf],
                                             errors="coerce").to_numpy(float)
                        if np.isfinite(_arr).any():        # drop components Prime left empty
                            comp[_lab] = _arr
                _fin = dG[np.isfinite(dG)]
                if _fin.size:
                    mmgbsa_median = float(np.median(_fin))
                    if nac.any():
                        _c, _cond_label = dG[nac], "NAC frames"
                    else:
                        _c, _cond_label = dG[np.argsort(attack)[:max(1, nfr // 10)]], "closest-approach 10%"
                    _c = _c[np.isfinite(_c)]
                    mmgbsa_cond = float(np.median(_c)) if _c.size else np.nan
        except Exception:
            pass

    # ── merged per-frame CSV ──────────────────────────────────────────────────────────
    hdr = ["frame", "time_ns", "attack_Od_Ca", "sn2_angle", "nac_competent"]
    cols = [np.array(idx), t, attack, angle, nac.astype(int)]
    for k in cradle: hdr.append(f"Fcradle_{k}"); cols.append(cradle_d[k])
    for k in clamp: hdr.append(f"clamp_{k}"); cols.append(clamp_d[k])
    for j, name in enumerate(res8): hdr.append("engage_" + name.split("·")[0]); cols.append(eng[j])
    hdr += ["protein_Rg", "lig_com_disp", "mmgbsa_dG_bind"]; cols += [rg, com, dG]
    for _lab in comp:
        hdr.append("mmgbsa_" + _lab); cols.append(comp[_lab])
    np.savetxt(out / "04_Defluorination_Geometry.csv", np.column_stack(cols),
               delimiter=",", header=",".join(hdr), comments="")

    # ── figures (black axis labels/values, neutral grey grid, CFG colours/DPI) ──────────
    _dpi = int(CFG.VIS_FIGURE_DPI)
    _LF, _FA = CFG.VIS_FONT_LEGEND, CFG.VIS_LEGEND_FRAME_ALPHA

    def _dgrid(ax):
        ax.grid(axis="x", color=_DEFL_GRID_X, lw=float(CFG.DEFLUOR_GRID_X_LW), alpha=float(CFG.DEFLUOR_GRID_X_ALPHA))
        ax.grid(axis="y", color=_DEFL_GRID_Y, lw=float(CFG.DEFLUOR_GRID_Y_LW), alpha=float(CFG.DEFLUOR_GRID_Y_ALPHA))
        ax.set_axisbelow(True)

    # 01 reactive-pose trajectory - three stacked panels sharing one Time axis: (A) SN2 attack geometry
    # (Oδ···Cα distance on the left axis, backside attack angle on the right), (B) leaving-F fluoride
    # cradle, (C) carboxylate clamp. The near-attack window (global NAC frames: attack < cut AND angle >
    # cut) is shaded green in every panel so cradle/clamp engagement can be read at the frames when the
    # ligand is attack-ready. It is the same global reaction-coordinate mask on all three panels - never a
    # per-residue claim - and is named in each panel's legend so it cannot be misread as residue-specific.
    # The attack angle is a property of the ligand carbon (Oδ-Cα-F), so its right-hand axis appears on
    # panel A only, not mirrored onto the cradle/clamp panels where it would be meaningless. Residue names
    # come from the trace keys (mapped per homolog), never hardcoded; tags min-separate in y so they never
    # overlap. Distance panels use unit (1 Å) y-ticks over a 0-based range spanning the data and the cut-off.
    _axlab = float(CFG.VIS_FONT_AXIS_LABEL); _lblf = float(CFG.VIS_FONT_TICK_DENSE)
    _band_a = float(CFG.DEFLUOR_NAC_BAND_ALPHA)
    _xleft = -total_ns * float(CFG.DEFLUOR_LEFT_MARGIN)     # left xlim: empty margin the residue tags sit inside
    _xlab = -total_ns * float(CFG.DEFLUOR_EDGE_LABEL_X)     # tag x anchor just left of t=0 (ha="right")

    def _start_level(d):
        return float(np.nanmean(d[:max(1, len(d) // 50)]))   # avg of the first ~2 % of frames

    def _edge_labels(ax, entries):
        """Write each (label, colour, start_y) horizontally in the left gap, right-aligned just left of the
        first frame and min-separated in y (falling back down if the stack overruns the top) so tags that
        start at nearly the same distance never overlap."""
        _lo, _hi = ax.get_ylim(); _span = _hi - _lo
        _gap = _span * float(CFG.DEFLUOR_EDGE_LABEL_GAP)
        _ent = sorted(entries, key=lambda e: e[2]); _ys = [e[2] for e in _ent]
        for _i in range(1, len(_ys)):
            if _ys[_i] < _ys[_i - 1] + _gap:
                _ys[_i] = _ys[_i - 1] + _gap
        _shift = max(0.0, _ys[-1] - (_hi - _span * 0.03))
        for (_lab, _col, _), _y in zip(_ent, _ys):
            ax.text(_xlab, min(_y - _shift, _hi - _span * 0.03), _lab, ha="right", va="center",
                    color=_col, fontsize=_lblf, fontweight="bold", clip_on=False)

    fig, (axA, axB, axC) = plt.subplots(3, 1, figsize=tuple(CFG.DEFLUOR_MERGED_FIGSIZE), sharex=True,
                                        gridspec_kw={"height_ratios": list(CFG.DEFLUOR_PANEL_HEIGHT_RATIOS)})

    # panel A - SN2 attack geometry (distance left, backside angle right)
    axA.plot(t, attack, color=_DEFL_OKABE[0], lw=0.9)
    _cut_d = axA.axhline(_DEFL_NAC_DIST, ls="--", color=_DEFL_OKABE[0], alpha=0.55)
    axA.set_ylabel("Oδ···Cα attack distance (Å)", color=_DEFL_OKABE[0], labelpad=3, fontsize=_axlab)
    _top = float(np.ceil(max(float(np.nanmax(attack)), _DEFL_NAC_DIST)))
    axA.set_ylim(0, _top); axA.set_yticks(np.arange(0, _top + 0.001, 1))
    axA.tick_params(axis="y", labelcolor=_DEFL_OKABE[0]); _dgrid(axA)
    _edge_labels(axA, [(f"ASP{R['nuc']}·Nu", _DEFL_OKABE[0], _start_level(attack))])
    _aA = axA.twinx()
    _aA.plot(t, angle, color=_DEFL_OKABE[1], lw=0.7, alpha=0.85)
    _cut_a = _aA.axhline(_DEFL_NAC_ANGLE, ls="--", color=_DEFL_OKABE[1], alpha=0.6)
    _aA.set_ylabel("Oδ–Cα–F backside attack angle (°)", color=_DEFL_OKABE[1], labelpad=3, fontsize=_axlab)
    _aA.set_ylim(0, 180); _aA.set_yticks(np.arange(0, 181, 30)); _aA.tick_params(axis="y", labelcolor=_DEFL_OKABE[1])
    _band = axA.fill_between(t, 0, axA.get_ylim()[1], where=nac, color=_DEFL_GREEN, alpha=_band_a, step="mid")
    _leg = _aA.legend([_cut_d, _cut_a, _band],
                      [f"NAC distance ≤ {_DEFL_NAC_DIST:.1f} Å", f"SN2 in-line attack ≥ {_DEFL_NAC_ANGLE:.0f}°",
                       f"NAC-competent ({nac_pct:.1f}%)"],
                      loc="lower right", ncol=3, framealpha=_FA, fontsize=_LF,
                      handlelength=1.8, columnspacing=1.2, borderpad=0.4)
    _leg.set_zorder(20)

    # panel B - fluoride cradle (leaving F to each stabiliser)
    _cmax = float(np.nanmax([np.nanmax(d) for d in cradle_d.values()])) * 1.05
    _lab_top = []
    for (k, d), c in zip(cradle_d.items(), _DEFL_OKABE):
        axB.plot(t, d, lw=0.8, color=c); _lab_top.append((f"F···{k}", c, _start_level(d)))
    _dlB = axB.axhline(_DEFL_NAC_DIST, ls="--", color=_DEFL_GREY, alpha=0.7)
    axB.set_ylabel("leaving F ··· donor distance (Å)", labelpad=2, color=_DEFL_AXTXT)
    axB.set_ylim(0, _cmax); axB.set_yticks(np.arange(0, _cmax + 0.001, 1))
    axB.tick_params(axis="y", labelcolor=_DEFL_AXTXT); _dgrid(axB)
    _bandB = axB.fill_between(t, 0, _cmax, where=nac, color=_DEFL_GREEN, alpha=_band_a, step="mid")
    _edge_labels(axB, _lab_top)
    _lB = axB.legend([_dlB, _bandB], [f"F cradled ≤ {_DEFL_NAC_DIST:.1f} Å", "near-attack window"],
                     loc="lower right", ncol=2, framealpha=_FA, fontsize=_LF,
                     handlelength=1.8, columnspacing=1.2, borderpad=0.4)
    _lB.set_zorder(20)

    # panel C - carboxylate clamp (each Arg NHx to the ligand carboxylate)
    _lmax = float(np.nanmax([np.nanmax(d) for d in clamp_d.values()])) * 1.05
    _lab_bot = []
    for (k, d), c in zip(clamp_d.items(), [_DEFL_OKABE[4], _DEFL_OKABE[2]]):
        axC.plot(t, d, lw=0.8, color=c); _lab_bot.append((k, c, _start_level(d)))
    _dlC = axC.axhline(_DEFL_NAC_DIST, ls="--", color=_DEFL_GREY, alpha=0.7)
    axC.set_ylabel("Arg NHx ··· carboxylate O (Å)", labelpad=2, color=_DEFL_AXTXT)
    axC.set_ylim(0, _lmax); axC.set_yticks(np.arange(0, _lmax + 0.001, 1))
    axC.tick_params(axis="y", labelcolor=_DEFL_AXTXT); _dgrid(axC)
    _bandC = axC.fill_between(t, 0, _lmax, where=nac, color=_DEFL_GREEN, alpha=_band_a, step="mid")
    _edge_labels(axC, _lab_bot)
    _lC = axC.legend([_dlC, _bandC], [f"clamp engaged ≤ {_DEFL_NAC_DIST:.1f} Å", "near-attack window"],
                     loc="upper right", ncol=2, framealpha=_FA, fontsize=_LF,
                     handlelength=1.8, columnspacing=1.2, borderpad=0.4)
    _lC.set_zorder(20)

    # one shared Time axis: label + tick VALUES on the bottom panel only
    axC.set_xlabel("Time (ns)", labelpad=2, color=_DEFL_AXTXT)
    axC.set_xticks(np.arange(0, total_ns + 1, 50)); axC.tick_params(axis="x", labelcolor=_DEFL_AXTXT)
    axA.set_xlim(_xleft, total_ns * 1.005)
    fig.tight_layout(h_pad=0.6); fig.savefig(out / "01_Reactive_Pose_Trajectory.svg", dpi=_dpi); plt.close(fig)

    # 02 reactive summary - each bar's y-tick label AND its value (written just right of the bar, never
    # inside) take that bar's colour, so residue/criterion, bar and number all read as one coloured unit.
    fig, (axa, axb2) = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1.7, 1]})
    names = list(occ.keys()); vals = [occ[k] for k in names]
    _cols_a = list(_DEFL_OKABE[:len(names)])
    _txt_a = [_readable(c) for c in _cols_a]   # bars keep the bright hue; label/value text is legible
    axa.barh(names, vals, color=_cols_a); axa.invert_yaxis()
    axa.set_xlabel(f"% of trajectory within {_DEFL_ENGAGE:.0f} Å of ligand")
    axa.grid(axis="x", alpha=0.3); axa.set_xlim(0, max(vals) * 1.28 + 1)
    for _lbl, _col in zip(axa.get_yticklabels(), _txt_a):
        _lbl.set_color(_col)
    for i, (v, _col) in enumerate(zip(vals, _txt_a)):
        axa.text(v + 0.3, i, f"{v:.1f}%", va="center", fontsize=CFG.VIS_FONT_LEGEND, color=_col)
    crit = [(f"attack distance < {_DEFL_NAC_DIST:.1f} Å", dist_pct, _DEFL_OKABE[0]),
            (f"SN2 angle > {_DEFL_NAC_ANGLE:.0f}°", angle_pct, _DEFL_OKABE[1]),
            ("NAC-competent (both)", nac_pct, _DEFL_GREEN)]
    _cols_b = [c[2] for c in crit]
    _txt_b = [_readable(c) for c in _cols_b]
    axb2.barh([c[0] for c in crit], [c[1] for c in crit], color=_cols_b); axb2.invert_yaxis()
    axb2.set_xlabel("% of frames satisfying criterion")
    axb2.grid(axis="x", alpha=0.3); axb2.set_xlim(0, max(max(c[1] for c in crit) * 1.35, 1) + 1)
    for _lbl, _col in zip(axb2.get_yticklabels(), _txt_b):
        _lbl.set_color(_col)
    for i, (c, _col) in enumerate(zip(crit, _txt_b)):
        axb2.text(c[1] + 0.3, i, f"{c[1]:.1f}%", va="center", fontsize=CFG.VIS_FONT_TICK, color=_col)
    fig.tight_layout(); fig.savefig(out / "02_Reactive_Summary.svg", dpi=_dpi); plt.close(fig)

    # 03 binding vs reactivity (only when MM-GBSA is available for this rank)
    if np.isfinite(dG).any():
        m = np.isfinite(dG)
        _dv = dG[m]
        # clip the y-axis to the physical spread: a few Prime blown-up minimisations reach
        # hundreds of kcal/mol. A Tukey fence (Q1 − 3·IQR) drops those so the real 0…−20 band
        # fills the panel instead of leaving it mostly empty.
        _q1, _q3 = np.percentile(_dv, [25, 75])
        _phys = _dv[_dv >= _q1 - 3.0 * (_q3 - _q1)]
        _lo = float(_phys.min()) if _phys.size else float(_dv.min())
        fig, ax = plt.subplots(figsize=(7.8, 5.2))
        ax.scatter(attack[m & ~nac], dG[m & ~nac], s=14, color=_DEFL_GREY, alpha=0.7, label="non-NAC frame")
        if (m & nac).any():
            ax.scatter(attack[m & nac], dG[m & nac], s=28, color=_DEFL_GREEN, label="NAC-competent frame")
        ax.axvline(_DEFL_NAC_DIST, ls="--", color=_DEFL_OKABE[0], alpha=0.6, label=f"attack cutoff {_DEFL_NAC_DIST} Å")
        ax.set_ylim(_lo * 1.1, 5.0)
        _xa = attack[m]
        # 0.5-Å ticks, horizontal: a step of 2 left only "2" and "4" on the axis, and the 90° rotation
        # was needless for short numbers.
        _x0 = np.floor(float(_xa.min()) * 2) / 2
        _x1 = np.ceil(float(_xa.max()) * 2) / 2
        ax.set_xticks(np.arange(_x0, _x1 + 0.01, 0.5))
        ax.set_xlabel("attack distance (Å)")
        ax.set_ylabel("MM-GBSA ΔG$_{bind}$ (kcal/mol)  ·  relative-only")
        ax.grid(alpha=0.3); ax.legend(loc="upper right", ncol=3, fontsize=_LF, framealpha=_FA)
        fig.tight_layout(); fig.savefig(out / "03_Binding_vs_Reactivity.svg", dpi=_dpi); plt.close(fig)

    # The MM-GBSA energy-component decomposition (reactive vs rest) is drawn once, in Step 07's
    # per-rank MMGBSA_NAC_Decomposition figure, from the same per-frame Prime terms; the per-component
    # columns are still written to the geometry CSV below so 07 (or any reader) can decompose them.

    # ── run log + figure-description file ─────────────────────────────────────────────
    log = [f"FAcD defluorination MD analysis - Rank {rank}",
           f"ligand atoms: alpha-C={fs.atom[aC_aid].pdbname.strip()} "
           f"alpha-C F={[fs.atom[i].pdbname.strip() for i in F_aids]} (leaving F chosen per frame) "
           f"carboxylate-O={[fs.atom[i].pdbname.strip() for i in cox_aids]}",
           f"frames {nfr} (stride {stride}) · {total_ns:.0f} ns",
           f"NAC-competent (dist<{_DEFL_NAC_DIST} Å & angle>{_DEFL_NAC_ANGLE:.0f}°): {nac_pct:.2f}%",
           f"attack distance: mean {attack.mean():.2f} Å  min {attack.min():.2f} Å",
           f"SN2 attack angle: mean {angle.mean():.1f}°",
           f"ligand COM displacement: start {com[:5].mean():.1f} Å  end {com[-20:].mean():.1f} Å",
           f"protein Rg: mean {rg.mean():.2f} Å"]
    if np.isfinite(mmgbsa_median):
        log.append(f"MM-GBSA ΔG_bind: median {mmgbsa_median:.1f} kcal/mol "
                   f"(relative-only; GB overstabilises anionic PFAS)")
        log.append(f"conditioned ΔG_bind ({_cond_label}): {mmgbsa_cond:.1f} kcal/mol")
    else:
        log.append("MM-GBSA ΔG_bind: not available for this rank")
    log.append("per-residue engagement occupancy:")
    log += [f"  {name:16s} {v:5.1f}%" for name, v in occ.items()]
    (out / "05_Analysis_Log.log").write_text("\n".join(log) + "\n")
    (out / "06_Figure_Descriptions.txt").write_text(
        f"FAcD Defluorination MD Analysis - Rank {rank} · {nfr} frames · {total_ns:.0f} ns\n"
        f"01 reactive-pose trajectory - 3 panels sharing one Time axis: (A) SN2 attack geometry\n"
        f"   (Oδ···Cα distance + Oδ-Cα-F backside angle), (B) fluoride cradle (F···stabilisers),\n"
        f"   (C) carboxylate clamp (Arg···carboxylate); green = near-attack window (global NAC frames)\n"
        f"02 reactive summary (per-residue engagement + NAC criterion decomposition)\n"
        f"03 binding vs reactivity (MM-GBSA ΔG_bind vs attack distance; only if MM-GBSA present)\n"
        f"04 per-frame geometry CSV (geometry + mmgbsa_dG_bind + per-component terms; the MM-GBSA\n"
        f"   component decomposition itself is Step 07's MMGBSA_NAC_Decomposition figure)\n"
        f"05 run log · 06 this figure-description file\n")
    _dg_note = f" · ΔG {mmgbsa_median:.0f}" if np.isfinite(mmgbsa_median) else ""
    _echo(f"    ✔ Defluor      : NAC {nac_pct:.1f}% · attack min {attack.min():.2f} Å{_dg_note} → {out.name}/")
    return out


def plot_defluor_combined(md_dir: Path, ligands: "dict | None" = None,
                          controls: "set | None" = None) -> None:
    """Cross-rank defluorination comparison (drawn at Finalise), each bar tagged with its ligand
    short name and the CFG control label (e.g. 3R3U-FA) for easy comparative reading."""
    ligands, controls = ligands or {}, controls or set()
    analysis = _analysis_dir(md_dir.parent)                                       # 06_Analysis
    root = analysis / getattr(CFG, "DEFLUOR_OUTPUT_SUBDIR", "Defluorination")     # per-rank subfolders
    if not root.is_dir():
        return
    rows = []
    # Rank from the 'Defluorination_R{N}' folder name (note: _natural_rank expects the '_R_N' spelling
    # and does NOT match this one). Non-controls by rank, the 3R3U control LAST - consistent with every
    # other combined figure.
    def _defl_rank(p: Path) -> int:
        m = re.search(r"_R(\d+)", p.name)
        return int(m.group(1)) if m else 0
    for d in sorted(root.glob("Defluorination_R*"),
                    key=lambda p: (_defl_rank(p) in controls, _defl_rank(p))):
        csv = d / "04_Defluorination_Geometry.csv"
        if not csv.is_file():
            continue
        try:
            df = pd.read_csv(csv)
        except Exception:
            continue
        rk = d.name.replace("Defluorination_", "")
        _dg = np.nan
        if "mmgbsa_dG_bind" in df.columns:
            _v = pd.to_numeric(df["mmgbsa_dG_bind"], errors="coerce").dropna()
            _dg = float(_v.median()) if len(_v) else np.nan
        rows.append((rk, 100 * df["nac_competent"].mean(),
                     df["attack_Od_Ca"].min(), df["attack_Od_Ca"].mean(), df["sn2_angle"].mean(), _dg))
    if not rows:
        return
    rks = [r[0] for r in rows]

    def _rank_lig_label(rk):
        _d = re.sub(r"\D", "", str(rk))
        _lig = _ctrl_label(int(_d) if _d else rk, ligands, controls)
        return f"{_lig}\n{rk}"
    _labs = [_rank_lig_label(rk) for rk in rks]

    def _labels(a, vals, fmt):
        top = max([v for v in vals if np.isfinite(v)] or [1])
        for i, v in enumerate(vals):
            if np.isfinite(v):
                a.text(i, v + top * 0.02, fmt.format(v), ha="center", va="bottom", fontsize=CFG.VIS_FONT_LEGEND)

    _dgs = [r[5] for r in rows]
    _have_dg = any(np.isfinite(v) for v in _dgs)
    fig, _axg = plt.subplots(2, 2, figsize=(11, 8.6))   # 2×2, balanced
    ax = _axg.flatten()
    nac_vals = [r[1] for r in rows]
    # Per-rank colours (one colour per rank, consistent across all four panels; the control distinct) -
    # mirrors 01_Physics_Build_Solvation_QC so a rank reads the same colour everywhere.
    _palette = [CFG.VIS_ACCENT["green"], CFG.VIS_ACCENT["amber"], CFG.VIS_ACCENT["vermillion"],
                CFG.VIS_ACCENT["blue"], CFG.VIS_ACCENT["magenta"], CFG.VIS_ACCENT["sky"]]
    _rank_ints = [int(re.sub(r"\D", "", str(rk)) or 0) for rk in rks]
    cols = _ctrl_palette(_rank_ints, controls, _palette)
    ax[0].bar(rks, nac_vals, color=cols)
    ax[0].set_ylabel("NAC-competent (%)"); ax[0].set_ylim(0, max(max(nac_vals), 1.0) * 1.2)
    _labels(ax[0], nac_vals, "{:.2f}%")
    ax[1].bar(rks, [r[2] for r in rows], color=cols); ax[1].set_ylabel("min attack distance (Å)")
    ax[1].axhline(_DEFL_NAC_DIST, ls="--", color=_DEFL_GREY); _labels(ax[1], [r[2] for r in rows], "{:.2f}")
    ax[2].bar(rks, [r[4] for r in rows], color=cols); ax[2].set_ylabel("mean SN2 angle (°)")
    ax[2].axhline(_DEFL_NAC_ANGLE, ls="--", color=_DEFL_GREY); _labels(ax[2], [r[4] for r in rows], "{:.1f}°")
    _used = [ax[0], ax[1], ax[2]]
    if _have_dg:
        _fin = [v if np.isfinite(v) else 0.0 for v in _dgs]
        ax[3].bar(rks, _fin, color=cols)
        ax[3].set_ylabel("median MM-GBSA ΔG$_{bind}$ (kcal/mol)  ·  relative")
        ax[3].axhline(0, color=_DEFL_GREY, lw=0.8)
        for i, v in enumerate(_dgs):
            ax[3].text(i, _fin[i], "n/a" if not np.isfinite(v) else f"{v:.1f}", ha="center",
                       va="bottom" if _fin[i] >= 0 else "top", fontsize=CFG.VIS_FONT_LEGEND)
        _used.append(ax[3])
    else:
        ax[3].axis("off")
    # X-axis labels only on the BOTTOM row (top panels share the same categories); tick labels coloured
    # per rank to match the bars - the 01_Physics_Build_Solvation_QC convention.
    _bottom_axes = {ax[2], ax[3] if _have_dg else ax[1]}
    for a in _used:
        a.grid(axis="y", alpha=0.3)
        a.set_xticks(range(len(rks)))
        if a in _bottom_axes:
            a.set_xticklabels(_labs)
            for _t, _c in zip(a.get_xticklabels(), cols):
                _t.set_color(_c); _t.set_fontweight("bold")
            # rank tick labels already identify the complexes; no redundant axis title
        else:
            a.tick_params(labelbottom=False)
    fig.tight_layout()
    fig.savefig(analysis / "05_Defluorination_Combined_AllRanks.svg",
                dpi=int(CFG.VIS_FIGURE_DPI))
    plt.close(fig)
    _echo(f"  ✔ Defluorination combined figure → {analysis.name}/{_utils_mod.deflx_fig_name('05_Defluorination_Combined_AllRanks.svg')}")


# =============================================================================
# SECTION 9: ESP PHYSICS - MSJ BUILDERS & DETACHED JOB ENVIRONMENT
# =============================================================================

STAGES_ALL = ("merge", "watermap", "build", "md")   # execution order: WaterMap needs only the ESP complex

def _build_msj(lig_indices: str) -> str:
    """Desmond System Builder .msj - every physical setting from CFG (§17b PHYS_*).

    Explicit-solvent box (PHYS_SOLVENT_MODEL, PHYS_BOX_SHAPE, PHYS_BOX_BUFFER_A), auto-neutralise
    with PHYS_COUNTERION, PHYS_SALT_CONC_M background salt, and PHYS_FORCEFIELD. Ions/salt are kept
    ≥ PHYS_ION_EXCLUDE_A from the ligand so a counter-ion never lands on the Asp-Oδ and corrupts the
    NAC; ion_awayfrom takes the 1-based LIG atom indices.
    """
    return (
        f'task {{ task = "desmond:auto" }}\n'
        f'build_geometry {{\n'
        f'  add_counterion = {{ ion = "{CFG.PHYS_COUNTERION}" number = "neutralize_system" }}\n'
        f'  box = {{ shape = "{CFG.PHYS_BOX_SHAPE}" size = [{CFG.PHYS_BOX_BUFFER_A} {CFG.PHYS_BOX_BUFFER_A} {CFG.PHYS_BOX_BUFFER_A} ] size_type = "buffer" }}\n'
        f'  override_forcefield = "{CFG.PHYS_FORCEFIELD}"\n'
        f'  rezero_system = false\n'
        f'  salt = {{ concentration = {CFG.PHYS_SALT_CONC_M} negative_ion = "{CFG.PHYS_SALT_NEG_ION}" positive_ion = "{CFG.PHYS_SALT_POS_ION}" }}\n'
        f'  solvent = "{CFG.PHYS_SOLVENT_MODEL}"\n'
        f'  ion_awaydistance = {CFG.PHYS_ION_EXCLUDE_A}\n'
        f'  ion_awayfrom = [ {lig_indices} ]\n'
        f'}}\n'
        # The `water` key is REQUIRED: assign_forcefield defaults water to SPC and re-stamps every
        # solvent residue, so build_geometry's solvent="TIP3P" (box geometry only) is silently overridden.
        f'assign_forcefield {{ forcefield = "{CFG.PHYS_FORCEFIELD}" water = "{CFG.PHYS_SOLVENT_MODEL}" }}\n'
    )


def _restraint_block() -> str:
    """The ligand + backbone positional restraint (Desmond `restrain`), or '' when MD_RESTRAIN_LIGAND
    is off. One definition, used by BOTH the final relaxation stage and production so the near-attack
    docked pose is held continuously from equilibration into production.

    Positionally restrain the ligand heavy atoms so a small substrate cannot diffuse out of the pocket
    (the escape artefact) or relax out of the reactive geometry, plus a gentle backbone anchor so the
    protein does not translate/tumble in the box and drag the ligand restraint off the moving site.
    Proven msj syntax (cf. $SCHRODINGER .../data/desmond/kinetics_membrane_md.msj production stage).
    """
    if not getattr(CFG, "MD_RESTRAIN_LIGAND", False):
        return ""
    _lig_res = getattr(CFG, "LIGAND_RESNAME_ASSERT", "LIG")
    return (
        f'  restrain = [\n'
        f'    {{atom = "asl: (res.ptype {_lig_res}) and not (atom.ele H)" force_constant = {CFG.MD_RESTRAIN_LIG_FORCE_K}}}\n'
        f'    {{atom = "asl: (backbone) and not (atom.ele H)" force_constant = {CFG.MD_RESTRAIN_BB_FORCE_K}}}\n'
        f'  ]\n'
    )


def _md_production(time_ps: float, interval_ps: float) -> str:
    """MD production stage appended to the Desmond relaxation protocol (minimise + staged NVT/NPT
    equilibration with restraints, then production). NPT at MD_EQUIL_TARGET_T / PHYS_MD_PRESSURE_BAR,
    RESPA PHYS_MD_TIMESTEP_PS, energies every PHYS_MD_ENESEQ_PS ps - all from CFG §17b.

    The ligand + backbone restraint (`_restraint_block`) is ON for production; the same block is also
    injected into the final relaxation stage by `_md_msj`, so the docked catalytic pose is held from
    equilibration onward rather than being allowed to collapse just before production begins.
    """
    temp = CFG.MD_EQUIL_TARGET_T                                   # one temperature for MD + MM-GBSA
    dt = " ".join(str(x) for x in CFG.PHYS_MD_TIMESTEP_PS)
    _restrain = _restraint_block()
    return (
        f'\nsimulate {{\n'
        f'  title       = "Production MD"\n'
        f'  time        = {time_ps}\n'
        f'  timestep    = [{dt} ]\n'
        f'  temperature = {temp}\n'
        f'  pressure    = [{CFG.PHYS_MD_PRESSURE_BAR} isotropic ]\n'
        f'  ensemble = {{ class = NPT method = MTK thermostat.tau = {CFG.PHYS_MD_THERMOSTAT_TAU} barostat.tau = {CFG.PHYS_MD_BAROSTAT_TAU} }}\n'
        f'  randomize_velocity = {{ first = 0.0 interval = inf temperature = {temp} seed = {CFG.PHYS_MD_SEED} }}\n'
        f'  eneseq.interval    = {CFG.PHYS_MD_ENESEQ_PS}\n'
        f'  trajectory = {{ interval = {interval_ps} center = solute write_velocity = false }}\n'
        f'{_restrain}'
        f'}}\n'
    )


def _clean_job_env() -> dict:
    """A fresh login-shell environment for launching WaterMap, detached from the $SCHRODINGER/run session.

    The script runs under $SCHRODINGER/run, which sets PYTHONHOME/PYTHONPATH/LD_LIBRARY_PATH/
    SCHRODINGER_EXEC. A WaterMap job that inherits those cannot stage its GCMC ligand companion to the
    sub-job and dies at stage 8. A plain login-shell environment - critically with $SCHRODINGER on PATH
    so the sub-stages find the multisim/watermap utilities - lets the job run as it would from a
    terminal. Paired with the `bash -lc` launch in run_watermap.
    """
    keep = ("HOME", "USER", "LOGNAME", "DISPLAY", "LANG", "LC_ALL", "TERM", "SSH_AUTH_SOCK",
            "SCHRODINGER_LICENSE_FILE", "SCHRODINGER_CUSTOM_MONOMER_DB_PATH")
    env = {k: os.environ[k] for k in keep if k in os.environ}
    env["SCHRODINGER"] = SCHRO
    env["PATH"] = f"{SCHRO}:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
    env["TMPDIR"] = os.environ.get("TMPDIR", "/tmp")
    env.setdefault("USER", os.environ.get("USER", os.environ.get("LOGNAME", "user")))
    env.setdefault("LOGNAME", env["USER"])
    env.setdefault("LANG", "C.UTF-8")
    return env


def _cancel_launched_jobs() -> None:
    """Cancel every still-registered job on the job server. Idempotent (guarded), safe to call twice.

    `jsc list` defaults to active jobs (RUNNING/WAITING/PAUSED); the JobId is column 0 and the job name
    column 1. Only names THIS run registered are matched, so a concurrent user's jobs are never touched.
    """
    global _CLEANUP_DONE
    if _CLEANUP_DONE or not _LAUNCHED_JOBS:
        return
    _CLEANUP_DONE = True
    names = set(_LAUNCHED_JOBS)
    try:
        out = subprocess.run([f"{SCHRO}/jsc", "list", "-j"], capture_output=True, text=True, timeout=30).stdout
        # Match a registered job AND its stage subjobs: multisim runs the MD as stage subjobs named
        # "<jobname>_<stage>" (production is the last), so cancelling only the exact parent name misses
        # the detached stage that is actually on the GPU. Prefix-match catches "desmond_md_job_R_8_7".
        ids = [p[0] for ln in out.splitlines()
               if len(p := ln.split()) >= 2 and (p[1] in names or any(p[1].startswith(n + "_") for n in names))]
        if ids:
            _echo(f"\n  {_C.WARNING}[cleanup] script exiting - cancelling {len(ids)} running job(s): "
                  f"{', '.join(sorted(names))}{_C.ENDC}")
            subprocess.run([f"{SCHRO}/jsc", "cancel", *ids], timeout=90)
    except Exception:
        pass   # jsc unreachable → still fall through to the direct Desmond process-tree kill below
    _kill_desmond_jobs(names)


def _kill_desmond_jobs(names: "set[str]") -> None:
    """Force-terminate any Desmond MD process chain still on the GPU for a job THIS run launched.

    multisim hands the MD to jobserverd, which runs it under its own `job_supervisord` re-parented to
    init - so killing the 06 process (and its multisim parent) orphans `gdesmond` on the GPU, and a
    `jsc cancel` on the parent name does not always reach the detached stage subjob. Belt-and-braces:
    match running processes by one of THIS run's registered jobnames in their command line and kill the
    process group by explicit PID. Only names the run registered are matched (never a bare pattern that
    could hit a concurrent run's identically-named job), and only the Desmond executables - never an
    arbitrary process. NOTE: this runs from the SIGINT/SIGTERM/SIGHUP handler and atexit; a `kill -9`
    (SIGKILL) on 06 is uncatchable and will still orphan the job - use Ctrl-C or `kill` for a clean stop."""
    md_names = [n for n in names if n.startswith("desmond_md_job_")]
    if not md_names:
        return
    try:
        out = subprocess.run(["ps", "-eo", "pid,args"], capture_output=True, text=True, timeout=15).stdout
    except Exception:
        return
    _exes = ("gdesmond", "desmond_driver", "job_supervisord", "chorus_multijob")
    victims = []
    for ln in out.splitlines():
        parts = ln.strip().split(None, 1)
        if len(parts) == 2 and parts[0].isdigit() and any(e in parts[1] for e in _exes) \
                and any(n in parts[1] for n in md_names):
            victims.append(int(parts[0]))
    if not victims:
        return
    for _sig in (signal.SIGTERM, signal.SIGKILL):
        for pid in victims:
            try:
                os.killpg(os.getpgid(pid), _sig)
            except (ProcessLookupError, PermissionError):
                try:
                    os.kill(pid, _sig)
                except (ProcessLookupError, PermissionError):
                    pass
        if _sig is signal.SIGTERM:
            time.sleep(2)


def _run_tracked(cmd: list, *, cwd=None, stdout=None, stderr=None, check: bool = False):
    """subprocess.run for a LOCAL Schrödinger step (SID; the MM-GBSA thermal_mmgbsa driver), but launched
    in its own session and registered in _LAUNCHED_PROCS so a signal/atexit handler can terminate its whole
    process group. Same job-server jobs die via _cancel_launched_jobs; these local ones die here - together
    they are the MM-GBSA/SID equivalent of the WaterMap/MD kill-cleanup. Blocks like subprocess.run; with
    check=True a non-zero return raises CalledProcessError. start_new_session so one os.killpg reaches the
    tool AND every worker it forked, and so the terminal's own Ctrl-C does not race the handler."""
    proc = subprocess.Popen(cmd, cwd=cwd, stdout=stdout, stderr=stderr, start_new_session=True)
    _LAUNCHED_PROCS.add(proc)
    try:
        rc = proc.wait()
    finally:
        _LAUNCHED_PROCS.discard(proc)
    if check and rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)
    return rc


def _kill_launched_procs() -> None:
    """Terminate every still-registered LOCAL subprocess (SID, MM-GBSA drivers) and its process group, so a
    killed script leaves none running detached under init. Job-server Prime subjobs are cancelled separately
    by _cancel_launched_jobs. SIGTERM the groups, a short grace, then SIGKILL survivors."""
    procs = [p for p in list(_LAUNCHED_PROCS) if p.poll() is None]
    for _sig in (signal.SIGTERM, signal.SIGKILL):
        for p in procs:
            if p.poll() is not None:
                continue
            try:
                os.killpg(os.getpgid(p.pid), _sig)
            except Exception:
                try:
                    p.send_signal(_sig)
                except Exception:
                    pass
        if _sig is signal.SIGTERM:
            time.sleep(2)


_OOMD_GUARD: "OomdGuard | None" = None   # set in main once primed; restored on signal/atexit


def _install_job_cleanup() -> None:
    """Cancel the run's job-server jobs AND kill its local subprocesses on normal exit and on
    SIGINT/SIGTERM/SIGHUP, so a killed script never leaves WaterMap/MD/MM-GBSA jobs under jobserverd or
    SID/thermal_mmgbsa drivers detached under init. The same handlers restore systemd-oomd: the signal path
    re-raises with SIG_DFL and never unwinds the guard's `with`/atexit, so without this a Ctrl-C would leave
    oomd masked and stopped on the host permanently."""
    atexit.register(_kill_launched_procs)
    atexit.register(_cancel_launched_jobs)
    atexit.register(lambda: _OOMD_GUARD and _OOMD_GUARD.restore())

    def _handler(signum, _frame):
        _cancel_launched_jobs()                 # cancel job-server Prime/MD/WaterMap jobs
        _kill_launched_procs()                  # kill local SID / thermal_mmgbsa driver process groups
        if _OOMD_GUARD is not None:
            _OOMD_GUARD.restore()               # unmask systemd-oomd before the process dies on the signal
        signal.signal(signum, signal.SIG_DFL)   # restore default and re-raise for the correct exit status
        os.kill(os.getpid(), signum)

    for _s in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
        try:
            signal.signal(_s, _handler)
        except Exception:
            pass


# -----------------------------------------------------------------------------
# SECTION 10: INPUT DISCOVERY (handover complexes + ESP charges)
# -----------------------------------------------------------------------------
def _resnum(v) -> int:
    m = re.search(r"(\d+)", str(v or ""))
    return int(m.group(1)) if m else 0


def newest_ranked_csv(run: Path) -> Path:
    prod = run / "1_Boltz2_Production"
    cands = (sorted(prod.glob(CFG.GLOB_RANKED_CSV), key=lambda p: p.stat().st_mtime) or
             sorted(prod.glob("*Ranked*.csv"), key=lambda p: p.stat().st_mtime))
    if not cands:
        sys.exit(f"No ranked CSV ({CFG.GLOB_RANKED_CSV}) under {prod}")
    return cands[-1]


def ranked_rows_by_jobname(ranked_csv: Path) -> dict:
    """Map job_name → ranked-CSV row, so each handover complex can look up its Mapped_Base etc."""
    with open(ranked_csv, newline="") as fh:
        return {r.get("job_name", "").strip(): r for r in csv.DictReader(fh)}


def discover_handover(run: Path) -> list:
    """Discover the MD-selected complexes from 05's handover folder - the single source of truth.

    05 writes `06_<tier>_<count>hits_Molecular_Handover_Files/R{N}_<stem>.pdb`, where N is the complex's
    Scientific_Rank (so the third selected hit can be R8, not R3). This keys the whole run on that same
    N and stem - never on MD_Rank, which is only a 1..k position over the selected set. Returns a list of
    {rank, stem, prepared} sorted by rank, so the output folders (R_{rank}) match the handover exactly.
    """
    base = run / "5_TopN_and_Preparation" / "3_Comparative_Analysis"
    # Newest by mtime, NOT lexical: 05 names the folder 06_<tier>_<count>hits_… and never removes stale
    # ones, so a lexical [-1] sorts 06_…_3hits_… AFTER 06_…_12hits_… and would run the older, smaller cohort.
    hos = sorted(base.glob("06_*_Molecular_Handover_Files"), key=lambda p: p.stat().st_mtime)
    if not hos:
        sys.exit(f"No handover folder (06_*_Molecular_Handover_Files) under {base}")
    ho = hos[-1]
    out = []
    for pdb in sorted(ho.glob("R*_*.pdb")):
        m = re.match(r"R(\d+)_(.+)\.pdb$", pdb.name)
        if m:
            out.append({"rank": int(m.group(1)), "stem": m.group(2), "prepared": pdb})
    if not out:
        sys.exit(f"No R<N>_<stem>.pdb handover complexes in {ho.name}")
    ranks = [e["rank"] for e in out]
    if len(ranks) != len(set(ranks)):
        sys.exit(f"Duplicate handover ranks in {ho.name}: {ranks} - each R{{N}} must be unique")
    return sorted(out, key=lambda e: e["rank"])


def find_esp(esp_dir: Path, stem: str) -> Path:
    """The ESP .mae for one complex, matched exactly by stem (05 writes <stem>_ESP.mae)."""
    hits = sorted(esp_dir.glob(f"{stem}_ESP.mae"))
    if len(hits) != 1:
        raise RuntimeError(f"{len(hits)} matches for {stem}_ESP.mae in {esp_dir.name} (expected 1)")
    return hits[0]


# -----------------------------------------------------------------------------
# SECTION 11: STAGE 1 - ESP MERGE
# -----------------------------------------------------------------------------
def merge_esp(prepared_pdb: Path, esp_mae: Path, out_mae: Path, base_resnum: int) -> dict:
    """Write the ESP ligand charges onto the prepared holo complex, matched strictly by atom name.

    Every LIG atom name must be present in the ESP map; a missing name raises rather than silently
    taking a positional charge (a wrong-atom charge is invisible downstream). The catalytic-base
    protonation state at Mapped_Base is reported and warned on when it is not HID.
    """
    from schrodinger import structure
    cx = structure.StructureReader.read(str(prepared_pdb))
    _esp_st = structure.StructureReader.read(str(esp_mae))
    esp_by_name = {a.pdbname.strip(): a.partial_charge for a in _esp_st.atom}
    esp_formal = round(sum(a.formal_charge for a in _esp_st.atom))     # the ligand's true net charge (e.g. −1 for a carboxylate)
    lig = [a for a in cx.atom if a.pdbres.strip() == "LIG"]
    if not lig:
        raise RuntimeError(f"no LIG atoms in {prepared_pdb.name}")
    missing = sorted({a.pdbname.strip() for a in lig} - set(esp_by_name))
    if missing:
        raise RuntimeError(f"LIG atoms absent from ESP map (name mismatch): {', '.join(missing)}")

    tot = 0.0
    for a in lig:
        a.partial_charge = esp_by_name[a.pdbname.strip()]
        tot += a.partial_charge

    # A neutral/corrupt ESP file (charges never written, or run on the wrong protonation state) sums to
    # a net charge that does not match the ligand's formal charge. Building MD on a mis-charged ligand
    # silently changes the α-carbon electrophilicity the SN2 depends on, so refuse rather than proceed.
    if round(tot) != esp_formal:
        raise RuntimeError(
            f"ESP net charge {tot:+.3f} e (rounds to {round(tot):+d}) ≠ ligand formal charge "
            f"{esp_formal:+d} in {esp_mae.name} - the ESP charges look wrong/neutral; refusing to "
            f"build MD on a mis-charged ligand")

    base_state = "?"
    for res in cx.residue:
        if res.resnum == base_resnum and res.pdbres.strip().startswith(("HI", "HS")):
            hs = {a.pdbname.strip() for a in res.atom} & {"HD1", "HE2"}
            base_state = ("HID" if hs == {"HD1"} else "HIP" if hs == {"HD1", "HE2"}
                          else "HIE" if hs == {"HE2"} else "bare")
            break

    out_mae.parent.mkdir(parents=True, exist_ok=True)
    cx.write(str(out_mae))
    return {"lig_atoms": len(lig), "lig_charge_sum": round(tot, 4), "base_state": base_state}


# -----------------------------------------------------------------------------
# SECTION 12: STAGE 3 - SYSTEM BUILDER (minimise-volume + build + ESP into force field)
# -----------------------------------------------------------------------------
def _lig_atom_indices(complex_mae: Path) -> str:
    from schrodinger import structure
    st = structure.StructureReader.read(str(complex_mae))
    return " ".join(str(a.index) for a in st.atom if a.pdbres.strip() == "LIG")


def run_build(complex_mae: Path, jobname: str, wd: Path) -> Path:
    """Reorient the complex to the smallest orthorhombic box, then solvate/ionise it with Desmond.

    Minimise-volume (the System Builder GUI operation) is a rigid rotation applied before the buffer
    is added, so the box fits the complex and is not oversized (fewer waters, faster MD). ESP charges
    and the NAC geometry are rotation-invariant.
    """
    from schrodinger import structure
    from schrodinger.application.desmond.system_builder_util import DesmondBoxSize
    st = structure.StructureReader.read(str(complex_mae))
    DesmondBoxSize().minimizeVolume([st])
    oriented = wd / f"{jobname}_minvol.mae"
    st.write(str(oriented))
    msj = wd / "build.msj"
    msj.write_text(_build_msj(_lig_atom_indices(oriented)))
    out_cms = wd / f"{jobname}-out.cms"
    _LAUNCHED_JOBS.add(jobname)                       # cancelled on Ctrl-C/kill if it does not finish
    cmd = [f"{SCHRO}/utilities/multisim", "-JOBNAME", jobname, "-HOST", "localhost",
           "-maxjob", "1", "-m", str(msj), "-o", str(out_cms), str(oriented), "-WAIT"]
    """multisim inherits the terminal and writes its own lines at their own indent - the JobId lands
    at column 0, out of step with every other line this phase prints. Stream it instead and re-emit
    each line through _log at the step's indent, which also mirrors multisim's output into the step
    log file. check=True is reproduced explicitly so a build failure still raises."""
    proc = subprocess.Popen(cmd, cwd=str(wd), stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in proc.stdout:
        line = line.strip()
        if line:
            _log(line)
    rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)
    _LAUNCHED_JOBS.discard(jobname)
    return out_cms


def reapply_esp_to_cms(out_cms: Path, esp_mae: Path) -> float:
    """Write the ESP ligand charges into the .cms ffio_block - the charges the MD engine integrates.

    The ffio_block (not the m_atom partial_charge) is what Desmond reads; ffio.site is 1-indexed and
    matched to the ligand atoms by name. After writing, the .cms is re-read through msys (the engine's
    own charges) and any LIG atom still differing from ESP raises - a silent revert to OPLS4 would
    invalidate the defluorination result. Idempotent: re-running only re-verifies an already-ESP .cms.
    """
    from schrodinger import structure
    from schrodinger.application.desmond import cms as cmsmod
    from schrodinger.application.desmond.packages import topo
    esp = {a.pdbname.strip(): a.partial_charge for a in structure.StructureReader.read(str(esp_mae)).atom}
    model = cmsmod.Cms(str(out_cms))

    def _partition_state(m):
        """A healthy Desmond build partitions the full_system atoms across component CTs, so
        sum(comp_ct atoms) == fsys atoms and a solvent CT is present. A build whose partition
        collapsed into the full_system alone (R_1 shipped comp_ct sum 4702 vs fsys 28009, no
        solvent CT) fails both tests."""
        _sum = sum(ct.atom_total for ct in m.comp_ct)
        _has_solvent = any("water" in (ct.title or "").lower() for ct in m.comp_ct)
        return (_sum == m.atom_total and _has_solvent), _sum, _has_solvent

    # Guard the file AS READ, before any write: assert the absolute partition invariant so an already
    # collapsed build (sum(comp_ct) != fsys, or no solvent CT) is refused rather than re-written on resume.
    _ok0, _sum0, _solv0 = _partition_state(model)
    if not _ok0:
        raise RuntimeError(
            f"reapply_esp_to_cms: {out_cms.name} is not a healthy build "
            f"(sum(comp_ct)={_sum0:,} vs fsys={model.atom_total:,}, solvent_ct={_solv0}); "
            f"refusing to touch it. Rebuild the system (restore from the *_3-out.tgz archive).")

    n_lig = applied = 0
    for ct in model.comp_ct:
        if not any(a.pdbres.strip() == "LIG" for a in ct.atom):
            continue
        for a in ct.atom:
            if a.pdbres.strip() == "LIG":
                n_lig += 1
                if a.pdbname.strip() in esp:
                    ct.ffio.site[a.index].charge = esp[a.pdbname.strip()]
                    applied += 1
    # Atomic write, re-asserting the same invariant on the written .tmp. A Cms.write() that
    # flattens the partition leaves an unusable build (MM-GBSA and topo.read_cms iterate comp_ct);
    # the collapsed write is raised and never replaces the good file.
    _tmp = out_cms.with_suffix(out_cms.suffix + ".tmp")
    model.write(str(_tmp))
    _chk = cmsmod.Cms(str(_tmp))
    _ok1, _sum1, _solv1 = _partition_state(_chk)
    if not _ok1 or _chk.atom_total != model.atom_total:
        _tmp.unlink(missing_ok=True)
        raise RuntimeError(
            f"reapply_esp_to_cms: write corrupted the built system "
            f"(sum(comp_ct) {_sum0:,}→{_sum1:,}, fsys {model.atom_total:,}→{_chk.atom_total:,}, "
            f"solvent_ct {_solv0}→{_solv1}); refusing to overwrite {out_cms.name}.")
    _tmp.replace(out_cms)

    msys, cms = topo.read_cms(str(out_cms))
    bad = []
    verified = 0
    lig_sum = 0.0
    for i, a in enumerate(cms.atom):
        if a.pdbres.strip() == "LIG":
            lig_sum += msys.atom(i).charge
            want = esp.get(a.pdbname.strip())
            if want is None:
                continue
            if abs(msys.atom(i).charge - want) > 1e-3:
                bad.append(f"{a.pdbname.strip()} FF={msys.atom(i).charge:+.4f} ESP={want:+.4f}")
            else:
                verified += 1
    # Fail-CLOSED: the whole point of this function is to guarantee ESP reached the force field. A
    # ligand PDB-name mismatch after System Builder writes nothing and verifies nothing; a lenient
    # `if bad` test would pass silently while returning the ESP file's own sum. Require every LIG atom.
    if n_lig == 0:
        raise RuntimeError("reapply_esp_to_cms: no LIG atoms in the built cms.")
    if applied != n_lig or verified != n_lig or bad:
        raise RuntimeError(
            f"ESP charges did NOT reach the MD force field for all {n_lig} LIG atoms "
            f"(applied {applied}, verified {verified})"
            + ("; " + "; ".join(bad) if bad else "")
            + " - likely a ligand PDB-name mismatch after System Builder.")
    return round(lig_sum, 4)   # the sum READ BACK from the engine, not from the ESP file


# -----------------------------------------------------------------------------
# SECTION 13: STAGE 4 - MOLECULAR DYNAMICS
# -----------------------------------------------------------------------------
def _md_msj(time_ps: float, interval_ps: float) -> str:
    relax = Path(f"{SCHRO}/mmshare-v7.3/data/desmond/desmond_npt_relax.msj")
    if not relax.exists():
        # sorted() so a machine with several mmshare-* suites picks the same one every run.
        relax = next(iter(sorted(glob.glob(f"{SCHRO}/mmshare-*/data/desmond/desmond_npt_relax.msj"))), None)
    if not relax or not Path(relax).exists():
        raise RuntimeError("desmond_npt_relax.msj not found - refusing to run production on an "
                           "unequilibrated box (equilibration must precede production).")
    relax_text = Path(relax).read_text()

    # The stock relaxation releases ALL restraints in its final "NPT and no restraints" stage. With the
    # ligand unrestrained for that stage, the substrate relaxes out of the near-attack docked pose
    # (Boltz/PrepWizard place Oδ···Cα in-line at ~170°; the classical OPLS4 minimum is a side-on ~90°
    # contact) BEFORE production's restraint engages, so production then locks the collapsed pose and
    # every downstream frame - including the Step-07 QM/MM starting geometry - is non-reactive. Inject
    # the production restraint block into that final stage so the ligand pose is held continuously from
    # equilibration into production. The stage title is the marker; if a Schrödinger update renames it,
    # fail loudly rather than silently run the collapsing protocol.
    _rblock = _restraint_block()
    if _rblock:
        _title_re = re.compile(r'(title\s*=\s*")NPT and no restraints([^"]*")')
        if not _title_re.search(relax_text):
            raise RuntimeError("desmond_npt_relax.msj: final unrestrained stage not found by title - "
                               "cannot inject the ligand pose restraint; refusing to run the pose-"
                               "collapsing relaxation protocol.")
        relax_text = _title_re.sub(rf'\g<1>NPT with ligand + backbone pose restraint\g<2>\n{_rblock.rstrip()}',
                                   relax_text, count=1)
    return relax_text + _md_production(time_ps=time_ps, interval_ps=interval_ps)


def _extract_with_progress(tgz: Path, wd: Path, members: list, label: str) -> None:
    """Extract members from tgz into wd (dropping one leading path component), drawing the shared \r
    progress line from tar's read offset into the archive (/proc/PID/fdinfo pos), so a multi-GB
    trajectory unpack reports live GB-consumed / GB-total. Raises on a non-zero tar exit so a truncated
    or corrupt archive is caught exactly as `tar ... check=True` would."""
    total = max(1, tgz.stat().st_size)
    proc = subprocess.Popen(["tar", "xzf", str(tgz), "-C", str(wd), "--strip-components=1", *members])
    _tgz_real = os.path.realpath(str(tgz))

    def _pos_for_pid(_pid: str) -> "int | None":
        """Read offset into the archive from any fd of _pid that points at the tgz. `tar xzf` decodes
        through a gzip child that holds the compressed file, so the live offset lives on the child."""
        _fdd = f"/proc/{_pid}/fd"
        try:
            _fds = os.listdir(_fdd)
        except OSError:
            return None
        for _fd in _fds:
            try:
                if os.path.realpath(os.path.join(_fdd, _fd)) != _tgz_real:
                    continue
                with open(f"/proc/{_pid}/fdinfo/{_fd}") as _fh:
                    for _ln in _fh:
                        if _ln.startswith("pos:"):
                            return int(_ln.split()[1])
            except (OSError, ValueError):
                continue
        return None

    def _read_pos() -> "int | None":
        _pids = [str(proc.pid)]
        try:
            _pids += open(f"/proc/{proc.pid}/task/{proc.pid}/children").read().split()
        except OSError:
            pass
        for _pid in _pids:
            _p = _pos_for_pid(_pid)
            if _p is not None:
                return _p
        return None

    _worker_progress(f"Extracting the MD Simulation file for {label}: {total / 1e9:.1f} GB archive - unpacking…")
    _last_pct, _last_t = -100, 0.0
    while proc.poll() is None:
        _pos = _read_pos()
        _now = time.time()
        if _pos is not None:
            _pos = min(_pos, total)
            _pct = 100 * _pos // total
            if _pct - _last_pct >= 5 or _now - _last_t >= 20:   # throttle: every ~5% or ~20 s
                _worker_progress(f"Extracting the MD Simulation file for {label}: "
                                 f"{_pos / 1e9:.1f}/{total / 1e9:.1f} GB ({_pct}%)")
                _last_pct, _last_t = _pct, _now
        time.sleep(1.0)
    if proc.returncode != 0:
        raise RuntimeError(f"tar extract failed (rc={proc.returncode}) for {tgz.name}")
    _worker_progress(f"Extracting the MD Simulation file for {label}: {total / 1e9:.1f} GB unpacked - done.")


def _unpack_md_production(wd: Path, jobname: str, expected_time_ps: "float | None" = None) -> None:
    """Put the production trajectory at the job-dir root as {job}_trj + {job}.ene and repoint the cms.

    multisim leaves the production stage (the final `simulate`) in one of two shapes, both the same
    finished run: loose at the job-dir root as {job}_trj + {job}.ene, or - after its own 'Cleaning up
    files' pass - packed into the highest-numbered stage archive {job}_N-out.tgz. The earlier archives
    are equilibration stages (the last a restrained 24 ps NPT relax); installing one of those as
    production would silently score an unequilibrated trajectory. A packed archive is therefore accepted
    only when its own cfg reports last_time == the requested production length, which no relax stage
    does. 07 and thermal_mmgbsa read the trajectory flat, so the accepted stage is unpacked to the root
    and its staged name flattened. If neither a loose trajectory nor a matching production archive is
    present the MD did not finish, so raise rather than fall back to an equilibration stage.
    """
    cms = wd / f"{jobname}-out.cms"

    def _repoint_cms() -> None:
        """Repoint the cms's s_chorus_trajectory_file from the staged name ({job}_N_trj) to {job}_trj.
        multisim writes the production-stage trajectory name into the cms; thermal_mmgbsa (Prime
        MM-GBSA) and 07 both resolve the trajectory from this reference, so a stale name aborts MM-GBSA
        with 'No trajectory found associated with CMS'."""
        if not cms.exists():
            return
        try:
            b = cms.read_bytes()
            m = re.search(re.escape(jobname.encode()) + rb"_\d+_trj", b)
            if m and not (wd / m.group(0).decode()).exists():
                cms.write_bytes(b.replace(m.group(0), f"{jobname}_trj".encode()))
        except Exception:
            pass

    # Loose production trajectory already at the root - the common case, and the resume no-op.
    if (wd / f"{jobname}_trj").exists():
        _repoint_cms()
        _worker_progress(f"{jobname}: production trajectory already unpacked - skipping extraction.")
        return

    """Production packed into a stage archive. Take the highest-numbered stage (production is the last
    simulate), but install it only after its cfg confirms it is the full-length run - an equilibration
    tgz can never stand in for production. The cfg sits ahead of the frames in the stream, so reading it
    with --occurrence=1 stops tar early rather than decompressing the whole multi-GB archive."""
    def _stage_num(p: Path) -> int:
        return int(re.search(rf"{re.escape(jobname)}_(\d+)-out\.tgz$", p.name).group(1))
    for tgz in sorted(wd.glob(f"{jobname}_*-out.tgz"), key=_stage_num, reverse=True):
        stage = f"{jobname}_{_stage_num(tgz)}"
        try:
            cfg = subprocess.run(["tar", "xzf", str(tgz), f"{stage}/{stage}-out.cfg",
                                  "--occurrence=1", "-O"], capture_output=True, timeout=300).stdout
        except (subprocess.SubprocessError, OSError):
            continue
        m = re.search(rb"last_time\s*=\s*([0-9.]+)", cfg)
        if not m:
            continue
        if expected_time_ps is not None and abs(float(m.group(1)) - expected_time_ps) > 1.0:
            continue                                     # an equilibration stage, not production
        # --strip-components=1 drops the {stage}/ parent so the trajectory lands directly in the job dir
        # rather than a nested {stage}/ subdir; only the {stage}_N basename is then flattened to {job}.
        # A non-zero tar exit (gzip CRC / truncated member) raises inside _extract_with_progress, so a
        # corrupt archive is caught and the production tgz below is never removed.
        _extract_with_progress(tgz, wd, [f"{stage}/{stage}_trj", f"{stage}/{stage}.ene"], jobname)
        (wd / f"{stage}_trj").rename(wd / f"{jobname}_trj")
        ene = wd / f"{stage}.ene"
        if ene.exists():
            ene.rename(wd / f"{jobname}.ene")
        # Leave the job dir as Maestro does: production loose at the root, the archive gone. The small
        # equilibration archives ({job}_1..6-out.tgz) stay - Maestro keeps those too. Removed only after
        # the verified extract above, so the loose trajectory is never deleted without its replacement.
        tgz.unlink(missing_ok=True)
        _repoint_cms()
        return

    raise RuntimeError(
        f"{jobname}: no production trajectory - neither a loose {jobname}_trj at the job-dir root nor a "
        f"stage archive whose last_time matches the requested production length ({expected_time_ps} ps). "
        f"The MD did not finish; re-run it for this rank.")


def run_md(system_cms: Path, jobname: str, wd: Path, time_ns: float, frames: int) -> Path:
    time_ps = time_ns * 1000.0
    interval_ps = max(round(time_ps / max(frames, 1), 4), 1.0)
    msj = wd / "md.msj"
    msj.write_text(_md_msj(time_ps, interval_ps))
    out_cms = wd / f"{jobname}-out.cms"
    _LAUNCHED_JOBS.add(jobname)                           # cancelled on Ctrl-C/kill if it does not finish
    # multisim blocks on -WAIT; the heartbeat thread tails the live .ene alongside it and refreshes a
    # single \r line - current ns / total ns, %, ns/day and ETA - so a 1000 ns run reads as live progress.
    with MDHeartbeat(wd, wd / f"{jobname}_multisim.log", f"MD {jobname}", time_ns):
        subprocess.run([f"{SCHRO}/utilities/multisim", "-JOBNAME", jobname, "-HOST", "localhost",
                        "-SUBHOST", "localhost", "-maxjob", "1", "-m", str(msj), "-o", str(out_cms),
                        str(system_cms), "-WAIT"], cwd=str(wd), check=True)
    _LAUNCHED_JOBS.discard(jobname)
    return out_cms                                        # raw multisim output; extraction runs on the CPU worker


# -----------------------------------------------------------------------------
# SECTION 14: STAGE 2 - WATERMAP (+ CSV export for Step 07)
# -----------------------------------------------------------------------------
def run_watermap(complex_mae: Path, jobname: str, wd: Path, time_ns: float, lig_dist: float) -> Path:
    """Run WaterMap through Schrödinger's own WaterMapInput, on a local scratch, and copy the result back.

    WaterMapInput splits the complex into receptor + ligand, tags the ct types, truncates the protein
    to the active site and writes the GPU msj (S-OPLS ligand, TIP4P water, holo via retain_ligand). It
    is driven in-process with real floats because the watermap_inp CLI would multiply a string
    simulation_time. The job runs in a /tmp scratch under a minimal environment (see _clean_job_env and
    Critic's Corner note 3); a WaterMap-GPU licence shortfall is retried from the checkpoint.
    """
    from schrodinger import structure
    from schrodinger.application.watermap import watermap_inp
    cx = structure.StructureReader.read(str(complex_mae))
    lig_idx = [a.index for a in cx.atom if a.pdbres.strip() == "LIG"]
    if not lig_idx:
        raise RuntimeError(f"no LIG atoms in {complex_mae.name} for WaterMap")
    ligand = cx.extract(lig_idx)
    protein = cx.copy()
    protein.deleteAtoms(lig_idx)
    wm = watermap_inp.WaterMapInput(protein, ligand, retain_ligand=bool(CFG.PHYS_WM_RETAIN_LIGAND),
                                    ligand_distance=float(lig_dist), simulation_time=float(time_ns))

    scratch = Path(tempfile.mkdtemp(prefix=f"{jobname}_wm_"))
    cwd = os.getcwd()
    ok = False
    try:
        os.chdir(str(scratch))
        wm.write(jobname)                             # writes {jobname}-in.maegz + {jobname}.msj
        os.chdir(cwd)
        in_mae, msj = scratch / f"{jobname}-in.maegz", scratch / f"{jobname}.msj"
        src = scratch / f"{jobname}_wm.maegz"
        log = scratch / f"{jobname}_multisim.log"
        ckpt = scratch / f"{jobname}-multisim_checkpoint"
        last_tail = ""
        max_tries = CFG.PHYS_WM_MAX_TRIES
        _LAUNCHED_JOBS.add(jobname)                    # cancelled on Ctrl-C/kill if it does not finish
        for attempt in range(1, max_tries + 1):
            if attempt > 1 and ckpt.exists():
                cmd = [f"{SCHRO}/watermap", "-JOBNAME", jobname, "-HOST", "localhost",
                       "-RESTART", ckpt.name, "-WAIT"]
            else:
                cmd = [f"{SCHRO}/watermap", "-JOBNAME", jobname, "-HOST", "localhost",
                       "-m", msj.name, in_mae.name, "-WAIT"]   # basenames: cwd is the scratch
            # Launch through a bash login shell with a detached environment (see _clean_job_env), and
            # with the input given by BASENAME from the scratch cwd. Both matter: a job that inherits the
            # $SCHRODINGER/run session, or is handed an absolute input path, fails to stage its GCMC ligand
            # companion to the sub-job and dies at stage 8.
            try:
                subprocess.run(["bash", "-lc", " ".join(shlex.quote(c) for c in cmd)],
                               cwd=str(scratch), check=True, env=_clean_job_env())
            except subprocess.CalledProcessError:
                pass                                  # the driver may exit non-zero; the maegz is authoritative
            if src.exists():
                break
            # WaterMap is stochastic (GCMC) and licence-contended: the SAME input often succeeds on a
            # rerun. Retry up to 3 times, resuming from the checkpoint; licence shortfalls wait longer.
            lines = log.read_text().splitlines() if log.exists() else []
            last_tail = "\n".join(lines[-30:]) if lines else "(no multisim log)"
            if attempt < max_tries:
                wait = 240 if any(("not enough total licenses" in ln) or ("FAILED_PRECONDITION" in ln)
                                  for ln in lines) else 60
                _log(f"[watermap] {jobname}: attempt {attempt}/{max_tries} produced no output - "
                     f"retrying in {wait} s")
                time.sleep(wait)
        if log.exists():
            shutil.copy2(log, wd / f"{jobname}_multisim.log")
        if not src.exists():
            lines = log.read_text().splitlines() if log.exists() else []
            key = next((ln.strip() for ln in reversed(lines)
                        if any(t in ln for t in ("not enough total licenses", "not found", "ERROR:",
                                                 "Error:", "Exception", "FAILED"))),
                       "(cause not in log - see scratch)")
            raise RuntimeError(f"WaterMap produced no {src.name}.\n"
                               f"    cause  : {key}\n"
                               f"    scratch: {scratch}  (kept for debug)\n"
                               f"    log    :\n{last_tail}")
        _LAUNCHED_JOBS.discard(jobname)                # server job has finished (WAIT returned with a maegz)
        out = wd / f"{jobname}_wm.maegz"
        shutil.copy2(src, out)
        if in_mae.exists():
            # Retain the WaterMap INPUT frame alongside the output. Its Cα are the reference the static
            # hydration sites are aligned onto before use (Step 07 load_watermap_reference_ca + frame
            # superposition); the *_wm.maegz output carries only sites, no protein backbone, so without
            # this the sites cannot be superimposed onto any MD frame and the WaterMap guidance no-ops.
            shutil.copy2(in_mae, wd / in_mae.name)
        cluster = scratch / f"{jobname}-cluster.maegz"
        if cluster.exists():
            shutil.copy2(cluster, wd / cluster.name)
        ok = True
        return out
    finally:
        os.chdir(cwd)
        if ok:
            shutil.rmtree(scratch, ignore_errors=True)


def export_watermap_csv(wm_maegz: Path, csv_path: Path) -> int:
    """Export the hydration sites to the CSV Step 07 reads: Site, Occupancy, dH, -TdS, dG, #HB(WW/PW/LW).

    Property names are the actual keys on a completed *_wm.maegz site structure:
      dG=r_watermap_deltaG  dH=r_watermap_deltaH  -TdS=r_watermap_-TdeltaS  occupancy=r_watermap_occupancy
      H-bonds=r_watermap_hbond_ww/pw/lw  site=i_watermap_site_num  (dH + (-TdS) = dG holds).
    A site is any atom carrying dG; atoms without it are skipped.
    """
    from schrodinger import structure
    _DG = ("r_watermap_deltaG", "r_watermap_free_energy")
    st = None
    for s in structure.StructureReader(str(wm_maegz)):
        if any(any(a.property.get(k) is not None for k in _DG) for a in s.atom):
            st = s
            break
    if st is None:
        return 0

    def _p(a, *keys):
        for k in keys:
            v = a.property.get(k)
            if v is not None:
                return v
        return None

    rows, n = [], 0
    for a in st.atom:
        dg = _p(a, *_DG)
        if dg is None:
            continue
        n += 1
        rows.append({
            "Site": _p(a, "i_watermap_site_num", "i_watermap_sitenum") or n,
            "Occupancy": _p(a, "r_watermap_occupancy", "r_watermap_density"),
            "dH (kcal/mol)": _p(a, "r_watermap_deltaH", "r_watermap_potential_energy_relative"),
            "-TdS (kcal/mol)": _p(a, "r_watermap_-TdeltaS", "r_watermap_entropy"),
            "dG (kcal/mol)": dg,
            "#HB(WW)": _p(a, "r_watermap_hbond_ww"),
            "#HB(PW)": _p(a, "r_watermap_hbond_pw"),
            "#HB(LW)": _p(a, "r_watermap_hbond_lw"),
        })
    # Atomic write: a kill between copy2(maegz) and this open() must not leave a truncated CSV that the
    # existence-only resume gate then treats as done. tmp + replace, so the CSV appears whole or not at all.
    _tmp = csv_path.with_suffix(csv_path.suffix + ".tmp")
    with open(_tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["Site"])
        w.writeheader()
        w.writerows(rows)
    _tmp.replace(csv_path)
    return n



# =============================================================================
# SECTION 15: MAIN
# =============================================================================
def plot_watermap_combined(wm_dir: Path, analysis_dir: Path) -> None:
    """One compact multi-panel figure of every rank's active-site hydration landscape, drawn straight
    from the WaterMaps already exported (watermap_R_<rank>.csv). Each panel sorts the sites by free
    energy and colours displaceable (ΔG>0) vs stable (ΔG<0) waters, annotating the counts and the mean
    ΔG - so the pocket a substrate must dewet to bind reads at a glance across the whole lead set.
    Redrawn as each WaterMap lands (it reads only the CSVs present), and skipped without failing the
    stage when none are ready. Saved to 06_Analysis/02_WaterMap_Landscapes_AllRanks.svg."""
    try:
        _csvs = sorted(wm_dir.glob("watermap_R_*.csv"),
                       key=lambda p: int(re.search(r"watermap_R_(\d+)\.csv$", p.name).group(1)))
        _csvs = [c for c in _csvs if re.search(r"watermap_R_(\d+)\.csv$", c.name)]
        if not _csvs:
            return
        _stable = CFG.VIS_ACCENT["blue"]; _disp = CFG.VIS_ACCENT["amber"]
        _ft = CFG.VIS_FONT_TICK
        # rank -> ligand name (nice panel titles), from the ranked sheet when available
        _lig = {}
        try:
            _prod = wm_dir.parent.parent / "1_Boltz2_Production"
            _rk = (_utils_mod.latest_by_mtime(_prod.glob(CFG.GLOB_RANKED_CSV))
                   or _utils_mod.latest_by_mtime(_prod.glob("*Ranked*.csv"))) if _prod.is_dir() else None
            if _rk is not None:
                _rdf = pd.read_csv(_rk, low_memory=False)
                _nc = next((c for c in ("Ligand_Name", "ligand_name", "job_name") if c in _rdf.columns), None)
                if "Scientific_Rank" in _rdf.columns and _nc:
                    for _, _r in _rdf.iterrows():
                        try:
                            _lig[int(_r["Scientific_Rank"])] = str(_r[_nc]).split("_")[-1]
                        except Exception:
                            pass
        except Exception:
            pass
        _panels = []
        for _c in _csvs:
            try:
                _df = pd.read_csv(_c)
                _col = next((k for k in _df.columns if k.strip().lower().startswith("dg")), None)
                if _col is None:
                    continue
                _dg = np.sort(pd.to_numeric(_df[_col], errors="coerce").dropna().to_numpy())
                if _dg.size:
                    _panels.append((int(re.search(r"_R_(\d+)\.csv$", _c.name).group(1)), _dg))
            except Exception:
                continue
        if not _panels:
            return
        _n = len(_panels); _ncol = 2 if _n > 1 else 1; _nrow = (_n + _ncol - 1) // _ncol
        # share the y-axis within a row: one ΔG label + scale per row (left panel), x-axis stays per panel
        fig, axes = plt.subplots(_nrow, _ncol, figsize=(6.0 * _ncol, 3.2 * _nrow), squeeze=False, sharey="row")
        for _i, (_rank, _dg) in enumerate(_panels):
            ax = axes[_i // _ncol][_i % _ncol]
            _x = np.arange(_dg.size)
            ax.bar(_x, _dg, color=[_stable if v < 0 else _disp for v in _dg], width=0.9, linewidth=0)
            ax.axhline(0, color=CFG.VIS_INK["near_black"], lw=0.8)
            _nu = int((_dg > 0).sum())
            _head = f"R{_rank}" + (f" · {_lig[_rank]}" if _rank in _lig else "")
            ax.text(0.03, 0.96, f"{_head}\n{_nu} displaceable / {_dg.size} sites\nmean ΔG {float(_dg.mean()):.2f} kcal/mol",
                    transform=ax.transAxes, ha="left", va="top", fontsize=_ft, linespacing=1.35,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.7))
            ax.set_xlabel("hydration site (sorted by ΔG)", fontsize=_ft)     # x per panel
            if _i % _ncol == 0:
                ax.set_ylabel("water ΔG  (kcal/mol)", fontsize=_ft)          # one y label per row (leftmost)
            else:
                ax.tick_params(labelleft=False)
            for _s in ("top", "right"):
                ax.spines[_s].set_visible(False)
        for _j in range(_n, _nrow * _ncol):
            axes[_j // _ncol][_j % _ncol].axis("off")
        _lh = [Patch(fc=_disp, ec="none", label="displaceable (ΔG > 0)"),
               Patch(fc=_stable, ec="none", label="stable (ΔG < 0)")]
        axes[0][0].legend(handles=_lh, loc="upper right", frameon=True, framealpha=0.85,
                          edgecolor="#C8C8C8", fontsize=_ft - 1, ncol=1, borderpad=0.5)   # inside 1st panel
        fig.tight_layout()
        _out = analysis_dir / "02_WaterMap_Landscapes_AllRanks.svg"
        analysis_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(_out, dpi=int(CFG.VIS_FIGURE_DPI), bbox_inches="tight")
        plt.close(fig)
        _echo(f"  ✔ WaterMap combined figure → {analysis_dir.name}/{_utils_mod.deflx_fig_name('02_WaterMap_Landscapes_AllRanks.svg')} ({_n} ranks)")
    except Exception as _e:
        _warn(f"[watermap] combined landscape figure skipped ({_e}).")


def _parse_args_merged():
    ap = argparse.ArgumentParser(
        description=f"{CFG.PROJECT_NAME} Step 06 - ESP Physics: WaterMap → System Builder → MD → SID → MM-GBSA → Defluorination",
        epilog=("examples:\n"
                "  production : python 06_Physics_Validation_DeFluorX.py Boltz-2_Run_20260309T085406Z/\n"
                "  quick test : python 06_Physics_Validation_DeFluorX.py Boltz-2_Run_20260309T085406Z/ --test"),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_folder", nargs="?", default=None,
                    help="Boltz-2_Run_* directory (auto-detect latest if omitted)")
    ap.add_argument("--test", action="store_true",
                    help=f"quick run: WaterMap {CFG.PHYS_TEST_WM_NS:g} ns · MD {CFG.PHYS_TEST_MD_NS:g} ns/{CFG.PHYS_TEST_MD_FRAMES} frames")
    ap.add_argument("--md-ns", type=float, default=CFG.PHYS_MD_NS)
    ap.add_argument("--md-frames", type=int, default=CFG.PHYS_MD_FRAMES)
    ap.add_argument("--wm-ns", type=float, default=CFG.PHYS_WM_NS)
    ap.add_argument("--lig-dist", type=float, default=CFG.PHYS_WM_SITE_A)
    ap.add_argument("--stages", default=",".join(STAGES_ALL),
                    help="subset of {merge,watermap,build,md}; SID+MM-GBSA run sequentially inside the md stage")
    ap.add_argument("--out", default=None)
    ap.add_argument("--pipeline-mode", action="store_true",
                    help="called from 00_00_run_pipeline_DeFluorX.sh; delegates oomd masking to the runner.")
    a = ap.parse_args()
    if a.test:
        a.md_ns, a.md_frames, a.wm_ns = CFG.PHYS_TEST_MD_NS, CFG.PHYS_TEST_MD_FRAMES, CFG.PHYS_TEST_WM_NS
    stages = [s.strip() for s in a.stages.split(",") if s.strip()]
    unknown = [s for s in stages if s not in STAGES_ALL]
    if unknown:
        sys.exit(f"Unknown --stages token(s): {', '.join(unknown)}  (valid: {', '.join(STAGES_ALL)})")
    a.stage_set = set(stages)
    return a


def _complex_label(entry: dict, ranked_map: dict) -> str:
    row = ranked_map.get(entry["stem"], {})
    return (f"R_{entry['rank']}  {row.get('Ligand_Name', '') or entry['stem']} × "
            f"{row.get('Protein_Name', '')}").rstrip(" ×")


def _hydrate(entry: dict, ranked_map: dict, dirs: dict) -> None:
    """Populate the keys every phase reads - esp / base / complex_mae - from pure lookups, with NO side
    effects (no copy, no merge). _phase_merge sets these as a by-product, so a --stages run that skips
    'merge' would otherwise KeyError inside build/md/watermap; run this up front for every entry instead."""
    stem, rank = entry["stem"], entry["rank"]
    if entry.get("esp") is None:
        entry["esp"] = find_esp(dirs["esp"], stem)
    row = ranked_map.get(stem)
    if row is None:
        raise RuntimeError(f"no ranked-CSV row (job_name) matches handover stem {stem}")
    entry.setdefault("base", _resnum(row.get("Mapped_Base", "")))
    entry.setdefault("complex_mae", dirs["esp_cx"] / f"R_{rank}_{stem}_ESP_Complex.mae")


def _read_physics_qc(entry: dict, dirs: dict) -> "dict | None":
    """Read the built .cms + WaterMap .maegz for one rank into the numbers the QC figure plots.

    Two reads of the same built system: schrodinger.structure/Cms for composition, box and geometry
    (comp_ct, water/ion counts, per-atom coordinates), and topo.read_cms for the msys partial charges
    the MD engine actually integrates (the ESP charges written into the ffio force field). The α-carbon
    is the ligand carbon bonded to fluorine - the reactive centre - and the 'crucial' hydration site is
    the WaterMap site nearest that α-carbon. Returns None (never raises) if a file is missing/unreadable.
    """
    from schrodinger import structure
    from schrodinger.application.desmond import cms as _cmsmod
    from schrodinger.application.desmond.packages import topo

    rank = entry["rank"]
    cms_path = entry.get("setup_cms") or (dirs["sb"] / f"desmond_setup_R_{rank}" / f"desmond_setup_R_{rank}-out.cms")
    maegz = dirs["wm"] / f"watermap_R_{rank}" / f"watermap_R_{rank}_wm.maegz"
    if not Path(cms_path).exists():
        return None

    m = _cmsmod.Cms(str(cms_path))
    d: dict = {"rank": rank}
    titles = [ct.title for ct in m.comp_ct]
    d["comp_ct"] = len(m.comp_ct)
    d["atoms_total"] = m.atom_total
    d["n_water"] = sum(ct.atom_total for ct, t in zip(m.comp_ct, titles) if "water" in (t or "").lower()) // 3
    d["n_na"] = sum(1 for t in titles if t == "Na+")
    d["n_cl"] = sum(1 for t in titles if t == "Cl-")
    try:
        d["box_A"] = [round(m.box[i], 1) for i in (0, 4, 8)]
    except Exception:
        d["box_A"] = None

    # Ligand geometry (for the α-C coordinate) from the structure-level atoms.
    ligC_xyz, F_xyz = [], []
    for ct in m.comp_ct:
        for a in ct.atom:
            if a.pdbres.strip() == "LIG":
                if a.element == "C":
                    ligC_xyz.append((a.x, a.y, a.z))
                elif a.element == "F":
                    F_xyz.append((a.x, a.y, a.z))
    aC_xyz = (min(ligC_xyz, key=lambda cc: min(_math.dist(cc, f) for f in F_xyz))
              if (F_xyz and ligC_xyz) else None)

    # Partial charges the engine integrates (msys), per ligand atom.
    msys, cms = topo.read_cms(str(cms_path))
    lig_at = [(a.element, msys.atom(i).charge) for i, a in enumerate(cms.atom) if a.pdbres.strip() == "LIG"]
    d["lig_sumq"] = round(sum(q for _, q in lig_at), 4)
    d["lig_natoms"] = len(lig_at)
    Cs = sorted([q for e, q in lig_at if e == "C"], reverse=True)
    d["carboxyl_C"] = round(Cs[0], 4) if Cs else None       # most positive C = the carboxylate carbon
    d["alpha_C"] = round(Cs[1], 4) if len(Cs) > 1 else None  # next = the reactive α-carbon (bears the F)
    d["F_charges"] = [round(q, 4) for e, q in lig_at if e == "F"]
    d["O_charges"] = [round(q, 4) for e, q in lig_at if e == "O"]

    # WaterMap hydration sites (ΔG per site) + the crucial water nearest the α-carbon.
    d["dG_list"], d["n_sites"], d["n_unstable"] = [], 0, 0
    d["crucial_dG"], d["crucial_dist_A"] = None, None
    if maegz.exists():
        sites = []
        for st in structure.StructureReader(str(maegz)):
            for a in st.atom:
                dg = a.property.get("r_watermap_deltaG", a.property.get("r_watermap_free_energy"))
                if dg is not None:
                    sites.append((round(float(dg), 2), (a.x, a.y, a.z)))
        d["dG_list"] = [s[0] for s in sites]
        d["n_sites"] = len(sites)
        d["n_unstable"] = sum(1 for g, _ in sites if g > 0)
        if aC_xyz and sites:
            near = min(sites, key=lambda s: _math.dist(s[1], aC_xyz))
            d["crucial_dG"] = near[0]
            d["crucial_dist_A"] = round(_math.dist(near[1], aC_xyz), 2)
    return d


def _analysis_dir(physics_root: Path) -> Path:
    """The one folder every Step-06 figure is written to - build/solvation QC, MD trajectory QC, and the
    MM-GBSA plots - so all analysis figures for the run sit together rather than scattered across subdirs."""
    d = physics_root / "06_Analysis"
    d.mkdir(parents=True, exist_ok=True)
    return d


def make_physics_qc_figure(entries: list, dirs: dict, out_root: Path, ligands: dict, controls: set = None) -> None:
    """Draw the build/solvation QC figure - ESP charge gradient, per-atom charges, WaterMap ΔG, box size.

    One 4-panel snapshot of what the physics build produced for every MD-selected rank, from the numbers
    the MD engine integrates (ESP charges in the ffio force field) and the WaterMap thermodynamics.
    Best-effort: any missing/unreadable system is dropped, and the whole figure is skipped rather than
    failing the run. Written to 6_Physics_Validation/06_Analysis/.
    """
    recs = []
    for e in entries:
        if e.get("_skip"):
            continue
        try:
            r = _read_physics_qc(e, dirs)
        except Exception as exc:
            _warn(f"[qc] R_{e['rank']} physics read failed - {str(exc).splitlines()[0]}")
            r = None
        if r:
            recs.append(r)
    if not recs:
        _warn("[qc] no built systems readable - physics QC figure skipped.")
        return

    recs.sort(key=lambda r: r["rank"])
    ranks = [r["rank"] for r in recs]
    controls = controls or set()
    labs = [_ctrl_label(rk, ligands, controls) for rk in ranks]
    _A, _INK = CFG.VIS_ACCENT, CFG.VIS_INK
    palette = [_A["green"], _A["amber"], _A["vermillion"], _A["blue"], _A["magenta"], _A["sky"]]
    cols = _ctrl_palette(ranks, controls, palette)
    xp = np.arange(len(recs))
    fa = CFG.VIS_FONT_ANNOT

    def _colour_xticks(axobj):
        for t, c in zip(axobj.get_xticklabels(), cols):
            t.set_color(c); t.set_fontweight("bold")

    fig, ax = plt.subplots(2, 2, figsize=(12, 8.6))

    # A - α-carbon ESP charge gradient (the reactive centre, as the MD force field sees it).
    a = ax[0, 0]
    aC = [r["alpha_C"] for r in recs]
    if all(v is not None for v in aC):
        a.bar(xp, aC, color=cols, edgecolor=_INK["dark"], linewidth=0.9, width=0.6, zorder=3)
        for xi, v in zip(xp, aC):
            a.text(xi, v + 0.006, f"{v:+.3f}", ha="center", fontsize=fa, fontweight="bold")
        a.set_ylim(-0.02, max(aC) * 1.28)
    a.axhline(0, color=_INK["soft"], lw=0.8)
    a.set_xticks(xp); a.tick_params(labelbottom=False)
    _qc_group_seps(a, len(recs))
    a.set_ylabel("α-carbon ESP charge (e)\nreactive centre, in the MD force field")

    # B - ligand per-atom ESP charges (α-C, carboxyl-C, F, O); Σq = −1.000.
    bx = ax[0, 1]
    for i, r in enumerate(recs):
        pts = [v for v in (r["alpha_C"], r["carboxyl_C"]) if v is not None] + r["F_charges"] + r["O_charges"]
        xs = [i + (_j - len(pts) / 2 + 0.5) * 0.12 for _j in range(len(pts))]
        bx.scatter(xs, pts, color=cols[i], s=48, edgecolor=_INK["dark"], linewidth=0.6, zorder=3)
    bx.set_xticks(xp); bx.tick_params(labelbottom=False)
    bx.axhline(0, color=_INK["soft"], lw=0.8)
    _qc_group_seps(bx, len(recs))
    bx.set_ylabel("Ligand per-atom ESP charge (e)\nα-C · carboxyl-C (+) · F · O (−);  Σq = −1.000")

    # C - WaterMap hydration-site ΔG (>0 = displaceable water); ★ = the water nearest the α-carbon.
    cx = ax[1, 0]
    data = [r["dG_list"] for r in recs]
    if any(data):
        vp = cx.violinplot(data, positions=xp, showmedians=True, widths=0.72)
        for i, pc in enumerate(vp["bodies"]):
            pc.set_facecolor(cols[i]); pc.set_alpha(0.5); pc.set_edgecolor(_INK["dark"])
        for k in ("cbars", "cmins", "cmaxes", "cmedians"):
            if k in vp:
                vp[k].set_color(_INK["dark"]); vp[k].set_linewidth(1)
        ymax = max((max(dl) for dl in data if dl), default=1.0)
        ymin = min((min(dl) for dl in data if dl), default=-1.0)
        for i, r in enumerate(recs):
            cx.text(i, ymax + 1.3, f"{r['n_sites']} sites\n{r['n_unstable']} displaceable",
                    ha="center", va="bottom", fontsize=fa - 0.5, fontweight="bold", color=cols[i])
            cw = r.get("crucial_dG")
            if cw is not None:
                cx.scatter([i], [cw], marker="*", s=240, color=_A["star"],
                           edgecolor=_INK["dark"], linewidth=1.1, zorder=6)
                cx.annotate(f"crucial H₂O\n{cw:+.1f}, {r['crucial_dist_A']} Å", (i, cw),
                            xytext=(i + 0.28, cw + 0.4), fontsize=fa - 1.5, color=_INK["dark"])
        cx.set_ylim(ymin - 1.2, ymax + 3)
        cx.axhline(0, color=_INK["dark"], lw=1.0, ls="--")
        cx.text(len(recs) - 0.52, 0.15, "ΔG = 0", ha="right", va="bottom",
                fontsize=fa - 1, color=_INK["soft"])
        cx.text(len(recs) - 0.52, -0.15, "stable ↓ · displaceable ↑", ha="right", va="top",
                fontsize=fa - 1.5, color=_INK["soft"])
        cx.text(0.015, 0.985, "★ crucial H₂O (nearest reactive α-C)", transform=cx.transAxes,
                ha="left", va="top", fontsize=fa - 1.5, color=_INK["soft"])
    cx.set_xticks(xp); cx.set_xticklabels(labs); _colour_xticks(cx)
    _qc_group_seps(cx, len(recs))
    cx.set_ylabel("WaterMap hydration-site ΔG (kcal/mol)")

    # D - solvated-system size: water count vs total atoms, box dims annotated.
    dx = ax[1, 1]
    w = 0.36
    dx.bar(xp - w / 2, [r["n_water"] for r in recs], w, label="water molecules",
           color=_A["blue"], edgecolor=_INK["dark"], linewidth=0.6, zorder=3)
    dx.bar(xp + w / 2, [r["atoms_total"] for r in recs], w, label="total system atoms",
           color=_INK["soft"], edgecolor=_INK["dark"], linewidth=0.6, zorder=3)
    for i, r in enumerate(recs):
        dx.text(i - w / 2, r["n_water"] + 400, f"{r['n_water']:,}", ha="center", fontsize=fa - 1)
        _box = (f"\n{r['box_A'][0]:.0f}×{r['box_A'][1]:.0f}×{r['box_A'][2]:.0f} Å" if r["box_A"] else "")
        dx.text(i + w / 2, r["atoms_total"] + 400, f"{r['atoms_total']:,}{_box}",
                ha="center", fontsize=fa - 1.5)
    dx.set_xticks(xp); dx.set_xticklabels(labs); _colour_xticks(dx)
    _qc_group_seps(dx, len(recs))
    dx.set_ylim(0, max(r["atoms_total"] for r in recs) * 1.18)
    dx.set_ylabel("Solvated-system size (atom / water count)\ncomp_ct = 5 · 2 Na⁺ · 1 Cl⁻ (all)")
    dx.legend(loc="upper left", fontsize=CFG.VIS_FONT_LEGEND)

    fig.tight_layout()
    qc_dir = _analysis_dir(out_root)
    out_path = qc_dir / "01_Physics_Build_Solvation_QC.svg"
    plt.savefig(out_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    _ok(f"[qc] ✔ physics build/solvation QC → {qc_dir.name}/{_utils_mod.deflx_fig_name(out_path.name)}  ({len(recs)} rank(s))")


def _read_md_qc(job_dir: Path, jobname: str) -> "dict | None":
    """Read one finished MD job's SID .eaf + Desmond .ene into the trajectory-QC arrays.

    From the SID event-analysis file: protein Cα-RMSD and ligand RMSD (fit on protein) per frame, and
    per-residue Cα-RMSF. From the .ene: system temperature per step. Returns None (never raises) if the
    SID .eaf is absent (SID not yet run) - the .ene panel degrades gracefully to empty.
    """
    from schrodinger.utils import sea

    eaf = job_dir / f"{jobname}{CFG.SUFFIX_SID_EAF}"
    ene = job_dir / f"{jobname}.ene"
    if not eaf.exists():
        return None

    m = sea.Map(eaf.read_text(errors="ignore"))
    dt_ns = float(m["TrajectoryInterval_ps"].val) / 1000.0        # ns per SID frame

    def _sel(name: str, seltype: str):
        for e in m["Keywords"]:
            k = list(e.keys())[0]
            if k != name:
                continue
            b = e[k]
            try:
                if b["SelectionType"].val == seltype:
                    return b
            except Exception:
                continue
        return None

    def _arr(blk) -> list:
        return [float(x.val) for x in blk["Result"]] if blk is not None else []

    ca = _sel("RMSD", "C-Alpha")
    lg = _sel("RMSD", "Ligand_wrt_protein")
    rf = _sel("RMSF", "C-Alpha")
    d: dict = {"jobname": jobname}
    d["ca"] = _arr(ca)
    d["lg"] = _arr(lg)
    d["rf"] = _arr(rf)
    d["t"] = [i * dt_ns for i in range(len(d["ca"]))]
    # Residue numbers for the RMSF x-axis: 'A:MET_1' → 1; fall back to a 1..N ordinal if unparsable.
    resids = []
    if rf is not None and "ProteinResidues" in rf:
        for x in rf["ProteinResidues"]:
            mm = re.search(r"(\d+)\s*$", str(x.val))
            resids.append(int(mm.group(1)) if mm else None)
    if not resids or any(v is None for v in resids):
        resids = list(range(1, len(d["rf"]) + 1))
    d["res"] = resids

    # Temperature trace (col 9 of the .ene; col 0 = time ps). Drop the leading ramp point.
    te, T = [], []
    if ene.exists():
        for ln in ene.read_text(errors="ignore").splitlines():
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            c = ln.split()
            if len(c) >= 10:
                try:
                    te.append(float(c[0]) / 1000.0); T.append(float(c[9]))
                except ValueError:
                    continue
    d["te"], d["T"] = te[1:], T[1:]
    return d


def make_md_qc_figure(md_dir: Path, out_root: Path, ligands: dict, controls: set = None) -> None:
    """Draw the MD trajectory-QC figure - Cα-RMSD, ligand RMSD, temperature, Cα-RMSF - across the ranks.

    A single stability snapshot proving the production runs are trustworthy before their MM-GBSA / Step-07
    numbers are believed: did the protein equilibrate (Cα-RMSD), did the ligand stay in the pocket (ligand
    RMSD, fit on protein), was the thermostat stable (T), and which regions stayed rigid (Cα-RMSF). Reads
    the SID .eaf + .ene already on disk; best-effort - a rank without SID output is dropped, and the whole
    figure is skipped rather than failing the run. Written to 6_Physics_Validation/06_Analysis/.

    NB: when CFG.MD_RESTRAIN_LIGAND is set the production runs under a positional restraint (ligand heavy
    atoms + backbone), so the Cα-RMSD and ligand-RMSD panels here show RESTRAINT-ENFORCED stability, not
    spontaneous retention - a disclosing footnote is stamped on the figure. The unrestrained reaction
    barrier is the QM/MM ΔE‡ in Step 07, not these positional-restraint panels.
    """
    job_dirs = sorted((d for d in md_dir.iterdir()
                       if d.is_dir() and re.match(r"desmond_md_job_R(?:ank)?_\d", d.name)),
                      key=_natural_rank)
    recs = []
    for jd in job_dirs:
        try:
            r = _read_md_qc(jd, jd.name)
        except Exception as exc:
            _warn(f"[qc] {jd.name} MD read failed - {str(exc).splitlines()[0]}")
            r = None
        if r and r["ca"]:
            _rk = _rank_of(jd.name)
            r["rank"] = int(_rk) if str(_rk).isdigit() else _rk   # int keys the ligand-name map (_lookup_ligands)
            recs.append(r)
    if not recs:
        _warn("[qc] no SID .eaf found - MD trajectory-QC figure skipped.")
        return

    controls = controls or set()
    _rk_list = [r["rank"] for r in recs]
    labs = [_ctrl_label(rk, ligands, controls) for rk in _rk_list]
    _A, _INK = CFG.VIS_ACCENT, CFG.VIS_INK
    palette = [_A["green"], _A["amber"], _A["vermillion"], _A["blue"], _A["magenta"], _A["sky"]]
    cols = _ctrl_palette(_rk_list, controls, palette)
    fl = CFG.VIS_FONT_LEGEND
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.6))

    # A - protein Cα-RMSD vs time (equilibration / drift).
    a = ax[0, 0]
    for i, r in enumerate(recs):
        a.plot(r["t"], r["ca"], color=cols[i], lw=1.3, label=labs[i])
    a.set_xlabel("time (ns)")
    a.set_ylabel("Protein Cα-RMSD (Å)\nvs the minimised start - lower = more stable")

    # B - ligand RMSD after fitting on the protein (did the PFAS stay in the pocket).
    b = ax[0, 1]
    for i, r in enumerate(recs):
        b.plot(r["t"], r["lg"], color=cols[i], lw=1.3, label=labs[i])
    b.set_xlabel("time (ns)")
    _lg_ylab = ("Ligand RMSD (Å), fit on protein\nhigher = drift within the pocket restraint"
                if getattr(CFG, "MD_RESTRAIN_LIGAND", False)
                else "Ligand RMSD (Å), fit on protein\nhigher = drifting out of the pocket")
    b.set_ylabel(_lg_ylab)

    # C - system temperature vs the target (thermostat stability), zoomed to a tight band.
    c = ax[1, 0]
    has_T = False
    for i, r in enumerate(recs):
        if r["te"]:
            c.plot(r["te"], r["T"], color=cols[i], lw=0.7, alpha=0.85, label=labs[i])
            has_T = True
    tgt = CFG.MD_EQUIL_TARGET_T
    c.axhline(tgt, color=_INK["dark"], lw=1.0, ls="--")
    c.text(0.99, 0.53, f"{tgt:g} K target", transform=c.transAxes, ha="right", va="bottom",
           fontsize=CFG.VIS_FONT_ANNOT - 1, color=_INK["soft"])
    if has_T:
        c.set_ylim(tgt - 12, tgt + 10)
    c.set_xlabel("time (ns)")
    c.set_ylabel("System temperature (K)\nthermostat stability around the target")

    # D - per-residue Cα-RMSF (which regions stayed rigid; termini are expectedly mobile).
    dd = ax[1, 1]
    for i, r in enumerate(recs):
        dd.plot(r["res"], r["rf"], color=cols[i], lw=1.0, label=labs[i])
    dd.set_xlabel("residue number")
    dd.set_ylabel("Protein Cα-RMSF (Å)\nper-residue flexibility over the run")

    # one legend for all four panels (identical FA/DFA/TFA/control series), placed inside panel A
    ax[0, 0].legend(loc="lower right", ncol=max(1, len(labs)), fontsize=fl,
                    frameon=True, framealpha=0.9)
    fig.tight_layout()
    qc_dir = _analysis_dir(out_root)
    out_path = qc_dir / "03_MD_Trajectory_QC.svg"
    plt.savefig(out_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    _ok(f"[qc] ✔ MD trajectory QC → {qc_dir.name}/{_utils_mod.deflx_fig_name(out_path.name)}  ({len(recs)} rank(s))")


def _phase_merge(entry: dict, ranked_map: dict, dirs: dict) -> None:
    """Import the prepared complex + write the ESP charges → 01_Prepared_Proteins / 02_ESP_Charged_Complexes."""
    rank, stem, prepared = entry["rank"], entry["stem"], entry["prepared"]
    esp = find_esp(dirs["esp"], stem)
    row = ranked_map.get(stem)
    if row is None:
        raise RuntimeError(f"no ranked-CSV row (job_name) matches handover stem {stem}")
    base = _resnum(row.get("Mapped_Base", ""))
    entry["esp"], entry["base"] = esp, base
    entry["complex_mae"] = dirs["esp_cx"] / f"R_{rank}_{stem}_ESP_Complex.mae"
    shutil.copy2(prepared, dirs["prot"] / prepared.name)
    rep = merge_esp(prepared, esp, entry["complex_mae"], base)
    _base_ok = rep["base_state"] == "HID"
    _base_c = _C.OKGREEN if _base_ok else _C.WARNING
    _echo(f"  {_C.BOLD}{_complex_label(entry, ranked_map)}{_C.ENDC}")
    _echo(f"       {_C.OKGREEN}✔{_C.ENDC} imported     → {dirs['prot'].name}/")
    _echo(f"       {_C.OKGREEN}✔{_C.ENDC} ESP merged   → {dirs['esp_cx'].name}/")
    _echo(f"         ligand    {rep['lig_atoms']} atoms · charge {rep['lig_charge_sum']:+.3f} e")
    _echo(f"         base      His{base} {_base_c}{rep['base_state']}{_C.ENDC}")
    if not _base_ok:
        _warn(f"catalytic base His{base} is {rep['base_state']}, not HID - geometry may be corrupted")


def _phase_watermap(entry: dict, dirs: dict, a) -> None:
    """WaterMap around the ESP complex → 03_WaterMaps (holo). 3 tries, then SKIP (never fails the complex)."""
    rank = entry["rank"]
    wm_dir = dirs["wm"] / f"watermap_R_{rank}"; wm_dir.mkdir(parents=True, exist_ok=True)
    csv_out = dirs["wm"] / f"watermap_R_{rank}.csv"
    if list(wm_dir.glob("*_wm.maegz")) and csv_out.exists():
        _ok(f"[watermap] R_{rank} already done ({csv_out.name}) - skipping"); return
    _log(f"[watermap] R_{rank} {a.wm_ns} ns (holo, {a.lig_dist} Å active site, S-OPLS/TIP4P), "
         f"up to {CFG.PHYS_WM_MAX_TRIES} tries…")
    try:
        with _timed("watermap", rank):
            wmout = run_watermap(entry["complex_mae"], f"watermap_R_{rank}", wm_dir, a.wm_ns, a.lig_dist)
        n = export_watermap_csv(wmout, csv_out)
        if n <= 0:
            # export writes NO csv when the maegz carries no dG site, so claiming ✔ here would print a
            # success for a file that does not exist and the next run's existence-gate would re-run the
            # whole multi-hour WaterMap anyway. Report it as a skip so it is visible and retried.
            _warn(f"[watermap] R_{rank} produced 0 hydration sites - no CSV written; will retry next run.")
        else:
            _ok(f"[watermap] R_{rank} ✔ {wmout.name} · {n} sites → {csv_out.name}")
            # Redraw the combined landscape now this CSV has landed - it reads only the CSVs present,
            # so the cross-rank figure fills in rank by rank as the WaterMaps complete.
            plot_watermap_combined(dirs["wm"], _analysis_dir(dirs["wm"].parent))
    except Exception as exc:
        _fail(f"[watermap] R_{rank} ✘ FAILED after 3 tries - SKIPPED, continuing. "
              f"{(str(exc).splitlines() or ['<no message>'])[0]}")


def _phase_build(entry: dict, dirs: dict) -> None:
    """Desmond System Builder (minimise-volume) + write ESP into the .cms force field → 04_System_Builder."""
    rank = entry["rank"]
    sb_dir = dirs["sb"] / f"desmond_setup_R_{rank}"; sb_dir.mkdir(parents=True, exist_ok=True)
    setup_cms = sb_dir / f"desmond_setup_R_{rank}-out.cms"
    entry["setup_cms"] = setup_cms
    if not entry["complex_mae"].exists():
        raise RuntimeError(f"no ESP complex {entry['complex_mae'].name} - merge first")
    if setup_cms.exists():
        # Resume over an existing build. reapply_esp_to_cms REFUSES a collapsed/corrupt build
        # (the comp_ct invariant), so a resume must not treat that refusal as a permanent failure -
        # reapply is the only path here and run_build lives in the rebuild branch below. Quarantine the
        # bad file and fall through to a fresh build, rather than bricking the rank on every re-run
        # when its on-disk setup_cms is collapsed (comp_ct=1). A genuine ESP name-mismatch will still
        # fail on the fresh build's own reapply, which is correct.
        try:
            q = reapply_esp_to_cms(setup_cms, entry["esp"])
            _ok(f"[build] R_{rank} ✔ already built - ESP re-verified (sum {q:+.3f} e)")
            return
        except RuntimeError as _e:
            _bad = setup_cms.with_suffix(setup_cms.suffix + f".corrupt.{time.strftime('%Y%m%d_%H%M%S')}")
            setup_cms.rename(_bad)
            _log(f"[build] R_{rank} ⚠ existing build unusable - {_e}; moved to {_bad.name}, rebuilding.")
    _log(f"[build] R_{rank} System Builder (minimize-volume, {CFG.PHYS_SOLVENT_MODEL}, "
         f"{CFG.PHYS_FORCEFIELD}, {CFG.PHYS_SALT_CONC_M} M {CFG.PHYS_SALT_POS_ION}{CFG.PHYS_SALT_NEG_ION})…")
    with _timed("build", rank):
        setup_cms = run_build(entry["complex_mae"], f"desmond_setup_R_{rank}", sb_dir)
    entry["setup_cms"] = setup_cms
    q = reapply_esp_to_cms(setup_cms, entry["esp"])
    _ok(f"[build] R_{rank} ✔ {setup_cms.name} · ESP applied to force field (sum {q:+.3f} e)")


def _phase_md(entry: dict, dirs: dict, a) -> "Path | None":
    """Desmond MD (relax + production) → 05_MD_Simulations; unpack _trj/.ene. Returns the job dir for SID/MM-GBSA."""
    rank = entry["rank"]
    setup_cms = entry.get("setup_cms") or (dirs["sb"] / f"desmond_setup_R_{rank}" / f"desmond_setup_R_{rank}-out.cms")
    md_dir = dirs["md"] / f"desmond_md_job_R_{rank}"; md_dir.mkdir(parents=True, exist_ok=True)
    md_cms = md_dir / f"desmond_md_job_R_{rank}-out.cms"
    if md_cms.exists():
        _ok(f"[md] R_{rank} ✔ already done ({md_cms.name}) - skipping"); return md_dir
    if not setup_cms.exists():
        raise RuntimeError(f"no built system {setup_cms.name} - build first")
    reapply_esp_to_cms(setup_cms, entry["esp"])          # re-verify ESP reached the FF before integrating
    _log(f"[md] R_{rank} {a.md_ns} ns production ({a.md_frames} frames, "
         f"NPT {CFG.MD_EQUIL_TARGET_T:g} K, relax + production)…")
    if getattr(CFG, "MD_RESTRAIN_LIGAND", False):
        _log(f"[md] R_{rank} NOTE: positional restraint ON (ligand k={CFG.MD_RESTRAIN_LIG_FORCE_K}, "
             f"backbone k={CFG.MD_RESTRAIN_BB_FORCE_K} kcal/mol/Å²) - pocket retention will be "
             f"restraint-enforced, not spontaneous; unrestrained reactivity is the Step-07 QM/MM barrier.")
    with _timed("md", rank):
        run_md(setup_cms, f"desmond_md_job_R_{rank}", md_dir, a.md_ns, a.md_frames)
    _ok(f"[md] R_{rank} ✔ {md_cms.name} - extraction + SID + MM-GBSA next (CPU)")
    return md_dir


def main() -> int:
    _utils_mod.install_console_rule_filter()   # collapse stacked separator rules
    _install_job_cleanup()                     # Ctrl-C/kill → cancel the run's job-server jobs, no orphans
    t0 = time.perf_counter()
    a = _parse_args_merged()
    stages = a.stage_set

    if not os.path.isfile(SCHROD_RUN):
        print(f"ERROR: Schrödinger 'run' not found at {SCHROD_RUN}."); return 1

    run = (Path(a.run_folder).resolve() if a.run_folder
           else (_SCRIPT_DIR / resolve_run_dir(None)).resolve())
    out_root = Path(a.out) if a.out else run / "6_Physics_Validation"
    dirs = {
        "esp":    run / "5_TopN_and_Preparation" / "4_Ligand_ESP_Charges",
        "prot":   out_root / "01_Prepared_Proteins",
        "esp_cx": out_root / "02_ESP_Charged_Complexes",
        "wm":     out_root / "03_WaterMaps",
        "sb":     out_root / "04_System_Builder",
        "md":     out_root / "05_MD_Simulations",
    }
    for d in (dirs["prot"], dirs["esp_cx"], dirs["wm"], dirs["sb"], dirs["md"]):
        d.mkdir(parents=True, exist_ok=True)

    _open_step_log(out_root)   # 6_Physics_Validation/00_Physics_Validation.log (colour-preserving, fresh)
    print_script_banner("06_Physics_Validation_DeFluorX.py",
                        "ESP Physics - WaterMap → System Builder → MD → SID → MM-GBSA → Defluorination")

    ranked = newest_ranked_csv(run)
    ranked_map = ranked_rows_by_jobname(ranked)
    entries = discover_handover(run)
    _echo(f"  Run Name    : {run.name}")
    _echo(f"  Ranked CSV  : {ranked.name}")
    _echo(f"  Handover    : {len(entries)} complex(es) - ranks {[e['rank'] for e in entries]}")
    _echo(f"  Output root : {out_root}")
    _echo(f"  Settings    : MD {a.md_ns:g} ns/{a.md_frames} fr · WaterMap {a.wm_ns:g} ns · "
          f"site {a.lig_dist:g} Å · stages {sorted(stages)}")

    ok, failed = [], []

    # Hydrate esp/base/complex_mae for every entry up front (pure lookups), so any --stages subset -
    # not just a run that includes 'merge' - has the keys the later phases read. Without this,
    # --stages build/md/watermap KeyError inside the phase and (for watermap) the retry handler swallows
    # it, so main returns 0 = PASS having produced nothing.
    for e in entries:
        try:
            _hydrate(e, ranked_map, dirs)
        except Exception as exc:
            _fail(f"[hydrate] R_{e['rank']} FAILED - {(str(exc).splitlines() or ['<no message>'])[0]}")
            e["_skip"] = True
            failed.append((_complex_label(e, ranked_map),
                           f"hydrate: {(str(exc).splitlines() or ['<no message>'])[0]}"))

    # Prompt for the (optional) sudo password NOW, up front, so the run is fully unattended
    # afterwards - the user can walk away and SID / MM-GBSA stay protected from systemd-oomd.
    global _OOMD_GUARD
    _guard = OomdGuard(active=not a.pipeline_mode)
    _guard.__enter__()
    _OOMD_GUARD = _guard        # reachable from the signal/atexit handlers so a kill still unmasks oomd

    if "merge" in stages:
        _section(f"Step 1/4 - Import + ESP merge  ({len(entries)} complex)")
        for _i, e in enumerate(entries):
            if _i:
                _echo("")
            try:
                _phase_merge(e, ranked_map, dirs)
            except Exception as exc:
                _fail(f"[merge] R_{e['rank']} FAILED - {str(exc).splitlines()[0]}")
                e["_skip"] = True; failed.append((_complex_label(e, ranked_map), f"merge: {str(exc).splitlines()[0]}"))

    if "watermap" in stages:
        _section(f"Step 2/4 - WaterMap  ({a.wm_ns:g} ns)")
        _first = True
        for e in entries:
            if e.get("_skip"):
                continue
            if not _first:
                _echo("")
            _first = False
            _phase_watermap(e, dirs, a)

    if "build" in stages:
        _section("Step 3/4 - System Builder  (minimise-volume)")
        _first = True
        for e in entries:
            if e.get("_skip"):
                continue
            if not _first:
                _echo("")
            _first = False
            try:
                _phase_build(e, dirs)
            except Exception as exc:
                _fail(f"[build] R_{e['rank']} FAILED - {str(exc).splitlines()[0]}")
                e["_skip"] = True; failed.append((_complex_label(e, ranked_map), f"build: {str(exc).splitlines()[0]}"))

    # Build + WaterMap are both done for every rank by here → draw the physics build/solvation QC figure
    # (ESP charge gradient, per-atom charges, WaterMap ΔG, box size). Best-effort; never fails the run.
    if {"build", "watermap"} & stages:
        _section("Physics QC - build & solvation snapshot")
        try:
            make_physics_qc_figure(entries, dirs, out_root, _lookup_ligands(run), _lookup_controls(run))
        except Exception as exc:
            _warn(f"[qc] physics QC figure skipped - {str(exc).splitlines()[0]}")

    _mmgbsa_status = "ok"
    if "md" in stages:
        _section(f"Step 4/4 - MD → SID → MM-GBSA → Defluorination  ({a.md_ns:g} ns · MD on GPU, post-processing pipelined on CPU)")
        run_root = run
        md_dir = dirs["md"]
        with _guard:                                     # reuse the guard primed at start (idempotent)
            # Pipelined per rank (GPU produces, CPU consumes). Each rank's MD runs on the GPU (blocking,
            # sequential - there is one GPU); the instant it lands, that rank's extract → SID → MM-GBSA →
            # Defluorination is put on a queue and the GPU immediately starts the NEXT rank's MD. A single
            # CPU worker drains the queue one rank at a time: if a later MD finishes while an earlier
            # rank's post-processing is still running, the new rank waits in the queue. Exactly ONE
            # Prime MM-GBSA batch touches the scratch disk at a time, while the GPU runs the next rank's
            # MD through the hours of CPU post-processing instead of sitting idle.
            # NB: MD (main thread) and SID/MM-GBSA (worker thread) log concurrently - lines interleave.
            _active_entries = [e for e in entries if not e.get("_skip")]
            _job_total = len(_active_entries)
            _post_q: "queue.Queue" = queue.Queue()
            _post_state = {"warn": False}

            def _post_worker():
                """CPU consumer: Extraction → SID → MM-GBSA → Defluorination, one queued rank at a time (the waiting list)."""
                while True:
                    _item = _post_q.get()
                    if _item is None:
                        _post_q.task_done(); break
                    _jd, _rank, _e, _jp = _item
                    try:
                        with _timed("extract", _rank):                             # unpack production {job}_trj/ + {job}.ene (CPU)
                            _unpack_md_production(_jd, _jd.name, a.md_ns * 1000.0)
                        with _timed("sid", _rank):
                            _sid_todo = scan_jobs([_jd])
                            if not _sid_todo:
                                _echo(f"  {_C.CYAN}▸{_C.ENDC} [Job {_jp}/{_job_total}] Rank {_rank}: SID already done - advancing to MM-GBSA.")
                            process_jobs(_sid_todo, _jp, _job_total)          # SID  (blocking, CPU)
                        with _timed("mmgbsa", _rank):
                            _mg = run_mmgbsa(_jd, _jd.name, _rank_of(_jd.name))   # MM-GBSA (blocking, CPU)
                        if _mg is None:
                            _post_state["warn"] = True
                        else:
                            # This rank is complete - draw its own MM-GBSA profile now rather than
                            # waiting for the cross-rank Finalise pass.
                            try:
                                _plot_rank_mmgbsa(md_dir, _jd, _jd.name, _rank, _mg, quiet=True)
                            except Exception as _exc:
                                _warn(f"[mmgbsa] R_{_rank} per-rank figure skipped - "
                                      f"{str(_exc).splitlines()[0]}")
                        with _timed("defluor", _rank):                            # defluorination geometry
                            run_defluorination(_jd, _jd.name, _rank, md_dir, run_root, a.md_ns)
                    except Exception as _exc:
                        _warn(f"[sid/mmgbsa/defluor] R_{_rank}: {str(_exc).splitlines()[0]}")
                        _post_state["warn"] = True
                    finally:
                        _wait_settled(_jd)           # let this rank's post-processing outputs settle before the next
                        _post_q.task_done()

            _worker = threading.Thread(target=_post_worker, name="md-post-worker", daemon=True)
            _worker.start()

            _job_pos = 0
            for e in entries:
                if e.get("_skip"):
                    continue
                _job_pos += 1
                if _job_pos > 1:
                    _echo("")
                rank = e["rank"]
                try:
                    jd = _phase_md(e, dirs, a)                          # MD  (GPU, blocking)
                except Exception as exc:
                    _fail(f"[md] R_{rank} FAILED - {str(exc).splitlines()[0]}")
                    failed.append((_complex_label(e, ranked_map), f"md: {str(exc).splitlines()[0]}"))
                    continue
                if jd is None:
                    continue
                # Let the just-landed MD outputs settle before post-processing reads them and before the
                # next rank's MD starts, so an in-flight multisim write never collides with a fresh run.
                _log(f"[md] R_{rank} landed - waiting until job files stop changing before post-processing / next MD…")
                if not _wait_settled(jd):
                    _warn(f"[md] R_{rank} settle timed out after {_SETTLE_TIMEOUT_SEC}s - proceeding anyway")
                ok.append(_complex_label(e, ranked_map))               # MD landed; SID/MM-GBSA issues are WARN
                _post_q.put((jd, rank, e, _job_pos))                   # queue CPU post-processing; GPU moves to next rank

            _post_q.put(None)                                          # end of MD stream - drain the waiting list
            _worker.join()                                            # wait for all queued SID → MM-GBSA → Defluorination
            if _post_state["warn"]:
                _mmgbsa_status = "warn"

            # Finalise: read the per-rank MM-GBSA CSVs and draw the combined cross-rank plots (no re-run).
            _section("Finalise - MM-GBSA combined plots")
            job_dirs = sorted((d for d in md_dir.iterdir()
                               if d.is_dir() and re.match(r"desmond_md_job_R(?:ank)?_\d", d.name)),
                              key=_natural_rank)
            if job_dirs:
                _st = run_mmgbsa_phase(md_dir, run_root, plots_only=True)   # read existing CSVs → plots only, never re-run
                if _st != "ok":
                    _mmgbsa_status = _st
                try:
                    plot_defluor_combined(md_dir, _lookup_ligands(run_root), _lookup_controls(run_root))
                except Exception as exc:
                    _warn(f"[defluor] combined figure skipped - {str(exc).splitlines()[0]}")
                # MD trajectory QC (Cα-RMSD · ligand RMSD · temperature · Cα-RMSF) from the SID .eaf + .ene.
                try:
                    make_md_qc_figure(md_dir, out_root, _lookup_ligands(run_root), _lookup_controls(run_root))
                except Exception as exc:
                    _warn(f"[qc] MD trajectory-QC figure skipped - {str(exc).splitlines()[0]}")

    # Combined WaterMap landscape - drawn every run from the WaterMaps already on disk (never gated on a
    # fresh export, so a resume where WaterMap is skipped still refreshes it), like the other cross-rank
    # figures in 06_Analysis. Reads only the CSVs present and no-ops when none exist.
    try:
        plot_watermap_combined(dirs["wm"], _analysis_dir(out_root))
    except Exception as exc:
        _warn(f"[watermap] combined landscape figure skipped - {str(exc).splitlines()[0]}")

    _section(f"Summary - {len(ok)} ok, {len(failed)} failed  (stages {sorted(stages)})")
    for t in ok:
        _echo(f"  {_C.OKGREEN}✔{_C.ENDC} {t}")
    for t, why in failed:
        _echo(f"  {_C.FAIL}✘{_C.ENDC} {t}  →  {why}")
    _emit_timings(out_root)                               # per-job + per-phase wall-clock → log + 00_Phase_Timings.csv
    _guard.__exit__(None, None, None)                    # restore systemd-oomd (idempotent if Step 4 already did)
    print_elapsed(t0, "06_Physics_Validation_DeFluorX.py")
    return EXIT_WARN if (failed or _mmgbsa_status == "warn") else 0


if __name__ == "__main__":
    sys.exit(main())
