#!/usr/bin/env python3
"""
===============================================================================
FAcDs Pipeline  |  Step 06  |  Desmond SID + Prime MM-GBSA Post-Processing
===============================================================================
Runs Schrödinger Event Analysis (event_analysis.py) and Simulation Interaction
Diagram analysis (analyze_simulation.py) on completed Desmond molecular
dynamics trajectories, producing the *_SID-out.eaf files that are consumed
downstream by 07_MD_Thermodynamics_QMMM_Engine_FAcDs.py.

It then runs Prime MM-GBSA (thermal_mmgbsa.py <job>-out.cms) on every completed
MD job — the end-state ligand binding free energy over the MD ensemble — and
plots per-job + combined ΔG_bind. MM-GBSA is complementary to the QSite QM/MM
reaction barrier (Step 07): it scores BINDING, not C–F bond cleavage.

Uses the central CFG / ProjectUtils modules for logging, console styling, and
conventions shared across the pipeline. Runs under the project 'PFAS' conda
environment and shells out to $SCHRODINGER/run for the Schrödinger interpreter
— it does NOT need to be launched with $SCHRODINGER/run.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 05 July 2026
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
                  .../desmond_md_job_R_N/*-prime-mmgbsa.csv  (per-frame ΔG_bind)
                  .../MolecularDynamics/Prime_MMGBSA/00_MMGBSA_Summary.csv
                  .../MolecularDynamics/Prime_MMGBSA/01_MMGBSA_Combined_AllRanks.png
                  .../MolecularDynamics/Prime_MMGBSA/Rank_NN_MMGBSA_Profile_*.png
  Upstream      : Desmond molecular dynamics simulations (manual Maestro step).
  Downstream    : 07_MD_Thermodynamics_QMMM_Engine_FAcDs.py (Step 07; consumes EAF output).
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
  5. MM-GBSA Cost & Validity: thermal_mmgbsa.py runs Prime per frame. On a
     100k-frame trajectory this is intractable — set CFG.MMGBSA_STEP_SIZE to
     subsample. GB implicit solvent overstabilises anionic PFAS, so ΔG_bind is
     a RELATIVE ranking only; it scores binding, not the QSite reaction barrier.
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

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

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

# Auto-set SCHRODINGER if the env var is absent (mirrors the 08 engine default).
os.environ.setdefault("SCHRODINGER", "/opt/schrodinger")
SCHRODINGER = os.environ["SCHRODINGER"]
SCHROD_RUN = os.path.join(SCHRODINGER, "run")

_SEP = "============================================================================="
_DEFAULT_FRAME_TOTAL = 100_000   # heartbeat fallback when the trajectory length is unreadable
EXIT_WARN = 3   # step completed but a complementary part (MM-GBSA) was deferred/failed; the pipeline runner renders WARN and continues (0=PASS, 1=hard error, 3=warn)

# Tokens the SID-out.eaf Result vector carries per trajectory frame. Governs the
# completion gate (is_eaf_complete): threshold = EAF_TOKENS_PER_FRAME · traj_frames.
# 1 is the standard per-frame scalar series; set to the measured ratio once read
# off a known-good, fully-analysed EAF so a partial high-k EAF cannot pass.
EAF_TOKENS_PER_FRAME = getattr(CFG, "EAF_TOKENS_PER_FRAME", 1)


# =============================================================================
# SECTION 2: SMALL HELPERS
# =============================================================================
def _echo(msg: str = "") -> None:
    """Print to terminal immediately (flush) — keeps live progress visible."""
    print(msg, flush=True)


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
        _echo("Optional: sudo masks systemd-oomd to prevent Out-Of-Memory kills during long SID runs.")
        _echo("Enter your sudo password to enable masking, or press Enter/Ctrl-D to skip (NORMAL mode).")
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

        _echo("Masking systemd-oomd...")
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
# SECTION 4: HEARTBEAT (progress of a long-running analyze_simulation step)
# =============================================================================
class Heartbeat:
    """Periodically report a long Schrödinger step's progress by tailing its log.

    Emits newline-terminated "[PROGRESS] frame X/total" lines (NOT an in-place
    \\r ticker) so the progress is captured in the redirected pipeline log.
    Each step prints its own per-frame marker; pass the matching regex(es):
      • analyze_simulation.py (SID)  → "analyzing frame# N"   (default)
      • thermal_mmgbsa.py (MM-GBSA)  → "Reading frame N" / "Structure N"
    The latest captured N is shown as a percentage of the known trajectory length.
    """

    def __init__(self, log_file: Path, label: str, total: int,
                 patterns: list[str] | None = None, interval: int = 120):
        self.log_file = Path(log_file)
        self.label = label
        self.total = max(int(total), 1)
        self.patterns = [re.compile(p) for p in (patterns or [r"analyzing frame# (\d+)"])]
        self.interval = max(5, int(interval))
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._start = 0.0
        self._prime_announced = False   # print the Prime-phase transition banner once

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

    def _run(self):
        while not self._stop.wait(self.interval):
            elapsed_min = int((time.time() - self._start) // 60)
            try:
                text = self.log_file.read_text(errors="ignore")
            except Exception:
                text = ""
            # Prime MM-GBSA scoring phase: trajectory read is complete and the main
            # log is silent while Prime minimises every frame in a separate job.
            if self._PRIME_HANDOFF.search(text):
                if not self._prime_announced:
                    self._prime_announced = True
                    _echo(f"    [PHASE 2/2] {self.label}: trajectory read complete — "
                          f"Prime MM-GBSA minimisation started on {self.total:,} frames "
                          f"(the long phase; Prime writes no per-frame log, so progress "
                          f"below is by subjob or elapsed time).")
                prog = self._prime_progress()
                if prog is not None:
                    _done, _tot = prog
                    _pct = min(100, _done * 100 // max(_tot, 1))
                    _echo(f"    [PROGRESS] {self.label}: Prime scoring "
                          f"{_done}/{_tot} subjobs done ({_pct}%, {elapsed_min}m elapsed)")
                else:
                    _echo(f"    [PROGRESS] {self.label}: Prime minimising "
                          f"{self.total:,} frames — working, {elapsed_min}m elapsed")
                continue
            frame = self._latest_frame(text)
            pct = min(100, frame * 100 // self.total)
            _echo(f"    [PROGRESS] {self.label}: frame {frame:,}/{self.total:,} "
                  f"({pct}%, {elapsed_min}m elapsed)")

    def __enter__(self):
        self._start = time.time()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1)
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
    _echo("Scanning job directories...")
    to_run: list[Path] = []
    completed = running = pending = 0

    for d in job_dirs:
        job_name = d.name
        rank = _rank_of(job_name)
        out_eaf = d / f"{job_name}_SID-out.eaf"
        cms_file = d / f"{job_name}-out.cms"
        trj_dir = d / f"{job_name}_trj"

        """
        1. No -out.cms => MD simulation not finished. A live process here means
        the SIMULATION is running (job-name match avoids SID processes).
        """
        if not cms_file.is_file():
            if _proc_alive(job_name):
                running += 1
                _echo(f"  - Rank {rank}: RUNNING (simulation active)")
            else:
                pending += 1
                _echo(f"  - Rank {rank}: PENDING (simulation not finished yet)")
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
                _echo(f"  - Rank {rank}: COMPLETED (SID-out.eaf has {of}/{tf} frames)")
            else:
                to_run.append(d)
                _echo(f"  - Rank {rank}: INCOMPLETE (SID-out.eaf has {of}/{tf} frames; will re-run)")
        else:
            to_run.append(d)
            _echo(f"  - Rank {rank}: READY (simulation finished, SID pending)")

    _echo("")
    _echo("Scan Summary:")
    _echo(f"  Total Jobs Found : {len(job_dirs)}")
    _echo(f"  Completed Jobs   : {completed}")
    _echo(f"  Running Jobs     : {running}")
    _echo(f"  Pending Jobs     : {pending}")
    _echo(f"  Jobs to Process  : {len(to_run)}")
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
# =============================================================================
"""
Mirrors the manual workflow:
    cd <desmond_md_job_R_N> ; $SCHRODINGER/run thermal_mmgbsa.py <job>-out.cms
run folder-wise and idempotent (skip when a valid results CSV already exists),
then plots per-job + combined ΔG_bind. MM-GBSA is complementary to the QSite
QM/MM reaction barrier (Step 07), not a replacement: it scores binding, not
bond cleavage. GB implicit solvent overstabilises anionic PFAS → relative only.
"""

def _mmgbsa_csv(job_dir: Path, job_name: str):
    """Locate a thermal_mmgbsa results CSV in `job_dir` (name varies by version).
    Returns the newest match or None."""
    hits = []
    for pat in (f"{job_name}-out-prime-mmgbsa.csv", f"{job_name}*prime*mmgbsa*.csv",
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


def run_mmgbsa(job_dir: Path, job_name: str, rank: str) -> Path | None:
    """Run thermal_mmgbsa.py on <job>-out.cms inside its MD folder. Idempotent:
    returns the existing CSV when already computed. Returns the results CSV path
    or None on failure."""
    cms_file = job_dir / f"{job_name}-out.cms"
    if not cms_file.is_file():
        _echo(f"  [Rank {rank}] MM-GBSA skipped — no {cms_file.name}.")
        return None
    existing = _mmgbsa_csv(job_dir, job_name)
    if existing is not None:
        _echo(f"  [Rank {rank}] MM-GBSA already done — reusing {existing.name}.")
        return existing

    """
    Parallelise frame subjobs across cores: total cores − reserve (same cap as
    the rest of the pipeline, CFG.GLOBAL_MAX_WORKERS = cpu_count − PREP_CPU_RESERVE).
    """
    _ncpu = max(1, int(getattr(CFG, "GLOBAL_MAX_WORKERS", max(1, (os.cpu_count() or 4) - 2))))
    cmd = [SCHROD_RUN, "thermal_mmgbsa.py", cms_file.name,
           "-j", f"{job_name}_mmgbsa", "-HOST", f"localhost:{_ncpu}"]
    if getattr(CFG, "MMGBSA_STEP_SIZE", 0) and CFG.MMGBSA_STEP_SIZE > 0:
        cmd += ["-step_size", str(CFG.MMGBSA_STEP_SIZE)]
    _echo(f"  [Rank {rank}] Running MM-GBSA on {_ncpu} cores: {' '.join(cmd[1:])}")
    _echo(f"  [Rank {rank}] [PHASE 1/2] reading trajectory frames "
          f"(thermal_mmgbsa), then hands off to a {_ncpu}-subjob Prime minimisation.")
    log = job_dir / f"{job_name}_mmgbsa.log"
    _timeout = getattr(CFG, "MMGBSA_TIMEOUT_SEC", 0) or None
    _interval = max(5, int(getattr(CFG, "MMGBSA_PROGRESS_INTERVAL_SEC", 30)))
    # Ground-truth denominator (same source the SID phase uses).
    _total = traj_frame_count(job_dir / f"{job_name}_trj") or _DEFAULT_FRAME_TOTAL
    """
    thermal_mmgbsa logs "Reading frame N..." during the (serial) trajectory
    extraction phase and "Structure N" during the parallel Prime phase. The
    Heartbeat tails the log and emits newline "[PROGRESS] frame X/total" lines —
    captured in the redirected pipeline log (an in-place \\r ticker is not).
    """
    _hb_patterns = [r"Reading frame (\d+)", r"[Ss]tructure[:\s]+(\d+)",
                    r"[Ff]rame\s+(\d+)", r"Processing\s+(\d+)"]
    try:
        with open(log, "w") as fh:
            proc = subprocess.Popen(cmd, cwd=str(job_dir), stdout=fh,
                                    stderr=subprocess.STDOUT)
            with Heartbeat(log, f"MM-GBSA Rank {rank}", _total,
                           patterns=_hb_patterns, interval=_interval):
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
    return _mmgbsa_csv(job_dir, job_name)


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
    return pd.to_numeric(df[col], errors="coerce").dropna()


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


def plot_mmgbsa_individual(out_dir: Path, job_name: str, rank: str, dg: "pd.Series") -> None:
    """Per-job MM-GBSA: ΔG_bind time series + distribution."""
    if dg.empty:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [2, 1]})
    ax1.plot(range(len(dg)), dg.values, color="#0072B2", linewidth=1.0, alpha=0.5)
    ax1.plot(range(len(dg)), dg.rolling(max(1, len(dg) // 20), min_periods=1).mean().values,
             color="#222222", linewidth=2.2)
    # Central line uses the configured headline estimator (CFG.MMGBSA_AVERAGING),
    # so this figure, the violin and the combined bar all report the same value.
    _central = _avg_dg(dg)
    _est_label = str(getattr(CFG, "MMGBSA_AVERAGING", "mean")).lower()
    ax1.axhline(_central, color="#D55E00", linestyle="--", linewidth=1.5,
                label=f"{_est_label} = {_central:.1f} kcal/mol")
    ax1.set_xlabel("MM-GBSA frame", fontweight="bold")
    ax1.set_ylabel("ΔG_bind (kcal/mol)", fontweight="bold")
    ax1.set_title(f"MM-GBSA binding free energy — Rank {rank}", fontsize=12, fontweight="bold")
    ax1.legend(frameon=True, fontsize=9)
    ax2.hist(dg.values, bins=25, color="#1B9E77", alpha=0.85, orientation="horizontal")
    ax2.axhline(_central, color="#D55E00", linestyle="--", linewidth=1.5)
    ax2.set_xlabel("Frames", fontweight="bold")
    ax2.set_title("Distribution", fontsize=11, fontweight="bold")
    _rk = f"{int(rank):02d}" if str(rank).isdigit() else str(rank)
    out_path = out_dir / f"Rank_{_rk}_MMGBSA_Profile_{job_name}.png"
    plt.savefig(out_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    _echo(f"  [Rank {rank}] MM-GBSA profile saved: {out_path.name}")


def plot_mmgbsa_combined(out_dir: Path, per_job: list, tiers: dict) -> None:
    """Combined MM-GBSA across ranks: ΔG_bind distribution (violin) + mean bar."""
    per_job = [(r, dg) for r, dg in per_job if not dg.empty]
    if not per_job:
        _echo("  [!] MM-GBSA combined skipped — no parsed ΔG_bind series.")
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    per_job.sort(key=lambda x: int(x[0]))
    labels = [f"#{r}" for r, _ in per_job]
    data = [dg.values for _, dg in per_job]
    means = [_avg_dg(dg) for _, dg in per_job]
    _mode = str(getattr(CFG, "MMGBSA_AVERAGING", "mean")).lower()
    _avg_label = {"mean": "Arithmetic-mean", "median": "Median"}.get(_mode, "Boltzmann-weighted")
    cols = [CFG.TIER_COLOUR.get(tiers.get(int(r), ""), "#888888") for r, _ in per_job]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 6.5), gridspec_kw={"wspace": 0.18})
    # showmeans draws the arithmetic mean tick; only correct when that IS the configured
    # estimator. For median/boltzmann, overlay the configured central value instead.
    parts = ax1.violinplot(data, showmeans=(_mode == "mean"), showextrema=False)
    if _mode != "mean":
        ax1.scatter(range(1, len(means) + 1), means, marker="_", s=400,
                    color="#111111", zorder=5, label=f"{_avg_label} ΔG")
    for i, b in enumerate(parts["bodies"]):
        b.set_facecolor(cols[i]); b.set_alpha(0.75)
    ax1.set_xticks(range(1, len(labels) + 1)); ax1.set_xticklabels(labels)
    ax1.set_ylabel("ΔG_bind (kcal/mol)", fontweight="bold")
    ax1.set_title("MM-GBSA binding free-energy distribution per rank", fontsize=12, fontweight="bold")
    ax2.bar(labels, means, color=cols, alpha=0.88, edgecolor="#222222", linewidth=0.6)
    for i, v in enumerate(means):
        ax2.text(i, v, f"{v:.1f}", ha="center",
                 va="bottom" if v >= 0 else "top", fontsize=9, fontweight="bold")
    ax2.axhline(0.0, color="#222222", linewidth=0.8)
    ax2.set_ylabel(f"{_avg_label} ⟨ΔG_bind⟩ (kcal/mol)", fontweight="bold")
    ax2.set_title(f"{_avg_label} MM-GBSA ΔG_bind — more negative = tighter binding", fontsize=12, fontweight="bold")
    out_path = out_dir / "01_MMGBSA_Combined_AllRanks.png"
    plt.savefig(out_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    _echo(f"  MM-GBSA combined figure saved: {out_path.resolve()}")


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
    _echo(f"Prime MM-GBSA — {len(job_dirs)} completed MD job(s)")
    _echo(_SEP)
    # Figures folder derived from the resolved MD dir (not a hardcoded path).
    out_dir = md_dir / getattr(CFG, "MMGBSA_OUTPUT_SUBDIR", "Prime_MMGBSA")
    tiers = _lookup_tiers(run_root)
    per_job = []
    summary_rows = []
    _failed = 0
    for d in job_dirs:
        job_name = d.name
        rank = _rank_of(job_name)
        csv = run_mmgbsa(d, job_name, rank)
        if csv is None:
            _failed += 1
            continue
        dg = _mmgbsa_dg_series(csv)
        try:
            plot_mmgbsa_individual(out_dir, job_name, rank, dg)
        except Exception as e:
            _echo(f"  [Rank {rank}] MM-GBSA per-job plot failed ({e}) — skipped.")
        per_job.append((rank, dg))
        if not dg.empty:
            summary_rows.append({
                "Scientific_Rank": rank, "Job_Name": job_name,
                "degrader_tier": tiers.get(int(rank), "") if str(rank).isdigit() else "",
                "N_Frames": int(dg.size),
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
        plot_mmgbsa_combined(out_dir, per_job, tiers)
    except Exception as e:
        _echo(f"  [!] MM-GBSA combined plot failed ({e}) — skipped.")
    if _failed:
        _echo(f"  [!] MM-GBSA: {_failed} of {len(job_dirs)} job(s) failed "
              f"(see per-rank rc/cause above) — step will report WARN, not PASS.")
        return "warn"
    return "ok"


# =============================================================================
# SECTION 9: MAIN
# =============================================================================
def main() -> int:
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
         (no stop, no killed job), so MM-GBSA is no longer deferred when another job
         happens to be running. This is what actually prevents the /tmp overflow.
      2. Belt-and-braces — point SCHRODINGER_TMPDIR/TMPDIR at the working disk for
         any tool that still honours them (driver-side temp, non-server steps).
    """
    _scratch_ok = _utils_mod.ensure_jobserver_on_working_disk(
        SCHRODINGER, _SCRIPT_DIR,
        getattr(CFG, "SCHRODINGER_JOBSERVER_SUBDIR", "_Schrodinger_JobServer"),
        echo=lambda m="": _echo(f"  {m}"))
    _scratch = md_dir / getattr(CFG, "SCHRODINGER_SCRATCH_SUBDIR", "_Schrodinger_Scratch")
    _scratch.mkdir(parents=True, exist_ok=True)
    os.environ["SCHRODINGER_TMPDIR"] = str(_scratch)
    os.environ["TMPDIR"] = str(_scratch)

    _echo(_SEP)
    _echo("Starting Desmond SID Analysis Post-processing")
    _echo(f"Run Directory      : {run_dir}")
    _echo(f"Molecular Dynamics : {md_dir}")
    _echo(f"Schrödinger scratch : {_scratch}  (on the run's working disk, not /tmp)")
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
            _echo("No jobs ready for SID analysis (all complete or pending).")
        _mmgbsa_status = run_mmgbsa_phase(md_dir, run_root, scratch_ok=_scratch_ok)

    _echo("")
    # Remove Schrödinger scratch now that all jobs finished and outputs are in the
    # working folders (mirrors the temp-thumbnail cleanup in 03). Guarded on the dir
    # name so only the dedicated scratch folder can ever be removed.
    if _scratch.name == getattr(CFG, "SCHRODINGER_SCRATCH_SUBDIR", "_Schrodinger_Scratch") and _scratch.exists():
        time.sleep(getattr(CFG, "SCHRODINGER_SCRATCH_COOLDOWN_SEC", 5))  # let outputs settle first
        shutil.rmtree(_scratch, ignore_errors=True)
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
