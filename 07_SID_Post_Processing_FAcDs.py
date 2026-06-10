#!/usr/bin/env python3
# =============================================================================
# FAcDs Pipeline  |  Step 07  |  Desmond SID Post-Processing
# =============================================================================
# Runs Schrödinger Event Analysis (event_analysis.py) and Simulation Interaction
# Diagram analysis (analyze_simulation.py) on completed Desmond molecular
# dynamics trajectories, producing the *_SID-out.eaf files that are consumed
# downstream by 08_MD_Thermodynamics_QMMM_Engine_FAcDs.py.
#
# This is the Python port of the former run_sid_analysis_FAcDs.sh /
# 07_SID_Post_Processing_FAcDs.sh shell orchestrator. The logic is unchanged;
# the rewrite adopts the central CFG / ProjectUtils modules so Step 07 shares
# the same logging, console styling, and conventions as the rest of the
# pipeline. The script itself runs under the project 'PFAS' conda environment
# and shells out to $SCHRODINGER/run for the Schrödinger interpreter — it does
# NOT need to be launched with $SCHRODINGER/run.
#
# Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
# Date   : 10 June 2026
# =============================================================================
# Usage (standalone):
#   python 07_SID_Post_Processing_FAcDs.py [Boltz-2_Run_Directory]
#
# Usage (from pipeline runner — oomd already managed by 00_00_run_pipeline):
#   python 07_SID_Post_Processing_FAcDs.py [Boltz-2_Run_Directory] --pipeline-mode
#
# ── Dependency Map ────────────────────────────────────────────────────────────
#   Script        : 07_SID_Post_Processing_FAcDs.py
#   Role          : Step 07 — Desmond post-simulation post-processing.
#                   Produces EAF interaction files for the Step 08 MD engine.
#   Imports from  : 00_02_Project_Config_FAcDs.py  (CFG — project metadata)
#                   00_03_Project_Utils_FAcDs.py   (ConsoleColours, logging, banners)
#   Reads         : Boltz-2_Run_X/7_Physics_Validation/MolecularDynamics/desmond_md_job_Rank_N/*-out.cms
#                   Boltz-2_Run_X/7_Physics_Validation/MolecularDynamics/desmond_md_job_Rank_N/*_trj
#   Writes        : Boltz-2_Run_X/7_Physics_Validation/MolecularDynamics/desmond_md_job_Rank_N/*_SID-in.eaf
#                   Boltz-2_Run_X/7_Physics_Validation/MolecularDynamics/desmond_md_job_Rank_N/*_SID-out.eaf
#                   Boltz-2_Run_X/7_Physics_Validation/MolecularDynamics/desmond_md_job_Rank_N/*.log
#   Upstream      : Desmond molecular dynamics simulations (manual Maestro step).
#   Downstream    : 08_MD_Thermodynamics_QMMM_Engine_FAcDs.py (Step 08; consumes EAF output).
# ─────────────────────────────────────────────────────────────────────────────
# ── The Critic's Corner: Known Limitations & Failure Points ───────────────────
#   1. Sudo Requirement: Masking systemd-oomd requires root privileges. When
#      run standalone the script prompts for sudo; when run from the pipeline
#      runner (--pipeline-mode) the credential is already cached and oomd
#      masking is owned by the runner across Steps 07 and 08.
#   2. Trajectory Volume: Processing 100,000 frames sequentially is time-intensive
#      (typically 10–15 hours). Do not run multiple instances concurrently.
#   3. Local I/O Capping: Relies on -LOCAL to prevent huge tmp partition writes;
#      requires sufficient disk space in the destination filesystem.
#   4. No Desmond Jobs: If no completed Desmond jobs are found, the script exits
#      with status 0 (non-fatal). The pipeline runner continues to Step 08,
#      which will then report that no EAF files are present.
# =============================================================================

import argparse
import importlib.util as _ilu
import os
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

# -------------------------------------------------------------------------------
# Step 1: Pipeline modules (00_02 CFG, 00_03 utils) via importlib
# -------------------------------------------------------------------------------
def _load_module(name: str, path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Required module not found: {path}")
    spec = _ilu.spec_from_file_location(name, str(path))
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_SCRIPT_DIR = Path(__file__).resolve().parent
_utils_mod = _load_module("ProjectUtils", _SCRIPT_DIR / "00_03_Project_Utils_FAcDs.py")
_cfg_mod = _load_module("ProjectConfig", _SCRIPT_DIR / "00_02_Project_Config_FAcDs.py")

CFG = _cfg_mod.CFG()
print_script_banner = _utils_mod.print_script_banner
print_elapsed = _utils_mod.print_elapsed

# Auto-set SCHRODINGER if the env var is absent (mirrors the 08 engine default).
os.environ.setdefault("SCHRODINGER", "/opt/schrodinger")
SCHRODINGER = os.environ["SCHRODINGER"]
SCHROD_RUN = os.path.join(SCHRODINGER, "run")

_SEP = "============================================================================="
_DEFAULT_FRAME_TOTAL = 100_000   # heartbeat fallback when the trajectory length is unreadable


# ===============================================================================
# SECTION 2: SMALL HELPERS
# ===============================================================================
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


def is_eaf_complete(out_frames: int, traj_frames: int) -> bool:
    """A SID-out.eaf is complete when it holds at least as many frames as the
    trajectory. If the trajectory count cannot be read (traj_frames == 0), an
    existing non-empty EAF is treated as complete — never destructively re-run
    on doubt.
    """
    if traj_frames > 0:
        return out_frames >= traj_frames
    return out_frames > 0


def _natural_rank(job_dir: Path) -> int:
    """Sort key for desmond_md_job_Rank_N directories (numeric, -V style)."""
    m = re.search(r"_Rank_(\d+)", job_dir.name)
    return int(m.group(1)) if m else 0


def _rank_of(job_name: str) -> str:
    """Extract the Rank token after '_md_job_Rank_' (matches bash ${JOBNAME##*_md_job_Rank_})."""
    return job_name.split("_md_job_Rank_")[-1]


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


# ===============================================================================
# SECTION 3: SYSTEMD-OOMD MANAGEMENT (standalone mode only)
# ===============================================================================
class OomdGuard:
    """Mask systemd-oomd to prevent Out-Of-Memory kills during long SID runs.

    Standalone mode only. When called via the pipeline runner the runner owns
    oomd masking across both Step 07 (this script) and Step 08 (MD engine), so
    this guard is a no-op (`active=False`) to avoid premature unmask between the
    two steps. A background thread keeps the sudo credential warm.
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
        _echo("This script requires sudo to mask systemd-oomd to prevent Out-Of-Memory kills.")
        # Prime the sudo credential interactively against the controlling tty.
        try:
            with open("/dev/tty") as _tty:
                subprocess.run(["sudo", "-v"], stdin=_tty, check=True)
        except Exception:
            subprocess.run(["sudo", "-v"], check=True)

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


# ===============================================================================
# SECTION 4: HEARTBEAT (progress of a long-running analyze_simulation step)
# ===============================================================================
class Heartbeat:
    """Periodically report the analyse-simulation progress by tailing its log.

    analyze_simulation.py prints "analyzing frame# N..." per frame; the
    heartbeat reads the latest N and converts it to a percentage against the
    known trajectory length.
    """

    def __init__(self, log_file: Path, label: str, total: int):
        self.log_file = Path(log_file)
        self.label = label
        self.total = max(int(total), 1)
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._start = 0.0

    def _frame_pat(self) -> re.Pattern:
        return re.compile(r"analyzing frame# (\d+)")

    def _run(self):
        pat = self._frame_pat()
        while not self._stop.wait(120):
            elapsed_min = int((time.time() - self._start) // 60)
            frame = 0
            try:
                text = self.log_file.read_text(errors="ignore")
                matches = pat.findall(text)
                if matches:
                    frame = int(matches[-1])
            except Exception:
                frame = 0
            pct = frame * 100 // self.total
            _echo(f"    [PROGRESS] {self.label}: frame {frame}/{self.total} "
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


# ===============================================================================
# SECTION 5: RUN-DIRECTORY RESOLUTION
# ===============================================================================
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
    _echo("Usage: python 07_SID_Post_Processing_FAcDs.py [Boltz-2_Run_Directory]")
    sys.exit(1)


# ===============================================================================
# SECTION 6: SCAN PHASE
# ===============================================================================
def scan_jobs(job_dirs: list[Path], md_dir: Path) -> list[Path]:
    """Classify every desmond_md_job_Rank_* directory and return the subset
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

        # 1. No -out.cms => MD simulation not finished. A live process here means
        #    the SIMULATION is running (job-name match avoids SID processes).
        if not cms_file.is_file():
            if _proc_alive(job_name):
                running += 1
                _echo(f"  - Rank {rank}: RUNNING (simulation active)")
            else:
                pending += 1
                _echo(f"  - Rank {rank}: PENDING (simulation not finished yet)")
            continue

        # 2. -out.cms present => MD done. Classify by SID-out.eaf completeness,
        #    measured against the actual trajectory length.
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


# ===============================================================================
# SECTION 7: SID ANALYSIS PER JOB
# ===============================================================================
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
    """Iterate the work list, generating SID-in then SID-out for each job."""
    total = len(to_run)
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
        run_event_analysis(d, job_name, cms_file, in_eaf)
        total_frames = run_analyze_simulation(d, job_name, cms_file, trj_dir, out_eaf, label)

        if out_eaf.is_file() and is_eaf_complete(out_eaf_frames(out_eaf), total_frames):
            of = out_eaf_frames(out_eaf)
            _echo(f"  [SUCCESS] Created output: {out_eaf} ({of}/{total_frames} frames)")
            _echo(f"  [LOGS] Logs written to {job_name}_event_analysis.log "
                  f"and {job_name}_analyze_simulation.log")
        else:
            _echo(f"  [ERROR] Failed to generate complete output: {out_eaf}")
            sys.exit(1)


# ===============================================================================
# SECTION 8: MAIN
# ===============================================================================
def main() -> int:
    parser = argparse.ArgumentParser(
        description=f"{CFG.PROJECT_NAME} Step 07 — Desmond SID Post-Processing")
    parser.add_argument("run_dir", nargs="?", default=None,
                        help="Boltz-2_Run_* directory (auto-detect latest if omitted)")
    parser.add_argument("--pipeline-mode", action="store_true",
                        help="Set when called from 00_00_run_pipeline_FAcDs.sh; "
                             "delegates systemd-oomd masking to the pipeline runner.")
    args = parser.parse_args()

    t0 = time.time()
    print_script_banner("07_SID_Post_Processing_FAcDs.py",
                        "Step 07 — Desmond SID Post-Processing")

    if not os.path.isfile(SCHROD_RUN):
        _echo(f"ERROR: Schrödinger 'run' not found at {SCHROD_RUN}. "
              f"Set the SCHRODINGER environment variable.")
        return 1

    run_dir = resolve_run_dir(args.run_dir)
    md_dir = _SCRIPT_DIR / run_dir / "7_Physics_Validation" / "MolecularDynamics"

    if not md_dir.is_dir():
        _echo(f"ERROR: MolecularDynamics directory not found: {md_dir}")
        return 1

    _echo(_SEP)
    _echo("Starting Desmond SID Analysis Post-processing")
    _echo(f"Run Directory      : {run_dir}")
    _echo(f"Molecular Dynamics : {md_dir}")
    _echo(_SEP)

    job_dirs = sorted(
        (d for d in md_dir.iterdir()
         if d.is_dir() and d.name.startswith("desmond_md_job_Rank_")),
        key=_natural_rank,
    )
    if not job_dirs:
        _echo("No job directories found matching pattern: desmond_md_job_Rank_*")
        return 0

    to_run = scan_jobs(job_dirs, md_dir)
    if not to_run:
        _echo("No jobs ready for SID analysis. Exiting.")
        return 0

    # Mask oomd (standalone only); the context manager restores it on exit.
    with OomdGuard(active=not args.pipeline_mode):
        process_jobs(to_run)

    _echo("")
    _echo(_SEP)
    _echo("All available Desmond SID Analysis jobs completed successfully.")
    _echo(_SEP)
    print_elapsed(t0, "07_SID_Post_Processing_FAcDs.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
