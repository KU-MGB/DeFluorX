#!/usr/bin/env bash
# =============================================================================
# FAcDs Pipeline  |  Post-Processing  |  Desmond SID Analysis Orchestrator
# =============================================================================
# Post-simulation Post-Processing script that runs Schrödinger's Event Analysis
# and Simulation Interaction Diagram (SID) post-processing on molecular dynamics
# trajectories sequentially.
#
# Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
# Date   : 10 June 2026
# =============================================================================
# Usage:
#   bash run_sid_analysis_FAcDs.sh [Boltz-2_Run_Directory]
#
# ── Dependency Map ───────────────────────────────────────────────────────────
#   Script        : run_sid_analysis_FAcDs.sh
#   Role          : Post-processing — Desmond post-simulation post-processing.
#   Imports from  : None.
#   Reads         : Boltz-2_Run_X/7_Physics_Validation/MolecularDynamics/desmond_md_job_Rank_N/*-out.cms
#                   Boltz-2_Run_X/7_Physics_Validation/MolecularDynamics/desmond_md_job_Rank_N/*_trj
#   Writes        : Boltz-2_Run_X/7_Physics_Validation/MolecularDynamics/desmond_md_job_Rank_N/*_SID-in.eaf
#                   Boltz-2_Run_X/7_Physics_Validation/MolecularDynamics/desmond_md_job_Rank_N/*_SID-out.eaf
#                   Boltz-2_Run_X/7_Physics_Validation/MolecularDynamics/desmond_md_job_Rank_N/*.log
#   Upstream      : Desmond molecular dynamics simulations (Step 07 preparation).
#   Downstream    : 07_MD_Thermodynamics_QMMM_Engine_FAcDs.py (consumes EAF output).
# ─────────────────────────────────────────────────────────────────────────────
# ── The Critic's Corner: Known Limitations & Failure Points ──────────────────
#   1. Sudo Requirement: Masking systemd-oomd requires root privileges via sudo.
#      The script will prompt for sudo access at startup.
#   2. Trajectory Volume: Processing 100,000 frames sequentially is time-intensive
#      (typically 10–15 hours). Do not run multiple instances concurrently.
#   3. Local I/O Capping: Relies on -LOCAL to prevent huge tmp partition writes;
#      requires sufficient disk space in the destination filesystem.
# ─────────────────────────────────────────────────────────────────────────────
# =============================================================================

set -euo pipefail

export SCHRODINGER=/opt/schrodinger

# Dynamically resolve root directory location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

RUN_DIR="${1:-}"

# Auto-detect latest Boltz-2_Run_* directory if not provided
if [ -z "$RUN_DIR" ]; then
    _auto=$(find . -maxdepth 1 -type d -name 'Boltz-2_Run_*' 2>/dev/null | sort | tail -1)
    if [ -n "$_auto" ]; then
        RUN_DIR=$(basename "$_auto")
        echo "Auto-detected run directory: $RUN_DIR"
    else
        echo "ERROR: Run directory not specified and no Boltz-2_Run_* found."
        echo "Usage: $0 [Boltz-2_Run_Directory]"
        exit 1
    fi
fi

MD_DIR="${SCRIPT_DIR}/${RUN_DIR}/7_Physics_Validation/MolecularDynamics"

if [ ! -d "$MD_DIR" ]; then
    echo "ERROR: MolecularDynamics directory not found: $MD_DIR"
    exit 1
fi

echo "============================================================================="
echo "Starting Desmond SID Analysis Post-processing"
echo "Run Directory      : $RUN_DIR"
echo "Molecular Dynamics : $MD_DIR"
echo "============================================================================="

# Find all desmond_md_job_Rank_X directories (numerically by Rank)
JOB_DIRS=$(find "$MD_DIR" -maxdepth 1 -type d -name 'desmond_md_job_Rank_*' | sort -V)

if [ -z "$JOB_DIRS" ]; then
    echo "No job directories found matching pattern: desmond_md_job_Rank_*"
    exit 0
fi

# ── Helpers ──────────────────────────────────────────────────────────────────
# Count frames actually written into a SID-out.eaf (the analysed result vector).
out_eaf_frames() {
    python3 -c "
import sys, re
try:
    with open(sys.argv[1]) as f:
        content = f.read()
    m = re.search(r'Result\s*=\s*\[([\d\.\-\s+eE]+)\]', content)
    print(len(m.group(1).split()) if m else 0)
except Exception:
    print(0)
" "$1" 2>/dev/null || echo 0
}

# Ground-truth frame count of a Desmond trajectory, read via the Schrodinger
# traj API. This is the denominator for completion — it adapts automatically
# whether the trajectory holds 1,000 or 100,000 frames. Returns 0 if unreadable.
traj_frame_count() {
    local _trj="$1"
    [ -d "$_trj" ] || { echo 0; return; }
    "$SCHRODINGER/run" python3 -c "
import sys
try:
    from schrodinger.application.desmond.packages import traj
    print(len(traj.read_traj(sys.argv[1])))
except Exception:
    print(0)
" "$_trj" 2>/dev/null || echo 0
}

# A SID-out.eaf is complete when it holds at least as many frames as the
# trajectory. If the trajectory count cannot be read (tf=0), an existing
# non-empty EAF is treated as complete — never destructively re-run on doubt.
is_eaf_complete() {
    local _of="$1" _tf="$2"
    if [ "$_tf" -gt 0 ]; then
        [ "$_of" -ge "$_tf" ]
    else
        [ "$_of" -gt 0 ]
    fi
}

# ── Scanning Phase ───────────────────────────────────────────────────────────
echo "Scanning job directories..."
TOTAL_JOBS=0
COMPLETED_JOBS=0
RUNNING_JOBS=0
PENDING_JOBS=0
TO_RUN_JOBS=0

for dir in $JOB_DIRS; do
    JOBNAME=$(basename "$dir")
    RANK="${JOBNAME##*_md_job_Rank_}"
    TOTAL_JOBS=$((TOTAL_JOBS + 1))

    JOBDIR="${MD_DIR}/${JOBNAME}"
    OUT_EAF="${JOBDIR}/${JOBNAME}_SID-out.eaf"
    CMS_FILE="${JOBDIR}/${JOBNAME}-out.cms"
    TRJ_DIR="${JOBDIR}/${JOBNAME}_trj"

    # 1. No -out.cms => MD simulation not finished. Only here does a live
    #    process mean the SIMULATION is running (avoids matching SID processes).
    if [ ! -f "$CMS_FILE" ]; then
        if pgrep -f "$JOBNAME" >/dev/null; then
            RUNNING_JOBS=$((RUNNING_JOBS + 1))
            echo "  - Rank $RANK: RUNNING (simulation active)"
        else
            PENDING_JOBS=$((PENDING_JOBS + 1))
            echo "  - Rank $RANK: PENDING (simulation not finished yet)"
        fi
        continue
    fi

    # 2. -out.cms present => MD done. Classify by SID-out.eaf completeness,
    #    measured against the actual trajectory length.
    if [ -f "$OUT_EAF" ]; then
        _of=$(out_eaf_frames "$OUT_EAF")
        _tf=$(traj_frame_count "$TRJ_DIR")
        if is_eaf_complete "$_of" "$_tf"; then
            COMPLETED_JOBS=$((COMPLETED_JOBS + 1))
            echo "  - Rank $RANK: COMPLETED (SID-out.eaf has $_of/$_tf frames)"
        else
            TO_RUN_JOBS=$((TO_RUN_JOBS + 1))
            echo "  - Rank $RANK: INCOMPLETE (SID-out.eaf has $_of/$_tf frames; will re-run)"
        fi
    else
        TO_RUN_JOBS=$((TO_RUN_JOBS + 1))
        echo "  - Rank $RANK: READY (simulation finished, SID pending)"
    fi
done

echo ""
echo "Scan Summary:"
echo "  Total Jobs Found : $TOTAL_JOBS"
echo "  Completed Jobs   : $COMPLETED_JOBS"
echo "  Running Jobs     : $RUNNING_JOBS"
echo "  Pending Jobs     : $PENDING_JOBS"
echo "  Jobs to Process  : $TO_RUN_JOBS"
echo "============================================================================="

if [ "$TO_RUN_JOBS" -eq 0 ]; then
    echo "No jobs ready for SID analysis. Exiting."
    exit 0
fi

# ── Mask systemd-oomd to prevent early termination ───────────────────────────
_OOMD_RESTORED=0
restore_oomd() {
    if [ "$_OOMD_RESTORED" -eq 0 ]; then
        echo "Restoring systemd-oomd services..."
        sudo systemctl unmask systemd-oomd.socket 2>/dev/null || true
        sudo systemctl start systemd-oomd 2>/dev/null || true
        _OOMD_RESTORED=1
    fi
}

# Prompt for sudo password at the beginning (matching run_pipeline behaviour)
echo "This script requires sudo to mask systemd-oomd to prevent Out-Of-Memory kills."
sudo -v
( while kill -0 $$ 2>/dev/null; do sudo -vn 2>/dev/null; sleep 60; done ) &
_SUDO_KEEPALIVE_PID=$!
trap 'kill "$_SUDO_KEEPALIVE_PID" 2>/dev/null || true; restore_oomd' EXIT ERR INT TERM

echo "Masking systemd-oomd..."
sudo systemctl stop systemd-oomd 2>/dev/null || true
sudo systemctl mask systemd-oomd.socket 2>/dev/null || true

# ── Heartbeat: print progress of a long-running step every interval ──────────
# Usage: heartbeat <log_file> <label> & HB_PID=$!  ... then: kill "$HB_PID"
heartbeat() {
    local _log="$1" _label="$2" _total="${3:-100000}" _start _elapsed _frame _pct
    _start=$(date +%s)
    while true; do
        sleep 120
        _elapsed=$(( ($(date +%s) - _start) / 60 ))
        # analyze_simulation.py prints "analyzing frame# N..." per frame
        _frame=$(grep -oE 'analyzing frame# [0-9]+' "$_log" 2>/dev/null | tail -1 | grep -oE '[0-9]+' || echo 0)
        _frame=${_frame:-0}
        _pct=$(( _frame * 100 / _total ))
        echo "    [PROGRESS] ${_label}: frame ${_frame}/${_total} (${_pct}%, ${_elapsed}m elapsed)"
    done
}

# ── Run Loop ─────────────────────────────────────────────────────────────────
JOB_INDEX=0
for dir in $JOB_DIRS; do
    JOBNAME=$(basename "$dir")
    RANK="${JOBNAME##*_md_job_Rank_}"
    
    JOBDIR="${MD_DIR}/${JOBNAME}"
    OUT_EAF="${JOBDIR}/${JOBNAME}_SID-out.eaf"
    IN_EAF="${JOBDIR}/${JOBNAME}_SID-in.eaf"
    CMS_FILE="${JOBDIR}/${JOBNAME}-out.cms"
    TRJ_DIR="${JOBDIR}/${JOBNAME}_trj"

    # 1. Skip if already finished (EAF holds >= trajectory frame count).
    if [ -f "$OUT_EAF" ]; then
        _of=$(out_eaf_frames "$OUT_EAF")
        _tf=$(traj_frame_count "$TRJ_DIR")

        if is_eaf_complete "$_of" "$_tf"; then
            continue
        else
            JOB_INDEX=$((JOB_INDEX + 1))
            echo ""
            echo "============================================================================="
            echo "[Job $JOB_INDEX/$TO_RUN_JOBS] Re-processing (Incomplete): $JOBNAME (Rank $RANK)"
            echo "============================================================================="
            echo "  [RE-RUN] Output EAF is incomplete ($_of/$_tf frames). Re-running analysis..."
            # Preserve the partial EAF instead of deleting it — a misjudged
            # "incomplete" must never destroy a good result.
            _bak="${OUT_EAF}.incomplete.$(date +%Y%m%d%H%M%S).bak"
            mv -f "$OUT_EAF" "$_bak"
            echo "  [BACKUP] Previous EAF moved to: $_bak"
        fi
    else
        # Verify readiness before echoing and entering directory
        if [ ! -f "$CMS_FILE" ] || pgrep -f "$JOBNAME" >/dev/null; then
            continue
        fi

        JOB_INDEX=$((JOB_INDEX + 1))
        echo ""
        echo "============================================================================="
        echo "[Job $JOB_INDEX/$TO_RUN_JOBS] Processing: $JOBNAME (Rank $RANK)"
        echo "============================================================================="
    fi

    # Verify trajectory folder is present
    if [ ! -d "$TRJ_DIR" ] && [ ! -f "${TRJ_DIR}.xtc" ]; then
        echo "  [WARNING] Trajectory folder not found: $TRJ_DIR. Skipping Rank $RANK."
        continue
    fi

    cd "$JOBDIR"

    # Step 1: event_analysis.py (generates SID-in.eaf)
    if [ ! -f "$IN_EAF" ]; then
        echo "  Running event_analysis.py to generate SID-in.eaf..."
        $SCHRODINGER/run event_analysis.py analyze \
            "$CMS_FILE" \
            -out "${JOBNAME}_SID" > "${JOBNAME}_event_analysis.log" 2>&1
    else
        echo "  [SKIP] SID-in.eaf already exists: $IN_EAF"
    fi

    # Step 2: analyze_simulation.py (generates SID-out.eaf)
    # Using local execution (-LOCAL) to avoid remote server overhead and cap memory
    # Resolve the trajectory length up front so the heartbeat reports a true
    # percentage and completion can be validated afterwards.
    TOTAL_FRAMES=$(traj_frame_count "$TRJ_DIR")
    [ "$TOTAL_FRAMES" -gt 0 ] || TOTAL_FRAMES=100000
    echo "  Running analyze_simulation.py (long step, ~10-15h; $TOTAL_FRAMES frames)..."
    heartbeat "${JOBDIR}/${JOBNAME}_analyze_simulation.log" "Job $JOB_INDEX/$TO_RUN_JOBS Rank $RANK" "$TOTAL_FRAMES" &
    HB_PID=$!
    $SCHRODINGER/run analyze_simulation.py -NOJOBID -LOCAL \
        "$CMS_FILE" \
        "$TRJ_DIR" \
        "${JOBNAME}_SID-out.eaf" \
        "${JOBNAME}_SID-in.eaf" > "${JOBNAME}_analyze_simulation.log" 2>&1
    kill "$HB_PID" 2>/dev/null || true
    wait "$HB_PID" 2>/dev/null || true

    if [ -f "$OUT_EAF" ] && is_eaf_complete "$(out_eaf_frames "$OUT_EAF")" "$TOTAL_FRAMES"; then
        echo "  [SUCCESS] Created output: $OUT_EAF ($(out_eaf_frames "$OUT_EAF")/$TOTAL_FRAMES frames)"
        echo "  [LOGS] Logs written to ${JOBNAME}_event_analysis.log and ${JOBNAME}_analyze_simulation.log"
    else
        echo "  [ERROR] Failed to generate complete output: $OUT_EAF"
        exit 1
    fi

    # Return to root directory
    cd "$SCRIPT_DIR"
done

echo ""
echo "============================================================================="
echo "All available Desmond SID Analysis jobs completed successfully."
echo "============================================================================="
