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
    
    # 1. Check if process is active in system
    if pgrep -f "$JOBNAME" >/dev/null; then
        RUNNING_JOBS=$((RUNNING_JOBS + 1))
        echo "  - Rank $RANK: RUNNING (simulation active)"
        continue
    fi
    
    # 2. Check if output EAF exists and is complete
    if [ -f "$OUT_EAF" ]; then
        _frames=$(python3 -c "
import sys, re
try:
    with open(sys.argv[1]) as f:
        content = f.read()
    m = re.search(r'Result\s*=\s*\[([\d\.\-\s+eE]+)\]', content)
    print(len(m.group(1).split()) if m else 0)
except Exception:
    print(0)
" "$OUT_EAF" 2>/dev/null || echo 0)
        
        if [ "$_frames" -ge 99990 ]; then
            COMPLETED_JOBS=$((COMPLETED_JOBS + 1))
            echo "  - Rank $RANK: COMPLETED (EAF exists with $_frames frames)"
        else
            TO_RUN_JOBS=$((TO_RUN_JOBS + 1))
            echo "  - Rank $RANK: INCOMPLETE (EAF exists but has $_frames/100000 frames; will re-run)"
        fi
        continue
    fi
    
    # 3. Check if CMS exists (simulation finished but SID not run)
    if [ -f "$CMS_FILE" ]; then
        TO_RUN_JOBS=$((TO_RUN_JOBS + 1))
        echo "  - Rank $RANK: READY (simulation finished, SID pending)"
    else
        PENDING_JOBS=$((PENDING_JOBS + 1))
        echo "  - Rank $RANK: PENDING (simulation not finished yet)"
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

# ── Run Loop ─────────────────────────────────────────────────────────────────
for dir in $JOB_DIRS; do
    JOBNAME=$(basename "$dir")
    RANK="${JOBNAME##*_md_job_Rank_}"
    
    JOBDIR="${MD_DIR}/${JOBNAME}"
    OUT_EAF="${JOBDIR}/${JOBNAME}_SID-out.eaf"
    IN_EAF="${JOBDIR}/${JOBNAME}_SID-in.eaf"
    CMS_FILE="${JOBDIR}/${JOBNAME}-out.cms"
    TRJ_DIR="${JOBDIR}/${JOBNAME}_trj"

    # 1. Skip if already finished (output exists with >= 100000 frames)
    if [ -f "$OUT_EAF" ]; then
        _frames=$(python3 -c "
import sys, re
try:
    with open(sys.argv[1]) as f:
        content = f.read()
    m = re.search(r'Result\s*=\s*\[([\d\.\-\s+eE]+)\]', content)
    print(len(m.group(1).split()) if m else 0)
except Exception:
    print(0)
" "$OUT_EAF" 2>/dev/null || echo 0)
        
        if [ "$_frames" -ge 99990 ]; then
            continue
        else
            echo ""
            echo "============================================================================="
            echo "Re-processing Job (Incomplete): $JOBNAME (Rank $RANK)"
            echo "============================================================================="
            echo "  [RE-RUN] Output EAF is incomplete ($_frames/100000 frames). Re-running analysis..."
            rm -f "$OUT_EAF"
        fi
    else
        # Verify readiness before echoing and entering directory
        if pgrep -f "$JOBNAME" >/dev/null || [ ! -f "$CMS_FILE" ]; then
            continue
        fi
        
        echo ""
        echo "============================================================================="
        echo "Processing Job: $JOBNAME (Rank $RANK)"
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
    echo "  Running analyze_simulation.py..."
    $SCHRODINGER/run analyze_simulation.py -NOJOBID -LOCAL \
        "$CMS_FILE" \
        "$TRJ_DIR" \
        "${JOBNAME}_SID-out.eaf" \
        "${JOBNAME}_SID-in.eaf" > "${JOBNAME}_analyze_simulation.log" 2>&1

    if [ -f "$OUT_EAF" ]; then
        echo "  [SUCCESS] Created output: $OUT_EAF"
        echo "  [LOGS] Logs written to ${JOBNAME}_event_analysis.log and ${JOBNAME}_analyze_simulation.log"
    else
        echo "  [ERROR] Failed to generate output: $OUT_EAF"
        exit 1
    fi

    # Return to root directory
    cd "$SCRIPT_DIR"
done

echo ""
echo "============================================================================="
echo "All available Desmond SID Analysis jobs completed successfully."
echo "============================================================================="
