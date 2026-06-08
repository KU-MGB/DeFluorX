#!/usr/bin/env bash
# =============================================================================
# FAcDs Pipeline Runner
# Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
# Date   : 10 June 2026
# =============================================================================
# Usage:
#   bash 00_00_run_pipeline_FAcDs.sh [--run-id=<name>] [--dry-run] [--resume-from=<N>]
#
# Run mode (auto-detected when --run-id is omitted):
#   Resume  : Boltz-2_<name>/ directory found   → python 02 --resume <run_id>
#   Fresh   : no run directory found             → python 02 --fasta <fasta> --smi <smi>
#             FASTA = C_Final_Merged_for_Boltz-2.fasta (step-01 output)
#             SMI   = first *.smi file found in repo root
#
# Monitor progress in a second terminal:
#   tail -f <log file printed at start>
#
# All stdout + stderr from every step is written to the log file.
# Step headers are echoed to both terminal and log for live progress tracking.
# Pipeline halts immediately on any step failure.
# A timing summary table is printed at the end.
#
# ── The Critic's Corner: Known Limitations & Failure Points ──────────────────
#   1. Sequential Execution: Steps are strictly ordered; if Step 02 fails,
#      downstream analysis (03-07) cannot be launched until fixed.
#   2. Log Interleaving: Concurrent runs in the same directory will interleave
#      output in the same log file unless unique --run-id is provided.
#   3. Environment: Assumes the 'PFAS' conda environment is correctly configured
#      via Step 00; lacks internal dependency verification.
#   4. Path Assumptions: Relies on `C_Final_Merged_for_Boltz-2.fasta` from Step 01
#      being present for fresh runs.
# ─────────────────────────────────────────────────────────────────────────────
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Argument parsing ──────────────────────────────────────────────────────────

DRY_RUN=0
RESUME_FROM=0
RUN_ID=""

for _arg in "$@"; do
    case "$_arg" in
        --dry-run)       DRY_RUN=1 ;;
        --resume-from=*) RESUME_FROM="${_arg#*=}" ;;
        --run-id=*)      RUN_ID="${_arg#*=}" ;;
    esac
done

# ── RUN_ID resolution ──────────────────────────────────────────────────────────
# Priority: (1) explicit --run-id dir exists → resume
#           (2) any Boltz-2_* dir exists     → resume latest found
#           (3) neither                       → fresh run (needs *.smi; FASTA from step 01)

_PIPELINE_MODE="resume"
_FRESH_FASTA="C_Final_Merged_for_Boltz-2.fasta"
_FRESH_SMI=""

if [[ -n "$RUN_ID" && ! "$RUN_ID" =~ ^[A-Za-z0-9_.-]+$ ]]; then
    echo "ERROR: --run-id may only contain letters, digits, underscores, hyphens, and dots." >&2
    exit 1
fi

if [[ -n "$RUN_ID" && -d "${SCRIPT_DIR}/${RUN_ID}" ]]; then
    _PIPELINE_MODE="resume"
else
    _auto=$(find "${SCRIPT_DIR}" -maxdepth 1 -type d -name 'Boltz-2_*' 2>/dev/null | sort | tail -1)
    if [[ -n "$_auto" ]]; then
        RUN_ID="$(basename "$_auto")"
        _PIPELINE_MODE="resume"
    else
        _PIPELINE_MODE="fresh"
        _FRESH_SMI=$(find "${SCRIPT_DIR}" -maxdepth 1 -name '*.smi' 2>/dev/null | sort | tail -1)
        if [[ -z "$_FRESH_SMI" ]]; then
            echo "ERROR: Fresh run requires a .smi ligand file in ${SCRIPT_DIR}" >&2
            exit 1
        fi
    fi
fi

LOG_DIR="${SCRIPT_DIR}/${RUN_ID}/0_FAcDs_Pipeline_Logs"
[[ "$_PIPELINE_MODE" == "fresh" ]] && LOG_DIR="${SCRIPT_DIR}/0_FAcDs_Pipeline_Staging_Logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"

# Redirect all stdout+stderr through a single tee — one write to log, one to terminal.
exec > >(tee >(sed 's/\x1B\[[0-9;]*[mKABCDEFGHJKSTfhilmnprsu]//g' >> "$LOG_FILE")) 2>&1

# ── Helpers ──────────────────────────────────────────────────────────────────

_sep="════════════════════════════════════════════════════════════════════════════════"

# Write to log + terminal (exec tee above handles the split)
_log() { echo "$@"; }
_tee() { echo "$@"; }

# Timing registry — parallel arrays
STEP_NAMES=()
STEP_TIMES=()
STEP_STATUS=()

_fmt_elapsed() {
    local s=$1
    if   (( s < 60   )); then echo "${s}s"
    elif (( s < 3600 )); then printf "%dm%02ds" $(( s/60 )) $(( s%60 ))
    else                       printf "%dh%02dm%02ds" $(( s/3600 )) $(( s%3600/60 )) $(( s%60 ))
    fi
}

run_step() {
    local name="$1"; shift
    # Extract the leading numeric prefix (e.g. "07" from "07  MD …") as the step number.
    # 10# forces base-10 so "08"/"09" never trip octal parsing; :-0 guards a non-numeric prefix.
    local step_num _digits="${name%%[^0-9]*}"
    step_num=$(( 10#${_digits:-0} ))

    _tee ""
    _tee "$_sep"
    _tee "  STEP : ${name}"

    # --dry-run: print what would run, skip execution.
    if [[ $DRY_RUN -eq 1 ]]; then
        _tee "  [DRY-RUN] Would execute: $*"
        _tee "$_sep"
        return 0
    fi

    # --resume-from N: skip steps whose number is strictly less than N.
    if (( step_num < RESUME_FROM )); then
        _tee "  [SKIP] Step ${step_num} (--resume-from=${RESUME_FROM})"
        _tee "$_sep"
        STEP_NAMES+=("$name"); STEP_TIMES+=(0); STEP_STATUS+=("SKIP")
        return 0
    fi

    _tee "  CMD  : $*"
    _tee "  START: $(date '+%Y-%m-%d %H:%M:%S')"
    _tee "$_sep"

    local t0=$SECONDS
    local elapsed=0
    local status="PASS"
    if "$@"; then
        elapsed=$(( SECONDS - t0 ))
        STEP_NAMES+=("$name"); STEP_TIMES+=("$elapsed"); STEP_STATUS+=("$status")
        _tee "  [PASS] ${name}  ($(_fmt_elapsed $elapsed))"
    else
        local exit_code=$?
        elapsed=$(( SECONDS - t0 ))
        status="FAIL"
        STEP_NAMES+=("$name"); STEP_TIMES+=("$elapsed"); STEP_STATUS+=("$status")
        _tee "  [FAIL] ${name}  (exit ${exit_code})"
        _tee ""
        _tee "  Pipeline halted. Tail the log for details:"
        _tee "    tail -n 50 ${LOG_FILE}"
        _print_timing_table
        exit "$exit_code"
    fi
}

_print_timing_table() {
    local total=0
    _tee ""
    _tee "$_sep"
    _tee "  TIMING SUMMARY"
    _tee "$_sep"
    _tee "  $(printf '%-36s  %10s  %s' 'STEP' 'ELAPSED' 'STATUS')"
    _tee "  $(printf '%-36s  %10s  %s' '----' '-------' '------')"
    for i in "${!STEP_NAMES[@]}"; do
        local t="${STEP_TIMES[$i]}"
        local st="${STEP_STATUS[$i]}"
        _tee "  $(printf '%-36s  %10s  %s' "${STEP_NAMES[$i]}" "$(_fmt_elapsed $t)" "$st")"
        total=$(( total + t ))
    done
    _tee "  $(printf '%-36s  %10s' '──────────────────────────────────' '──────────')"
    _tee "  $(printf '%-36s  %10s' 'TOTAL WALL TIME' "$(_fmt_elapsed $total)")"
    _tee "$_sep"
}

# ── Conda activation ──────────────────────────────────────────────────────────

if [[ -z "${CONDA_BASE:-}" ]]; then
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
fi
if [[ -z "${CONDA_BASE:-}" ]]; then
    CONDA_BASE="$HOME/miniconda3"
fi

if [[ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate PFAS
else
    _tee "WARNING: conda not found at ${CONDA_BASE}. Set CONDA_BASE env var or install miniconda."
fi

PIPELINE_START=$SECONDS

# ── Sudo credential cache ─────────────────────────────────────────────────────
# Prompt for sudo password NOW (while user is present) so the cached credential
# is available for the systemd-oomd mask/unmask commands that wrap step 07.
# The four systemctl commands also have NOPASSWD in sudoers as a belt-and-braces
# guarantee — no expiry risk even if steps 01–06 run longer than 2–3 hours.
echo ""
echo "  Step 07 (MD thermodynamics) requires sudo to mask systemd-oomd."
echo "  systemctl commands are covered by NOPASSWD in sudoers."
# Non-interactive credential refresh — succeeds silently if already cached;
# falls through harmlessly if no TTY (nohup/background context).
sudo -vn 2>/dev/null || true
( while kill -0 $$ 2>/dev/null; do sudo -vn 2>/dev/null; sleep 60; done ) &
_SUDO_KEEPALIVE_PID=$!
trap 'kill "$_SUDO_KEEPALIVE_PID" 2>/dev/null || true; sudo systemctl unmask systemd-oomd.socket 2>/dev/null || true; sudo systemctl start systemd-oomd 2>/dev/null || true; true' EXIT
echo ""

# ── Header ────────────────────────────────────────────────────────────────────

_tee "$_sep"
_tee "  FAcDs Pipeline Run"
_tee "  Started : $(date '+%Y-%m-%d %H:%M:%S')"
_tee "  Env     : ${CONDA_DEFAULT_ENV:-unknown}"
_tee "  Python  : $(python --version 2>&1)"
_tee "  Mode    : ${_PIPELINE_MODE^^}"
if [[ "$_PIPELINE_MODE" == "resume" ]]; then
    _tee "  Run ID  : ${RUN_ID}"
else
    _tee "  FASTA   : ${_FRESH_FASTA}"
    _tee "  SMI     : ${_FRESH_SMI}"
fi
_tee "  Log     : ${LOG_FILE}"
_tee "$_sep"
echo ""
echo "  Monitor: tail -f ${LOG_FILE}"
echo ""

# ── Pipeline steps ────────────────────────────────────────────────────────────

run_step "00  Environment check" \
    python 00_01_Environment_Installation_FAcDs.py --export

run_step "01  Merge sequences" \
    python 01_Merge_FAcDs.py \
        --master    A_Labelled_15-Seq.fasta \
        --secondary B_Downloaded-Blast_Uniprot_NCBI.fasta \
        --output    C_Final_Merged_for_Boltz-2.fasta

if [[ "$_PIPELINE_MODE" == "resume" ]]; then
    run_step "02  Production (Boltz-2 scoring)" \
        python 02_Production_FAcDs.py --resume "$RUN_ID"
else
    run_step "02  Production (Boltz-2 scoring — fresh)" \
        python 02_Production_FAcDs.py --fasta "$_FRESH_FASTA" --smi "$_FRESH_SMI"
    _new_run=$(find "${SCRIPT_DIR}" -maxdepth 1 -type d -name 'Boltz-2_*' 2>/dev/null | sort | tail -1)
    if [[ -z "$_new_run" ]]; then
        _tee "  ERROR: Step 02 did not create a Boltz-2_* run directory."
        exit 1
    fi
    RUN_ID="$(basename "$_new_run")"
    _tee "  Fresh run directory: ${RUN_ID}"
fi

run_step "03  Validation figures" \
    python 03_Validation_Figures_FAcDs.py "$RUN_ID"

run_step "04  Phylogeny" \
    python 04_Phylogeny_FAcDs.py "$RUN_ID"

run_step "05  CIF/PDB preparation" \
    python 05_CIF-PDB_Preparation_FAcDs.py "$RUN_ID"

run_step "06  Top-N extraction" \
    python 06_Top-N_Extraction_FAcDs.py "$RUN_ID"

if [[ $DRY_RUN -eq 0 && 7 -ge $RESUME_FROM ]]; then
    sudo systemctl stop systemd-oomd && sudo systemctl mask systemd-oomd.socket
fi
run_step "07  MD thermodynamics + QM/MM engine" \
    python 07_MD_Thermodynamics_QMMM_Engine_FAcDs.py "$RUN_ID"
if [[ $DRY_RUN -eq 0 && 7 -ge $RESUME_FROM ]]; then
    sudo systemctl unmask systemd-oomd.socket && sudo systemctl start systemd-oomd
fi

# ── Footer ────────────────────────────────────────────────────────────────────

_tee ""
_tee "$_sep"
_tee "  ALL STEPS COMPLETE"
_tee "  Finished: $(date '+%Y-%m-%d %H:%M:%S')"
_tee "  Full log: ${LOG_FILE}"

_print_timing_table
