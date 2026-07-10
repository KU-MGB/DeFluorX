#!/usr/bin/env bash
# =============================================================================
# FAcDs Pipeline Runner
# Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
# Date   : 10 July 2026
# =============================================================================
# Usage (non-interactive / scripted mode):
#   bash 00_00_run_pipeline_FAcDs.sh [--run-id=<name>] [--dry-run] [--resume-from=<N>]
#
# Usage (interactive — default):
#   bash 00_00_run_pipeline_FAcDs.sh
#   → prompted: Fresh or Resume?  (+ optional sudo password — press Enter to skip)
#   → runs fully unattended thereafter
#
# Run modes:
#   Interactive : Script asks "Fresh or Resume" and (if Resume) which run dir.
#   Non-interactive (--run-id supplied):
#     --run-id directory exists → resume that run
#     ERROR if directory not found (fresh runs require the interactive prompt)
#
# Fresh run:
#   Step 01 merges sequences → C_INP_Merged_for_Boltz-2.fasta
#   Step 02 creates a new Boltz-2_<timestamp>/ directory; its name is
#   auto-captured and passed to all downstream steps (03–07).
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
#   4. Fresh mode RUN_ID: The run directory created by Step 02 is auto-detected
#      after that step completes; if Step 02 creates no directory, the pipeline
#      halts with a diagnostic message.
#   5. --resume-from with fresh mode: Not supported (no existing run to skip
#      into). Use --resume-from only together with --run-id.
# ─────────────────────────────────────────────────────────────────────────────
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Argument parsing ──────────────────────────────────────────────────────────

DRY_RUN=0
RESUME_FROM=0
RUN_ID=""
# Exit code a step may return to mean "completed, but a complementary part was
# deferred/failed" (e.g. 06 when MM-GBSA is skipped or a Prime job fails). The
# runner renders this as WARN and continues rather than a false PASS or a hard
# abort. Keep in sync with EXIT_WARN in 06_SID_Prime-MMGBSA_FAcDs.py.
readonly WARN_EXIT_CODE=3

# Terminal colours for the interactive prelude (mode/sudo prompts). Only emitted
# to a real TTY; the log copy is ANSI-stripped downstream, so these never pollute
# the file. Empty on non-interactive stdout (pipe/redirect).
if [[ -t 1 ]]; then
    _C_GREEN=$'\033[0;32m'; _C_RED=$'\033[0;31m'; _C_YELLOW=$'\033[0;33m'
    _C_BOLD=$'\033[1m';     _C_RESET=$'\033[0m'
else
    _C_GREEN=""; _C_RED=""; _C_YELLOW=""; _C_BOLD=""; _C_RESET=""
fi

for _arg in "$@"; do
    case "$_arg" in
        --dry-run)       DRY_RUN=1 ;;
        --resume-from=*) RESUME_FROM="${_arg#*=}" ;;
        --run-id=*)      RUN_ID="${_arg#*=}" ;;
    esac
done

if [[ -n "$RUN_ID" && ! "$RUN_ID" =~ ^[A-Za-z0-9_.-]+$ ]]; then
    echo "ERROR: --run-id may only contain letters, digits, underscores, hyphens, and dots." >&2
    exit 1
fi

# Validate --resume-from is a non-negative integer at parse time; otherwise the later
# integer comparisons ([[ 7 -ge $RESUME_FROM ]]) throw "integer expression expected"
# and abort under set -e.
if [[ ! "$RESUME_FROM" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --resume-from must be a non-negative integer." >&2
    exit 1
fi

# ── Canonical filenames ───────────────────────────────────────────────────────
_FRESH_FASTA="C_INP_Merged_for_Boltz-2.fasta"
_FRESH_SMI="D_INP_PFAS-27_Ligands.smi"

# ── Interactive run-mode selection ────────────────────────────────────────────
# Skipped when --run-id is supplied (non-interactive / scripted context).
_PIPELINE_MODE="fresh"   # default; overridden below

if [[ -n "$RUN_ID" ]]; then
    # Non-interactive path: directory must exist for a resume
    if [[ -d "${SCRIPT_DIR}/${RUN_ID}" ]]; then
        _PIPELINE_MODE="resume"
    else
        echo "ERROR: --run-id '${RUN_ID}' directory not found." >&2
        echo "       For a fresh run, omit --run-id and use the interactive prompt." >&2
        exit 1
    fi
else
    # Collect available Boltz-2 run directories (sorted oldest → newest)
    _available_runs=()
    while IFS= read -r _d; do
        _available_runs+=("$(basename "$_d")")
    done < <(find "${SCRIPT_DIR}" -maxdepth 1 -type d -name 'Boltz-2_*' 2>/dev/null | sort)

    echo ""
    echo "  ══════════════════════════════════════════════════════════════════════════════"
    echo "  FAcDs Pipeline — Run Mode Selection"
    echo "  ══════════════════════════════════════════════════════════════════════════════"
    echo ""

    if [[ ${#_available_runs[@]} -gt 0 ]]; then
        echo "  Existing Boltz-2 run directories:"
        for i in "${!_available_runs[@]}"; do
            printf "    [%d]  %s\n" "$((i+1))" "${_available_runs[$i]}"
        done
        echo ""
        echo "  Select mode:"
        echo "    [F]  Fresh   — start a new prediction run from scratch"
        echo "    [R]  Resume  — continue from an existing run (default: latest)"
        echo ""
        read -r -p "  Your choice [F/R, default=R]: " _mode_choice </dev/tty
        _mode_choice="${_mode_choice:-R}"
    else
        echo "  No existing Boltz-2 run directories found."
        echo "  A fresh run is required."
        echo ""
        _mode_choice="F"
    fi

    case "${_mode_choice^^}" in
        F|FRESH)
            _PIPELINE_MODE="fresh"
            echo "  Mode: FRESH — a new run directory will be created by Step 02."
            ;;
        R|RESUME)
            _PIPELINE_MODE="resume"
            if [[ ${#_available_runs[@]} -eq 0 ]]; then
                echo "  ERROR: Resume selected but no Boltz-2 run directories exist." >&2
                echo "         Re-run without --run-id and select Fresh." >&2
                exit 1
            elif [[ ${#_available_runs[@]} -eq 1 ]]; then
                RUN_ID="${_available_runs[0]}"
                echo "  Mode: RESUME — auto-selected: ${RUN_ID}"
            else
                echo ""
                echo "  Enter the number of the run to resume"
                printf "  [1–%d, default=%d for latest]: " \
                    "${#_available_runs[@]}" "${#_available_runs[@]}"
                read -r _run_idx </dev/tty
                _run_idx="${_run_idx:-${#_available_runs[@]}}"
                if [[ "$_run_idx" =~ ^[0-9]+$ \
                   && "$_run_idx" -ge 1 \
                   && "$_run_idx" -le "${#_available_runs[@]}" ]]; then
                    RUN_ID="${_available_runs[$((_run_idx-1))]}"
                else
                    echo "  Invalid selection — defaulting to latest."
                    RUN_ID="${_available_runs[-1]}"
                fi
                echo "  Mode: RESUME — selected: ${RUN_ID}"
            fi
            ;;
        *)
            echo "  ERROR: Unrecognised choice '${_mode_choice}'. Enter F or R." >&2
            exit 1
            ;;
    esac
fi

# ── Validate prerequisites ────────────────────────────────────────────────────
if [[ "$_PIPELINE_MODE" == "fresh" ]]; then
    if [[ ! -f "${SCRIPT_DIR}/${_FRESH_SMI}" ]]; then
        echo "" >&2
        echo "ERROR: Ligand SMILES file not found: ${_FRESH_SMI}" >&2
        echo "       Expected at: ${SCRIPT_DIR}/${_FRESH_SMI}" >&2
        exit 1
    fi
fi

if [[ "$_PIPELINE_MODE" == "fresh" && "$RESUME_FROM" -gt 0 ]]; then
    echo "" >&2
    echo "WARNING: --resume-from=${RESUME_FROM} is not meaningful for a fresh run" >&2
    echo "         (no existing run to skip into). Flag will be ignored." >&2
    RESUME_FROM=0
fi

# ── Resume from a chosen step (interactive; resume mode only) ──────────────────
# Lets the user skip already-completed early steps (e.g. jump straight to 03 or 06
# without re-running the hours-long Step 02). Steps below the chosen number are
# marked SKIP. Prompted only when resuming interactively and --resume-from was not
# already given on the command line; Enter (default) runs the whole pipeline.
if [[ "$_PIPELINE_MODE" == "resume" && "$RESUME_FROM" -eq 0 && -t 0 ]]; then
    echo ""
    echo "  ${_C_BOLD}── Resume from which step? ──${_C_RESET}"
    echo "    [1] 01  Merge sequences"
    echo "    [2] 02  Production (Boltz-2 scoring)      ${_C_YELLOW}← heaviest${_C_RESET}"
    echo "    [3] 03  Validation figures"
    echo "    [4] 04  Dendrogram"
    echo "    [5] 05  Top-N selection + PDB preparation"
    echo "    [6] 06  SID + Prime MM-GBSA"
    echo "    [7] 07  MD thermodynamics + QM/MM engine"
    read -r -p "  Start from step [1-7, default=Enter = run everything]: " _from_step </dev/tty || _from_step=""
    if [[ -z "$_from_step" ]]; then
        echo "  ${_C_GREEN}→ Running the full pipeline (all steps).${_C_RESET}"
    elif [[ "$_from_step" =~ ^[1-7]$ ]]; then
        RESUME_FROM="$_from_step"
        printf "  ${_C_GREEN}→ Resuming from Step %02d onwards; Steps below %02d will be skipped.${_C_RESET}\n" \
            "$_from_step" "$_from_step"
    else
        echo "  ${_C_YELLOW}→ Invalid selection '${_from_step}' — running the full pipeline.${_C_RESET}"
    fi
fi

# ── Sudo credential cache (OPTIONAL) ──────────────────────────────────────────
# systemd-oomd masking stops the Linux OOM-killer from terminating the
# memory-heavy Steps 06-07. It needs sudo, which is OPTIONAL: not every user has
# sudo on every machine. Enter the password to enable masking (recommended for
# long unattended runs); press Enter to skip and run in NORMAL mode (no sudo).
#
# Scripts/operations that use sudo (only when enabled):
#   06_SID_Prime-MMGBSA_FAcDs.py              → systemctl stop / mask systemd-oomd
#   07_MD_QMMM_Defluorination_FAcDs.py → systemctl mask / unmask / start systemd-oomd
SUDO_ENABLED=0
_SUDO_KEEPALIVE_PID=""
echo ""
echo "  ── Optional sudo: systemd-oomd masking for Steps 06-07 ──"
echo "    Sudo is used only by: 06_SID_Prime-MMGBSA_FAcDs.py and 07_MD_QMMM_Defluorination_FAcDs.py"
echo "    (systemctl stop/mask/unmask/start systemd-oomd around the OOM-prone phase)."
echo ""
# Pressing Enter (empty) skips silently → NORMAL mode. A NON-EMPTY entry is treated as a
# password attempt: a wrong password is flagged and re-prompted (up to 3 attempts) so a typo
# does not silently drop a long unattended run into NORMAL mode. Three failures (or no sudo
# rights) fall back to NORMAL mode.
_sudo_attempt=0
while true; do
    read -r -s -p "  Enter sudo password to enable oomd masking, or press Enter to skip (NORMAL mode): " _sudo_pw </dev/tty || _sudo_pw=""
    echo ""
    if [[ -z "$_sudo_pw" ]]; then
        echo "  ${_C_YELLOW}→ No password entered — running in NORMAL mode (systemd-oomd NOT masked).${_C_RESET}"
        break
    elif echo "$_sudo_pw" | sudo -S -v 2>/dev/null; then
        SUDO_ENABLED=1
        echo "  ${_C_GREEN}→ sudo enabled — systemd-oomd will be masked during Steps 06-07.${_C_RESET}"
        # Keepalive: refresh credential every 60 s for long-running pipelines
        ( while kill -0 $$ 2>/dev/null; do sudo -vn 2>/dev/null; sleep 60; done ) &
        _SUDO_KEEPALIVE_PID=$!
        break
    else
        _sudo_attempt=$(( _sudo_attempt + 1 ))
        if (( _sudo_attempt >= 3 )); then
            echo "  ${_C_RED}→ Wrong password 3 times (or no sudo access) — running in NORMAL mode (systemd-oomd NOT masked).${_C_RESET}"
            break
        fi
        echo "  ${_C_RED}→ Wrong password. Try again, or press Enter to skip (NORMAL mode).  [attempt ${_sudo_attempt}/3]${_C_RESET}"
    fi
done
unset _sudo_pw
echo ""

# ── Log directory setup ───────────────────────────────────────────────────────
if [[ "$_PIPELINE_MODE" == "resume" ]]; then
    LOG_DIR="${SCRIPT_DIR}/${RUN_ID}/0_FAcDs_Pipeline_Logs"
else
    # Fresh: staging directory until Step 02 creates the real run directory
    LOG_DIR="${SCRIPT_DIR}/0_FAcDs_Pipeline_Staging_Logs"
fi
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"

# Fresh runs stream to a staging log until Step 02 creates the real run directory;
# _STAGING_SYNC_DEST is then set to the run-dir copy so it can be re-synced after
# every downstream step (and on exit), keeping the run-dir log complete — not
# truncated at Step 02. Resume runs already log straight into the run dir.
_STAGING_SYNC_DEST=""
_sync_staging_log() {
    [[ -n "$_STAGING_SYNC_DEST" && -f "$LOG_FILE" ]] && cp "$LOG_FILE" "$_STAGING_SYNC_DEST" 2>/dev/null || true
}

# ── Background / foreground selection ─────────────────────────────────────────
# Foreground (default): live output streamed to terminal + log via tee.
# Background: detach the run, write to the log only, and return the shell so the
#             terminal can be closed. Steps 06-07 sudo calls rely on the NOPASSWD
#             sudoers entries, so the detached (tty-less) process still works.
_BG_MODE=0
echo ""
read -r -p "  Run unattended in background (detach; log only, default = N)? [y/N]: " _bg_choice </dev/tty || _bg_choice=""
case "${_bg_choice^^}" in
    Y|YES) _BG_MODE=1 ;;
esac

# ── Helpers ──────────────────────────────────────────────────────────────────

_sep="════════════════════════════════════════════════════════════════════════════════"

# _tee is a thin wrapper around echo: it does NOT itself fork to a file.
# Terminal/log duplication is handled once, globally, by the `exec > >(tee ...)`
# redirect installed above — so every echo (and all child-process output) is
# already mirrored to both the terminal and ${LOG_FILE}. The wrapper exists only
# to give pipeline-progress lines a single, greppable call site.
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
    # Extract leading numeric prefix as step number (base-10 safe)
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
    # Capture the exit code without tripping `set -e` (the `||` guards it). A step
    # may return 0 (PASS), WARN_EXIT_CODE (completed with a deferred/failed
    # complementary part → WARN, continue), or any other non-zero (hard FAIL → halt).
    local exit_code=0
    "$@" || exit_code=$?
    elapsed=$(( SECONDS - t0 ))
    if (( exit_code == 0 )); then
        status="PASS"
        STEP_NAMES+=("$name"); STEP_TIMES+=("$elapsed"); STEP_STATUS+=("$status")
        _tee "  [PASS] ${name}  ($(_fmt_elapsed $elapsed))"
        _sync_staging_log   # keep the run-dir log copy current after each completed step
    elif (( exit_code == WARN_EXIT_CODE )); then
        status="WARN"
        STEP_NAMES+=("$name"); STEP_TIMES+=("$elapsed"); STEP_STATUS+=("$status")
        _tee "  [WARN] ${name}  ($(_fmt_elapsed $elapsed)) — completed with warnings (exit ${exit_code}); pipeline continues."
        _sync_staging_log
    else
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
    local line_c1="$(printf '─%.0s' {1..60})"
    local line_c2="$(printf '─%.0s' {1..12})"
    local line_c3="$(printf '─%.0s' {1..10})"

    _tee ""
    _tee "  TIMING SUMMARY"
    _tee "  ┌${line_c1}┬${line_c2}┬${line_c3}┐"
    _tee "  │ $(printf '%-58s │ %10s │ %-8s' 'STEP' 'ELAPSED' 'STATUS') │"
    _tee "  ├${line_c1}┼${line_c2}┼${line_c3}┤"
    for i in "${!STEP_NAMES[@]}"; do
        local t="${STEP_TIMES[$i]}"
        local st="${STEP_STATUS[$i]}"
        local name="${STEP_NAMES[$i]}"
        name="${name//—/-}"
        name="${name//–/-}"
        # Cap to the column width (58) so a long step label cannot overflow and
        # push the box borders out of alignment.
        (( ${#name} > 58 )) && name="${name:0:55}..."
        _tee "  │ $(printf '%-58s │ %10s │ %-8s' "$name" "$(_fmt_elapsed $t)" "$st") │"
        total=$(( total + t ))
    done
    _tee "  ├${line_c1}┼${line_c2}┼${line_c3}┤"
    _tee "  │ $(printf '%-58s │ %10s │ %-8s' 'TOTAL WALL TIME' "$(_fmt_elapsed $total)" "") │"
    _tee "  └${line_c1}┴${line_c2}┴${line_c3}┘"
    # Legend only when a step is anything other than a plain PASS, so a clean run
    # stays uncluttered.
    if printf '%s\n' "${STEP_STATUS[@]}" | grep -qvx 'PASS'; then
        _tee "    PASS = completed · WARN = completed, a complementary part deferred/failed (e.g. MM-GBSA) · SKIP = --resume-from · FAIL = halted"
    fi
}

# ── Conda activation ──────────────────────────────────────────────────────────
# If the caller already has the PFAS environment active, trust it and skip the
# activation dance entirely — this avoids aborting on hosts where conda is not a
# shell function or lives in a non-standard directory.
if [[ "${CONDA_DEFAULT_ENV:-}" == "PFAS" ]]; then
    _tee "  PFAS conda environment already active — activation bypassed."
elif true; then

if [[ -z "${CONDA_BASE:-}" ]]; then
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
fi
if [[ -z "${CONDA_BASE:-}" ]]; then
    CONDA_BASE="$HOME/miniconda3"
fi

if [[ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate PFAS || { _tee "ERROR: failed to activate conda env 'PFAS'. Aborting (would otherwise fall back to system Python without gemmi/rdkit/boltz)."; exit 1; }
    if [[ "${CONDA_DEFAULT_ENV:-}" != "PFAS" ]]; then
        _tee "ERROR: 'PFAS' env is not active after activation (got '${CONDA_DEFAULT_ENV:-none}'). Aborting."
        exit 1
    fi
else
    _tee "ERROR: conda not found at ${CONDA_BASE}. Set CONDA_BASE env var or install miniconda. Aborting."
    exit 1
fi

fi   # end conda activation (bypassed when PFAS already active)

# ── Pipeline body ─────────────────────────────────────────────────────────────
# Wrapped in a function so it can run either in the foreground (output tee'd to
# terminal + log) or detached in the background (output to log only). The trap is
# registered here so cleanup fires in whichever context actually runs the body.
_run_all_steps() {
trap '_sync_staging_log
      kill "$_SUDO_KEEPALIVE_PID" 2>/dev/null || true
      if [[ "$SUDO_ENABLED" == "1" ]]; then
        sudo systemctl unmask systemd-oomd.socket 2>/dev/null || true
        sudo systemctl start systemd-oomd 2>/dev/null || true
      fi
      true' EXIT

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
    _tee "  FASTA   : ${_FRESH_FASTA}  (produced by Step 01)"
    _tee "  SMI     : ${_FRESH_SMI}"
    _tee "  Run ID  : <assigned by Step 02>"
fi
_tee "  Log     : ${LOG_FILE}"
_tee "$_sep"
echo ""
echo "  Monitor: tail -f ${LOG_FILE}"
echo ""

# ── Pipeline steps ────────────────────────────────────────────────────────────

# Verify the environment is complete.
run_step "00a  Environment check" \
    python 00_03_Environment_FAcDs.py

# Refresh the canonical root PFAS.yml + requirements.txt every run (current host versions,
# export timestamp in the header) so they are always present and up to date.
run_step "00b  Environment export" \
    python 00_03_Environment_FAcDs.py --export || true

run_step "01  Merge sequences" \
    python 01_Merge_FAcDs.py \
        --master    A_Labelled_15-Seq.fasta \
        --secondary B_Downloaded-Blast_Uniprot_NCBI.fasta \
        --output    "${_FRESH_FASTA}"

# ── Step 02: Production (Boltz-2 co-folding) ──────────────────────────────────
if [[ "$_PIPELINE_MODE" == "resume" ]]; then
    run_step "02  Production (Boltz-2 scoring — resume)" \
        python 02_Production_FAcDs.py --resume "$RUN_ID"
else
    run_step "02  Production (Boltz-2 scoring — fresh)" \
        python 02_Production_FAcDs.py \
            --fasta "${_FRESH_FASTA}" \
            --smi   "${_FRESH_SMI}"

    # Auto-detect the run directory just created by Step 02 — pick the most-recently
    # MODIFIED Boltz-2_* dir (not the lexically last), so a pre-existing dir that sorts
    # later by name is not mistaken for the one Step 02 just made. Still best-effort:
    # a concurrent run in the same folder can race this.
    _new_run=$(find "${SCRIPT_DIR}" -maxdepth 1 -type d -name 'Boltz-2_*' -printf '%T@ %p\n' 2>/dev/null \
               | sort -n | tail -1 | cut -d' ' -f2-)
    if [[ -z "$_new_run" ]]; then
        _tee ""
        _tee "  ERROR: Step 02 did not create a Boltz-2_* run directory."
        _tee "         Check the Step 02 log above for errors."
        exit 1
    fi
    RUN_ID="$(basename "$_new_run")"
    _tee "  Fresh run directory: ${RUN_ID}"

    # Move / copy the staging log into the real run directory
    _real_log_dir="${SCRIPT_DIR}/${RUN_ID}/0_FAcDs_Pipeline_Logs"
    mkdir -p "$_real_log_dir"
    _STAGING_SYNC_DEST="${_real_log_dir}/$(basename "$LOG_FILE")"
    # The EXIT trap (registered above) already calls _sync_staging_log on exit — do
    # NOT register a second EXIT trap here, as it would override the oomd/keepalive
    # cleanup. Per-step re-syncs run via run_step; sync once now for the run dir.
    _sync_staging_log
    _tee "  Pipeline log copied to: ${_real_log_dir}/ (re-synced after each step)"
fi

# ── Steps 03–07: downstream analysis (all modes use RUN_ID) ──────────────────

run_step "03  Validation figures" \
    python 03_Validation_Figures_FAcDs.py "$RUN_ID"

run_step "04  Dendrogram" \
    python 04_Dendrogram_FAcDs.py "$RUN_ID"

run_step "05  Top-N selection + CIF/PDB generation & preparation (MD-ready cohort)" \
    python 05_TopN_and_PDB_Preparation_FAcDs.py "$RUN_ID"

# Mask systemd-oomd before Steps 06 (SID) and 07 (MD), the OOM-prone phase.
# The SID script is called with --pipeline-mode so it does NOT unmask on exit,
# leaving the mask in place for the subsequent MD thermodynamics engine.
# Guard on step 07 (the higher step number): mask whenever the MD engine runs.
if [[ $DRY_RUN -eq 0 && 7 -ge $RESUME_FROM ]]; then
    if [[ "$SUDO_ENABLED" == "1" ]]; then
        sudo systemctl stop systemd-oomd 2>/dev/null || true
        sudo systemctl mask systemd-oomd.socket 2>/dev/null || true
    else
        _tee "  [NORMAL MODE] systemd-oomd not masked (no sudo) — Steps 06-07 run unprotected from the OOM-killer."
    fi
fi

run_step "06  SID + Prime MM-GBSA post-processing (*_SID-out.eaf + MM-GBSA)" \
    python 06_SID_Prime-MMGBSA_FAcDs.py "$RUN_ID" --pipeline-mode

run_step "07  MD thermodynamics + QM/MM engine" \
    python 07_MD_QMMM_Defluorination_FAcDs.py "$RUN_ID"

# Restore systemd-oomd after both Steps 06 and 07 have completed (only if masked).
if [[ $DRY_RUN -eq 0 && "$SUDO_ENABLED" == "1" && 7 -ge $RESUME_FROM ]]; then
    sudo systemctl unmask systemd-oomd.socket && sudo systemctl start systemd-oomd
fi

# ── Footer ────────────────────────────────────────────────────────────────────

_tee ""
_tee "$_sep"
_tee "  ALL STEPS COMPLETE"
_tee "  Finished: $(date '+%Y-%m-%d %H:%M:%S')"
_tee "  Run ID  : ${RUN_ID}"
_tee "  Full log: ${LOG_FILE}"

_print_timing_table
}

# ── Dispatch: foreground (tee to terminal + log) or background (detached) ──────
_ANSI_STRIP='s/\x1B\[[0-9;]*[mKABCDEFGHJKSTfhilmnprsu]//g'
if [[ "$_BG_MODE" -eq 1 ]]; then
    # Detach: log-only output, stdin from /dev/null, removed from the job table
    # so closing the terminal (SIGHUP) does not kill the run. 'set -m' puts the
    # subshell in its own process group (PGID == PID) so every step + child can
    # be killed together with a single negative-PID signal.
    set -m
    ( _run_all_steps ) > >(sed -u "$_ANSI_STRIP" >> "$LOG_FILE") 2>&1 </dev/null &
    _BG_PID=$!
    set +m
    disown 2>/dev/null || true
    # The original sudo keepalive monitors the parent shell ($$), which exits on
    # detach below — re-anchor it to the detached run so credentials stay fresh for
    # the whole background run (only relevant when sudo was enabled).
    if [[ "$SUDO_ENABLED" == "1" ]]; then
        kill "$_SUDO_KEEPALIVE_PID" 2>/dev/null || true
        ( while kill -0 "$_BG_PID" 2>/dev/null; do sudo -vn 2>/dev/null; sleep 60; done ) &
        disown 2>/dev/null || true
    fi
    echo ""
    echo "  Pipeline detached — PID ${_BG_PID}. Safe to close this terminal."
    echo "  Monitor:  tail -f ${LOG_FILE}"
    echo "  Stop it :  kill -- -${_BG_PID}      # kills the run and its current step"
    echo ""
    exit 0
else
    # Foreground: original behaviour — one write to log (ANSI-stripped), one to terminal.
    exec > >(tee >(sed -u "$_ANSI_STRIP" >> "$LOG_FILE")) 2>&1
    _run_all_steps
fi
