#!/usr/bin/env bash
# ==============================================================================
#  disk_guard_FAcDs.sh — protect running Desmond MD from a scratch-disk blow-out
# ==============================================================================
#
#  WHY THIS EXISTS
#  ---------------
#  Prime MM-GBSA (and QSite QM/MM) write their per-subjob scratch through the
#  Schrodinger job-server daemon. When that daemon was started with
#  `-tmpdir /tmp`, every subjob lands on the OS disk, not the WD working disk.
#  A single MM-GBSA can spill up to ~600 GB. If it fills the disk on the last
#  frame it not only crashes itself but ALSO kills any co-resident Desmond MD
#  (which needs a little /tmp to checkpoint) — losing days of trajectory.
#
#  This guard treats free disk space as the thing to protect. MM-GBSA and QSite
#  are idempotent (06/07 re-run them from scratch on the next pass), so when the
#  disk is about to fill it sacrifices THOSE jobs and leaves Desmond untouched.
#
#  WHAT IT DOES
#  ------------
#    every INTERVAL seconds, read free bytes on the scratch filesystem:
#      free > WARN_GB   → silent (healthy)
#      free <= WARN_GB  → log + desktop notce (no action) — you have time to react
#      free <= CRIT_GB  → EMERGENCY: kill only the MM-GBSA / QSite job trees,
#                         never Desmond, then keep watching so a second hog can
#                         also be caught. This stops disk growth and saves the MD.
#
#  It never touches: desmond, gdesmond, mdsim, multisim, cmj_startup (Desmond MD),
#  nor the jobserver daemon itself (killing that orphans everything).
#
#  USAGE
#  -----
#    ./disk_guard_FAcDs.sh                 # defaults, foreground (Ctrl-C to stop)
#    nohup ./disk_guard_FAcDs.sh &         # run detached, survives logout
#    WARN_GB=150 CRIT_GB=90 ./disk_guard_FAcDs.sh
#    DRY_RUN=1 ./disk_guard_FAcDs.sh       # log what it WOULD kill, kill nothing
#    SCRATCH_PATH=/tmp ./disk_guard_FAcDs.sh
#
#  All output is logged to the WD disk (GUARD_LOG), never to the disk it guards.
# ==============================================================================
set -u

# --- configuration (override via environment) --------------------------------
SCRATCH_PATH="${SCRATCH_PATH:-/tmp}"                       # filesystem to guard
WARN_GB="${WARN_GB:-120}"                                  # notify below this
CRIT_GB="${CRIT_GB:-70}"                                   # kill hogs below this
INTERVAL="${INTERVAL:-20}"                                 # poll seconds
GUARD_LOG="${GUARD_LOG:-/mnt/wdpassport/FAcDs/disk_guard.log}"
DRY_RUN="${DRY_RUN:-0}"                                    # 1 = never kill

# Process patterns to SACRIFICE when critical (restartable Schrodinger scoring).
# thermal_mmgbsa runs on a *-out.cms whose name contains "desmond_md_job..."; match
# the SCRIPT name, never the data filename, so PROTECT below can't shield it.
KILL_PATTERNS='thermal_mmgbsa|run_qsite|run_jaguar_backend|jexec|prime_mmgbsa|/prime\b'
# Process patterns to ALWAYS PROTECT (never kill, even if matched above). These are
# Desmond MD engine binaries + the job-server daemon — matched by EXECUTABLE name,
# not by the bare word "desmond" (which appears inside MM-GBSA .cms filenames).
PROTECT_PATTERNS='desmond_driver|gdesmond|/mdsim|cmj_startup|chorus_multijob|multisim|jobserverd|job_supervisor'

# --- single-instance lock ----------------------------------------------------
LOCK="/mnt/wdpassport/FAcDs/.disk_guard.lock"
if [ -e "$LOCK" ] && kill -0 "$(cat "$LOCK" 2>/dev/null)" 2>/dev/null; then
    echo "disk_guard already running (PID $(cat "$LOCK")). Exiting." >&2
    exit 1
fi
echo $$ > "$LOCK"
trap 'rm -f "$LOCK"' EXIT INT TERM

log() { printf '%s  %s\n' "$(date '+%F %T')" "$*" | tee -a "$GUARD_LOG"; }

notify() {
    # best-effort desktop toast; silent if headless
    command -v notify-send >/dev/null 2>&1 && \
        DISPLAY="${DISPLAY:-:0}" notify-send -u critical "FAcDs disk guard" "$1" 2>/dev/null || true
}

free_gb() {
    # integer GiB free on the guarded filesystem
    df -PB1G "$SCRATCH_PATH" 2>/dev/null | awk 'NR==2 {print $4}'
}

# Kill the restartable scoring job trees; never Desmond. Returns count killed.
sacrifice_hogs() {
    local killed=0 pid args
    while read -r pid args; do
        [ -z "$pid" ] && continue
        # skip anything protected
        echo "$args" | grep -Eq "$PROTECT_PATTERNS" && continue
        echo "$args" | grep -Eq "$KILL_PATTERNS"    || continue
        if [ "$DRY_RUN" = "1" ]; then
            log "  [DRY_RUN] would kill PID $pid : ${args:0:90}"
        else
            log "  KILL PID $pid : ${args:0:90}"
            kill -TERM "$pid" 2>/dev/null
        fi
        killed=$((killed+1))
    done < <(ps -eo pid=,args= | grep -Ev "$PROTECT_PATTERNS" | grep -E "$KILL_PATTERNS")
    # give TERM a moment, then hard-kill stragglers (still never Desmond)
    if [ "$DRY_RUN" != "1" ] && [ "$killed" -gt 0 ]; then
        sleep 5
        while read -r pid args; do
            [ -z "$pid" ] && continue
            echo "$args" | grep -Eq "$PROTECT_PATTERNS" && continue
            echo "$args" | grep -Eq "$KILL_PATTERNS"    || continue
            log "  KILL -9 PID $pid (did not exit on TERM)"
            kill -9 "$pid" 2>/dev/null
        done < <(ps -eo pid=,args= | grep -Ev "$PROTECT_PATTERNS" | grep -E "$KILL_PATTERNS")
    fi
    return 0
}

# --- main loop ---------------------------------------------------------------
log "=== disk_guard started : guard=$SCRATCH_PATH  WARN=${WARN_GB}G  CRIT=${CRIT_GB}G  interval=${INTERVAL}s  dry_run=$DRY_RUN ==="
armed_warn=0     # so WARN logs once per crossing, not every tick
while true; do
    fg="$(free_gb)"
    if [ -z "$fg" ]; then
        log "WARN: cannot read free space on $SCRATCH_PATH — retrying"
        sleep "$INTERVAL"; continue
    fi

    if [ "$fg" -le "$CRIT_GB" ]; then
        log "CRITICAL: only ${fg}G free on $SCRATCH_PATH (<= ${CRIT_GB}G). Sacrificing MM-GBSA/QSite to save Desmond."
        notify "CRITICAL ${fg}G free — killing MM-GBSA/QSite to protect Desmond MD."
        sacrifice_hogs
        log "  post-kill free: $(free_gb)G. Continuing to watch."
        armed_warn=0
    elif [ "$fg" -le "$WARN_GB" ]; then
        if [ "$armed_warn" -eq 0 ]; then
            log "WARN: ${fg}G free on $SCRATCH_PATH (<= ${WARN_GB}G). No action yet; will kill hogs at ${CRIT_GB}G."
            notify "Low disk: ${fg}G free. Will protect Desmond at ${CRIT_GB}G."
            armed_warn=1
        fi
    else
        armed_warn=0
    fi
    sleep "$INTERVAL"
done
