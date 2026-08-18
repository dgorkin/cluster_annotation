#!/usr/bin/env bash
# Cluster Annotation Browser launcher.
#
# Binds to localhost only; you reach it from your laptop over an SSH tunnel. `start` runs the app
# detached so it survives logging out — a 34-cluster generation takes tens of minutes and used to
# die with the SSH session.
#
#   ./run_app.sh [start] [PORT]   start in the background; prints the ssh tunnel command
#   ./run_app.sh fg [PORT]        run in the foreground (Ctrl-C to quit)
#   ./run_app.sh status           running? which port? with the tunnel command to paste
#   ./run_app.sh stop             stop the background instance
#   ./run_app.sh restart [PORT]
#   ./run_app.sh logs             follow the app log
#   ./run_app.sh doctor           preflight: env, Rscript, config inputs, secrets permissions
#   ./run_app.sh counts [CONFIG]  count cells per cluster from the Seurat object (--force to redo)
#
# PORT is a starting point, not a requirement: if it is taken (a stale instance, another user on
# this shared box) the next free port is used and reported. Default 8501.
set -euo pipefail
cd "$(dirname "$0")"

# Runs from the dedicated conda env (override with CLUSTER_ANNOTATION_ENV=/path/to/env).
# Uses the env's interpreter directly so no `conda activate` is needed.
ENV_DIR="${CLUSTER_ANNOTATION_ENV:-$HOME/.conda/envs/cluster_annotation}"
PY="$ENV_DIR/bin/python"
STREAMLIT="$ENV_DIR/bin/streamlit"

RUN_DIR=".run"
PID_FILE="$RUN_DIR/app.pid"
PORT_FILE="$RUN_DIR/app.port"
LOG_FILE="logs/app.log"
DEFAULT_PORT=8501

die() { printf 'error: %s\n' "$*" >&2; exit 1; }

require_env() {
    [[ -x "$STREAMLIT" ]] || die "streamlit not found at $STREAMLIT
  The app runs from a dedicated conda env. Create it with:
    conda create -y -n cluster_annotation python=3.11
    $ENV_DIR/bin/pip install -r requirements.txt
  or point CLUSTER_ANNOTATION_ENV at an existing env. See ./run_app.sh doctor."
}

# First free TCP port at or above $1, so a busy 8501 is routed around instead of fatal.
pick_port() {
    "$PY" - "$1" <<'PYEOF'
import socket, sys
start = int(sys.argv[1])
for port in range(start, start + 50):
    s = socket.socket()
    try:
        s.bind(("127.0.0.1", port))
    except OSError:
        continue
    finally:
        s.close()
    print(port)
    break
else:
    sys.exit(f"no free port in {start}-{start + 49}")
PYEOF
}

app_pid() {  # echo the live pid, or nothing
    [[ -f "$PID_FILE" ]] || return 0
    local pid
    pid="$(cat "$PID_FILE" 2>/dev/null || true)"
    [[ -n "$pid" ]] || return 0
    # Confirm it is still OUR process: a recycled pid belonging to something else must not be
    # reported as the app, and must never be killed by `stop`.
    if kill -0 "$pid" 2>/dev/null && grep -qa streamlit "/proc/$pid/cmdline" 2>/dev/null; then
        printf '%s' "$pid"
    fi
}

tunnel_hint() {
    local port="$1" host
    host="$(hostname -f 2>/dev/null || hostname)"
    cat <<EOF

  On your laptop, open the tunnel then browse:
    ssh -N -L ${port}:localhost:${port} ${USER}@${host}
    http://localhost:${port}
EOF
}

cmd_start() {
    require_env
    local pid; pid="$(app_pid)"
    if [[ -n "$pid" ]]; then
        echo "already running (pid $pid, port $(cat "$PORT_FILE" 2>/dev/null || echo '?'))"
        tunnel_hint "$(cat "$PORT_FILE" 2>/dev/null || echo "$DEFAULT_PORT")"
        echo "  (./run_app.sh restart to bounce it)"
        return 0
    fi
    mkdir -p "$RUN_DIR" logs
    local port; port="$(pick_port "${1:-$DEFAULT_PORT}")"
    setsid nohup "$STREAMLIT" run app/app.py \
        --server.address 127.0.0.1 \
        --server.port "$port" \
        --server.headless true \
        --browser.gatherUsageStats false \
        >>"$LOG_FILE" 2>&1 &
    local new_pid=$!
    echo "$new_pid" > "$PID_FILE"
    echo "$port" > "$PORT_FILE"
    # Give it a moment to bind, and surface an immediate crash rather than a bare "started".
    "$PY" - "$new_pid" "$port" "$LOG_FILE" <<'PYEOF'
import socket, sys, time
pid, port, log = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
import os
for _ in range(60):
    try:
        os.kill(pid, 0)
    except OSError:
        sys.exit(f"process exited during startup — see {log}")
    s = socket.socket()
    if s.connect_ex(("127.0.0.1", port)) == 0:
        s.close()
        sys.exit(0)
    s.close()
    time.sleep(0.25)
sys.exit(f"did not start listening on {port} within 15s — see {log}")
PYEOF
    echo "started (pid $new_pid, port $port), logging to $LOG_FILE"
    tunnel_hint "$port"
}

cmd_stop() {
    local pid; pid="$(app_pid)"
    if [[ -z "$pid" ]]; then
        echo "not running"
        rm -f "$PID_FILE" "$PORT_FILE"
        return 0
    fi
    kill "$pid"
    for _ in $(seq 40); do
        kill -0 "$pid" 2>/dev/null || break
        sleep 0.25
    done
    if kill -0 "$pid" 2>/dev/null; then
        echo "did not exit on TERM, sending KILL"
        kill -9 "$pid" 2>/dev/null || true
    fi
    rm -f "$PID_FILE" "$PORT_FILE"
    echo "stopped (was pid $pid)"
}

cmd_status() {
    local pid; pid="$(app_pid)"
    if [[ -z "$pid" ]]; then
        echo "not running"
        [[ -f "$PID_FILE" ]] && echo "  (stale pidfile $PID_FILE — ./run_app.sh stop clears it)"
        return 1
    fi
    local port; port="$(cat "$PORT_FILE" 2>/dev/null || echo '?')"
    echo "running: pid $pid, port $port, up $(ps -o etime= -p "$pid" | tr -d ' ')"
    tunnel_hint "$port"
}

case "${1:-start}" in
    start)          shift || true; cmd_start "${1:-$DEFAULT_PORT}" ;;
    [0-9]*)         cmd_start "$1" ;;                 # backwards compatible: ./run_app.sh 8502
    fg|foreground)  shift || true
                    require_env
                    port="$(pick_port "${1:-$DEFAULT_PORT}")"
                    tunnel_hint "$port"
                    exec "$STREAMLIT" run app/app.py \
                        --server.address 127.0.0.1 --server.port "$port" \
                        --server.headless true --browser.gatherUsageStats false ;;
    stop)           cmd_stop ;;
    restart)        shift || true; cmd_stop; cmd_start "${1:-$DEFAULT_PORT}" ;;
    status)         cmd_status ;;
    logs)           tail -f "$LOG_FILE" ;;
    doctor)         require_env; exec "$PY" scripts/doctor.py "${@:2}" ;;
    counts)         require_env; exec "$PY" scripts/cell_counts.py "${@:2}" ;;
    # Print the header comment block as the help text, stopping at the first non-comment line so
    # it cannot drift out of sync when the header grows.
    -h|--help|help) awk 'NR>1 && /^#/ {sub(/^# ?/, ""); print; next} NR>1 {exit}' "$0" ;;
    *)              die "unknown command '$1' (try: ./run_app.sh --help)" ;;
esac
