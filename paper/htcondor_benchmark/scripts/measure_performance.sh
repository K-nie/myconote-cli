#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# measure_performance.sh
# Wrap an arbitrary command with `/usr/bin/time -v` and convert its output
# into a performance.json file that aggregate_metrics.py can consume.
#
# Usage:
#     bash measure_performance.sh STAGE_NAME OUT_DIR -- <command...>
#
# Example:
#     bash measure_performance.sh predict /tmp/out -- myconote-cli predict \
#         --genome genome.fa --out /tmp/out
#
# Produces:
#     OUT_DIR/time_STAGE_NAME.log    (raw `time -v` output)
#     OUT_DIR/perf_STAGE_NAME.json   (parsed metrics)
#
# Downstream, run_*.sh wrappers concatenate per-stage perf_*.json files into
# a single performance.json that aggregate_metrics.py reads.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

if [[ "$#" -lt 4 ]]; then
    echo "Usage: $0 STAGE_NAME OUT_DIR -- <command...>" >&2
    exit 2
fi

STAGE="$1"; shift
OUT_DIR="$1"; shift
SEP="$1"; shift
if [[ "$SEP" != "--" ]]; then
    echo "ERROR: expected '--' before command, got '$SEP'" >&2
    exit 2
fi

mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/time_${STAGE}.log"
JSON="$OUT_DIR/perf_${STAGE}.json"

# ── Pick the GNU `time` binary. BSD time on macOS does not support -v. ────
TIME_BIN=""
if [[ -x /usr/bin/time ]]; then
    if /usr/bin/time -v true 2>/dev/null; then
        TIME_BIN="/usr/bin/time -v"
    fi
fi
if [[ -z "$TIME_BIN" ]] && command -v gtime &>/dev/null; then
    TIME_BIN="gtime -v"
fi
if [[ -z "$TIME_BIN" ]]; then
    echo "ERROR: GNU time (-v) not found. Install with 'brew install gnu-time' on macOS." >&2
    exit 3
fi

# ── Run the command under time -v ────────────────────────────────────────
set +e
$TIME_BIN -o "$LOG" "$@"
EXIT_CODE=$?
set -e

# ── Parse time -v output into JSON ───────────────────────────────────────
python3 - "$LOG" "$JSON" "$STAGE" "$EXIT_CODE" << 'PYEOF'
import json
import re
import sys

log_path, json_path, stage, exit_code = sys.argv[1:5]

data = {
    'stage': stage,
    'exit_code': int(exit_code),
    'wall_seconds': 0.0,
    'user_seconds': 0.0,
    'sys_seconds': 0.0,
    'cpu_percent': 0,
    'peak_rss_kb': 0,
    'peak_rss_mb': 0.0,
    'major_page_faults': 0,
    'minor_page_faults': 0,
    'fs_inputs': 0,
    'fs_outputs': 0,
    'voluntary_ctx_switches': 0,
    'involuntary_ctx_switches': 0,
}

try:
    with open(log_path) as f:
        for line in f:
            line = line.strip()
            if 'Elapsed (wall clock) time' in line:
                m = re.search(r'(\d+):?(\d+)?:(\d+(?:\.\d+)?)', line)
                if m:
                    h = int(m.group(1)) if m.group(2) else 0
                    mm = int(m.group(2)) if m.group(2) else int(m.group(1))
                    ss = float(m.group(3))
                    data['wall_seconds'] = round(h * 3600 + mm * 60 + ss, 2)
            elif line.startswith('User time (seconds):'):
                data['user_seconds'] = float(line.split(':')[1].strip())
            elif line.startswith('System time (seconds):'):
                data['sys_seconds'] = float(line.split(':')[1].strip())
            elif line.startswith('Percent of CPU this job got:'):
                v = line.split(':')[1].strip().rstrip('%')
                try:
                    data['cpu_percent'] = int(v)
                except ValueError:
                    pass
            elif line.startswith('Maximum resident set size (kbytes):'):
                kb = int(line.split(':')[1].strip())
                data['peak_rss_kb'] = kb
                data['peak_rss_mb'] = round(kb / 1024, 1)
            elif line.startswith('Major (requiring I/O) page faults:'):
                data['major_page_faults'] = int(line.split(':')[1].strip())
            elif line.startswith('Minor (reclaiming a frame) page faults:'):
                data['minor_page_faults'] = int(line.split(':')[1].strip())
            elif line.startswith('File system inputs:'):
                data['fs_inputs'] = int(line.split(':')[1].strip())
            elif line.startswith('File system outputs:'):
                data['fs_outputs'] = int(line.split(':')[1].strip())
            elif line.startswith('Voluntary context switches:'):
                data['voluntary_ctx_switches'] = int(line.split(':')[1].strip())
            elif line.startswith('Involuntary context switches:'):
                data['involuntary_ctx_switches'] = int(line.split(':')[1].strip())
except FileNotFoundError:
    pass

with open(json_path, 'w') as f:
    json.dump(data, f, indent=2)
PYEOF

echo "  [perf] stage=$STAGE  log=$LOG  json=$JSON  exit=$EXIT_CODE"
exit $EXIT_CODE
