#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

make Release >/dev/null

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

threads_list="${THREADS_LIST:-1 2 4}"
n_list="${N_LIST:-50 512 1024 2048}"

printf '%-8s %-8s %-12s\n' "N" "threads" "calc_seconds"
for n in $n_list; do
    for threads in $threads_list; do
        cfg="$tmpdir/config-N${n}-T${threads}.yaml"
        cat >"$cfg" <<EOF
TimeIntegration: RungeKuttaClassic
EnablePlotting: false
OpenMPMinCells: 0
OpenMPThreads: $threads
OpenMPLargeMinCells: 0
OpenMPLargeThreads: $threads
CFLReserve: 0.45
N: $n
order: 5
t0: 0.0
t1: 2.0e-4
x0: 0.0
x1: 0.5
c_0: 300.0
rho_0: 1.5
freq_i: 3000.0
mu: 0.1
sigma: 0.02
height: 1.0
EOF
        output="$(cd "$tmpdir" && OPENBLAS_NUM_THREADS=1 "$repo_root/bin/Release/DiscontinousGalerkin" "$cfg" </dev/null)"
        seconds="$(printf '%s\n' "$output" | awk '/^elapsed time:/ {print $3; exit}')"
        printf '%-8s %-8s %-12s\n' "$n" "$threads" "${seconds:-n/a}"
    done
done
