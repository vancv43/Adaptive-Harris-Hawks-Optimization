#!/bin/bash
# Phayung Meesad(1), Van Van Chau(1,2,*), and Tuan Minh Nguyen(2,*)
# (1) Department of Information Technology Management, King Mongkut's University of Technology North Bangkok,
#     Bangkok, Thailand, phayung.m@itd.kmutnb.ac.th
# (2) Faculty of Information Technology, Posts and Telecommunications Institute of Technology,
#     Ho Chi Minh City, Vietnam, vancv@ptit.edu.vn, minhtuan@ptit.edu.vn
# (*) Corresponding author
#
# Reruns every optimization of the paper, then the analysis. 22,320 runs in total, several hundred
# core-hours; each experiment is split into seed shards that can run on separate machines.
set -e
cd "$(dirname "$0")"
PY=${PYTHON:-python}
W=${WORKERS:-4}
$PY test_r1.py
for k in 0 1 2 3; do mkdir -p ../runs/main-s$k; $PY run_r1.py main --shard $k --nshards 4 --workers $W --out ../runs/main-s$k; done
mkdir -p ../runs/abl-col
$PY run_r1.py ablation --workers $W --out ../runs/abl-col
$PY run_r1.py collapse --workers $W --out ../runs/abl-col
for k in 0 1 2 3; do mkdir -p ../runs/sens-s$k; $PY run_r1.py sensitivity --shard $k --nshards 4 --workers $W --out ../runs/sens-s$k; done
for k in 0 1; do mkdir -p ../runs/conf-s$k; $PY run_r1.py confirm --shard $k --nshards 2 --workers $W --out ../runs/conf-s$k; done
bash run_analysis.sh
