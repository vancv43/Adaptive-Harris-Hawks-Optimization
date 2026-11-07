#!/bin/bash
# Phayung Meesad(1), Van Van Chau(1,2,*), and Tuan Minh Nguyen(2)
# (1) Department of Information Technology Management, King Mongkut's University of Technology North Bangkok,
#     Bangkok, Thailand, phayung.m@itd.kmutnb.ac.th
# (2) Faculty of Information Technology, Posts and Telecommunications Institute of Technology,
#     Ho Chi Minh City, Vietnam, vancv@ptit.edu.vn, minhtuan@ptit.edu.vn
# (*) Corresponding author
#
# Recomputes every table of the paper from the saved runs in ../runs; no optimization is rerun.
set -e
cd "$(dirname "$0")"
PY=${PYTHON:-python}
$PY analyze_r1.py ../runs ../results/analysis
$PY make_r1_tables.py ../runs ../results/tables   # LaTeX tables, not kept in the repository
$PY best_feasible_all.py ../runs ../results/analysis
$PY analyze_confirm.py ../runs ../results/analysis > ../results/analysis/confirm_summary.txt
echo 'statistics in results/analysis, tables in results/tables'
