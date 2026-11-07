# Phayung Meesad(1), Van Van Chau(1,2,*), and Tuan Minh Nguyen(2)
# (1) Department of Information Technology Management, King Mongkut's University of Technology North Bangkok,
#     Bangkok, Thailand, phayung.m@itd.kmutnb.ac.th
# (2) Faculty of Information Technology, Posts and Telecommunications Institute of Technology,
#     Ho Chi Minh City, Vietnam, vancv@ptit.edu.vn, minhtuan@ptit.edu.vn
# (*) Corresponding author
#
"""Best feasible design of every optimizer on the engineering problems, Appendix B of the paper.

For each problem and optimizer, the run with the lowest objective among designs within the tolerance 1e-4 is taken.
The objective and the largest constraint violation are recomputed from its stored variables, which are
written with six significant digits. Usage: python3 best_feasible_all.py RESULTS_DIR OUT_DIR"""
import glob, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engineering_problems import PROBLEMS

RES, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
ORDER = ["GA", "PSO", "ACO", "ABC", "BAT", "GWO", "AOA", "SCSO", "HHO", "GCHHO", "m-HHO", "AHHO"]
d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(RES, "main-s*", "r1_main_*.csv")))], ignore_index=True)
e = d[d.suite == "eng"].copy()
e["bf"] = pd.to_numeric(e.best_feas, errors="coerce")
rows = []
for p in PROBLEMS:
    for c in ORDER:
        h = e[(e.problem == p) & (e.config == c) & np.isfinite(e.bf)]
        r = h.loc[h.bf.idxmin()]
        x = np.array([float(v) for v in r.best_feas_x.split(";")])
        obj, g = PROBLEMS[p]["raw"](x)
        inb = bool(np.all(x >= PROBLEMS[p]["lb"]) and np.all(x <= PROBLEMS[p]["ub"]))
        rows.append([p, c, obj, max(0.0, max(g)), ";".join(f"{v:.6g}" for v in x), len(h), inb, r.bf])
pd.DataFrame(rows, columns=["p", "c", "obj", "maxviol", "x", "nruns_feas", "inbounds", "logged"]).to_csv(
    f"{OUT}/best_feasible_all.csv", index=False)
