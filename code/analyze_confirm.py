# Phayung Meesad(1), Van Van Chau(1,2,*), and Tuan Minh Nguyen(2)
# (1) Department of Information Technology Management, King Mongkut's University of Technology North Bangkok,
#     Bangkok, Thailand, phayung.m@itd.kmutnb.ac.th
# (2) Faculty of Information Technology, Posts and Telecommunications Institute of Technology,
#     Ho Chi Minh City, Vietnam, vancv@ptit.edu.vn, minhtuan@ptit.edu.vn
# (*) Corresponding author
#
"""Locked analysis of the confirmation run: frac=1.0 against the default on fresh seeds.
Usage: python3 analyze_confirm.py RESULTS_DIR OUT_DIR   (expects RESULTS_DIR/conf-s0, conf-s1 with r1_confirm_s*.csv)"""
import glob, os, sys
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
EPS = 1e-8
d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(sys.argv[1], "conf-s*", "r1_confirm_*.csv")))], ignore_index=True)
assert (d.exp == "confirm").all(), "error rows present"
assert len(d) == 1440 and not d.duplicated(["suite", "problem", "config", "run"]).any(), len(d)
assert d.seed.min() >= 500000
def holm(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); a = np.empty(m); run = 0
    for i, k in enumerate(o): run = max(run, (m - i) * p[k]); a[k] = min(1, run)
    return a
c = d[d.suite != "eng"].copy(); c["le"] = np.log10(c.error.astype(float).clip(lower=0) + EPS)
c["fn"] = [("F" if s == "cec2017" else "G") + str(p) for s, p in zip(c.suite, c.problem)]
M = c.groupby(["fn", "config"])["le"].mean().unstack()
better = int((M["frac=1.0"] < M["default"]).sum()); worse = int((M["frac=1.0"] > M["default"]).sum())
p = wilcoxon(M["frac=1.0"], M["default"], zero_method="zsplit").pvalue
print(f"PRIMARY: frac=1.0 lower mean log error on {better}/20, higher on {worse}/20, Wilcoxon p={p:.4g}")
print("CONFIRMED" if (p < 0.05 and better > worse) else "NOT CONFIRMED")
rows = []
for f, g in c.groupby("fn"):
    a = g[g.config == "frac=1.0"].sort_values("run")["le"].values; b = g[g.config == "default"].sort_values("run")["le"].values
    rows.append([f, np.median(b - a), wilcoxon(a, b, zero_method="zsplit").pvalue if np.any(a != b) else 1.0])
r = pd.DataFrame(rows, columns=["fn", "median_gain", "p"]); r["p_holm"] = holm(r.p)
r["res"] = np.where(r.p_holm >= .05, "tie", np.where(r.median_gain > 0, "win", "loss"))
print("Run-level (frac=1.0 vs default):", r.res.value_counts().to_dict())
e = d[d.suite == "eng"]; er = []
for pr, g in e.groupby("problem"):
    a = g[g.config == "frac=1.0"].sort_values("run").best.values; b = g[g.config == "default"].sort_values("run").best.values
    er.append([pr, b.mean(), a.mean(), wilcoxon(a, b, zero_method="zsplit").pvalue])
er = pd.DataFrame(er, columns=["problem", "default_mean", "frac1_mean", "p"]); er["p_holm"] = holm(er.p)
print(M.round(4).to_string()); print(r.round(4).to_string(index=False)); print(er.to_string(index=False))
OUT = sys.argv[2] if len(sys.argv) > 2 else "analysis"; os.makedirs(OUT, exist_ok=True)
M.to_csv(f"{OUT}/confirm_meanlog.csv"); r.to_csv(f"{OUT}/confirm_runlevel.csv", index=False); er.to_csv(f"{OUT}/confirm_eng.csv", index=False)
