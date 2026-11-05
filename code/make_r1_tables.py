# Phayung Meesad(1), Van Van Chau(1,2,*), and Tuan Minh Nguyen(2,*)
# (1) Department of Information Technology Management, King Mongkut's University of Technology North Bangkok,
#     Bangkok, Thailand, phayung.m@itd.kmutnb.ac.th
# (2) Faculty of Information Technology, Posts and Telecommunications Institute of Technology,
#     Ho Chi Minh City, Vietnam, vancv@ptit.edu.vn, minhtuan@ptit.edu.vn
# (*) Corresponding author
#
"""LaTeX tables of the paper, computed from the per-run results.

Usage: python3 make_r1_tables.py RESULTS_DIR OUT_TEX_DIR
Statistics: mean log10(error+1e-8), Friedman ranks, Holm-corrected Wilcoxon.
"""
import glob, os, sys
import numpy as np, pandas as pd
from scipy.stats import friedmanchisquare, wilcoxon

RES, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
EPS = 1e-8
ORDER = ["GA", "PSO", "ACO", "ABC", "BAT", "GWO", "AOA", "SCSO", "HHO", "GCHHO", "m-HHO", "AHHO"]
BASE = ORDER[:-1]


def load(pat):
    return pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(RES, pat)))], ignore_index=True)


def holm(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); a = np.empty(m); r = 0
    for i, k in enumerate(o):
        r = max(r, (m - i) * p[k]); a[k] = min(1.0, r)
    return a


def pw(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return 1.0 if np.all(a == b) else wilcoxon(a, b, zero_method="zsplit").pvalue


def E(x):
    if pd.isna(x): return "--"
    if x == 0: return "0"
    s = f"{x:.3e}"; m, e = s.split("e"); return f"{m}E{int(e):+03d}"


def P(p):
    return f"{p:.3f}" if p >= 0.001 else "$<$0.001"


def lab(s, p):
    p = int(p)
    if s == "cec2017":
        return f"F{p}" if p < 2 else f"F{p} [F{p + 1}]"
    return f"F{p}"


main, abl, sens, col = load("main-s*/r1_main_*.csv"), load("abl-col/r1_ablation_*.csv"), \
    load("sens-s*/r1_sensitivity_*.csv"), load("abl-col/r1_collapse_*.csv")
for d in (main, abl, sens, col):
    d["key"] = [f"{s}:{p}" for s, p in zip(d.suite, d.problem)]
    num = d.suite != "eng"
    d.loc[num, "le"] = np.log10(d.loc[num, "error"].astype(float).clip(lower=0) + EPS)
probs = {s: sorted(main[main.suite == s].problem.astype(int).unique()) for s in ("cec2017", "cec2020")}
keys = {s: [f"{s}:{p}" for p in probs[s]] for s in probs}
allk = keys["cec2017"] + keys["cec2020"]
log = []

# ---------- descriptive longtables
for s, name in [("cec2017", "CEC2017"), ("cec2020", "CEC2020")]:
    L = [r"\begin{longtable}{llrrrrr}", rf"\caption{{Error statistics on {name} over 30 runs}}\label{{tab:{s}-desc}}\\",
         r"\toprule", r"Function & Optimizer & Min & Median & Mean & Max & STD \\", r"\midrule", r"\endfirsthead",
         r"\toprule", r"Function & Optimizer & Min & Median & Mean & Max & STD \\", r"\midrule", r"\endhead",
         r"\bottomrule", r"\endfoot"]
    for j, p in enumerate(probs[s]):
        g = main[(main.suite == s) & (main.problem.astype(str) == str(p))]
        for i, c in enumerate(ORDER):
            e = g[g.config == c].error.astype(float)
            L.append(f"{lab(s, p) if i == 0 else ''} & {c} & {E(e.min())} & {E(e.median())} & {E(e.mean())} & {E(e.max())} & {E(e.std())} \\\\")
        if j < len(probs[s]) - 1: L.append(r"\addlinespace")
    L.append(r"\end{longtable}")
    open(f"{OUT}/desc_{s}.tex", "w").write("\n".join(L) + "\n")

# ---------- Friedman ranks
cec = main[main.suite != "eng"]
M = cec.groupby(["key", "config"])["le"].mean().unstack()[ORDER]
MR = cec.groupby(["key", "config"]).error.mean().unstack()[ORDER]
rk = {s: M.loc[keys[s]].rank(axis=1).mean() for s in keys}; rk["all"] = M.loc[allk].rank(axis=1).mean()
rraw = MR.loc[allk].rank(axis=1).mean()
fr = {s: friedmanchisquare(*[M.loc[keys[s], c] for c in ORDER]) for s in keys}
fr["all"] = friedmanchisquare(*[M.loc[allk, c] for c in ORDER]); fr_raw = friedmanchisquare(*[MR.loc[allk, c] for c in ORDER])
eng = main[main.suite == "eng"]
EM = eng.groupby(["key", "config"]).best.mean().unstack()[ORDER]
rk_e = EM.rank(axis=1).mean(); fr_e = friedmanchisquare(*[EM[c] for c in ORDER])
L = [r"\begin{table}[htbp]", r"\caption{Average Friedman ranks, lower is better}\label{tab:ranks}", r"\centering",
     r"\begin{tabular}{lccccc}", r"\toprule",
     r"Optimizer & CEC2017 & CEC2020 & CEC Suites & CEC Suites, raw mean & Engineering \\", r"\midrule"]
for c in ORDER:
    L.append(f"{c} & {rk['cec2017'][c]:.2f} & {rk['cec2020'][c]:.2f} & {rk['all'][c]:.2f} & {rraw[c]:.2f} & {rk_e[c]:.2f} \\\\")
L += [r"\midrule", "Friedman $p$ & " + " & ".join(f"{x.pvalue:.1e}" for x in (fr["cec2017"], fr["cec2020"], fr["all"], fr_raw, fr_e)) + r" \\",
      r"\bottomrule", r"\end{tabular}", r"\end{table}"]
open(f"{OUT}/ranks.tex", "w").write("\n".join(L) + "\n")
log.append(("friedman", {k: (round(v.statistic, 2), v.pvalue) for k, v in fr.items()}, fr_raw.pvalue, fr_e.pvalue))


# ---------- AHHO versus each optimizer
def run_wtl(df, ks, col_):
    out = {c: [0, 0, 0] for c in BASE}
    for k in ks:
        g = df[df.key == k]; a = g[g.config == "AHHO"].sort_values("run")[col_].values
        ps, sg = [], []
        for c in BASE:
            b = g[g.config == c].sort_values("run")[col_].values
            ps.append(pw(a, b)); sg.append(np.sign(np.median(b - a)))
        for c, p, s_ in zip(BASE, holm(ps), sg):
            out[c][0 if (p < .05 and s_ > 0) else 2 if (p < .05 and s_ < 0) else 1] += 1
    return out


acr = {}
for s in keys:
    ps = [pw(M.loc[keys[s], "AHHO"], M.loc[keys[s], c]) for c in BASE]
    acr[s] = dict(zip(BASE, zip([(M.loc[keys[s], "AHHO"] < M.loc[keys[s], c]).sum() for c in BASE], holm(ps))))
w = {s: run_wtl(cec, keys[s], "le") for s in keys}
we = run_wtl(eng, list(eng.key.unique()), "best")
L = [r"\begin{table}[htbp]", r"\caption{AHHO against each optimizer}\label{tab:vs}", r"\centering", r"\small",
     r"\begin{tabular}{lcccccccc}", r"\toprule",
     r" & \multicolumn{3}{c}{CEC2017} & \multicolumn{3}{c}{CEC2020} & Engineering \\",
     r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-8}",
     r"Optimizer & AHHO better & $p$ Holm & W/T/L & AHHO better & $p$ Holm & W/T/L & W/T/L \\", r"\midrule"]
for c in BASE:
    a, b = acr["cec2017"][c], acr["cec2020"][c]
    L.append(f"{c} & {a[0]}/10 & {P(a[1])} & {'/'.join(map(str, w['cec2017'][c]))} & {b[0]}/10 & {P(b[1])} & "
             f"{'/'.join(map(str, w['cec2020'][c]))} & {'/'.join(map(str, we[c]))} \\\\")
L += [r"\bottomrule", r"\end{tabular}", r"\vspace{2pt}", "",
      r"{\small AHHO better counts the functions on which AHHO has the lower mean log error. W/T/L counts functions on which AHHO wins, ties or loses in run-level Wilcoxon tests, Holm-corrected within each function.}",
      r"\end{table}"]
open(f"{OUT}/vs.tex", "w").write("\n".join(L) + "\n")
log.append(("vs", acr, w, we))

# ---------- engineering
NAMES = {"PressureVessel": "Pressure vessel", "WeldedBeam": "Welded beam", "SpringDesign": "Spring", "SpeedReducer": "Speed reducer"}
EP = ["PressureVessel", "WeldedBeam", "SpringDesign", "SpeedReducer"]
L = [r"\begin{longtable}{llrrrrr}", r"\caption{Penalized fitness on the engineering problems over 30 runs}\label{tab:eng-desc}\\",
     r"\toprule", r"Problem & Optimizer & Min & Median & Mean & Max & STD \\", r"\midrule", r"\endfirsthead",
     r"\toprule", r"Problem & Optimizer & Min & Median & Mean & Max & STD \\", r"\midrule", r"\endhead", r"\bottomrule", r"\endfoot"]
L2 = [r"\begin{longtable}{llrrrrr}", r"\caption{Original objective, violation and feasibility on the engineering problems}\label{tab:eng-feas}\\",
      r"\toprule", r"Problem & Optimizer & Mean $f$ & Max violation & Feasible $10^{-4}$ & Feasible $10^{-6}$ & Best feasible $f$ \\",
      r"\midrule", r"\endfirsthead", r"\toprule",
      r"Problem & Optimizer & Mean $f$ & Max violation & Feasible $10^{-4}$ & Feasible $10^{-6}$ & Best feasible $f$ \\",
      r"\midrule", r"\endhead", r"\bottomrule", r"\endfoot"]
for j, pb in enumerate(EP):
    g = eng[eng.problem == pb]
    for i, c in enumerate(ORDER):
        h = g[g.config == c]; b = h.best
        nm = NAMES[pb] if i == 0 else ""
        L.append(f"{nm} & {c} & {E(b.min())} & {E(b.median())} & {E(b.mean())} & {E(b.max())} & {E(b.std())} \\\\")
        L2.append(f"{nm} & {c} & {E(h.f_raw.mean())} & {E(h.max_viol.max())} & {(h.max_viol <= 1e-4).mean():.2f} & "
                  f"{(h.max_viol <= 1e-6).mean():.2f} & {E(pd.to_numeric(h.best_feas, errors='coerce').min())} \\\\")
    if j < 3: L.append(r"\addlinespace"); L2.append(r"\addlinespace")
L.append(r"\end{longtable}"); L2.append(r"\end{longtable}")
open(f"{OUT}/eng_desc.tex", "w").write("\n".join(L) + "\n"); open(f"{OUT}/eng_feas.tex", "w").write("\n".join(L2) + "\n")
# best designs of AHHO
L = [r"\begin{table}[htbp]", r"\caption{Best AHHO designs on the engineering problems}\label{tab:eng-design}", r"\centering", r"\small",
     r"\begin{tabularx}{\textwidth}{l l X r r}", r"\toprule", r"Problem & Design & Variables & $f$ & Max violation \\", r"\midrule"]
for pb in EP:
    h = eng[(eng.problem == pb) & (eng.config == "AHHO")]
    r = h.loc[h.best.idxmin()]
    xs = ", ".join(f"{float(v):.5g}" for v in r.best_x.split(";"))
    L.append(f"{NAMES[pb]} & Lowest $F$ & {xs} & {r.f_raw:.6g} & {E(r.max_viol)} \\\\")
    hf = h.copy(); hf["bf"] = pd.to_numeric(hf.best_feas, errors="coerce"); rf = hf.loc[hf.bf.idxmin()]
    xs = ", ".join(f"{float(v):.5g}" for v in rf.best_feas_x.split(";"))
    from engineering_problems import PROBLEMS as _P
    _g = _P[pb]["raw"](np.array([float(v) for v in rf.best_feas_x.split(";")]))[1]
    L.append(f" & Best feasible & {xs} & {rf.bf:.6g} & {E(max(0.0, max(_g)))} \\\\")
L += [r"\bottomrule", r"\end{tabularx}", r"\end{table}"]
open(f"{OUT}/eng_design.tex", "w").write("\n".join(L) + "\n")

# ---------- ablation
A = pd.concat([main[main.config.isin(["HHO", "AHHO"])], abl])
VAR = ["HHO", "AHHO-EnergyOnly", "AHHO-PerturbOnly", "AHHO"]
Ac = A[A.suite != "eng"]; AM = Ac.groupby(["key", "config"])["le"].mean().unstack()[VAR]
AE = A[A.suite == "eng"].groupby(["key", "config"]).best.mean().unstack()[VAR]
ra = {s: AM.loc[keys[s]].rank(axis=1).mean() for s in keys}; rae = AE.rank(axis=1).mean()
pa = holm([pw(AM.loc[allk, "AHHO"], AM.loc[allk, c]) for c in VAR[:-1]]) .tolist() + [np.nan]
fa = friedmanchisquare(*[AM.loc[allk, c] for c in VAR])
L = [r"\begin{table}[htbp]", r"\caption{Factorial ablation: average ranks among the four variants}\label{tab:ablation-results}", r"\centering",
     r"\begin{tabular}{lcccccc}", r"\toprule",
     r"Variant & Energy & Perturbation & CEC2017 & CEC2020 & Engineering & $p$ Holm vs AHHO \\", r"\midrule"]
for c, en, pe, p in zip(VAR, ["No", "Yes", "No", "Yes"], ["No", "No", "Yes", "Yes"], pa):
    L.append(f"{'AHHO-Full' if c == 'AHHO' else c} & {en} & {pe} & {ra['cec2017'][c]:.2f} & {ra['cec2020'][c]:.2f} & {rae[c]:.2f} & {'--' if np.isnan(p) else P(p)} \\\\")
L += [r"\bottomrule", r"\end{tabular}", r"\vspace{2pt}", "",
      r"{\small The $p$ values compare mean log errors over the 20 CEC functions.}", r"\end{table}"]
open(f"{OUT}/ablation.tex", "w").write("\n".join(L) + "\n")
log.append(("ablation", fa.pvalue, ra, rae, pa))

# ---------- sensitivity
S = pd.concat([main[main.config == "AHHO"].assign(config="default"), sens])
SM = S[S.suite != "eng"].groupby(["key", "config"])["le"].mean().unstack()
SE = S[S.suite == "eng"].groupby(["key", "config"]).best.mean().unstack()
LV = [("alpha=0.25", r"$\alpha=0.25$"), ("alpha=1.0", r"$\alpha=1$"), ("p=0.15", r"$p=0.15$"), ("p=0.6", r"$p=0.6$"),
      ("bmax=0.25", r"$\beta_{max}=0.25$"), ("bmax=1.0", r"$\beta_{max}=1$"), ("bmin=0.01", r"$\beta_{min}=0.01$"),
      ("bmin=0.1", r"$\beta_{min}=0.1$"), ("delta=0", r"$\delta=0$"), ("delta=1e-4", r"$\delta=10^{-4}$"),
      ("delta=1e-2", r"$\delta=10^{-2}$"), ("frac=0.25", "fraction 0.25"), ("frac=0.75", "fraction 0.75"),
      ("frac=1.0", "fraction 1"), ("default", "Default")]
cols = [k for k, _ in LV]
rs = SM.loc[allk, cols].rank(axis=1).mean(); rse = SE[cols].rank(axis=1).mean()
fs = friedmanchisquare(*[SM.loc[allk, c] for c in cols])
psens = holm([pw(SM.loc[allk, k], SM.loc[allk, "default"]) for k in cols[:-1]])
L = [r"\begin{table}[htbp]", r"\caption{One-factor-at-a-time sensitivity of AHHO}\label{tab:sensitivity}", r"\centering",
     r"\begin{tabular}{lccccc}", r"\toprule",
     r"Setting & Rank, CEC Suites & Better than default & $p$ Holm & Rank, engineering \\", r"\midrule"]
for (k, t), p in zip(LV, list(psens) + [np.nan]):
    better = "--" if k == "default" else f"{int((SM.loc[allk, k] < SM.loc[allk, 'default']).sum())}/20"
    L.append(f"{t} & {rs[k]:.2f} & {better} & {'--' if np.isnan(p) else P(p)} & {rse[k]:.2f} \\\\")
L += [r"\bottomrule", r"\end{tabular}", r"\vspace{2pt}", "",
      r"{\small Ranks are among the 15 settings. Default: $\alpha=0.5$, $p=0.3$, $\beta_{max}=0.5$, $\beta_{min}=0.05$, $\delta=10^{-3}$, fraction 0.5.}",
      r"\end{table}"]
open(f"{OUT}/sensitivity.tex", "w").write("\n".join(L) + "\n")
log.append(("sens", fs.pvalue, rs.to_dict(), psens))

# ---------- collapse
CF = ["HHO", "AHHO delta=0", "AHHO delta=1e-2", "AHHO delta=1e-3"]
ps_all, rows = [], []
for k in sorted(col.key.unique(), key=lambda x: (x.split(":")[0], int(x.split(":")[1]))):
    g = col[col.key == k]; h = g[g.config == "HHO"].sort_values("run")["le"].values
    row = [k]
    for c in CF:
        row.append(g[g.config == c].error.median())
        if c != "HHO":
            ps_all.append(pw(g[g.config == c].sort_values("run")["le"].values, h))
    rows.append(row)
ph = holm(ps_all); mind = {}
for c in CF[1:]:
    t = col[col.config == c].trace.dropna()
    mind[c] = min(float(x.split("|")[1]) for s_ in t for x in s_.split(";"))
L = [r"\begin{table}[htbp]", r"\caption{Median final error after a forced collapse at 70\% of the budget}\label{tab:collapse}", r"\centering", r"\small",
     r"\begin{tabular}{lccccccc}", r"\toprule",
     r"Function & HHO & $\delta=0$ & $p$ Holm & $\delta=10^{-2}$ & $p$ Holm & $\delta=10^{-3}$ & $p$ Holm \\", r"\midrule"]
for i, r in enumerate(rows):
    s, p = r[0].split(":")
    L.append(f"{s.upper()} {lab(s, p)} & {E(r[1])} & {E(r[2])} & {P(ph[3 * i])} & {E(r[3])} & {P(ph[3 * i + 1])} & {E(r[4])} & {P(ph[3 * i + 2])} \\\\")
L += [r"\midrule", f"Minimum logged $D(t)$ & -- & {E(mind['AHHO delta=0'])} & & {E(mind['AHHO delta=1e-2'])} & & {E(mind['AHHO delta=1e-3'])} & \\\\",
      r"\bottomrule", r"\end{tabular}", r"\vspace{2pt}", "",
      r"{\small $\delta=10^{-3}$ is the default. The $p$ values compare each AHHO setting with HHO over paired runs.}", r"\end{table}"]
open(f"{OUT}/collapse.tex", "w").write("\n".join(L) + "\n")
log.append(("collapse", rows, ph, mind))

import pprint; pprint.pprint(log, width=180)
