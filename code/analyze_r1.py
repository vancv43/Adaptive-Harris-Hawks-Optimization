# Phayung Meesad(1), Van Van Chau(1,2,*), and Tuan Minh Nguyen(2,*)
# (1) Department of Information Technology Management, King Mongkut's University of Technology North Bangkok,
#     Bangkok, Thailand, phayung.m@itd.kmutnb.ac.th
# (2) Faculty of Information Technology, Posts and Telecommunications Institute of Technology,
#     Ho Chi Minh City, Vietnam, vancv@ptit.edu.vn, minhtuan@ptit.edu.vn
# (*) Corresponding author
#
"""Analysis of the runs: mean log10(error+1e-8), Friedman ranks, Holm-corrected Wilcoxon tests, feasibility.

Usage: python3 analyze_r1.py RESULTS_DIR OUT_DIR
"""
import glob, os, sys
import numpy as np, pandas as pd
from scipy.stats import friedmanchisquare, rankdata, wilcoxon

RES, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
EPS = 1e-8
ORDER = ["GA", "PSO", "ACO", "ABC", "BAT", "GWO", "AOA", "SCSO", "HHO", "GCHHO", "m-HHO", "AHHO"]


def load(pat):
    return pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(RES, pat)))], ignore_index=True)


def holm(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); adj = np.empty(m); run = 0
    for i, k in enumerate(o):
        run = max(run, (m - i) * p[k]); adj[k] = min(1.0, run)
    return adj


def pw(a, b):
    d = np.asarray(a) - np.asarray(b)
    if np.all(d == 0):
        return 1.0
    return wilcoxon(a, b, zero_method="zsplit").pvalue


def fid(s, p):
    return ("F" if s == "cec2017" else "G") + str(p)


main = load("main-s*/r1_main_*.csv")
abl = load("abl-col/r1_ablation_*.csv")
sens = load("sens-s*/r1_sensitivity_*.csv")
col = load("abl-col/r1_collapse_*.csv")
for d in (main, abl, sens, col):
    d["fn"] = [fid(s, p) if s != "eng" else p for s, p in zip(d.suite, d.problem)]
cec = main[main.suite != "eng"].copy()
cec["le"] = np.log10(cec.error.astype(float).clip(lower=0) + EPS)
fns = list(dict.fromkeys(cec.sort_values(["suite", "problem"]).fn))

# ---- descriptive table: Min Median Mean Max STD of error, plus mean log10 error
desc = cec.groupby(["fn", "config"]).agg(Min=("error", "min"), Median=("error", "median"), Mean=("error", "mean"),
                                         Max=("error", "max"), STD=("error", "std"), MeanLog=("le", "mean")).reset_index()
desc.to_csv(f"{OUT}/cec_descriptive.csv", index=False)

# ---- primary: Friedman over functions on mean log10 error
M = desc.pivot(index="fn", columns="config", values="MeanLog").loc[fns, ORDER]
R = M.rank(axis=1, method="average")
fr = friedmanchisquare(*[M[c] for c in ORDER])
rank = R.mean().sort_values()
pd.DataFrame({"mean_rank": rank}).to_csv(f"{OUT}/friedman_primary.csv")
# Holm Wilcoxon across functions, AHHO vs each
rows = []
for c in ORDER[:-1]:
    rows.append([c, (M["AHHO"] < M[c]).sum(), (M["AHHO"] > M[c]).sum(), pw(M["AHHO"], M[c])])
fw = pd.DataFrame(rows, columns=["vs", "AHHO_better_fns", "AHHO_worse_fns", "p"]); fw["p_holm"] = holm(fw.p)
fw.to_csv(f"{OUT}/wilcoxon_across_functions.csv", index=False)
# secondary raw-mean Friedman
MR = desc.pivot(index="fn", columns="config", values="Mean").loc[fns, ORDER]
rank_raw = MR.rank(axis=1).mean().sort_values()
fr_raw = friedmanchisquare(*[MR[c] for c in ORDER])

# ---- per-function run-level paired Wilcoxon (paired by run seed), Holm over the 11 comparisons
wtl = []
for f in fns:
    g = cec[cec.fn == f]; a = g[g.config == "AHHO"].sort_values("run")["le"].values
    ps, sg = [], []
    for c in ORDER[:-1]:
        b = g[g.config == c].sort_values("run")["le"].values
        ps.append(pw(a, b)); sg.append(np.sign(np.median(b - a)))
    for c, p, s in zip(ORDER[:-1], holm(ps), sg):
        res = "tie" if p >= 0.05 or s == 0 else ("win" if s > 0 else "loss")
        wtl.append([f, c, p, res])
wtl = pd.DataFrame(wtl, columns=["fn", "vs", "p_holm", "result"])
wtl.to_csv(f"{OUT}/per_function_wtl.csv", index=False)
wtl_sum = wtl.groupby("vs").result.value_counts().unstack(fill_value=0).reindex(ORDER[:-1])[["win", "tie", "loss"]]
wtl_sum.to_csv(f"{OUT}/wtl_summary.csv")

# ---- engineering
eng = main[main.suite == "eng"].copy()
eng["feas4"] = eng.max_viol <= 1e-4; eng["feas6"] = eng.max_viol <= 1e-6
et = eng.groupby(["fn", "config"]).agg(Min=("best", "min"), Median=("best", "median"), Mean=("best", "mean"),
                                       Max=("best", "max"), STD=("best", "std"), RawMean=("f_raw", "mean"),
                                       MaxViol=("max_viol", "max"), Feas4=("feas4", "mean"), Feas6=("feas6", "mean"),
                                       BestFeas=("best_feas", "min")).reset_index()
bx = eng.loc[eng.groupby(["fn", "config"]).best.idxmin(), ["fn", "config", "best_x", "f_raw", "max_viol"]]
et = et.merge(bx.rename(columns={"f_raw": "bx_raw", "max_viol": "bx_viol"}), on=["fn", "config"])
et.to_csv(f"{OUT}/engineering.csv", index=False)
eng_wtl = []
for f in eng.fn.unique():
    g = eng[eng.fn == f]; a = g[g.config == "AHHO"].sort_values("run").best.values; ps, sg = [], []
    for c in ORDER[:-1]:
        b = g[g.config == c].sort_values("run").best.values; ps.append(pw(a, b)); sg.append(np.sign(np.median(b - a)))
    for c, p, s in zip(ORDER[:-1], holm(ps), sg):
        eng_wtl.append([f, c, p, "tie" if p >= 0.05 or s == 0 else ("win" if s > 0 else "loss")])
eng_wtl = pd.DataFrame(eng_wtl, columns=["fn", "vs", "p_holm", "result"]); eng_wtl.to_csv(f"{OUT}/engineering_wtl.csv", index=False)
eng_rank = et.pivot(index="fn", columns="config", values="Mean")[ORDER].rank(axis=1).mean().sort_values()

# ---- ablation (HHO and AHHO from main)
A = pd.concat([main[main.config.isin(["HHO", "AHHO"])], abl])
Ac = A[A.suite != "eng"].copy(); Ac["le"] = np.log10(Ac.error.astype(float).clip(lower=0) + EPS)
AM = Ac.groupby(["fn", "config"])["le"].mean().unstack().loc[fns, ["HHO", "AHHO-EnergyOnly", "AHHO-PerturbOnly", "AHHO"]]
abl_rank = AM.rank(axis=1).mean(); abl_fr = friedmanchisquare(*[AM[c] for c in AM.columns])
abl_pair = {c: (int((AM["AHHO"] < AM[c]).sum()), int((AM["AHHO"] > AM[c]).sum()), pw(AM["AHHO"], AM[c])) for c in AM.columns[:-1]}
AE = A[A.suite == "eng"].groupby(["fn", "config"]).best.mean().unstack()[AM.columns]
AM.to_csv(f"{OUT}/ablation_meanlog.csv"); AE.to_csv(f"{OUT}/ablation_eng_mean.csv")

# ---- sensitivity: default = AHHO from main
S = pd.concat([main[main.config == "AHHO"].assign(config="default"), sens])
Sc = S[S.suite != "eng"].copy(); Sc["le"] = np.log10(Sc.error.astype(float).clip(lower=0) + EPS)
SM = Sc.groupby(["fn", "config"])["le"].mean().unstack().loc[fns]
sens_rank = SM.rank(axis=1).mean().sort_values(); sens_fr = friedmanchisquare(*[SM[c] for c in SM.columns])
sp = []
for c in SM.columns:
    if c == "default": continue
    sp.append([c, int((SM[c] < SM["default"]).sum()), int((SM[c] > SM["default"]).sum()), pw(SM[c], SM["default"])])
sp = pd.DataFrame(sp, columns=["level", "better_than_default", "worse", "p"]); sp["p_holm"] = holm(sp.p)
sp.to_csv(f"{OUT}/sensitivity_vs_default.csv", index=False)
SE = S[S.suite == "eng"].groupby(["fn", "config"]).best.mean().unstack()

# ---- collapse
col["le"] = np.log10(col.error.astype(float).clip(lower=0) + EPS)
CM = col.groupby(["fn", "config"]).agg(mean_err=("error", "mean"), med_err=("error", "median"), meanlog=("le", "mean")).reset_index()
col_t = []
for f in col.fn.unique():
    g = col[col.fn == f]; h = g[g.config == "HHO"].sort_values("run")["le"].values
    for c in ["AHHO delta=0", "AHHO delta=1e-3", "AHHO delta=1e-2"]:
        a = g[g.config == c].sort_values("run")["le"].values
        col_t.append([f, c, np.median(h - a), pw(a, h)])
col_t = pd.DataFrame(col_t, columns=["fn", "cfg", "median_log_gain_vs_HHO", "p"]); col_t["p_holm"] = holm(col_t.p)
CM.to_csv(f"{OUT}/collapse.csv", index=False); col_t.to_csv(f"{OUT}/collapse_tests.csv", index=False)
# sigma floor evidence from AHHO traces in collapse runs: min sigma/range field
sig = {}
for c in ["AHHO delta=0", "AHHO delta=1e-3", "AHHO delta=1e-2"]:
    t = col[col.config == c].trace.dropna()
    vals = [float(x.split("|")[2]) for s in t for x in s.split(";")]
    sig[c] = (min(vals), np.median(vals))

pd.set_option("display.width", 200); pd.set_option("display.max_columns", 30)
print("PRIMARY Friedman chi2=%.2f p=%.2e" % (fr.statistic, fr.pvalue)); print(rank.round(2).to_string())
print("\nRaw-mean Friedman p=%.2e" % fr_raw.pvalue); print(rank_raw.round(2).to_string())
print("\nAcross-function Wilcoxon AHHO vs"); print(fw.round(4).to_string(index=False))
print("\nRun-level W/T/L (AHHO perspective, 20 CEC functions)"); print(wtl_sum.to_string())
print("\nPer-function where AHHO loses:"); print(wtl[wtl.result == "loss"].groupby("fn").vs.apply(list).to_string())
print("\nAHHO rank per function:"); print(R["AHHO"].to_string())
print("\nEngineering mean-rank:"); print(eng_rank.round(2).to_string())
print(et[["fn", "config", "Mean", "STD", "MaxViol", "Feas4", "Feas6", "BestFeas"]].round(5).to_string(index=False))
print("\nEngineering W/T/L"); print(eng_wtl.groupby("vs").result.value_counts().unstack(fill_value=0).reindex(ORDER[:-1]).to_string())
print("\nABLATION Friedman p=%.3g" % abl_fr.pvalue); print(abl_rank.round(2).to_string()); print(abl_pair); print(AE.round(4).to_string())
print("\nSENSITIVITY Friedman p=%.3g" % sens_fr.pvalue); print(sens_rank.round(2).to_string()); print(sp.round(4).to_string(index=False))
print(SE.round(4).to_string())
print("\nCOLLAPSE"); print(CM.round(5).to_string(index=False)); print(col_t.round(4).to_string(index=False)); print("sigma/range min,median:", sig)
