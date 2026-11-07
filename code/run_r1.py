# Phayung Meesad(1), Van Van Chau(1,2,*), and Tuan Minh Nguyen(2)
# (1) Department of Information Technology Management, King Mongkut's University of Technology North Bangkok,
#     Bangkok, Thailand, phayung.m@itd.kmutnb.ac.th
# (2) Faculty of Information Technology, Posts and Telecommunications Institute of Technology,
#     Ho Chi Minh City, Vietnam, vancv@ptit.edu.vn, minhtuan@ptit.edu.vn
# (*) Corresponding author
#
"""Task runner. Usage: run_r1.py EXP --shard k --nshards m --workers w --out DIR"""
import argparse, csv, os, sys, time, traceback, warnings
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"  # one BLAS thread per worker
from multiprocessing import Pool
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from optimizers_r1 import Budget, run, BASE
from engineering_problems import PROBLEMS

CEC17 = [1, 3, 5, 7, 9, 13, 14, 21, 23, 24]
CEC20 = list(range(1, 11))
ALL = [("cec2017", f) for f in CEC17] + [("cec2020", f) for f in CEC20] + [("eng", k) for k in PROBLEMS]
RUNS, POP, DIM = 30, 30, 30
AB = {"AHHO-EnergyOnly": dict(use_perturb=False), "AHHO-PerturbOnly": dict(use_energy=False)}
SENS = {"alpha=0.25": dict(alpha=.25), "alpha=1.0": dict(alpha=1.), "p=0.15": dict(p=.15), "p=0.6": dict(p=.6),
        "bmax=0.25": dict(beta_max=.25), "bmax=1.0": dict(beta_max=1.), "bmin=0.01": dict(beta_min=.01),
        "bmin=0.1": dict(beta_min=.1), "delta=0": dict(delta=0.), "delta=1e-4": dict(delta=1e-4),
        "delta=1e-2": dict(delta=1e-2), "frac=0.25": dict(frac=.25), "frac=0.75": dict(frac=.75), "frac=1.0": dict(frac=1.)}
# Confirmation of the sensitivity finding on fresh seeds (test fixed before the runs).
CONFIRM = {"default": dict(), "frac=1.0": dict(frac=1.)}
CONFIRM_OFFSET = 500000  # disjoint from every seed used in main, ablation, sensitivity and collapse
COLLAPSE_FUNCS = [("cec2017", 5), ("cec2017", 7), ("cec2017", 9), ("cec2017", 21), ("cec2020", 2), ("cec2020", 3)]
COLLAPSE = {"HHO": dict(use_energy=False, use_perturb=False), "AHHO delta=0": dict(delta=0.),
            "AHHO delta=1e-3": dict(delta=1e-3), "AHHO delta=1e-2": dict(delta=1e-2)}
FIELDS = ["exp", "suite", "problem", "config", "run", "seed", "best", "error", "n_evals", "rt", "curve", "trace",
          "f_raw", "max_viol", "best_x", "best_feas", "best_feas_x"]


def problem(suite, pid):
    if suite == "eng":
        p = PROBLEMS[pid]
        return p["fitness"], p["lb"], p["ub"], p["dim"], 30000, None, p["raw"]
    from opfunu.cec_based import cec2017, cec2020
    mod, tag = (cec2017, "2017") if suite == "cec2017" else (cec2020, "2020")
    f = getattr(mod, f"F{pid}{tag}")(ndim=DIM)
    return f.evaluate, np.array(f.lb, float), np.array(f.ub, float), DIM, 150000, f.f_global, None


def seed(suite, pid, r, offset=0):
    return offset + {"cec2017": 1000 * pid if suite == "cec2017" else 0, "cec2020": 2000 * pid if suite == "cec2020" else 0,
            "eng": 2000}[suite] + r


def tasks(exp):
    if exp == "main":
        return [(exp, s, p, a, a, {}, r) for s, p in ALL for a in BASE for r in range(RUNS)]
    if exp == "ablation":
        return [(exp, s, p, c, "AHHO", k, r) for s, p in ALL for c, k in AB.items() for r in range(RUNS)]
    if exp == "sensitivity":
        return [(exp, s, p, c, "AHHO", k, r) for s, p in ALL for c, k in SENS.items() for r in range(RUNS)]
    if exp == "confirm":
        return [(exp, s, p, c, "AHHO", k, r) for s, p in ALL for c, k in CONFIRM.items() for r in range(RUNS)]
    if exp == "collapse":
        return [(exp, s, p, c, "AHHO", k, r) for s, p in COLLAPSE_FUNCS for c, k in COLLAPSE.items() for r in range(RUNS)]
    raise SystemExit(exp)


def safe(t):
    try:
        return one(t)
    except Exception:
        return {"exp": "ERROR", "suite": t[1], "problem": t[2], "config": t[3], "run": t[6], "trace": traceback.format_exc()}


def one(t):
    exp, s, p, cfg, algo, kw, r = t
    func, lb, ub, d, maxfe, fopt, raw = problem(s, p)
    sd = seed(s, p, r, CONFIRM_OFFSET if exp == "confirm" else 0)
    rng = np.random.default_rng(sd)
    kw = dict(kw)
    t0 = time.time()
    tr = []
    if algo == "AHHO":
        kw["trace"] = tr
    b = Budget(func, maxfe, raw=raw)
    if exp == "collapse":
        # Controlled test: start at 70% of the budget with every hawk within
        # 1e-8 of the range around one random point, so |E|<1 throughout.
        b.n = int(0.7 * maxfe)
        x0 = lb + rng.random(d) * (ub - lb)
        kw["init_X"] = np.clip(x0 + 1e-8 * (ub - lb) * rng.normal(0, 1, (POP, d)), lb, ub)
    run(algo, b, lb, ub, d, POP, rng, **kw)
    row = dict(exp=exp, suite=s, problem=p, config=cfg, run=r, seed=sd, best=b.best,
               error=(b.best - fopt) if fopt is not None else "", n_evals=b.n, rt=round(time.time() - t0, 2),
               curve=";".join(f"{v:.8g}" for v in b.curve))
    if tr:
        idx = np.unique(np.linspace(0, len(tr) - 1, 21).astype(int))
        row["trace"] = ";".join("|".join(f"{v:.5g}" for v in tr[i]) for i in idx)
    if raw is not None:
        fr, g = raw(b.best_x)
        row.update(f_raw=fr, max_viol=max(0.0, max(g)), best_x=";".join(f"{v:.10g}" for v in b.best_x),
                   best_feas=b.best_feas if np.isfinite(b.best_feas) else "",
                   best_feas_x=";".join(f"{v:.10g}" for v in b.best_feas_x) if b.best_feas_x is not None else "")
    return row


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("exp"); ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--out", default="."); ap.add_argument("--limit", type=int, default=0); ap.add_argument("--max-hours", type=float, default=11.0)
    a = ap.parse_args()
    all_t = tasks(a.exp)
    ts = all_t[a.shard::a.nshards]
    if a.limit:
        ts = ts[:a.limit]
    path = os.path.join(a.out, f"r1_{a.exp}_s{a.shard}.csv")
    done = set()
    if os.path.exists(path):
        with open(path) as fh:
            done = {(x["suite"], x["problem"], x["config"], x["run"]) for x in csv.DictReader(fh) if x["exp"] == a.exp and x["rt"]}
    ts = [t for t in ts if (t[1], str(t[2]), t[3], str(t[6])) not in done]
    print(a.exp, "shard", a.shard, "todo", len(ts), "of", len(all_t), flush=True)
    new = not os.path.exists(path)
    t0 = time.time()
    with open(path, "a", newline="") as fh, Pool(a.workers) as pool:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if new:
            w.writeheader()
        errors = 0
        for i, row in enumerate(pool.imap_unordered(safe, ts, chunksize=1)):
            w.writerow(row); fh.flush()
            if row["exp"] == "ERROR":
                errors += 1
                print(row["trace"], flush=True)
            if (time.time() - t0) / 3600 > a.max_hours:
                # Stop before an external time limit so that the rows written so far are kept.
                pool.terminate()
                print("STOPPED_AT_TIME_LIMIT", i + 1, "of", len(ts), flush=True)
                sys.exit(3)
            if i % 200 == 0:
                print(i, f"{(time.time()-t0)/60:.1f}min", flush=True)
    print("DONE" if not errors else f"FINISHED_WITH_{errors}_ERRORS", a.exp, a.shard, f"{(time.time()-t0)/60:.1f}min")
    sys.exit(1 if errors else 0)
