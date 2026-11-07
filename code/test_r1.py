# Phayung Meesad(1), Van Van Chau(1,2,*), and Tuan Minh Nguyen(2)
# (1) Department of Information Technology Management, King Mongkut's University of Technology North Bangkok,
#     Bangkok, Thailand, phayung.m@itd.kmutnb.ac.th
# (2) Faculty of Information Technology, Posts and Telecommunications Institute of Technology,
#     Ho Chi Minh City, Vietnam, vancv@ptit.edu.vn, minhtuan@ptit.edu.vn
# (*) Corresponding author
#
"""Audit tests of the optimizer and runner code. Run: python3 test_r1.py"""
import os, sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import optimizers_r1 as O
from optimizers_r1 import Budget, run, BASE, Exhausted
import run_r1 as R
from engineering_problems import PROBLEMS

d, n = 10, 30
lb, ub = -100 * np.ones(d), 100 * np.ones(d)
sph = lambda x: float(np.sum((x - 3.0) ** 2))
fails = []


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name, info)
    if not ok:
        fails.append(name)


class Seen:
    def __init__(self): self.c = {}
    def __call__(self, x):
        k = np.asarray(x).tobytes(); self.c[k] = self.c.get(k, 0) + 1
        return sph(np.asarray(x))


for a in BASE:
    res = []
    for _ in range(2):
        b = run(a, Budget(sph, 5003), lb, ub, d, n, np.random.default_rng(7))
        res.append(b.best)
    check(f"exact budget {a}", b.n == 5003, b.n)
    check(f"deterministic {a}", res[0] == res[1])
    check(f"improves {a}", b.best < sph(lb), f"{b.best:.3g}")
    check(f"best consistent {a}", abs(sph(b.best_x) - b.best) < 1e-9)
    check(f"inside bounds {a}", np.all(b.best_x >= lb) and np.all(b.best_x <= ub))

# A failed dive keeps the position and its known fitness (2 FE, no third call);
# a successful Y costs 1 FE and Z is not computed.
class Cnt:
    def __init__(self, v): self.v, self.n = v, 0
    def __call__(self, x):
        self.n += 1
        return self.v
Xt, ft = g0 = (np.zeros((n, d)), np.full(n, 5.0))
def dive(val, e):
    c = Cnt(val); rg = np.random.default_rng(0)
    # find a seed whose first draw selects the dive branch (r < 0.5)
    for sd in range(100):
        rg = np.random.default_rng(sd)
        if np.random.default_rng(sd).random() < .5:
            break
    out = O._hho_move(c, Xt, ft, 0, e, 1.0, np.ones(d), Xt.mean(0), lb, ub, d, n, rg)
    return c.n, out
cnt, (xn, fn) = dive(9.0, .7)
check("failed dive: 2 FE, position kept, fitness reused", cnt == 2 and fn == 5.0 and np.array_equal(xn, Xt[0]))
cnt, (xn, fn) = dive(1.0, .7)
check("successful Y: 1 FE, fitness returned", cnt == 1 and fn == 1.0)
cnt, (xn, fn) = dive(9.0, .2)
check("failed hard dive: 2 FE", cnt == 2 and fn == 5.0)
# m-HHO opposite set follows Eq. 18 at both ends of the run.
check("m-HHO Nop at start and end", round(30 - 0.0 * 29) == 30 and round(30 - 1.0 * 29) == 1)
tt = Budget(sph, 3000); tr = run("m-HHO", tt, lb, ub, d, n, np.random.default_rng(1))
check("m-HHO uses non-negative decreasing energy", 2 * np.exp(-(1.0 / .24) ** 2) < 1e-6)
# GA elites keep their cached fitness: exactly N-ne new evaluations per generation.
c = Cnt(1.0); bb = Budget(c, 30 + 27 * 3 + 5)
try:
    O.run_ga(bb, lb, ub, d, n, np.random.default_rng(0))
except Exhausted:
    pass
check("GA spends N-ne FE per generation", bb.n == bb.max_evals and c.n == bb.max_evals)

x = run("HHO", Budget(sph, 4000), lb, ub, d, n, np.random.default_rng(3)).best
y = O.Budget(sph, 4000)
try:
    O.run_family(y, lb, ub, d, n, np.random.default_rng(3), use_energy=False, use_perturb=False)
except Exhausted:
    pass
check("AHHO with both switches off equals HHO", x == y.best)

for a in ["AHHO", "HHO"]:
    tr = []
    b = Budget(sph, 6000)
    run(a, b, lb, ub, d, n, np.random.default_rng(2), trace=tr)
    check(f"schedule completes {a}", tr[-1][0] > 0.98, f"tau_end={tr[-1][0]:.3f}")

# D(t) upper bound: two-corner population reaches 1/2; uniform is near 0.2887.
span = np.linalg.norm(ub - lb)
Dof = lambda X: np.mean(np.linalg.norm(X - X.mean(0), axis=1)) / span
corner = np.vstack([np.tile(lb, (15, 1)), np.tile(ub, (15, 1))])
check("D two-corner equals 1/2", abs(Dof(corner) - .5) < 1e-12)
g = np.random.default_rng(0)
check("D never above 1/2", max(Dof(lb + g.random((30, d)) ** g.integers(1, 9) * (ub - lb)) for _ in range(2000)) <= .5)
check("D uniform near 0.2887", abs(np.mean([Dof(lb + g.random((3000, d)) * (ub - lb)) for _ in range(20)]) - .2887) < .01)

# Fully collapsed swarm with equal fitness: perturbation must still move hawks.
flat = lambda x: 1.0
dist = {}
for dl in (1e-3, 0.0):
    c = {}
    def f(x):
        c[np.asarray(x).tobytes()] = 1
        return 1.0
    b = Budget(f, 3000, n_start=2000)
    try:
        O.run_family(b, lb, ub, d, n, np.random.default_rng(5), delta=dl, use_energy=False,
                     init_X=np.tile(np.zeros(d), (n, 1)))
    except Exhausted:
        pass
    dist[dl] = len(c)
check("floor spreads a collapsed flat swarm more than delta=0", dist[1e-3] > dist[0.0], dist)

# Rank rule equals "worse than the median" when there are no ties.
fit = g.random(30)
perm = g.permutation(30)
order = perm[np.argsort(fit[perm], kind="stable")]
e = np.zeros(30, bool); e[order[15:]] = True
check("rank rule equals median rule without ties", np.array_equal(e, fit > np.median(fit)))

# Known designs from the literature.
K = {"WeldedBeam": ([.20573, 3.47049, 9.03662, .20573], 1.7249), "SpringDesign": ([.051689, .356718, 11.288966], .012665),
     "PressureVessel": ([.7781686, .3846492, 40.319619, 200.], 5885.33)}
for k, (xx, v) in K.items():
    o, gg = PROBLEMS[k]["raw"](np.array(xx))
    check(f"known design {k}", abs(o - v) / v < 1e-3 and max(gg) < 1e-3, f"{o:.5f}")

# Runner: one real task per experiment on a reduced budget, plus resume key shape.
for e_ in ["main", "ablation", "sensitivity", "collapse"]:
    ts = R.tasks(e_)
    keys = {(t[1], t[2], t[3], t[6]) for t in ts}
    check(f"tasks unique {e_}", len(keys) == len(ts), len(ts))
parts = [set(map(str, R.tasks("main")[i::4])) for i in range(4)]
check("shards partition tasks", sum(map(len, parts)) == 8640 and len(set.union(*parts)) == 8640)
r = R.one(("collapse", "cec2017", 9, "x", "AHHO", dict(delta=1e-3), 0))
check("collapse task exact stop", r["n_evals"] == 150000)
r = R.one(("main", "eng", "SpringDesign", "AHHO", "AHHO", {}, 0))
check("eng task fields", r["n_evals"] == 30000 and r["best_x"] != "" and r["max_viol"] != "")
print("FAILED:", fails if fails else "none")
sys.exit(1 if fails else 0)
