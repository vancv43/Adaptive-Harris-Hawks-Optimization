# Phayung Meesad(1), Van Van Chau(1,2,*), and Tuan Minh Nguyen(2)
# (1) Department of Information Technology Management, King Mongkut's University of Technology North Bangkok,
#     Bangkok, Thailand, phayung.m@itd.kmutnb.ac.th
# (2) Faculty of Information Technology, Posts and Telecommunications Institute of Technology,
#     Ho Chi Minh City, Vietnam, vancv@ptit.edu.vn, minhtuan@ptit.edu.vn
# (*) Corresponding author
#
"""Optimizers of the AHHO study.

Unified function-evaluation rule: every call of the objective is one FE
(initial population, rapid-dive candidates, onlookers, scouts, perturbed
hawks). `Budget` stops a run at exactly maxFEs by raising `Exhausted`, and it
records the best-so-far value itself, so no algorithm can overspend or lose its
last improvement. Schedules use progress tau = FE/maxFE instead of t/T_max, so
an algorithm that spends more than N FE per iteration still completes its
schedule. An unchanged or already evaluated point is never evaluated again.
"""
from math import gamma, pi, sin
import numpy as np


class Exhausted(Exception):
    pass


class Budget:
    def __init__(self, func, max_evals, n_start=0, n_points=100, raw=None, tol=1e-4):
        self.func, self.max_evals, self.n = func, int(max_evals), int(n_start)
        self.best, self.best_x = np.inf, None
        self.raw, self.tol = raw, tol
        self.best_feas, self.best_feas_x = np.inf, None
        self.every = max(1, self.max_evals // n_points)
        self.curve = []

    @property
    def tau(self):
        return min(1.0, self.n / self.max_evals)

    def __call__(self, x):
        if self.n >= self.max_evals:
            raise Exhausted
        f = float(self.func(x))
        self.n += 1
        if f < self.best:
            self.best, self.best_x = f, np.array(x, dtype=float)
        if self.raw is not None:
            obj, g = self.raw(x)
            if max(g) <= self.tol and obj < self.best_feas:
                self.best_feas, self.best_feas_x = float(obj), np.array(x, dtype=float)
        if self.n % self.every == 0:
            self.curve.append(self.best)
        return f


def _init(n, d, lb, ub, rng):
    return lb + rng.random((n, d)) * (ub - lb)


def _ev(obj, X):
    return np.array([obj(x) for x in X])


def _levy(d, rng, b=1.5, scale=0.01):
    s = (gamma(1 + b) * sin(pi * b / 2) / (gamma((1 + b) / 2) * b * 2 ** ((b - 1) / 2))) ** (1 / b)
    # Factor 0.01 as in the Levy flight of Heidari et al. (2019).
    return scale * rng.normal(0, s, d) / (np.abs(rng.normal(0, 1, d)) ** (1 / b) + 1e-12)


def _hho_move(obj, X, fit, i, e, j, rabbit, Xm, lb, ub, d, n, rng, ls=0.01):
    """One HHO move. Returns (new position, its fitness or None if not evaluated)."""
    xi = X[i]
    if abs(e) >= 1:
        if rng.random() >= 0.5:
            xr = X[rng.integers(n)]
            r1, r2 = rng.random(), rng.random()
            return xr - r1 * np.abs(xr - 2 * r2 * xi), None
        r3, r4 = rng.random(), rng.random()
        return (rabbit - Xm) - r3 * (lb + r4 * (ub - lb)), None
    r = rng.random()
    if r >= 0.5:
        if abs(e) >= 0.5:
            return (rabbit - xi) - e * np.abs(j * rabbit - xi), None
        return rabbit - e * np.abs(rabbit - xi), None
    base = xi if abs(e) >= 0.5 else Xm
    Y = np.clip(rabbit - e * np.abs(j * rabbit - base), lb, ub)
    fy = obj(Y)
    if fy < fit[i]:
        return Y, fy
    Z = np.clip(Y + rng.random(d) * _levy(d, rng, scale=ls), lb, ub)
    fz = obj(Z)
    if fz < fit[i]:
        return Z, fz
    return xi.copy(), fit[i]


def run_family(obj, lb, ub, d, n, rng, alpha=0.5, beta_max=0.5, beta_min=0.05, p=0.3,
               delta=1e-3, frac=0.5, use_energy=True, use_perturb=True,
               init_X=None, trace=None, post=None, levy_scale=0.01):
    """HHO (both switches off) and AHHO share this one code path.

    D(t), sigma(t) and the fitness ranking are computed from the snapshot at the
    start of iteration t, before any position update of that iteration. Hawks
    are ranked with ties broken at random; the worst ceil(frac*N) are eligible,
    and at least the best one is always protected. Without ties and with
    frac=0.5 this equals "worse than the median". Unlike a strict comparison
    with the median, it still selects hawks when a collapsed swarm has equal
    fitness values, which is the case the floor delta exists for.
    """
    X = _init(n, d, lb, ub, rng) if init_X is None else init_X.copy()
    fit = _ev(obj, X)
    rng_span = ub - lb
    norm = np.linalg.norm(rng_span) + 1e-12
    rabbit = X[np.argmin(fit)].copy()
    rf = fit.min()
    while True:
        tau = obj.tau
        Dt = float(np.mean(np.linalg.norm(X - X.mean(axis=0), axis=1)) / norm)
        sig = X.std(axis=0)
        sig_eff = np.maximum(sig, delta * rng_span)
        beta = beta_max - (beta_max - beta_min) * tau
        perm = rng.permutation(n)
        order = perm[np.argsort(fit[perm], kind="stable")]
        eligible = np.zeros(n, dtype=bool)
        eligible[order[max(1, n - int(np.ceil(frac * n))):]] = True
        mod = (1 + alpha * Dt) if use_energy else 1.0
        E = 2 * (2 * rng.random(n) - 1) * (1 - tau) * mod
        J = 2 * (1 - rng.random(n))
        Xm = X.mean(axis=0)
        if trace is not None:
            trace.append((tau, Dt, float(np.min(sig / rng_span)), float(np.mean(np.abs(E) >= 1))))
        for i in range(n):
            new, fnew = _hho_move(obj, X, fit, i, E[i], J[i], rabbit, Xm, lb, ub, d, n, rng, levy_scale)
            if fnew is not None and fnew < rf:
                rf, rabbit = fnew, new.copy()
            if use_perturb and eligible[i] and rng.random() < p * (1 - Dt):
                new = new + beta * rng.normal(0, 1, d) * sig_eff
                fnew = None
            X[i] = np.clip(new, lb, ub)
            fit[i] = obj(X[i]) if fnew is None else fnew
            if fit[i] < rf:
                rf, rabbit = fit[i], X[i].copy()
        if post is not None:
            rf, rabbit = post(obj, X, fit, rf, rabbit, lb, ub, d, n, rng)


def run_hho(obj, lb, ub, d, n, rng, **kw):
    kw.setdefault("use_energy", False)
    kw.setdefault("use_perturb", False)
    return run_family(obj, lb, ub, d, n, rng, **kw)


def run_ahho(obj, lb, ub, d, n, rng, **kw):
    return run_family(obj, lb, ub, d, n, rng, **kw)


def _post_gchho(obj, X, fit, rf, rabbit, lb, ub, d, n, rng, theta=0.0):
    # Song et al. (2021) Eqs. 16-20. O = G * (X_r1 - X_r2), C = k * LF(D), k = Judge > theta,
    # X_new = X_raw + O * C, greedy. theta = 0 is their final setting; LF is their Eq. 10 (no factor 0.01).
    for i in range(n):
        r1, r2 = rng.integers(n, size=2)
        O_ = rng.normal(0, 1, d) * (X[r1] - X[r2])
        C = (rng.random(d) > theta) * _levy(d, rng, scale=1.0)
        c = np.clip(X[i] + O_ * C, lb, ub)
        fc = obj(c)
        if fc < fit[i]:
            X[i], fit[i] = c, fc
            if fc < rf:
                rf, rabbit = fc, c.copy()
    return rf, rabbit


def run_gchho(obj, lb, ub, d, n, rng, **kw):
    return run_family(obj, lb, ub, d, n, rng, use_energy=False, use_perturb=False,
                      post=_post_gchho, levy_scale=1.0, **kw)


def run_mhho(obj, lb, ub, d, n, rng, k=0.24, e_start=2.0, e_end=0.0):
    """m-HHO of Gupta et al. (2020): Eq. 14 energy, Eqs. 15-16 dives, greedy selection,
    Eqs. 17-18 opposition on the top Nop hawks. The source gives no numeric values for
    k, E_start and E_end; E_start=2, E_end=0 and k=0.24 are read from its Fig. 4. The
    step factor a and jump strength J follow the same decreasing formula, as it states.
    """
    X = _init(n, d, lb, ub, rng)
    fit = _ev(obj, X)
    rabbit = X[np.argmin(fit)].copy()
    rf = fit.min()
    while True:
        tau = obj.tau
        E = (e_start - e_end) * np.exp(-(tau ** 2) / k ** 2) + e_end
        a_, J = E, E
        Xm = X.mean(axis=0)
        for i in range(n):
            xi = X[i].copy()
            if abs(E) >= 1:
                if rng.random() >= 0.5:
                    xr = X[rng.integers(n)]
                    new = xr - rng.random() * np.abs(xr - 2 * rng.random() * xi)
                else:
                    new = (rabbit - Xm) - rng.random() * (lb + rng.random() * (ub - lb))
            else:
                r = rng.random()
                if r >= 0.5:
                    new = (rabbit - xi) - E * np.abs(J * rabbit - xi) if abs(E) >= 0.5 else rabbit - E * np.abs(rabbit - xi)
                else:
                    base = xi if abs(E) >= 0.5 else Xm
                    Y = np.clip(rabbit - E * np.abs(J * rabbit - base), lb, ub)
                    fy = obj(Y)
                    if fy < fit[i]:
                        X[i], fit[i] = Y, fy
                    else:
                        Z = np.clip(rabbit - E * np.abs(J * rabbit - base) + a_ * _levy(d, rng), lb, ub)
                        fz = obj(Z)
                        if fz < fit[i]:
                            X[i], fit[i] = Z, fz
                    if fit[i] < rf:
                        rf, rabbit = fit[i], X[i].copy()
                    continue
            c = np.clip(new, lb, ub)
            fc = obj(c)
            if fc < fit[i]:  # greedy selection
                X[i], fit[i] = c, fc
                if fc < rf:
                    rf, rabbit = fc, c.copy()
        nop = int(round(n - obj.tau * (n - 1)))  # Eq. 18
        for i in np.argsort(fit)[:nop]:  # opposites of the top-fitted hawks, Eq. 17
            c = lb + ub - X[i]
            fc = obj(c)
            if fc < fit[i]:
                X[i], fit[i] = c, fc
                if fc < rf:
                    rf, rabbit = fc, c.copy()


def run_pso(obj, lb, ub, d, n, rng, w_max=0.9, w_min=0.4, c1=2.0, c2=2.0):
    X = _init(n, d, lb, ub, rng)
    V = np.zeros((n, d))
    vmax = 0.2 * (ub - lb)
    fit = _ev(obj, X)
    pb, pf = X.copy(), fit.copy()
    while True:
        g = pb[np.argmin(pf)].copy()
        w = w_max - (w_max - w_min) * obj.tau
        V = np.clip(w * V + c1 * rng.random((n, d)) * (pb - X) + c2 * rng.random((n, d)) * (g - X), -vmax, vmax)
        X = np.clip(X + V, lb, ub)
        fit = _ev(obj, X)
        m = fit < pf
        pb[m], pf[m] = X[m], fit[m]


def run_gwo(obj, lb, ub, d, n, rng):
    # Alpha, beta and delta persist across iterations, as in Mirjalili et al. (2014).
    X = _init(n, d, lb, ub, rng)
    fit = _ev(obj, X)
    Lf, Lx = [np.inf] * 3, [np.zeros(d) for _ in range(3)]
    while True:
        for i in range(n):
            f = fit[i]
            for k in range(3):
                if f < Lf[k]:
                    Lf[k + 1:], Lx[k + 1:] = Lf[k:2], Lx[k:2]
                    Lf[k], Lx[k] = f, X[i].copy()
                    break
        a = 2 - 2 * obj.tau
        for i in range(n):
            s = 0
            for l in Lx:
                A, C = 2 * a * rng.random(d) - a, 2 * rng.random(d)
                s = s + (l - A * np.abs(C * l - X[i]))
            X[i] = np.clip(s / 3, lb, ub)
        fit = _ev(obj, X)


def run_ga(obj, lb, ub, d, n, rng, pc=0.8, pm=0.1, elite=0.1):
    X = _init(n, d, lb, ub, rng)
    fit = _ev(obj, X)
    ne = max(1, int(elite * n))
    while True:
        o = np.argsort(fit)
        X, fit = X[o], fit[o]
        new = []
        while len(new) < n - ne:
            a, b = rng.integers(0, n, 2)
            p1 = X[a] if fit[a] < fit[b] else X[b]
            a, b = rng.integers(0, n, 2)
            p2 = X[a] if fit[a] < fit[b] else X[b]
            if rng.random() < pc:
                w = rng.random(d)
                c = w * p1 + (1 - w) * p2
            else:
                c = p1.copy()
            m = rng.random(d) < pm
            c[m] += rng.normal(0, 0.1, d)[m] * (ub - lb)[m]
            new.append(np.clip(c, lb, ub))
        new = np.array(new)
        fn = _ev(obj, new)
        X = np.vstack([X[:ne], new])
        fit = np.concatenate([fit[:ne], fn])


def run_abc(obj, lb, ub, d, n, rng, limit=None):
    limit = limit or n * d // 2
    X = _init(n, d, lb, ub, rng)
    fit = _ev(obj, X)
    trial = np.zeros(n)

    def step(i):
        k = rng.integers(n)
        while k == i:
            k = rng.integers(n)
        j = rng.integers(d)
        c = X[i].copy()
        c[j] = X[i, j] + rng.uniform(-1, 1) * (X[i, j] - X[k, j])
        c = np.clip(c, lb, ub)
        f = obj(c)
        if f < fit[i]:
            X[i], fit[i], trial[i] = c, f, 0
        else:
            trial[i] += 1

    while True:
        for i in range(n):
            step(i)
        pr = (1 / (1 + fit - fit.min())) if np.ptp(fit) > 0 else np.ones(n)
        pr = pr / pr.sum()
        for _ in range(n):
            step(rng.choice(n, p=pr))
        for i in np.where(trial > limit)[0]:
            X[i] = lb + rng.random(d) * (ub - lb)
            fit[i] = obj(X[i])
            trial[i] = 0


def run_bat(obj, lb, ub, d, n, rng, fmin=0, fmax=2, A0=1.0, r0=0.5, alpha=0.9, gam=0.9):
    X = _init(n, d, lb, ub, rng)
    V = np.zeros((n, d))
    vmax = 0.2 * (ub - lb)
    fit = _ev(obj, X)
    A, r = np.full(n, A0), np.full(n, r0)
    best = X[np.argmin(fit)].copy()
    bf = fit.min()
    t = 0
    while True:
        for i in range(n):
            V[i] = np.clip(V[i] + (X[i] - best) * (fmin + (fmax - fmin) * rng.random()), -vmax, vmax)
            c = np.clip(X[i] + V[i], lb, ub)
            if rng.random() > r[i]:
                c = np.clip(best + rng.uniform(-1, 1, d) * A.mean(), lb, ub)  # Yang (2010)
            f = obj(c)
            if f <= fit[i] and rng.random() < A[i]:
                X[i], fit[i] = c, f
                A[i] *= alpha
                r[i] = r0 * (1 - np.exp(-gam * (t + 1)))
            if fit[i] < bf:
                bf, best = fit[i], X[i].copy()
        t += 1


def run_aco(obj, lb, ub, d, n, rng, q=0.5, xi=0.85):
    A = _init(n, d, lb, ub, rng)
    fit = _ev(obj, A)
    w = np.exp(-(np.arange(n) ** 2) / (2 * (q * n) ** 2))
    w = w / w.sum()
    while True:
        o = np.argsort(fit)
        A, fit = A[o], fit[o]
        S = np.zeros((n, d))
        for s in range(n):
            l = rng.choice(n, p=w)
            S[s] = np.clip(A[l] + xi * np.sum(np.abs(A - A[l]), axis=0) / (n - 1) * rng.normal(0, 1, d), lb, ub)
        fs = _ev(obj, S)
        MX, MF = np.vstack([A, S]), np.concatenate([fit, fs])
        k = np.argsort(MF)[:n]
        A, fit = MX[k], MF[k]


def run_aoa(obj, lb, ub, d, n, rng, alpha=5, mu=0.499):
    X = _init(n, d, lb, ub, rng)
    fit = _ev(obj, X)
    best = X[np.argmin(fit)].copy()
    bf = fit.min()
    scale = (ub - lb) * mu + lb
    while True:
        tau = obj.tau
        MOA, MOP = 0.2 + 0.8 * tau, 1 - tau ** (1 / alpha)
        r1, r2 = rng.random((n, d)), rng.random((n, d))
        div = best / (MOP + 1e-12) * scale
        mul = best * MOP * scale
        sub, add = best - MOP * scale, best + MOP * scale
        X = np.where(r1 > MOA, np.where(r2 < 0.5, div, mul), np.where(r2 < 0.5, sub, add))
        X = np.clip(X, lb, ub)
        fit = _ev(obj, X)
        if fit.min() < bf:
            bf, best = fit.min(), X[np.argmin(fit)].copy()


def run_scso(obj, lb, ub, d, n, rng, S=2.0):
    X = _init(n, d, lb, ub, rng)
    fit = _ev(obj, X)
    best = X[np.argmin(fit)].copy()
    bf = fit.min()
    while True:
        rg = S - S * obj.tau
        for i in range(n):
            r, R = rg * rng.random(), 2 * rg * rng.random() - rg
            if abs(R) <= 1:
                X[i] = best - r * np.abs(rng.random(d) * best - X[i]) * np.cos(2 * pi * rng.random())  # Seyyedabbasi and Kiani (2022)
            else:
                X[i] = r * (X[rng.integers(n)] - rng.random(d) * X[i])
        X = np.clip(X, lb, ub)
        fit = _ev(obj, X)
        if fit.min() < bf:
            bf, best = fit.min(), X[np.argmin(fit)].copy()


BASE = {"GA": run_ga, "PSO": run_pso, "GWO": run_gwo, "ACO": run_aco, "ABC": run_abc,
        "BAT": run_bat, "AOA": run_aoa, "SCSO": run_scso, "HHO": run_hho,
        "GCHHO": run_gchho, "m-HHO": run_mhho, "AHHO": run_ahho}


def run(algo, budget, lb, ub, d, n, rng, **kw):
    try:
        BASE[algo](budget, lb, ub, d, n, rng, **kw)
    except Exhausted:
        pass
    return budget
