# Phayung Meesad(1), Van Van Chau(1,2,*), and Tuan Minh Nguyen(2,*)
# (1) Department of Information Technology Management, King Mongkut's University of Technology North Bangkok,
#     Bangkok, Thailand, phayung.m@itd.kmutnb.ac.th
# (2) Faculty of Information Technology, Posts and Telecommunications Institute of Technology,
#     Ho Chi Minh City, Vietnam, vancv@ptit.edu.vn, minhtuan@ptit.edu.vn
# (*) Corresponding author
#
"""
Classic constrained engineering design benchmark problems, used as Part B
real-world validation alongside the CEC2017 numerical suite. Constraints are
handled with a static quadratic penalty on inequality violations.

Each problem exposes: dim, lb, ub, fitness(x), and a literature-reported
best objective value (best_known) used only for reference/sanity checks in
the reporting stage (not used by the optimizers).
"""
import numpy as np

PENALTY = 1.0e7


def _penalty(gs):
    return PENALTY * sum(max(0.0, g) ** 2 for g in gs)


# ---------------------------------------------------------- Pressure vessel
# NOTE: Kannan & Kramer (1994) define x1, x2
# (shell/head thickness) as discrete multiples of 0.0625in rolled-plate
# stock; this implementation uses the fully continuous relaxation of all
# four variables (bounds below), not the mixed discrete-continuous
# formulation of the original source. This is a deliberate simplification,
# stated here explicitly, not an oversight.
def pressure_vessel_raw(x):
    x1, x2, x3, x4 = x
    obj = (0.6224 * x1 * x3 * x4 + 1.7781 * x2 * x3 ** 2
           + 3.1661 * x1 ** 2 * x4 + 19.84 * x1 ** 2 * x3)
    g = [
        -x1 + 0.0193 * x3,
        -x2 + 0.00954 * x3,
        -np.pi * x3 ** 2 * x4 - (4 / 3) * np.pi * x3 ** 3 + 1296000,
        x4 - 240,
    ]
    return obj, g


def pressure_vessel(x):
    obj, g = pressure_vessel_raw(x)
    return obj + _penalty(g)


PRESSURE_VESSEL = dict(
    name="Pressure Vessel Design", dim=4,
    lb=np.array([0.0625, 0.0625, 10.0, 10.0]),
    ub=np.array([6.1875, 6.1875, 200.0, 200.0]),
    fitness=pressure_vessel, raw=pressure_vessel_raw, best_known=6059.7143,
)


# ------------------------------------------------------------- Welded beam
def welded_beam_raw(x):
    x1, x2, x3, x4 = x
    P, L, E, G = 6000.0, 14.0, 30e6, 12e6
    tau_max, sigma_max, delta_max = 13600.0, 30000.0, 0.25

    M = P * (L + x2 / 2)
    R = np.sqrt(x2 ** 2 / 4 + ((x1 + x3) / 2) ** 2)
    J = 2 * (np.sqrt(2) * x1 * x2 * (x2 ** 2 / 12 + ((x1 + x3) / 2) ** 2))
    tau1 = P / (np.sqrt(2) * x1 * x2)
    tau2 = M * R / J
    tau = np.sqrt(tau1 ** 2 + 2 * tau1 * tau2 * (x2 / (2 * R)) + tau2 ** 2)
    sigma = 6 * P * L / (x4 * x3 ** 2)
    delta = 4 * P * L ** 3 / (E * x3 ** 3 * x4)
    Pc = (4.013 * E * np.sqrt(x3 ** 2 * x4 ** 6 / 36) / L ** 2
          * (1 - x3 / (2 * L) * np.sqrt(E / (4 * G))))

    obj = 1.10471 * x1 ** 2 * x2 + 0.04811 * x3 * x4 * (14.0 + x2)
    g = [
        tau - tau_max,
        sigma - sigma_max,
        x1 - x4,
        0.10471 * x1 ** 2 + 0.04811 * x3 * x4 * (14 + x2) - 5.0,
        0.125 - x1,
        delta - delta_max,
        P - Pc,
    ]
    return obj, g


def welded_beam(x):
    obj, g = welded_beam_raw(x)
    return obj + _penalty(g)


WELDED_BEAM = dict(
    name="Welded Beam Design", dim=4,
    lb=np.array([0.1, 0.1, 0.1, 0.1]),
    ub=np.array([2.0, 10.0, 10.0, 2.0]),
    fitness=welded_beam, raw=welded_beam_raw, best_known=1.724852,
)


# --------------------------------------------------- Tension/compression spring
def spring_design_raw(x):
    x1, x2, x3 = x
    obj = (x3 + 2) * x2 * x1 ** 2
    g = [
        1 - (x2 ** 3 * x3) / (71785 * x1 ** 4),
        (4 * x2 ** 2 - x1 * x2) / (12566 * (x2 * x1 ** 3 - x1 ** 4)) + 1 / (5108 * x1 ** 2) - 1,
        1 - 140.45 * x1 / (x2 ** 2 * x3),
        (x1 + x2) / 1.5 - 1,
    ]
    return obj, g


def spring_design(x):
    obj, g = spring_design_raw(x)
    return obj + _penalty(g)


SPRING_DESIGN = dict(
    name="Tension/Compression Spring Design", dim=3,
    lb=np.array([0.05, 0.25, 2.0]),
    ub=np.array([2.0, 1.3, 15.0]),
    fitness=spring_design, raw=spring_design_raw, best_known=0.012665,
)


# ------------------------------------------------------------- Speed reducer
def speed_reducer_raw(x):
    x1, x2, x3, x4, x5, x6, x7 = x
    x3 = round(x3)  # number of teeth is integer

    obj = (0.7854 * x1 * x2 ** 2 * (3.3333 * x3 ** 2 + 14.9334 * x3 - 43.0934)
           - 1.508 * x1 * (x6 ** 2 + x7 ** 2)
           + 7.4777 * (x6 ** 3 + x7 ** 3)
           + 0.7854 * (x4 * x6 ** 2 + x5 * x7 ** 2))

    g = [
        27 / (x1 * x2 ** 2 * x3) - 1,
        397.5 / (x1 * x2 ** 2 * x3 ** 2) - 1,
        1.93 * x4 ** 3 / (x2 * x3 * x6 ** 4) - 1,
        1.93 * x5 ** 3 / (x2 * x3 * x7 ** 4) - 1,
        np.sqrt((745 * x4 / (x2 * x3)) ** 2 + 16.9e6) / (110 * x6 ** 3) - 1,
        np.sqrt((745 * x5 / (x2 * x3)) ** 2 + 157.5e6) / (85 * x7 ** 3) - 1,
        x2 * x3 / 40 - 1,
        5 * x2 / x1 - 1,
        x1 / (12 * x2) - 1,
        (1.5 * x6 + 1.9) / x4 - 1,
        (1.1 * x7 + 1.9) / x5 - 1,
    ]
    return obj, g


def speed_reducer(x):
    obj, g = speed_reducer_raw(x)
    return obj + _penalty(g)


SPEED_REDUCER = dict(
    name="Speed Reducer Design", dim=7,
    lb=np.array([2.6, 0.7, 17.0, 7.3, 7.8, 2.9, 5.0]),
    ub=np.array([3.6, 0.8, 28.0, 8.3, 8.3, 3.9, 5.5]),
    fitness=speed_reducer, raw=speed_reducer_raw, best_known=2994.4711,
)


PROBLEMS = {
    "PressureVessel": PRESSURE_VESSEL,
    "WeldedBeam": WELDED_BEAM,
    "SpringDesign": SPRING_DESIGN,
    "SpeedReducer": SPEED_REDUCER,
}
