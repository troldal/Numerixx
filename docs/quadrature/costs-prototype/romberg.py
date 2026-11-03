# Romberg prototype (Phase 1; kept as proto3.py of its session), aligned with docs/quadrature/PHASE1_API.md sections A to E
# Changes from proto.py: the scale per term, (b - a) * (f / 2^k) where b - a is finite, else hw * (f / 2^(k-1));
# T_0 = (b - a) (fa / 2) + (b - a) (fb / 2); the cap at 16 levels; tolerance validation and the clamp of rel_tol.
# References from mpmath at 60 digits. Run: python -I romberg.py
import math, sys
import mpmath as mp
mp.mp.dps = 60
EPS = 2.0 ** -52
DMIN = 5e-324
DMAX = sys.float_info.max


class Sum:  # Neumaier
    def __init__(s): s.s = 0.0; s.c = 0.0
    def add(s, x):
        t = s.s + x
        s.c += ((s.s - t) + x) if abs(s.s) >= abs(x) else ((x - t) + s.s)
        s.s = t
    def get(s): return s.s + s.c


def isfinite(x): return math.isfinite(x)


def romberg(f, a, b, abs_tol=0.0, rel_tol=math.sqrt(EPS), kmin=4, kmax=16):
    n = [0]
    def fn(x):
        n[0] += 1
        return f(x)
    if not (isfinite(abs_tol) and isfinite(rel_tol)) or abs_tol < 0 or rel_tol < 0: return ('non_finite', None, None, 0, 0)
    rel_tol = min(rel_tol, 1.0)
    if not (isfinite(a) and isfinite(b)): return ('non_finite', None, None, 0, 0)
    if a == b: return ('ok', 0.0, 0.0, 0, 0)
    # Mirrors quadrature.hpp's detail::Interval::term and detail::romberg operation for operation.
    width = b - a; finite = isfinite(width); wide = finite and abs(width) >= 1
    hw = b / 2 - a / 2; c = a / 2 + b / 2
    lo_, hi_ = min(a, b), max(a, b)
    def point(x): return lo_ if x < lo_ else (hi_ if hi_ < x else x)
    def term(w, y):  # hw w y
        if not finite: return (hw * y) * w if abs(y) <= 1 else hw * (w * y)
        half = w / 2
        return (width * half) * y if wide else width * (half * y)
    w = 1.0
    fa = fn(a); fb = fn(b)
    ta, tb = term(w, fa), term(w, fb)
    L = abs(ta) + abs(tb)
    if not isfinite(L): return ('non_finite', None, None, 0, n[0])
    T = ta + tb
    prev = [T]; best = None
    for k in range(1, kmax + 1):
        if k > 1: w = w / 2
        m = 2 ** (k - 1)
        s = Sum(); sl = 0.0
        for j in range(1, m + 1):
            u = (2 * j - 1 - m) / m
            t = term(w, fn(point(c + hw * u)))
            s.add(t); sl += abs(t)
        T = prev[0] / 2 + s.get(); L = L / 2 + sl
        if not isfinite(L): return ('non_finite', None, None, k, n[0])
        row = [T]; p4 = 1.0
        for j in range(1, k + 1):
            p4 *= 4
            row.append(row[j - 1] + (row[j - 1] - prev[j - 1]) / (p4 - 1))
        est = abs(row[k] - prev[k - 1])
        if not isfinite(est): return ('non_finite', None, None, k, n[0])
        floor = 50.0 * max(EPS * L, DMIN)
        err = max(est, floor)
        best = (row[k], err)
        prev = row
        if k >= kmin and est <= max(abs_tol, rel_tol * abs(row[k]), floor):
            return ('ok', row[k], err, k, n[0])
    return ('not_converged', best[0], best[1], kmax, n[0])


def zero_samples(x):  # (x (1 - x) prod_(j=1..15) (16 x - j))^2
    p = x * (1 - x)
    for j in range(1, 16): p *= 16 * x - j
    return p * p


pi = math.pi
cases = [
    ('e^x [0,1]', math.exp, 0.0, 1.0, {}, mp.e - 1),
    ('e^x [0,1] 1e-14', math.exp, 0.0, 1.0, {'rel_tol': 1e-14}, mp.e - 1),
    ('e^x [1,0]', math.exp, 1.0, 0.0, {}, 1 - mp.e),
    ('1/(1+x^2) [0,1]', lambda x: 1 / (1 + x * x), 0.0, 1.0, {}, mp.pi / 4),
    ('1/(1+x^2) [0,1] 1e-14', lambda x: 1 / (1 + x * x), 0.0, 1.0, {'rel_tol': 1e-14}, mp.pi / 4),
    ('sin [0,pi]', math.sin, 0.0, pi, {}, 1 - mp.cos(mp.mpf(pi))),
    ('sin^2 [0,2pi]', lambda x: math.sin(x) ** 2, 0.0, 2 * pi, {}, mp.mpf(pi) - mp.sin(4 * mp.mpf(pi)) / 4),
    ('sin^2 [0,2pi] 1e-14', lambda x: math.sin(x) ** 2, 0.0, 2 * pi, {'rel_tol': 1e-14}, mp.mpf(pi) - mp.sin(4 * mp.mpf(pi)) / 4),
    ('x^2 [0,3]', lambda x: x * x, 0.0, 3.0, {}, mp.mpf(9)),
    ('x^2 [0,3] 1e-14', lambda x: x * x, 0.0, 3.0, {'rel_tol': 1e-14}, mp.mpf(9)),
    ('x^3 [-1,1]', lambda x: x * x * x, -1.0, 1.0, {}, mp.mpf(0)),
    ('|x-1/5| [0,1]', lambda x: abs(x - 0.2), 0.0, 1.0, {}, mp.mpf(17) / 50),
    ('|x-1/5| [0,1] 1e-14', lambda x: abs(x - 0.2), 0.0, 1.0, {'rel_tol': 1e-14}, mp.mpf(17) / 50),
    ('peak [0,1]', lambda x: math.exp(-1e4 * (x - 0.5) ** 2), 0.0, 1.0, {}, mp.sqrt(mp.pi) * mp.erf(50) / 100),
    ('cos(1000x) [0,1]', lambda x: math.cos(1000 * x), 0.0, 1.0, {}, mp.sin(1000) / 1000),
    ('1e-300 e^x [0,1]', lambda x: 1e-300 * math.exp(x), 0.0, 1.0, {}, mp.mpf('1e-300') * (mp.e - 1)),
    ('sqrt x [0,1]', math.sqrt, 0.0, 1.0, {}, mp.mpf(2) / 3),
    ('sqrt x [0,1] 1e-14', math.sqrt, 0.0, 1.0, {'rel_tol': 1e-14}, mp.mpf(2) / 3),
    ('1/sqrt x [0,1]', lambda x: 1 / math.sqrt(x) if x > 0 else math.inf, 0.0, 1.0, {}, None),
    # edges
    ('1e308 [0,1e-10]', lambda x: 1e308, 0.0, 1e-10, {}, mp.mpf('1e308') * mp.mpf(1e-10)),
    ('0 [-max,max]', lambda x: 0.0, -DMAX, DMAX, {}, mp.mpf(0)),
    ('1e-300 [-max,max]', lambda x: 1e-300, -DMAX, DMAX, {}, mp.mpf('1e-300') * 2 * mp.mpf(DMAX)),
    ('1e300 [0,3dmin]', lambda x: 1e300, 0.0, 3 * DMIN, {}, mp.mpf('1e300') * 3 * mp.mpf(DMIN)),
    ('denorm_min [-max,max]', lambda x: DMIN, -DMAX, DMAX, {}, 2 * mp.mpf(DMAX) * mp.mpf(DMIN)),
    ('outside [D,3D]', lambda x: 1e300 if (x < DMIN or x > 3 * DMIN) else 0.0, DMIN, 3 * DMIN, {}, mp.mpf(0)),
    ('1e-310 [0,1e10]', lambda x: 1e-310, 0.0, 1e10, {}, mp.mpf(1e-310) * mp.mpf(1e10)),
    ('estimate overflow', lambda x: {0.0: -0.9 * DMAX, 1.0: -0.9 * DMAX, 0.5: 0.9 * DMAX}.get(x, 0.0), 0.0, 1.0, {}, None),
    ('R(2,2) overflow', lambda x: {0.0: -0.3 * DMAX, 1.0: -0.3 * DMAX, 0.5: -0.7 * DMAX, 0.25: DMAX, 0.75: DMAX}.get(x, 0.0), 0.0, 1.0, {}, None),
    ('1/(1+x^2) reversed', lambda x: 1 / (1 + x * x), 1.0, 0.0, {}, -mp.pi / 4),
    ('1/(1+x^2) rel 1e300', lambda x: 1 / (1 + x * x), 0.0, 1.0, {'rel_tol': 1e300}, mp.pi / 4),
    ('1/(1+x^2) rel 0', lambda x: 1 / (1 + x * x), 0.0, 1.0, {'rel_tol': 0.0}, mp.pi / 4),
    ('1e-320 e^x [0,1]', lambda x: 1e-320 * math.exp(x), 0.0, 1.0, {}, mp.mpf(1e-320) * (mp.e - 1)),
    ('1e308 sin x [0,2pi]', lambda x: 1e308 * math.sin(x), 0.0, 2 * pi, {}, None),
    ('zero samples [0,1]', zero_samples, 0.0, 1.0, {}, None),
    ('x^3 [-1,1] abs 0', lambda x: x * x * x, -1.0, 1.0, {'abs_tol': 0.0}, mp.mpf(0)),
    ('nan inside', lambda x: math.nan if x == 0.5 else x, 0.0, 1.0, {}, None),
    ('inf at end', lambda x: math.inf if x == 1.0 else x, 0.0, 1.0, {}, None),
]

if __name__ == '__main__':
    print(f"{'case':26s} {'status':14s} {'k':>3s} {'evals':>7s} {'estimate':>10s} {'actual':>10s} {'value':>24s}")
    for name, f, a, b, kw, ref in cases:
        st, v, e, k, n = romberg(f, a, b, **kw)
        act = abs(mp.mpf(v) - ref) if (v is not None and ref is not None) else None
        print(f"{name:26s} {st:14s} {k:3d} {n:7d} {'' if e is None else f'{e:.2e}':>10s} "
              f"{'' if act is None else mp.nstr(act, 3):>10s} {'' if v is None else repr(v):>24s}")
    # the not_converged best values, for the tests
    for name, f, a, b, kw in [('|x-1/5| 1e-14', lambda x: abs(x - 0.2), 0.0, 1.0, {'rel_tol': 1e-14}),
                              ('sqrt 1e-14', math.sqrt, 0.0, 1.0, {'rel_tol': 1e-14})]:
        print(name, romberg(f, a, b, **kw))
