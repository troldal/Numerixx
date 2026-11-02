# Double-only prototype of TanhSinh (Phase 3b): PHASE3_API.md A to E, as PHASE3B_API.md pins them (C and D),
# through the method protocol of PHASE3A_API.md (start, advance, and integrate's loop). Python floats are IEEE double
# with round to nearest, as the C++ on x86-64. The nodes are read from tab_double.txt, the C++ double table dumped by
# dump.cpp, so that the evaluation counts compare the same computation.
import math
import os
import sys

EPS = 2.0 ** -52
DMIN = 5e-324
MAX = sys.float_info.max
HALF_MAX = MAX / 2
CAP = 10       # levels 0 to 9
MIN_LEVEL = 4  # ready from level 4
DEFAULT_REL = 2.0 ** -26

class Table:
    def __init__(self, path):
        self.levels = [[] for _ in range(CAP)]
        for line in open(path):
            k, j, comp, w = line.split()
            self.levels[int(k)].append((int(j), float(comp), float(w)))
        for lev in self.levels:
            lev.sort()


HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = Table(os.path.join(HERE, 'tab_double.txt'))


class Sum:  # Phase 1's detail::Sum (Neumaier)
    def __init__(s):
        s.sum = 0.0
        s.correction = 0.0

    def add(s, x):
        t = s.sum + x
        if abs(s.sum) >= abs(x):
            s.correction += (s.sum - t) + x
        else:
            s.correction += (x - t) + s.sum
        s.sum = t

    def value(s):
        return s.sum + s.correction


def larger(a, b):
    return b if a < b else a


def small_product(a, b):
    m = abs(a)
    if m <= 1.0:
        return m * abs(b) <= HALF_MAX
    return abs(b) <= HALF_MAX / m


class Interval:  # Phase 1's detail::Interval, for a < b
    def __init__(s, a, b):
        s.lo, s.hi = a, b
        s.c = a / 2 + b / 2
        s.hw = b / 2 - a / 2
        s.width = b - a
        s.finite = math.isfinite(s.width)
        s.wide = s.finite and abs(s.width) >= 1.0

    def term(s, w, y):
        if not s.finite:
            if abs(y) <= 1.0:
                return (s.hw * y) * w
            return s.hw * (w * y)
        half = w / 2
        if s.wide:
            return (s.width * half) * y
        return s.width * (half * y)


class Failure(Exception):
    def __init__(s, reason, best=None):
        s.reason, s.best = reason, best


class Kept:
    __slots__ = ('j', 'point', 'value', 'size')

    def __init__(s, j, point, value, size):
        s.j, s.point, s.value, s.size = j, point, value, size


class Side:
    def __init__(s, end, centre):
        s.end = end
        s.cut = None          # nodes at or beyond it are not sampled
        s.kind = 'table'      # 'table' (uncut), 'round', 'prune', 'tail'
        s.outer = centre      # the outermost kept point (the centre counts, j = 0)
        s.below = None        # the kept point with the largest t below outer
        s.distinct = None     # the kept point with the largest t farther from the end than outer
        s.W = 0.0             # the table's weights beyond outer, every level so far
        s.history = ['table']

    def set_cut(s, j, kind):
        s.cut, s.kind = j, kind
        if s.history[-1] != kind:
            s.history.append(kind)

    def keep(s, k):
        d = lambda p: abs(s.end - p)
        if k.j > s.outer.j:
            s.below = s.outer
            if d(s.outer.point) > d(k.point):
                s.distinct = s.outer
            s.outer = k
            return
        if s.below is None or k.j > s.below.j:
            s.below = k
        if d(k.point) > d(s.outer.point) and (s.distinct is None or k.j > s.distinct.j):
            s.distinct = k


UNBOUNDED = None


def geometric_end(h, te, tp):
    if not (te < tp):
        return UNBOUNDED
    r = te / tp
    room = 1.0 - r
    scaled = h * te
    if scaled > room * HALF_MAX:
        return UNBOUNDED
    return scaled * (1.0 / room)


def power_end(side):
    if side.distinct is None:
        return UNBOUNDED
    de, dp = abs(side.end - side.outer.point), abs(side.end - side.distinct.point)
    ye, yp = abs(side.outer.value), abs(side.distinct.value)
    if not small_product(de, ye) or not small_product(dp, yp):
        return UNBOUNDED
    re, rp = de * ye, dp * yp
    if rp == 0.0:
        return 0.0 if re == 0.0 else UNBOUNDED
    rho = de / dp
    if not (re < rp):
        return UNBOUNDED
    q = re / rp
    if q >= 1.0:
        return UNBOUNDED
    num = (1.0 - rho) * q
    den = rho * (1.0 - q)
    if not (den > 0.0) or num > den * HALF_MAX:
        return UNBOUNDED
    ratio = num / den - 1.0
    factor = 0.0 if ratio < 0.0 else ratio
    if not small_product(re, factor):
        return UNBOUNDED
    return re * factor


class State:
    pass


COUNT = [0]  # the evaluations of the current run, for a failure raised inside start


def _fn(s, fn, p):
    assert s.iv.lo < p < s.iv.hi, 'fn evaluated at or outside an end'
    s.evaluations += 1
    COUNT[0] = s.evaluations
    try:
        return float(fn(p))
    except (ValueError, ZeroDivisionError, OverflowError):
        return math.nan


def _point(s, side_index, comp):
    return s.iv.hi - s.iv.hw * comp if side_index == 0 else s.iv.lo + s.iv.hw * comp


def _inside(s, p):
    return s.iv.lo < p < s.iv.hi


def _add_w(s, level):
    # W_k = W_(k-1) + level k's weights beyond the outermost kept node, added one by one in increasing t, while the
    # side is uncut or cut by rounding (a pruned or tail side never returns to rounding).
    for side in s.sides:
        if side.kind not in ('table', 'round'):
            continue
        for j, comp, w in TABLE.levels[level]:
            if j > side.outer.j:
                side.W += w


def _report(s):
    # The step of level k: the slivers, the end terms, the reported value, the floor and the estimate.
    if not (math.isfinite(s.R) and math.isfinite(s.L1)):
        raise Failure('non_finite')
    floor = 50.0 * larger(EPS * s.L1, DMIN)
    sliver, sliver_size, ends, unbounded = 0.0, 0.0, 0.0, False
    for side in s.sides:
        if side.kind == 'round':
            piece = s.iv.term(s.h * side.W, side.outer.value)
            sliver += piece
            sliver_size += abs(piece)
            end = power_end(side)
        elif side.kind == 'prune' or side.outer.size == 0.0:
            end = 0.0
        elif side.below is None:
            end = s.h * side.outer.size  # only the centre kept
        else:
            end = geometric_end(s.h, side.outer.size, side.below.size)
        if end is UNBOUNDED:
            unbounded = True
        else:
            ends += end
    total = s.R + sliver
    total_floor = floor + sliver_size
    if not (math.isfinite(total) and math.isfinite(total_floor)):
        raise Failure('non_finite')
    current = None
    if s.previous is not None:
        diff = abs(total - s.previous)
        base = larger(diff, total_floor)
        estimate = base if unbounded else base + ends
        if not math.isfinite(estimate):
            raise Failure('non_finite')
        if not unbounded:
            current = (total, larger(estimate, total_floor))
    s.step = dict(current=current, attainable=2.0 * total_floor, ready=s.k >= MIN_LEVEL,
                  exhausted=s.k == CAP - 1, level=s.k, evaluations=s.evaluations,
                  kinds=tuple(side.kind for side in s.sides), ends=None if unbounded else ends, sliver=sliver)
    s.previous = total


def start(fn, a, b):
    s = State()
    s.iv = Interval(a, b)
    s.evaluations = 0
    s.k, s.h = 0, 1.0
    s.previous = None
    m = s.iv.c
    if not _inside(s, m):  # no point of T strictly inside: nothing can be sampled
        s.step = dict(current=None, attainable=0.0, ready=False, exhausted=True, level=0, evaluations=0)
        s.sides = []
        return s
    y0 = _fn(s, fn, m)
    if not math.isfinite(y0):
        raise Failure('non_finite')
    tc = s.iv.term(1.0, y0)  # w(0) = 1
    if not math.isfinite(tc):
        raise Failure('non_finite')
    centre = Kept(0, m, y0, abs(tc))
    s.sides = [Side(s.iv.hi, centre), Side(s.iv.lo, centre)]
    # Level 0's sweep: each side outward at t = 1, 2, ..., right side first.
    found = [[], []]
    for i, side in enumerate(s.sides):
        for j, comp, w in TABLE.levels[0][1:]:
            p = _point(s, i, comp)
            if not _inside(s, p):
                side.set_cut(j, 'round')
                break
            y = _fn(s, fn, p)
            if not math.isfinite(y):
                if comp < EPS:
                    side.set_cut(j, 'tail')
                    break
                raise Failure('non_finite')
            t = s.iv.term(w, y)
            if not math.isfinite(t):
                raise Failure('non_finite')
            found[i].append((j, p, y, t))
    # Pruning: the threshold uses the L1 of the whole sweep, the centre included (h = 1).
    sweep = Sum()
    sweep.add(abs(tc))
    for i in range(2):
        for x in found[i]:
            sweep.add(abs(x[3]))
    l1_sweep = sweep.value()
    if not math.isfinite(l1_sweep):
        raise Failure('non_finite')
    # Each side drops its run of outermost small terms, from the outside in (D2).
    for i, side in enumerate(s.sides):
        while found[i] and abs(found[i][-1][3]) < EPS * l1_sweep:
            side.set_cut(found[i][-1][0], 'prune')
            found[i].pop()
    values, sizes = Sum(), Sum()
    values.add(1.0 * tc)
    sizes.add(1.0 * abs(tc))
    for i, side in enumerate(s.sides):
        for j, p, y, t in found[i]:
            values.add(1.0 * t)
            sizes.add(1.0 * abs(t))
            side.keep(Kept(j, p, y, abs(t)))
    s.R, s.L1 = values.value(), sizes.value()
    _add_w(s, 0)
    _report(s)
    return s


def advance(fn, s):
    s.k += 1
    s.h = s.h / 2
    l1_prev = s.L1
    threshold = (EPS * l1_prev) / (2.0 * s.h)  # |h term| < eps L1_(k-1) / 2, in term scale (D2)
    new_values, new_sizes = Sum(), Sum()
    for i, side in enumerate(s.sides):
        for j, comp, w in TABLE.levels[s.k]:
            if side.cut is not None and j >= side.cut:
                break
            p = _point(s, i, comp)
            outer = j > side.outer.j
            if not _inside(s, p):
                side.set_cut(j, 'round')
                break
            y = _fn(s, fn, p)
            if not math.isfinite(y):
                if outer and comp < EPS:
                    side.set_cut(j, 'tail')
                    break
                raise Failure('non_finite')
            t = s.iv.term(w, y)
            if not math.isfinite(t):
                raise Failure('non_finite')
            if outer and abs(t) < threshold:
                side.set_cut(j, 'prune')
                break
            new_values.add(s.h * t)
            new_sizes.add(s.h * abs(t))
            side.keep(Kept(j, p, y, abs(t)))
    s.R = s.R / 2 + new_values.value()
    s.L1 = s.L1 / 2 + new_sizes.value()
    _add_w(s, s.k)
    _report(s)


def integrate(fn, a, b, rel_tol=None, abs_tol=0.0, observer=None):
    """Phase 3a's integrate on TanhSinh: validation, the ends in order, the loop of B. Returns a dict."""
    rel_tol = DEFAULT_REL if rel_tol is None else min(rel_tol, 1.0)
    if a == b:
        return dict(status='ok', value=0.0, abs_error=0.0, evaluations=0, steps=[])
    sign = 1.0
    if b < a:
        a, b, sign = b, a, -1.0
    steps = []
    s = None
    COUNT[0] = 0
    try:
        s = start(fn, a, b)
        while True:
            steps.append(s.step)
            if observer:
                observer(s.step)
            cur = s.step['current']
            passed = s.step['ready'] and cur is not None and \
                cur[1] <= larger(larger(abs_tol, rel_tol * abs(cur[0])), s.step['attainable'])
            if passed:
                return dict(status='ok', value=sign * cur[0], abs_error=cur[1], evaluations=s.evaluations,
                            steps=steps, state=s)
            if s.step['exhausted']:
                best = None if cur is None else (sign * cur[0], cur[1])
                return dict(status='not_converged', best=best, evaluations=s.evaluations, steps=steps, state=s)
            advance(fn, s)
    except Failure as f:
        return dict(status=f.reason, best=None, evaluations=COUNT[0], steps=steps, state=s)


def show(name, r, exact=None):
    if r['status'] == 'ok':
        v, e = r['value'], r['abs_error']
        out = f"ok {v!r} +- {e:.3g}"
    elif r.get('best'):
        v, e = r['best']
        out = f"{r['status']}, best {v!r} +- {e:.3g}"
    else:
        v = e = None
        out = f"{r['status']}, no best"
    err = ''
    if exact is not None and v is not None:
        err = f"  error {abs(v - exact):.3g}" + ('  FALSE SUCCESS' if r['status'] == 'ok' and abs(v - exact) > e
                                                 else ('  (estimate below error)' if abs(v - exact) > e else ''))
    s = r.get('state')
    kinds = '' if s is None or not getattr(s, 'sides', None) else \
        '  sides ' + ' / '.join('>'.join(side.history) for side in s.sides)
    print(f"{name:44s} {out}  [{r['evaluations']} evaluations]{err}{kinds}")
