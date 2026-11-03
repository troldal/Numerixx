# Prototype of TanhSinh on infinite intervals (nxx::quadrature Phase 4, PHASE4_API.md): exp-sinh on a half line,
# sinh-sinh on the line, on the machinery of Phase 3b's ts.py (beside this file) unchanged where it carries over:
# levels (t = j / 512), the level-0 sweep and pruning from the outside in, the level-k pruning threshold
# eps L1_(k-1) / 2, cuts by kind (table, round, prune, tail), the end terms by kind, the slivers and the power-law
# excess, the raw sum, the floor 50 max(eps L1, denorm_min), the estimate and the protocol loop (ready from level 4,
# 10 levels). Generic over the arithmetic: IEEE double through Python floats, or mpmath at 113 / 168 bits with eps,
# min, denorm_min and max emulated (fn's operations overflow to inf and underflow to 0 through Ops).
#
# Maps (c > 0), nodes from exp alone (e = exp(t), sinh t and cosh t by series below t = 1/2, E = exp(c sinh t)):
#   'es' [a, inf):  x = a + E (side 0, toward inf, weight c cosh t E) and x = a + 1/E (side 1, toward a, weight
#                   c cosh t / E); the centre a + 1, weight c. (-inf, b] is fn(-u) on [-b, inf), which evaluates fn
#                   at b - E and b - 1/E exactly (round to nearest is symmetric), as the C++ does directly.
#   'ss' (-inf, inf): x = +-sinh(c sinh t) (sinh, cosh of u from E and 1/E, series below u = 1/2), weight
#                   c cosh t cosh(c sinh t); the centre 0, weight c.
#   'ts' [a, b]:    Phase 3's map (c = 1), from tab_double.txt in double (bit for bit ts.py) or generated.
# A side toward infinity ends at the last node with E and weight <= max/2, decided before dividing (c cosh t <=
# (max/2) q; the line's node exists where that exp-sinh entry does); one toward a finite end at the last node
# whose distance (comp, 1/E) is >= min(). A tail toward a finite end: distance < eps, as Phase 3.
#
# The default is the approved design (PHASE4_API.md draft 3, Q1 to Q8):
#   c = 3/2 on half lines and the line (finite intervals keep c = 1 and Phase 3's rules);
#   form='one_q': q = exp(-u) on exp's present domain, E = 1/q, and the line's node from the exp-sinh entries at the
#     same t (x = (E - q)/2, sinh's series below u = 1/2; w = (w_right + w_left)/2): one table, exp not widened;
#   uf_fix: toward infinity a node is pruned only if w max(|y|, min()) is below the threshold, and a node with
#     |y| < min() that is not so ends the side as a table end (geometric term over the last two kept values), once
#     the side (the centre included) has kept a normal value (uf_fix='r2': without that condition, round 2);
#   tail_inf='wide': toward infinity a value not finite is a tail where E (sinh u) > 1/eps, or where the last kept
#     term before it is already below the pruning threshold (level 0 decides after both sweeps); else non_finite.
# cpp=True (double, c = 3/2, form 'one_q' only): the exp-sinh and sinh-sinh nodes from es_double.txt, the C++ table
#   dumped by dump.cpp (lines "level j q f [E x]"), with the weights formed as the C++ forms them at use (f E, f q,
#   (f E + f q) / 2), so that evaluation counts compare the same computation; without it Python's exp builds the
#   table, an ulp or two from the C++ one.
# Round 1: c=1, form='base', uf_fix=False, tail_inf='inv_eps'. Round 2's other variants: osc=True|'sticky',
# tail_inf='any'|'never', inf_end='zero'|'last', form='q'|'one'|'one_raw', c='pi/2' or 'n/d', tcap, fcap, nlev.
import math
import os
import sys
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
NLEV = 10


# ---------------------------------------------------------------- arithmetic
class Double:
    name = 'double'
    prec = 53
    mp_ = False
    eps = 2.0 ** -52
    nmin = 2.0 ** -1022
    dmin = 5e-324
    max = sys.float_info.max
    half_max = max / 2
    zero, one, inf, nan = 0.0, 1.0, math.inf, math.nan
    ln_max = math.log(sys.float_info.max)
    default_rel = 2.0 ** -26

    @staticmethod
    def num(x): return float(x)
    @staticmethod
    def isfinite(x): return math.isfinite(x)
    @staticmethod
    def exp(x): return math.exp(x)
    @staticmethod
    def pi(): return math.pi
    def enter(self): pass


class MP:
    def __init__(s, prec, emin, emax, subnormal, name):
        s.name, s.prec, s.mp_ = name, prec, True
        with mp.workprec(prec):
            s.eps = mp.mpf(2) ** (1 - prec)
            s.nmin = mp.mpf(2) ** emin
            s.dmin = s.nmin * s.eps if subnormal else s.nmin
            s.max = (2 - s.eps) * mp.mpf(2) ** emax
            s.half_max = s.max / 2
            s.zero, s.one = mp.mpf(0), mp.mpf(1)
            s.inf, s.nan = mp.inf, mp.nan
            s.ln_max = mp.log(s.max)
            s.ln_dmin = mp.log(s.dmin)
            s.lg_max = mp.log(s.max, 2)
            s.lg_dmin = mp.log(s.dmin, 2)
            s.default_rel = mp.sqrt(s.eps)
    def enter(s): mp.mp.prec = s.prec
    def num(s, x): return mp.mpf(x)
    def isfinite(s, x): return mp.isfinite(x) and abs(x) <= s.max
    def exp(s, x): return mp.exp(x)
    def pi(s): return +mp.pi


DOUBLE = Double()
QUAD = MP(113, -16382, 16383, True, 'binary128')
FLOAT = MP(24, -126, 127, True, 'float')
BF50 = MP(168, -2147483312, 2147483311, False, 'mp168')


class Ops:
    """fn's arithmetic as a C++ fn would see it: overflow to inf, underflow to 0 (mp emulated)."""
    def __init__(s, A):
        s.A = A

    def clamp(s, v):
        A = s.A
        if not A.mp_:
            return v
        if mp.isnan(v):
            return v
        if abs(v) > A.max:
            return A.inf if v > 0 else -A.inf
        if v != 0 and abs(v) < A.dmin:
            return A.zero
        return v

    def exp(s, z):
        A = s.A
        if not A.mp_:
            try:
                return math.exp(z)
            except OverflowError:
                return math.inf
        if mp.isnan(z): return z
        if z > A.ln_max: return A.inf
        if z < A.ln_dmin: return A.zero
        return s.clamp(mp.exp(z))

    def pow(s, x, k):  # x >= 0
        A = s.A
        if not A.mp_:
            try:
                return x ** k
            except (OverflowError, ZeroDivisionError):
                return math.inf
        if x == 0:
            return A.inf if k < 0 else A.zero
        if mp.isinf(x):
            return A.inf if k > 0 else A.zero
        lg = k * mp.log(abs(x), 2)
        if lg > A.lg_max: return A.inf
        if lg < A.lg_dmin: return A.zero
        return s.clamp(x ** k)

    def mul(s, a, b):
        return s.clamp(a * b)

    def log(s, x):
        return math.log(x) if not s.A.mp_ else mp.log(x)

    def sqrt(s, x):
        return math.sqrt(x) if not s.A.mp_ else mp.sqrt(x)

    def cos(s, x):
        return math.cos(x) if not s.A.mp_ else mp.cos(x)

    def sin(s, x):
        return math.sin(x) if not s.A.mp_ else mp.sin(x)

    def cosh(s, x):
        A = s.A
        if not A.mp_:
            try:
                return math.cosh(x)
            except OverflowError:
                return math.inf
        if abs(x) > A.ln_max + 0.69: return A.inf
        return s.clamp(mp.cosh(x))


# ---------------------------------------------------------------- tables
def sinh_cosh(A, t):
    """sinh t, cosh t from exp alone; series below 1/2 (e - 1/e cancels there)."""
    if t < 0.5:
        t2 = t * t
        sh, term, k = t, t, 1
        thr = A.eps / 8
        while True:
            term = term * t2 / ((2 * k) * (2 * k + 1))
            if term <= thr * sh:
                break
            sh += term
            k += 1
        ch, term, k = A.one, A.one, 1
        while True:
            term = term * t2 / ((2 * k - 1) * (2 * k))
            if term <= thr * ch:
                break
            ch += term
            k += 1
        return sh, ch
    e = A.exp(t)
    ie = A.one / e
    return (e - ie) / 2, (e + ie) / 2


class Tables:
    """levels[name][k] = [(j, d, w)], j in units of 2^-9, sorted; d the table's value (comp, E, 1/E, sinh u)."""
    def __init__(s, A, kind, c, tcap=None, fcap=None, nlev=NLEV, form='base', cpp=False):
        # form: 'base' E = exp(u) and 1/E; 'q' q = exp(-u) on exp's present domain (u <= -ln min), E = 1/q (a side
        # toward infinity then ends where q < min); for 'ss' also 'one' (the exp-sinh entries at the same t:
        # x = (E - 1/E)/2, sinh's series below u = 1/2, w = (w_right + w_left)/2), 'one_raw' (no series) and
        # 'one_q' ('one' on the 'q' entries).
        s.A, s.kind, s.nlev, s.form = A, kind, nlev, form
        unit = 2 ** (nlev - 1)
        A.enter()
        if not isinstance(c, str):
            cc = A.num(c)
        elif c == 'pi/2':
            cc = A.pi() / 2
        else:
            n_, d_ = c.split('/')
            cc = A.num(int(n_)) / int(d_)
        s.c = cc
        s.names = {'ts': ('comp',), 'es': ('E', 'iE'), 'ss': ('S',)}[kind]
        s.levels = {n: [[] for _ in range(nlev)] for n in s.names}
        s.tmax = {}
        if kind == 'ts' and A is DOUBLE and c == 1 and tcap is None and fcap is None and nlev == NLEV:
            for line in open(os.path.join(HERE, 'tab_double.txt')):
                k, j, comp, w = line.split()
                if int(j) > 0:
                    s.levels['comp'][int(k)].append((int(j), float(comp), float(w)))
            s.tmax['comp'] = max(x[0] for lev in s.levels['comp'] for x in lev) / 512
            return
        if cpp:
            if not (A is DOUBLE and kind in ('es', 'ss') and c == '3/2' and form == ('q' if kind == 'es' else 'one_q')
                    and tcap is None and fcap is None and nlev == NLEV):
                raise ValueError('cpp=True: double, c = 3/2 and the default form only')
            s.load_cpp()
            return
        alive = {n: True for n in s.names}
        j = 1
        ln_inv_nmin = -math.log(float(A.nmin)) if not A.mp_ else -mp.log(A.nmin)
        while any(alive.values()):
            k = 0 if j % unit == 0 else nlev - 1 - ((j & -j).bit_length() - 1)
            t = A.num(j) / unit
            nodes = s.node(t, cc, ln_inv_nmin)
            for n in s.names:
                if not alive[n]:
                    continue
                cap = tcap if n in ('E', 'S') else fcap
                if nodes[n] is None or (cap is not None and t > cap):
                    alive[n] = False
                    continue
                s.levels[n][k].append((j, nodes[n][0], nodes[n][1]))
                s.tmax[n] = float(t)
            j += 1

    def load_cpp(s):
        """The C++ double table, es_double.txt: "level j q f" toward the finite end only, "level j q f E x" toward
        infinity and on the line too; t = 0 (j = 0) is the centre, not a side's node."""
        for line in open(os.path.join(HERE, 'es_double.txt')):
            v = line.split()
            k, j = int(v[0]), int(v[1])
            if j == 0:
                continue
            q, f = float(v[2]), float(v[3])
            outward = len(v) == 6
            if s.kind == 'es':
                s.levels['iE'][k].append((j, q, f * q))
                if outward:
                    E = float(v[4])
                    s.levels['E'][k].append((j, E, f * E))
            elif outward:
                E, x = float(v[4]), float(v[5])
                s.levels['S'][k].append((j, x, (f * E + f * q) / 2))
        for n in s.names:
            for lev in s.levels[n]:
                lev.sort()
            s.tmax[n] = max(x[0] for lev in s.levels[n] for x in lev) / 512

    def node(s, t, c, ln_inv_nmin):
        A = s.A
        sh, ch = sinh_cosh(A, t)
        u = c * sh
        out = {}
        if s.kind == 'ts':
            if -2 * u < (math.log(float(A.nmin)) if not A.mp_ else mp.log(A.nmin)):
                return {'comp': None}
            q = A.exp(-2 * u)
            d = A.one + q
            out['comp'] = None if q < A.nmin else (2 * q / d, c * ch * (4 * q / (d * d)))
            return out
        if s.form != 'base':
            return s.node_q(c, ln_inv_nmin, ch, u)
        E = A.exp(u) if u <= A.ln_max else None
        if s.kind == 'es':
            if E is None or not A.isfinite(E) or not (c * ch <= A.half_max / E):
                out['E'] = None
            else:
                out['E'] = (E, c * ch * E)
            if u > ln_inv_nmin:
                out['iE'] = None
            else:
                iE = A.one / E
                out['iE'] = (iE, c * ch * iE) if iE >= A.nmin else None
            return out
        # ss
        if E is None or not A.isfinite(E) or not (c * ch <= A.half_max / E):
            return {'S': None}
        iE = A.one / E
        if u < 0.5:
            S, C = sinh_cosh(A, u)
        else:
            S, C = (E - iE) / 2, (E + iE) / 2
        w = c * ch * C
        return {'S': (S, w)}

    def node_q(s, c, ln_inv_nmin, ch, u):
        A = s.A
        none = {'E': None, 'iE': None, 'S': None}
        if s.form in ('q', 'one_q'):
            if u > ln_inv_nmin:
                return none
            q = A.exp(-u)
            if q < A.nmin:
                return none
            iE = q
            toward_inf = c * ch <= A.half_max * q   # decided before dividing: E = 1/q, wE <= max/2
            E = A.one / q if toward_inf else None
        else:
            if u > A.ln_max:
                return none
            E = A.exp(u)
            if not A.isfinite(E):
                return none
            iE = A.one / E
            toward_inf = c * ch <= A.half_max / E
        wiE = c * ch * iE
        wE = c * ch * E if toward_inf else None
        if s.kind == 'es':
            return {'E': (E, wE) if toward_inf else None,
                    'iE': (iE, wiE) if iE >= A.nmin else None}
        if not toward_inf:
            return none
        if s.form == 'q':
            if u < 0.5:
                S, C = sinh_cosh(A, u)
            else:
                S, C = (E - iE) / 2, (E + iE) / 2
            w = c * ch * C
            return {'S': (S, w) if A.isfinite(w) else None}
        S = sinh_cosh(A, u)[0] if (u < 0.5 and s.form != 'one_raw') else (E - iE) / 2
        return {'S': (S, (wE + wiE) / 2)}

    def count(s):
        return {n: sum(len(l) for l in s.levels[n]) for n in s.names}


_TABLES = {}


def tables(A, kind, c, tcap=None, fcap=None, nlev=NLEV, form='base', cpp=False):
    if kind == 'es' and form.startswith('one'):
        form = 'q' if form == 'one_q' else 'base'
    if kind == 'ts':
        form = 'base'
        cpp = False  # the finite table is tab_double.txt in double anyway
    key = (A.name, kind, str(c), tcap, fcap, nlev, form, cpp)
    if key not in _TABLES:
        _TABLES[key] = Tables(A, kind, c, tcap, fcap, nlev, form, cpp)
    return _TABLES[key]


# ---------------------------------------------------------------- the machinery of ts.py, generic
class Sum:
    def __init__(s, z):
        s.sum = z
        s.correction = z

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


class Failure(Exception):
    def __init__(s, reason, best=None):
        s.reason, s.best = reason, best


class Kept:
    __slots__ = ('j', 'point', 'value', 'size')

    def __init__(s, j, point, value, size):
        s.j, s.point, s.value, s.size = j, point, value, size


class Side:
    def __init__(s, levels, end, centre, point, infinite):
        s.levels = levels          # the side's table, per level
        s.end = end                # the finite end, or None toward infinity
        s.point = point            # d -> x
        s.infinite = infinite
        s.cut = None
        s.kind = 'table'
        s.outer = centre
        s.below = None
        s.distinct = None
        s.W = None
        s.history = ['table']

    def set_cut(s, j, kind):
        s.cut, s.kind = j, kind
        if s.history[-1] != kind:
            s.history.append(kind)

    def keep(s, k):
        if s.end is None:
            if k.j > s.outer.j:
                s.below = s.outer
                s.outer = k
            elif s.below is None or k.j > s.below.j:
                s.below = k
            return
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


class Engine:
    """One run's arithmetic and options."""
    def __init__(s, A, opts):
        s.A = A
        s.min_level = opts.get('min_level', 4)
        s.cap = opts.get('nlev', NLEV)
        s.tail_inf = opts.get('tail_inf', 'wide')   # 'wide', 'inv_eps' (E > 1/eps), 'any' (beyond kept), 'never'
        s.prune = opts.get('prune', True)
        s.inf_end = opts.get('inf_end', 'geometric')   # end term toward infinity at the table's end / a tail
        s.uf_fix = opts.get('uf_fix', True)   # toward infinity: prune on w max(|y|, min()); an underflowed value
                                               # that is not negligible so ends the side as a table end
        s.osc = opts.get('osc', False)         # toward infinity: unbounded when the last two values differ in sign

    def small_product(s, a, b):
        m = abs(a)
        if m <= 1:
            return m * abs(b) <= s.A.half_max
        return abs(b) <= s.A.half_max / m

    def geometric_end(s, h, te, tp):
        A = s.A
        if not (te < tp):
            return UNBOUNDED
        r = te / tp
        room = 1 - r
        scaled = h * te
        if scaled > room * A.half_max:
            return UNBOUNDED
        return scaled * (1 / room)

    def power_end(s, side):
        A = s.A
        if side.distinct is None:
            return UNBOUNDED
        de, dp = abs(side.end - side.outer.point), abs(side.end - side.distinct.point)
        ye, yp = abs(side.outer.value), abs(side.distinct.value)
        if not s.small_product(de, ye) or not s.small_product(dp, yp):
            return UNBOUNDED
        re, rp = de * ye, dp * yp
        lost = (re < A.nmin and ye > 0) or (rp < A.nmin and yp > 0)
        if (yp == 0) if lost else (rp == 0):
            return A.zero if (not lost and re == 0) else UNBOUNDED
        rho = de / dp
        scaled = rho * ye if lost else re
        base = yp if lost else rp
        if not (scaled < base):
            return UNBOUNDED
        q = scaled / base
        if q >= 1:
            return UNBOUNDED
        num = (1 - rho) * q
        den = rho * (1 - q)
        if not (den > 0) or num > den * A.half_max:
            return UNBOUNDED
        ratio = num / den - 1
        factor = A.zero if ratio < 0 else ratio
        if lost:
            if not s.small_product(ye, factor):
                return UNBOUNDED
            excess = ye * factor
            if not s.small_product(de, excess):
                return UNBOUNDED
            return de * excess
        if not s.small_product(re, factor):
            return UNBOUNDED
        return re * factor


class Map:
    """The interval and its map: the centre, the sides and the term."""
    def __init__(s, A, kind, a, b, tab):
        s.A, s.kind, s.a, s.b, s.tab = A, kind, a, b, tab
        if kind == 'ts':
            s.lo, s.hi = a, b
            s.m = a / 2 + b / 2
            s.hw = b / 2 - a / 2
            s.width = b - a
            s.finite = A.isfinite(s.width)
            s.wide = s.finite and abs(s.width) >= 1
            s.centre, s.w0 = s.m, A.one
            hi, lo, hw = b, a, s.hw
            s.specs = [('comp', b, lambda d: hi - hw * d, False), ('comp', a, lambda d: lo + hw * d, False)]
        elif kind == 'es':
            s.centre, s.w0 = a + 1, tab.c
            s.specs = [('E', None, lambda d: a + d, True), ('iE', a, lambda d: a + d, False)]
        else:
            s.centre, s.w0 = A.zero, tab.c
            s.specs = [('S', None, lambda d: d, True), ('S', None, lambda d: -d, True)]

    def inside(s, p):
        A = s.A
        if s.kind == 'ts':
            return s.lo < p < s.hi
        if s.kind == 'es':
            return s.a < p and A.isfinite(p)
        return A.isfinite(p)

    def term(s, w, y):
        if s.kind != 'ts':
            return w * y
        if not s.finite:
            if abs(y) <= 1:
                return (s.hw * y) * w
            return s.hw * (w * y)
        half = w / 2
        if s.wide:
            return (s.width * half) * y
        return s.width * (half * y)


class State:
    pass


def _fn(s, fn, p):
    A = s.A
    s.evaluations += 1
    s.count[0] = s.evaluations
    try:
        y = fn(p)
    except (ValueError, ZeroDivisionError, OverflowError):
        return A.nan
    if A.mp_:
        y = s.ops.clamp(mp.mpf(y))
    else:
        y = float(y)
    return y


def _finite(s, y):
    return s.A.isfinite(y)


def _is_tail(s, side, d):
    A = s.A
    if not side.infinite:
        return d < A.eps
    if s.eng.tail_inf == 'never':
        return False
    if s.eng.tail_inf == 'any':
        return True
    return d > 1 / A.eps  # 'inv_eps', and 'wide' (whose second clause is in start and advance)


def _add_w(s, level):
    for side in s.sides:
        if side.infinite or side.kind not in ('table', 'round'):
            continue
        for j, d, w in side.levels[level]:
            if j > side.outer.j:
                side.W += w


def _report(s):
    A, eng = s.A, s.eng
    if not (_finite(s, s.R) and _finite(s, s.L1)):
        raise Failure('non_finite')
    floor = 50 * larger(A.eps * s.L1, A.dmin)
    sliver, sliver_size, ends, unbounded = A.zero, A.zero, A.zero, False
    for side in s.sides:
        if side.kind == 'round':
            piece = s.map.term(s.h * side.W, side.outer.value)
            sliver += piece
            sliver_size += abs(piece)
            end = eng.power_end(side)
        elif side.kind == 'prune' or side.outer.size == 0:
            end = A.zero
        elif side.below is None:
            end = s.h * side.outer.size
        elif side.infinite and eng.osc and (side.outer.value * side.below.value < 0 or
                                            (eng.osc == 'sticky' and getattr(side, 'osc_seen', False))):
            side.osc_seen = True
            end = UNBOUNDED
        elif side.infinite and eng.inf_end == 'zero':
            end = A.zero
        elif side.infinite and eng.inf_end == 'last':
            end = s.h * side.outer.size
        else:
            end = eng.geometric_end(s.h, side.outer.size, side.below.size)
        if end is UNBOUNDED:
            unbounded = True
        else:
            ends += end
    total = s.R + sliver
    total_floor = floor + sliver_size
    if not (_finite(s, total) and _finite(s, total_floor)):
        raise Failure('non_finite')
    current = None
    if s.previous is not None:
        diff = abs(total - s.previous)
        base = larger(diff, total_floor)
        estimate = base if unbounded else base + ends
        if not _finite(s, estimate):
            raise Failure('non_finite')
        if not unbounded:
            current = (total, larger(estimate, total_floor))
    s.step = dict(current=current, attainable=2 * total_floor, ready=s.k >= eng.min_level,
                  exhausted=s.k == eng.cap - 1, level=s.k, evaluations=s.evaluations,
                  kinds=tuple(side.kind for side in s.sides), ends=None if unbounded else ends, sliver=sliver,
                  diff=None if s.previous is None else abs(total - s.previous))
    s.previous = total


def start(s, fn):
    A = s.A
    s.k, s.h = 0, A.one
    s.previous = None
    mp_ = s.map
    m = mp_.centre
    if not (mp_.inside(m) if mp_.kind != 'es' else (mp_.a < m and A.isfinite(m))):
        s.step = dict(current=None, attainable=A.zero, ready=False, exhausted=True, level=0, evaluations=0)
        s.sides = []
        return
    y0 = _fn(s, fn, m)
    if not _finite(s, y0):
        raise Failure('non_finite')
    tc = mp_.term(mp_.w0, y0)
    if not _finite(s, tc):
        raise Failure('non_finite')
    centre = Kept(0, m, y0, abs(tc))
    s.sides = []
    for name, end, point, infinite in mp_.specs:
        side = Side(s.tab.levels[name], end, centre, point, infinite)
        side.W = A.zero
        s.sides.append(side)
    found = [[] for _ in s.sides]
    pending = {}
    for i, side in enumerate(s.sides):
        for j, d, w in side.levels[0]:
            p = side.point(d)
            if not mp_.inside(p):
                side.set_cut(j, 'table' if side.infinite else 'round')
                break
            y = _fn(s, fn, p)
            if not _finite(s, y):
                if _is_tail(s, side, d):
                    side.set_cut(j, 'tail')
                    break
                if side.infinite and s.eng.tail_inf == 'wide':
                    pending[i] = j
                    break
                raise Failure('non_finite')
            t = mp_.term(w, y)
            if not _finite(s, t):
                raise Failure('non_finite')
            found[i].append((j, p, y, t, w))
    sweep = Sum(A.zero)
    sweep.add(abs(tc))
    for f in found:
        for x in f:
            sweep.add(abs(x[3]))
    l1_sweep = sweep.value()
    if not _finite(s, l1_sweep):
        raise Failure('non_finite')
    thr0 = A.eps * l1_sweep
    for i, j in pending.items():
        # the wide tail rule: the last term before the value that is not finite is already below the threshold
        last = abs(found[i][-1][3]) if found[i] else abs(tc)
        if last < thr0:
            s.sides[i].set_cut(j, 'tail')
            s.log.append(('wide0', i, 'tail'))
        else:
            s.log.append(('wide0', i, 'non_finite'))
            raise Failure('non_finite')
    size = lambda side, x: abs(x[3])
    if s.eng.uf_fix:
        size = lambda side, x: abs(x[3])  # after the truncation, |term| decides as in round 1
        for i, side in enumerate(s.sides):
            if not side.infinite:
                continue
            for idx, x in enumerate(found[i]):
                prev = found[i][idx - 1][2] if idx > 0 else y0
                if s.eng.uf_fix == 'zero':
                    trigger, evidence = x[2] == 0, prev != 0
                else:
                    trigger, evidence = abs(x[2]) < A.nmin, s.eng.uf_fix == 'r2' or abs(prev) >= A.nmin
                if trigger and x[4] * A.nmin >= thr0 and evidence:
                    side.set_cut(x[0], 'table')
                    s.log.append(('uf0', i, x[0] / 512))
                    del found[i][idx:]
                    break
    if s.eng.prune:
        for i, side in enumerate(s.sides):
            while found[i] and size(side, found[i][-1]) < thr0:
                side.set_cut(found[i][-1][0], 'prune')
                found[i].pop()
    values, sizes = Sum(A.zero), Sum(A.zero)
    values.add(tc)
    sizes.add(abs(tc))
    for i, side in enumerate(s.sides):
        for j, p, y, t, w in found[i]:
            values.add(t)
            sizes.add(abs(t))
            side.keep(Kept(j, p, y, abs(t)))
    s.R, s.L1 = values.value(), sizes.value()
    _add_w(s, 0)
    _report(s)


def advance(s, fn):
    A = s.A
    mp_ = s.map
    s.k += 1
    s.h = s.h / 2
    threshold = (A.eps * s.L1) / (2 * s.h)
    new_values, new_sizes = Sum(A.zero), Sum(A.zero)
    for side in s.sides:
        for j, d, w in side.levels[s.k]:
            if side.cut is not None and j >= side.cut:
                break
            p = side.point(d)
            outer = j > side.outer.j
            if not mp_.inside(p):
                side.set_cut(j, 'table' if side.infinite else 'round')
                break
            y = _fn(s, fn, p)
            if not _finite(s, y):
                if outer and _is_tail(s, side, d):
                    side.set_cut(j, 'tail')
                    break
                if outer and side.infinite and s.eng.tail_inf == 'wide':
                    if side.outer.size < threshold:
                        side.set_cut(j, 'tail')
                        s.log.append(('widek', s.k, 'tail'))
                        break
                    s.log.append(('widek', s.k, 'non_finite'))
                raise Failure('non_finite')
            t = mp_.term(w, y)
            if not _finite(s, t):
                raise Failure('non_finite')
            if s.eng.uf_fix == 'zero':
                uf_trigger = y == 0 and side.outer.value != 0
            else:
                uf_trigger = abs(y) < A.nmin and (s.eng.uf_fix == 'r2' or abs(side.outer.value) >= A.nmin)
            if s.eng.uf_fix and outer and side.infinite and uf_trigger:
                kind = 'prune' if w * A.nmin < threshold else 'table'
                side.set_cut(j, kind)
                s.log.append(('ufk', s.k, kind))
                break
            if s.eng.prune and outer and abs(t) < threshold:
                side.set_cut(j, 'prune')
                break
            new_values.add(s.h * t)
            new_sizes.add(s.h * abs(t))
            side.keep(Kept(j, p, y, abs(t)))
    s.R = s.R / 2 + new_values.value()
    s.L1 = s.L1 / 2 + new_sizes.value()
    _add_w(s, s.k)
    _report(s)


def integrate(fn, a, b, A=DOUBLE, rel_tol=None, abs_tol=0, c='3/2', observer=None, **opts):
    """fn(x) in A's arithmetic; a, b may be +-inf (float('inf') or mp.inf). Returns a dict as ts.py's."""
    A.enter()
    rel_tol = A.default_rel if rel_tol is None else min(A.num(rel_tol), A.one)
    abs_tol = A.num(abs_tol)
    a, b = A.num(a), A.num(b)
    if a == b:
        return dict(status='ok', value=A.zero, abs_error=A.zero, evaluations=0, steps=[])
    sign = 1
    if b < a:
        a, b, sign = b, a, -1
    ia, ib = not mp.isfinite(a), not mp.isfinite(b)
    g = fn
    if ia and ib:
        kind = 'ss'
    elif ib:
        kind = 'es'
    elif ia:
        kind = 'es'
        a, b = -b, A.inf
        g = lambda u: fn(-u)
    else:
        kind = 'ts'
    tab = tables(A, kind, 1 if kind == 'ts' else c, opts.get('tcap'), opts.get('fcap'), opts.get('nlev', NLEV),
                 opts.get('form', 'one_q'), opts.get('cpp', False))
    s = State()
    s.A = A
    s.ops = Ops(A)
    s.eng = Engine(A, opts)
    s.tab = tab
    s.map = Map(A, kind, a, b, tab)
    s.evaluations = 0
    s.count = [0]
    s.log = []
    steps = []
    try:
        start(s, g)
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
            advance(s, g)
    except Failure as f:
        return dict(status=f.reason, best=None, evaluations=s.count[0], steps=steps, state=s)
