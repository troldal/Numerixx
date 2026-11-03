# Globally adaptive G(n)-K(2n+1) prototype for nxx::quadrature Phase 2 (BRIEF.md 2.3, 3.2, 3.3, 3.6, 3.7).
# One loop, generic over the arithmetic (IEEE double through Python floats, or mpmath mpf at a fixed precision).
# Mirrors Phase 1's detail::Interval (point clamp, per-term scale order) and detail::Sum (Neumaier).
import heapq, math, sys
import mpmath as mp
import rules

DBL_MAX = sys.float_info.max
sys.setrecursionlimit(20000)
TRACE_K = (4, 8, 16, 32, 64, 128)
RESET_R = (2, 10, 100)


class Double:
    name = 'double'
    eps = 2.0 ** -52
    dmin = 5e-324
    dmax = DBL_MAX
    digits = 53
    min_exp = -1021  # numeric_limits<double>::min_exponent
    max_exp = 1024
    isfinite = staticmethod(math.isfinite)
    num = staticmethod(float)
    zero = 0.0


class MPF:
    """mpmath at a fixed binary precision; mpmath has no overflow or subnormals, so dmin / dmax are symbolic."""
    def __init__(self, prec=168):
        self.name = f'mpf{prec}'
        self.prec = prec
        self.digits = prec
        with mp.workprec(prec):
            self.eps = mp.mpf(2) ** (1 - prec)
            self.dmin = mp.mpf(2) ** -100000
            self.dmax = mp.mpf(2) ** 100000
        self.min_exp, self.max_exp = -100000, 100000
        self.zero = mp.mpf(0)
    def isfinite(self, x): return mp.isfinite(x)
    def num(self, x): return mp.mpf(x)


class Sum:  # Neumaier, as Phase 1's detail::Sum
    __slots__ = ('s', 'c')
    def __init__(s, z): s.s = z; s.c = z
    def add(s, x):
        t = s.s + x
        if abs(s.s) >= abs(x): s.c += (s.s - t) + x
        else: s.c += (x - t) + s.s
        s.s = t
    def get(s): return s.s + s.c


class Interval:
    """Phase 1's detail::Interval: nodes c + hw u clamped into [lo, hi]; term(w, y) = hw w y in the safe order."""
    __slots__ = ('a', 'b', 'lo', 'hi', 'c', 'hw', 'width', 'finite', 'wide')
    def __init__(iv, ctx, a, b):
        iv.a, iv.b = a, b
        iv.lo, iv.hi = (a, b) if a <= b else (b, a)
        iv.c = a / 2 + b / 2
        iv.hw = b / 2 - a / 2
        # whether b - a overflows is decided before forming it (Phase 1, finding 3)
        iv.finite = not ((a < 0) != (b < 0) and abs(a) / 2 + abs(b) / 2 > ctx.dmax / 2)
        iv.width = b - a if iv.finite else None
        iv.wide = iv.finite and abs(iv.width) >= 1
    def point(iv, u):
        x = iv.c + iv.hw * u
        return iv.lo if x < iv.lo else (iv.hi if iv.hi < x else x)
    def raw(iv, u):
        return iv.c + iv.hw * u
    def term(iv, w, y):
        if not iv.finite:
            return (iv.hw * y) * w if abs(y) <= 1 else iv.hw * (w * y)
        half = w / 2
        return (iv.width * half) * y if iv.wide else iv.width * (half * y)


class Rule:
    def __init__(self, ctx, n=7, kind='double', force=True):
        """kind: 'double' = the formulas evaluated in double (emulates the plain-T generator);
        'rounded' = correctly rounded from the 120-digit reference; for MPF always rounded at ctx precision."""
        ref = rules.ref_rule(n)
        if isinstance(ctx, MPF):
            with mp.workprec(ctx.prec):
                r = [(+ref['x'][i], +ref['wk'][i], +ref['wg'][i], +ref['d'][i]) for i in range(len(ref['x']))]
                if force:
                    ci = len(r) // 2; s = mp.mpf(0)
                    for i in range(len(r)):
                        if i != ci: s += r[i][3]
                    r[ci] = (r[ci][0], r[ci][1], r[ci][2], -s)
        else:
            r = rules.double_rule(n, ref) if kind == 'double' else rules.rounded_rule(ref)
            if force: r = rules.force_centre(r)
        self.n = n
        self.x = [t[0] for t in r]
        self.wk = [t[1] for t in r]
        self.wg = [t[2] for t in r]
        self.d = [t[3] for t in r]
        self.gauss = [i for i in range(len(r)) if self.wg[i] != 0]
        self.xmax = self.x[-1]
        self.m = len(r)


class Result:
    __slots__ = ('status', 'value', 'abs_error', 'evals', 'panels', 'bisections', 'max_run', 'max_run_prev',
                 'heap_max', 'final_unsplit', 'd14', 'max_ratio', 'reason', 'trips', 'trips_prev', 'undone', 'trips_r')
    def __init__(r, **k):
        for s in Result.__slots__: setattr(r, s, k.get(s))
    def __repr__(r):
        return (f"{r.status:13s} v={r.value!r} e={r.abs_error!r} evals={r.evals} panels={r.panels} "
                f"bis={r.bisections} run={r.max_run} reason={r.reason}")


def integrate(f, a, b, ctx=Double, rule=None, abs_tol=None, rel_tol=None, mesh=10, stag_k=None, stag_rule='min',
              cap=None, floor_c=50, diffw=True, panels_out=None, d14_undo=True, stag_reset=None, tol_floor_mult=1):
    """Returns Result. stag_k None: no stagnation rule. cap: max bisections (None: unlimited)."""
    if rule is None: rule = Rule(ctx)
    Z = ctx.zero
    eps, dmin = ctx.eps, ctx.dmin
    if abs_tol is None: abs_tol = Z
    if rel_tol is None: rel_tol = ctx.num(math.sqrt(ctx.eps)) if ctx is Double else mp.sqrt(eps)
    n_ev = [0]
    def fn(x):
        n_ev[0] += 1
        try:
            return f(x)
        except (ZeroDivisionError, ValueError, OverflowError):
            return math.inf if ctx is Double else mp.inf
    # A. validation
    if not (ctx.isfinite(abs_tol) and ctx.isfinite(rel_tol)) or abs_tol < 0 or rel_tol < 0:
        return Result(status='non_finite', reason='tolerance', evals=0)
    if rel_tol > 1: rel_tol = ctx.num(1)
    if not (ctx.isfinite(a) and ctx.isfinite(b)): return Result(status='non_finite', reason='end', evals=0)
    if a == b: return Result(status='ok', value=Z, abs_error=Z, evals=0, panels=0, bisections=0)
    sign = 1
    if b < a: a, b, sign = b, a, -1  # orientation: integrate over [min, max] and negate (a choice, see report)

    xs, wk, wg, dw, gauss = rule.x, rule.wk, rule.wg, rule.d, rule.gauss
    xmax = rule.xmax
    m = rule.m

    def splittable(lo, hi):
        iv = Interval(ctx, lo, hi)
        mid = iv.c
        if not (lo < mid < hi): return False
        for p, q in ((lo, mid), (mid, hi)):
            ch = Interval(ctx, p, q)
            if not (p < ch.raw(-xmax) and ch.raw(xmax) < q): return False
        return True

    stats = dict(d14=0, unsplit=0, max_ratio=0.0, undone=0)

    def evaluate(lo, hi):
        """Returns (K, err, floor, L1) or None if L1 is not finite (some f_i not finite, or overflow)."""
        iv = Interval(ctx, lo, hi)
        K = Sum(Z); D = Sum(Z); L1 = Z
        G = None if diffw else Sum(Z)
        for i in range(m):
            y = fn(iv.point(xs[i]))
            t = iv.term(wk[i], y)
            K.add(t); L1 += abs(t)
            if diffw: D.add(iv.term(dw[i], y))
            elif wg[i] != 0: G.add(iv.term(wg[i], y))
        if not ctx.isfinite(L1): return None
        Kv = K.get()
        Dv = abs(D.get()) if diffw else abs(Kv - G.get())
        if not (ctx.isfinite(Kv) and ctx.isfinite(Dv)): return 'overflow'
        floor = floor_c * max(eps * L1, dmin)
        if L1 > 0:
            r = float(Dv / (eps * L1)) if eps * L1 > 0 else 0.0
            stats['max_ratio'] = max(stats['max_ratio'], r)
        return Kv, Dv, floor, L1

    heap = []   # (-err, seq, lo, hi, K, err, floor)
    final = []  # (lo, hi, K, err, floor)
    seq = [0]
    Ksum, Esum, Fsum, Xsum = Sum(Z), Sum(Z), Sum(Z), Sum(Z)

    def collect(lo, hi, out):
        """Evaluate [lo, hi] into out; D14: a non-finite value on a splittable panel bisects it. Returns status:
        'ok', 'non_finite' (overflow with finite fn, or D14 reached an unsplittable panel)."""
        r = evaluate(lo, hi)
        if r == 'overflow': return 'non_finite'
        if r is None:
            if splittable(lo, hi):
                stats['d14'] += 1
                mid = Interval(ctx, lo, hi).c
                s1 = collect(lo, mid, out)
                if s1 != 'ok': return s1
                return collect(mid, hi, out)
            return 'non_finite_fn'
        out.append((lo, hi) + r)
        return 'ok'

    def commit(out):
        for lo, hi, K, err, floor, L1 in out:
            e = max(err, floor)
            Ksum.add(K); Esum.add(e); Fsum.add(floor); Xsum.add(e - floor)
            if err <= floor:
                final.append((lo, hi, K, e, floor))
            elif not splittable(lo, hi):
                stats['unsplit'] += 1
                final.append((lo, hi, K, e, floor))
            else:
                seq[0] += 1
                heapq.heappush(heap, (-err, seq[0], lo, hi, K, e, floor))

    def add_panel(lo, hi):
        out = []
        s = collect(lo, hi, out)
        if s != 'ok': return 'non_finite'
        commit(out)
        return 'ok'

    # initial mesh
    G0 = Interval(ctx, a, b)
    pts = [a]
    for i in range(1, mesh):
        u = ctx.num(-1) + ctx.num(2 * i) / mesh
        pts.append(G0.point(u))
    pts.append(b)
    for i in range(mesh):
        if pts[i] < pts[i + 1]:
            s = add_panel(pts[i], pts[i + 1])
            if s != 'ok': return Result(status=s, reason='fn or overflow (initial)', evals=n_ev[0])

    bis = 0; run = 0; run_prev = 0; max_run = 0; max_run_prev = 0
    Emin = None; Eprev = None; heap_max = len(heap)
    status = None; reason = None
    trips = {}; trips_prev = {}; trips_r = {}; rstate = {R: [None, 0] for R in RESET_R + tuple(-r for r in RESET_R)}
    while True:
        I = Ksum.get(); E = Esum.get(); F = Fsum.get()
        if not (ctx.isfinite(I) and ctx.isfinite(E) and ctx.isfinite(F)):
            status, reason = 'non_finite', 'sum overflow'; break
        tol = max(abs_tol, rel_tol * abs(I), tol_floor_mult * F)
        if E <= tol: status = 'ok'; break
        if not heap: status, reason = 'not_converged', 'final panels'; break
        if cap is not None and bis >= cap: status, reason = 'not_converged', 'cap'; break
        if stag_k is not None:
            cnt = run if stag_rule == 'min' else run_prev
            if cnt >= stag_k: status, reason = 'not_converged', 'stagnation'; break
        _, _, lo, hi, K, e, floor = heapq.heappop(heap)
        mid = Interval(ctx, lo, hi).c
        out = []
        s = collect(lo, mid, out)
        if s == 'ok': s = collect(mid, hi, out)
        if s == 'non_finite_fn' and d14_undo:
            # the bisection is undone: the parent becomes final with its own value and estimate
            stats['undone'] += 1
            final.append((lo, hi, K, e, floor))
            continue
        if s != 'ok': status, reason = 'non_finite', 'fn or overflow'; break
        Ksum.add(-K); Esum.add(-e); Fsum.add(-floor); Xsum.add(-(e - floor))
        commit(out)
        bis += 1
        heap_max = max(heap_max, len(heap))
        En = Esum.get()
        if Emin is None or En < Emin: Emin = En; run = 0
        else: run += 1
        Xn = Xsum.get()
        for R in RESET_R + tuple(-r for r in RESET_R):
            st = rstate[R]
            V = En if R > 0 else Xn
            if st[0] is None or V < st[0] or V > abs(R) * st[0]: st[0] = V; st[1] = 0
            else: st[1] += 1
            for k in TRACE_K:
                if st[1] >= k and (R, k) not in trips_r: trips_r[(R, k)] = (bis, n_ev[0])
                if st[1] >= k and stag_k == k and stag_reset == R and status is None:
                    pass
        if stag_reset is not None: run = rstate[stag_reset][1]
        if Eprev is None or En < Eprev: run_prev = 0
        else: run_prev += 1
        Eprev = En
        max_run = max(max_run, run); max_run_prev = max(max_run_prev, run_prev)
        for k in TRACE_K:
            if run >= k and k not in trips: trips[k] = (bis, n_ev[0])
            if run_prev >= k and k not in trips_prev: trips_prev[k] = (bis, n_ev[0])

    if status == 'non_finite':
        return Result(status=status, reason=reason, evals=n_ev[0], bisections=bis, trips=trips, trips_prev=trips_prev,
                      max_run=max_run, trips_r=trips_r)
    # the result: a fresh compensated sum over every panel
    V = Sum(Z); Er = Sum(Z)
    for t in heap: V.add(t[4]); Er.add(t[5])
    for t in final: V.add(t[2]); Er.add(t[3])
    value = V.get(); err = Er.get()
    if panels_out is not None:
        panels_out.extend(sorted([(t[2], t[3], t[4], t[5], t[6], 'heap') for t in heap] + [(t[0], t[1], t[2], t[3], t[4], 'final') for t in final]))
    if not (ctx.isfinite(value) and ctx.isfinite(err)):
        return Result(status='non_finite', reason='sum overflow (final)', evals=n_ev[0], bisections=bis)
    return Result(status=status, value=sign * value, abs_error=err, evals=n_ev[0], panels=len(heap) + len(final),
                  bisections=bis, max_run=max_run, max_run_prev=max_run_prev, heap_max=heap_max,
                  final_unsplit=stats['unsplit'], d14=stats['d14'], max_ratio=stats['max_ratio'], reason=reason,
                  trips=trips, trips_prev=trips_prev, undone=stats['undone'], trips_r=trips_r)


if __name__ == '__main__':
    R = Rule(Double)
    print(integrate(math.exp, 0.0, 1.0, rule=R), float(abs(mp.e - 1 - integrate(math.exp, 0.0, 1.0, rule=R).value)))
    print(integrate(lambda x: x ** -0.9, 0.0, 1.0, rule=R, mesh=1))
    print(integrate(lambda x: math.log(abs(x)), -1.0, 1.0, rule=R, mesh=1))
    print(integrate(lambda x: 1.0, 0.0, 1.0, rule=R, rel_tol=0.0))
