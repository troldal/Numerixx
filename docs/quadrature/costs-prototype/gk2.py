# The Phase 2 loop as specified after review (PHASE2_API.md), with the review's alternatives as options:
#   descent L    D14's descent goes at most L levels below the non-finite panel that started it (None: unbounded),
#                iterative, every split counted against the cap together with the bisections.
#   nf 'undo'    (a) a bisection whose descent leaves a non-finite leaf is undone: the parent becomes final with
#                its last finite value and estimate.
#   nf 'carry'   (b) every non-finite leaf becomes final with value 0 and the estimate of its last finite ancestor
#                (the bisected parent); the finite pieces of the descent stay.
#   taint        once a run has undone or dropped such a piece it cannot end ok: not_converged with best. A mesh
#                panel with no finite ancestor is non_finite in both modes.
#   stag 'i'     reset (m = X, count = 0) when X < m or X > 10 m; else count += 1      (the spec)
#   stag 'ii'    X < m: m = X, count = 0; X > 10 m: m = X, count += 1; else count += 1   (numerics)
#   stag 'iii'   X < m: m = X, count = 0; else count += 1                                 (simplicity)
#                X = sum over panels of (e_p - floor_p); trip at count >= stag_k.
#   test         E <= max(abs_tol, rel_tol |I|, 2F); floor_p = 50 max(eps L1_p, denorm_min); mesh 10;
#                cap = 64 digits bisections plus descent splits.
#   diffw        True: D from the difference weights (centre forced, if the Rule was built so); False: D = K - G
#                from two Neumaier sums.
import heapq, math, sys
import mpmath as mp
from gk import Double, MPF, Sum, Interval, Rule


class Res:
    def __init__(r, **k): r.__dict__.update(k)
    def get(r, k, d=None): return r.__dict__.get(k, d)
    def __repr__(r):
        return ' '.join(f"{k}={v}" for k, v in r.__dict__.items() if k not in ('trace',))


def integrate(f, a, b, ctx=Double, rule=None, abs_tol=None, rel_tol=None, mesh=10, cap=None, descent=4,
              nf='undo', taint=True, stag='i', stag_k=32, floor_c=50, diffw=True, floor_mult=2):
    if rule is None: rule = Rule(ctx)
    Z = ctx.zero; eps, dmin = ctx.eps, ctx.dmin
    if cap is None: cap = 64 * ctx.digits
    if abs_tol is None: abs_tol = Z
    if rel_tol is None: rel_tol = ctx.num(math.sqrt(ctx.eps)) if ctx is Double else mp.sqrt(eps)
    nev = [0]
    INF = math.inf if ctx is Double else mp.inf
    def fn(x):
        nev[0] += 1
        try: return f(x)
        except (ZeroDivisionError, ValueError, OverflowError): return INF  # C++ gives inf or NaN: non-finite either way
    if not (ctx.isfinite(abs_tol) and ctx.isfinite(rel_tol)) or abs_tol < 0 or rel_tol < 0:
        return Res(status='non_finite', reason='tolerance', evals=0)
    if rel_tol > 1: rel_tol = ctx.num(1)
    if not (ctx.isfinite(a) and ctx.isfinite(b)): return Res(status='non_finite', reason='end', evals=0)
    if a == b: return Res(status='ok', value=Z, abs_error=Z, evals=0)
    sign = 1
    if b < a: a, b, sign = b, a, -1
    xs, wk, wg, dw, M = rule.x, rule.wk, rule.wg, rule.d, rule.m
    xmax = rule.xmax

    def splittable(lo, hi):
        mid = Interval(ctx, lo, hi).c
        if not (lo < mid < hi): return False
        for p, q in ((lo, mid), (mid, hi)):
            ch = Interval(ctx, p, q)
            if not (p < ch.raw(-xmax) and ch.raw(xmax) < q): return False
        return True

    st = dict(splits=0, descents=0, maxdepth=0, undone=0, dropped=0, maxratio=0.0)

    def evaluate(lo, hi):
        iv = Interval(ctx, lo, hi)
        K = Sum(Z); D = Sum(Z); G = Sum(Z); L1 = Z
        for i in range(M):
            y = fn(iv.point(xs[i]))
            t = iv.term(wk[i], y)
            K.add(t); L1 += abs(t)
            if diffw: D.add(iv.term(dw[i], y))
            elif wg[i] != 0: G.add(iv.term(wg[i], y))
        if not ctx.isfinite(L1): return None
        Kv = K.get()
        Dv = abs(D.get()) if diffw else abs(Kv - G.get())
        if not (ctx.isfinite(Kv) and ctx.isfinite(Dv)): return 'overflow'
        if L1 > 0 and eps * L1 > 0: st['maxratio'] = max(st['maxratio'], float(Dv / (eps * L1)))
        return Kv, Dv, floor_c * max(eps * L1, dmin)

    work = [0]  # bisections + descent splits, against the cap

    def resolve(lo, hi, r):
        """r = evaluate(lo, hi). Returns ('overflow',), or (pieces, leaves): the finite pieces and the non-finite
        leaves of a descent of at most `descent` levels (iterative)."""
        if r == 'overflow': return ('overflow',)
        if r is not None: return ([(lo, hi) + r], [])
        st['descents'] += 1
        pieces, leaves = [], []
        stack = [(lo, hi, 0, None)]
        while stack:
            p, q, d, rr = stack.pop()
            if d > 0:
                rr = evaluate(p, q)
                if rr == 'overflow': return ('overflow',)
                if rr is not None: pieces.append((p, q) + rr); continue
            st['maxdepth'] = max(st['maxdepth'], d)
            if (descent is None or d < descent) and splittable(p, q) and work[0] < cap:
                work[0] += 1; st['splits'] += 1
                mid = Interval(ctx, p, q).c
                stack.append((mid, q, d + 1, None)); stack.append((p, mid, d + 1, None))
            else:
                leaves.append((p, q))
        return (pieces, leaves)

    heap, final = [], []
    seq = [0]
    S = dict(K=Sum(Z), E=Sum(Z), F=Sum(Z), X=Sum(Z))

    def add(K, e, fl, sgn=1):
        S['K'].add(sgn * K); S['E'].add(sgn * e); S['F'].add(sgn * fl); S['X'].add(sgn * (e - fl))

    def commit(pieces):
        for lo, hi, K, err, fl in pieces:
            e = max(err, fl)
            add(K, e, fl)
            if err <= fl or not splittable(lo, hi): final.append((lo, hi, K, e, fl))
            else:
                seq[0] += 1; heapq.heappush(heap, (-err, seq[0], lo, hi, K, e, fl))

    tainted = False
    G0 = Interval(ctx, a, b)
    pts = [a] + [G0.point(ctx.num(-1) + ctx.num(2 * i) / mesh) for i in range(1, mesh)] + [b]
    for i in range(mesh):
        if pts[i] < pts[i + 1]:
            res = resolve(pts[i], pts[i + 1], evaluate(pts[i], pts[i + 1]))
            if res[0] == 'overflow': return Res(status='non_finite', reason='overflow', evals=nev[0], **st)
            pieces, leaves = res
            if leaves: return Res(status='non_finite', reason='mesh panel', evals=nev[0], **st)
            commit(pieces)

    bis = 0; m = None; count = 0; heap_max = len(heap)
    status = reason = None
    while True:
        I, E, F = S['K'].get(), S['E'].get(), S['F'].get()
        if not (ctx.isfinite(I) and ctx.isfinite(E) and ctx.isfinite(F)): status, reason = 'non_finite', 'sum'; break
        if E <= max(abs_tol, rel_tol * abs(I), floor_mult * F):
            status, reason = ('not_converged', 'non-finite piece') if tainted else ('ok', None); break
        if not heap: status, reason = 'not_converged', 'final panels'; break
        if work[0] >= cap: status, reason = 'not_converged', 'cap'; break
        if stag_k is not None and count >= stag_k: status, reason = 'not_converged', 'stagnation'; break
        _, _, lo, hi, K, e, fl = heapq.heappop(heap)
        mid = Interval(ctx, lo, hi).c
        work[0] += 1; bis += 1
        pieces, leaves, ovf = [], [], False
        for p, q in ((lo, mid), (mid, hi)):
            res = resolve(p, q, evaluate(p, q))
            if res[0] == 'overflow': ovf = True; break
            pieces += res[0]; leaves += res[1]
        if ovf: status, reason = 'non_finite', 'overflow'; break
        if leaves and nf == 'undo':
            st['undone'] += 1; tainted = taint
            final.append((lo, hi, K, e, fl))           # stays in the sums as it was
        else:
            add(K, e, fl, -1)
            commit(pieces)
            for p, q in leaves:                        # 'carry'
                st['dropped'] += 1; tainted = taint
                add(Z, e, Z); final.append((p, q, Z, e, Z))
        heap_max = max(heap_max, len(heap))
        X = S['X'].get()
        if m is None or X < m: m = X; count = 0
        elif stag == 'i' and X > 10 * m: m = X; count = 0
        elif stag == 'ii' and X > 10 * m: m = X; count += 1
        else: count += 1

    if status == 'non_finite': return Res(status=status, reason=reason, evals=nev[0], bis=bis, **st)
    V, Er = Sum(Z), Sum(Z)
    for t in heap: V.add(t[4]); Er.add(t[5])
    for t in final: V.add(t[2]); Er.add(t[3])
    v, er = V.get(), Er.get()
    if not (ctx.isfinite(v) and ctx.isfinite(er)): return Res(status='non_finite', reason='sum', evals=nev[0], bis=bis, **st)
    return Res(status=status, reason=reason, value=sign * v, abs_error=er, evals=nev[0], bis=bis, work=work[0],
               heap_max=heap_max, tainted=tainted, **st)
