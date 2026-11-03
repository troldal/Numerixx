# G(n)-K(2n+1) rules, computed independently of any published table.
#  - Gauss: Golub-Welsch (mpmath eigsy on the Jacobi matrix).
#  - Kronrod: the Stieltjes polynomial E_{n+1} = P_{n+1} + sum c_j P_j (j = n-1, n-3, ...), its coefficients exact
#    (Fractions) from the orthogonality conditions int E P_n P_k = 0, k odd <= n, a triangular system; its zeros by
#    mpmath polyroots at high precision; weights by the interpolatory formulas
#       Kronrod node xi:  w = 2 / ((n+1) P_n(xi) E'(xi))
#       Gauss node x_i:   w = w_G + 2 / ((n+1) P_n'(x_i) E(x_i))      (difference weight d_i = the second term)
#    derived from the Lagrange basis (int P_n x^n = 2 / ((2n+1) k_n), k_{n+1} / k_n = (2n+1)/(n+1)); cross-checked
#    against the moment equations int x^m = sum w x^m, m <= 3n+1.
#  - "plain double": the same formulas evaluated in double (Newton on P_n and on E in double, the recurrences in
#    double), emulating the arithmetic-only generator of BRIEF 3.4 in T = double.
# Usage: python -I rules.py   (checks against Boost's 115-digit 15-point tables and prints the double errors)
from fractions import Fraction as Fr
import math
import mpmath as mp

EPS = 2.0 ** -52


def legendre_polys(m):
    P = [[Fr(1)], [Fr(0), Fr(1)]]
    for k in range(1, m):
        a = [Fr(0)] + [Fr(2 * k + 1, k + 1) * c for c in P[k]]
        b = [Fr(k, k + 1) * c for c in P[k - 1]] + [Fr(0)] * 2
        P.append([a[i] - b[i] for i in range(k + 2)])
    return P


def pmul(p, q):
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                r[i + j] += a * b
    return r


def pint(p):  # integral over [-1, 1]
    return sum(c * Fr(2, i + 1) for i, c in enumerate(p) if i % 2 == 0)


def stieltjes(n):
    """Exact Legendre coefficients c[j] of E_{n+1} with c[n+1] = 1."""
    P = legendre_polys(n + 2)
    c = {n + 1: Fr(1)}
    # condition k (odd, k <= n) involves c_j for j >= n - k: solve for c_{n-k}
    for k in range(1, n + 1, 2):
        j0 = n - k
        PnPk = pmul(P[n], P[k])
        rhs = sum(cj * pint(pmul(P[j], PnPk)) for j, cj in c.items())
        coef = pint(pmul(P[j0], PnPk))
        c[j0] = -rhs / coef
    return c, P


def ref_rule(n, dps=120):
    """High-precision (mpf) nodes/weights: returns dict with lists sorted ascending over [-1, 1]."""
    with mp.workdps(dps):
        c, P = stieltjes(n)
        # Gauss by Golub-Welsch
        J = mp.zeros(n, n)
        for k in range(1, n):
            b = k / mp.sqrt(4 * mp.mpf(k) ** 2 - 1)
            J[k - 1, k] = b; J[k, k - 1] = b
        Ev, Q = mp.eigsy(J)
        g = sorted((Ev[i] if abs(Ev[i]) > mp.mpf(10) ** (-dps + 10) else mp.mpf(0), 2 * Q[0, i] ** 2) for i in range(n))
        # E in the monomial basis
        Em = [Fr(0)] * (n + 2)
        for j, cj in c.items():
            for i, a in enumerate(P[j]):
                Em[i] += cj * a
        coeffs = [mp.mpf(Em[i].numerator) / Em[i].denominator for i in range(n + 1, -1, -1)]
        xi = sorted(mp.re(r) for r in mp.polyroots(coeffs, maxsteps=500, extraprec=4 * dps))

        def Pn_and_d(x, m):
            p0, p1 = mp.mpf(1), x
            d0, d1 = mp.mpf(0), mp.mpf(1)
            if m == 0: return p0, d0
            for k in range(1, m):
                p2 = ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
                d2 = d0 + (2 * k + 1) * p1  # P'_{k+1} = P'_{k-1} + (2k+1) P_k
                p0, p1, d0, d1 = p1, p2, d1, d2
            return p1, d1

        def E_and_d(x):
            s = mp.mpf(0); sd = mp.mpf(0)
            for j, cj in c.items():
                p, d = Pn_and_d(x, j)
                q = mp.mpf(cj.numerator) / cj.denominator
                s += q * p; sd += q * d
            return s, sd

        # polish Kronrod zeros by Newton
        xi = [x - E_and_d(x)[0] / E_and_d(x)[1] for x in xi]
        xi = [x - E_and_d(x)[0] / E_and_d(x)[1] for x in xi]
        nodes = []
        for x, w in g:
            p, dp = Pn_and_d(x, n)
            e, _ = E_and_d(x)
            d = 2 / ((n + 1) * dp * e)
            nodes.append((x, w + d, w, d))
        for x in xi:
            p, _ = Pn_and_d(x, n)
            _, de = E_and_d(x)
            wk = 2 / ((n + 1) * p * de)
            nodes.append((x, wk, mp.mpf(0), wk))
        nodes.sort()
        # moment check
        worst = mp.mpf(0)
        deg = 3 * n + 1 if n % 2 == 0 else 3 * n + 2
        for m in range(deg + 1):
            ex = mp.mpf(2) / (m + 1) if m % 2 == 0 else mp.mpf(0)
            worst = max(worst, abs(sum(w * x ** m for x, w, _, _ in nodes) - ex))
        worst_g = max(abs(sum(w * x ** m for x, w in g) - (mp.mpf(2) / (m + 1) if m % 2 == 0 else 0)) for m in range(2 * n))
        return dict(n=n, x=[t[0] for t in nodes], wk=[t[1] for t in nodes], wg=[t[2] for t in nodes],
                    d=[t[3] for t in nodes], c=c, moment_err=worst, moment_err_g=worst_g)


def double_rule(n, ref):
    """The same formulas evaluated in IEEE double, Newton from the rounded reference node (as a safeguarded Newton
    in T converges to the same fixed point)."""
    c = {j: float(v) for j, v in ref['c'].items()}  # the coefficients rounded once (a triangular solve in T adds a few eps)

    def Pn_and_d(x, m):
        p0, p1 = 1.0, x
        d0, d1 = 0.0, 1.0
        if m == 0: return p0, d0
        for k in range(1, m):
            p2 = ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
            d2 = d0 + (2 * k + 1) * p1
            p0, p1, d0, d1 = p1, p2, d1, d2
        return p1, d1

    def E_and_d(x):
        s = 0.0; sd = 0.0
        for j in sorted(c):
            p, d = Pn_and_d(x, j)
            s += c[j] * p; sd += c[j] * d
        return s, sd

    def newton(fd, x):
        for _ in range(60):
            v, d = fd(x)
            if d == 0: break
            nx = x - v / d
            if nx == x: break
            if abs(nx - x) <= 2 * EPS * abs(x):
                x = nx; break
            x = nx
        return x

    out = []
    for i, X in enumerate(ref['x']):
        is_gauss = ref['wg'][i] != 0
        x0 = float(X)
        if x0 == 0.0:
            x = 0.0
        elif is_gauss:
            x = newton(lambda t: Pn_and_d(t, n), x0)
        else:
            x = newton(E_and_d, x0)
        if is_gauss:
            _, dp = Pn_and_d(x, n)
            wg = 2.0 / ((1.0 - x) * (1.0 + x) * dp * dp)
            e, _ = E_and_d(x)
            d = 2.0 / ((n + 1) * dp * e)
            wk = wg + d
        else:
            p, _ = Pn_and_d(x, n)
            _, de = E_and_d(x)
            wk = 2.0 / ((n + 1) * p * de)
            wg = 0.0; d = wk
        out.append((x, wk, wg, d))
    # symmetrise as a generator storing the non-negative half would
    m = len(out)
    for i in range(m // 2):
        x, wk, wg, d = out[m - 1 - i]
        out[i] = (-x, wk, wg, d)
    return out


def rounded_rule(ref):
    return [(float(ref['x'][i]), float(ref['wk'][i]), float(ref['wg'][i]), float(ref['d'][i])) for i in range(len(ref['x']))]


def force_centre(rule):
    """d_i with their sum forced to 0 through the centre weight (BRIEF 3.6)."""
    m = len(rule); ci = m // 2
    assert rule[ci][0] == 0.0
    s = 0.0
    for i in range(m):
        if i != ci: s += rule[i][3]
    out = list(rule)
    x, wk, wg, d = out[ci]
    out[ci] = (x, wk, wg, -s)
    return out


BOOST_K15_X = ['0', '2.0778495500789846760068940377324491347978440714517064971384573461986693844943520226910343227183698530560857645062738e-01',
               '4.0584515137739716690660641207696146334738201409937012638704325179466381322612565532831268972774658776528675866604802e-01',
               '5.8608723546769113029414483825872959843678075060436095130499289319880373607444407464511674498935942098956811555121368e-01',
               '7.4153118559939443986386477328078840707414764714139026011995535196742987467218051379282683236686324705969251809311201e-01',
               '8.6486442335976907278971278864092620121097230707408814860145771276706770813259572103585847859604590541475281326027862e-01',
               '9.4910791234275852452618968404785126240077093767061778354876910391306333035484014080573077002792572414430073966699522e-01',
               '9.9145537112081263920685469752632851664204433837033470129108741357244173934653407235924503509626841760744349505339308e-01']
BOOST_K15_W = ['2.0948214108472782801299917489171426369776208022370431671299800656137515132325648616816908211675949102392971459688215e-01',
               '2.0443294007529889241416199923464908471651760418071835742447095312045467698546598879348374292009347554167803659293064e-01',
               '1.9035057806478540991325640242101368282607807545535835588544088036744058072410212679605964605106377593834568683551139e-01',
               '1.6900472663926790282658342659855028410624490030294424149734006755695680921619029112936702403855359908156070095656537e-01',
               '1.4065325971552591874518959051023792039988975724799857556174546893312708093090950408097379122415555910759700350860143e-01',
               '1.0479001032225018383987632254151801744375665421383061189339065133963746321576289524167571627509311333949422518201492e-01',
               '6.3092092629978553290700663189204286665071157211550707113605545146983997477964874928199170264504441995865872491871943e-02',
               '2.2935322010529224963732008058969591993560811275746992267507430254711815787976075946156368168156289483493617134063245e-02']
BOOST_G7_W = ['4.1795918367346938775510204081632653061224489795918367346938775510204081632653061224489795918367346938775510204081633e-01',
              '3.8183005050511894495036977548897513387836508353386273475108345103070554643412970834868465934404480145031467176458536e-01',
              '2.7970539148927666790146777142377958248692506522659876453701403269361881043056267681324094290119761876632337521337205e-01',
              '1.2948496616886969327061143267908201832858740225994666397720863872465523497204230871562541816292084508948440200163443e-01']


def check_boost():
    ref = ref_rule(7)
    with mp.workdps(120):
        half = list(zip(ref['x'], ref['wk'], ref['wg']))[7:]
        ex = max(abs(half[i][0] - mp.mpf(BOOST_K15_X[i])) for i in range(8))
        ew = max(abs(half[i][1] - mp.mpf(BOOST_K15_W[i])) for i in range(8))
        gw = [half[i][2] for i in range(0, 8, 2)]
        eg = max(abs(gw[i] - mp.mpf(BOOST_G7_W[i])) for i in range(4))
        print(f"K15 vs Boost 115-digit: nodes {mp.nstr(ex, 3)}, K weights {mp.nstr(ew, 3)}, G weights {mp.nstr(eg, 3)}")
        print(f"moment residual K15 (deg 23) {mp.nstr(ref['moment_err'], 3)}, G7 (deg 13) {mp.nstr(ref['moment_err_g'], 3)}")
    return ref


def errors_in_eps(rule, ref):
    with mp.workdps(60):
        dx = [abs(mp.mpf(rule[i][0]) - ref['x'][i]) / EPS for i in range(len(rule))]
        dk = [abs(mp.mpf(rule[i][1]) - ref['wk'][i]) / EPS for i in range(len(rule))]
        dg = [abs(mp.mpf(rule[i][2]) - ref['wg'][i]) / EPS for i in range(len(rule))]
        dd = [abs(mp.mpf(rule[i][3]) - ref['d'][i]) / EPS for i in range(len(rule))]
        sumd = abs(sum(mp.mpf(r[3]) for r in rule)) / EPS
        sumk = abs(sum(mp.mpf(r[1]) for r in rule) - 2) / EPS
        sumg = abs(sum(mp.mpf(r[2]) for r in rule) - 2) / EPS
        return dict(max_dx=float(max(dx)), sum_dwk=float(sum(dk)), sum_dwg=float(sum(dg)), sum_dd=float(sum(dd)),
                    sum_d=float(sumd), sum_wk=float(sumk), sum_wg=float(sumg))


if __name__ == '__main__':
    ref = check_boost()
    for n in (7, 10, 15):
        r = ref_rule(n)
        print(f"n={n}: K{2*n+1} moment residual {mp.nstr(r['moment_err'], 3)}, G{n} {mp.nstr(r['moment_err_g'], 3)}")
    dbl = double_rule(7, ref)
    rnd = rounded_rule(ref)
    for name, rule in (('plain double', dbl), ('plain double, centre forced', force_centre(dbl)),
                       ('correctly rounded', rnd), ('correctly rounded, centre forced', force_centre(rnd))):
        print(f"{name:34s}", {k: round(v, 3) for k, v in errors_in_eps(rule, ref).items()})
    print("sum |d_i| =", sum(abs(r[3]) for r in rnd), " (K - G = sum d_i f_i)")
