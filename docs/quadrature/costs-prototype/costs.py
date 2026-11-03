# The counts of docs/sphinx/quadrature/costs.md from independent implementations: Phase 1's Romberg prototype
# (romberg.py, proto3.py of its session), Phase 2's Gauss-Kronrod prototype with its final configuration (gk.py, gk2.py,
# rules.py; the configuration of final2.py) and the tanh-sinh prototype (../tanh-sinh-prototype, ts.py and es.py with
# the C++ exp-sinh table). Python floats are IEEE doubles; the Gauss-Kronrod prototype uses its own double rule, so its
# values may differ from the C++ in the last bits while its counts agree.
# Usage: python -I costs.py (needs mpmath)
import math, os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', 'tanh-sinh-prototype'))
sys.path.insert(0, here)
import es, gk, gk2, romberg, ts

CFG = dict(mesh=10, descent=2, nf='undo', taint=True, stag='iii', stag_k=32, diffw=False)

if __name__ == '__main__':
    os.chdir(os.path.join(here, '..', 'tanh-sinh-prototype'))  # es.py reads es_double.txt from there
    rule = gk.Rule(gk.Double, kind='double', force=False)
    for name, f in [('4/(1+x^2)', lambda x: 4.0 / (1.0 + x * x)), ('sqrt x', math.sqrt),
                    ('1/sqrt x', lambda x: 1.0 / math.sqrt(x))]:
        r = ts.integrate(f, 0.0, 1.0)
        print(f"{name:10} TanhSinh  {r['status']:14} {r.get('value')!r:22} {r.get('abs_error')!r:24}",
              r['evaluations'])
        n = [0]
        def counted(x, f=f):
            n[0] += 1
            return f(x)
        r = gk2.integrate(counted, 0.0, 1.0, rule=rule, rel_tol=2.0 ** -26, **CFG)
        print(f"{name:10} Kronrod   {r.status:14} {r.get('value')!r:22} {r.get('abs_error')!r:24} {n[0]}")
        n[0] = 0
        try:
            print(f"{name:10} Romberg   {romberg.romberg(counted, 0.0, 1.0)!r} {n[0]}")
        except Exception as e:  # 1 / sqrt x: fn(0) is infinite, which the prototype raises on
            print(f"{name:10} Romberg   raised {e!r} after {n[0]}")
    r = es.integrate(lambda x: 1.0 / (1.0 + x * x), -math.inf, math.inf, cpp=True)
    print(f"{'line':10} TanhSinh  {r['status']:14} {r.get('value')!r:22} {r.get('abs_error')!r:24} {r['evaluations']}")
