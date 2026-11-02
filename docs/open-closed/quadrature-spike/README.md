# Spike: the quadrature method protocol (`docs/quadrature/PHASE3A_API.md`, draft 3)

Evidence for the proposal, not part of the build. From the project root, with a scratch directory `$S`:

```sh
mkdir -p $S/include/numerixx && cp include/numerixx/math.hpp $S/include/numerixx/
# normalize: integrate() orders the ends (Q4's recommendation); omit it for the alternative.
python3 docs/open-closed/quadrature-spike/quadrature.patch.py include/numerixx/quadrature.hpp $S/draft2.hpp normalize
# Draft 3 on top: Step, the stepped Romberg and Gauss-Kronrod, the loop in integrate(), observed().
python3 docs/open-closed/quadrature-spike/steps.py $S/draft2.hpp $S/include/numerixx/quadrature.hpp

# Bit for bit: the same battery against today's header and the patched one (Boost.Multiprecision on the include path).
g++ -std=c++23 -O1 -DBOOST_MP_STANDALONE -Iinclude   docs/open-closed/quadrature-spike/dump.cpp -o $S/dump_old
g++ -std=c++23 -O1 -DBOOST_MP_STANDALONE -I$S/include docs/open-closed/quadrature-spike/dump.cpp -o $S/dump_new
$S/dump_old > $S/old.txt && $S/dump_new > $S/new.txt && cmp $S/old.txt $S/new.txt

# Methods written outside nxx, on the public protocol alone.
g++ -std=c++23 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -I$S/include docs/open-closed/quadrature-spike/user_methods.cpp

# Progress: observed(), the library's methods stepped, a user's iterative method.
g++ -std=c++23 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -I$S/include docs/open-closed/quadrature-spike/progress.cpp

# Printing each step: the example of PHASE3A_API.md, "Example: printing each step".
g++ -std=c++23 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -I$S/include docs/open-closed/quadrature-spike/print_steps.cpp

# Orientation: integrate(fn, b, a) against -integrate(fn, a, b), the whole result, with today's header.
g++ -std=c++23 -O1 -Iinclude docs/open-closed/quadrature-spike/orient.cpp
```

The existing tests `tests/test_quadrature_{gauss,romberg,kronrod}.cpp` are built against `$S/include` unchanged.
