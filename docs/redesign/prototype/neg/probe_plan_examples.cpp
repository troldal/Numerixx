// Plan examples verbatim (6.10/6.11/7.2) that are not in test_core: brent{x_tol{1e-12}, 60}, runtime tolerance.
#include "prelude.hpp"
#include <cstdio>
int main(int argc, char**) {
    auto a = r::brent{nxx::x_tol{1e-12}, 60}(f, nxx::bracket{1.0, 2.0});
    auto tol = nxx::tolerance<double>::make(argc * 1e-9);
    auto b = r::bisection{nxx::x_tol{*tol}}(f, nxx::bracket{1.0, 2.0});
    std::printf("brent{x_tol{1e-12},60}: %s x=%.17g; bisection{x_tol{runtime 1e-9}}: x=%.17g after %u it (true root 1.4142135623730951)\n",
                a ? "ok" : "fail", a ? a->x : 0.0, b ? b->x : 0.0, b ? b->used.iterations : 0u);
}
