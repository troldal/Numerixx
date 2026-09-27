// Not a negative test: a user helper that takes a tolerance as a plain double.
// P2564 (consteval propagation) makes the TEMPLATE version an immediate function on GCC/Clang; MSVC lacks P2564.
#include "prelude.hpp"
#include <cstdio>
template<class T> constexpr auto make_solver(T t) { return r::bisection{nxx::x_tol{t}}; }   // escalates to consteval
int main() {
    constexpr auto s = make_solver(1e-6);                      // constant argument: OK only with P2564
    auto res = s(f, nxx::bracket{0.0, 2.0});
    std::printf("escalation probe: %s\n", res ? "ok" : "fail");
}
