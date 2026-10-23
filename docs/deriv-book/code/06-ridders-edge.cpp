// E10 (Chapter 6): Ridders near the edge of a domain. log x at x = 0.01 has s = max(|x|, 1) = 1, so the default
// initial step 1/32 reaches below 0, where log is NaN: Ridders restarts below the failing step and still returns
// 100 with an estimate that covers the error, at the price of more evaluations (measured 24, against 14 in the
// interior). Central<2> with its default step eps^(1/3) (about 6e-6, no restart possible) does not fail, but its
// relative error is the truncation h^2 / (3 x^2), about 1e-7: the step is scaled with s = 1, not with x. Quoted
// with \snippet{06-ridders-edge}{edge-restart}.

#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

int main()
{
    using namespace nxx::deriv;
    // clang-format off
// <edge-restart>
// log x at 0.01: the default initial step 1/32 crosses 0,
// where log is NaN. Ridders restarts below the failing
// step; a stencil with its default step (about 6e-6) never
// reaches 0.
int calls = 0;
auto const f = [&calls](double x) {
    ++calls;
    return std::log(x);
};
auto const r = derivative<Ridders>(f, 0.01);
auto const d = derivative<Central<2>>(f, 0.01);
// </edge-restart>
    // clang-format on
    if (!r || !d) return 1;
    double const exact = 100.0; // 1 / x
    double const actual = std::abs(r->value - exact);
    double const stencil_error = std::abs(d->value - exact) / exact;
    std::println("Ridders     {:.12f} +- {:.1e}, error {:.1e}, {} evaluations in all",
        r->value, *r->abs_error, actual, calls - 2);
    std::println("Central<2>  {:.12f}  rel. error {:.1e}", d->value, stencil_error);
    // Ridders: covered, and more than the 16 evaluations an interior point of f' takes at most (the restart), within
    // the cap of 64. The stencil: the measured relative error is 1.2e-7 (truncation at a step scaled with s = 1).
    int const ridders_calls = calls - 2; // the stencil made the last two calls
    bool const ok = actual <= *r->abs_error && ridders_calls > 16 && ridders_calls <= 64 && stencil_error <= 1e-6;
    return ok ? 0 : 1;
}
