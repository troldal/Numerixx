// E17 (Chapter 3, section 3.6): a function less accurate than machine epsilon. exp(x) with its value perturbed by a
// deterministic relative noise eta = 1e-8 (book::noisy from code/noise.hpp, the same perturbation as DS19), as the
// output of a solver converged to 1e-8 would be. The noise level replaces epsilon in the rounding error, so the
// step that balances it is about eta^(1/3) s and the accuracy about eta^(2/3). Ridders' estimate assumes epsilon:
// at x = 1 it happens to cover the error, which is why the example prints rather than asserts it, and the text
// quotes the coverage over 300 points from DS19 instead.

#include "noise.hpp"

#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

int main()
{
    using nxx::deriv::derivative;
    using nxx::deriv::Ridders;

    // clang-format off
// <noisy-function>
// exp(x), accurate to a relative eta = 1e-8 only
double const eta = 1e-8;
auto const f = [eta](double x) {
    return book::noisy(std::exp(x), x, eta);
};
double const x = 1.0; // so s = max(|x|, 1) = 1

auto const d = derivative(f, x);          // eps^(1/3) s
auto const d_eta = derivative(f, x, std::cbrt(eta));
auto const r = derivative<Ridders>(f, x);
// </noisy-function>
    // clang-format on

    if (!d || !d_eta || !r) return 1;
    double const exact = std::exp(x);
    double const rel_d = std::abs(d->value - exact) / exact;
    double const rel_eta = std::abs(d_eta->value - exact) / exact;
    double const r_actual = std::abs(r->value - exact);
    std::println("default step {:.1e}: relative error {:.1e}", nxx::math::root_epsilon<3, double>(), rel_d);
    std::println("step eta^(1/3) = {:.1e}: relative error {:.1e}", std::cbrt(eta), rel_eta);
    std::println("Ridders: actual error {:.1e}, estimate {:.1e}, ratio {:.2f}", r_actual, *r->abs_error,
        r_actual / *r->abs_error);

    // clang-format off
// <noisy-coverage>
// does Ridders' estimate cover the actual error?
int covered = 0, cases = 0;
for (int k = 0; k <= 20; ++k) {
    double const xk = 0.5 + 0.1 * k;
    auto const rk = derivative<Ridders>(f, xk);
    if (!rk) continue;
    ++cases;
    double const err = std::abs(rk->value - std::exp(xk));
    if (err <= *rk->abs_error) ++covered;
}
std::println("covered in {} of {} points", covered, cases);
// </noisy-coverage>
    // clang-format on

    // The fitted step beats the default (DS19: 1.2e-6 against 9.4e-4); the model gives about eta^(2/3) = 5e-6 for
    // the fitted step, asserted with a factor 20 of room for the luck of the noise at two points.
    bool const ok = rel_eta < rel_d && rel_eta <= 1e-4 && cases > 0;
    return ok ? 0 : 1;
}
