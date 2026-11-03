// Numerixx demo: integrals of a function of one variable with nxx::quadrature.

#include <numerixx/quadrature.hpp>

#include <cmath>
#include <expected>
#include <format>
#include <limits>
#include <print>
#include <string_view>

#ifdef NUMERIXX_DEMO_MULTIPRECISION
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif

namespace {

    namespace qd = nxx::quadrature;

    // fn with a count of its evaluations, which the methods do not report.
    template<typename F>
    struct Counted
    {
        F fn;
        int *count;

        double operator()(double x) const
        {
            ++*count;
            return fn(x);
        }
    };

    template<typename M = qd::GaussKronrod<15>, typename F>
    void show(std::string_view label, F const &fn, double a, double b, double exact, qd::Criteria<double> c = {})
    {
        int n = 0;
        auto const r = qd::integrate<M>(Counted<F>{ fn, &n }, a, b, c);
        if (r)
            std::println("   {:<34} {:.17g} +- {:.1e}  (error {:.1e}, {} evaluations)", label, r->value, r->abs_error,
                         std::abs(r->value - exact), n);
        else // not_converged here, which carries the best estimate
            std::println("   {:<34} {}: best {:.17g} +- {:.1e}  ({} evaluations)", label, r.error(),
                         r.error().best->value, r.error().best->abs_error, n);
    }

} // namespace

int main()
{
    // 1. A smooth integral with the default method: adaptive G7-K15 from 10 panels, 150 evaluations at least. The
    //    estimate is |K15 - G7| summed over the panels: it estimates the 7-point rule's error, so it overstates the
    //    returned 15-point value's.
    std::println("A smooth integral, the default method (GaussKronrod<15>):");
    show("e^x on [0, 1]", [](double x) { return std::exp(x); }, 0.0, 1.0, std::exp(1.0) - 1.0);
    show("1 / (1 + x^2) on [0, 1], 1e-14", [](double x) { return 1.0 / (1.0 + x * x); }, 0.0, 1.0,
         0.78539816339744831, { .rel_tol = 1e-14 });

    // 2. A strong endpoint singularity: x^-0.9 cos x on [0, 1]. Beyond about x^-0.65, |K - G| understates the end
    //    panel's error by a ratio bisection does not change. Substituting x = u^k with k (1 + alpha) >= 1, here
    //    k = 10, turns it into 10 cos(u^10), smooth. TanhSinh, whose points crowd double exponentially towards the
    //    ends, takes it directly.
    std::println("\nAn endpoint singularity stronger than x^-0.65:");
    double const exact = 9.7718427116392900; // sum_n (-1)^n / ((2n)! (2n + 0.1))
    auto const singular = [](double x) { return std::pow(x, -0.9) * std::cos(x); };
    show("x^-0.9 cos x, directly", singular, 0.0, 1.0, exact);
    show("10 cos(u^10), x = u^10", [](double u) { return 10.0 * std::cos(std::pow(u, 10.0)); }, 0.0, 1.0, exact);
    show<qd::TanhSinh>("x^-0.9 cos x by TanhSinh", singular, 0.0, 1.0, exact);

    // 3. A failure with its best estimate: e^x with noise at 1e-10 never resolves at rel_tol 1e-12, and the
    //    bisections stop lowering the estimate (stagnation). Failure::best holds where it got.
    std::println("\nA failure, with the best estimate reached:");
    show("e^x (1 + 1e-10 sin 1e7 x), 1e-12", [](double x) { return std::exp(x) * (1.0 + 1e-10 * std::sin(1e7 * x)); },
         0.0, 1.0, std::exp(1.0) - 1.0, { .rel_tol = 1e-12 });

    // 4. Infinite intervals: TanhSinh alone takes an infinite end, by exp-sinh on a half line and sinh-sinh on the
    //    line; every other method returns non_finite. The maps centre at 0 or at the finite end at scale 1: a peak
    //    far from there is split at the peak, or shifted into fn.
    std::println("\nInfinite intervals, by TanhSinh:");
    double const inf = std::numeric_limits<double>::infinity();
    show<qd::TanhSinh>("x^2 e^-x on [0, inf)", [](double x) { return x * x * std::exp(-x); }, 0.0, inf, 2.0);
    show<qd::TanhSinh>("e^-x^2 on (-inf, inf)", [](double x) { return std::exp(-x * x); }, -inf, inf,
                       1.7724538509055160);
    show<qd::TanhSinh>("x^-1.5 on [1, inf)", [](double x) { return 1.0 / (x * std::sqrt(x)); }, 1.0, inf, 2.0);
    std::println("   the default on [0, inf): {}", qd::integrate([](double x) { return std::exp(-x); }, 0.0, inf).error());

    // 5. A multiprecision integral: in a 50-digit type the default rel_tol is sqrt(eps), about 1e-25; a caller who
    //    chose the type for its digits sets rel_tol. The nodes are computed once, at the first use.
#ifdef NUMERIXX_DEMO_MULTIPRECISION
    using mp = boost::multiprecision::cpp_bin_float_50;
    std::println("\nMultiprecision (cpp_bin_float_50), 4 atan(1) = pi:");
    for (char const *tol : { "1e-25", "1e-45" }) {
        auto const r =
            qd::integrate([](mp const &x) { return mp(4 / (1 + x * x)); }, mp(0), mp(1), { .rel_tol = mp(tol) });
        if (r) std::println("   rel_tol {:<6} {} +- {}", tol, r->value.str(48), r->abs_error.str(2));
    }
#else
    std::println("\nMultiprecision: built without Boost.Multiprecision (it comes with the tests).");
#endif
}
