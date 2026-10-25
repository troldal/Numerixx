// Numerixx demo: roots of a function of one variable in a bracket with nxx::roots.

#include <numerixx/roots.hpp>

#include <cmath>
#include <expected>
#include <format>
#include <numbers>
#include <print>
#include <string_view>

namespace {

    namespace rt = nxx::roots;

    template<typename T>
    void show(std::string_view label, std::expected<rt::Root<T>, rt::error> const &r)
    {
        if (r)
            std::println("   {:<30} {:.17g} +- {:.1e}", label, r->value, r->abs_error);
        else
            std::println("   {:<30} error: {}", label, r.error());
    }

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

    template<typename M, typename F>
    void count(std::string_view label, F const &fn, rt::Bracket<double> const &bracket, rt::Criteria<double> c = {})
    {
        int n = 0;
        auto const r = rt::root<M>(Counted<F>{ fn, &n }, bracket, c);
        if (r)
            std::println("   {:<30} {:.17g} +- {:.1e}  ({} evaluations)", label, r->value, r->abs_error, n);
        else
            std::println("   {:<30} error: {}  ({} evaluations)", label, r.error(), n);
    }

    void first_root()
    {
        std::println("1. sqrt(2) as the root of x^2 - 2 in [1, 2], by Brent's method");
        auto const f = [](double x) { return x * x - 2.0; };
        auto const r = rt::root<rt::Brent>(f, rt::Bracket{ 1.0, 2.0 });
        show("root", r);
        std::println("   {:<30} {:.17g}", "exact", std::numbers::sqrt2);
        // abs_error bounds the distance to the point where fn, as evaluated, changes sign.
        if (r) std::println("   {:<30} {:.1e}", "actual error", std::abs(r->value - std::numbers::sqrt2));
    }

    void methods()
    {
        // Bisection gains one bit per evaluation, whatever fn is. Brent interpolates and is much faster on a smooth fn
        // with a simple root, but slower on a multiple root, where the interpolation converges only linearly.
        std::println("\n2. Bisection against Brent");
        auto const dottie = [](double x) { return std::cos(x) - x; };
        auto const wallis = [](double x) { return (x * x - 2.0) * x - 5.0; };
        auto const triple = [](double x) { return (x - 1.0) * (x - 1.0) * (x - 1.0); };
        count<rt::Bisection>("cos x = x: Bisection", dottie, { 0.0, 1.0 });
        count<rt::Brent>("cos x = x: Brent", dottie, { 0.0, 1.0 });
        count<rt::Bisection>("x^3 - 2x - 5: Bisection", wallis, { 2.0, 3.0 });
        count<rt::Brent>("x^3 - 2x - 5: Brent", wallis, { 2.0, 3.0 });
        count<rt::Bisection>("(x - 1)^3: Bisection", triple, { 0.0, 3.0 });
        count<rt::Brent>("(x - 1)^3: Brent", triple, { 0.0, 3.0 });
    }

    void kepler()
    {
        // Kepler's equation E - e sin E = M gives a planet's eccentric anomaly E from its mean anomaly M. For e < 1,
        // f(0) = -M and f(pi) = pi - M bracket the root for every M in (0, pi).
        std::println("\n3. Kepler's equation E - e sin E = M for e = 0.9, a table over M");
        double const e = 0.9;
        for (double const m : { 0.1, 0.5, 1.0, 2.0, 3.0 }) {
            auto const fn = [e, m](double big_e) { return big_e - e * std::sin(big_e) - m; };
            show(std::format("M = {:.1f}", m), rt::root<rt::Brent>(fn, rt::Bracket{ 0.0, std::numbers::pi }));
        }
    }

    void tolerances()
    {
        // The default criteria give the full precision of the type at every magnitude: a root at 1e-12 keeps its 16
        // digits, with no tolerance to set. A looser tolerance saves evaluations.
        std::println("\n4. Tolerances");
        auto const trace = [](double x) { return std::log(x / 1e-12); };
        count<rt::Brent>("ln(x / 1e-12): default", trace, { 1e-20, 1.0 });
        auto const f = [](double x) { return x * x - 2.0; };
        count<rt::Brent>("x^2 - 2: default", f, { 1.0, 2.0 });
        count<rt::Brent>("x^2 - 2: rel_tol = 1e-6", f, { 1.0, 2.0 }, { .rel_tol = 1e-6 });

        // A sign change near 0 that fn misses is approached down to the type's smallest numbers: set abs_tol there.
        auto const tiny = [](double x) { return x - 1e-300; };
        count<rt::Bisection>("x - 1e-300: default", tiny, { 0.0, 1.0 });
        count<rt::Bisection>("x - 1e-300: abs_tol = 1e-12", tiny, { 0.0, 1.0 }, { .abs_tol = 1e-12 });
    }

    void around_zero()
    {
        // A bracket around 0 is first split at 0, so a root at 0 costs one evaluation more.
        std::println("\n5. A bracket around 0");
        auto const sin = [](double x) { return std::sin(x); };
        count<rt::Brent>("sin x in [-1, 2]", sin, { -1.0, 2.0 });
        // Where fn is not finite at 0 (a removable singularity), that split is skipped.
        auto const sinc = [](double x) { return std::sin(x) / x - 0.5; };
        count<rt::Brent>("sin(x)/x = 0.5 in [-1, 3]", sinc, { -1.0, 3.0 });
    }

    template<typename T>
    void precision(std::string_view type)
    {
        auto const f = [](T x) { return x * x - T(2); };
        if (auto const r = rt::root<rt::Brent>(f, rt::Bracket{ T(1), T(2) }))
            std::println("   {:<12} abs_error {:.1e}, error {:.1e}",
                type,
                static_cast<double>(r->abs_error),
                static_cast<double>(std::abs(r->value - std::sqrt(T(2)))));
    }

    void types()
    {
        // The error is a few units of the type's precision. MSVC makes long double the same as double.
        std::println("\n6. sqrt(2) in float, double and long double");
        precision<float>("float");
        precision<double>("double");
        precision<long double>("long double");
    }

    void compile_time()
    {
        std::println("\n7. At compile time");
        constexpr auto f = [](double x) { return x * x - 2.0; };
        constexpr auto r = rt::root<rt::Brent>(f, rt::Bracket{ 1.0, 2.0 });
        static_assert(r.has_value() && nxx::math::abs(r->value - 1.4142135623730951) <= r->abs_error);
        show("constexpr sqrt(2)", r);
        // Bisection reaches a dyadic root exactly: 0.5, 0.25, then 0.375.
        static_assert(
            rt::root<rt::Bisection>([](double x) { return x - 0.375; }, rt::Bracket{ 0.0, 1.0 })->value == 0.375);
    }

    void errors()
    {
        std::println("\n8. Failures are values, not exceptions");
        auto const f = [](double x) { return x * x - 2.0; };
        show("x^2 - 2 in [2, 3]", rt::root<rt::Brent>(f, rt::Bracket{ 2.0, 3.0 })); // the same sign at both ends
        auto const log = [](double x) { return std::log(x) + 1.0; };
        show("log x + 1 in [0, 1]", rt::root<rt::Brent>(log, rt::Bracket{ 0.0, 1.0 })); // log 0 is -inf
        show("x^2 - 2, rel_tol = -1", rt::root<rt::Brent>(f, rt::Bracket{ 1.0, 2.0 }, { .rel_tol = -1.0 }));

        // A pole changes sign like a root, and is returned as one: check |fn| at the result where fn may have poles.
        auto const pole = [](double x) { return 1.0 / (x - 1.0 / 3.0); };
        auto const p = rt::root<rt::Brent>(pole, rt::Bracket{ 0.0, 1.0 });
        show("1 / (x - 1/3) in [0, 1]", p);
        if (p) std::println("   {:<30} {:.1e}", "|fn| there", std::abs(pole(p->value)));
    }

    void composition()
    {
        // Results compose through std::expected: the root of one solve sets up the next, and an error at any stage
        // stops the chain.
        std::println("\n9. Composition: the x where x e^x = 1, then the y where y^2 = that x");
        auto const omega = [](double x) { return x * std::exp(x) - 1.0; };
        auto const both = rt::root<rt::Brent>(omega, rt::Bracket{ 0.0, 1.0 }).and_then([](rt::Root<double> const &w) {
            auto const g = [x = w.value](double y) { return y * y - x; };
            return rt::root<rt::Brent>(g, rt::Bracket{ 0.0, 1.0 });
        });
        show("sqrt(W(1))", both);
        std::println("   {:<30} {:.17g}", "via std::sqrt", std::sqrt(0.56714329040978387300));
    }

} // namespace

int main()
{
    first_root();
    methods();
    kepler();
    tolerances();
    around_zero();
    types();
    compile_time();
    errors();
    composition();
}
