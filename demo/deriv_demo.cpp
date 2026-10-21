// Numerixx demo: numerical derivatives of a scalar function with nxx::deriv.

#include <numerixx/deriv.hpp>

#include <cmath>
#include <expected>
#include <limits>
#include <print>
#include <string_view>

namespace {

    using nxx::deriv::Central;
    using nxx::deriv::derivative;
    using nxx::deriv::Forward;
    using nxx::deriv::Ridders;

    std::string_view to_string(nxx::deriv::error e)
    {
        switch (e) {
        case nxx::deriv::error::non_finite:
            return "non_finite";
        }
        return "unknown";
    }

    template<typename T>
    void show(std::string_view label, std::expected<T, nxx::deriv::error> const &d)
    {
        if (d)
            std::println("   {:<28} {:.17g}", label, *d);
        else
            std::println("   {:<28} error: {}", label, to_string(d.error()));
    }

    // A derivative estimate, with its error estimate when the method provides one (Ridders). std::format does not
    // consider the conversion to T, so the value is passed explicitly.
    template<typename T>
    void show(std::string_view label, std::expected<nxx::deriv::Estimate<T>, nxx::deriv::error> const &d)
    {
        if (!d)
            std::println("   {:<28} error: {}", label, to_string(d.error()));
        else if (d->abs_error)
            std::println("   {:<28} {:.17g} +- {:.1e}", label, d->value, *d->abs_error);
        else
            std::println("   {:<28} {:.17g}", label, d->value);
    }

    void default_step()
    {
        std::println("1. d/dx exp(x) at x = 1 with the default step");
        auto const exp = [](double x) { return std::exp(x); };
        auto const d = derivative(exp, 1.0);
        show("estimate", d);
        std::println("   {:<28} {:.17g}", "exact", std::exp(1.0));
        if (d) std::println("   {:<28} {:.1e}", "relative error", std::abs(*d - std::exp(1.0)) / std::exp(1.0));
    }

    void step_size()
    {
        // Too large a step leaves truncation error (h^2 / 6 * f'''), too small a step leaves rounding error
        // (eps * f / h); the default step eps^(1/3) * max(|x|, 1) sits near the minimum of their sum.
        std::println("\n2. The error of d/dx sin(x) at x = 1 against the step h");
        auto const sin = [](double x) { return std::sin(x); };
        for (int k = 1; k <= 13; ++k) {
            double const h = std::pow(10.0, -k);
            if (auto const d = derivative(sin, 1.0, h))
                std::println("   h = 1e-{:<2}  error {:.1e}", k, std::abs(*d - std::cos(1.0)));
        }
        if (auto const d = derivative(sin, 1.0))
            std::println("   default h = {:.1e}  error {:.1e}",
                nxx::math::root_epsilon<3, double>(),
                std::abs(*d - std::cos(1.0)));
    }

    template<typename T>
    void precision(std::string_view type)
    {
        auto const exp = [](T x) { return std::exp(x); };
        if (auto const d = derivative(exp, T(1)))
            std::println("   {:<12} relative error {:.1e}",
                type,
                static_cast<double>(std::abs(*d - std::exp(T(1))) / std::exp(T(1))));
    }

    void types()
    {
        // The error with the default step is about eps^(2/3) of the type. MSVC makes long double the same as double.
        std::println("\n3. The same derivative in float, double and long double");
        precision<float>("float");
        precision<double>("double");
        precision<long double>("long double");
    }

    void compile_time()
    {
        // The central difference is exact for x^3 up to h^2: 3 x^2 + h^2, exact in binary for h = 1/4.
        std::println("\n4. At compile time: d/dx x^3 at x = 2 with h = 0.25");
        constexpr auto cube = [](double x) { return x * x * x; };
        constexpr auto d = derivative(cube, 2.0, 0.25);
        static_assert(d.has_value() && *d == 12.0625);
        show("constexpr estimate", d);
    }

    void errors()
    {
        std::println("\n5. Failures are values, not exceptions");
        auto const log = [](double x) { return std::log(x); };
        show("d/dx log(x) at x = 1", derivative(log, 1.0));
        show("d/dx log(x) at x = 0", derivative(log, 0.0)); // fn(-h) is NaN
    }

    void composition()
    {
        // Newton's iteration for x^3 = 2, with the derivative estimated numerically; an error would stop the chain.
        std::println("\n6. Composition: Newton's method for the cube root of 2");
        auto const f = [](double x) { return x * x * x - 2.0; };
        auto const newton_step = [&f](double x) {
            return derivative(f, x).transform([&f, x](double df) { return x - f(x) / df; });
        };
        std::expected<double, nxx::deriv::error> x = 1.0;
        for (int i = 1; i <= 6; ++i) {
            x = x.and_then(newton_step);
            show(std::format("step {}", i), x);
        }
        std::println("   {:<28} {:.17g}", "exact", std::cbrt(2.0));
    }

    template<typename S>
    void stencil_error(std::string_view name)
    {
        auto const exp = [](double x) { return std::exp(x); };
        if (auto const d = derivative<S>(exp, 1.0))
            std::println("   {:<12} relative error {:.1e}", name, std::abs(*d - std::exp(1.0)) / std::exp(1.0));
    }

    void stencils()
    {
        // A central rule of order P costs P evaluations and reaches about eps^(P/(P+1)) with its default step.
        std::println("\n7. Higher-order stencils: d/dx exp(x) at x = 1");
        stencil_error<Central<2>>("Central<2>");
        stencil_error<Central<4>>("Central<4>");
        stencil_error<Central<6>>("Central<6>");
        stencil_error<Central<10>>("Central<10>");

        // One-sided rules never evaluate fn on the other side of x, so they work at the edge of fn's domain.
        std::println("\n8. One-sided stencils: fn defined only for x >= 0 (exp, NaN below 0), at x = 0");
        auto const right = [](double x) { return x < 0.0 ? std::numeric_limits<double>::quiet_NaN() : std::exp(x); };
        show("Central<2>", derivative<Central<2>>(right, 0.0));
        show("Forward<2>", derivative<Forward<2>>(right, 0.0));
        show("Forward<4>", derivative<Forward<4>>(right, 0.0));
    }

    void ridders()
    {
        // Ridders adapts the step and estimates its own error: compare the estimate with the actual error.
        std::println("\n9. Ridders' extrapolation, with an error estimate");
        auto const exp = [](double x) { return std::exp(x); };
        auto const d = derivative<Ridders>(exp, 1.0);
        show("d/dx exp(x) at x = 1", d);
        if (d) std::println("   {:<28} {:.1e}", "actual error", std::abs(d->value - std::exp(1.0)));

        // Near the edge of fn's domain the default initial step is too large; Ridders restarts below it.
        auto const log = [](double x) { return std::log(x); };
        auto const e = derivative<Ridders>(log, 0.01);
        show("d/dx log(x) at x = 0.01", e);
        if (e) std::println("   {:<28} {:.1e}", "actual error", std::abs(e->value - 100.0));
    }

    void ridders_compile_time()
    {
        std::println("\n10. Ridders at compile time, error estimate included");

        // The tableau is exact for a cubic once two levels agree: 12 exactly, with the rounding floor as the estimate.
        constexpr auto cube = [](double x) { return x * x * x; };
        constexpr auto d = derivative<Ridders>(cube, 2.0, 1.0 / 16);
        static_assert(d.has_value() && d->value == 12.0 && d->abs_error.has_value());
        show("d/dx x^3 at x = 2", d);

        // Runge's function with the default step. The compiler runs the whole method, restart rules included, and
        // checks the actual error against the estimate: f'(x) = -50 x / (1 + 25 x^2)^2, -25 / 7.25^2 at x = 1/2.
        constexpr auto runge = [](double x) { return 1.0 / (1.0 + 25.0 * x * x); };
        constexpr auto r = derivative<Ridders>(runge, 0.5);
        constexpr double exact = -25.0 / (7.25 * 7.25);
        static_assert(r.has_value() && nxx::math::abs(r->value - exact) <= *r->abs_error);
        show("d/dx 1/(1 + 25x^2) at 1/2", r);
        std::println("   {:<28} {:.1e}", "actual error", nxx::math::abs(r->value - exact));

        // Failures are values in constant expressions too: fn is NaN below 0, so every level around 0 fails.
        constexpr auto half_line = [](double x) { return x < 0.0 ? std::numeric_limits<double>::quiet_NaN() : x; };
        static_assert(derivative<Ridders>(half_line, 0.0).error() == nxx::deriv::error::non_finite);
        show("d/dx x on x >= 0, at 0", derivative<Ridders>(half_line, 0.0));
    }

} // namespace

int main()
{
    default_step();
    step_size();
    types();
    compile_time();
    errors();
    composition();
    stencils();
    ridders();
    ridders_compile_time();
}
