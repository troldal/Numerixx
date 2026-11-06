// Sketch 7: a method with a parameter, given as a value: the bracket cut at a configurable fraction from the end
// with the smaller |f| (0.5 is bisection). And a method whose parameter has no default (named form refused).
#include <numerixx/roots.hpp>
#include <cmath>
#include <print>

namespace rt = nxx::roots;

namespace user {

    template<nxx::real P = double>
    struct Section
    {
        P fraction = P(1) / P(3); // in (0, 1); a double parameter used with any T: T(fraction)

        template<nxx::real T>
        struct State { T a, fa, b, fb; rt::Step<T> step; };

        template<nxx::real T>
        static constexpr rt::Step<T> report(T const &a, T const &fa, T const &b, T const &fb)
        { return rt::Step<T>{ a, fa, b, fb, nxx::math::abs(fa) <= nxx::math::abs(fb) ? a : b }; }

        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
        { return State<T>{ s.lower, s.at_lower[0], s.upper, s.at_upper[0], report(s.lower, s.at_lower[0], s.upper, s.at_upper[0]) }; }

        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &) const
        {
            T const w = T(T(s.b / T(2)) - T(s.a / T(2))); // half the width, no overflow
            T const r = T(T(2) * T(fraction));
            bool const from_a = nxx::math::abs(s.fa) <= nxx::math::abs(s.fb);
            T const x = from_a ? T(s.a + T(r * w)) : T(s.b - T(r * w));
            if (!(s.a < x && x < s.b)) { s.step.exhausted = true; return {}; }
            T const fx = T(std::invoke(fn, x));
            if (!nxx::math::isfinite(fx)) return std::unexpected(rt::error::non_finite);
            if (fx == T(0)) { s.step = rt::Step<T>{ x, fx, x, fx, x }; return {}; }
            if ((fx < T(0)) == (s.fa < T(0))) { s.a = x; s.fa = fx; } else { s.b = x; s.fb = fx; }
            s.step = report(s.a, s.fa, s.b, s.fb);
            return {};
        }
    };

    // The same with a required parameter (no default constructor).
    struct Fixed : Section<>
    {
        explicit constexpr Fixed(double f) : Section<>{ f } {}
    };

} // namespace user

int main()
{
    auto const cubic = [](double x) { return x * x * x - 2.0 * x - 5.0; };
    double const exact = 2.0945514815423265;
    for (double fr : { 0.5, 1.0 / 3.0, 0.1, 0.01 }) {
        int calls = 0;
        auto counted = [&](double x) { ++calls; return cubic(x); };
        auto const r = rt::solve(user::Section<>{ .fraction = fr }, counted, rt::Bracket{ 2.0, 3.0 });
        if (r) std::println("fraction {:<6.3g} {:.17g} +- {:.3g}  |x - exact| {:.3g}  calls {}", fr, r->value, r->abs_error,
                            std::abs(r->value - exact), calls);
        else std::println("fraction {:<6.3g} error {} (calls {})", fr, r.error(), calls);
    }
    std::println("named, default fraction: {:.17g}", rt::solve<user::Section<>>(cubic, rt::Bracket{ 2.0, 3.0 })->value);
    std::println("value Fixed(0.25): {:.17g}", rt::solve(user::Fixed(0.25), cubic, rt::Bracket{ 2.0, 3.0 })->value);
    int n = 0;
    (void)rt::solve(user::Fixed(0.25), rt::observer([&](rt::Step<double> const &) { ++n; }), cubic, rt::Bracket{ 2.0, 3.0 });
    std::println("value + observer: {} reports", n);
    // A bad parameter: fraction 0 means no progress; the method reports exhausted -> not_converged.
    auto const bad = rt::solve(user::Section<>{ .fraction = 0.0 }, cubic, rt::Bracket{ 2.0, 3.0 });
    std::println("fraction 0: {}", bad ? "root" : std::format("error {}", bad.error()));
    // fraction 1.5: outside the bracket on the first step -> exhausted
    auto const bad2 = rt::solve(user::Section<>{ .fraction = 1.5 }, cubic, rt::Bracket{ 2.0, 3.0 });
    std::println("fraction 1.5: {}", bad2 ? "root" : std::format("error {}", bad2.error()));
#ifdef NAMED_FIXED
    (void)rt::solve<user::Fixed>(cubic, rt::Bracket{ 2.0, 3.0 });
#endif
}
