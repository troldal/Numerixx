// Spike (draft 2): methods written outside nxx, against the public protocol only (no nxx::detail), used through
// integrate(). Built against the header that quadrature.patch.py writes with `normalize` (integrate() passes a < b).
#include <numerixx/quadrature.hpp>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>

namespace qd = nxx::quadrature;

namespace user {

    // A fixed rule with a run-time parameter: the composite midpoint rule. integrate() passes a < b, both finite.
    // fn is evaluated only inside [a, b]: a width that overflows is refused, and a point is clamped (it can round
    // past b by an ulp). The parameter is the method's precondition, asserted here.
    struct Midpoint
    {
        int panels = 100;

        template<nxx::real T, typename F>
        constexpr std::expected<T, qd::error> apply(F &fn, T const &a, T const &b) const
        {
            assert(panels >= 1);
            T const h = T(T(b - a) / T(panels));
            if (!nxx::math::isfinite(h)) return std::unexpected(qd::error::non_finite);
            T sum = T(0);
            for (int i = 0; i < panels; ++i) {
                T x = T(a + T(h * T(T(i) + T(1) / T(2))));
                if (b < x) x = b;
                sum = T(sum + T(std::invoke(fn, x)));
            }
            T const r = T(sum * h);
            if (!nxx::math::isfinite(r)) return std::unexpected(qd::error::non_finite);
            return r;
        }
    };
    static_assert(qd::fixed_rule<Midpoint, double>); // beside the definition: the protocol's diagnostics are here

    // A combinator over any estimating method: x = a + (b - a) u^2 at the lower end a (integrate() passes a < b).
    // Smooth after the substitution for fn written in x - a, or with a = 0; elsewhere x is quantised to ulp(a).
    template<typename Inner>
    struct SquareAtLowerEnd
    {
        Inner inner{};

        template<nxx::real T, typename F>
        constexpr std::expected<qd::Integral<T>, qd::Failure<T>>
            estimate(F &fn, T const &a, T const &b, qd::Criteria<T> const &criteria) const
        {
            T const width = T(b - a);
            if (!nxx::math::isfinite(width)) return std::unexpected(qd::Failure<T>{ qd::error::non_finite, std::nullopt });
            auto const g = [&](T const &u) -> T {
                T x = T(a + T(width * T(u * u)));
                if (b < x) x = b;
                return T(T(T(2) * width) * T(u * T(std::invoke(fn, x))));
            };
            return qd::integrate(inner, g, T(0), T(1), criteria);
        }
    };

    // What a method receives, recorded through a pointer (the idiom for state: the member is const).
    struct Seen
    {
        double a = 0, b = 0, abs_tol = -1, rel_tol = -1;
        int calls = 0;
    };
    struct Spy
    {
        Seen *seen = nullptr;

        template<nxx::real T, typename F>
        constexpr std::expected<qd::Integral<T>, qd::Failure<T>>
            estimate(F &, T const &a, T const &b, qd::Criteria<T> const &c) const
        {
            *seen = Seen{ a, b, c.abs_tol, c.rel_tol, seen->calls + 1 };
            return qd::Integral<T>{ T(1), T(0) };
        }
    };

    // Broken methods: integrate() turns their results into error::non_finite in every build.
    struct NaNValue
    {
        template<nxx::real T, typename F>
        constexpr std::expected<qd::Integral<T>, qd::Failure<T>> estimate(F &, T const &, T const &,
            qd::Criteria<T> const &) const
        { return qd::Integral<T>{ std::numeric_limits<T>::quiet_NaN(), T(0) }; }
    };
    struct NegativeError
    {
        template<nxx::real T, typename F>
        constexpr std::expected<qd::Integral<T>, qd::Failure<T>> estimate(F &, T const &, T const &,
            qd::Criteria<T> const &) const
        { return qd::Integral<T>{ T(1), T(-1) }; }
    };

    // Not methods: the wrong return type, a member that is not const; and a type of both kinds.
    struct WrongReturn
    {
        template<nxx::real T, typename F>
        constexpr T apply(F &, T const &, T const &) const { return T(0); }
    };
    struct NotConst
    {
        template<nxx::real T, typename F>
        constexpr std::expected<T, qd::error> apply(F &, T const &, T const &) { return T(0); }
    };
    struct Both
    {
        template<nxx::real T, typename F>
        constexpr std::expected<T, qd::error> apply(F &, T const &, T const &) const { return T(0); }
        template<nxx::real T, typename F>
        constexpr std::expected<qd::Integral<T>, qd::Failure<T>> estimate(F &, T const &, T const &,
            qd::Criteria<T> const &) const
        { return qd::Integral<T>{ T(0), T(0) }; }
    };

} // namespace user

constexpr double square(double x) { return x * x; }

// A user fixed rule, in constant expressions, named and by value; the midpoint rule errs by (b - a) h^2 f'' / 24,
// exactly 1/4 for x^2 on [0, 3] with 3 panels (h = 1, every node dyadic, so 8.75 is exact).
static_assert(qd::integrate(user::Midpoint{ .panels = 3 }, square, 0.0, 3.0).value() == 8.75);
static_assert(qd::integrate(user::Midpoint{ .panels = 3 }, square, 3.0, 0.0).value() == -8.75); // integrate() negates
static_assert(qd::integrate<user::Midpoint>(square, 0.0, 3.0).value() > 8.999);
static_assert(qd::integrate<user::Midpoint>(square, 0.0, std::numeric_limits<double>::infinity())
              == std::unexpected(qd::error::non_finite)); // the library's validation, before the method
// A combinator over a library method, in constant expressions: g(u) = 54 u^5, exact for K15.
static_assert(qd::integrate(user::SquareAtLowerEnd<qd::GaussKronrod<15>>{}, square, 0.0, 3.0)->value > 8.999999);
static_assert(qd::integrate(user::SquareAtLowerEnd<qd::Romberg>{}, square, 3.0, 0.0)->value < -8.999999);

// The concepts (F defaults to T(T const &)), and the entry points' constraints.
template<typename M>
concept integrable = requires { qd::integrate<M>(square, 0.0, 1.0); };
static_assert(qd::estimating_method<user::SquareAtLowerEnd<qd::Romberg>, double>);
static_assert(!qd::fixed_rule<user::WrongReturn, double> && !integrable<user::WrongReturn>);
static_assert(!qd::fixed_rule<user::NotConst, double> && !integrable<user::NotConst>);
static_assert(!qd::fixed_rule<qd::Romberg, double> && !qd::estimating_method<qd::GaussLegendre<7>, double>);
static_assert(!integrable<qd::GaussLegendre<0>> && !integrable<qd::GaussLegendre<65>> && !integrable<int>);
static_assert(!integrable<user::Both>); // a method is of one kind: named, a type of both kinds is ambiguous
// A braced list is never a method: the named forms take no method argument.
template<typename M>
concept braced_method = requires { qd::integrate<M>(square, 0.0, 1.0, { 3 }); };
static_assert(!braced_method<user::Midpoint>);

int main()
{
    // What the method receives: the ends in order, rel_tol clamped to 1, abs_tol unclamped.
    user::Seen seen;
    auto const spy = qd::integrate(user::Spy{ &seen }, square, 2.0, -1.0, { .abs_tol = 1e300, .rel_tol = 5 });
    std::printf("Spy: a %g, b %g, abs_tol %g, rel_tol %g, value %g\n", seen.a, seen.b, seen.abs_tol, seen.rel_tol, spy->value);
    // The library's validation stops before the method: the spy is not called.
    seen.calls = 0;
    bool const refused = !qd::integrate(user::Spy{ &seen }, square, 0.0, std::nan(""))
                         && !qd::integrate(user::Spy{ &seen }, square, 0.0, 1.0, { .abs_tol = -1 })
                         && qd::integrate(user::Spy{ &seen }, square, 1.0, 1.0)->value == 0.0;
    std::printf("validation before the method: %s, spy calls %d\n", refused ? "yes" : "no", seen.calls);
    // Broken results become non_finite, in every build.
    std::printf("NaN value -> %s; negative abs_error -> %s\n",
        qd::integrate(user::NaNValue{}, square, 0.0, 1.0) ? "accepted" : "non_finite",
        qd::integrate(user::NegativeError{}, square, 0.0, 1.0) ? "accepted" : "non_finite");

    int calls = 0;
    auto const counted = [&calls](auto const &f) {
        return [&calls, f](double x) {
            ++calls;
            return f(x);
        };
    };
    auto const report = [&](char const *name, auto const &r, double exact) {
        if (r)
            std::printf("%-46s calls %5d  error %.2g  abs_error %.2g\n", name, calls, std::fabs(r->value - exact), r->abs_error);
        else
            std::printf("%-46s calls %5d  failed (reason %d)\n", name, calls, static_cast<int>(r.error().reason));
        calls = 0;
    };
    // 1/sqrt x on [0, 1] = 2: the substitution makes it the constant 2, the best case for it.
    auto const inv_sqrt = counted([](double x) { return 1.0 / std::sqrt(x); });
    report("1/sqrt x  GaussKronrod<15>", qd::integrate(inv_sqrt, 0.0, 1.0, { .rel_tol = 1e-12 }), 2.0);
    report("1/sqrt x  SquareAtLowerEnd<GaussKronrod<15>>",
        qd::integrate(user::SquareAtLowerEnd<qd::GaussKronrod<15>>{}, inv_sqrt, 0.0, 1.0, { .rel_tol = 1e-12 }), 2.0);
    report("1/sqrt x  reversed, SquareAtLowerEnd<GaussKronrod<15>>",
        qd::integrate(user::SquareAtLowerEnd<qd::GaussKronrod<15>>{}, inv_sqrt, 1.0, 0.0, { .rel_tol = 1e-12 }), -2.0);
    // cos x / sqrt x on [0, 1] = 2 int_0^1 cos(u^2) du = 1.8090484758005441629 (mpmath, 30 digits): g = 2 cos(u^2).
    double const fresnel = 1.8090484758005441629;
    auto const cos_sqrt = counted([](double x) { return std::cos(x) / std::sqrt(x); });
    report("cos x/sqrt x  GaussKronrod<15>", qd::integrate(cos_sqrt, 0.0, 1.0, { .rel_tol = 1e-12 }), fresnel);
    report("cos x/sqrt x  SquareAtLowerEnd<GaussKronrod<15>>",
        qd::integrate(user::SquareAtLowerEnd<qd::GaussKronrod<15>>{}, cos_sqrt, 0.0, 1.0, { .rel_tol = 1e-12 }), fresnel);
    report("1/sqrt x  SquareAtLowerEnd<Romberg> (samples u = 0)",
        qd::integrate(user::SquareAtLowerEnd<qd::Romberg>{}, inv_sqrt, 0.0, 1.0, { .rel_tol = 1e-12 }), 2.0);
    // fn evaluated only inside [a, b]: the constant 1/4 on [-max, max].
    double const max = std::numeric_limits<double>::max();
    bool inside = true;
    auto const quarter = [&](double x) {
        inside = inside && -max <= x && x <= max;
        return 0.25;
    };
    auto const m = qd::integrate(user::Midpoint{}, quarter, -max, max);
    std::printf("Midpoint on [-max, max]: %s, fn only inside: %s\n", m ? "value" : "non_finite", inside ? "yes" : "no");
}
