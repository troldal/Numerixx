// Spike: methods written outside nxx, against the public protocol only (no nxx::detail), used through integrate().
#include <numerixx/quadrature.hpp>

#include <cmath>
#include <cstdio>

namespace qd = nxx::quadrature;

namespace user {

    // A fixed rule with a run-time parameter: the composite midpoint rule on `panels` equal panels. Orientation comes
    // for free: h = (b - a) / panels is negative for b < a.
    struct Midpoint
    {
        int panels = 100;

        template<nxx::real T, typename F>
        constexpr std::expected<T, qd::error> apply(F &fn, T const &a, T const &b) const
        {
            T const h = T(T(b - a) / T(panels));
            T sum = T(0);
            for (int i = 0; i < panels; ++i) sum = T(sum + T(std::invoke(fn, T(a + T(h * T(T(i) + T(1) / T(2)))))));
            T const r = T(sum * h);
            if (!nxx::math::isfinite(r)) return std::unexpected(qd::error::non_finite);
            return r;
        }
    };

    // A method combinator: x = a + (b - a) u^2 moves an endpoint singularity at a of type (x - a)^-1/2 into a smooth
    // integrand, then any estimating method integrates over u in [0, 1] through the public integrate(). Composition,
    // not a new algorithm.
    template<typename Inner>
    struct SquareAtLowerEnd
    {
        Inner inner{};

        template<nxx::real T, typename F>
        constexpr std::expected<qd::Integral<T>, qd::Failure<T>>
            estimate(F &fn, T const &a, T const &b, qd::Criteria<T> const &criteria) const
        {
            T const width = T(b - a);
            auto const g = [&](T const &u) -> T {
                T const x = T(a + T(width * T(u * u)));
                return T(T(T(2) * width) * T(u * T(std::invoke(fn, x))));
            };
            return qd::integrate(g, T(0), T(1), criteria, inner);
        }
    };

    // Not methods: the wrong return type, and a member that is not const.
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

} // namespace user

constexpr double square(double x) { return x * x; }

// A user fixed rule, in constant expressions, by name and by value: the midpoint rule errs by (b - a) h^2 f'' / 24.
static_assert(qd::integrate<user::Midpoint>(square, 0.0, 3.0).value() > 8.999
              && qd::integrate<user::Midpoint>(square, 0.0, 3.0).value() < 9.0);
static_assert(qd::integrate(square, 0.0, 3.0, user::Midpoint{ .panels = 3 }).value() == 8.75); // 9 - 3 * 1 * 2 / 24
static_assert(qd::integrate(square, 3.0, 0.0, user::Midpoint{ .panels = 3 }).value() == -8.75); // oriented
static_assert(qd::integrate<user::Midpoint>(square, 0.0, std::numeric_limits<double>::infinity())
              == std::unexpected(qd::error::non_finite)); // the library's validation, before the method
// A user combinator over a library method, in constant expressions: sqrt-free, so x^2 on [0, 3] stays exact enough.
static_assert(qd::integrate(square, 0.0, 3.0, {}, user::SquareAtLowerEnd<qd::GaussKronrod<15>>{})->value > 8.999999);
static_assert(qd::integrate<user::SquareAtLowerEnd<qd::Romberg>>(square, 0.0, 3.0)->value < 9.000001);

// The protocol's concepts, and the entry points' constraints.
template<typename M>
concept integrable = requires { qd::integrate<M>(square, 0.0, 1.0); };
static_assert(qd::fixed_rule<user::Midpoint, double(double), double>);
static_assert(qd::estimating_method<user::SquareAtLowerEnd<qd::Romberg>, double(double), double>);
static_assert(!qd::fixed_rule<user::WrongReturn, double(double), double> && !integrable<user::WrongReturn>);
static_assert(!qd::fixed_rule<user::NotConst, double(double), double> && !integrable<user::NotConst>);
static_assert(!qd::fixed_rule<qd::Romberg, double(double), double>); // the two kinds stay apart
static_assert(!qd::estimating_method<qd::GaussLegendre<7>, double(double), double>);
static_assert(!integrable<qd::GaussLegendre<0>> && !integrable<qd::GaussLegendre<65>> && !integrable<int>);

int main()
{
    int calls = 0;
    auto const inv_sqrt = [&calls](double x) {
        ++calls;
        return 1.0 / std::sqrt(x);
    };
    auto const report = [&](char const *name, auto const &r) {
        if (r)
            std::printf("%-34s calls %5d  value %.17g  abs_error %.3g\n", name, calls, r->value, r->abs_error);
        else
            std::printf("%-34s calls %5d  failed (reason %d)\n", name, calls, static_cast<int>(r.error().reason));
        calls = 0;
    };
    // 1/sqrt x on [0, 1] = 2: GaussKronrod<15> directly, and through the user's substitution (smooth: 2 u / u).
    report("GaussKronrod<15>", qd::integrate(inv_sqrt, 0.0, 1.0, { .rel_tol = 1e-12 }));
    report("SquareAtLowerEnd<GaussKronrod<15>>",
        qd::integrate(inv_sqrt, 0.0, 1.0, { .rel_tol = 1e-12 }, user::SquareAtLowerEnd<qd::GaussKronrod<15>>{}));
    report("SquareAtLowerEnd<Romberg>",
        qd::integrate(inv_sqrt, 0.0, 1.0, { .rel_tol = 1e-12 }, user::SquareAtLowerEnd<qd::Romberg>{}));
    auto const value = qd::integrate(inv_sqrt, 0.0, 1.0, user::Midpoint{ .panels = 1000 });
    std::printf("%-34s calls %5d  value %.17g\n", "Midpoint{1000}", calls, *value);
}
