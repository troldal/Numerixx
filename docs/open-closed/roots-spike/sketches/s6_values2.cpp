// Sketch 6: a user method with values = 2, through the separate-derivative form solve<M>(fn, dfn, bracket).
// A safeguarded Newton on a bracket: Newton's step from the end with the smaller |f| where it lands strictly inside
// (and tol away from the ends), else the midpoint.
#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>
#include <cmath>
#include <print>

namespace rt = nxx::roots;

namespace user {

    struct NewtonBisect
    {
        static constexpr std::size_t values = 2;

        template<nxx::real T>
        struct State
        {
            T a, fa, da, b, fb, db;
            rt::Step<T> step;
        };

        template<nxx::real T>
        static constexpr rt::Step<T> report(State<T> const &s)
        { return rt::Step<T>{ s.a, s.fa, s.b, s.fb, nxx::math::abs(s.fa) <= nxx::math::abs(s.fb) ? s.a : s.b }; }

        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &, rt::Start<T, 2> const &s, rt::Tolerance<T> const &) const
        {
            State<T> st{ s.lower, s.at_lower[0], s.at_lower[1], s.upper, s.at_upper[0], s.at_upper[1], {} };
            st.step = report(st);
            return st;
        }

        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &t) const
        {
            bool const from_a = nxx::math::abs(s.fa) <= nxx::math::abs(s.fb);
            T const x0 = from_a ? s.a : s.b, f0 = from_a ? s.fa : s.fb, d0 = from_a ? s.da : s.db;
            T const mid = T(T(s.a / T(2)) + T(s.b / T(2)));
            T const margin = t.at(x0);
            T x = mid;
            if (d0 != T(0) && nxx::math::isfinite(d0) && T(T(s.b - s.a) / T(4)) > margin) {
                T const n = T(x0 - T(f0 / d0));
                if (T(s.a + margin) <= n && n <= T(s.b - margin)) x = n;
                else if (s.a <= n && n < T(s.a + margin)) x = T(s.a + margin);
                else if (T(s.b - margin) < n && n <= s.b) x = T(s.b - margin);
            }
            if (!(s.a < x && x < s.b)) { s.step.exhausted = true; return {}; }
            auto const v = std::invoke(fn, x); // std::array<T, 2> from the separate form; any tuple-like otherwise
            T const fx = T(std::get<0>(v)), dx = T(std::get<1>(v));
            if (!nxx::math::isfinite(fx)) return std::unexpected(rt::error::non_finite);
            if (fx == T(0)) { s.step = rt::Step<T>{ x, fx, x, fx, x }; return {}; }
            if ((fx < T(0)) == (s.fa < T(0))) { s.a = x; s.fa = fx; s.da = dx; }
            else { s.b = x; s.fb = fx; s.db = dx; }
            s.step = report(s);
            return {};
        }
    };

} // namespace user

int main()
{
    auto const f = [](double x) { return x * x * x - 2.0 * x - 5.0; };
    auto const df = [](double x) { return 3.0 * x * x - 2.0; };
    auto const fdf = [](double x) { return std::pair{ x * x * x - 2.0 * x - 5.0, 3.0 * x * x - 2.0 }; };
    auto const fdf3 = [](double x) { return std::tuple{ x * x * x - 2.0 * x - 5.0, 3.0 * x * x - 2.0, 6.0 * x }; };
    double const exact = 2.0945514815423265;
    int calls = 0;
    auto const counted_f = [&](double x) { ++calls; return f(x); };
    auto const show = [&](char const *name, std::expected<rt::Root<double>, rt::error> const &r) {
        if (r) std::println("{:34} {:.17g} +- {:.3g}  |x - exact| = {:.3g}  calls {}", name, r->value, r->abs_error, std::abs(r->value - exact), calls);
        else std::println("{:34} error {}", name, r.error());
        calls = 0;
    };
    show("NewtonBisect, separate df", rt::solve<user::NewtonBisect>(counted_f, df, rt::Bracket{ 2.0, 3.0 }));
    show("NewtonBisect, derivative_of", rt::solve<user::NewtonBisect>(counted_f, nxx::deriv::derivative_of(counted_f), rt::Bracket{ 2.0, 3.0 }));
    show("Newton (lib), separate df", rt::solve<rt::Newton>(counted_f, df, rt::Bracket{ 2.0, 3.0 }));
    show("NewtonBisect, pair fn", rt::solve<user::NewtonBisect>(fdf, rt::Bracket{ 2.0, 3.0 }));
    show("NewtonBisect, 3-tuple fn", rt::solve<user::NewtonBisect>(fdf3, rt::Bracket{ 2.0, 3.0 }));
    show("NewtonBisect as value, pair fn", rt::solve(user::NewtonBisect{}, fdf, rt::Bracket{ 2.0, 3.0 }));
#ifdef VALUE_SEPARATE
    show("value + separate df", rt::solve(user::NewtonBisect{}, f, df, rt::Bracket{ 2.0, 3.0 }));
#endif
#ifdef PLAIN_FN
    show("values=2 with a plain fn", rt::solve<user::NewtonBisect>(f, rt::Bracket{ 2.0, 3.0 }));
#endif
}
