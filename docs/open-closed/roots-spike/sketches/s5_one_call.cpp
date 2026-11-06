// Sketch 5: one-call methods. (a) bisection by hand, on a bracket; (b) a combinator from a guess: search() for a
// bracket, then any bracketing method; (c) the observer given to a one-call method (expected: no overload).
#include <numerixx/roots.hpp>
#include <cmath>
#include <print>

namespace rt = nxx::roots;

namespace user {

    // (a) Bisection by hand, whole run in one call.
    struct HandBisection
    {
        template<nxx::real T, typename F>
        constexpr std::expected<rt::Root<T>, rt::error> solve(F &fn, rt::Start<T> const &s, rt::Tolerance<T> const &t) const
        {
            T a = s.lower, fa = s.at_lower[0], b = s.upper;
            while (true) {
                T const m = T(T(a / T(2)) + T(b / T(2)));
                T const half = T(T(b / T(2)) - T(a / T(2)));
                if (half <= t.at(m) || !(a < m && m < b)) return rt::Root<T>{ m, half };
                T const fm = T(std::invoke(fn, m));
                if (!nxx::math::isfinite(fm)) return std::unexpected(rt::error::non_finite);
                if (fm == T(0)) return rt::Root<T>{ m, t.at(m) };
                if ((fm < T(0)) == (fa < T(0))) { a = m; fa = fm; }
                else b = m;
            }
        }
    };

    // (b) From a guess: search() for a bracket, then a bracketing method. Tolerance has to go back to Criteria.
    template<typename Finisher = rt::Brent>
    struct SearchThen
    {
        Finisher finisher{};

        template<nxx::real T, typename F>
        constexpr std::expected<rt::Root<T>, rt::error> solve(F &fn, T const &guess, rt::Tolerance<T> const &t) const
        {
            rt::Criteria<T> const c{ .abs_tol = t.absolute, .rel_tol = t.relative };
            return rt::search(fn, guess).and_then([&](rt::Bracket<T> const &b) { return rt::solve(finisher, fn, b, c); });
        }
    };

} // namespace user

static_assert(rt::bracketing_solver<user::HandBisection, double>);
static_assert(!rt::bracketing_method<user::HandBisection, double>);
static_assert(rt::open_solver<user::SearchThen<>, double>);

int main()
{
    auto const cubic = [](double x) { return x * x * x - 2.0 * x - 5.0; };
    double const exact = 2.0945514815423265;
    auto const show = [&](char const *name, std::expected<rt::Root<double>, rt::error> const &r) {
        if (r) std::println("{:28} {:.17g} +- {:.3g}  |x - exact| = {:.3g}", name, r->value, r->abs_error, std::abs(r->value - exact));
        else std::println("{:28} error {}", name, r.error());
    };
    show("HandBisection [2, 3]", rt::solve<user::HandBisection>(cubic, rt::Bracket{ 2.0, 3.0 }));
    show("Bisection (lib) [2, 3]", rt::solve<rt::Bisection>(cubic, rt::Bracket{ 2.0, 3.0 }));
    show("HandBisection rel 1e-6", rt::solve<user::HandBisection>(cubic, rt::Bracket{ 2.0, 3.0 }, { .rel_tol = 1e-6 }));
    show("SearchThen<Brent> from 10", rt::solve<user::SearchThen<>>(cubic, 10.0));
    show("SearchThen{Toms748} from -3", rt::solve(user::SearchThen<rt::Toms748>{}, cubic, -3.0));
    show("SearchThen<HandBisection> 10", rt::solve<user::SearchThen<user::HandBisection>>(cubic, 10.0));
    // Tolerance round trip into Criteria: does solve() accept the same request back?
    show("SearchThen, abs 1e-3", rt::solve<user::SearchThen<>>(cubic, 10.0, { .abs_tol = 1e-3, .rel_tol = 0.0 }));
#ifdef OBSERVE
    (void)rt::solve<user::HandBisection>(rt::observer([](rt::Step<double> const &) {}), cubic, rt::Bracket{ 2.0, 3.0 });
#endif
}
