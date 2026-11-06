// Sketch 4: an open method of one's own (Chebyshev's third-order step, values = 3), handing over to the library's
// Brent once two iterates bracket the root. Public API only. MinStep = false: the step as the textbook writes it.
#include <numerixx/roots.hpp>
#include <cmath>
#include <optional>
#include <print>
#include <tuple>

namespace rt = nxx::roots;

namespace user {

    template<bool MinStep = true>
    struct Chebyshev
    {
        static constexpr std::size_t values = 3;

        template<nxx::real T>
        struct State
        {
            T x, f, d, d2;   // the latest iterate and fn's three values there
            T prev, fprev;   // the iterate before
            std::optional<rt::bracketing_state_t<rt::Brent, T>> brent; // after the hand-over
            rt::Step<T> step;
        };

        template<nxx::real T>
        static constexpr rt::Step<T> iterates(T const &a, T const &fa, T const &b, T const &fb)
        { return a < b ? rt::Step<T>{ a, fa, b, fb, b } : rt::Step<T>{ b, fb, a, fa, b }; } // point: b, the latest

        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &fn, T const &guess, rt::Tolerance<T> const &) const
        {
            auto const [f, d, d2] = std::invoke(fn, guess); // what if fn gives 4 values, or a Joined callable?
            if (!nxx::math::isfinite(T(f))) return std::unexpected(rt::error::non_finite);
            return State<T>{ guess, T(f), T(d), T(d2), guess, T(f), std::nullopt, iterates(guess, T(f), guess, T(f)) };
        }

        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &t) const
        {
            if (s.brent) {
                auto const r = rt::Brent{}.advance(fn, *s.brent, t);
                s.step = s.brent->step; // copy the finisher's report into ours
                return r;
            }
            if (s.d == T(0) || !nxx::math::isfinite(s.d) || !nxx::math::isfinite(s.d2)) {
                s.step.exhausted = true;
                return {};
            }
            T const u = T(s.f / s.d);
            T dx = T(u * T(T(1) + T(T(u * s.d2) / T(T(2) * s.d))));
            if constexpr (MinStep) { // at least tol(x), so that an iterate converging from one side crosses the root
                T const tol = t.at(s.x);
                if (nxx::math::abs(dx) < tol) dx = dx < T(0) ? T(-tol) : tol;
            }
            T const x = T(s.x - dx);
            if (!nxx::math::isfinite(x)) { s.step.exhausted = true; return {}; }
            auto const [f, d, d2] = std::invoke(fn, x);
            if (!nxx::math::isfinite(T(f))) return std::unexpected(rt::error::non_finite);
            s.prev = s.x; s.fprev = s.f;
            s.x = x; s.f = T(f); s.d = T(d); s.d2 = T(d2);
            s.step = iterates(s.prev, s.fprev, s.x, s.f);
            if (s.f == T(0)) { s.step = rt::Step<T>{ s.x, s.f, s.x, s.f, s.x }; return {}; }
            if ((s.f < T(0)) != (s.fprev < T(0))) { // a sign change: hand over to Brent on the two iterates
                rt::Start<T> const b{ s.step.lower, { s.step.f_lower }, s.step.upper, { s.step.f_upper } };
                auto started = rt::Brent{}.start(fn, b, t);
                if (!started) return std::unexpected(started.error());
                s.brent = *started;
                s.step = s.brent->step;
            }
            return {};
        }
    };

} // namespace user

static_assert(rt::open_method<user::Chebyshev<>, double>);
static_assert(rt::values_of<user::Chebyshev<>> == 3);

template<typename M, typename F>
void run(char const *name, F const &f, double guess, double exact, bool trace = false)
{
    int calls = 0, reports = 0;
    auto counted = [&](double x) { ++calls; return f(x); };
    auto const print = [&](rt::Step<double> const &s) {
        ++reports;
        if (trace) std::println("    [{:.17g}, {:.17g}] f=({:.3g}, {:.3g}) point {:.17g}", s.lower, s.upper, s.f_lower, s.f_upper, s.point);
    };
    auto const r = rt::solve<M>(rt::observer(print), counted, guess);
    if (r) std::println("  {:26} {:.17g} +- {:.3g}  |x - exact| = {:.3g}  calls {:3}  {}", name, r->value, r->abs_error,
                        std::abs(r->value - exact), calls, std::abs(r->value - exact) <= r->abs_error ? "ok" : "BOUND FAILS");
    else std::println("  {:26} error {} (calls {}, reports {})", name, r.error(), calls, reports);
}

int main()
{
    auto const sq3 = [](double x) { return std::tuple{ x * x - 2.0, 2.0 * x, 2.0 }; };
    auto const cubic3 = [](double x) { return std::array{ x * x * x - 2.0 * x - 5.0, 3.0 * x * x - 2.0, 6.0 * x }; };
    auto const cos3 = [](double x) { return std::tuple{ std::cos(x) - x, -std::sin(x) - 1.0, -std::cos(x) }; };
    std::println("x^2 - 2 from 3 (Newton-like iterates converge from above):");
    run<user::Chebyshev<false>>("Chebyshev, textbook step", sq3, 3.0, std::sqrt(2.0));
    run<user::Chebyshev<true>>("Chebyshev, min step tol", sq3, 3.0, std::sqrt(2.0), true);
    run<rt::Halley>("Halley (lib)", sq3, 3.0, std::sqrt(2.0));
    std::println("cubic from 2:");
    run<user::Chebyshev<false>>("Chebyshev, textbook step", cubic3, 2.0, 2.0945514815423265);
    run<user::Chebyshev<true>>("Chebyshev, min step tol", cubic3, 2.0, 2.0945514815423265, true);
    run<rt::Halley>("Halley (lib)", cubic3, 2.0, 2.0945514815423265);
    std::println("cos x - x from 0:");
    run<user::Chebyshev<true>>("Chebyshev, min step tol", cos3, 0.0, 0.73908513321516064);
    run<rt::Halley>("Halley (lib)", cos3, 0.0, 0.73908513321516064);
    std::println("on a bracket? {}", rt::bracketing_method<user::Chebyshev<>, double>);
}
