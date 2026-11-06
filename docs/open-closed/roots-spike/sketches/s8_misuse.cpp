// Sketch 8: misuse. Build with -DCASE=n; CASE=0 runs the run-time cases.
#include <numerixx/roots.hpp>
#include <print>

namespace rt = nxx::roots;

namespace user {

    struct Bisect // a correct base, varied below
    {
        template<nxx::real T> struct State { T a, fa, b, fb; rt::Step<T> step; };
        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
        { return State<T>{ s.lower, s.at_lower[0], s.upper, s.at_upper[0], { s.lower, s.at_lower[0], s.upper, s.at_upper[0], s.lower } }; }
        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &) const
        {
            T const m = T(T(s.a / T(2)) + T(s.b / T(2)));
            T const fm = T(std::invoke(fn, m));
            if ((fm < T(0)) == (s.fa < T(0))) { s.a = m; s.fa = fm; } else { s.b = m; s.fb = fm; }
            s.step = rt::Step<T>{ s.a, s.fa, s.b, s.fb, s.a };
            return {};
        }
    };

    // Run time: keeps the wrong half (loses the sign change).
    struct LosesSign : Bisect
    {
        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &) const
        {
            T const m = T(T(s.a / T(2)) + T(s.b / T(2)));
            T const fm = T(std::invoke(fn, m));
            if ((fm < T(0)) != (s.fa < T(0))) { s.a = m; s.fa = fm; } else { s.b = m; s.fb = fm; } // wrong half
            s.step = rt::Step<T>{ s.a, s.fa, s.b, s.fb, s.a };
            return {};
        }
    };

    // Run time: reports its ends in the wrong order.
    struct Unordered : Bisect
    {
        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
        { return State<T>{ s.lower, s.at_lower[0], s.upper, s.at_upper[0], { s.upper, s.at_upper[0], s.lower, s.at_lower[0], s.lower } }; }
    };

    // Run time: point outside its points (the midpoint forgotten: point 0).
    struct PointOutside : Bisect
    {
        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
        { return State<T>{ s.lower, s.at_lower[0], s.upper, s.at_upper[0], { s.lower, s.at_lower[0], s.upper, s.at_upper[0], T(0) } }; }
    };

    // Run time: evaluates outside the caller's bracket (extrapolates).
    struct Extrapolates : Bisect
    {
        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &) const
        {
            T const x = T(s.b + T(1));
            T const fx = T(std::invoke(fn, x));
            s.b = x; s.fb = fx;
            s.step = rt::Step<T>{ s.a, s.fa, s.b, s.fb, s.a };
            return {};
        }
    };

    // Run time: never narrows its bracket and never reports exhausted (hangs). Guarded by an observer in main.
    struct Stuck : Bisect
    {
        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &, State<T> &, rt::Tolerance<T> const &) const { return {}; }
    };

    // CASE 1: advance not const.
    struct NonConst : Bisect
    {
        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &t)
        { return Bisect::advance(fn, s, t); }
    };

    // CASE 5: no step member.
    struct NoStep
    {
        template<nxx::real T> struct State { T a, fa, b, fb; };
        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
        { return State<T>{ s.lower, s.at_lower[0], s.upper, s.at_upper[0] }; }
        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &, State<T> &, rt::Tolerance<T> const &) const { return {}; }
    };

    // CASE 6: start returns the state, not an expected.
    struct PlainStart : Bisect
    {
        template<nxx::real T, typename F>
        constexpr State<T> start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
        { return State<T>{ s.lower, s.at_lower[0], s.upper, s.at_upper[0], { s.lower, s.at_lower[0], s.upper, s.at_upper[0], s.lower } }; }
    };

    // CASE 7: advance returns std::expected<void, rt::error> but takes the state by value.
    struct ByValue : Bisect
    {
        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &, State<T>, rt::Tolerance<T> const &) const { return {}; }
    };

    // CASE 8: start written for double, not a template on T (works for double only?).
    struct DoubleOnly : Bisect
    {
        template<typename F>
        constexpr std::expected<State<double>, rt::error> start(F &, rt::Start<double> const &s, rt::Tolerance<double> const &) const
        { return State<double>{ s.lower, s.at_lower[0], s.upper, s.at_upper[0], { s.lower, s.at_lower[0], s.upper, s.at_upper[0], s.lower } }; }
    };
} // namespace user

auto const f = [](double x) { return x * x - 2.0; };

int main()
{
#if CASE == 0
    auto show = [](char const *name, std::expected<rt::Root<double>, rt::error> const &r) {
        if (r) std::println("{:34} {:.17g} +- {:.3g}", name, r->value, r->abs_error);
        else std::println("{:34} error {}", name, r.error());
    };
    show("Bisect (correct)", rt::solve<user::Bisect>(f, rt::Bracket{ 1.0, 2.0 }));
    show("LosesSign", rt::solve<user::LosesSign>(f, rt::Bracket{ 1.0, 2.0 }));
    show("Unordered", rt::solve<user::Unordered>(f, rt::Bracket{ 1.0, 2.0 }));
    show("PointOutside", rt::solve<user::PointOutside>(f, rt::Bracket{ 1.0, 2.0 }));
    int evals_outside = 0;
    auto watched = [&](double x) { if (x < 1.0 || x > 2.0) ++evals_outside; return f(x); };
    show("Extrapolates", rt::solve<user::Extrapolates>(watched, rt::Bracket{ 1.0, 2.0 }));
    std::println("  evaluations outside [1, 2] before the error: {}", evals_outside);
    int n = 0;
    show("Stuck (observer stops at 1e6)", rt::solve<user::Stuck>(rt::observer([&](rt::Step<double> const &) {
        if (++n == 1000000) { std::println("  Stuck: 1e6 steps, no end in sight"); std::exit(0); }
    }), f, rt::Bracket{ 1.0, 2.0 }));
#elif CASE == 1
    (void)rt::solve<user::NonConst>(f, rt::Bracket{ 1.0, 2.0 });
#elif CASE == 2
    (void)rt::solve<rt::Brent>(rt::observer([](double x) { (void)x; }), f, rt::Bracket{ 1.0, 2.0 });
#elif CASE == 3
    (void)rt::solve<rt::Brent>(rt::observer([](rt::Step<double> const &) { return true; }), f, rt::Bracket{ 1.0, 2.0 });
#elif CASE == 4
    (void)rt::solve<rt::Brent>(f, 1.0);
#elif CASE == 5
    (void)rt::solve<user::NoStep>(f, rt::Bracket{ 1.0, 2.0 });
#elif CASE == 6
    (void)rt::solve<user::PlainStart>(f, rt::Bracket{ 1.0, 2.0 });
#elif CASE == 7
    (void)rt::solve<user::ByValue>(f, rt::Bracket{ 1.0, 2.0 });
#elif CASE == 8
    (void)rt::solve<user::DoubleOnly>(f, rt::Bracket{ 1.0, 2.0 });
    (void)rt::solve<user::DoubleOnly>([](float x) { return x * x - 2.0f; }, rt::Bracket{ 1.0f, 2.0f });
#elif CASE == 9
    auto const print = [](rt::Step<double> const &) {};
    (void)rt::solve<rt::Brent>(print, f, rt::Bracket{ 1.0, 2.0 }); // forgot observer()
#elif CASE == 10
    static_assert(rt::bracketing_method<user::NonConst, double>);
#elif CASE == 11
    (void)rt::solve<rt::Brent>(rt::observer([](rt::Step<float> const &) {}), f, rt::Bracket{ 1.0, 2.0 }); // wrong T
#endif
}
