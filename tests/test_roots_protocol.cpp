// The method protocol of nxx::roots on a bracket (Phase 7, step 1, docs/roots/PHASE7_API.md): the loop's checks with
// a scripted method, observing, the methods as values, Bisection's stop test, the fix for types without subnormals,
// methods of one's own and the concepts.
//
// References are closed forms; the scripted cases are derived by hand from the loop's rules (PHASE7_API.md B).

#include "roots_check.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <vector>

using namespace roots_check;

namespace {

    // ---- A scripted method: its steps come from a list ------------------------------------------------------------

    // advance past the list's end returns no_sign_change, an error no loop branch produces, so a script too short
    // shows in a test; fail_at makes that advance return not_converged instead.
    template<std::size_t N>
    struct Script
    {
        std::array<rt::Step<double>, N> steps;
        int fail_at = -1;

        struct State
        {
            std::array<rt::Step<double>, N> steps;
            std::size_t next;
            int fail_at;
            rt::Step<double> step;
        };

        template<typename F>
        constexpr std::expected<State, rt::error> start(F &, rt::Start<double> const &, rt::Tolerance<double> const &) const
        { return State{ steps, 1, fail_at, steps[0] }; }

        template<typename F>
        constexpr std::expected<void, rt::error> advance(F &, State &s, rt::Tolerance<double> const &) const
        {
            if (static_cast<int>(s.next) == s.fail_at) return std::unexpected(rt::error::not_converged);
            if (s.next >= N) return std::unexpected(rt::error::no_sign_change);
            s.step = s.steps[s.next++];
            return {};
        }
    };

    constexpr auto half_off = [](double x) { return x - 0.5; }; // [0, 1] is a valid bracket, 0 not inside it
    constexpr rt::Bracket<double> unit{ 0.0, 1.0 };

    // abs_tol 1e-3, rel_tol 2 eps: tol(x) = 1e-3 + 2 eps |x|, so a bracket 1e-3 wide meets the stop test.
    template<std::size_t N>
    constexpr std::expected<rt::Root<double>, rt::error> on_unit(Script<N> const &m)
    { return rt::solve(m, half_off, unit, { .abs_tol = 1e-3 }); }

    constexpr rt::Step<double> s0{ 0.0, -0.5, 1.0, 0.5, 0.5 };
    constexpr rt::Step<double> s1{ 0.25, -0.25, 0.75, 0.25, 0.75 };
    constexpr rt::Step<double> s2{ 0.4995, -0.0005, 0.5005, 0.0005, 0.5005 }; // about 1e-3 wide: meets the test

    // Each case of the loop (PHASE7_API.md B), true where the loop does what the comment says; checked in constant
    // evaluation and again at run time.
    constexpr std::array<bool, 20> loop_cases()
    {
        using rt::error;
        using rt::Step;
        return {
            // The first step that meets the request is returned, abs_error the larger distance to an end: the width
            // at an end point.
            on_unit(Script<3>{ { s0, s1, s2 } })->value == 0.5005,
            on_unit(Script<3>{ { s0, s1, s2 } })->abs_error == rt::detail::distance_bound(0.4995, 0.5005),
            // An interior point: the larger of its two distances to the ends.
            on_unit(Script<1>{ { Step<double>{ 0.4995, -0.0005, 0.5005, 0.0005, 0.5 } } })->abs_error
                == std::max(rt::detail::distance_bound(0.5, 0.4995), rt::detail::distance_bound(0.5005, 0.5)),
            // A point that is an end where fn is exactly 0: that point, abs_error tol(point).
            on_unit(Script<2>{ { s0, Step<double>{ 0.25, -0.25, 0.5, 0.0, 0.5 } } })
                == rt::Root<double>{ 0.5, 1e-3 + 2.0 * deps * 0.5 },
            // A NaN or an infinity anywhere in a step: non_finite.
            on_unit(Script<2>{ { s0, Step<double>{ 0.25, -0.25, dinf, 0.25, 0.5 } } }).error() == error::non_finite,
            on_unit(Script<2>{ { s0, Step<double>{ 0.25, dnan, 0.75, 0.25, 0.5 } } }).error() == error::non_finite,
            // The point outside its points; a step outside the caller's bracket; a lost sign change: invalid_step.
            on_unit(Script<2>{ { s0, Step<double>{ 0.25, -0.25, 0.75, 0.25, 0.8 } } }).error() == error::invalid_step,
            on_unit(Script<2>{ { s0, Step<double>{ 0.25, -0.25, 1.5, 0.25, 0.75 } } }).error() == error::invalid_step,
            on_unit(Script<2>{ { s0, Step<double>{ 0.25, 0.25, 0.75, 0.25, 0.75 } } }).error() == error::invalid_step,
            // The start is checked as every step is.
            on_unit(Script<1>{ { Step<double>{ -0.5, -1.0, 1.0, 0.5, 0.5 } } }).error() == error::invalid_step,
            // Opposite signs at a single point (lower == upper).
            on_unit(Script<2>{ { s0, Step<double>{ 0.5, -0.25, 0.5, 0.25, 0.5 } } }).error() == error::invalid_step,
            // The ends in the wrong order: no point lies between them.
            on_unit(Script<2>{ { s0, Step<double>{ 0.75, 0.25, 0.25, -0.25, 0.5 } } }).error() == error::invalid_step,
            // fn exactly 0 at an end that is not the point: no sign change shown, so invalid_step.
            on_unit(Script<2>{ { s0, Step<double>{ 0.25, -0.25, 0.5, 0.0, 0.25 } } }).error() == error::invalid_step,
            // Progress: the same bracket twice, or one not nested in the last, is invalid_step (numerics review of
            // draft 1, finding 3); one end kept and the other moved inside is nested.
            on_unit(Script<3>{ { s0, s1, s1 } }).error() == error::invalid_step,
            on_unit(Script<3>{ { s0, s1, Step<double>{ 0.2, -0.3, 0.6, 0.1, 0.6 } } }).error() == error::invalid_step,
            on_unit(Script<3>{ { s0, Step<double>{ 0.0, -0.5, 0.75, 0.25, 0.75 }, s2 } })->value == 0.5005,
            // A point moved without moving the bracket is the same bracket: invalid_step.
            on_unit(Script<2>{ { s0, Step<double>{ 0.0, -0.5, 1.0, 0.5, 0.0 } } }).error() == error::invalid_step,
            // An error from advance is returned as it is.
            on_unit(Script<3>{ { s0, s1, s2 }, 1 }).error() == error::not_converged,
            on_unit(Script<1>{ { s0 } }).error() == error::no_sign_change,
            // The caller's bracket is the validated one: on [-1, 1] zero first halves it to [0, 1] (fn(0) = -0.5),
            // so a step reaching below 0 is outside it.
            rt::solve(Script<2>{ { s0, Step<double>{ -0.5, -1.0, 0.75, 0.25, 0.75 } } }, half_off,
                rt::Bracket{ -1.0, 1.0 }, { .abs_tol = 1e-3 })
                    .error()
                == error::invalid_step,
        };
    }

    static_assert(std::ranges::all_of(loop_cases(), [](bool b) { return b; }));

    // ---- Observing changes nothing, for every method on a bracket, named and as a value -----------------------------

    // Lambdas, not functions: MSVC warns that const on a function type (F const & with F a function) is ignored.
    constexpr auto g = [](double x) { return x * x - 2.0; };
    constexpr auto wal_fn = [](double x) { return (x * x - 2.0) * x - 5.0; }; // Wallis's x^3 - 2x - 5
    constexpr auto g1d = [](double x) { return std::array<double, 2>{ x * x - 2.0, 2.0 * x }; };
    constexpr auto g2d = [](double x) { return std::array<double, 3>{ x * x - 2.0, 2.0 * x, 2.0 }; };
    constexpr rt::Bracket<double> one_two{ 1.0, 2.0 };

    // The named form, the value form and both observed forms give the same result; reports counts the observer's
    // calls of one observed run.
    template<typename M, typename F>
    constexpr bool same_everywhere(F const &fn)
    {
        int reports = 0;
        auto count = [&reports](rt::Step<double> const &) { ++reports; };
        auto const named = rt::solve<M>(fn, one_two);
        auto const value = rt::solve(M{}, fn, one_two);
        auto const observed = rt::solve<M>(rt::observer(count), fn, one_two);
        int const once = reports;
        auto const observed_value = rt::solve(M{}, rt::observer(count), fn, one_two);
        return named.has_value() && named == value && named == observed && named == observed_value && once > 0
               && reports == 2 * once;
    }

    static_assert(same_everywhere<rt::Bisection>(g));
    static_assert(same_everywhere<rt::Brent>(g));
    static_assert(same_everywhere<rt::Toms748>(g));
    static_assert(same_everywhere<rt::AndersonBjorck>(g));
    static_assert(same_everywhere<rt::Ridders>(g));
    static_assert(same_everywhere<rt::RegulaFalsi>(g));
    static_assert(same_everywhere<rt::Illinois>(g));
    static_assert(same_everywhere<rt::Newton>(g1d));
    static_assert(same_everywhere<rt::Halley>(g2d));

    // A report after the start and one after each advance, each advance one evaluation where fn is never exactly 0:
    // reports = evaluations - 2 + 1 (the two ends are evaluated by solve()).
    template<typename M, typename F>
    constexpr bool one_report_per_evaluation(F const &fn)
    {
        int reports = 0;
        int evaluations = 0;
        auto counted = [&evaluations, &fn](double const &x) {
            ++evaluations;
            return fn(x);
        };
        auto count = [&reports](rt::Step<double> const &) { ++reports; };
        (void)rt::solve<M>(rt::observer(count), counted, one_two);
        return reports == evaluations - 1;
    }

    static_assert(one_report_per_evaluation<rt::Bisection>(g));
    static_assert(one_report_per_evaluation<rt::Brent>(g));
    static_assert(one_report_per_evaluation<rt::Toms748>(g));
    static_assert(one_report_per_evaluation<rt::AndersonBjorck>(g));
    static_assert(one_report_per_evaluation<rt::Ridders>(g));
    static_assert(one_report_per_evaluation<rt::Illinois>(g));

    // An observer with a non-const call operator and a count the caller reads.
    struct Tally
    {
        int n = 0;
        constexpr void operator()(rt::Step<double> const &) { ++n; }
    };

    // The number of steps of M on x^2 - 2 over [1, 2], from an lvalue recorder.
    template<typename M>
    std::size_t steps_of()
    {
        std::vector<rt::Step<double>> steps;
        auto record = [&steps](rt::Step<double> const &s) { steps.push_back(s); };
        (void)rt::solve<M>(rt::observer(record), g, one_two);
        return steps.size();
    }

    // ---- Methods of one's own ---------------------------------------------------------------------------------------

    // Pegasus (Dowell and Jarratt 1972): false position whose retained value is scaled by f1 / (f1 + f2), each new
    // point kept tol(point) inside the ends, as bracketing_method advises for convergence.
    struct Pegasus
    {
        template<nxx::real T>
        struct State
        {
            T x0, f0, g0; // the retained point, fn there, its scaled value
            T x1, f1;     // the latest point
            rt::Step<T> step;
        };

        template<nxx::real T>
        static constexpr rt::Step<T> report(State<T> const &s)
        {
            T const point = nxx::math::abs(s.f1) <= nxx::math::abs(s.f0) ? s.x1 : s.x0;
            return s.x0 < s.x1 ? rt::Step<T>{ s.x0, s.f0, s.x1, s.f1, point } : rt::Step<T>{ s.x1, s.f1, s.x0, s.f0, point };
        }

        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
        {
            State<T> st{ s.lower, s.at_lower[0], s.at_lower[0], s.upper, s.at_upper[0], {} };
            st.step = report(st);
            return st;
        }

        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &t) const
        {
            T const lo = s.step.lower;
            T const hi = s.step.upper;
            T const mid = T(T(lo / T(2)) + T(hi / T(2)));
            T const q = T(s.f1 / T(s.f1 - s.g0)); // in (0, 1): f1 and g0 have opposite signs
            T x = T(s.x1 - T(q * T(s.x1 - s.x0)));
            T const margin = t.at(s.step.point);
            if (T(T(hi - lo) / T(4)) <= margin)
                x = mid;
            else if (x < T(lo + margin))
                x = T(lo + margin);
            else if (x > T(hi - margin))
                x = T(hi - margin);
            if (!(lo < x && x < hi)) return std::unexpected(rt::error::not_converged);
            T const fx = T(std::invoke(fn, x));
            if (!nxx::math::isfinite(fx)) return std::unexpected(rt::error::non_finite);
            if (fx == T(0)) {
                s.step = rt::Step<T>{ lo, s.step.f_lower, x, fx, x };
                return {};
            }
            if ((fx < T(0)) != (s.f1 < T(0))) {
                s.x0 = s.x1;
                s.f0 = s.f1;
                s.g0 = s.f1;
            } else {
                s.g0 = T(s.g0 * T(s.f1 / T(s.f1 + fx)));
            }
            s.x1 = x;
            s.f1 = fx;
            s.step = report(s);
            return {};
        }
    };

    // A section method with a parameter: the bracket cut at a fraction of its width from the end with the smaller
    // |f|. No default: a method with a required parameter is given as a value.
    struct Section
    {
        double fraction;

        explicit constexpr Section(double f) : fraction(f) {}

        struct State
        {
            double a, fa, b, fb;
            rt::Step<double> step;
        };

        static constexpr rt::Step<double> report(double a, double fa, double b, double fb)
        { return rt::Step<double>{ a, fa, b, fb, nxx::math::abs(fa) < nxx::math::abs(fb) ? a : b }; }

        template<typename F>
        constexpr std::expected<State, rt::error> start(F &, rt::Start<double> const &s, rt::Tolerance<double> const &) const
        {
            return State{ s.lower, s.at_lower[0], s.upper, s.at_upper[0],
                report(s.lower, s.at_lower[0], s.upper, s.at_upper[0]) };
        }

        template<typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State &s, rt::Tolerance<double> const &t) const
        {
            double const width = s.b - s.a;
            double x = s.step.point == s.a ? s.a + fraction * width : s.b - fraction * width;
            double const margin = t.at(s.step.point);
            if (width / 4 <= margin) x = s.a / 2 + s.b / 2;
            else if (x < s.a + margin) x = s.a + margin;
            else if (x > s.b - margin) x = s.b - margin;
            if (!(s.a < x && x < s.b)) return std::unexpected(rt::error::not_converged);
            double const fx = std::invoke(fn, x);
            if (!nxx::math::isfinite(fx)) return std::unexpected(rt::error::non_finite);
            if (fx == 0.0) {
                s.step = rt::Step<double>{ s.a, s.fa, x, fx, x };
                return {};
            }
            if ((fx < 0.0) == (s.fa < 0.0)) {
                s.a = x;
                s.fa = fx;
            } else {
                s.b = x;
                s.fb = fx;
            }
            s.step = report(s.a, s.fa, s.b, s.fb);
            return {};
        }
    };

    // A user's safeguarded Newton (values = 2): Newton's step from the end with the smaller |f| where it lands
    // tol inside the bracket, else the midpoint. Reached through the separate-derivative form of solve().
    struct SafeNewton
    {
        static constexpr int values = 2;

        struct State
        {
            double a;
            std::array<double, 2> va;
            double b;
            std::array<double, 2> vb;
            rt::Step<double> step;
        };

        static constexpr rt::Step<double> report(State const &s)
        {
            return rt::Step<double>{ s.a, s.va[0], s.b, s.vb[0],
                nxx::math::abs(s.va[0]) < nxx::math::abs(s.vb[0]) ? s.a : s.b };
        }

        template<typename F>
        constexpr std::expected<State, rt::error> start(F &, rt::Start<double, 2> const &s, rt::Tolerance<double> const &) const
        {
            State st{ s.lower, s.at_lower, s.upper, s.at_upper, {} };
            st.step = report(st);
            return st;
        }

        template<typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State &s, rt::Tolerance<double> const &t) const
        {
            bool const at_a = s.step.point == s.a;
            std::array<double, 2> const &v = at_a ? s.va : s.vb;
            double const from = at_a ? s.a : s.b;
            double const margin = t.at(from);
            double x = s.a / 2 + s.b / 2;
            if (v[1] != 0.0) {
                double const newton = from - v[0] / v[1];
                if (s.a + margin < newton && newton < s.b - margin) x = newton;
            }
            if (!(s.a < x && x < s.b)) return std::unexpected(rt::error::not_converged);
            auto const r = std::invoke(fn, x);
            std::array<double, 2> const vx{ std::get<0>(r), std::get<1>(r) };
            if (!nxx::math::isfinite(vx[0])) return std::unexpected(rt::error::non_finite);
            if (vx[0] == 0.0) {
                s.step = rt::Step<double>{ s.a, s.va[0], x, vx[0], x };
                return {};
            }
            if ((vx[0] < 0.0) == (s.va[0] < 0.0)) {
                s.a = x;
                s.va = vx;
            } else {
                s.b = x;
                s.vb = vx;
            }
            s.step = report(s);
            return {};
        }
    };

    static_assert(rt::bracketing_method<Pegasus, double> && rt::bracketing_method<Pegasus, bin50>);
    static_assert(rt::bracketing_method<Section, double> && !std::default_initializable<Section>);
    static_assert(rt::bracketing_method<SafeNewton, double> && rt::values_of<SafeNewton> == 2);
    // The named form needs a default-constructible method.
    template<typename M>
    concept named_on_bracket = requires { rt::solve<M>(g, one_two); };
    static_assert(named_on_bracket<Pegasus> && !named_on_bracket<Section>);
    static_assert(rt::solve(Section{ 0.25 }, g, one_two).has_value());

    // ---- The concepts --------------------------------------------------------------------------------------------

    static_assert(rt::bracketing_method<rt::Bisection, double> && rt::bracketing_method<rt::Brent, double>
                  && rt::bracketing_method<rt::Toms748, double> && rt::bracketing_method<rt::AndersonBjorck, double>
                  && rt::bracketing_method<rt::Ridders, double> && rt::bracketing_method<rt::RegulaFalsi, double>
                  && rt::bracketing_method<rt::Illinois, double> && rt::bracketing_method<rt::Newton, double>
                  && rt::bracketing_method<rt::Halley, double>);
    static_assert(rt::bracketing_method<rt::Brent, bin50> && rt::bracketing_method<rt::Halley, dec50>
                  && rt::bracketing_method<rt::Toms748, float>);
    static_assert(!rt::bracketing_method<rt::Secant, double> && !rt::bracketing_method<rt::Steffensen, double>);
    static_assert(!rt::bracketing_method<int, double>);
    struct Empty
    {
    };
    static_assert(!rt::bracketing_method<Empty, double>);
    static_assert(rt::values_of<rt::Newton> == 2 && rt::values_of<rt::Halley> == 3 && rt::values_of<rt::Brent> == 1
                  && rt::values_of<Empty> == 1);
    // Only a static integral constant counts (numerics review of draft 1, finding 7).
    struct WithTable
    {
        std::array<double, 3> values;
    };
    struct WithStaticTable
    {
        static constexpr std::array<double, 3> values{};
    };
    struct WithCount
    {
        int values;
    };
    static_assert(rt::values_of<WithTable> == 1 && rt::values_of<WithStaticTable> == 1 && rt::values_of<WithCount> == 1);
    // At least 1: fn's value is always read (review of step 1, finding 4).
    struct NoValues
    {
        static constexpr int values = 0;
    };
    struct NegativeValues
    {
        static constexpr int values = -2;
    };
    static_assert(rt::values_of<NoValues> == 1 && rt::values_of<NegativeValues> == 1);

    // Brent's members around a state taken another way, to show the concept rejects exactly that.
    template<typename Pass>
    struct Reshaped
    {
        template<typename F>
        constexpr std::expected<rt::bracketing_state_t<rt::Brent, double>, rt::error>
            start(F &fn, rt::Start<double> const &s, rt::Tolerance<double> const &t) const
        { return rt::Brent{}.start(fn, s, t); }

        template<typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, Pass s, rt::Tolerance<double> const &t) const
        {
            auto copy = s;
            return rt::Brent{}.advance(fn, copy, t);
        }
    };
    using BrentState = rt::bracketing_state_t<rt::Brent, double>;
    static_assert(rt::bracketing_method<Reshaped<BrentState &>, double>);
    // A state taken by value or by const reference would be stepped as a copy, for ever (usage review of draft 1,
    // finding 2).
    static_assert(!rt::bracketing_method<Reshaped<BrentState>, double>);
    static_assert(!rt::bracketing_method<Reshaped<BrentState const &>, double>);
    // A non-const advance cannot be called on the const method solve() holds.
    struct NonConstAdvance
    {
        template<typename F>
        constexpr std::expected<BrentState, rt::error>
            start(F &fn, rt::Start<double> const &s, rt::Tolerance<double> const &t) const
        { return rt::Brent{}.start(fn, s, t); }

        template<typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, BrentState &s, rt::Tolerance<double> const &t)
        { return rt::Brent{}.advance(fn, s, t); }
    };
    static_assert(!rt::bracketing_method<NonConstAdvance, double>);

    // Today's calls with a generic lambda as fn resolve, named, as a value and observed.
    constexpr auto generic = [](auto const &x) { return x * x - 2.0; };
    static_assert(rt::solve<rt::Brent>(generic, one_two).has_value());
    static_assert(rt::solve<rt::Brent>(generic, one_two, { .rel_tol = 1e-8 }).has_value());
    static_assert(rt::solve(rt::Brent{}, generic, one_two).has_value());
    static_assert(rt::solve<rt::Secant>(generic, 1.0).has_value());
    constexpr bool generic_observed()
    {
        auto ignore = [](auto const &) {};
        return rt::solve<rt::Toms748>(rt::observer(ignore), generic, one_two)
               == rt::solve<rt::Toms748>(generic, one_two);
    }
    static_assert(generic_observed());

    // ---- Bisection's stop test (PHASE7_API.md B, Q3) ---------------------------------------------------------------

    // x - 2.4 on [1, 3] at rel_tol 0.5: the start's midpoint 2 is within tol(2) = 1 of both ends, so 2 with
    // abs_error 1, and no evaluation inside (numerics review of draft 1, finding 2: the width alone, raised by a unit,
    // would have gone on).
    constexpr auto shifted = [](double x) { return x - 2.4; };
    // 1/x on [-1, 1] at abs_tol 1 (at run time only: zero first evaluates 1/0): fn(0) is infinite, so zero first is
    // skipped; the midpoint 0 is within tol(0) = 1 of both ends, though the width 2, raised by a unit, exceeds 2 tol,
    // so the pole at 0 is returned with abs_error 1, as documented.
    constexpr auto pole_at_0 = [](double x) { return 1.0 / x; };
    // An odd-unit bracket [1, 1 + 3 eps]: 1.0 + 1.5 * deps rounds to even, 1 + 2 eps, so the sign change lies between
    // 1 + eps and 1 + 2 eps. 3 eps <= 2 tol(point) at the start, and the midpoint 1 + 1.5 eps rounds to even too,
    // 1 + 2 eps, which is 2 eps from the lower end.
    constexpr auto step_at = [](double x) { return x < 1.0 + 1.5 * deps ? -1.0 : 1.0; };
    // A saturated tolerance: x - 1 on [-max, max] at abs_tol max / 4, rel_tol 1. Zero first gives [0, max] (fn(0)
    // = -1); its midpoint max / 2 has tol = max / 4 + max / 2, and 2 tol saturates at max, so the start meets the
    // width test: max / 2 with abs_error max / 2, its distance to either end.
    constexpr auto minus_one = [](double x) { return x - 1.0; };

    // The cases above, true where Bisection returns what the comments derive.
    constexpr std::array<bool, 3> bisection_cases()
    {
        return {
            rt::solve<rt::Bisection>(shifted, rt::Bracket{ 1.0, 3.0 }, { .rel_tol = 0.5 }) == rt::Root<double>{ 2.0, 1.0 },
            rt::solve<rt::Bisection>(step_at, rt::Bracket{ 1.0, 1.0 + 3.0 * deps }) == rt::Root<double>{ 1.0 + 2.0 * deps, 2.0 * deps },
            rt::solve<rt::Bisection>(minus_one, rt::Bracket{ -dmax, dmax }, { .abs_tol = dmax / 4, .rel_tol = 1.0 }) == rt::Root<double>{ dmax / 2, dmax / 2 },
        };
    }

    static_assert(std::ranges::all_of(bisection_cases(), [](bool b) { return b; }));

    // ---- distance_bound is symmetric (PHASE7_API.md B relies on it) --------------------------------------------------

    template<typename T>
    constexpr bool symmetric(T const &u, T const &v)
    { return rt::detail::distance_bound(u, v) == rt::detail::distance_bound(v, u); }

    template<typename T>
    constexpr bool symmetric_cases()
    {
        T const max = (std::numeric_limits<T>::max)();
        T const tiny = std::numeric_limits<T>::denorm_min();
        T const least = (std::numeric_limits<T>::min)();
        return symmetric(T(1), T(1.5)) && symmetric(T(1), T(3)) && symmetric(T(-2), T(-0.75))
               && symmetric(T(-1), T(2)) && symmetric(T(max / T(-3)), T(max / T(2.5)))
               && symmetric(T(-max), T(max)) && symmetric(T(0), T(3)) && symmetric(T(-0.0), T(-5))
               && symmetric(tiny, T(T(3) * tiny)) && symmetric(least, T(T(2.5) * least))
               && symmetric(T(-least), T(T(1.5) * least));
    }

    static_assert(symmetric_cases<float>() && symmetric_cases<double>() && symmetric_cases<long double>());

    // ---- The fix for types without subnormals (PHASE7_API.md F) -----------------------------------------------------

    // [min, 2.5 min] around a root r = c min, default criteria, in a type whose denorm_min is min: before the fix every
    // method returned min with abs_error 2.5 min > 2 tol(min) = 2 min. Each result must now bracket r within
    // abs_error, and abs_error be at most 2 tol(value), with tol(x) = min + 2 eps |x| (2 eps |x| flushes to 0 here).
    template<typename T, typename M>
    bool tiny_bracket_root(T const &c)
    {
        T const least = (std::numeric_limits<T>::min)();
        T const r = T(c * least);
        rt::Bracket<T> const b{ least, T(T(2.5) * least) };
        std::expected<rt::Root<T>, rt::error> result;
        if constexpr (rt::values_of<M> == 1) {
            result = rt::solve<M>([&](T const &x) { return T(T(x / least) - c); }, b);
        } else if constexpr (rt::values_of<M> == 2) {
            result = rt::solve<M>([&](T const &x) { return std::array<T, 2>{ T(T(x / least) - c), T(T(1) / least) }; }, b);
        } else {
            result = rt::solve<M>(
                [&](T const &x) { return std::array<T, 3>{ T(T(x / least) - c), T(T(1) / least), T(0) }; }, b);
        }
        if (!result) return false;
        T const tol = T(least + T(T(2) * std::numeric_limits<T>::epsilon() * abs(result->value)));
        return abs(T(result->value - r)) <= result->abs_error && result->abs_error <= T(T(2) * tol);
    }

    template<typename T>
    bool tiny_brackets_all()
    {
        bool ok = true;
        for (T const &c : { T(1.25), T(2) }) {
            ok = ok && tiny_bracket_root<T, rt::Bisection>(c) && tiny_bracket_root<T, rt::Brent>(c)
                 && tiny_bracket_root<T, rt::Toms748>(c) && tiny_bracket_root<T, rt::AndersonBjorck>(c)
                 && tiny_bracket_root<T, rt::Ridders>(c) && tiny_bracket_root<T, rt::RegulaFalsi>(c)
                 && tiny_bracket_root<T, rt::Illinois>(c) && tiny_bracket_root<T, rt::Newton>(c)
                 && tiny_bracket_root<T, rt::Halley>(c);
        }
        return ok;
    }

} // namespace

TEST_CASE("roots:protocol:loop")
{
    // The constant-evaluation cases again at run time.
    auto const cases = loop_cases();
    for (std::size_t i = 0; i < cases.size(); ++i) {
        INFO("loop case ", i);
        CHECK(cases[i]);
    }
    CHECK(std::format("{}", rt::error::invalid_step) == "invalid_step");
}

TEST_CASE("roots:protocol:observing")
{
    CHECK(same_everywhere<rt::Bisection>(g));
    CHECK(same_everywhere<rt::Toms748>(g));
    CHECK(same_everywhere<rt::Halley>(g2d));
    CHECK(one_report_per_evaluation<rt::Ridders>(g));

    // An lvalue recorder sees every step, in order: the start's bracket first, then nested brackets, the last the
    // result's.
    std::vector<rt::Step<double>> steps;
    auto record = [&steps](rt::Step<double> const &s) { steps.push_back(s); };
    auto const r = rt::solve<rt::Brent>(rt::observer(record), g, one_two);
    REQUIRE(r.has_value());
    REQUIRE(steps.size() >= 2);
    CHECK(steps.front().lower == 1.0);
    CHECK(steps.front().upper == 2.0);
    CHECK(steps.back().point == r->value);
    for (std::size_t i = 1; i < steps.size(); ++i) {
        CHECK(steps[i - 1].lower <= steps[i].lower);
        CHECK(steps[i].upper <= steps[i - 1].upper);
    }

    // A callable with a non-const call operator: an lvalue is kept by reference, so its count is the caller's; an
    // rvalue is kept by value, and the run is the same.
    Tally tally;
    auto const s = rt::solve<rt::Bisection>(rt::observer(tally), g, one_two);
    REQUIRE(s.has_value());
    CHECK(tally.n == static_cast<int>(steps_of<rt::Bisection>()));
    CHECK(rt::solve<rt::Bisection>(rt::observer(Tally{}), g, one_two) == s);
    // A mutable lambda as an observer.
    int calls = 0;
    auto mutable_counter = [&calls, local = 0](rt::Step<double> const &) mutable { calls = ++local; };
    CHECK(rt::solve<rt::Bisection>(rt::observer(mutable_counter), g, one_two) == s);
    CHECK(calls == tally.n);
}

TEST_CASE("roots:protocol:methods-of-ones-own")
{
    // Pegasus on x^3 - 2x - 5 over [2, 3] (Wallis's root) and x^2 - 2 over [1, 2], double and cpp_bin_float_50. The
    // conditioning terms: for the cubic, fn's rounding near the root (about 10 eps of its terms of up to 9) over
    // |f'(root)| = 11.2, taken as 30 eps / 11.2; for x^2 - 2, x * x rounds by at most eps near sqrt 2 and subtracting
    // 2 is exact (Sterbenz), eps / 2.83, plus half a unit from rounding the reference: 0.85 eps, taken as eps (in
    // cpp_bin_float_50 as 2 eps, the reference parsed from 70 digits).
    CHECK(brackets(rt::solve<Pegasus>(wal_fn, rt::Bracket{ 2.0, 3.0 }), wallis, 30.0 * deps / 11.2));
    CHECK(brackets(rt::solve<Pegasus>(g, one_two), sqrt2, deps));
    bin50 const bin_sqrt2 = from_string<bin50>(sqrt2_70);
    auto const big = rt::solve<Pegasus>([](bin50 const &x) { return bin50(x * x - 2); }, rt::Bracket<bin50>{ 1, 2 });
    REQUIRE(big.has_value());
    bin50 const bin_eps = std::numeric_limits<bin50>::epsilon();
    bin50 const bin_tol = bin50(std::numeric_limits<bin50>::denorm_min() + bin50(2) * bin_eps * abs(big->value));
    CHECK(brackets(big, bin_sqrt2, bin50(bin50(2) * bin_eps), bin50(bin50(2) * bin_tol)));

    // A section method given as a value, with its parameter, against the same references.
    CHECK(brackets(rt::solve(Section{ 0.25 }, g, one_two), sqrt2, deps));
    CHECK(brackets(rt::solve(Section{ 0.1 }, wal_fn, rt::Bracket{ 2.0, 3.0 }), wallis, 30.0 * deps / 11.2));

    // A user's method with values = 2 through the separate-derivative form, and with fn returning both values.
    auto const df = [](double x) { return 2.0 * x; };
    CHECK(brackets(rt::solve<SafeNewton>(g, df, one_two), sqrt2, deps));
    CHECK(brackets(rt::solve<SafeNewton>(g1d, one_two), sqrt2, deps));
}

TEST_CASE("roots:protocol:bisection")
{
    auto const cases = bisection_cases();
    for (std::size_t i = 0; i < cases.size(); ++i) {
        INFO("bisection case ", i);
        CHECK(cases[i]);
    }
    CHECK(rt::solve<rt::Bisection>(pole_at_0, rt::Bracket{ -1.0, 1.0 }, { .abs_tol = 1.0 })
          == rt::Root<double>{ 0.0, 1.0 });

    // The battery's changed results (PHASE7_API.md, Spike): each within 2 tol of a sign change, one evaluation
    // earlier than before this phase. x^2 - 2 on [0, 10]: after k halvings the bracket is 10 2^-k wide, and at the
    // stop 5 units of T's spacing at sqrt 2, so its ends are exact. With 2 tol(sqrt 2) about 5.66 eps, the width
    // test first passes at k = p + 1 for a significand of p bits (10 2^-(p+1) = 5 2^-p < 5.66 eps = 5.66 2^(1-p) / 2,
    // while at k = p the width 10 2^-p fails, and so does the half width 5 2^-p > tol): 2 + 24 = 26 evaluations in
    // float, 2 + 53 = 55 in double, 2 + 64 = 66 in x87's long double.
    // The conditioning term: near sqrt 2, x * x rounds by at most half a unit of 2, eps, and subtracting 2 is exact
    // (Sterbenz), so the sign change moves by at most eps / 2.83; the double literal sqrt2 adds half a unit, 0.5 eps.
    int n_float = 0;
    auto const sq_float = [&n_float](float x) { ++n_float; return x * x - 2.0f; };
    auto const rf = rt::solve<rt::Bisection>(sq_float, rt::Bracket{ 0.0f, 10.0f });
    REQUIRE(rf.has_value());
    float const feps = std::numeric_limits<float>::epsilon();
    CHECK(nxx::math::abs(rf->value - static_cast<float>(sqrt2)) <= rf->abs_error + feps);
    CHECK(rf->abs_error <= 2.0f * (std::numeric_limits<float>::denorm_min() + 2.0f * feps * rf->value));
    CHECK(n_float == 26);

    auto const [rd, n_double] = counted<rt::Bisection>(g, rt::Bracket{ 0.0, 10.0 });
    CHECK(brackets(rd, sqrt2, deps));
    CHECK(n_double == 55);

    int n_long = 0;
    auto const sq_long = [&n_long](long double x) { ++n_long; return x * x - 2.0L; };
    auto const rl = rt::solve<rt::Bisection>(sq_long, rt::Bracket{ 0.0L, 10.0L });
    REQUIRE(rl.has_value());
    long double const leps = std::numeric_limits<long double>::epsilon();
    long double const lsqrt2 = from_string<long double>(sqrt2_70);
    CHECK(nxx::math::abs(rl->value - lsqrt2) <= rl->abs_error + leps);
    CHECK(rl->abs_error <= 2.0L * (std::numeric_limits<long double>::denorm_min() + 2.0L * leps * rl->value));

    // The count's derivation holds for 53 and 64 bits (MSVC's long double is double); in binary128 (Emscripten) the
    // stop comes at a width of 2.5 units, where midpoints round. The e^x case is a regression baseline recorded with
    // x87's 64-bit significand: its ends round and expl's last bits decide the path.
    if constexpr (std::numeric_limits<long double>::digits == 53) CHECK(n_long == 55);
    if constexpr (std::numeric_limits<long double>::digits == 64) {
        CHECK(n_long == 66);

        // e^x - 2 on [-700, 700]: ln 2 within abs_error plus e^x's 2 ulps over f' = 2. Never later than the 75
        // evaluations of the header at 2dd3228 (regression baseline).
        int n_exp = 0;
        auto const exp2 = [&n_exp](long double x) { ++n_exp; return std::exp(x) - 2.0L; };
        auto const re = rt::solve<rt::Bisection>(exp2, rt::Bracket{ -700.0L, 700.0L });
        REQUIRE(re.has_value());
        long double const lln2 = 0.693147180559945309417232121458176568L;
        CHECK(nxx::math::abs(re->value - lln2) <= re->abs_error + 2.0L * leps);
        CHECK(re->abs_error <= 2.0L * (std::numeric_limits<long double>::denorm_min() + 2.0L * leps * re->value));
        CHECK(n_exp <= 75);
    }
}

TEST_CASE("roots:protocol:tiny-brackets")
{
    CHECK(tiny_brackets_all<bin50>());
    CHECK(tiny_brackets_all<bin50_et>());
    CHECK(tiny_brackets_all<dec50>());
}

TEST_CASE("roots:protocol:distance-bound")
{
    CHECK(symmetric_cases<bin50>());
    CHECK(symmetric_cases<dec50>());
}

TEST_CASE("roots:protocol:negative-zero")
{
    // Ends at -0.0: not around 0, so no zero first; the validated bracket starts at -0.0 and every step lies inside
    // it, 0.0 included (0.0 < -0.0 is false). x - 0.25 is exact near 0.25 (Sterbenz), so the sign change is at 0.25: c = 0.
    auto const quarter = [](double x) { return x - 0.25; };
    auto const neg_quarter = [](double x) { return x + 0.25; };
    CHECK(brackets(rt::solve<rt::Brent>(quarter, rt::Bracket{ -0.0, 1.0 }), 0.25, 0.0));
    CHECK(brackets(rt::solve<rt::Bisection>(quarter, rt::Bracket{ -0.0, 1.0 }), 0.25, 0.0));
    CHECK(brackets(rt::solve<rt::Toms748>(quarter, rt::Bracket{ 1.0, -0.0 }), 0.25, 0.0));
    CHECK(brackets(rt::solve<rt::AndersonBjorck>(neg_quarter, rt::Bracket{ -1.0, -0.0 }), -0.25, 0.0));
    CHECK(brackets(rt::solve<rt::Ridders>(neg_quarter, rt::Bracket{ -0.0, -1.0 }), -0.25, 0.0));
    // fn exactly 0 at -0.0: the end is the root.
    auto const identity = [](double x) { return x; };
    auto const at_zero = rt::solve<rt::Illinois>(identity, rt::Bracket{ -0.0, 1.0 });
    REQUIRE(at_zero.has_value());
    CHECK(at_zero->value == 0.0);
}
