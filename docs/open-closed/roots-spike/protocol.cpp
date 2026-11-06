// The loop of PHASE7_API.md B, branch by branch, with a scripted method whose steps come from a list; observing
// changes nothing; the concepts accept and reject what they should. Everything here is checked in constant evaluation
// (static_assert) and again at run time.
#include <numerixx/roots.hpp>

#include <array>
#include <cstdio>
#include <limits>

namespace rt = nxx::roots;

// ---- A scripted method: its steps come from a list; an `advance` past the list's end is an error. ----------------

template<std::size_t N>
struct Script
{
    std::array<rt::Step<double>, N> steps;
    int fail_at = -1; // the step whose advance returns not_converged, if any

    struct State
    {
        std::array<rt::Step<double>, N> steps;
        std::size_t next;
        int fail_at;
        int advances;
        rt::Step<double> step;
    };

    template<typename F>
    constexpr std::expected<State, rt::error> start(F &, rt::Start<double> const &, rt::Tolerance<double> const &) const
    { return State{ steps, 1, fail_at, 0, steps[0] }; }

    template<typename F>
    constexpr std::expected<State, rt::error> start(F &, double const &, rt::Tolerance<double> const &) const
    { return State{ steps, 1, fail_at, 0, steps[0] }; }

    template<typename F>
    constexpr std::expected<void, rt::error> advance(F &, State &s, rt::Tolerance<double> const &) const
    {
        if (static_cast<int>(s.next) == s.fail_at) return std::unexpected(rt::error::not_converged);
        if (s.next >= N) return std::unexpected(rt::error::no_sign_change); // a script too short: visible in a test
        s.step = s.steps[s.next++];
        ++s.advances;
        return {};
    }
};

constexpr double f(double const &x) { return x - 0.5; } // fn: [0, 1] is a valid bracket around 0.5
constexpr rt::Bracket<double> unit{ 0.0, 1.0 };
constexpr double inf = std::numeric_limits<double>::infinity();

template<std::size_t N>
constexpr std::expected<rt::Root<double>, rt::error> on_unit(Script<N> const &m)
{ return rt::solve(m, f, unit, { .abs_tol = 1e-3 }); }

template<std::size_t N>
constexpr std::expected<rt::Root<double>, rt::error> from(Script<N> const &m)
{ return rt::solve(m, f, 0.25, { .abs_tol = 1e-3 }); }

// The bracket [0, 1], then [0.25, 0.75], then [0.4995, 0.5005], which meets 2 tol = 2e-3 + ...: the third step stops.
constexpr rt::Step<double> s0{ 0.0, -0.5, 1.0, 0.5, 0.5 };
constexpr rt::Step<double> s1{ 0.25, -0.25, 0.75, 0.25, 0.75 };
constexpr rt::Step<double> s2{ 0.4995, -0.0005, 0.5005, 0.0005, 0.5005 };

// Stops at the first step that meets the request, with abs_error the larger distance to an end (here the width).
static_assert(on_unit(Script<3>{ { s0, s1, s2 } })->value == 0.5005);
static_assert(on_unit(Script<3>{ { s0, s1, s2 } })->abs_error > 0.0009); // the width, about 1e-3
// An interior point: abs_error is the larger of the two distances.
static_assert(on_unit(Script<1>{ { rt::Step<double>{ 0.4995, -0.0005, 0.5005, 0.0005, 0.5 } } })->abs_error > 0.00049); // half the width
// A point that is an end where fn is exactly 0: that point, abs_error tol(point).
static_assert(on_unit(Script<2>{ { s0, rt::Step<double>{ 0.25, -0.25, 0.5, 0.0, 0.5 } } })->value == 0.5);
// A NaN or an infinity anywhere in a step: non_finite.
static_assert(on_unit(Script<2>{ { s0, rt::Step<double>{ 0.25, -0.25, inf, 0.25, 0.5 } } }).error()
              == rt::error::non_finite);
// The point outside its points; a step outside the caller's bracket; a lost sign change: invalid_step.
static_assert(on_unit(Script<2>{ { s0, rt::Step<double>{ 0.25, -0.25, 0.75, 0.25, 0.8 } } }).error()
              == rt::error::invalid_step);
static_assert(on_unit(Script<2>{ { s0, rt::Step<double>{ 0.25, -0.25, 1.5, 0.25, 0.75 } } }).error()
              == rt::error::invalid_step);
static_assert(on_unit(Script<2>{ { s0, rt::Step<double>{ 0.25, 0.25, 0.75, 0.25, 0.75 } } }).error()
              == rt::error::invalid_step);
// Opposite signs at a single point: invalid_step.
static_assert(on_unit(Script<2>{ { s0, rt::Step<double>{ 0.5, -0.25, 0.5, 0.25, 0.5 } } }).error()
              == rt::error::invalid_step);
// fn exactly 0 at an end that is not the point: no sign change shown, so invalid_step on a bracket.
static_assert(on_unit(Script<2>{ { s0, rt::Step<double>{ 0.25, -0.25, 0.5, 0.0, 0.25 } } }).error()
              == rt::error::invalid_step);
// Progress: the same bracket again, or one not nested in the last, is invalid_step (numerics review, finding 3).
static_assert(on_unit(Script<3>{ { s0, s1, s1 } }).error() == rt::error::invalid_step);
static_assert(on_unit(Script<3>{ { s0, s1, rt::Step<double>{ 0.2, -0.3, 0.6, 0.1, 0.6 } } }).error()
              == rt::error::invalid_step);
// One end kept, the other moved inside: nested.
static_assert(on_unit(Script<3>{ { s0, rt::Step<double>{ 0.0, -0.5, 0.75, 0.25, 0.75 }, s2 } })->value == 0.5005);
// A method that cannot go on returns not_converged from advance, passed through as it is (any error is).
static_assert(on_unit(Script<3>{ { s0, s1, s2 }, 1 }).error() == rt::error::not_converged);

// From a guess: unbracketed steps (same signs) count against the cap; bracketed ones do not; once bracketed, an
// unbracketed step is invalid.
constexpr rt::Step<double> g0{ 0.25, -0.25, 0.25, -0.25, 0.25 };
constexpr rt::Step<double> g1{ 0.25, -0.25, 0.3, -0.2, 0.3 };
static_assert(from(Script<3>{ { g0, g1, s2 } })->value == 0.5005);
// A sign change found and lost again (simplicity review, finding 3): allowed; the step without one counts to the cap.
static_assert(from(Script<5>{ { g0, g1, s1, g1, s2 } })->value == 0.5005);
// From a guess, a bracketed step after a bracketed one must nest too.
static_assert(from(Script<4>{ { g0, s1, s1, s2 } }).error() == rt::error::invalid_step);
// The cap: 4 times the bits of double, 212 unbracketed steps after the start, then not_converged.
constexpr auto capped = [] {
    Script<214> m{};
    for (auto &s : m.steps) s = g0;
    m.steps[213] = s2;
    return m;
}();
static_assert(from(capped).error() == rt::error::not_converged);
constexpr auto under_cap = [] {
    Script<214> m{};
    for (auto &s : m.steps) s = g0;
    m.steps[212] = s2;
    return m;
}();
static_assert(from(under_cap)->value == 0.5005);
constexpr auto well_under_cap = [] {
    Script<214> m{};
    for (auto &s : m.steps) s = g0;
    m.steps[211] = s2;
    return m;
}();
static_assert(from(well_under_cap)->value == 0.5005);

// ---- Observing changes nothing, for every library method, in both forms. ------------------------------------------

constexpr double g(double const &x) { return x * x - 2.0; }
constexpr std::array<double, 2> g1d(double const &x) { return { x * x - 2.0, 2.0 * x }; }
constexpr std::array<double, 3> g2d(double const &x) { return { x * x - 2.0, 2.0 * x, 2.0 }; }

template<typename M, typename F, typename S>
constexpr bool same_observed(F const &fn, S const &start)
{
    int reports = 0;
    auto count = [&reports](rt::Step<double> const &) { ++reports; };
    auto const plain = rt::solve<M>(fn, start);
    auto const named = rt::solve<M>(rt::observer(count), fn, start);
    auto const value = rt::solve(M{}, rt::observer(count), fn, start);
    return plain == named && plain == value && reports > 0;
}

static_assert(same_observed<rt::Bisection>(g, rt::Bracket{ 1.0, 2.0 }));
static_assert(same_observed<rt::Brent>(g, rt::Bracket{ 1.0, 2.0 }));
static_assert(same_observed<rt::Toms748>(g, rt::Bracket{ 1.0, 2.0 }));
static_assert(same_observed<rt::AndersonBjorck>(g, rt::Bracket{ 1.0, 2.0 }));
static_assert(same_observed<rt::Ridders>(g, rt::Bracket{ 1.0, 2.0 }));
static_assert(same_observed<rt::RegulaFalsi>(g, rt::Bracket{ 1.0, 2.0 }));
static_assert(same_observed<rt::Illinois>(g, rt::Bracket{ 1.0, 2.0 }));
static_assert(same_observed<rt::Newton>(g1d, rt::Bracket{ 1.0, 2.0 }));
static_assert(same_observed<rt::Halley>(g2d, rt::Bracket{ 1.0, 2.0 }));
static_assert(same_observed<rt::Newton>(g1d, 1.0));
static_assert(same_observed<rt::Halley>(g2d, 1.0));
static_assert(same_observed<rt::Secant>(g, 1.0));
static_assert(same_observed<rt::Steffensen>(g, 1.0));

// Bisection reports once after its start and once per evaluation; its last report is the result's bracket.
constexpr int bisection_reports()
{
    int reports = 0;
    int evaluations = 0;
    auto counted = [&evaluations](double const &x) { ++evaluations; return x * x - 2.0; };
    auto count = [&reports](rt::Step<double> const &) { ++reports; };
    (void)rt::solve<rt::Bisection>(rt::observer(count), counted, rt::Bracket{ 1.0, 2.0 });
    return evaluations - 2 + 1 == reports ? reports : -1; // 2 end evaluations, 1 report after the start
}
static_assert(bisection_reports() > 0);

// ---- The concepts. ------------------------------------------------------------------------------------------------

static_assert(rt::bracketing_method<rt::Brent, double>);
static_assert(rt::bracketing_method<rt::Newton, double>);
static_assert(rt::open_method<rt::Newton, double>);
static_assert(rt::open_method<rt::Secant, double>);
static_assert(!rt::open_method<rt::Brent, double>);
static_assert(!rt::bracketing_method<rt::Secant, double>);
static_assert(!rt::bracketing_method<int, double>);
static_assert(rt::values_of<rt::Halley> == 3 && rt::values_of<rt::Brent> == 1);
struct WithTable { std::array<double, 3> values; };
struct WithStaticTable { static constexpr std::array<double, 3> values{}; };
static_assert(rt::values_of<WithTable> == 1 && rt::values_of<WithStaticTable> == 1); // numerics review, finding 7
// A state taken by value in advance (a forgotten &) is not a method (usage review, finding 2).
struct ByValue
{
    struct State { rt::Step<double> step; };
    template<typename F>
    constexpr std::expected<State, rt::error> start(F &, rt::Start<double> const &s, rt::Tolerance<double> const &) const
    { return State{ { s.lower, s.at_lower[0], s.upper, s.at_upper[0], s.lower } }; }
    template<typename F>
    constexpr std::expected<void, rt::error> advance(F &, State, rt::Tolerance<double> const &) const { return {}; }
};
static_assert(!rt::bracketing_method<ByValue, double>);
static_assert(rt::open_method<rt::Secant, double> && rt::open_method<rt::Steffensen, double>
              && rt::open_method<rt::Halley, double>);

// Today's calls with a generic lambda as fn still resolve.
constexpr auto generic = [](auto const &x) { return x * x - 2.0; };
static_assert(rt::solve<rt::Brent>(generic, rt::Bracket{ 1.0, 2.0 }).has_value());
static_assert(rt::solve<rt::Brent>(generic, rt::Bracket{ 1.0, 2.0 }, { .rel_tol = 1e-8 }).has_value());
static_assert(rt::solve(rt::Brent{}, generic, rt::Bracket{ 1.0, 2.0 }).has_value());
static_assert(rt::solve<rt::Secant>(generic, 1.0).has_value());

int main()
{
    // The same at run time, for one branch of each kind.
    bool ok = on_unit(Script<3>{ { s0, s1, s2 } })->value == 0.5005
              && on_unit(Script<2>{ { s0, rt::Step<double>{ 0.25, 0.25, 0.75, 0.25, 0.75 } } }).error()
                     == rt::error::invalid_step
              && on_unit(Script<3>{ { s0, s1, s1 } }).error() == rt::error::invalid_step
              && from(capped).error() == rt::error::not_converged && same_observed<rt::Toms748>(g, rt::Bracket{ 1.0, 2.0 })
              && same_observed<rt::Steffensen>(g, 1.0);
    std::printf("%s\n", ok ? "protocol: all checks passed" : "protocol: FAILED");
    return ok ? 0 : 1;
}
