// The method protocol of nxx::quadrature (Phase 3a, docs/quadrature/PHASE3A_API.md): methods as types, the loop that
// integrate() runs for an iterative method, the ends put in order, the checked results, and the observer.
//
// The methods of one's own here are the tests' own, written to the contract of the concepts' comments. References are
// closed forms, or mpmath at 80 digits for cos x / sqrt x; every tolerance is derived where it is used. The scripted
// method drives each branch of integrate()'s loop with steps from a list, so those tests check the loop's decisions,
// not a numerical method.

#include "quadrature_check.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

using namespace quadrature_check;

namespace {

    using gk15 = qd::GaussKronrod<15>;

    // --- Methods of one's own -------------------------------------------------------------------------------------

    // A fixed rule with a run-time parameter: the composite midpoint rule on `panels` panels of [a, b], a < b. fn is
    // evaluated only inside [a, b]: a width that overflows is refused, and a point that rounds past b is clamped.
    struct Midpoint
    {
        int panels = 100;

        template<nxx::real T, typename F>
        constexpr std::expected<T, qd::error> apply(F &fn, T const &a, T const &b) const
        {
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
    static_assert(qd::fixed_rule<Midpoint, double>);

    // A combinator over any estimating method: x = a + (b - a) u^2 at the lower end a (integrate() passes a < b), then
    // the inner method over u in [0, 1] through the public integrate(). It removes a singularity (x - a)^(-1/2) at a.
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
    static_assert(qd::estimating_method<SquareAtLowerEnd<gk15>, double>);

    // An iterative method: the trapezoidal rule, halving its step. Its own limits: ready from level 3, exhausted at
    // level 20, the attainable error a floor of 50 eps times its sum of |terms|. fn only in [a, b]: the midpoint as
    // a/2 + b/2, a width that overflows refused, a point a + i h clamped to b. It refines only while its next step
    // (b - a) / 2^k is exact, so that its sums are the trapezoidal sums of [a, b]: a step that rounds (in the subnormal
    // range, 1.5 denorm_min to 2) would sample past b and sum the wrong width.
    struct Trapezoid
    {
        template<nxx::real T>
        struct State
        {
            T a, b, h, sum, abs_sum, previous;
            int level;
            qd::Step<T> step;
        };

        // Whether level k's step (b - a) / 2^k is exact: in radix 2 it is until it is subnormal.
        template<nxx::real T>
        static constexpr bool exact_step(T const &width, int k)
        {
            T const n = T(1 << k);
            T const h = T(width / n);
            return h > T(0) && T(h * n) == width;
        }

        template<nxx::real T>
        static constexpr void report(State<T> &s)
        {
            T const value = T(s.h * s.sum);
            T const floor = T(T(50) * T(std::numeric_limits<T>::epsilon()) * T(s.h * s.abs_sum));
            T const change = nxx::math::abs(T(value - s.previous));
            s.step = qd::Step<T>{ .current = qd::Integral<T>{ value, change > floor ? change : floor },
                .attainable = floor,
                .ready = s.level >= 3,
                .exhausted = s.level >= 20 || !exact_step(T(s.b - s.a), s.level + 1) };
            s.previous = value;
        }

        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, qd::Failure<T>> start(F &fn, T const &a, T const &b) const
        {
            qd::Failure<T> const invalid{ qd::error::non_finite, std::nullopt };
            T const width = T(b - a);
            if (!nxx::math::isfinite(width)) return std::unexpected(invalid);
            if (!exact_step(width, 1)) return std::unexpected(qd::Failure<T>{ qd::error::not_converged, std::nullopt });
            T const fa = T(std::invoke(fn, a));
            T const fb = T(std::invoke(fn, b));
            T const mid = T(std::invoke(fn, T(T(a / T(2)) + T(b / T(2)))));
            State<T> s{ a, b, T(width / T(2)), T(T(T(fa + fb) / T(2)) + mid),
                T(T(T(nxx::math::abs(fa) + nxx::math::abs(fb)) / T(2)) + nxx::math::abs(mid)),
                T(T(width / T(2)) * T(fa + fb)), 1, {} };
            report(s);
            if (!nxx::math::isfinite(s.step.current->value)) return std::unexpected(invalid);
            return s;
        }

        template<nxx::real T, typename F>
        constexpr std::expected<void, qd::Failure<T>> advance(F &fn, State<T> &s) const
        {
            ++s.level;
            int const n = 1 << s.level;
            s.h = T(T(s.b - s.a) / T(n));
            for (int i = 1; i < n; i += 2) {
                T x = T(s.a + T(T(i) * s.h));
                if (s.b < x) x = s.b;
                T const y = T(std::invoke(fn, x));
                s.sum = T(s.sum + y);
                s.abs_sum = T(s.abs_sum + nxx::math::abs(y));
            }
            report(s);
            if (!nxx::math::isfinite(s.step.current->value))
                return std::unexpected(qd::Failure<T>{ qd::error::non_finite, std::nullopt });
            return {};
        }

        template<nxx::real T>
        constexpr std::expected<void, qd::Failure<T>> settle(State<T> &) const
        { return {}; }
    };
    static_assert(qd::iterative_method<Trapezoid, double>);

    // An iterative method that cannot be copied: Romberg's steps, through state_t. It can be observed, since the
    // observed method holds it by reference.
    struct MoveOnly
    {
        constexpr MoveOnly() = default;
        MoveOnly(MoveOnly const &) = delete;
        constexpr MoveOnly(MoveOnly &&) = default;

        template<nxx::real T, typename F>
        constexpr std::expected<qd::state_t<qd::Romberg, T, F>, qd::Failure<T>> start(F &fn, T const &a, T const &b) const
        { return qd::Romberg{}.start(fn, a, b); }

        template<typename F, typename S>
        constexpr auto advance(F &fn, S &s) const -> decltype(qd::Romberg{}.advance(fn, s))
        { return qd::Romberg{}.advance(fn, s); }

        template<typename S>
        constexpr auto settle(S &s) const -> decltype(qd::Romberg{}.settle(s))
        { return qd::Romberg{}.settle(s); }
    };
    static_assert(qd::iterative_method<MoveOnly, double> && !std::copy_constructible<MoveOnly>);

    // Methods that take fn only as a plain function of double: as a pointer, which integrate() passes for a function,
    // or as a reference, which it never passes. Each returns (b - a) fn(a) at its start.
    struct OneStep
    {
        qd::Step<double> step;
    };

    constexpr std::expected<OneStep, qd::Failure<double>> one_step(double fa, double const &a, double const &b)
    { return OneStep{ qd::Step<double>{ .current = qd::Integral<double>{ (b - a) * fa, 0.0 }, .exhausted = true } }; }

    struct PointerOnly
    {
        constexpr std::expected<OneStep, qd::Failure<double>> start(double (*fn)(double const &), double const &a,
            double const &b) const
        { return one_step(fn(a), a, b); }

        constexpr std::expected<void, qd::Failure<double>> advance(double (*)(double const &), OneStep &) const
        { return {}; }

        constexpr std::expected<void, qd::Failure<double>> settle(OneStep &) const
        { return {}; }
    };

    struct ReferenceOnly
    {
        constexpr std::expected<OneStep, qd::Failure<double>> start(double (&fn)(double const &), double const &a,
            double const &b) const
        { return one_step(fn(a), a, b); }

        constexpr std::expected<void, qd::Failure<double>> advance(double (&)(double const &), OneStep &) const
        { return {}; }

        constexpr std::expected<void, qd::Failure<double>> settle(OneStep &) const
        { return {}; }
    };

    // What an estimating method receives, recorded through a pointer (the member is const).
    struct Seen
    {
        double a = 0, b = 0, abs_tol = -1, rel_tol = -1;
        int calls = 0;
    };

    struct Spy
    {
        Seen *seen = nullptr;

        template<nxx::real T, typename F>
        constexpr std::expected<qd::Integral<T>, qd::Failure<T>> estimate(F &, T const &a, T const &b,
            qd::Criteria<T> const &c) const
        {
            *seen = Seen{ a, b, c.abs_tol, c.rel_tol, seen->calls + 1 };
            return qd::Integral<T>{ T(1), T(0) };
        }
    };

    // An estimating method that takes infinite ends (PHASE4_API.md A): the spy, opted in.
    struct UnboundedSpy
    {
        static constexpr bool infinite_ends = true;

        Seen *seen = nullptr;

        template<nxx::real T, typename F>
        constexpr std::expected<qd::Integral<T>, qd::Failure<T>> estimate(F &, T const &a, T const &b,
            qd::Criteria<T> const &c) const
        {
            *seen = Seen{ a, b, c.abs_tol, c.rel_tol, seen->calls + 1 };
            return qd::Integral<T>{ T(1), T(0) };
        }
    };

    // Methods that break the contract's checkable part: integrate() returns non_finite for each.
    template<typename R>
    struct Returning
    {
        R result;

        template<typename F>
        constexpr R estimate(F &, double const &, double const &, qd::Criteria<double> const &) const
        { return result; }
    };

    struct NaNRule
    {
        template<typename F>
        constexpr std::expected<double, qd::error> apply(F &, double const &, double const &) const
        { return std::numeric_limits<double>::quiet_NaN(); }
    };

    // Not methods.
    struct WrongReturn
    {
        template<nxx::real T, typename F>
        constexpr T apply(F &, T const &, T const &) const
        { return T(0); }
    };

    struct NotConst
    {
        template<nxx::real T, typename F>
        constexpr std::expected<T, qd::error> apply(F &, T const &, T const &)
        { return T(0); }
    };

    // advance or settle taking the state by value or by const reference: they would step or settle a copy, and the
    // loop would go on for ever (docs/roots/PHASE7_API.md, Q17).
    template<bool AdvanceByValue, bool SettleByConstReference>
    struct CopiedState
    {
        struct State
        {
            qd::Step<double> step;
        };

        template<typename F>
        constexpr std::expected<State, qd::Failure<double>> start(F &, double const &, double const &) const
        { return State{}; }

        template<typename F>
        requires AdvanceByValue
        constexpr std::expected<void, qd::Failure<double>> advance(F &, State) const
        { return {}; }

        template<typename F>
        requires(!AdvanceByValue)
        constexpr std::expected<void, qd::Failure<double>> advance(F &, State &) const
        { return {}; }

        constexpr std::expected<void, qd::Failure<double>> settle(State const &) const
        requires SettleByConstReference
        { return {}; }

        constexpr std::expected<void, qd::Failure<double>> settle(State &) const
        requires(!SettleByConstReference)
        { return {}; }
    };

    struct NoSettle
    {
        struct State
        {
            qd::Step<double> step;
        };

        template<typename F>
        constexpr std::expected<State, qd::Failure<double>> start(F &, double const &, double const &) const
        { return State{}; }

        template<typename F>
        constexpr std::expected<void, qd::Failure<double>> advance(F &, State &) const
        { return {}; }
    };

    // A scripted iterative method: its steps come from a list, so that every branch of integrate()'s loop can be
    // driven. At index settle_at, settle replaces the step by `settled`; the failures come where they are set. It never
    // calls fn.
    struct Scripted
    {
        std::array<qd::Step<double>, 4> steps;
        int settle_at = -1;
        qd::Step<double> settled{};
        std::optional<qd::Failure<double>> fail_start{};
        int fail_advance_at = -1;
        int fail_settle_at = -1;
        int *settles = nullptr;  // counted where given
        int *advances = nullptr; // counted where given

        struct State
        {
            int i;
            qd::Step<double> step;
        };

        static constexpr qd::Failure<double> scripted_failure{ qd::error::not_converged, qd::Integral<double>{ 7, 1 } };

        template<typename F>
        constexpr std::expected<State, qd::Failure<double>> start(F &, double const &, double const &) const
        {
            if (fail_start) return std::unexpected(*fail_start);
            return State{ 0, steps[0] };
        }

        template<typename F>
        constexpr std::expected<void, qd::Failure<double>> advance(F &, State &s) const
        {
            // The loop never advances an exhausted step: a call there would be reported as non_finite.
            if (s.step.exhausted) return std::unexpected(qd::Failure<double>{ qd::error::non_finite, std::nullopt });
            if (advances != nullptr) ++*advances;
            if (s.i == fail_advance_at) return std::unexpected(scripted_failure);
            ++s.i;
            s.step = steps[static_cast<std::size_t>(s.i)];
            return {};
        }

        constexpr std::expected<void, qd::Failure<double>> settle(State &s) const
        {
            if (settles != nullptr) ++*settles;
            if (s.i == fail_settle_at) return std::unexpected(scripted_failure);
            if (s.i == settle_at) s.step = settled;
            return {};
        }
    };
    static_assert(qd::iterative_method<Scripted, double>);

    constexpr double cubic(double x) { return x * x * x - x; }
    constexpr double one(double) { return 1.0; }
    constexpr double runge(double x) { return 1.0 / (1.0 + x * x); }
    constexpr double unit(double const &) { return 1.0; }
    constexpr auto quiet = [](qd::Step<double> const &) {};
    constexpr auto generic_square = [](auto const &x) { return x * x; };

    template<typename M>
    concept integrable = requires { qd::integrate<M>(square, 0.0, 1.0); };
    template<typename M>
    concept observable = requires { qd::integrate<M>(qd::observer(quiet), square, 0.0, 1.0); }
                         && requires { qd::integrate(M{}, qd::observer(quiet), square, 0.0, 1.0); };
    template<typename M>
    concept integrable_unit = requires { qd::integrate<M>(unit, 0.0, 1.0); };
    template<typename M>
    concept braced_method = requires { qd::integrate<M>(square, 0.0, 1.0, { 3 }); };

    // The scripted steps (meets rel_tol 1e-2 on 2, or not).
    constexpr qd::Integral<double> two{ 2.0, 1e-3 };
    constexpr qd::Integral<double> wide_two{ 2.0, 1.0 };
    constexpr qd::Step<double> pass{ .current = two };
    constexpr qd::Step<double> miss{ .current = wide_two };
    constexpr qd::Criteria<double> loose{ .rel_tol = 1e-2 };
    constexpr qd::Failure<double> non_finite{ qd::error::non_finite, std::nullopt };

    template<typename R>
    constexpr R negated(R r)
    {
        auto const flip = [](qd::Integral<double> i) {
            i.value = -i.value;
            return i;
        };
        if (r) return flip(*r);
        qd::Failure<double> f = r.error();
        if (f.best) f.best = flip(*f.best);
        return std::unexpected(f);
    }

} // namespace

// --- The concepts ----------------------------------------------------------------------------------------------------

static_assert(qd::fixed_rule<qd::GaussLegendre<7>, double> && !qd::fixed_rule<qd::GaussLegendre<0>, double>);
static_assert(qd::iterative_method<qd::Romberg, double> && qd::iterative_method<gk15, double>);
static_assert(qd::iterative_method<qd::Romberg, bin50> && qd::iterative_method<gk15, float>);
static_assert(qd::iterative_method<qd::TanhSinh, double>);
// The kinds stay apart.
static_assert(!qd::fixed_rule<qd::Romberg, double> && !qd::estimating_method<qd::Romberg, double>);
static_assert(!qd::iterative_method<qd::GaussLegendre<7>, double> && !qd::estimating_method<gk15, double>);
// What is not a method: a wrong return type, a member that is not const, an iterative method without settle, int.
static_assert(!qd::fixed_rule<WrongReturn, double> && !integrable<WrongReturn>);
static_assert(!qd::fixed_rule<NotConst, double> && !integrable<NotConst>);
static_assert(!qd::iterative_method<NoSettle, double> && !integrable<NoSettle>);
// The state by non-const reference in advance and settle; the same type with both so is a method.
static_assert(qd::iterative_method<CopiedState<false, false>, double>);
static_assert(!qd::iterative_method<CopiedState<true, false>, double> && !integrable<CopiedState<true, false>>);
static_assert(!qd::iterative_method<CopiedState<false, true>, double> && !integrable<CopiedState<false, true>>);
// The observer is nxx::observer, shared with the other modules (numerixx/math.hpp).
static_assert(std::same_as<decltype(qd::observer(quiet)), decltype(nxx::observer(quiet))>);
static_assert(std::same_as<qd::Observer<int>, nxx::Observer<int>>);
static_assert(!integrable<int> && !integrable<qd::GaussLegendre<0>> && !integrable<qd::GaussLegendre<65>>);
// A braced list is never deduced as a method: with M named, a braced list after the ends is criteria, which a fixed rule
// does not take. (With M named, a braced list in the method's place is the method: integrate<M>({...}, fn, a, b).)
static_assert(!braced_method<Midpoint>);
// An observer only with an iterative method, in both forms.
static_assert(observable<qd::Romberg> && observable<gk15> && observable<qd::TanhSinh> && observable<Trapezoid>);
static_assert(!observable<qd::GaussLegendre<7>> && !observable<Midpoint>);
static_assert(!observable<SquareAtLowerEnd<gk15>>); // an estimating method has no steps to report
// A method that cannot be copied is integrated and observed, named and as a value: Romberg's computation.
static_assert(integrable<MoveOnly> && observable<MoveOnly>);
static_assert(qd::integrate<MoveOnly>(qd::observer(quiet), runge, 0.0, 1.0, { .rel_tol = 1e-10 })
              == qd::integrate<qd::Romberg>(runge, 0.0, 1.0, { .rel_tol = 1e-10 }));
static_assert(qd::integrate(MoveOnly{}, qd::observer(quiet), runge, 1.0, 0.0, { .rel_tol = 1e-10 })
              == qd::integrate<qd::Romberg>(runge, 1.0, 0.0, { .rel_tol = 1e-10 }));
// The concepts' default F is what integrate() passes for a plain function of T, a pointer to it, so the concepts and
// integrate() agree on a method that takes fn only as such: one that takes the pointer is both a method and
// accepted, and one that takes a reference to the function is neither. (The rule is exact: (2 - 0) 1.)
static_assert(qd::iterative_method<PointerOnly, double> && integrable_unit<PointerOnly>);
static_assert(qd::integrate<PointerOnly>(unit, 0.0, 2.0) == qd::Integral<double>{ 2.0, 0.0 });
static_assert(!qd::iterative_method<ReferenceOnly, double> && !integrable_unit<ReferenceOnly>);
// Today's calls with a generic lambda as fn still resolve: no overload calls an argument on trial.
static_assert(qd::integrate<qd::Romberg>(generic_square, 0.0, 1.0, { .rel_tol = 1e-12 }).has_value());
static_assert(qd::integrate(qd::Romberg{}, generic_square, 0.0, 1.0).has_value());
static_assert(qd::integrate(generic_square, 0.0, 1.0).has_value());
static_assert(qd::integrate<qd::GaussLegendre<3>>(generic_square, 0.0, 1.0).has_value());

// --- The forms ---------------------------------------------------------------------------------------------------------

// A method named and the same method as a value are one computation.
static_assert(qd::integrate<qd::GaussLegendre<7>>(runge, 0.0, 1.0) == qd::integrate(qd::GaussLegendre<7>{}, runge, 0.0, 1.0));
static_assert(qd::integrate<qd::Romberg>(runge, 0.0, 1.0, { .rel_tol = 1e-10 })
              == qd::integrate(qd::Romberg{}, runge, 0.0, 1.0, { .rel_tol = 1e-10 }));
static_assert(qd::integrate(runge, 0.0, 1.0) == qd::integrate(gk15{}, runge, 0.0, 1.0));
static_assert(qd::integrate<qd::TanhSinh>(runge, 0.0, 1.0) == qd::integrate(qd::TanhSinh{}, runge, 0.0, 1.0));

TEST_CASE("quadrature: protocol: a method named and as a value")
{
    CHECK(qd::integrate<qd::TanhSinh>(runge, 0.0, 1.0) == qd::integrate(qd::TanhSinh{}, runge, 0.0, 1.0));
    CHECK(qd::integrate<qd::GaussLegendre<7>>(runge, 0.0, 1.0) == qd::integrate(qd::GaussLegendre<7>{}, runge, 0.0, 1.0));
    CHECK(qd::integrate<qd::Romberg>(runge, 0.0, 1.0, { .rel_tol = 1e-10 })
          == qd::integrate(qd::Romberg{}, runge, 0.0, 1.0, { .rel_tol = 1e-10 }));
    CHECK(qd::integrate(runge, 0.0, 1.0) == qd::integrate(gk15{}, runge, 0.0, 1.0));
}

// --- The ends in order -------------------------------------------------------------------------------------------------

// integrate(fn, b, a) is exactly -integrate(fn, a, b), the whole result, for every method. Romberg at rel_tol 0 on
// 1 / (1 + x^2) over [0, 1.25] is a case where it was not before Phase 3a: the reversed run's abs_error differed in
// its last bits (the plain sum of |terms| behind the floor, summed in the other order).
static_assert(qd::integrate<qd::Romberg>(runge, 1.25, 0.0, { .rel_tol = 0.0 })
              == negated(qd::integrate<qd::Romberg>(runge, 0.0, 1.25, { .rel_tol = 0.0 })));
static_assert(qd::integrate<qd::TanhSinh>(runge, 1.25, 0.0, { .rel_tol = 0.0 })
              == negated(qd::integrate<qd::TanhSinh>(runge, 0.0, 1.25, { .rel_tol = 0.0 })));
// Infinite ends too, for the method that takes them (PHASE4_API.md A): ordered as (0, inf) and (-inf, inf).
static_assert(qd::integrate<qd::TanhSinh>(runge, dinf, 0.0) == negated(qd::integrate<qd::TanhSinh>(runge, 0.0, dinf)));
static_assert(qd::integrate<qd::TanhSinh>(runge, dinf, -dinf)
              == negated(qd::integrate<qd::TanhSinh>(runge, -dinf, dinf)));

TEST_CASE("quadrature: protocol: reversed ends give exactly the negated result")
{
    auto const peak = [](double x) { return 1.0 / ((x - 0.3) * (x - 0.3) + 1e-4); };
    for (double tol : { 1e-6, 1e-10, 0.0 }) {
        qd::Criteria<double> const c{ .rel_tol = tol };
        CHECK(qd::integrate<qd::Romberg>(runge, 1.25, 0.0, c) == negated(qd::integrate<qd::Romberg>(runge, 0.0, 1.25, c)));
        CHECK(qd::integrate<qd::Romberg>(peak, 1.0, -0.5, c) == negated(qd::integrate<qd::Romberg>(peak, -0.5, 1.0, c)));
        CHECK(qd::integrate(peak, 1.0, -0.5, c) == negated(qd::integrate(peak, -0.5, 1.0, c)));
        CHECK(qd::integrate<qd::TanhSinh>(peak, 1.0, -0.5, c)
              == negated(qd::integrate<qd::TanhSinh>(peak, -0.5, 1.0, c)));
        CHECK(qd::integrate<Trapezoid>(runge, 2.0, -1.0, c) == negated(qd::integrate<Trapezoid>(runge, -1.0, 2.0, c)));
    }
    auto const g = qd::integrate<qd::GaussLegendre<20>>(peak, 1.0, -0.5);
    REQUIRE(g.has_value());
    CHECK(*g == -*qd::integrate<qd::GaussLegendre<20>>(peak, -0.5, 1.0));
    CHECK(*qd::integrate(Midpoint{ .panels = 3 }, square, 3.0, 0.0) == -8.75);

    // integrate()'s comment for a fixed rule as a value: Midpoint{.panels = 1000} on e^x over [0, 1]. The composite
    // midpoint rule errs by -(b - a) h^2 f''(xi) / 24 (h = 1/1000, f'' = e^x in [1, e]), so it falls short of e - 1 by
    // between h^2 / 24 = 4.2e-8 and e h^2 / 24 = 1.2e-7; its plain sum of 1000 terms rounds by at most 1000 eps (e - 1)
    // = 3.8e-13, which the margin 1e-12 covers.
    auto const m = qd::integrate(Midpoint{ .panels = 1000 }, [](double x) { return std::exp(x); }, 0.0, 1.0);
    REQUIRE(m.has_value());
    double const shortfall = static_cast<double>(bin100(bin100(e_minus_1) - bin100(*m)));
    CHECK(shortfall >= 1e-6 / 24.0 - 1e-12);
    CHECK(shortfall <= 2.7183 * 1e-6 / 24.0 + 1e-12);
}

// --- What a method receives, and the validation before it --------------------------------------------------------------

TEST_CASE("quadrature: protocol: an estimating method receives the ends in order and the validated criteria")
{
    Seen seen;
    auto const r = qd::integrate(Spy{ &seen }, square, 2.0, -1.0, { .abs_tol = 1e300, .rel_tol = 5.0 });
    REQUIRE(r.has_value());
    CHECK(seen.calls == 1);
    CHECK(seen.a == -1.0);
    CHECK(seen.b == 2.0);
    CHECK(seen.abs_tol == 1e300); // unclamped: it enters only a maximum
    CHECK(seen.rel_tol == 1.0);   // clamped to 1
    CHECK(r->value == -1.0);      // the spy's 1, negated for the reversed ends

    // The validation comes first: the spy is never called.
    seen = Seen{};
    CHECK(qd::integrate(Spy{ &seen }, square, 0.0, dnan) == std::unexpected(non_finite));
    CHECK(qd::integrate(Spy{ &seen }, square, -dinf, 1.0) == std::unexpected(non_finite));
    CHECK(qd::integrate(Spy{ &seen }, square, 0.0, 1.0, { .abs_tol = -1.0 }) == std::unexpected(non_finite));
    CHECK(qd::integrate(Spy{ &seen }, square, 0.0, 1.0, { .rel_tol = dnan }) == std::unexpected(non_finite));
    CHECK(qd::integrate(Spy{ &seen }, square, 2.0, 2.0) == qd::Integral<double>{ 0.0, 0.0 });
    CHECK(seen.calls == 0);

    // A method that opts in to infinite ends receives them in order (b - a and a/2 + b/2 are its to avoid), the
    // result negated for reversed ends; a NaN end and the same infinity twice never reach it.
    CHECK(qd::integrate(UnboundedSpy{ &seen }, square, dinf, 0.0) == qd::Integral<double>{ -1.0, 0.0 });
    CHECK(seen.calls == 1);
    CHECK(seen.a == 0.0);
    CHECK(seen.b == dinf);
    CHECK(qd::integrate(UnboundedSpy{ &seen }, square, dinf, -dinf) == qd::Integral<double>{ -1.0, 0.0 });
    CHECK(seen.calls == 2);
    CHECK(seen.a == -dinf);
    CHECK(seen.b == dinf);
    CHECK(qd::integrate(UnboundedSpy{ &seen }, square, -dinf, 1.0) == qd::Integral<double>{ 1.0, 0.0 });
    CHECK(seen.calls == 3);
    CHECK(seen.a == -dinf);
    CHECK(seen.b == 1.0);
    CHECK(qd::integrate(UnboundedSpy{ &seen }, square, dnan, dinf) == std::unexpected(non_finite));
    CHECK(qd::integrate(UnboundedSpy{ &seen }, square, -dinf, -dinf) == qd::Integral<double>{ 0.0, 0.0 });
    CHECK(seen.calls == 3);
}

// --- The checked results ----------------------------------------------------------------------------------------------

using broken = std::expected<qd::Integral<double>, qd::Failure<double>>;
static_assert(qd::integrate(Returning<broken>{ qd::Integral<double>{ dnan, 0.0 } }, one, 0.0, 1.0)
              == std::unexpected(non_finite));
static_assert(qd::integrate(Returning<broken>{ qd::Integral<double>{ 1.0, -1.0 } }, one, 0.0, 1.0)
              == std::unexpected(non_finite));
static_assert(qd::integrate(Returning<broken>{ qd::Integral<double>{ 1.0, dinf } }, one, 0.0, 1.0)
              == std::unexpected(non_finite));
// A best only with not_converged, and only a sound one.
static_assert(qd::integrate(Returning<broken>{ std::unexpected(qd::Failure<double>{ qd::error::non_finite, two }) }, one,
                  0.0, 1.0)
              == std::unexpected(non_finite));
static_assert(qd::integrate(Returning<broken>{ std::unexpected(
                                qd::Failure<double>{ qd::error::not_converged, qd::Integral<double>{ dnan, 0.0 } }) },
                  one, 0.0, 1.0)
              == std::unexpected(non_finite));
static_assert(qd::integrate(Returning<broken>{ std::unexpected(qd::Failure<double>{ qd::error::not_converged, two }) }, one,
                  0.0, 1.0)
              == std::unexpected(qd::Failure<double>{ qd::error::not_converged, two }));
static_assert(qd::integrate(NaNRule{}, one, 0.0, 1.0) == std::unexpected(qd::error::non_finite));

TEST_CASE("quadrature: protocol: a result that breaks the contract is non_finite")
{
    CHECK(qd::integrate(Returning<broken>{ qd::Integral<double>{ dnan, 0.0 } }, one, 0.0, 1.0) == std::unexpected(non_finite));
    CHECK(qd::integrate(Returning<broken>{ qd::Integral<double>{ 1.0, -1.0 } }, one, 0.0, 1.0)
          == std::unexpected(non_finite));
    CHECK(qd::integrate(NaNRule{}, one, 0.0, 1.0) == std::unexpected(qd::error::non_finite));
    CHECK(qd::integrate(Scripted{ { qd::Step<double>{ .current = qd::Integral<double>{ dnan, 0.0 } }, miss, miss, miss } },
              one, 0.0, 1.0, { .rel_tol = 0.0 })
          == std::unexpected(non_finite));
}

// --- The loop of an iterative method ----------------------------------------------------------------------------------

// Each branch, by the scripted method, in constant evaluation (the run-time twins are in the test case below).
// Not ready: it goes on; ready and passing: success.
static_assert(qd::integrate(Scripted{ { qd::Step<double>{ .current = two, .ready = false }, pass, miss, miss } }, one, 0.0,
                  1.0, loose)
              == two);
// Ready with no current, not exhausted: it goes on.
static_assert(qd::integrate(Scripted{ { qd::Step<double>{}, pass, miss, miss } }, one, 0.0, 1.0, loose) == two);
// Passing but untrusted: not_converged with it as the best.
static_assert(qd::integrate(Scripted{ { qd::Step<double>{ .current = two, .trusted = false }, miss, miss, miss } }, one,
                  0.0, 1.0, loose)
              == std::unexpected(qd::Failure<double>{ qd::error::not_converged, two }));
// Exhausted and failing: not_converged with the best; without one where current is empty; exhausted before ready.
static_assert(qd::integrate(Scripted{ { miss, qd::Step<double>{ .current = wide_two, .exhausted = true }, miss, miss } },
                  one, 0.0, 1.0, loose)
              == std::unexpected(qd::Failure<double>{ qd::error::not_converged, wide_two }));
static_assert(qd::integrate(Scripted{ { miss, qd::Step<double>{ .exhausted = true }, miss, miss } }, one, 0.0, 1.0, loose)
              == std::unexpected(qd::Failure<double>{ qd::error::not_converged, std::nullopt }));
static_assert(qd::integrate(Scripted{ { qd::Step<double>{ .current = two, .ready = false, .exhausted = true }, miss, miss,
                                  miss } },
                  one, 0.0, 1.0, loose)
              == std::unexpected(qd::Failure<double>{ qd::error::not_converged, two }));
// Passing before settle and failing after: it goes on from the settled step, here to a success at the next step.
static_assert(qd::integrate(Scripted{ { pass, pass, miss, miss }, 0, miss }, one, 0.0, 1.0, loose) == two);
// Passing, then failing after settle while exhausted (Phase 2's double settle): not_converged with the settled best.
static_assert(qd::integrate(Scripted{ { qd::Step<double>{ .current = two, .exhausted = true }, miss, miss, miss }, 0,
                                  qd::Step<double>{ .current = wide_two, .exhausted = true } },
                  one, 0.0, 1.0, loose)
              == std::unexpected(qd::Failure<double>{ qd::error::not_converged, wide_two }));
// A Failure from start, advance or settle is returned as it is.
static_assert(qd::integrate(Scripted{ .steps = { pass, miss, miss, miss }, .fail_start = Scripted::scripted_failure },
                  one, 0.0, 1.0, loose)
              == std::unexpected(Scripted::scripted_failure));
static_assert(qd::integrate(Scripted{ .steps = { miss, miss, miss, miss }, .fail_advance_at = 1 }, one, 0.0, 1.0, loose)
              == std::unexpected(Scripted::scripted_failure));
static_assert(qd::integrate(Scripted{ .steps = { miss, pass, miss, miss }, .fail_settle_at = 1 }, one, 0.0, 1.0, loose)
              == std::unexpected(Scripted::scripted_failure));
// A step that cannot be tested is non_finite, before any arithmetic on it (a NaN would not be a constant expression):
// a current that is not finite, an attainable error that is infinite, NaN or negative.
static_assert(qd::integrate(Scripted{ { qd::Step<double>{ .current = qd::Integral<double>{ dnan, 0.0 } }, miss, miss,
                                  miss } },
                  one, 0.0, 1.0, { .rel_tol = 0.0 })
              == std::unexpected(non_finite));
static_assert(qd::integrate(Scripted{ { qd::Step<double>{ .current = two, .attainable = dinf }, miss, miss, miss } }, one,
                  0.0, 1.0, loose)
              == std::unexpected(non_finite));
static_assert(qd::integrate(Scripted{ { qd::Step<double>{ .current = two, .attainable = dnan }, miss, miss, miss } }, one,
                  0.0, 1.0, loose)
              == std::unexpected(non_finite));
static_assert(qd::integrate(Scripted{ { qd::Step<double>{ .current = two, .attainable = -1.0 }, miss, miss, miss } }, one,
                  0.0, 1.0, loose)
              == std::unexpected(non_finite));
// The attainable error lets a ready step stop: an estimate of 1 meets max(0, 1e-2 * 2, attainable 1).
static_assert(qd::integrate(Scripted{ { qd::Step<double>{ .current = wide_two, .attainable = 1.0 }, miss, miss, miss } },
                  one, 0.0, 1.0, loose)
              == wide_two);

TEST_CASE("quadrature: protocol: the loop's branches")
{
    int settles = 0;
    int advances = 0;
    auto const counted = [&](Scripted s) {
        s.settles = &settles;
        s.advances = &advances;
        return s;
    };
    auto const reset = [&] { settles = advances = 0; };

    CHECK(qd::integrate(counted(Scripted{ { qd::Step<double>{ .current = two, .ready = false }, pass, miss, miss } }), one,
              0.0, 1.0, loose)
          == two);
    CHECK(settles == 1);
    CHECK(advances == 1);

    reset();
    CHECK(qd::integrate(counted(Scripted{ { miss, qd::Step<double>{ .exhausted = true }, miss, miss } }), one, 0.0, 1.0,
              loose)
          == std::unexpected(qd::Failure<double>{ qd::error::not_converged, std::nullopt }));
    CHECK(advances == 1); // never advanced from the exhausted step

    // Passing before settle and failing after it: settle, go on, then pass and settle again.
    reset();
    CHECK(qd::integrate(counted(Scripted{ { pass, pass, miss, miss }, 0, miss }), one, 0.0, 1.0, loose) == two);
    CHECK(settles == 2);
    CHECK(advances == 1);

    // Phase 2's double settle: passing, failing after settle, exhausted: settled twice, the settled best returned.
    reset();
    CHECK(qd::integrate(counted(Scripted{ { qd::Step<double>{ .current = two, .exhausted = true }, miss, miss, miss }, 0,
              qd::Step<double>{ .current = wide_two, .exhausted = true } }),
              one, 0.0, 1.0, loose)
          == std::unexpected(qd::Failure<double>{ qd::error::not_converged, wide_two }));
    CHECK(settles == 2);
    CHECK(advances == 0);

    // The same branches with an observer: the same decisions, settle forwarded, a report per start and advance.
    reset();
    int reports = 0;
    auto const count_reports = [&reports](qd::Step<double> const &) { ++reports; };
    CHECK(qd::integrate(counted(Scripted{ { pass, pass, miss, miss }, 0, miss }), qd::observer(count_reports), one, 0.0,
              1.0, loose)
          == two);
    CHECK(settles == 2);
    CHECK(reports == 2);
    reset();
    reports = 0;
    CHECK(qd::integrate(counted(Scripted{ { qd::Step<double>{ .current = two, .exhausted = true }, miss, miss, miss }, 0,
                            qd::Step<double>{ .current = wide_two, .exhausted = true } }),
              qd::observer(count_reports), one, 0.0, 1.0, loose)
          == std::unexpected(qd::Failure<double>{ qd::error::not_converged, wide_two }));
    CHECK(settles == 2);
    CHECK(reports == 1);
}

TEST_CASE("quadrature: protocol: the loop's branches at run time")
{
    // The run-time twins of the static_asserts above that have none in the test case before.
    CHECK(qd::integrate(Scripted{ { qd::Step<double>{}, pass, miss, miss } }, one, 0.0, 1.0, loose) == two);
    CHECK(qd::integrate(Scripted{ { qd::Step<double>{ .current = two, .trusted = false }, miss, miss, miss } }, one, 0.0,
              1.0, loose)
          == std::unexpected(qd::Failure<double>{ qd::error::not_converged, two }));
    CHECK(qd::integrate(Scripted{ { miss, qd::Step<double>{ .current = wide_two, .exhausted = true }, miss, miss } }, one,
              0.0, 1.0, loose)
          == std::unexpected(qd::Failure<double>{ qd::error::not_converged, wide_two }));
    CHECK(qd::integrate(Scripted{ { qd::Step<double>{ .current = two, .ready = false, .exhausted = true }, miss, miss,
                            miss } },
              one, 0.0, 1.0, loose)
          == std::unexpected(qd::Failure<double>{ qd::error::not_converged, two }));
    CHECK(qd::integrate(Scripted{ .steps = { pass, miss, miss, miss }, .fail_start = Scripted::scripted_failure }, one, 0.0,
              1.0, loose)
          == std::unexpected(Scripted::scripted_failure));
    CHECK(qd::integrate(Scripted{ .steps = { miss, miss, miss, miss }, .fail_advance_at = 1 }, one, 0.0, 1.0, loose)
          == std::unexpected(Scripted::scripted_failure));
    CHECK(qd::integrate(Scripted{ .steps = { miss, pass, miss, miss }, .fail_settle_at = 1 }, one, 0.0, 1.0, loose)
          == std::unexpected(Scripted::scripted_failure));
    for (double attainable : { dinf, dnan, -1.0 }) {
        CAPTURE(attainable);
        CHECK(qd::integrate(Scripted{ { qd::Step<double>{ .current = two, .attainable = attainable }, miss, miss, miss } },
                  one, 0.0, 1.0, loose)
              == std::unexpected(non_finite));
    }
    CHECK(qd::integrate(Scripted{ { qd::Step<double>{ .current = wide_two, .attainable = 1.0 }, miss, miss, miss } }, one,
              0.0, 1.0, loose)
          == wide_two);
}

// --- The observer -----------------------------------------------------------------------------------------------------

// An observer changes nothing, in both forms, reversed ends included.
static_assert(qd::integrate(qd::observer(quiet), runge, 0.0, 1.0) == qd::integrate(runge, 0.0, 1.0));
static_assert(qd::integrate<qd::Romberg>(qd::observer(quiet), runge, 1.0, 0.0, { .rel_tol = 1e-10 })
              == qd::integrate<qd::Romberg>(runge, 1.0, 0.0, { .rel_tol = 1e-10 }));
static_assert(qd::integrate(qd::Romberg{}, qd::observer(quiet), runge, 1.0, 0.0, { .rel_tol = 1e-10 })
              == qd::integrate(qd::Romberg{}, runge, 1.0, 0.0, { .rel_tol = 1e-10 }));
static_assert(qd::integrate(gk15{}, qd::observer(quiet), cubic, 2.0, -1.0) == qd::integrate(cubic, 2.0, -1.0));
static_assert(qd::integrate<qd::TanhSinh>(qd::observer(quiet), runge, 1.0, 0.0)
              == qd::integrate<qd::TanhSinh>(runge, 1.0, 0.0));
static_assert(qd::integrate(qd::TanhSinh{}, qd::observer(quiet), runge, 1.0, 0.0)
              == qd::integrate(qd::TanhSinh{}, runge, 1.0, 0.0));
// An observed TanhSinh takes infinite ends: Observed forwards the opt-in.
static_assert(qd::integrate<qd::TanhSinh>(qd::observer(quiet), runge, dinf, -dinf)
              == qd::integrate<qd::TanhSinh>(runge, dinf, -dinf));
// In constant evaluation an observer may record, here into a transient vector.
static_assert([] {
    std::vector<qd::Step<double>> steps;
    auto record = [&steps](qd::Step<double> const &s) { steps.push_back(s); };
    auto const r = qd::integrate<qd::Romberg>(qd::observer(record), runge, 0.0, 1.0, { .rel_tol = 1e-10 });
    return r && !steps.empty() && steps.back().current == *r;
}());
static_assert([] {
    std::vector<qd::Step<double>> steps;
    auto record = [&steps](qd::Step<double> const &s) { steps.push_back(s); };
    auto const r = qd::integrate<qd::TanhSinh>(qd::observer(record), runge, 0.0, 1.0);
    return r && steps.size() == 5 && !steps.front().current && steps.back().current == *r;
}());

TEST_CASE("quadrature: protocol: an observer reports each step and changes nothing")
{
    SUBCASE("Romberg: a step per level, level k at 2^k + 1 evaluations, the last step the result")
    {
        int count = 0;
        auto const e = Counted<double (*)(double)>{ [](double x) { return std::exp(x); }, &count };
        std::vector<std::pair<int, qd::Step<double>>> steps;
        auto const record = [&](qd::Step<double> const &s) { steps.emplace_back(count, s); };
        auto const r = qd::integrate<qd::Romberg>(qd::observer(record), e, 0.0, 1.0, { .rel_tol = 1e-12 });
        int const total = count;
        count = 0;
        REQUIRE(r.has_value());
        CHECK(r == qd::integrate<qd::Romberg>(e, 0.0, 1.0, { .rel_tol = 1e-12 }));
        CHECK(count == total);
        REQUIRE(!steps.empty());
        for (std::size_t i = 0; i < steps.size(); ++i) {
            CAPTURE(i);
            CHECK(steps[i].first == (1 << (i + 1)) + 1);
            CHECK(steps[i].second.ready == (i + 1 >= 4));
            CHECK(!steps[i].second.exhausted);
        }
        CHECK(steps.back().second.current == *r);
    }
    SUBCASE("Gauss-Kronrod: the mesh, then 30 evaluations per bisection; a mesh-only run's step is the result")
    {
        int count = 0;
        auto const peak = Counted<double (*)(double)>{ [](double x) { return 1.0 / ((x - 0.3) * (x - 0.3) + 1e-6); },
            &count };
        std::vector<int> evaluations;
        auto const record = [&](qd::Step<double> const &) { evaluations.push_back(count); };
        auto const r = qd::integrate<gk15>(qd::observer(record), peak, 0.0, 1.0, { .rel_tol = 1e-10 });
        REQUIRE(r.has_value());
        REQUIRE(evaluations.size() > 1);
        for (std::size_t i = 0; i < evaluations.size(); ++i) {
            CAPTURE(i);
            CHECK(evaluations[i] == 150 + 30 * static_cast<int>(i)); // no descent: fn is finite everywhere
        }
        // 1000 (atan 700 + atan 300) (mpmath, 80 digits); the estimate is honest.
        CHECK(error_of(r->value, "3136.830762145301293399296335784189135890259424004257320130132129794127") <= r->abs_error);

        // x^3 - x is integrated exactly by the mesh, so the run stops at its first step; the running and fresh sums add
        // the same panels in the same order, so the reported step is the result.
        std::vector<qd::Step<double>> steps;
        auto const keep = [&steps](qd::Step<double> const &s) { steps.push_back(s); };
        auto const c = qd::integrate(gk15{}, qd::observer(keep), cubic, -1.0, 2.0);
        REQUIRE(c.has_value());
        REQUIRE(steps.size() == 1);
        CHECK(steps.front().current == *c);
    }
    SUBCASE("an lvalue observer is kept by reference; a mutable lambda counts")
    {
        struct Recorder
        {
            std::vector<qd::Step<double>> steps;
            void operator()(qd::Step<double> const &s) { steps.push_back(s); }
        } recorder;
        auto const r = qd::integrate<qd::Romberg>(qd::observer(recorder), runge, 0.0, 1.0, { .rel_tol = 1e-10 });
        REQUIRE(r.has_value());
        CHECK(!recorder.steps.empty());
        CHECK(recorder.steps.back().current == *r);

        int calls = 0;
        auto counting = [&calls, n = 0](qd::Step<double> const &) mutable { calls = ++n; };
        (void)qd::integrate(qd::Romberg{}, qd::observer(counting), runge, 0.0, 1.0, { .rel_tol = 1e-10 });
        CHECK(calls == static_cast<int>(recorder.steps.size()));
    }
    SUBCASE("a method that cannot be copied is observed, and reports Romberg's steps")
    {
        std::vector<qd::Step<double>> steps, romberg;
        auto const keep = [&steps](qd::Step<double> const &s) { steps.push_back(s); };
        auto const keep_romberg = [&romberg](qd::Step<double> const &s) { romberg.push_back(s); };
        MoveOnly const m;
        auto const r = qd::integrate(m, qd::observer(keep), runge, 0.0, 1.0, { .rel_tol = 1e-10 });
        auto const q = qd::integrate<qd::Romberg>(qd::observer(keep_romberg), runge, 0.0, 1.0, { .rel_tol = 1e-10 });
        REQUIRE(r.has_value());
        CHECK(r == q);
        REQUIRE(steps.size() == romberg.size());
        CHECK(std::ranges::equal(steps, romberg, [](auto const &x, auto const &y) {
            return x.current == y.current && x.attainable == y.attainable && x.ready == y.ready
                   && x.exhausted == y.exhausted && x.trusted == y.trusted;
        }));
    }
    SUBCASE("the reports are for the ends in order")
    {
        std::vector<qd::Step<double>> steps;
        auto const keep = [&steps](qd::Step<double> const &s) { steps.push_back(s); };
        auto const r = qd::integrate<qd::Romberg>(qd::observer(keep), runge, 1.0, 0.0, { .rel_tol = 1e-10 });
        REQUIRE(r.has_value());
        REQUIRE(!steps.empty());
        CHECK(steps.back().current->value == -r->value);
    }
}

// --- Methods of one's own -------------------------------------------------------------------------------------------

// The midpoint rule errs by (b - a) h^2 f'' / 24, exact for a quadratic: x^2 over [0, 3] with 3 panels is 9 - 1/4, and
// every node (0.5, 1.5, 2.5) and sum is exact in double.
static_assert(*qd::integrate(Midpoint{ .panels = 3 }, square, 0.0, 3.0) == 8.75);
static_assert(*qd::integrate(Midpoint{ .panels = 3 }, square, 3.0, 0.0) == -8.75);
static_assert(qd::integrate<Midpoint>(square, 0.0, dinf) == std::unexpected(qd::error::non_finite));

// The trapezoidal rule on x^3 - x over [-1, 2] errs by exactly (h^2 / 12) (f'(2) - f'(-1)) = 3 h^2 / 4
// (Euler-Maclaurin ends after one term: f''' is constant). Its change between levels is 9 h^2 / 4, at most
// 1e-6 * 9/4 first at level 12 (h = 3 / 4096), where every node, product and sum is exact in double: 9/4 + 27 / 2^26.
static_assert(*qd::integrate<Trapezoid>(cubic, -1.0, 2.0, { .rel_tol = 1e-6 }) == qd::Integral<double>{ 2.25 + 27.0 / 67108864.0, 81.0 / 67108864.0 });

TEST_CASE("quadrature: protocol: methods of one's own")
{
    CHECK(*qd::integrate(Midpoint{ .panels = 3 }, square, 0.0, 3.0) == 8.75);

    SUBCASE("a fixed rule evaluates fn only inside the interval")
    {
        int outside = 0;
        auto const quarter = Confined{ [](double) { return 0.25; }, -dmax, dmax, &outside };
        CHECK(qd::integrate(Midpoint{}, quarter, -dmax, dmax) == std::unexpected(qd::error::non_finite));
        auto const linear = Confined{ [](double x) { return x; }, 0.0, 1.0, &outside };
        CHECK(qd::integrate(Midpoint{ .panels = 7 }, linear, 1.0, 0.0).has_value());
        CHECK(outside == 0);
    }
    SUBCASE("an iterative method: e^x over [0, 1]")
    {
        // The trapezoidal error (e - 1)(h^2 / 12 - h^4 / 720 + ...) is 3.3e-11 at level 16, h = 2^-16, where the change
        // 3 times that first meets 1e-10 * (e - 1); the plain sum's rounding adds at most about n u sum |f| h = 1.2e-11.
        auto const e = [](double x) { return std::exp(x); };
        int levels = 0;
        auto const r = qd::integrate(Trapezoid{}, qd::observer([&levels](qd::Step<double> const &) { ++levels; }), e, 0.0,
            1.0, { .rel_tol = 1e-10 });
        REQUIRE(r.has_value());
        CHECK(levels == 16);
        CHECK(error_of(r->value, e_minus_1) <= 4.5e-11);
        CHECK(error_of(r->value, e_minus_1) <= r->abs_error);
        CHECK(r->abs_error <= 1e-10 * 1.7182818284590452);
    }
    SUBCASE("an iterative method stops where its step cannot be represented")
    {
        // On [0, 12 D], D = denorm_min, the trapezoidal rule's step 12 D / 2^k is exact to level 2 (3 D); at level 3,
        // 1.5 D would round to 2 D and sample 14 D. So the method reports exhausted at level 2, before it is ready:
        // not_converged, with the level-2 sum of the constant 1, exactly 12 D, as the best, and fn never outside.
        int outside = 0;
        auto const constant = Confined{ [](double) { return 1.0; }, 0.0, 12.0 * dmin, &outside };
        auto const r = qd::integrate(Trapezoid{}, constant, 0.0, 12.0 * dmin);
        REQUIRE(!r.has_value());
        CHECK(r.error().reason == qd::error::not_converged);
        REQUIRE(r.error().best.has_value());
        CHECK(r.error().best->value == 12.0 * dmin);
        CHECK(outside == 0);
    }
    SUBCASE("a combinator over a library method")
    {
        // x = u^2 turns 1/sqrt x into g(u) = 2 u / fl(sqrt(fl(u^2))) = 2 within about 3 ulps per evaluation (u^2,
        // sqrt, the division; the product by 2 is exact): K15's weights sum to 2 per panel, so the compensated sums
        // keep the result within 8 eps of 2.
        auto const inv_sqrt = [](double x) { return 1.0 / std::sqrt(x); };
        auto const r = qd::integrate(SquareAtLowerEnd<gk15>{}, inv_sqrt, 0.0, 1.0, { .rel_tol = 1e-12 });
        REQUIRE(r.has_value());
        CHECK(std::fabs(r->value - 2.0) <= 8.0 * deps);
        CHECK(qd::integrate(SquareAtLowerEnd<gk15>{}, inv_sqrt, 1.0, 0.0, { .rel_tol = 1e-12 }) == negated(r));
        // cos x / sqrt x becomes 2 cos(u^2), smooth: the estimate covers the error and meets the request.
        auto const cos_sqrt = [](double x) { return std::cos(x) / std::sqrt(x); };
        auto const c = qd::integrate(SquareAtLowerEnd<gk15>{}, cos_sqrt, 0.0, 1.0, { .rel_tol = 1e-12 });
        REQUIRE(c.has_value());
        CHECK(error_of(c->value, "1.809048475800544162949576733665114291615983641191003479196640074902159") <= c->abs_error);
        CHECK(c->abs_error <= 1e-12 * 1.81);
        // Over Romberg, which samples u = 0 where the substituted integrand is 0 times inf: non_finite, as Romberg is
        // documented to be on a singular end.
        CHECK(!qd::integrate(SquareAtLowerEnd<qd::Romberg>{}, inv_sqrt, 0.0, 1.0).has_value());
    }
    SUBCASE("multiprecision, with and without expression templates")
    {
        auto const inv_sqrt50 = [](bin50 const &x) -> bin50 { return bin50(1) / sqrt(x); };
        auto const r = qd::integrate(SquareAtLowerEnd<gk15>{}, inv_sqrt50, bin50(0), bin50(1), { .rel_tol = 1e-40 });
        REQUIRE(r.has_value());
        CHECK(error_of(r->value, "2") <= 8.0 * eps_of<bin50>());
        auto const cos_sqrt_et = [](bin50_et const &x) -> bin50_et { return bin50_et(cos(x) / sqrt(x)); };
        auto const c = qd::integrate(SquareAtLowerEnd<gk15>{}, cos_sqrt_et, bin50_et(0), bin50_et(1), { .rel_tol = 1e-30 });
        REQUIRE(c.has_value());
        CHECK(error_of(c->value, "1.809048475800544162949576733665114291615983641191003479196640074902159") <= c->abs_error);
        CHECK(c->abs_error <= 1e-30 * 1.81);
    }
}
