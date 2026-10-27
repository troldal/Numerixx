// Numerical derivatives in nxx::roots (Phase 3, docs/roots/PHASE3_API.md): derivative callables that return a
// std::expected, nxx::deriv::derivative_of among them, and a failed derivative counting as one that is not finite.
//
// The counts of the numerical derivatives are calls of the user's f (the value and every stencil node), against an
// independent Python implementation of the specification and of nxx::deriv's Central<2> with its default step
// (.scratch/roots3/proto3.py), for + - * / and sqrt (correctly rounded) in double, assuming no FMA contraction, as in
// Phase 2. Every other case checks accuracy (roots_check.hpp).

#include "roots_check.hpp"

#include <numerixx/deriv.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <expected>
#include <limits>
#include <optional>
#include <string>
#include <utility>

using namespace roots_check;
namespace dv = nxx::deriv;

namespace {

    // The prototype's expressions, in the same order.
    constexpr auto sq = [](double x) { return x * x - 2.0; };
    constexpr auto wal = [](double x) { return (x * x - 2.0) * x - 5.0; };
    constexpr auto triple = [](double x) { return (x - 1.0) * (x - 1.0) * (x - 1.0); };
    constexpr auto cube_e3 = [](double x) { return x * x * x - 1e-3; };
    constexpr auto reciprocal = [](double x) { return 1.0 / x - 1.0; };
    constexpr auto x20 = [](double x) {
        double const x2 = x * x;
        double const x4 = x2 * x2;
        double const x8 = x4 * x4;
        return x8 * x8 * x4 - 1.0;
    };
    auto const sqrt_fn = [](double x) { return std::sqrt(x) - 0.5; };
    auto const ln = [](double x) { return std::log(x / 1e-12); };

    constexpr auto dsq = [](double x) { return 2.0 * x; };
    constexpr auto d2sq = [](double) { return 2.0; };
    constexpr auto dwal = [](double x) { return 3.0 * x * x - 2.0; };

    // fn with numerical derivatives (Central<2>: f', and f'' for Halley), every call of f counted.
    template<typename M, typename F, typename S>
    std::pair<std::expected<rt::Root<double>, rt::error>, int> numerical(F fn, S start)
    {
        int n = 0;
        Counted<F> const f{ fn, &n };
        if constexpr (std::same_as<M, rt::Newton>) {
            auto const r = rt::solve<M>(f, dv::derivative_of(f), start);
            return { r, n };
        } else {
            auto const r = rt::solve<M>(f, dv::derivative_of(f), dv::derivative_of<dv::Central<2>, 2>(f), start);
            return { r, n };
        }
    }

    // fn with the given separate derivatives, the calls of fn counted.
    template<typename M, typename F, typename S, typename... Ds>
    std::pair<std::expected<rt::Root<double>, rt::error>, int> separate(F fn, S start, Ds const &...ds)
    {
        int n = 0;
        auto const r = rt::solve<M>(Counted<F>{ fn, &n }, ds..., start);
        return { r, n };
    }

    template<typename M, typename F, typename D, typename S>
    concept solvable_with = requires(F f, D d, S s) { rt::solve<M>(f, d, s); };

    template<typename M, typename F, typename D, typename D2, typename S>
    concept solvable_with_two = requires(F f, D d, D2 d2, S s) { rt::solve<M>(f, d, d2, s); };

    template<typename M, typename F, typename S>
    concept solvable = requires(F f, S s) { rt::solve<M>(f, s); };

    using Value = double (*)(double);
    using B = rt::Bracket<double>;

} // namespace

// Constant evaluation: derivative_of, a derivative returning an expected, and a failed one (C). The conditioning term
// of x^2 - 2 is 7.9e-17 (test_roots_open.cpp).
static_assert([] {
    auto const r = rt::solve<rt::Newton>(sq, dv::derivative_of(sq), 1.0);
    return r && nxx::math::abs(r->value - sqrt2) <= r->abs_error + 8e-17;
}());
static_assert([] {
    auto const r = rt::solve<rt::Halley>(sq, dv::derivative_of(sq), dv::derivative_of<dv::Central<2>, 2>(sq), 1.0);
    return r && nxx::math::abs(r->value - sqrt2) <= r->abs_error + 8e-17;
}());
static_assert([] {
    auto const r = rt::solve<rt::Newton>(sq, dv::derivative_of(sq), rt::Bracket{ 1.0, 2.0 });
    return r && nxx::math::abs(r->value - sqrt2) <= r->abs_error + 8e-17;
}());
static_assert([] {
    auto const sqf = [](float x) { return x * x - 2.0f; };
    auto const r = rt::solve<rt::Halley>(sqf, dv::derivative_of(sqf), dv::derivative_of<dv::Central<2>, 2>(sqf), 1.0f);
    return r && nxx::math::abs(r->value - 1.41421356f) <= r->abs_error + 1.5e-7f; // c = eps |ref| in float
}());
static_assert([] {
    auto const sql = [](long double x) { return x * x - 2.0L; };
    auto const r = rt::solve<rt::Newton>(sql, dv::derivative_of(sql), 1.0L);
    return r && nxx::math::abs(r->value - 1.41421356237309504880L) <= r->abs_error + 2.0L * 1.0842e-19L;
}());
// A derivative returning an expected gives Newton's result with the plain derivative, exactly; a failed one (always)
// leaves Newton from a guess without a step.
static_assert(rt::solve<rt::Newton>(sq, [](double x) -> std::expected<double, int> { return 2.0 * x; }, 1.0)
              == rt::solve<rt::Newton>(sq, dsq, 1.0));
static_assert(rt::solve<rt::Newton>(sq, [](double) -> std::expected<double, int> { return std::unexpected(1); }, 1.0)
              == std::unexpected(rt::error::not_converged));
static_assert([] { // on a bracket every step is a bisection step, and the bracket still finishes the root
    auto const r = rt::solve<rt::Newton>(sq, [](double) -> std::expected<double, int> { return std::unexpected(1); },
        rt::Bracket{ 1.0, 2.0 });
    return r && nxx::math::abs(r->value - sqrt2) <= r->abs_error + 8e-17;
}());

// What a separate derivative may return (P3): a value, or a std::expected whose value converts to T; fn itself not.
using ExpectedValue = std::expected<double, int> (*)(double);
using ExpectedEstimate = std::expected<dv::Estimate<double>, dv::error> (*)(double);
static_assert(solvable_with<rt::Newton, Value, ExpectedValue, double> && solvable_with<rt::Newton, Value, ExpectedValue, B>);
static_assert(solvable_with<rt::Newton, Value, ExpectedEstimate, double>);
static_assert(solvable_with_two<rt::Halley, Value, Value, ExpectedValue, double>); // mixed: a plain f', an expected f''
static_assert(solvable_with_two<rt::Halley, Value, ExpectedEstimate, ExpectedValue, B>);
static_assert(!solvable_with<rt::Newton, Value, std::expected<void, int> (*)(double), double>);
static_assert(!solvable_with<rt::Newton, Value, std::expected<std::string, int> (*)(double), double>);
static_assert(!solvable_with<rt::Newton, ExpectedValue, Value, double>); // fn returning an expected: deferred
static_assert(!solvable<rt::Brent, ExpectedValue, B>);
static_assert(!solvable_with<rt::Brent, Value, ExpectedValue, B>); // a derivative only for Newton and Halley

// The independent review of Phase 3: the constraint admitted conversions the body could not perform.
namespace {
    // Converts to double only as an rvalue: the body converted a const copy.
    struct RvalueOnly
    {
        constexpr operator double() && { return 2.0; }
    };

    // Not copyable: an expected of it returned by reference was copied.
    struct NoCopy
    {
        constexpr NoCopy() = default;
        NoCopy(NoCopy const &) = delete;
        constexpr NoCopy(NoCopy &&) = default;
    };

    struct ByReference
    {
        std::expected<double, NoCopy> e;
        constexpr std::expected<double, NoCopy> const &operator()(double const &) const { return e; }
    };
} // namespace
// f' = 2 everywhere: Newton on x^2 - 2 from 1 with that derivative, against the same constant derivative as a value.
static_assert(rt::solve<rt::Newton>(sq, [](double) -> std::expected<RvalueOnly, int> { return RvalueOnly{}; }, 1.0)
              == rt::solve<rt::Newton>(sq, [](double) { return 2.0; }, 1.0));
static_assert(rt::solve<rt::Newton>(sq, ByReference{ 2.0 }, 1.0) == rt::solve<rt::Newton>(sq, [](double) { return 2.0; }, 1.0));
// A value that converts only as an rvalue, returned by const reference: rejected by the constraint, as *r is const.
static_assert(!solvable_with<rt::Newton, Value, std::expected<RvalueOnly, int> const &(*)(double), double>);

TEST_CASE("roots:numerical: derivative_of, calls of f against the prototype")
{
    auto const check = [](auto const &fn, auto start, int newton, int halley) {
        CHECK(numerical<rt::Newton>(fn, start).second == newton);
        CHECK(numerical<rt::Halley>(fn, start).second == halley);
    };
    // A point costs 3 calls with f' (f, two nodes) and 6 with f' and f'' (f, two nodes, f again and two nodes).
    check(sq, 1.0, 21, 30);
    check(sq, 10.0, 27, 42);
    check(sq, B{ 1.0, 2.0 }, 24, 36);
    check(wal, 2.0, 18, 30);
    check(wal, 10.0, 33, 48);
    check(wal, B{ 2.0, 3.0 }, 21, 36);
    check(reciprocal, 0.5, 21, 18); // Halley: one point more than with the analytic derivatives (2)
    check(reciprocal, B{ 0.5, 3.0 }, 27, 36);
    check(x20, 0.5, 66, 66);
    check(x20, 2.0, 60, 66);
    check(x20, B{ 0.0, 5.0 }, 33, 66);
    check(sqrt_fn, 1.0, 15, 30); // Newton: one point more than with the analytic derivative (4)
    check(sqrt_fn, B{ 0.0, 1.0 }, 12, 24); // the failed derivative at 0 is evaluated, and the solve goes on
    // Where |f'| is as small as the derivative's own error (B): the outcome changes, either way.
    check(triple, B{ 0.0, 3.0 }, 549, 1008); // analytic: 91 and 54
    check(cube_e3, 0.0, 105, 204); // analytic: not_converged, f'(0) = 0
    check(cube_e3, B{ 0.0, 1.0 }, 33, 54);
}

TEST_CASE("roots:numerical: derivative_of, accuracy")
{
    // The conditioning terms of test_roots_open.cpp. x^3 - 1e-3 near 0.1: x^3's rounding, about 2 half-ulps of 1e-3
    // (2.2e-19), over f' = 0.03, 7.2e-18; and fl(1e-3)'s cube root lies 6.9e-19 above 0.1: 8e-18 in all.
    for (auto const &r : { numerical<rt::Newton>(sq, 1.0).first, numerical<rt::Halley>(sq, 10.0).first,
             numerical<rt::Newton>(sq, B{ 1.0, 2.0 }).first })
        CHECK(brackets(r, sqrt2, 8e-17));
    CHECK(brackets(numerical<rt::Newton>(wal, 10.0).first, wallis, 2e-16));
    CHECK(brackets(numerical<rt::Halley>(wal, B{ 2.0, 3.0 }).first, wallis, 2e-16));
    CHECK(brackets(numerical<rt::Newton>(reciprocal, 0.5).first, 1.0, 1.2e-16));
    CHECK(brackets(numerical<rt::Halley>(x20, 2.0).first, 1.0, 1.1e-16));
    CHECK(brackets(numerical<rt::Newton>(triple, B{ 0.0, 3.0 }).first, 1.0, 0.0));
    CHECK(brackets(numerical<rt::Newton>(cube_e3, 0.0).first, 0.1, 8e-18));
    CHECK(brackets(numerical<rt::Halley>(cube_e3, B{ 0.0, 1.0 }).first, 0.1, 8e-18));
    CHECK(brackets(numerical<rt::Newton>(sqrt_fn, B{ 0.0, 1.0 }).first, 0.25, 0.0));

    // libm: cos x - x and Kepler, with Ridders' derivative too (1.4e-16 and 2.9e-16, as in Phase 2).
    auto const dottie_fn = [](double x) { return std::cos(x) - x; };
    constexpr double kepler_root = 1.862086686874532271834067; // E - 0.9 sin E = 1, 0.9 as a double (mpmath)
    auto const kepler_fn = [](double x) { return x - 0.9 * std::sin(x) - 1.0; };
    CHECK(brackets(numerical<rt::Newton>(dottie_fn, 1.0).first, dottie, 1.4e-16));
    CHECK(brackets(numerical<rt::Halley>(kepler_fn, B{ 0.0, 3.141592653589793 }).first, kepler_root, 2.9e-16));
    CHECK(brackets(rt::solve<rt::Newton>(dottie_fn, dv::derivative_of<dv::Ridders>(dottie_fn), 1.0), dottie, 1.4e-16));
    CHECK(brackets(rt::solve<rt::Newton>(kepler_fn, dv::derivative_of<dv::Ridders>(kepler_fn), 1.0), kepler_root,
        2.9e-16));
    // A step of one's own: a lambda around derivative(fn, x, h).
    CHECK(brackets(
        rt::solve<rt::Newton>(dottie_fn, [&](double x) { return dv::derivative(dottie_fn, x, 1e-4); }, 1.0), dottie,
        1.4e-16));

    // ln(x / 1e-12) on [1e-20, 1]: every point below 6e-6 has a failed f', the root's neighbourhood included, so Newton
    // bisects throughout; Halley also passes through (6e-6, 1.2e-4), where only f'' failed. 2.5e-28: Phase 2's term.
    CHECK(brackets(numerical<rt::Newton>(ln, B{ 1e-20, 1.0 }).first, 1e-12, 2.5e-28));
    CHECK(brackets(numerical<rt::Halley>(ln, B{ 1e-20, 1.0 }).first, 1e-12, 2.5e-28));
}

TEST_CASE("roots:numerical: documented failures of a numerical derivative from a guess")
{
    auto const fails = [](auto const &fn, double x0, int newton, int halley) {
        auto const n = numerical<rt::Newton>(fn, x0);
        auto const h = numerical<rt::Halley>(fn, x0);
        CHECK(n.first == std::unexpected(rt::error::not_converged));
        CHECK(h.first == std::unexpected(rt::error::not_converged));
        CHECK(n.second == newton);
        CHECK(h.second == halley);
    };
    // (x - 1)^3 from 2: f' = 3 e^2 + h^2 stalls the steps below e = h, and the cap ends it (1 + 212 points).
    fails(triple, 2.0, 639, 1278);
    // x^2 - 2 from 1e-12: f(x - h) and f(x + h) round to the same value, f' = 0 at the guess.
    fails(sq, 1e-12, 3, 6);
    // ln(x / 1e-12) from 1e-13: the step 6e-6 reaches x < 0, so f' fails at the guess (its two nodes are evaluated).
    fails(ln, 1e-13, 3, 6);
}

TEST_CASE("roots:numerical: a failed derivative counts as one that is not finite")
{
    auto const dsq_expected = [](double x) -> std::expected<double, int> { return 2.0 * x; };
    auto const dsq_estimate = [](double x) -> std::expected<dv::Estimate<double>, dv::error> {
        return dv::Estimate<double>{ 2.0 * x, std::nullopt };
    };
    auto const failed = [](double) -> std::expected<double, int> { return std::unexpected(0); };
    auto const dwal_expected = [](double x) -> std::expected<double, int> { return 3.0 * x * x - 2.0; };

    // An expected's value is the value: the same results and counts as the plain derivative, bit for bit.
    for (double const x0 : { 1.0, 10.0 }) {
        CHECK(separate<rt::Newton>(sq, x0, dsq_expected) == separate<rt::Newton>(sq, x0, dsq));
        CHECK(separate<rt::Newton>(sq, x0, dsq_estimate) == separate<rt::Newton>(sq, x0, dsq));
        CHECK(separate<rt::Halley>(sq, x0, dsq_estimate, d2sq) == separate<rt::Halley>(sq, x0, dsq, d2sq));
    }
    CHECK(separate<rt::Newton>(wal, B{ 2.0, 3.0 }, dwal_expected) == separate<rt::Newton>(wal, B{ 2.0, 3.0 }, dwal));

    // Halley with an f'' that always fails takes Newton's step everywhere: Newton's result and count of fn.
    for (double const x0 : { 2.0, 10.0 })
        CHECK(separate<rt::Halley>(wal, x0, dwal, failed) == separate<rt::Newton>(wal, x0, dwal));
    CHECK(separate<rt::Halley>(wal, B{ 2.0, 3.0 }, dwal, failed) == separate<rt::Newton>(wal, B{ 2.0, 3.0 }, dwal));

    // A derivative that fails only below 1.42, against one returning infinity there: bit for bit the same. On [1, 2]
    // Newton's step from 1.4167 would be accepted, and becomes a bisection step; from 10 the iterate 1.4145 has no
    // step. Each differs from the always-analytic run, so the failure is exercised.
    auto const dsq_partial = [](double x) -> std::expected<double, int> {
        if (x < 1.42) return std::unexpected(0);
        return 2.0 * x;
    };
    auto const dsq_inf = [](double x) { return x < 1.42 ? dinf : 2.0 * x; };
    for (auto const &start : { B{ 1.0, 2.0 } }) {
        auto const partial = separate<rt::Newton>(sq, start, dsq_partial);
        CHECK(partial == separate<rt::Newton>(sq, start, dsq_inf));
        CHECK(partial.second != separate<rt::Newton>(sq, start, dsq).second);
        CHECK(brackets(partial.first, sqrt2, 8e-17));
    }
    auto const from_ten = separate<rt::Newton>(sq, 10.0, dsq_partial);
    CHECK(from_ten == separate<rt::Newton>(sq, 10.0, dsq_inf));
    CHECK(from_ten.first == std::unexpected(rt::error::not_converged));

    // The same for f'': below 1.42 Halley takes Newton's step (on [1, 2] from 1.4, where Halley's would be accepted).
    auto const d2_partial = [](double x) -> std::expected<double, int> {
        if (x < 1.42) return std::unexpected(0);
        return 2.0;
    };
    auto const d2_inf = [](double x) { return x < 1.42 ? dinf : 2.0; };
    for (auto const &start : { B{ 1.0, 2.0 }, B{ 1.3, 2.0 } }) {
        auto const partial = separate<rt::Halley>(sq, start, dsq, d2_partial);
        CHECK(partial == separate<rt::Halley>(sq, start, dsq, d2_inf));
        CHECK(partial != separate<rt::Halley>(sq, start, dsq, d2sq));
        CHECK(brackets(partial.first, sqrt2, 8e-17));
    }
    auto const halley_ten = separate<rt::Halley>(sq, 10.0, dsq, d2_partial);
    CHECK(halley_ten == separate<rt::Halley>(sq, 10.0, dsq, d2_inf));
    CHECK(halley_ten != separate<rt::Halley>(sq, 10.0, dsq, d2sq));
}

TEST_CASE("roots:numerical: float, long double and multiprecision types")
{
    // x^2 - 2 with c = eps |ref| (fn's rounding), as in Phase 2; Newton with derivative_of and Halley with two.
    auto const check = [](auto zero) {
        using T = decltype(zero);
        T const eps = std::numeric_limits<T>::epsilon();
        T const tiny = T(std::numeric_limits<T>::denorm_min());
        T const s2 = from_string<T>(sqrt2_70);
        auto const sq_t = [](T const &x) { return T(x * x - 2); };
        auto const d = dv::derivative_of(sq_t);
        auto const d2 = dv::derivative_of<dv::Central<2>, 2>(sq_t);
        for (auto const &r : { rt::solve<rt::Newton>(sq_t, d, T(1)), rt::solve<rt::Halley>(sq_t, d, d2, T(10)),
                 rt::solve<rt::Newton>(sq_t, d, rt::Bracket<T>{ T(1), T(2) }),
                 rt::solve<rt::Halley>(sq_t, d, d2, rt::Bracket<T>{ T(1), T(2) }) })
            CHECK(brackets(r, s2, T(eps * s2), T(2 * T(tiny + T(2 * eps * s2)))));
    };
    check(float{});
    check(static_cast<long double>(0));
    check(bin50{});
    check(bin50_et{});
}
