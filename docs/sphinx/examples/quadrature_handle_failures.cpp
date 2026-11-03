// The program of docs/sphinx/quadrature/handle-failures.md: what each error means, a failure's best, recovering by
// another method or a looser request, and watching the steps of a run that fails. The page shows the parts between the
// markers and the output in the raw strings; main returns non-zero unless the program prints that output and what the
// page says holds. The integrands use only + - * /, std::sqrt, std::abs and integer operations, which are exact or
// correctly rounded, and the methods only arithmetic, so the output is the same on every target that does not contract
// a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/quadrature.hpp>

#include <cmath>
#include <limits>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <bit>
#include <cstdint>
#include <string>

namespace {

    namespace qd = nxx::quadrature;

    // Noise in [-1, 1] from a hash of x's bits, as in the tests: it never resolves, however close the points.
    double hash(double x)
    {
        std::uint64_t h = std::bit_cast<std::uint64_t>(x) * 0x9E3779B97F4A7C15ull;
        h ^= h >> 29;
        h *= 0xBF58476D1CE4E5B9ull;
        h ^= h >> 32;
        return static_cast<double>(h >> 11) / 4503599627370496.0 - 1.0;
    }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <reasons>
    auto const show = [](char const *what, auto const &r) {
        if (r) std::println("{:<38} ok", what);
        else std::println("{:<38} {}", what, r.error());
    };
    double const inf = std::numeric_limits<double>::infinity();
    auto const inverse = [](double x) { return 1.0 / x; };
    show("the default on [1, inf)", qd::integrate(inverse, 1.0, inf));
    show("a negative rel_tol", qd::integrate(inverse, 1.0, 2.0, {.rel_tol = -1.0}));
    show("sqrt(x - 0.5) on [0, 1]", qd::integrate([](double x) { return std::sqrt(x - 0.5); }, 0.0, 1.0));
    show("Romberg on 1/x over [0, 1]", qd::integrate<qd::Romberg>(inverse, 0.0, 1.0));
    show("the default on 1/x over [0, 1]", qd::integrate(inverse, 0.0, 1.0));
    // </reasons>

    // <isolated>
    double const nan = std::numeric_limits<double>::quiet_NaN();
    int evaluations = 0;
    auto const spike = [&evaluations, inf](double x) {
        ++evaluations;
        return x == 0.35 ? inf : 1.0; // 0.35 is one of the default's points
    };
    auto const s = qd::integrate(spike, 0.0, 1.0);
    if (s) std::println("inf at 0.35: {:.17g} +- {:.1e}, {} evaluations", s->value, s->abs_error, evaluations);
    show("NaN on (0.3, 0.4)", qd::integrate([nan](double x) { return x > 0.3 && x < 0.4 ? nan : 1.0; }, 0.0, 1.0));
    // </isolated>

    // <best>
    auto const noisy = [](double x) { return 1.0 / (1.0 + x) * (1.0 + 1e-10 * hash(x)); };
    auto const n = qd::integrate(noisy, 0.0, 1.0, {.rel_tol = 1e-12});
    if (!n && n.error().best)
        std::println("noisy: {}, best {:.15f} +- {:.1e}", n.error(), n.error().best->value, n.error().best->abs_error);
    auto const d = qd::integrate(inverse, 0.0, 1.0);
    if (!d && d.error().best)
        std::println("1/x:   {}, best {:.15f} +- {:.1e}", d.error(), d.error().best->value, d.error().best->abs_error);
    // </best>

    // <fallback>
    auto const g = [](double x) { return 1.0 / std::sqrt(x) + std::abs(x - 1.0 / 3.0); };
    double const a = 0.0;
    double const b = 1.0;
    auto const t = qd::integrate<qd::TanhSinh>(g, a, b);
    auto const r = t.or_else([&](qd::Failure<double> const &f) -> decltype(t) {
        if (f.reason != qd::error::not_converged || !f.best) return t; // non_finite, or no estimate at all
        auto const k = qd::integrate(g, a, b);                          // interior features: the default
        if (k && std::abs(k->value - f.best->value) <= f.best->abs_error + k->abs_error) return k;
        return t;
    });
    if (!t) std::println("TanhSinh: {}", t.error());
    if (r) std::println("fallback: {:.15f} +- {:.1e}", r->value, r->abs_error);
    // </fallback>

    // <looser>
    auto const looser = qd::integrate(noisy, 0.0, 1.0, {.rel_tol = 1e-10});
    if (looser) std::println("{:.15f} +- {:.1e}", looser->value, looser->abs_error);
    // </looser>

    // <diverge>
    auto const loose = qd::integrate(inverse, 0.0, 1.0, {.rel_tol = 0.2});
    if (loose) std::println("1/x on [0, 1], rel_tol 0.2: ok, {:.4f} +- {:.1e}", loose->value, loose->abs_error);
    for (double const lower : {1e-2, 1e-4, 1e-8, 1e-16}) {
        auto const moved = qd::integrate(inverse, lower, 1.0);
        if (moved) std::println("1/x on [{:.0e}, 1]: {:.4f}", lower, moved->value);
    }
    // </diverge>

    // <observe>
    int steps = 0;
    evaluations = 0;
    auto const watch = [&steps](qd::Step<double> const &step) {
        if (step.current) std::print("step {}: {:.6f} +- {:.1e}", steps, step.current->value, step.current->abs_error);
        else std::print("step {}: no estimate", steps);
        std::println("{}{}", step.ready ? ", ready" : "", step.exhausted ? ", exhausted" : "");
        ++steps;
    };
    auto const counted = [&evaluations](double x) {
        ++evaluations;
        return 1.0 / x;
    };
    auto const w = qd::integrate<qd::TanhSinh>(qd::observer(watch), counted, 0.0, 1.0);
    if (!w) std::println("{} after {} steps and {} evaluations", w.error(), steps, evaluations);
    // </observe>
    // clang-format on

    constexpr std::string_view reasons_output = R"(
the default on [1, inf)                non_finite
a negative rel_tol                     non_finite
sqrt(x - 0.5) on [0, 1]                non_finite
Romberg on 1/x over [0, 1]             non_finite
the default on 1/x over [0, 1]         not_converged
)";
    constexpr std::string_view isolated_output = R"(
inf at 0.35: 1 +- 1.1e-14, 180 evaluations
NaN on (0.3, 0.4)                      non_finite
)";
    constexpr std::string_view best_output = R"(
noisy: not_converged, best 0.693147180561986 +- 3.4e-12
1/x:   not_converged, best 31.515095842644442 +- 1.8e+00
)";
    constexpr std::string_view fallback_output = R"(
TanhSinh: not_converged
fallback: 2.277777770619851 +- 3.1e-08
)";
    constexpr std::string_view looser_output = R"(
0.693147180560495 +- 1.1e-11
)";
    constexpr std::string_view diverge_output = R"(
1/x on [0, 1], rel_tol 0.2: ok, 9.3344 +- 1.8e+00
1/x on [1e-02, 1]: 4.6052
1/x on [1e-04, 1]: 9.2103
1/x on [1e-08, 1]: 18.4207
1/x on [1e-16, 1]: 36.8414
)";
    constexpr std::string_view observe_output = R"(
step 0: no estimate
step 1: no estimate
step 2: no estimate
step 3: no estimate
step 4: no estimate, ready
step 5: no estimate, ready
step 6: no estimate, ready
step 7: no estimate, ready
step 8: no estimate, ready
step 9: no estimate, ready, exhausted
not_converged after 10 steps and 5159 evaluations
)";
    std::string const expected = std::string(docs::quoted(reasons_output)) + std::string(docs::quoted(isolated_output))
                                 + std::string(docs::quoted(best_output)) + std::string(docs::quoted(fallback_output))
                                 + std::string(docs::quoted(looser_output)) + std::string(docs::quoted(diverge_output))
                                 + std::string(docs::quoted(observe_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    double const ln2 = 0.69314718055994530941723212145817656807550013436026;
    // The noisy run's best: "here within its abs_error alone" (2.0e-12 against 3.4e-12); "what you can count on is
    // abs_error plus the noise's bound, 1e-10 times the integral of |f|, about 7e-11 in all".
    if (!(!n && n.error().best && std::abs(n.error().best->value - ln2) <= n.error().best->abs_error)) return 2;
    if (!(n.error().best->abs_error + 1e-10 * ln2 < 7.5e-11 && n.error().best->abs_error + 1e-10 * ln2 > 6.5e-11))
        return 12;
    // "TanhSinh gives no best on 1/x" (the observed run too).
    auto const ts = qd::integrate<qd::TanhSinh>(inverse, 0.0, 1.0);
    if (!(!ts && !ts.error().best && !w && !w.error().best)) return 3;
    // "The same best whatever the request": the default on 1/x at 1e-4 and 1e-12 gives the same failure.
    if (!(qd::integrate(inverse, 0.0, 1.0, {.rel_tol = 1e-4}) == d
          && qd::integrate(inverse, 0.0, 1.0, {.rel_tol = 1e-12}) == d))
        return 4;
    // The fallback: TanhSinh failed with a best; the default's result was accepted; it is within its abs_error of the
    // exact 2 + 5/18 (the integral of 1 / sqrt x is 2, of |x - 1/3| over [0, 1] (1/3)^2 / 2 + (2/3)^2 / 2 = 5/18).
    double const exact = 2.0 + 5.0 / 18.0;
    if (!(!t && t.error().best && r && std::abs(r->value - exact) <= r->abs_error)) return 5;
    // "The default alone would have succeeded too": the fallback's result is the default's.
    if (!(qd::integrate(g, a, b) == r)) return 6;
    // The looser request on the noisy function: ok, and within its abs_error of ln 2.
    if (!(looser && std::abs(looser->value - ln2) <= looser->abs_error)) return 7;
    // "32 bisections in a row fail to improve its estimate, which is what happened on 1/x": the mesh's 150 evaluations
    // and 32 bisections of 30.
    int stagnant = 0;
    auto const s32 = qd::integrate([&stagnant](double x) { ++stagnant; return 1.0 / x; }, 0.0, 1.0);
    if (!(s32 == d && stagnant == 150 + 32 * 30)) return 8;
    // "At a strong singularity the default can succeed falsely": x^-0.875 on [0, 1], whose integral is 8, ends ok with
    // an abs_error (1.2e-7) below its actual error (4.5e-7); x^-0.875 = 1 / (x^(1/2) x^(1/4) x^(1/8)) by square roots.
    auto const strong = qd::integrate([](double x) {
        double const r2 = std::sqrt(x);
        double const r4 = std::sqrt(r2);
        return 1.0 / (r2 * r4 * std::sqrt(r4));
    }, 0.0, 1.0);
    if (!(strong && std::abs(strong->value - 8.0) > strong->abs_error)) return 9;
    // "The observer changes nothing": the observed run's failure is the unobserved one's.
    if (!(w == ts)) return 10;
    // "On an integral that converges, the steps show the estimate falling level by level": 1 / sqrt x on [0, 1], the
    // last estimate at least a million times below the first.
    double first = -1.0;
    double last = -1.0;
    auto const falling = qd::integrate<qd::TanhSinh>(qd::observer([&](qd::Step<double> const &step) {
        if (!step.current) return;
        if (first < 0.0) first = step.current->abs_error;
        last = step.current->abs_error;
    }), [](double x) { return 1.0 / std::sqrt(x); }, 0.0, 1.0);
    if (!(falling && first > 0.0 && last * 1e6 < first)) return 11;
    // The isolated value: 0.35 is one of the default's points, and one bisection, 30 evaluations beyond the mesh's 150
    // (the output), stepped around it; the result is exact (the integral of 1). "A value no point lands on is never
    // seen": TanhSinh does not evaluate at 0.35.
    if (!(s && s->value == 1.0)) return 13;
    int at = 0;
    auto const ts_spike = qd::integrate<qd::TanhSinh>([&at](double x) { at += x == 0.35; return 1.0; }, 0.0, 1.0);
    if (!(ts_spike && at == 0)) return 14;
    // The divergent integral: "at rel_tol 0.2 the default reports ok, about 9.3, after its first 150 evaluations";
    // "the integral over [10^-k, 1] grows like k ln 10" (exactly k ln 10 = ln(1 / lower)): within the result's own
    // estimate plus the reference's rounding, std::log within an ulp and the negation exact, 2 eps relative (the
    // independent review: the first version allowed an unjustified 1e-6).
    int loose_count = 0;
    auto const loose_counted = qd::integrate([&loose_count](double x) { ++loose_count; return 1.0 / x; }, 0.0, 1.0,
                                             {.rel_tol = 0.2});
    if (!(loose_counted == loose && loose_count == 150)) return 15;
    for (double const lower : {1e-2, 1e-4, 1e-8, 1e-16}) {
        auto const moved = qd::integrate(inverse, lower, 1.0);
        double const log_lower = -std::log(lower);
        double const ulp = std::numeric_limits<double>::epsilon();
        if (!(moved && std::abs(moved->value - log_lower) <= moved->abs_error + 2.0 * ulp * log_lower)) return 16;
    }
    return 0;
}
