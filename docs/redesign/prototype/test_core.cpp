// Feasibility prototype for PLAN_v1: 1-D core, combinators, derivative_of, steps_view, FXT at the edge.
// Build with -DNXX_FIX_UNWRAP_FAILURE (see neg/neg_deriv_unwrap.cpp for what happens without it).
#include "nxx/roots.hpp"
#include "nxx/deriv.hpp"
#ifndef NXX_NO_PIPES
#  include "nxx/pipes.hpp"
#endif
#include <cstdio>
#include <limits>
#include <ranges>
#include <string_view>

namespace r = nxx::roots;
namespace d = nxx::deriv;

constexpr auto f  = [](double x) { return x * x - 2.0; };
constexpr auto df = [](double x) { return 2.0 * x; };

// ------------------------------------------------ refined types (tier A/B)
constexpr nxx::tolerance tol_ctad{1e-8};                          // alias-template CTAD
static_assert(std::is_same_v<decltype(tol_ctad), const nxx::tolerance<double>>);
static_assert(nxx::tolerance<double>::make(-1.0).error() == nxx::errc::invalid_input);
static_assert(nxx::tolerance<double>::make(1e-3).value().value() == 1e-3);
static_assert(!std::is_constructible_v<nxx::tolerance<double>, bool>);
static_assert(!std::is_convertible_v<nxx::tolerance<double>, nxx::rel_tolerance<double>>);   // roles not interchangeable
static_assert(std::is_constructible_v<nxx::tolerance<double>, double>);   // TRUE although a run-time double cannot construct it
static_assert(nxx::max_iterations::make(0).error() == nxx::errc::invalid_input);
static_assert(!std::is_constructible_v<nxx::max_iterations, bool>);
static_assert(nxx::bracket<double>::make(2.0, 1.0)->lo() == 1.0);   // make() re-orders
static_assert(!nxx::bracket<double>::make(1.0, 1.0));               // equal endpoints: error
static_assert(!nxx::tolerance<double>::make(std::numeric_limits<double>::quiet_NaN()));   // consteval-safe isfinite
static_assert(!nxx::bracket<double>::make(0.0, std::numeric_limits<double>::infinity()));
static_assert(std::is_trivially_copyable_v<nxx::failure<r::root_estimate<double>>>);
constexpr bool result_trivially_copyable = std::is_trivially_copyable_v<nxx::result<r::root_estimate<double>>>;   // false on libc++
static_assert(std::regular<nxx::tolerance<double>> || !std::default_initializable<nxx::tolerance<double>>);

// ------------------------------------------------ solver concept
using prob_t = nxx::problem<std::reference_wrapper<const decltype(f)>, r::sign_bracket<double>>;
static_assert(nxx::iterative_solver_for<r::bisection<>, prob_t>);
static_assert(nxx::iterative_solver_for<r::brent<>, prob_t>);
static_assert(!std::is_invocable_v<r::newton<>, decltype(f), double>);        // no derivative source: not callable
static_assert(std::is_invocable_v<decltype(r::newton{}.with_derivative(df)), decltype(f), double>);
static_assert(!std::is_invocable_v<r::bisection<>, decltype(f), double>);      // bisection on a guess

// ------------------------------------------------ the headline chain (plan 6.11), evaluated at compile time
constexpr auto chain = nxx::first_of(
    r::newton{}.with_derivative(df).on(0.0),          // fails: f'(0) = 0 -> zero_derivative
    r::secant{nxx::default_step{}, 5}.on(0.0),        // fails: 5 iterations are not enough
    r::bisection{}.on(nxx::bracket{0.0, 2.0}));       // succeeds
static_assert(chain(f).has_value() && chain(f)->by == nxx::algo::bisection);

// search -> coarse bisection -> Newton polish (Kleisli "then")
constexpr auto pipeline = nxx::then(r::expand_out{}.on(nxx::bracket{2.0, 2.5}),
                                    r::bisection{nxx::x_tol{1e-4}},
                                    r::newton{}.with_derivative(df));
static_assert(pipeline(f).has_value() && pipeline(f)->by == nxx::algo::newton);

// Newton with a NUMERIC derivative policy inside a curried chain (f is not known when the chain is built)
constexpr auto chain_fd = nxx::first_of(
    r::newton{}.with_derivative(d::numeric{}).on(1.0),
    r::brent{}.on(nxx::bracket{0.0, 2.0}));
static_assert(chain_fd(f).has_value() && chain_fd(f)->by == nxx::algo::newton);

// brent at compile time
static_assert(r::brent{}(f, nxx::bracket{1.0, 2.0})->used.iterations < 12);

// derivative_of: a function-returning API, constexpr
constexpr auto dsq = d::derivative_of([](double x) { return x * x * x; });
static_assert(dsq(2.0).has_value());
static_assert(nxx::math::abs(*dsq(2.0) - 12.0) < 1e-8);

// criteria algebra
constexpr auto crit = nxx::x_tol{1e-6} || nxx::max_evaluations{20};
constexpr auto crit2 = nxx::default_step{} && nxx::min_iterations{3};

template<class T> constexpr std::size_t type_name_len() {
#if defined(_MSC_VER) && !defined(__clang__)
    return std::string_view{__FUNCSIG__}.size();
#else
    return std::string_view{__PRETTY_FUNCTION__}.size();
#endif
}

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { std::printf("  CHECK FAILED line %d: %s\n", __LINE__, #cond); ++failures; } } while (0)

const char* name(nxx::errc e) {
    switch (e) {
        case nxx::errc::no_sign_change: return "no_sign_change";
        case nxx::errc::budget_exhausted: return "budget_exhausted";
        case nxx::errc::evaluations_exhausted: return "evaluations_exhausted";
        case nxx::errc::stalled: return "stalled";
        case nxx::errc::zero_derivative: return "zero_derivative";
        case nxx::errc::callback_failed: return "callback_failed";
        case nxx::errc::non_finite_value: return "non_finite_value";
        case nxx::errc::singular: return "singular";
        default: return "other";
    }
}

int main(int argc, char**) {
    std::printf("[1] headline chain: ");
    auto res = chain(f);
    CHECK(res && res->by == nxx::algo::bisection);
    if (res) std::printf("x=%.17g iters=%u evals=%u\n", res->x, res->used.iterations, res->used.evaluations);

    std::printf("[2] pipeline expand->bisection->newton: ");
    auto pr = pipeline(f);
    CHECK(pr && nxx::math::abs(pr->x - 1.4142135623730951) < 1e-15);
    if (pr) std::printf("x=%.17g iters=%u evals=%u\n", pr->x, pr->used.iterations, pr->used.evaluations);

    std::printf("[3] newton with numeric derivative policy in a chain: ");
    auto rf = chain_fd(f);
    CHECK(rf && rf->by == nxx::algo::newton);
    if (rf) std::printf("x=%.17g iters=%u evals=%u\n", rf->x, rf->used.iterations, rf->used.evaluations);

    std::printf("[4] brent on x^2-2 over [1,2]: ");
    auto rb = r::brent{}(f, nxx::bracket{1.0, 2.0});
    CHECK(rb && nxx::math::abs(rb->x - 1.4142135623730951) < 4e-16);
    if (rb) std::printf("x=%.17g iters=%u evals=%u how=%d\n", rb->x, rb->used.iterations, rb->used.evaluations, int(rb->how));

    std::printf("[5] no sign change is an error with the best endpoint: ");
    auto ns = r::bisection{}([](double x) { return x * x - 5.0; }, nxx::bracket{3.0, 4.0});
    CHECK(!ns && ns.error().code == nxx::errc::no_sign_change && ns.error().best && ns.error().best->x == 3.0);
    if (!ns) std::printf("%s best x=%g\n", name(ns.error().code), ns.error().best->x);

    std::printf("[6] fallible callback keeps the user's error: ");
    enum class callback_error { diverged };
    auto g = [](double x) -> std::expected<double, callback_error> {
        if (x < 0) return std::unexpected(callback_error::diverged);
        return x * x - 2.0;
    };
    auto bad = r::bisection{}(g, nxx::bracket{-1.0, 2.0});
    CHECK(!bad && bad.error().code == nxx::errc::callback_failed && bad.error().cause && *bad.error().cause == callback_error::diverged);
    auto ok = nxx::first_of(r::bisection{}.on(nxx::bracket{-1.0, 2.0}), r::brent{}.on(nxx::bracket{0.0, 2.0}))(g);
    CHECK(ok && ok->by == nxx::algo::brent);
    std::printf("%s, then first_of -> %s x=%.17g\n", name(bad.error().code), ok ? "brent" : "?", ok ? ok->x : 0.0);

    std::printf("[7] newton on x^2+1 fails with best iterate: ");
    auto nn = r::newton{}.with_derivative(df)([](double x) { return x * x + 1.0; }, 0.5);
    CHECK(!nn && nn.error().code == nxx::errc::budget_exhausted && nn.error().best);
    if (!nn) std::printf("%s best x=%g |f|=%g\n", name(nn.error().code), nn.error().best->x, nn.error().best->fx);

    std::printf("[8] projection: clamped newton from 10 on [1,3]: ");
    auto pc = r::newton{}.with_derivative(df).with_projection(r::clamp_to{1.0, 3.0})(f, 10.0);
    CHECK(pc && nxx::math::abs(pc->x - 1.4142135623730951) < 1e-15);
    if (pc) std::printf("x=%.17g iters=%u\n", pc->x, pc->used.iterations);
    std::printf("    pinned at the edge -> stalled: ");
    auto pin = r::newton{}.with_derivative([](double) { return 1.0; }).with_projection(r::clamp_to{0.0, 3.0})(
        [](double x) { return x - 5.0; }, 1.0);
    CHECK(!pin && pin.error().code == nxx::errc::stalled && pin.error().best->x == 3.0);
    if (!pin) std::printf("%s best x=%g\n", name(pin.error().code), pin.error().best->x);

    std::printf("[9] warm_fallback: starved bisection -> newton from its best: ");
    auto wf = nxx::warm_fallback(r::bisection{nxx::never{}, 3}.on(nxx::bracket{0.0, 2.0}), r::newton{}.with_derivative(df))(f);
    CHECK(wf && wf->by == nxx::algo::newton);
    if (wf) std::printf("x=%.17g total iters=%u evals=%u\n", wf->x, wf->used.iterations, wf->used.evaluations);

    std::printf("[10] criteria algebra: ");
    auto c1 = r::bisection{crit}(f, nxx::bracket{0.0, 2.0});
    auto c2 = r::bisection{nxx::never{} || nxx::max_evaluations{10}}(f, nxx::bracket{0.0, 2.0});
    auto c3 = r::secant{crit2}(f, 1.0);
    CHECK(c1 && !c2 && c2.error().code == nxx::errc::evaluations_exhausted && c3 && c3->used.iterations >= 3);
    std::printf("x_tol||max_evals: %u it x=%.6g; never||max_evals(10): %s after %u evals; step&&min_it(3): %u it\n",
                c1 ? c1->used.iterations : 0, c1 ? c1->x : 0.0, c2 ? "ok" : name(c2.error().code), c2 ? 0 : c2.error().used.evaluations,
                c3 ? c3->used.iterations : 0);

    std::printf("[11] steps_view | take(8) over brent: ");
    auto p = r::brent{}.prepare(std::cref(f), nxx::bracket{1.0, 2.0});
    int n = 0;
    for (const auto& st : nxx::steps_view{r::brent{}, *p} | std::views::take(8)) {
        if (st) std::printf("%.10g ", st->b);
        ++n;
    }
    std::printf("(%d elements, ends at intrinsic stop)\n", n);
    std::printf("     steps_view | transform(estimate) | take(4) over newton: ");
    auto pn = r::newton{}.with_derivative(df).prepare(std::cref(f), 3.0);
    const auto nt = r::newton{}.with_derivative(df);
    for (auto e : nxx::steps_view{nt, *pn} | std::views::transform([&](const auto& st) {
                      return st.transform([&](const auto& s) { return nt.estimate(s); }); }) | std::views::take(4))
        if (e) std::printf("%.12g ", e->x);
    std::printf("\n");
    static_assert(std::ranges::input_range<decltype(nxx::steps_view{nt, *pn})>);
    static_assert(std::ranges::view<decltype(nxx::steps_view{nt, *pn})>);

    std::printf("[12] derivative_of at run time: ");
    auto dsin = d::derivative_of([](double x) { return std::sin(x); });
    auto dsin4 = d::derivative_of([](double x) { return std::sin(x); }, d::central_1_4);
    CHECK(dsin(1.0) && nxx::math::abs(*dsin(1.0) - std::cos(1.0)) < 1e-10);
    std::printf("err(central_1_2)=%.2e err(central_1_4)=%.2e\n", *dsin(1.0) - std::cos(1.0), *dsin4(1.0) - std::cos(1.0));

    std::printf("[13] run-time inputs through make(): ");
    const double lo = argc > 5 ? 9.0 : 0.0;
    auto rb2 = nxx::bracket<double>::make(lo, 2.0).transform_error([](nxx::errc e) {
                   return nxx::failure<r::root_estimate<double>>{e};
               }).and_then([&](nxx::bracket<double> b) { return r::brent{}(f, b); });
    auto budget = nxx::max_iterations::make(argc * 100);
    auto rb3 = r::bisection{nxx::floored_width{}, *budget}(f, *nxx::bracket<double>::make(lo, 2.0));
    CHECK(rb2 && rb3);
    std::printf("brent x=%.17g, bisection(budget=%u) x=%.17g\n", rb2 ? rb2->x : 0.0, budget->value(), rb3 ? rb3->x : 0.0);

#ifndef NXX_NO_PIPES
    std::printf("[14] FXT pipes at the edge: ");
    using fxt::operator|;
    const double x = chain(f) | fxt::transform([](const auto& s) { return s.x; })
                              | fxt::value_or(std::numeric_limits<double>::quiet_NaN());
    const double m = chain(f) | fxt::match([](const auto& s) { return s.x; }, [](const auto&) { return -1.0; });
    auto polished = r::bisection{nxx::x_tol{1e-3}}(f, nxx::bracket{0.0, 2.0})
                  | fxt::and_then([&](const auto& s) { return r::newton{}.with_derivative(df)(f, s); })
                  | fxt::tap([](const auto& s) { std::printf("(tap x=%.6g) ", s.x); });
    CHECK(x == m && polished);
    std::printf("x=%.17g match=%.17g polished=%.17g\n", x, m, polished ? polished->x : 0.0);
#endif

    std::printf("[15] sizes: failure<root_estimate<double>>=%zu result=%zu (result trivially copyable=%d) chain=%zu; type-name length of chain=%zu, pipeline=%zu\n",
                sizeof(nxx::failure<r::root_estimate<double>>), sizeof(nxx::result<r::root_estimate<double>>), int(result_trivially_copyable), sizeof(chain),
                type_name_len<decltype(chain)>(), type_name_len<decltype(pipeline)>());
    std::printf("%s (%d check failures)\n", failures ? "FAIL" : "PASS", failures);
    return failures ? 1 : 0;
}
