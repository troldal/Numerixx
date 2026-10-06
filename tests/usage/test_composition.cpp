// Composition across modules (DESIGN §10.2 exit criterion 11, §6.4, §6.12, §6.10): Newton with derivative_of(g) of a
// fallible g keeps g's error type and the derivative's cause and cost; the numeric derivative policy works in a
// curried chain; evaluation counts equal instrumented calls of f; cost-aware criteria see the derivative's evaluations.
//
// Mixing f and f' with different (non-none) error types is a compile error that names .transform_error; that case is
// covered by the compile-fail tests (DESIGN §9.1), not here.

#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <expected>
#include <functional>
#include <limits>
#include <optional>
#include <random>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

namespace r = nxx::roots;
namespace d = nxx::deriv;

namespace
{
    double rt(double v)
    {
        volatile double x = v;
        return x;
    }

    const double sqrt2 = std::sqrt(2.0);

    enum class eval_error { domain, too_big };

    // x^2 - 2, defined for x >= 0 only.
    constexpr auto g_fallible = [](double x) -> std::expected<double, eval_error> {
        if (x < 0.0) return std::unexpected(eval_error::domain);
        return x * x - 2.0;
    };

    constexpr auto f_plain = [](double x) { return x * x - 2.0; };

    template<class R>
    using cause_of = typename std::remove_cvref_t<R>::error_type::cause_type;

    // log(x) - log(c) with c = 1.79769e308, within 1.7e-6 of the largest double: near the root, a stencil point of the
    // numeric derivative overflows, and diff rejects it with invalid_input (DESIGN §7.1).
    constexpr double near_max     = 1.79769e308;
    const auto       log_near_max = [](double x) { return std::log(x) - std::log(near_max); };

    // A chain policy that continues after numerical errors and stops on input errors.
    struct stop_on_input_error
    {
        template<class Est, class UE>
        constexpr bool operator()(const nxx::failure<Est, UE>& err) const noexcept
        { return !nxx::is_input_error(err.code); }
    };
}    // namespace

TEST_SUITE("usage")
{
    TEST_CASE("newton with derivative_of(g) of a fallible g: g's error type, and success")
    {
        const auto res = r::newton {}.with_derivative(d::derivative_of(g_fallible))(g_fallible, rt(1.0));
        static_assert(std::is_same_v<cause_of<decltype(res)>, eval_error>);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(res)>, nxx::result<r::root_estimate<double>, eval_error>>);
        static_assert(std::is_same_v<nxx::callback_error_t<decltype(d::derivative_of(g_fallible)), double>, eval_error>);
        CHECK(res.has_value());
        if (res) {
            CHECK(res->by == r::algos::newton);
            CHECK(std::abs(res->x - sqrt2) < 1e-12);
            CHECK(res->used.evaluations == 1 + 3 * res->used.iterations);    // f(x0), then per step f' (2 calls) and f
        }

        // The common cause of f and f' is the one that is not none (DESIGN §6.4).
        using plain_f_fallible_df = decltype(r::newton {}.with_derivative(d::derivative_of(g_fallible))(f_plain, 1.0));
        static_assert(std::is_same_v<cause_of<plain_f_fallible_df>, eval_error>);
        using fallible_f_plain_df = decltype(r::newton {}.with_derivative(d::derivative_of(f_plain))(g_fallible, 1.0));
        static_assert(std::is_same_v<cause_of<fallible_f_plain_df>, eval_error>);
        using both_plain = decltype(r::newton {}.with_derivative(d::derivative_of(f_plain))(f_plain, 1.0));
        static_assert(std::is_same_v<cause_of<both_plain>, nxx::none>);
    }

    TEST_CASE("g fails inside the derivative's evaluation: callback_failed, g's cause, the cost counted")
    {
        std::uint32_t calls = 0;
        const auto    g     = nxx::fn::counted(
            [](double x) -> std::expected<double, eval_error> {
                if (x < 0.0) return std::unexpected(eval_error::domain);
                if (x > 1.0) return std::unexpected(eval_error::too_big);
                return x * x - 0.5;
            },
            calls);

        SUBCASE("the first stencil point")
        {
            // f(0.25) is fine; the derivative's first point, 0.25 - 0.5, is not.
            const auto res = r::newton {}.with_derivative(d::derivative_of(g, d::central_1_2, d::absolute { 0.5 }))(g, rt(0.25));
            static_assert(std::is_same_v<cause_of<decltype(res)>, eval_error>);
            CHECK_FALSE(res.has_value());
            if (!res) {
                CHECK(res.error().code == nxx::errc::callback_failed);
                CHECK(res.error().by == r::algos::newton);
                CHECK(res.error().cause == eval_error::domain);
                CHECK(res.error().used == nxx::counters { 1, 2 });
                CHECK(res.error().used.evaluations == calls);
                CHECK(nxx::best_x(res) == std::optional<double>(0.25));
            }
        }
        SUBCASE("the second stencil point")
        {
            // f(0.75) and f(0.25) are fine; f(1.25) is not: the failing derivative cost 2 calls.
            const auto res = r::newton {}.with_derivative(d::derivative_of(g, d::central_1_2, d::absolute { 0.5 }))(g, rt(0.75));
            CHECK_FALSE(res.has_value());
            if (!res) {
                CHECK(res.error().code == nxx::errc::callback_failed);
                CHECK(res.error().cause == eval_error::too_big);
                CHECK(res.error().used == nxx::counters { 1, 3 });
                CHECK(res.error().used.evaluations == calls);
            }
        }
    }

    TEST_CASE("g fails inside the derivative after successful steps: every evaluation is counted")
    {
        std::uint32_t calls = 0;
        const auto    g     = nxx::fn::counted(
            [](double x) -> std::expected<double, eval_error> {
                if (x < 0.0) return std::unexpected(eval_error::domain);
                return x * x - 0.5;
            },
            calls);
        // Central differences are exact on a quadratic: 3 -> 1.5833 -> 0.9496, whose stencil point 0.9496 - 1 < 0.
        const auto res = r::newton {}.with_derivative(d::derivative_of(g, d::central_1_2, d::absolute { 1.0 }))(g, rt(3.0));
        CHECK_FALSE(res.has_value());
        if (!res) {
            CHECK(res.error().code == nxx::errc::callback_failed);
            CHECK(res.error().cause == eval_error::domain);
            CHECK(res.error().used == nxx::counters { 3, 8 });    // 1 + 3 + 3 + 1
            CHECK(res.error().used.evaluations == calls);
            CHECK(std::abs(nxx::best_x(res).value_or(0.0) - 0.9495614) < 1e-6);    // the iterate with the smallest |f|
        }
    }

    TEST_CASE("the cause survives a chain: first_of over newton with derivative_of(g) and bisection")
    {
        const auto res =
            nxx::first_of(r::newton {}.with_derivative(d::derivative_of(g_fallible, d::central_1_2, d::absolute { 0.5 })).on(0.25),
                          r::bisection {}.on(nxx::bracket { -1.0, 2.0 }))(g_fallible);
        static_assert(std::is_same_v<cause_of<decltype(res)>, eval_error>);
        CHECK_FALSE(res.has_value());
        if (!res) {
            CHECK(res.error().code == nxx::errc::callback_failed);    // bisection's g(-1)
            CHECK(res.error().by == r::algos::bisection);
            CHECK(res.error().cause == eval_error::domain);
            CHECK(res.error().used == nxx::counters { 1, 3 });    // newton: 1 iteration, 2 calls; bisection: 1 call
        }
    }

    TEST_CASE("the numeric policy in a curried chain: the chain is built before f is known")
    {
        const auto chain =
            nxx::first_of(r::newton {}.with_derivative(nxx::deriv::numeric {}).on(1.0), r::brent {}.on(nxx::bracket { 0.0, 2.0 }));

        const auto res = chain(f_plain);
        CHECK(res.has_value());
        if (res) {
            CHECK(res->by == r::algos::newton);
            CHECK(std::abs(res->x - sqrt2) < 1e-12);
            CHECK(res->used.evaluations == 1 + 3 * res->used.iterations);
        }

        // The same chain value with another function type: an instrumented f ...
        std::uint32_t calls = 0;
        const auto    fc    = nxx::fn::counted(f_plain, calls);
        const auto    rc    = chain(fc);
        CHECK(rc.has_value());
        if (rc && res) {
            CHECK(rc->used.evaluations == calls);
            CHECK(rc->x == res->x);
        }

        // ... and a fallible g: the policy binds whichever function the chain is given.
        const auto rg = chain(g_fallible);
        static_assert(std::is_same_v<cause_of<decltype(rg)>, eval_error>);
        CHECK(rg.has_value());
        if (rg) {
            CHECK(rg->by == r::algos::newton);
            CHECK(std::abs(rg->x - sqrt2) < 1e-12);
        }

        // The whole chain runs at compile time.
        constexpr auto cchain =
            nxx::first_of(r::newton {}.with_derivative(nxx::deriv::numeric {}).on(1.0), r::brent {}.on(nxx::bracket { 0.0, 2.0 }));
        static_assert(cchain(f_plain).has_value() && cchain(f_plain)->by == r::algos::newton);
    }

    TEST_CASE("the numeric policy falls through: newton fails at a zero numeric derivative, brent succeeds")
    {
        const auto chain =
            nxx::first_of(r::newton {}.with_derivative(nxx::deriv::numeric {}).on(0.0), r::brent {}.on(nxx::bracket { 0.0, 2.0 }));
        std::uint32_t calls = 0;
        const auto    fc    = nxx::fn::counted(f_plain, calls);
        const auto    res   = chain(fc);
        CHECK(res.has_value());
        if (res) {
            CHECK(res->by == r::algos::brent);
            CHECK(std::abs(res->x - sqrt2) < 1e-14);
            CHECK(res->used.evaluations == calls);
        }
        const auto nt = r::newton {}.with_derivative(nxx::deriv::numeric {})(f_plain, 0.0);
        CHECK_FALSE(nt.has_value());
        if (!nt) {
            CHECK(nt.error().code == nxx::errc::zero_derivative);
            CHECK(nt.error().used == nxx::counters { 1, 3 });    // f(0), then the two points of f'(0)
        }
        const auto br = r::brent {}(f_plain, nxx::bracket { 0.0, 2.0 });
        if (res && !nt && br) {
            CHECK(res->used == nt.error().used + br->used);    // the success pays for the failed attempt
        }
    }

    TEST_CASE("an input code from a nested derivative inside a step is non_finite_value (regression)")
    {
        // DESIGN §6.3, §6.7: a step never reports invalid_input or non_finite_input. In the 5th step, Newton's numeric
        // derivative rejects its stencil at x4 (x4 + h overflows), and the driver reported that as invalid_input, an
        // input error, although the input 1e308 was valid. The counts and the best estimate are unchanged by the mapping.
        std::uint32_t calls = 0;
        const auto    fc    = nxx::fn::counted(log_near_max, calls);
        const auto    nt    = r::newton {}.with_derivative(d::numeric {});
        const auto    res   = nt(fc, rt(1e308));
        CHECK_FALSE(res.has_value());
        if (!res) {
            CHECK(res.error().code == nxx::errc::non_finite_value);
            CHECK_FALSE(nxx::is_input_error(res.error().code));
            CHECK(res.error().by == r::algos::newton);
            CHECK(res.error().used == nxx::counters { 5, 13 });    // f(x0), f' (2 calls) and f in 4 steps; the 5th calls nothing
            CHECK(res.error().used.evaluations == calls);
            CHECK(res.error().best.has_value());
            if (res.error().best) { CHECK(std::abs(res.error().best->x / near_max - 1.0) < 1e-8); }
        }

        // The same with derivative_of(g) as the secant's f: f itself fails mid-iteration.
        // g(x) = x/1e300 (log(x) - log(c) - 1), so g'(x) = (log(x) - log(c))/1e300.
        std::uint32_t gcalls = 0;
        const auto    g      = nxx::fn::counted([](double x) { return x / 1e300 * (std::log(x) - std::log(near_max) - 1.0); }, gcalls);
        const auto    sec    = r::secant {}(d::derivative_of(g), rt(1e308));
        CHECK_FALSE(sec.has_value());
        if (!sec) {
            CHECK(sec.error().code == nxx::errc::non_finite_value);
            CHECK(sec.error().used.iterations > 0u);
            CHECK(sec.error().used.evaluations == gcalls);
            CHECK(sec.error().best.has_value());
        }
    }

    TEST_CASE("first_of_with stopping on input errors falls through a nested input code to brent")
    {
        // Before the mapping, newton's failure was invalid_input, so this policy stopped the chain there.
        const auto newton_first = r::newton {}.with_derivative(d::numeric {}).on(rt(1e308));
        const auto chain =
            nxx::first_of_with(stop_on_input_error {}, newton_first, r::brent {}.on({ rt(1e308), rt(1.7976931348623157e308) }));
        std::uint32_t calls = 0;
        const auto    fc    = nxx::fn::counted(log_near_max, calls);
        const auto    res   = chain(fc);
        CHECK(res.has_value());
        if (res) {
            CHECK(res->by == r::algos::brent);
            CHECK(std::abs(res->x / near_max - 1.0) < 1e-12);
            CHECK(res->used.evaluations == calls);
        }
        // A genuine input error still stops the same chain.
        const auto stops = nxx::first_of_with(stop_on_input_error {},
                                              r::bisection {}.on({ rt(1.0), rt(1.0) }),
                                              r::brent {}.on({ rt(1e308), rt(1.7976931348623157e308) }));
        const auto res2  = stops(log_near_max);
        CHECK_FALSE(res2.has_value());
        if (!res2) { CHECK(res2.error().code == nxx::errc::invalid_input); }
    }

    TEST_CASE("steps_view maps a nested input code as the driver does")
    {
        // nxx::evaluate maps the numeric derivative's input code (DESIGN §6.4, §12 item 22), below both the view and the
        // driver: the view's last element is the driver's fault, non_finite_value, and the evaluations add up to the
        // driver's.
        const auto nt = r::newton {}.with_derivative(d::numeric {});
        const auto p  = nt.prepare(std::cref(log_near_max), rt(1e308));
        if (!p) {
            FAIL_CHECK("prepare only checks that the guess is finite");
            return;
        }
        std::vector<std::remove_cvref_t<decltype(*nxx::steps_view { nt, *p }.begin())>> elems;
        for (const auto& st : nxx::steps_view { nt, *p } | std::views::take(20)) elems.push_back(st);
        const auto driver = nt(log_near_max, rt(1e308));
        CHECK(elems.size() == 6u);    // init, 4 steps, the failing step
        if (elems.size() == 6 && !driver && elems[4] && !elems[5]) {
            CHECK(elems[5].error().code == nxx::errc::non_finite_value);
            CHECK(elems[5].error().code == driver.error().code);
            CHECK(elems[4]->nfev + elems[5].error().evaluations == driver.error().used.evaluations);
            CHECK(driver.error().best.has_value());
            if (driver.error().best) { CHECK(elems[4]->x == driver.error().best->x); }
        }
        else {
            FAIL_CHECK("the view and the driver disagree on where newton fails");
        }
    }

    TEST_CASE("an input code from a nested callable is non_finite_value at every evaluation site (DESIGN §12 item 22)")
    {
        // derivative_of(g) rejects a stencil that overflows (x + h at +-max) with invalid_input and no call of g. Inside
        // nxx::evaluate that fault becomes non_finite_value with its own evaluations; before the rule, prepare and init
        // passed invalid_input through, which looked like the solver's own rejection of its input. g'(x) is smooth and
        // bounded, with its root at 0.5.
        constexpr double mx     = (std::numeric_limits<double>::max)();
        std::uint32_t    gcalls = 0;
        const auto       g      = nxx::fn::counted([](double x) { return std::hypot(1.0, x - 0.5); }, gcalls);
        const auto       dg     = d::derivative_of(g);
        using est_t             = r::root_estimate<double>;

        // Checks the code, the cost and whether there is a best estimate (and its x).
        const auto check_failure = [&](const auto& res, nxx::counters used, std::optional<double> best_x) {
            CHECK_FALSE(res.has_value());
            if (res) return;
            CHECK(res.error().code == nxx::errc::non_finite_value);
            CHECK(res.error().used == used);
            CHECK(res.error().used.evaluations == gcalls);
            CHECK(res.error().best.has_value() == best_x.has_value());
            if (res.error().best && best_x) { CHECK(res.error().best->x == *best_x); }
        };

        // prepare (bracketing): the first sample, then the second (the first costs 2 calls of g).
        gcalls = 0;
        check_failure(r::brent {}(dg, { rt(-mx), rt(1.0) }), { 0, 0 }, std::nullopt);
        gcalls = 0;
        check_failure(r::brent {}(dg, { rt(1e308), rt(mx) }), { 0, 2 }, 1e308);
        gcalls = 0;
        check_failure(r::solve(dg, { rt(-mx), rt(1.0) }), { 0, 0 }, std::nullopt);

        // init: expand's first and second samples; the secant's f(x0), its f(x1), and f(x1) after a seed; newton's f(x0).
        gcalls = 0;
        check_failure(r::expand {}(dg, { rt(-mx), rt(1.0) }), { 0, 0 }, std::nullopt);
        gcalls = 0;
        check_failure(r::expand {}(dg, { rt(0.0), rt(mx) }), { 0, 2 }, 0.0);
        gcalls = 0;
        check_failure(r::secant {}(dg, rt(mx)), { 0, 0 }, std::nullopt);
        gcalls = 0;
        check_failure(r::secant {}(dg, rt(1.79593e308)), { 0, 2 }, 1.79593e308);    // x1 = x0 + h overflows the stencil
        gcalls = 0;
        check_failure(r::secant {}(dg, est_t { rt(1.79593e308), 1.0 }), { 0, 0 }, 1.79593e308);
        gcalls = 0;
        check_failure(r::newton {}.with_derivative(d::numeric {})(dg, rt(mx)), { 0, 0 }, std::nullopt);

        // steps_view element 0 is init, as the driver sees it.
        {
            const auto ex = r::expand {};
            const auto p  = ex.prepare(std::cref(dg), std::pair { rt(-mx), rt(1.0) });
            CHECK(p.has_value());
            if (p) {
                const auto first = *nxx::steps_view { ex, *p }.begin();
                CHECK_FALSE(first.has_value());
                if (!first) {
                    CHECK(first.error().code == nxx::errc::non_finite_value);
                    CHECK(first.error().evaluations == 0u);
                }
            }
        }

        // nxx::evaluate and diff over a nested derivative. diff's x is valid; its stencil reaches derivative_of(g)'s
        // overflow at the upper point, after 2 calls of g at the lower one.
        const auto ev = nxx::evaluate(dg, rt(mx));
        CHECK(ev == std::expected<double, nxx::fault<>> { std::unexpect, nxx::fault<> { nxx::errc::non_finite_value, 0, {} } });
        gcalls        = 0;
        const auto dd = d::diff(dg, rt(mx / ((1.0 + std::ldexp(1.0, -17)) * (1.0 + std::ldexp(1.0, -18)))));
        CHECK(dd == std::expected<double, nxx::fault<>> { std::unexpect, nxx::fault<> { nxx::errc::non_finite_value, 2, {} } });
        CHECK(gcalls == 2u);

        // Direct calls keep their codes: the input is the caller's own.
        CHECK(dg(rt(mx)) == std::expected<double, nxx::fault<>> { std::unexpect, nxx::fault<> { nxx::errc::invalid_input, 0, {} } });
        CHECK(d::diff(g, rt(std::numeric_limits<double>::quiet_NaN())) ==
              std::expected<double, nxx::fault<>> { std::unexpect, nxx::fault<> { nxx::errc::non_finite_input, 0, {} } });
        CHECK(d::diff(dg, rt(std::numeric_limits<double>::quiet_NaN())) ==
              std::expected<double, nxx::fault<>> { std::unexpect, nxx::fault<> { nxx::errc::non_finite_input, 0, {} } });

        // A fallible callback that returns a Numerixx fault itself: the evaluations and the cause are kept. Input codes
        // other than invalid_input and non_finite_input pass through (DESIGN §6.3).
        using F = nxx::fault<eval_error>;
        for (const auto& [given, mapped] : { std::pair { nxx::errc::invalid_input, nxx::errc::non_finite_value },
                                             std::pair { nxx::errc::non_finite_input, nxx::errc::non_finite_value },
                                             std::pair { nxx::errc::no_sign_change, nxx::errc::no_sign_change } })
        {
            const auto fallible = [given](double) -> std::expected<double, F> {
                return std::unexpected(F { given, 3, eval_error::too_big });
            };
            CHECK(nxx::evaluate(fallible, 1.0) == std::expected<double, F> { std::unexpect, F { mapped, 3, eval_error::too_big } });
            CHECK(nxx::evaluate_sample(fallible, 1.0) == std::expected<double, F> { std::unexpect, F { mapped, 3, eval_error::too_big } });
        }

        // A chain that stops on input errors no longer depends on the order of its alternatives: brent's failure at
        // -max is numerical, so the chain falls through to bisection either way.
        const auto brent_wide = r::brent {}.on({ -mx, 1.0 });
        const auto bisect01   = r::bisection {}.on({ 0.0, 1.0 });
        const auto forward    = nxx::first_of_with(stop_on_input_error {}, brent_wide, bisect01)(dg);
        const auto reverse    = nxx::first_of_with(stop_on_input_error {}, bisect01, brent_wide)(dg);
        CHECK(forward.has_value());
        CHECK(reverse.has_value());
        if (forward && reverse) {
            CHECK(forward->by == r::algos::bisection);
            CHECK(reverse->by == r::algos::bisection);
            CHECK(forward->x == reverse->x);
            CHECK(forward->used == reverse->used);    // brent's attempt cost nothing: the stencil failed before calling g
        }
    }

    TEST_CASE("cost_of(derivative_of(f)): the stencil's points")
    {
        static_assert(nxx::cost_of(d::derivative_of(f_plain)) == 2);
        static_assert(nxx::cost_of(d::derivative_of(f_plain, d::central_1_4)) == 4);
        static_assert(nxx::cost_of(d::numeric {}.bind(f_plain)) == 2);
        static_assert(nxx::cost_of(d::numeric { d::central_1_4 }.bind(f_plain)) == 4);
        CHECK(nxx::cost_of(d::derivative_of(g_fallible)) == 2);
        CHECK(nxx::cost_of(d::derivative_of(g_fallible, d::central_1_4)) == 4);
    }

    TEST_CASE("evaluation counts equal instrumented calls: newton with the numeric policy (property)")
    {
        std::mt19937                           rng(20260928u);
        std::uniform_real_distribution<double> guess(-10.0, 10.0);
        int                                    succeeded = 0;
        for (int i = 0; i < 200; ++i) {
            const double x0 = guess(rng);

            std::uint32_t calls2 = 0;
            const auto    r2     = r::newton {}.with_derivative(d::numeric {})(nxx::fn::counted(f_plain, calls2), x0);
            const auto    used2  = r2 ? r2->used.evaluations : r2.error().used.evaluations;
            CHECK(used2 == calls2);

            // A fallible g whose domain the guess may lie outside: failures count their evaluations too.
            std::uint32_t callsg = 0;
            const auto    rg     = r::newton {}.with_derivative(d::numeric {})(nxx::fn::counted(g_fallible, callsg), x0);
            const auto    usedg  = rg ? rg->used.evaluations : rg.error().used.evaluations;
            CHECK(usedg == callsg);

            // The chain: counts include the failed attempts.
            std::uint32_t callsc = 0;
            const auto    rc     = nxx::first_of(r::newton {}.with_derivative(d::numeric {}).on(x0),
                                                 r::brent {}.on(nxx::bracket { 0.0, 2.0 }))(nxx::fn::counted(g_fallible, callsc));
            const auto    usedc  = rc ? rc->used.evaluations : rc.error().used.evaluations;
            CHECK(usedc == callsc);

            if (r2 && rc) ++succeeded;
        }
        CHECK(succeeded > 150);
    }

    TEST_CASE("evaluation counts equal instrumented calls: the numeric policy with a stencil and a step (property)")
    {
        std::mt19937                           rng(20260929u);
        std::uniform_real_distribution<double> guess(-10.0, 10.0);
        int                                    succeeded = 0;
        for (int i = 0; i < 200; ++i) {
            const double x0 = guess(rng);

            std::uint32_t calls4 = 0;
            const auto    r4     = r::newton {}.with_derivative(d::numeric { d::central_1_4 })(nxx::fn::counted(f_plain, calls4), x0);
            const auto    used4  = r4 ? r4->used.evaluations : r4.error().used.evaluations;
            CHECK(used4 == calls4);
            if (r4) {
                CHECK(r4->used.evaluations == 1 + 5 * r4->used.iterations);    // f(x0), then per step f' (4 calls) and f
            }

            std::uint32_t callsr = 0;
            const auto    rr =
                r::newton {}.with_derivative(d::numeric { d::central_1_2, d::relative { 1e-6 } })(nxx::fn::counted(g_fallible, callsr), x0);
            const auto usedr = rr ? rr->used.evaluations : rr.error().used.evaluations;
            CHECK(usedr == callsr);

            if (r4) ++succeeded;
        }
        CHECK(succeeded > 150);
    }

    TEST_CASE("never || max_evaluations{10}: evaluations_exhausted after exactly 10 evaluations")
    {
        SUBCASE("bisection: two endpoint samples, then one evaluation per step")
        {
            std::uint32_t calls = 0;
            const auto    res =
                r::bisection { nxx::never {} || nxx::max_evaluations { 10 } }(nxx::fn::counted(f_plain, calls), nxx::bracket { 0.0, 2.0 });
            CHECK_FALSE(res.has_value());
            if (!res) {
                CHECK(res.error().code == nxx::errc::evaluations_exhausted);
                CHECK(res.error().by == r::algos::bisection);
                CHECK(res.error().used == nxx::counters { 8, 10 });
                CHECK(res.error().best.has_value());
                if (res.error().best) { CHECK(res.error().best->enclosure.has_value()); }
            }
            CHECK(calls == 10);
        }
        SUBCASE("newton with a numeric derivative: the criterion sees the derivative's evaluations")
        {
            std::uint32_t calls = 0;
            const auto    res   = r::newton { nxx::never {} || nxx::max_evaluations { 10 } }.with_derivative(
                d::numeric {})(nxx::fn::counted(f_plain, calls), rt(100.0));
            CHECK_FALSE(res.has_value());
            if (!res) {
                CHECK(res.error().code == nxx::errc::evaluations_exhausted);
                CHECK(res.error().used == nxx::counters { 3, 10 });    // 1 + 3 steps x (2 + 1)
            }
            CHECK(calls == 10);
        }
    }
}
