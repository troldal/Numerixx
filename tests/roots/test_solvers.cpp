// Unit tests of the spike's 1-D root solvers (DESIGN §7.2, §6.5, §9.2, §9.3): each solver's success, stop reason and
// algorithm id; the failures each one reports, with their code, place, cost and best estimate; run-time inputs in
// every accepted form; evaluation counts against instrumented calls; float and long double defaults; expand's growth
// rules; Newton's derivative sources; and which calls the family facades accept.
#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <expected>
#include <limits>
#include <numbers>
#include <optional>
#include <type_traits>
#include <utility>

// No floating-point contraction in this file's own code on Clang (the library's headers already have it off), so the
// test functions compute the same values, and the solvers take the same paths, on every compiler.
#if defined(__clang__)
#    pragma clang fp contract(off)
#endif

namespace
{
    namespace nr = nxx::roots;

    // A run-time value that no compiler can constant-fold.
    double opaque(double v)
    {
        volatile double sink = v;
        return sink;
    }

    constexpr auto sq2  = [](double x) { return x * x - 2.0; };
    constexpr auto dsq2 = [](double x) { return 2.0 * x; };

    const double root2 = std::sqrt(2.0);

    // The evaluations a result reports, on success and on failure.
    template<class R>
    std::uint32_t evaluations_of(const R& res)
    { return res ? res->used.evaluations : res.error().used.evaluations; }

    // solve rejects through its constraints, so std::is_invocable_v is false rather than a hard error (DESIGN §6.6): a
    // function that cannot take the bracket's scalar type, a list or array with other than two real ends, a 2-D array.
    // A requires-expression rather than a wrapper with a call operator nobody calls, which Clang's -Wunused-template
    // flags; a call that resolves to a deleted overload makes it false, as it makes std::is_invocable_v false.
    template<class... A>
    constexpr bool solve_accepts = requires(A&&... a) { nr::solve(std::forward<A>(a)...); };

    using takes_text = double (*)(const char*);
    using sq2_t      = decltype(sq2);
    static_assert(solve_accepts<sq2_t, const double (&)[2]>);
    static_assert(solve_accepts<sq2_t, std::pair<double, double>>);
    static_assert(!solve_accepts<takes_text, const double (&)[2]>);
    static_assert(!solve_accepts<takes_text, std::pair<double, double>>);
    static_assert(!solve_accepts<sq2_t, const double (&)[1]>);
    static_assert(!solve_accepts<sq2_t, const double (&)[3]>);
    static_assert(!solve_accepts<sq2_t, const int (&)[2]>);
    static_assert(!solve_accepts<sq2_t, const double (&)[2][2]>);

    // The error type of a fallible callback.
    enum class table_error : std::uint8_t { negative_argument = 3, outside_table = 7 };

    // A polynomial with a structural derivative (DESIGN D13 (c)).
    struct linear_poly
    {
        double slope;
        double operator()(double x) const { return slope * x; }
    };

    struct quadratic_poly
    {
        double      c0;
        double      operator()(double x) const { return x * x + c0; }
        linear_poly derivative() const { return linear_poly { 2.0 }; }
    };

    // Every solver's defaults on x^2 - 2 in T (DESIGN §3.5: defaults are achievable in every T).
    template<class T>
    void check_defaults_in()
    {
        const auto quad       = [](T x) { return x * x - T(2); };
        const auto dquad      = [](T x) { return T(2) * x; };
        const T    root       = std::sqrt(T(2));
        const T    tol        = T(16) * std::numeric_limits<T>::epsilon() * root;
        const auto check_root = [&](const auto& res, nxx::algo id) {
            if (res) {
                static_assert(std::is_same_v<std::remove_cvref_t<decltype(res->x)>, T>);
                CHECK(res->by == id);
                CHECK((res->how == nxx::stop_reason::criterion || res->how == nxx::stop_reason::exact_zero));
                CHECK(std::abs(res->x - root) <= tol);
                CHECK(res->fx == quad(res->x));
            }
            else
                FAIL_CHECK("the solver failed on x^2 - 2 with its defaults");
        };
        check_root(nr::bisection {}(quad, { T(1), T(2) }), nr::algos::bisection);
        check_root(nr::brent {}(quad, { T(1), T(2) }), nr::algos::brent);
        check_root(nr::secant {}(quad, T(1)), nr::algos::secant);
        check_root(nr::newton {}.with_derivative(dquad)(quad, T(1)), nr::algos::newton);
    }
}    // namespace

TEST_SUITE("roots")
{
    // ---- 1. Each solver on x^2 - 2 ----------------------------------------------------------------------------------

    TEST_CASE("solvers: bisection on x^2-2")
    {
        const auto res = nr::bisection {}(sq2, { 1.0, 2.0 });
        if (res) {
            CHECK(res->by == nr::algos::bisection);
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(std::abs(res->x - root2) <= res->uncertainty);
            CHECK(res->uncertainty <= 0x1p-50 * 1.5);    // floored_width{}: 2^-50 max(1, min(|lo|, |hi|))
            CHECK(res->fx == sq2(res->x));
            CHECK(res->used.evaluations == res->used.iterations + 2u);
            if (res->enclosure) {
                CHECK(res->enclosure->lo() <= root2);
                CHECK(root2 <= res->enclosure->hi());
            }
            else
                FAIL_CHECK("bisection returned no enclosure");
        }
        else
            FAIL_CHECK("bisection failed on x^2 - 2");
    }

    TEST_CASE("solvers: brent on x^2-2 costs 7 iterations and 9 evaluations")
    {
        const auto res = nr::brent {}(sq2, { 1.0, 2.0 });
        if (res) {
            CHECK(res->by == nr::algos::brent);
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(std::abs(res->x - root2) <= 16.0 * std::numeric_limits<double>::epsilon() * root2);
            CHECK(res->fx == sq2(res->x));
            CHECK(res->used.iterations == 7u);    // DESIGN §7.2
            CHECK(res->used.evaluations == 9u);
        }
        else
            FAIL_CHECK("brent failed on x^2 - 2");
    }

    TEST_CASE("solvers: secant on x^2-2")
    {
        const auto res = nr::secant {}(sq2, 1.0);
        if (res) {
            CHECK(res->by == nr::algos::secant);
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(std::abs(res->x - root2) <= 4.0 * std::numeric_limits<double>::epsilon() * root2);
            CHECK(res->fx == sq2(res->x));
            CHECK(res->used.evaluations == res->used.iterations + 2u);    // x0 and x0 + 2^-10, then one per step
            CHECK_FALSE(res->enclosure.has_value());
        }
        else
            FAIL_CHECK("secant failed on x^2 - 2");
    }

    TEST_CASE("solvers: newton with an analytic derivative on x^2-2")
    {
        const auto res = nr::newton {}.with_derivative(dsq2)(sq2, 1.0);
        if (res) {
            CHECK(res->by == nr::algos::newton);
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(std::abs(res->x - root2) <= 4.0 * std::numeric_limits<double>::epsilon() * root2);
            CHECK(res->fx == sq2(res->x));
            CHECK(res->used.evaluations == 2u * res->used.iterations + 1u);    // f(x0), then f' and f per step
        }
        else
            FAIL_CHECK("newton failed on x^2 - 2");
    }

    TEST_CASE("solvers: expand then bisection on x^2-2")
    {
        // The search: f(2) = 2 and f(2.5) = 4.25; lo moves (smaller |f|) to 2 / 1.6 = 1.25, where f = -0.4375.
        const auto search = nr::expand {}(sq2, nxx::bracket { 2.0, 2.5 });
        if (search) {
            CHECK(search->by == nr::algos::expand);
            CHECK(search->how == nxx::stop_reason::algorithm);
            CHECK(search->lo() == 2.0 / (8.0 / 5.0));
            CHECK(search->hi() == 2.0);    // the tightest pair: hi is the previous lo, not 2.5
            CHECK(search->flo() == sq2(search->lo()));
            CHECK(search->fhi() == 2.0);
            CHECK(search->used == nxx::counters { 1, 3 });
        }
        else
            FAIL_CHECK("expand failed on x^2 - 2 from [2, 2.5]");

        const auto chain = nxx::then(nr::expand {}.on(nxx::bracket { 2.0, 2.5 }), nr::bisection {});
        const auto res   = chain(sq2);
        if (res) {
            CHECK(res->by == nr::algos::bisection);
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(std::abs(res->x - root2) <= res->uncertainty);
            CHECK(res->uncertainty <= 0x1p-50 * 1.5);
            // 3 evaluations and 1 iteration for the search; the bisection stage starts from its samples, so it spends
            // one evaluation per step and never re-evaluates the ends.
            CHECK(res->used.evaluations == res->used.iterations + 2u);
            if (res->enclosure) {
                CHECK(res->enclosure->lo() >= 1.25);
                CHECK(res->enclosure->hi() <= 2.0);
            }
            else
                FAIL_CHECK("the chain returned no enclosure");
        }
        else
            FAIL_CHECK("then(expand, bisection) failed on x^2 - 2");
    }

    // ---- 2. No sign change ------------------------------------------------------------------------------------------

    TEST_CASE("solvers: no sign change fails before iterating")
    {
        const auto sq5            = [](double x) { return x * x - 5.0; };
        const auto check_no_roots = [&](const auto& res, nxx::algo id) {
            CHECK_FALSE(res.has_value());
            if (!res) {
                const auto& err = res.error();
                CHECK(err.code == nxx::errc::no_sign_change);
                CHECK(err.where == id);
                CHECK(err.used == nxx::counters { 0, 2 });
                if (err.best) {
                    CHECK(err.best->x == 3.0);    // the end with the smaller |f|
                    CHECK(err.best->fx == 4.0);
                    CHECK_FALSE(err.best->enclosure.has_value());
                }
                else
                    FAIL_CHECK("no best estimate");
            }
        };
        check_no_roots(nr::bisection {}(sq5, { 3.0, 4.0 }), nr::algos::bisection);
        check_no_roots(nr::brent {}(sq5, { 3.0, 4.0 }), nr::algos::brent);
    }

    // ---- 3. Fallible callbacks --------------------------------------------------------------------------------------

    TEST_CASE("solvers: a fallible callback's error is the failure's cause")
    {
        const auto guarded = [](double x) -> std::expected<double, table_error> {
            if (x < 0.0) return std::unexpected(table_error::negative_argument);
            return x - 0.5;
        };
        const auto check_first_sample = [](const auto& res, nxx::algo id) {
            static_assert(std::is_same_v<typename std::remove_cvref_t<decltype(res)>::error_type,
                                         nxx::failure<nr::root_estimate<double>, table_error>>);
            CHECK_FALSE(res.has_value());
            if (!res) {
                const auto& err = res.error();
                CHECK(err.code == nxx::errc::callback_failed);
                CHECK(err.where == id);
                CHECK(err.used == nxx::counters { 0, 1 });    // the failing call counts
                CHECK_FALSE(err.best.has_value());            // nothing was evaluated successfully
                if (err.cause)
                    CHECK(*err.cause == table_error::negative_argument);
                else
                    FAIL_CHECK("no cause");
            }
        };
        check_first_sample(nr::bisection {}(guarded, { -1.0, 2.0 }), nr::algos::bisection);
        check_first_sample(nr::brent {}(guarded, { -1.0, 2.0 }), nr::algos::brent);

        // A failure while iterating: bisection on [0, 2] evaluates 1 (fine), then 0.5 (refused).
        const auto gappy = [](double x) -> std::expected<double, table_error> {
            if (x > 0.4 && x < 0.6) return std::unexpected(table_error::outside_table);
            return x - 0.3;
        };
        const auto res = nr::bisection {}(gappy, { 0.0, 2.0 });
        CHECK_FALSE(res.has_value());
        if (!res) {
            const auto& err = res.error();
            CHECK(err.code == nxx::errc::callback_failed);
            CHECK(err.where == nr::algos::bisection);
            CHECK(err.used == nxx::counters { 2, 4 });
            if (err.cause)
                CHECK(*err.cause == table_error::outside_table);
            else
                FAIL_CHECK("no cause");
            if (err.best && err.best->enclosure) {
                CHECK(err.best->x == 0.0);
                CHECK(err.best->enclosure->lo() == 0.0);
                CHECK(err.best->enclosure->hi() == 1.0);    // the enclosure after the first step
            }
            else
                FAIL_CHECK("no best estimate with an enclosure");
        }
    }

    // ---- 4. Budget exhaustion ---------------------------------------------------------------------------------------

    TEST_CASE("solvers: newton on x^2+1 exhausts its default budget")
    {
        const auto no_real_root = [](double x) { return x * x + 1.0; };
        const auto res          = nr::newton {}.with_derivative(dsq2)(no_real_root, 0.5);
        CHECK_FALSE(res.has_value());
        if (!res) {
            const auto& err = res.error();
            CHECK(err.code == nxx::errc::budget_exhausted);
            CHECK(err.where == nr::algos::newton);
            CHECK(err.used.iterations == 30u);
            CHECK(err.used.evaluations == 61u);
            if (err.best) {
                CHECK(err.best->fx == no_real_root(err.best->x));
                CHECK(err.best->fx <= 1.25);    // no worse than the start
            }
            else
                FAIL_CHECK("no best estimate");
        }
    }

    // ---- 5. Projection ----------------------------------------------------------------------------------------------

    TEST_CASE("solvers: newton with a projection converges from outside the box")
    {
        double     first_arg = std::numeric_limits<double>::quiet_NaN();
        const auto recorded  = [&first_arg](double x) {
            if (std::isnan(first_arg)) first_arg = x;
            return x * x - 2.0;
        };
        const auto res = nr::newton {}.with_derivative(dsq2).with_projection(nr::clamp_to { 1.0, 3.0 })(recorded, 10.0);
        CHECK(first_arg == 3.0);    // the guess is projected before it is evaluated
        if (res) {
            CHECK(res->by == nr::algos::newton);
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(std::abs(res->x - root2) <= 4.0 * std::numeric_limits<double>::epsilon() * root2);
        }
        else
            FAIL_CHECK("projected newton failed on x^2 - 2");
    }

    TEST_CASE("solvers: an iterate pinned at the projection's edge is stalled")
    {
        const auto beyond     = [](double x) { return x - 5.0; };
        const auto unit_slope = [](double) { return 1.0; };

        // Newton: 1 -> 5, clamped to 3; then 3 -> 5, clamped to 3 again: pinned.
        const auto nres = nr::newton {}.with_derivative(unit_slope).with_projection(nr::clamp_to { 0.0, 3.0 })(beyond, 1.0);
        CHECK_FALSE(nres.has_value());
        if (!nres) {
            const auto& err = nres.error();
            CHECK(err.code == nxx::errc::stalled);
            CHECK(err.where == nr::algos::newton);
            CHECK(err.used == nxx::counters { 2, 4 });    // f(1); f'(1), f(3); f'(3)
            if (err.best) {
                CHECK(err.best->x == 3.0);
                CHECK(err.best->fx == -2.0);
            }
            else
                FAIL_CHECK("no best estimate");
        }
        CHECK(nxx::best_x(nres).value_or(0.0) == 3.0);

        // Secant: the same, from x0 = 1 and x1 = 1 + 2^-10.
        const auto sres = nr::secant {}.with_projection(nr::clamp_to { 0.0, 3.0 })(beyond, 1.0);
        CHECK_FALSE(sres.has_value());
        if (!sres) {
            const auto& err = sres.error();
            CHECK(err.code == nxx::errc::stalled);
            CHECK(err.where == nr::algos::secant);
            CHECK(err.used == nxx::counters { 2, 3 });    // f(x0), f(x1); f(3); pinned before evaluating
            if (err.best) {
                CHECK(err.best->x == 3.0);
                CHECK(err.best->fx == -2.0);
            }
            else
                FAIL_CHECK("no best estimate");
        }
        CHECK(nxx::best_x(sres).value_or(0.0) == 3.0);
    }

    TEST_CASE("solvers: an iterate projected off the reals is rejected before f is evaluated there")
    {
        constexpr double inf      = std::numeric_limits<double>::infinity();
        std::uint32_t    calls    = 0;
        std::uint32_t    df_calls = 0;
        const auto       reset    = [&calls, &df_calls] { calls = df_calls = 0; };

        // The start. clamp_to{inf, inf}, or reversed bounds with lo = inf, sends the guess to +inf, where 1/x is exactly
        // 0; a custom projection to NaN does the same for a three-way comparison, which is 0 at NaN. Both methods
        // reported an exact_zero success at the non-finite x after one evaluation there.
        const auto recip          = nxx::fn::counted([](double x) { return 1.0 / x; }, calls);
        const auto drecip         = nxx::fn::counted([](double x) { return -1.0 / (x * x); }, df_calls);
        const auto three_way      = nxx::fn::counted([](double x) { return x < 1.0 ? -1.0 : (x > 1.0 ? 1.0 : 0.0); }, calls);
        const auto to_nan         = [](double) { return std::numeric_limits<double>::quiet_NaN(); };
        const auto start_rejected = [&calls, &df_calls, &reset](const auto& res, nxx::algo where) {
            CHECK_FALSE(res.has_value());
            if (!res) {
                const auto& err = res.error();
                CHECK(err.code == nxx::errc::non_finite_input);
                CHECK(err.where == where);
                CHECK(err.used == nxx::counters {});
                CHECK_FALSE(err.best.has_value());
            }
            CHECK(calls + df_calls == 0u);    // never called at the non-finite start
            reset();
        };
        start_rejected(nr::secant {}.with_projection(nr::clamp_to { inf, inf })(recip, 5.0), nr::algos::secant);
        start_rejected(nr::secant {}.with_projection(nr::clamp_to { inf, 0.0 })(recip, 5.0), nr::algos::secant);
        start_rejected(nr::newton {}.with_derivative(drecip).with_projection(nr::clamp_to { inf, inf })(recip, 5.0), nr::algos::newton);
        start_rejected(nr::secant {}.with_projection(to_nan)(three_way, 3.0), nr::algos::secant);

        // A later iterate. A projection that marks x <= 0 as outside the domain with +inf: the first step of either
        // method from 10 on 1/x - 1/4 goes to about -5, projected to +inf, where f = -1/4 is finite. step_tol's threshold
        // 2^-k max(|x|, 1) is inf there, so both methods reported a criterion success at x = inf.
        const auto pos_or_inf    = [](double x) { return x > 0.0 ? x : std::numeric_limits<double>::infinity(); };
        const auto quarter       = nxx::fn::counted([](double x) { return 1.0 / x - 0.25; }, calls);
        const auto step_rejected = [&calls, &df_calls, &reset](const auto& res, nxx::algo where, nxx::counters used) {
            CHECK_FALSE(res.has_value());
            if (!res) {
                const auto& err = res.error();
                CHECK(err.code == nxx::errc::diverged);
                CHECK(err.where == where);
                CHECK(err.used == used);
                CHECK(err.used.evaluations == calls + df_calls);
                if (err.best) {
                    CHECK(err.best->x == 10.0);    // the start: |f| grows from 10 to 10 + 2^-10 10
                    CHECK(err.best->fx == 1.0 / 10.0 - 0.25);
                }
                else
                    FAIL_CHECK("no best estimate");
            }
            reset();
        };
        step_rejected(nr::secant {}.with_projection(pos_or_inf)(quarter, 10.0),
                      nr::algos::secant,
                      nxx::counters { 1, 2 });    // f(x0), f(x1)
        step_rejected(nr::newton {}.with_derivative(drecip).with_projection(pos_or_inf)(quarter, 10.0),
                      nr::algos::newton,
                      nxx::counters { 1, 2 });    // f(10); f'(10)

        // The same domain marked with a far finite value instead: the criteria saw only the proposed step (about 15),
        // and step_tol's threshold at x = 1e300 is huge, so both methods reported a criterion success at x = 1e300 or
        // max with f = -1/4. The step the criteria see is now the larger of the proposed and the actual one.
        for (const double remote : { 1e300, (std::numeric_limits<double>::max)() }) {
            CAPTURE(remote);
            const auto pos_or_remote = [remote](double x) { return x > 0.0 ? x : remote; };
            const auto sec           = nr::secant {}.with_projection(pos_or_remote)(quarter, 10.0);
            CHECK_FALSE(sec.has_value());
            if (!sec) CHECK(sec.error().used.evaluations == calls + df_calls);
            reset();
            const auto nwt = nr::newton {}.with_derivative(drecip).with_projection(pos_or_remote)(quarter, 10.0);
            CHECK_FALSE(nwt.has_value());
            if (!nwt) CHECK(nwt.error().used.evaluations == calls + df_calls);
            reset();
        }

        // Secant's second point: the guess 3 stays, and 3 + 2^-10 3 is projected to +inf. As where x0 + h overflows,
        // the second point is then 3 - 2^-10 3, and x - 2 solves.
        const auto above3   = [](double x) { return x > 3.0 ? std::numeric_limits<double>::infinity() : x; };
        const auto line     = nxx::fn::counted([](double x) { return x - 2.0; }, calls);
        const auto fallback = nr::secant {}.with_projection(above3)(line, 3.0);
        CHECK(fallback.has_value());
        if (fallback) CHECK(fallback->x == 2.0);
        reset();

        // A projection that pins one neighbour of the guess and sends the other off the reals: diverged, whichever side
        // it pins, after f(3) alone (the rule for both: off the reals beats pinned).
        const auto up_off   = [](double x) { return x > 3.0 ? std::numeric_limits<double>::infinity() : 3.0; };
        const auto down_off = [](double x) { return x < 3.0 ? std::numeric_limits<double>::infinity() : 3.0; };
        const auto r_up     = nr::secant {}.with_projection(up_off)(recip, 3.0);
        const auto r_down   = nr::secant {}.with_projection(down_off)(recip, 3.0);
        CHECK_FALSE(r_up.has_value());
        CHECK_FALSE(r_down.has_value());
        if (!r_up && !r_down) {
            CHECK(r_up.error().code == nxx::errc::diverged);
            CHECK(r_down.error().code == nxx::errc::diverged);
            CHECK(r_up.error().used == nxx::counters { 0, 1 });
            CHECK(r_down.error().used == nxx::counters { 0, 1 });
        }
        CHECK(calls == 2u);
        reset();

        // Reversed finite bounds send the iterate to the far bound: secant reported a criterion success at x = 1, with
        // f = 1 - 1e-15, from 0.5, because the criteria saw only the proposed step. Now the larger step is seen.
        const auto shifted  = nxx::fn::counted([](double x) { return x - 1e-15; }, calls);
        const auto reversed = nr::secant {}.with_projection(nr::clamp_to { 1.0, 0.0 })(shifted, 0.5);
        CHECK_FALSE(reversed.has_value());
        if (!reversed) {
            CHECK(reversed.error().code == nxx::errc::stalled);
            CHECK(reversed.error().used.evaluations == calls);
        }
        reset();

        // Only a projection that sends both neighbours of the guess off the reals fails: diverged, after f(3) alone.
        const auto only3  = [](double x) { return x == 3.0 ? x : std::numeric_limits<double>::infinity(); };
        const auto second = nr::secant {}.with_projection(only3)(recip, 3.0);
        CHECK_FALSE(second.has_value());
        if (!second) {
            const auto& err = second.error();
            CHECK(err.code == nxx::errc::diverged);
            CHECK(err.used == nxx::counters { 0, 1 });    // f(3) only
            if (err.best) {
                CHECK(err.best->x == 3.0);
                CHECK(err.best->fx == 1.0 / 3.0);
            }
            else
                FAIL_CHECK("no best estimate");
        }
        CHECK(calls == 1u);
    }

    // ---- 6. Open-method failures ------------------------------------------------------------------------------------

    TEST_CASE("solvers: a flat secant is stalled")
    {
        const auto level = [](double) { return 1.0; };
        const auto res   = nr::secant {}(level, 1.0);
        CHECK_FALSE(res.has_value());
        if (!res) {
            const auto& err = res.error();
            CHECK(err.code == nxx::errc::stalled);
            CHECK(err.where == nr::algos::secant);
            CHECK(err.used == nxx::counters { 1, 2 });
            if (err.best)
                CHECK(err.best->fx == 1.0);
            else
                FAIL_CHECK("no best estimate");
        }
    }

    TEST_CASE("solvers: newton at a zero derivative")
    {
        const auto res = nr::newton {}.with_derivative(dsq2)(sq2, 0.0);
        CHECK_FALSE(res.has_value());
        if (!res) {
            const auto& err = res.error();
            CHECK(err.code == nxx::errc::zero_derivative);
            CHECK(err.where == nr::algos::newton);
            CHECK(err.used == nxx::counters { 1, 2 });    // f(0), f'(0)
            if (err.best) {
                CHECK(err.best->x == 0.0);
                CHECK(err.best->fx == -2.0);
            }
            else
                FAIL_CHECK("no best estimate");
        }
    }

    TEST_CASE("solvers: a non-finite proposed step diverges")
    {
        // Newton: the step 1e300 / 1e-300 overflows before the iterate is evaluated.
        const auto huge       = [](double) { return 1e300; };
        const auto tiny_slope = [](double) { return 1e-300; };
        const auto nres       = nr::newton {}.with_derivative(tiny_slope)(huge, 1.0);
        CHECK_FALSE(nres.has_value());
        if (!nres) {
            const auto& err = nres.error();
            CHECK(err.code == nxx::errc::diverged);
            CHECK(err.where == nr::algos::newton);
            CHECK(err.used == nxx::counters { 1, 2 });    // f(1), f'(1); the proposed iterate is never evaluated
            if (err.best)
                CHECK(err.best->x == 1.0);
            else
                FAIL_CHECK("no best estimate");
        }

        // Secant: from x0 = 1e300 and x1 = x0 + 2^-10 x0, a secant slope of 2^-52 / 1e297 proposes a step of about
        // 4e312, which is not representable.
        const auto step_up = [](double x) { return x < 1.0005e300 ? 1.0 : 1.0 + 0x1p-52; };
        const auto sres    = nr::secant {}(step_up, 1e300);
        CHECK_FALSE(sres.has_value());
        if (!sres) {
            const auto& err = sres.error();
            CHECK(err.code == nxx::errc::diverged);
            CHECK(err.where == nr::algos::secant);
            CHECK(err.used == nxx::counters { 1, 2 });
            if (err.best)
                CHECK(err.best->x == 1e300);
            else
                FAIL_CHECK("no best estimate");
        }
    }

    // ---- 7. Poles ---------------------------------------------------------------------------------------------------

    TEST_CASE("solvers: a pole is not a root")
    {
        const auto tangent    = [](double x) { return std::tan(x); };
        const auto hyperbola  = [](double x) { return 1.0 / (x - 1.0 / 3.0); };
        const auto check_pole = [](const auto& res, nxx::algo id, double pole) {
            CHECK_FALSE(res.has_value());
            if (!res) {
                const auto& err = res.error();
                CHECK(err.code == nxx::errc::sign_change_not_root);
                CHECK(err.where == id);
                CHECK(err.used.iterations > 0u);
                CHECK(err.used.evaluations == err.used.iterations + 2u);
                if (err.best && err.best->enclosure) {    // the final enclosure, around the pole
                    CHECK(err.best->enclosure->lo() <= pole);
                    CHECK(pole <= err.best->enclosure->hi());
                    CHECK(std::abs(err.best->fx) > 1e10);
                }
                else
                    FAIL_CHECK("no best estimate with an enclosure");
            }
        };
        check_pole(nr::bisection {}(tangent, { 1.0, 2.0 }), nr::algos::bisection, std::numbers::pi / 2.0);
        check_pole(nr::brent {}(tangent, { 1.0, 2.0 }), nr::algos::brent, std::numbers::pi / 2.0);
        check_pole(nr::bisection {}(hyperbola, { 0.0, 1.0 }), nr::algos::bisection, 1.0 / 3.0);
        check_pole(nr::brent {}(hyperbola, { 0.0, 1.0 }), nr::algos::brent, 1.0 / 3.0);
    }

    TEST_CASE("solvers: a pole at an end with an infinite sample is not a root")
    {
        // The pole check's reference is the larger *finite* initial sample: an infinite one would switch it off.
        const auto reciprocal = [](double x) { return 1.0 / x; };
        const auto shifted    = [](double x) { return 1.0 / (x - 1.0); };
        const auto cubed      = [](double x) { return 1.0 / (x * x * x); };    // infinite at both ends of the bracket
        const auto check_pole = [](const auto& res, nxx::algo id) {
            CHECK_FALSE(res.has_value());
            if (!res) {
                CHECK(res.error().code == nxx::errc::sign_change_not_root);
                CHECK(res.error().where == id);
            }
        };
        check_pole(nr::bisection {}(reciprocal, { -1.0, 0.0 }), nr::algos::bisection);
        check_pole(nr::brent {}(reciprocal, { -1.0, 0.0 }), nr::algos::brent);
        check_pole(nr::solve(reciprocal, { -1.0, 0.0 }), nr::algos::brent);
        check_pole(nr::bisection {}(shifted, { 0.0, 1.0 }), nr::algos::bisection);
        check_pole(nr::brent {}(shifted, { 0.0, 1.0 }), nr::algos::brent);
        // Both samples infinite: only the finiteness test applies, and it rejects the infinite fx at the returned point.
        check_pole(nr::bisection {}(cubed, { -1e-200, 1e-200 }), nr::algos::bisection);
        check_pole(nr::brent {}(cubed, { -1e-200, 1e-200 }), nr::algos::brent);
    }

    // ---- 8. Infinite endpoint samples -------------------------------------------------------------------------------

    TEST_CASE("solvers: log from 0 to 2 starts from an infinite endpoint sample")
    {
        const auto logarithm = [](double x) { return std::log(x); };
        const auto check_one = [](const auto& res, nxx::algo id) {
            if (res) {
                CHECK(res->by == id);
                CHECK(std::abs(res->x - 1.0) <= 4.0 * std::numeric_limits<double>::epsilon());
                CHECK(std::isfinite(res->fx));
            }
            else
                FAIL_CHECK("log(x) on [0, 2] failed");
        };
        check_one(nr::bisection {}(logarithm, { 0.0, 2.0 }), nr::algos::bisection);
        check_one(nr::brent {}(logarithm, { 0.0, 2.0 }), nr::algos::brent);
    }

    // ---- 9. A root exactly at 0 -------------------------------------------------------------------------------------

    TEST_CASE("solvers: a root exactly at 0 terminates with the default tolerance")
    {
        const auto odd_cubic  = [](double x) { return x + x * x * x; };
        const auto check_zero = [](const auto& res, nxx::algo id) {
            if (res) {
                CHECK(res->by == id);
                CHECK((res->how == nxx::stop_reason::criterion || res->how == nxx::stop_reason::exact_zero));
                CHECK(std::abs(res->x) <= 0x1p-49);    // the absolute floor of floored_width{} at scale 1
                CHECK(res->used.iterations <= 100u);
            }
            else
                FAIL_CHECK("x + x^3 on [-1, 2] failed");
        };
        check_zero(nr::bisection {}(odd_cubic, { -1.0, 2.0 }), nr::algos::bisection);
        check_zero(nr::brent {}(odd_cubic, { -1.0, 2.0 }), nr::algos::brent);
    }

    // ---- 10. Extreme brackets ---------------------------------------------------------------------------------------

    TEST_CASE("solvers: bisection's midpoint is overflow-safe on an extreme bracket")
    {
        const auto linear = [](double x) { return x; };
        const auto res    = nr::bisection {}(linear, { -1.7e308, 1.7e308 });
        if (res) {
            CHECK(res->how == nxx::stop_reason::exact_zero);    // not a false resolution_limit at -1.7e308
            CHECK(res->x == 0.0);
            CHECK(res->used == nxx::counters { 1, 3 });
        }
        else
            FAIL_CHECK("bisection failed on [-1.7e308, 1.7e308]");
    }

    TEST_CASE("solvers: brent on an extreme bracket")
    {
        // DESIGN §7.2, common rules for bracketing solvers: the midpoint is overflow-safe.
        const auto linear = [](double x) { return x; };
        const auto res    = nr::brent {}(linear, { -1.7e308, 1.7e308 });
        if (res) {
            CHECK(std::isfinite(res->x));
            CHECK(std::abs(res->x) <= 0x1p-49);
        }
        else
            FAIL_CHECK("brent failed on [-1.7e308, 1.7e308]");
    }

    TEST_CASE("solvers: roots at 1e300")
    {
        // DESIGN §9.2: roots at 1e300.
        const auto shifted = [](double x) { return x - 1e300; };
        const auto unit    = [](double) { return 1.0; };
        const auto near    = [](const auto& res) {
            if (res)
                CHECK(std::abs(res->x - 1e300) <= 1e286);
            else
                FAIL_CHECK("no root found at 1e300");
        };
        near(nr::bisection {}(shifted, { 1e299, 1e301 }));
        near(nr::brent {}(shifted, { 1e299, 1e301 }));
        near(nr::newton {}.with_derivative(unit)(shifted, 1.1e300));
    }

    TEST_CASE("solvers: secant on a root at 1e300")
    {
        // DESIGN §9.2: roots at 1e300. The secant step through f1 (x1 - x0) / (f1 - f0) must not overflow here: the
        // true step is about -1e299.
        const auto shifted = [](double x) { return x - 1e300; };
        const auto res     = nr::secant {}(shifted, 1.1e300);
        if (res)
            CHECK(std::abs(res->x - 1e300) <= 1e286);
        else
            FAIL_CHECK("secant found no root at 1e300");
    }

    // ---- 11. Run-time inputs ----------------------------------------------------------------------------------------

    TEST_CASE("solvers: run-time brackets in every accepted form")
    {
        const double lo = opaque(1.0);
        const double hi = opaque(2.0);

        const auto check_forms = [&](const auto& solver) {
            const auto   fixed    = solver(sq2, nxx::bracket { 1.0, 2.0 });    // the literal, checked at compile time
            const auto   braced   = solver(sq2, { lo, hi });
            const auto   paired   = solver(sq2, std::pair { lo, hi });
            const auto   made     = solver(sq2, nxx::bracket<double>::make(lo, hi));
            const auto   reversed = solver(sq2, std::pair { hi, lo });    // make() re-orders
            const auto   made_rev = solver(sq2, nxx::bracket<double>::make(hi, lo));
            const auto   curried  = solver.on({ lo, hi })(sq2);
            const double arr[2]   = { lo, hi };    // an lvalue C array, by reference and by .on (cl once found both ambiguous)
            const auto   from_arr = solver(sq2, arr);
            const auto   on_arr   = solver.on(arr)(sq2);
            CHECK(from_arr.has_value());
            CHECK(on_arr.has_value());
            if (from_arr && on_arr && fixed) {
                CHECK(from_arr->x == fixed->x);
                CHECK(on_arr->x == fixed->x);
            }
            CHECK(fixed.has_value());
            CHECK(braced.has_value());
            CHECK(paired.has_value());
            CHECK(made.has_value());
            CHECK(reversed.has_value());
            CHECK(made_rev.has_value());
            CHECK(curried.has_value());
            if (fixed && braced && paired && made && reversed && made_rev && curried) {
                CHECK(std::abs(fixed->x - root2) <= fixed->uncertainty);
                for (const auto* other : { &*braced, &*paired, &*made, &*reversed, &*made_rev, &*curried }) {
                    CHECK(other->x == fixed->x);
                    CHECK(other->used == fixed->used);
                }
            }
        };
        check_forms(nr::bisection {});
        check_forms(nr::brent {});

        // solve(f, {lo, hi}) is brent (canonical call 1).
        const auto solved = nr::solve(sq2, { lo, hi });
        const auto direct = nr::brent {}(sq2, { lo, hi });
        if (solved && direct) {
            CHECK(solved->x == direct->x);
            CHECK(solved->used == direct->used);
            CHECK(solved->by == nr::algos::brent);
        }
        else
            FAIL_CHECK("solve or brent failed on a run-time bracket");

        // Run-time guesses.
        const double x0     = opaque(1.0);
        const auto   secres = nr::secant {}(sq2, x0);
        const auto   newres = nr::newton {}.with_derivative(dsq2)(sq2, x0);
        CHECK(secres.has_value());
        CHECK(newres.has_value());
    }

    TEST_CASE("solvers: an open method started at an exact root stops at once with an unknown uncertainty")
    {
        const auto check_start = [](const auto& res, nxx::algo id) {
            if (res) {
                CHECK(res->by == id);
                CHECK(res->how == nxx::stop_reason::exact_zero);
                CHECK(res->x == 2.0);
                CHECK(res->used.iterations == 0u);
                CHECK(std::isinf(res->uncertainty));    // no step was taken: unknown, not 0
            }
            else
                FAIL_CHECK("a solve from an exact root failed");
        };
        const auto shifted = [](double x) { return x - 2.0; };
        check_start(nr::secant {}(shifted, 2.0), nr::algos::secant);
        check_start(nr::newton {}.with_derivative([](double) { return 1.0; })(shifted, 2.0), nr::algos::newton);
    }

    TEST_CASE("solvers: invalid run-time inputs fail in-band at zero cost")
    {
        const double one          = opaque(1.0);
        const double not_a_number = opaque(std::numeric_limits<double>::quiet_NaN());
        const double infinite     = opaque(std::numeric_limits<double>::infinity());

        std::uint32_t calls       = 0;
        const auto    counted     = nxx::fn::counted(sq2, calls);
        const auto    check_input = [&calls](const auto& res, nxx::errc code, nxx::algo id) {
            CHECK_FALSE(res.has_value());
            if (!res) {
                const auto& err = res.error();
                CHECK(err.code == code);
                CHECK(err.where == id);
                CHECK(err.used == nxx::counters {});
                CHECK_FALSE(err.best.has_value());
            }
            CHECK(calls == 0u);
        };

        // Equal endpoints.
        check_input(nr::bisection {}(counted, { one, one }), nxx::errc::invalid_input, nr::algos::bisection);
        check_input(nr::brent {}(counted, { one, one }), nxx::errc::invalid_input, nr::algos::brent);
        check_input(nr::expand {}(counted, { one, one }), nxx::errc::invalid_input, nr::algos::expand);
        check_input(nr::bisection {}(counted, std::pair { one, one }), nxx::errc::invalid_input, nr::algos::bisection);
        // Non-finite endpoints.
        check_input(nr::bisection {}(counted, std::pair { not_a_number, one }), nxx::errc::invalid_input, nr::algos::bisection);
        check_input(nr::brent {}(counted, { one, infinite }), nxx::errc::invalid_input, nr::algos::brent);
        // make()'s error, forwarded unchanged.
        check_input(nr::bisection {}(counted, nxx::bracket<double>::make(one, one)), nxx::errc::invalid_input, nr::algos::bisection);
        const std::expected<nxx::bracket<double>, nxx::errc> refused { std::unexpect, nxx::errc::out_of_domain };
        check_input(nr::bisection {}(counted, refused), nxx::errc::out_of_domain, nr::algos::bisection);
        check_input(nr::brent {}(counted, refused), nxx::errc::out_of_domain, nr::algos::brent);
        check_input(nr::expand {}(counted, refused), nxx::errc::out_of_domain, nr::algos::expand);
        // Non-finite guesses.
        check_input(nr::secant {}(counted, infinite), nxx::errc::non_finite_input, nr::algos::secant);
        check_input(nr::secant {}(counted, not_a_number), nxx::errc::non_finite_input, nr::algos::secant);
        check_input(nr::newton {}.with_derivative(dsq2)(counted, infinite), nxx::errc::non_finite_input, nr::algos::newton);
        check_input(nr::newton {}.with_derivative(dsq2)(counted, not_a_number), nxx::errc::non_finite_input, nr::algos::newton);
    }

    // ---- 12. Evaluation counts --------------------------------------------------------------------------------------

    TEST_CASE("solvers: evaluation counts equal instrumented calls")
    {
        std::uint32_t calls    = 0;
        std::uint32_t df_calls = 0;
        const auto    cf       = nxx::fn::counted(sq2, calls);
        const auto    cdf      = nxx::fn::counted(dsq2, df_calls);
        const auto    count_is = [&calls, &df_calls](const auto& res) {
            CHECK(evaluations_of(res) == calls + df_calls);
            calls    = 0;
            df_calls = 0;
        };

        // Successes.
        count_is(nr::bisection {}(cf, { 1.0, 2.0 }));
        count_is(nr::brent {}(cf, { 1.0, 2.0 }));
        count_is(nr::secant {}(cf, 1.0));
        count_is(nr::newton {}.with_derivative(cdf)(cf, 1.0));
        count_is(nr::expand {}(cf, nxx::bracket { 2.0, 2.5 }));
        const auto chained = nxx::then(nr::expand {}.on(nxx::bracket { 2.0, 2.5 }), nr::bisection {})(cf);
        if (chained) CHECK(chained->used.evaluations == chained->used.iterations + 2u);    // no re-evaluation of the ends
        count_is(chained);

        // Failures.
        const auto c5 = nxx::fn::counted([](double x) { return x * x - 5.0; }, calls);
        count_is(nr::bisection {}(c5, { 3.0, 4.0 }));    // no sign change
        count_is(nr::brent {}(c5, { 3.0, 4.0 }));
        const auto cplus = nxx::fn::counted([](double x) { return x * x + 1.0; }, calls);
        count_is(nr::newton {}.with_derivative(cdf)(cplus, 0.5));    // budget exhausted
        count_is(nr::expand {}(cplus, { 1.0, 2.0 }));                // budget exhausted
        count_is(nr::newton {}.with_derivative(cdf)(cf, 0.0));       // zero derivative
        const auto clevel = nxx::fn::counted([](double) { return 1.0; }, calls);
        count_is(nr::secant {}(clevel, 1.0));    // flat
        const auto cbeyond = nxx::fn::counted([](double x) { return x - 5.0; }, calls);
        const auto cunit   = nxx::fn::counted([](double) { return 1.0; }, df_calls);
        count_is(nr::newton {}.with_derivative(cunit).with_projection(nr::clamp_to { 0.0, 3.0 })(cbeyond, 1.0));    // pinned
        count_is(nr::secant {}.with_projection(nr::clamp_to { 0.0, 3.0 })(cbeyond, 1.0));
        const auto ctan = nxx::fn::counted([](double x) { return std::tan(x); }, calls);
        count_is(nr::bisection {}(ctan, { 1.0, 2.0 }));    // a pole
        count_is(nr::brent {}(ctan, { 1.0, 2.0 }));
        const auto cgappy = nxx::fn::counted(
            [](double x) -> std::expected<double, table_error> {
                if (x > 0.4 && x < 0.6) return std::unexpected(table_error::outside_table);
                return x - 0.3;
            },
            calls);
        count_is(nr::bisection {}(cgappy, { 0.0, 2.0 }));    // the callback's own failure
        count_is(nr::brent {}(cgappy, { -1.0, 0.45 }));      // fails at the second sample
        const auto chuge = nxx::fn::counted([](double) { return 1e300; }, calls);
        const auto ctiny = nxx::fn::counted([](double) { return 1e-300; }, df_calls);
        count_is(nr::newton {}.with_derivative(ctiny)(chuge, 1.0));    // diverged
        const auto cstep = nxx::fn::counted([](double x) { return x < 1.0005e300 ? 1.0 : 1.0 + 0x1p-52; }, calls);
        count_is(nr::secant {}(cstep, 1e300));
    }

    // ---- 13. float and long double ----------------------------------------------------------------------------------

    TEST_CASE("solvers: float defaults converge on x^2-2") { check_defaults_in<float>(); }
    TEST_CASE("solvers: double defaults converge on x^2-2") { check_defaults_in<double>(); }
    TEST_CASE("solvers: long double defaults converge on x^2-2") { check_defaults_in<long double>(); }

    // ---- 14. expand -------------------------------------------------------------------------------------------------

    TEST_CASE("solvers: expand grows a positive window geometrically")
    {
        const auto far_right = [](double x) { return x - 1000.0; };
        const auto res       = nr::expand {}(far_right, { 1.0, 2.0 });
        static_assert(std::is_same_v<typename std::remove_cvref_t<decltype(res)>::value_type, nxx::solution<nr::sign_bracket<double>>>);
        if (res) {
            CHECK(res->by == nr::algos::expand);
            CHECK(res->how == nxx::stop_reason::algorithm);
            CHECK(res->lo() < 1000.0);
            CHECK(1000.0 < res->hi());
            CHECK(res->flo() < 0.0);
            CHECK(res->fhi() > 0.0);
            CHECK(res->flo() == far_right(res->lo()));
            CHECK(res->fhi() == far_right(res->hi()));
            CHECK(res->hi() == res->lo() * (8.0 / 5.0));     // the tightest pair: its lo is the previous hi
            CHECK(res->used == nxx::counters { 14, 16 });    // hi = 2 * 1.6^14 is the first beyond 1000
        }
        else
            FAIL_CHECK("expand failed on x - 1000 from [1, 2]");
    }

    TEST_CASE("solvers: expand mirrors its geometric growth for a negative window")
    {
        const auto far_left = [](double x) { return x + 1000.0; };
        const auto res      = nr::expand {}(far_left, { -2.0, -1.0 });
        if (res) {
            CHECK(res->by == nr::algos::expand);
            CHECK(res->how == nxx::stop_reason::algorithm);
            CHECK(res->lo() < -1000.0);
            CHECK(-1000.0 < res->hi());
            CHECK(res->flo() < 0.0);
            CHECK(res->fhi() > 0.0);
            CHECK(res->lo() == res->hi() * (8.0 / 5.0));    // the tightest pair: its hi is the previous lo
            CHECK(res->used == nxx::counters { 14, 16 });
        }
        else
            FAIL_CHECK("expand failed on x + 1000 from [-2, -1]");
    }

    TEST_CASE("solvers: expand without a root exhausts its budget with a best estimate")
    {
        const auto positive = [](double x) { return x * x + 1.0; };
        const auto res      = nr::expand {}(positive, { 1.0, 2.0 });
        static_assert(std::is_same_v<typename std::remove_cvref_t<decltype(res)>::error_type, nxx::failure<nr::root_estimate<double>>>);
        CHECK_FALSE(res.has_value());
        if (!res) {
            const auto& err = res.error();
            CHECK(err.code == nxx::errc::budget_exhausted);
            CHECK(err.where == nr::algos::expand);
            CHECK(err.used == nxx::counters { 60, 62 });
            if (err.best) {
                CHECK(err.best->fx == positive(err.best->x));
                CHECK(std::isinf(err.best->uncertainty));    // neither an enclosure nor a step: unknown, not the width
                CHECK(err.best->x > 0.0);                    // lo moved towards 0 by lo / 1.6, keeping the smaller |f|
                CHECK(err.best->x < 1.0);
            }
            else
                FAIL_CHECK("no best estimate");
        }
    }

    TEST_CASE("solvers: expand saturates at the largest finite values")
    {
        // Unclamped growth reached hi = inf, and bisection then called [1.25e308, inf] unsplittable: a false
        // resolution_limit success at |f| = 4.5e307.
        const auto near_max = [](double x) { return x - 1.7e308; };
        const auto window   = nxx::bracket { -1e307, 1e307 };
        const auto found    = nr::expand {}(near_max, window);
        if (found) {
            CHECK(std::isfinite(found->lo()));
            CHECK(std::isfinite(found->hi()));
            CHECK(found->hi() <= (std::numeric_limits<double>::max)());
        }
        else
            FAIL_CHECK("expand failed on x - 1.7e308");
        const auto solved = nxx::then(nr::expand {}.on(window), nr::bisection {})(near_max);
        if (solved) {
            CHECK(solved->how != nxx::stop_reason::resolution_limit);
            CHECK(std::abs(solved->x - 1.7e308) <= 1e-14 * 1.7e308);
        }
        else
            FAIL_CHECK("then(expand, bisection) failed on x - 1.7e308");

        // No root: once both ends are at +-max there is nothing left to search, so expand stops early with stalled.
        const auto below = [](double x) { return 1.0 / (1.0 + x * x) - 2.0; };
        const auto none  = nr::expand {}(below, { -1e300, 1e300 });
        CHECK_FALSE(none.has_value());
        if (!none) {
            CHECK(none.error().code == nxx::errc::stalled);
            CHECK(none.error().used.iterations < 60u);
            CHECK(none.error().best.has_value());
            if (none.error().best) CHECK(std::isinf(none.error().best->uncertainty));
        }
    }

    // ---- 15. Structural derivative ----------------------------------------------------------------------------------

    TEST_CASE("solvers: newton uses a structural derivative without with_derivative")
    {
        static_assert(nxx::has_derivative_source_v<quadratic_poly>);
        static_assert(!nxx::has_derivative_source_v<decltype(sq2)>);
        const auto res = nr::newton {}(quadratic_poly { -2.0 }, 1.0);
        if (res) {
            CHECK(res->by == nr::algos::newton);
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(std::abs(res->x - root2) <= 4.0 * std::numeric_limits<double>::epsilon() * root2);
            CHECK(res->used.evaluations == 2u * res->used.iterations + 1u);
        }
        else
            FAIL_CHECK("newton with a structural derivative failed");
    }

    // ---- 16. Invocability -------------------------------------------------------------------------------------------

    // A root_estimate needs x and f(x): an open method starts from its fx without evaluating f, so a value-initialised
    // fx (root_estimate{1.0}) would pass for an exact zero that was never evaluated.
    static_assert(!std::is_constructible_v<nr::root_estimate<double>, double>);
    static_assert(!std::is_default_constructible_v<nr::root_estimate<double>>);
    static_assert(std::is_constructible_v<nr::root_estimate<double>, double, double>);

    // min_iterations is a guard: every solver's constructor and rebuild reject it alone or under ||, and accept it under
    // && with a convergence test (the with_stop and CTAD paths are compile-fail cases).
    using guard_t    = nxx::min_iterations;
    using guard_or_t = decltype(nxx::x_tol { 1e-12 } || nxx::min_iterations { 2 });
    using guarded_t  = decltype(nxx::x_tol { 1e-12 } && nxx::min_iterations { 2 });
    using width_or_t = decltype(nxx::width_tol { 1e-12 } || nxx::min_iterations { 2 });
    static_assert(!std::is_constructible_v<nr::secant<nxx::options<guard_t>>, guard_t>);
    static_assert(!std::is_constructible_v<nr::secant<nxx::options<guard_or_t>>, guard_or_t>);
    static_assert(std::is_constructible_v<nr::secant<nxx::options<guarded_t>>, guarded_t>);
    static_assert(!std::is_constructible_v<nr::newton<nxx::options<guard_t>>, guard_t>);
    static_assert(!std::is_constructible_v<nr::newton<nxx::options<guard_or_t>>, guard_or_t>);
    static_assert(std::is_constructible_v<nr::newton<nxx::options<guarded_t>>, guarded_t>);
    static_assert(!std::is_constructible_v<nr::bisection<nxx::options<guard_t>>, guard_t>);
    static_assert(!std::is_constructible_v<nr::bisection<nxx::options<width_or_t>>, width_or_t>);
    template<class S, class O>
    concept rebuilds_with = requires(const S& s, O o) { s.rebuild(o); };
    static_assert(!rebuilds_with<nr::secant<>, nxx::options<guard_t>>);
    static_assert(!rebuilds_with<nr::newton<>, nxx::options<guard_or_t>>);
    static_assert(!rebuilds_with<nr::bisection<>, nxx::options<guard_t>>);
    static_assert(!rebuilds_with<nr::brent<>, nxx::options<guard_t>>);
    static_assert(rebuilds_with<nr::secant<>, nxx::options<guarded_t>>);
    static_assert(rebuilds_with<nr::brent<>, nxx::options<nxx::never>>);
    // A searcher has no configurable stop criterion, on the rebuild path too: any other criterion could stop on a state
    // without a sign change, and estimate() would forge a sign_bracket.
    static_assert(!rebuilds_with<nr::expand<>, nxx::options<nxx::x_tol<double>>>);
    static_assert(!rebuilds_with<nr::expand<>, nxx::options<guarded_t>>);
    static_assert(rebuilds_with<nr::expand<>, nxx::options<nxx::never>>);

    // A solver with its own tolerance (brent) takes a width criterion only in its constructor. Its intrinsic test runs
    // before the stop criterion, so a width criterion given to with_stop or rebuild, alone or at any depth of || and
    // &&, reported stop_reason::criterion once brent's own tolerance held: width 6.66e-16 for width_tol{1e-20, 0} on
    // x^2 - 2 over [1, 2] (DESIGN §6.8). with_stop deletes it with a reason and rebuild is constrained, so both are
    // false here rather than hard errors. Solvers without their own tolerance (bisection) keep taking it.
    template<class S, class C>
    constexpr bool with_stop_accepts = requires(const S& s, const C& c) { s.with_stop(c); };

    using wt_t           = nxx::width_tol<double>;
    using fw_t           = nxx::floored_width;
    using ft_t           = nxx::f_tol<double>;
    using me_t           = nxx::max_evaluations;
    using width_or_f_t   = decltype(std::declval<wt_t>() || std::declval<ft_t>());
    using width_and_f_t  = decltype(std::declval<wt_t>() && std::declval<ft_t>());
    using f_and_width_t  = decltype(std::declval<ft_t>() && std::declval<wt_t>());
    using fw_or_f_t      = decltype(std::declval<fw_t>() || std::declval<ft_t>());
    using nested_width_t = decltype((std::declval<ft_t>() || std::declval<me_t>()) &&
                                    (std::declval<me_t>() || (std::declval<ft_t>() && std::declval<fw_t>())));
    using f_or_budget_t  = decltype(std::declval<ft_t>() || std::declval<me_t>());
    using f_guarded_t    = decltype(std::declval<ft_t>() && std::declval<guard_t>());
    static_assert(nxx::detail::contains_width_v<wt_t> && nxx::detail::contains_width_v<fw_t>);
    static_assert(nxx::detail::contains_width_v<width_or_f_t> && nxx::detail::contains_width_v<f_and_width_t>);
    static_assert(nxx::detail::contains_width_v<nested_width_t> && nxx::detail::contains_width_v<const wt_t&>);
    static_assert(!nxx::detail::contains_width_v<ft_t> && !nxx::detail::contains_width_v<f_or_budget_t>);
    static_assert(!nxx::detail::contains_width_v<nxx::never> && !nxx::detail::contains_width_v<nxx::x_tol<double>>);
    static_assert(!nxx::detail::contains_width_v<f_guarded_t> && !nxx::detail::contains_width_v<double>);
    // Rejected on brent, through with_stop and through rebuild, whatever its tolerance.
    static_assert(!with_stop_accepts<nr::brent<>, wt_t>);
    static_assert(!with_stop_accepts<nr::brent<>, fw_t>);
    static_assert(!with_stop_accepts<nr::brent<>, width_or_f_t>);
    static_assert(!with_stop_accepts<nr::brent<>, width_and_f_t>);
    static_assert(!with_stop_accepts<nr::brent<>, f_and_width_t>);
    static_assert(!with_stop_accepts<nr::brent<>, nested_width_t>);
    static_assert(!with_stop_accepts<nr::brent<wt_t>, fw_t>);
    static_assert(!rebuilds_with<nr::brent<>, nxx::options<wt_t>>);
    static_assert(!rebuilds_with<nr::brent<>, nxx::options<width_or_f_t>>);
    static_assert(!rebuilds_with<nr::brent<>, nxx::options<nested_width_t>>);
    static_assert(!rebuilds_with<nr::brent<wt_t>, nxx::options<fw_t>>);
    // The detail key is closed too, and with_stop and every rebuild share one predicate (stop_allowed_v).
    static_assert(!std::is_constructible_v<nr::brent<fw_t, nxx::options<wt_t>>, nxx::detail::from_options_t, nxx::options<wt_t>, fw_t>);
    static_assert(std::is_constructible_v<nr::brent<fw_t, nxx::options<ft_t>>, nxx::detail::from_options_t, nxx::options<ft_t>, fw_t>);
    static_assert(!nxx::detail::stop_allowed_v<nr::brent<>, width_or_f_t> && nxx::detail::stop_allowed_v<nr::brent<>, ft_t>);
    static_assert(!nxx::detail::stop_allowed_v<nr::brent<>, guard_t> && !nxx::detail::stop_allowed_v<nr::bisection<>, guard_t>);
    static_assert(nxx::detail::stop_allowed_v<nr::bisection<>, width_and_f_t>);
    // Kept: brent's early exits and failures, its constructor, and every width criterion on bisection.
    static_assert(with_stop_accepts<nr::brent<>, ft_t>);
    static_assert(with_stop_accepts<nr::brent<>, me_t>);
    static_assert(with_stop_accepts<nr::brent<>, f_or_budget_t>);
    static_assert(with_stop_accepts<nr::brent<>, f_guarded_t>);
    static_assert(with_stop_accepts<nr::brent<wt_t>, ft_t>);
    static_assert(rebuilds_with<nr::brent<>, nxx::options<ft_t>>);
    static_assert(std::is_constructible_v<nr::brent<wt_t>, wt_t>);
    static_assert(std::is_constructible_v<nr::brent<>, fw_t>);
    static_assert(with_stop_accepts<nr::bisection<>, wt_t>);
    static_assert(with_stop_accepts<nr::bisection<>, fw_or_f_t>);
    static_assert(with_stop_accepts<nr::bisection<>, width_and_f_t>);
    static_assert(with_stop_accepts<nr::bisection<>, nested_width_t>);
    static_assert(rebuilds_with<nr::bisection<>, nxx::options<wt_t>>);

    // The route that works: the width criterion as brent's tolerance. Below brent's resolution floor (4 eps |b|) it
    // stops at the floor and reports resolution_limit, because the width does not meet the tolerance; with_stop gave
    // stop_reason::criterion at that same width. The tolerance is the smallest normal T, below the floor for every T: a
    // fixed 1e-20 is above it where long double is binary128 (wasm32, ε = 1.9e-34), and brent then meets it.
    TEST_CASE_TEMPLATE("solvers: brent with a width tolerance below its floor reports resolution_limit", T, float, double, long double)
    {
        constexpr T tiny = (std::numeric_limits<T>::min)();
        const auto  quad = [](T x) { return x * x - T(2); };
        const auto  res  = nr::brent { nxx::width_tol { tiny, T(0) } }(quad, { T(1), T(2) });
        if (res) {
            CHECK(res->how == nxx::stop_reason::resolution_limit);
            if (res->enclosure) {
                const T width = res->enclosure->hi() - res->enclosure->lo();
                CHECK(width > tiny);
                CHECK(width <= T(4) * std::numeric_limits<T>::epsilon() * T(2));
            }
            else
                FAIL_CHECK("brent returned no enclosure");
        }
        else
            FAIL_CHECK("brent failed on x^2 - 2 with a tolerance below its floor");
    }

    // A braced list is a bracket only with two ends of a real type: {x} was once taken as {x, 0} and solved on [0, x].
    template<class S, class F>
    concept takes_one_end = requires(const S& s, const F& fn) { s(fn, { 1.0 }); };
    template<class S, class F>
    concept takes_three_ends = requires(const S& s, const F& fn) { s(fn, { 1.0, 2.0, 3.0 }); };
    template<class S, class F>
    concept takes_int_ends = requires(const S& s, const F& fn) { s(fn, { 1, 2 }); };
    template<class S, class F>
    concept takes_two_ends = requires(const S& s, const F& fn) { s(fn, { 1.0, 2.0 }); };
    template<class S>
    concept binds_one_end = requires(const S& s) { s.on({ 1.0 }); };
    template<class S>
    concept binds_pointer = requires(const S& s, const double* p) { s.on(p); };
    struct line_t
    {
        double operator()(double x) const { return x; }
    };
    static_assert(!takes_one_end<nr::brent<>, line_t> && !takes_one_end<nr::bisection<>, line_t> && !takes_one_end<nr::expand<>, line_t>);
    static_assert(!takes_three_ends<nr::brent<>, line_t> && !takes_int_ends<nr::brent<>, line_t>);
    static_assert(takes_two_ends<nr::brent<>, line_t> && takes_two_ends<nr::expand<>, line_t>);
    static_assert(!binds_one_end<nr::bisection<>> && !binds_one_end<nr::expand<>>);
    static_assert(!binds_pointer<nr::brent<>> && !binds_pointer<nr::expand<>>);

    // Open methods take one guess, not a braced list: a deleted overload gives the reason (compile-fail cases
    // open_braced_bracket and newton_on_braced_bracket), for a named array too, and a pointer keeps "open methods take a
    // guess" (the .on catch-all takes a forwarding reference, so an array never decays to a pointer).
    template<class S>
    concept binds_two_ends = requires(const S& s) { s.on({ 1.0, 2.0 }); };
    template<class S>
    concept binds_array = requires(const S& s, const double (&a)[2]) { s.on(a); };
    template<class S, class F>
    concept takes_array = requires(const S& s, const F& fn, const double (&a)[2]) { s(fn, a); };
    using newton_line_t = decltype(nr::newton {}.with_derivative(line_t {}));
    static_assert(!takes_two_ends<nr::secant<>, line_t> && !takes_one_end<nr::secant<>, line_t> && !takes_array<nr::secant<>, line_t>);
    static_assert(!takes_two_ends<newton_line_t, line_t> && !takes_array<newton_line_t, line_t>);
    static_assert(!binds_two_ends<nr::secant<>> && !binds_array<nr::secant<>> && !binds_pointer<nr::secant<>>);
    static_assert(!binds_two_ends<newton_line_t> && !binds_array<newton_line_t>);
    static_assert(binds_two_ends<nr::brent<>> && binds_array<nr::brent<>>);

    // A bare number is not a tolerance (DESIGN §6.8): each solver deletes it with a reason (compile-fail cases
    // brent_number_tolerance, bisection_number_tolerance and secant_number_tolerance). brent's width_tolerance_v asked
    // W::applies_to of every W, so the brent probes were hard errors inside brent.hpp rather than false.
    template<class T>
    concept brent_from = requires(T t) { nr::brent { t }; };
    template<class T>
    concept bisection_from = requires(T t) { nr::bisection { t }; };
    template<class T>
    concept newton_from = requires(T t) { nr::newton { t }; };
    static_assert(!std::is_constructible_v<nr::brent<double>, double>);
    static_assert(!std::is_constructible_v<nr::brent<>, double> && !std::is_constructible_v<nr::brent<wt_t>, double>);
    static_assert(!brent_from<double> && !brent_from<float> && !brent_from<int>);
    static_assert(brent_from<wt_t> && brent_from<fw_t>);
    static_assert(!bisection_from<double> && !newton_from<double>);
    static_assert(bisection_from<wt_t> && newton_from<nxx::x_tol<double>>);
    // A number that the solver's own criterion type takes implicitly is still accepted: the deletion excludes it, so a
    // user-defined criterion with a converting constructor keeps its run-time spelling, brent<user_width>{tol}.
    struct user_width : nxx::criterion_base
    {
        static constexpr nxx::view_kind applies_to = nxx::view_kind::enclosure;
        double                          w;
        constexpr user_width(double a) noexcept : w(a) {}
        template<class U>
        constexpr U threshold(const U&, const U&) const noexcept
        { return U(w); }
        template<class V>
        constexpr nxx::verdict operator()(const V&, const V& next, nxx::counters) const
        {
            const auto e = next.enclosure();
            return e.hi() - e.lo() <= w ? nxx::verdict::converged : nxx::verdict::proceed;
        }
    };
    static_assert(std::is_constructible_v<nr::brent<user_width>, double>);
    static_assert(std::is_constructible_v<nr::bisection<nxx::options<user_width>>, double>);
    static_assert(brent_from<user_width>);

    // The facades check the solver protocol (DESIGN §6.6), so a solver that lacks part of it makes std::is_invocable_v
    // false, with a reasoned deletion (compile-fail case solver_incomplete), rather than a hard error inside
    // detail::run. Each type below hides one protocol member of bisection<>, which is what a solver without that
    // member looks like to the facade.
    struct hides_accepts : nr::bisection<>
    {
        static constexpr int accepts_v = 0;
    };
    struct hides_prepare : nr::bisection<>
    {
        void prepare() const = delete;
    };
    struct hides_id : nr::bisection<>
    {
        static constexpr int id = 0;
    };
    struct hides_options : nr::bisection<>
    {
        void options() const = delete;
    };
    struct hides_init : nr::bisection<>
    {
        void init() const = delete;
    };
    struct hides_step : nr::bisection<>
    {
        void step() const = delete;
    };
    struct hides_view : nr::bisection<>
    {
        void view() const = delete;
    };
    struct hides_estimate : nr::bisection<>
    {
        void estimate() const = delete;
    };
    struct hides_best : nr::bisection<>
    {
        void best() const = delete;
    };
    struct hides_intrinsic : nr::bisection<>
    {
        void intrinsic() const = delete;
    };
    // A complete solver: bisection's protocol under an id of its own.
    struct own_bisection : nr::bisection<>
    {
        static constexpr nxx::algo id = nxx::algo::user_first;
    };
    template<class S>
    constexpr bool runs_on_pair = std::is_invocable_v<const S&, const line_t&, const std::pair<double, double>&>;
    template<class S>
    constexpr bool runs_on_braced = std::is_invocable_v<const S&, const line_t&, const double (&)[2]>;
    static_assert(runs_on_pair<own_bisection> && runs_on_braced<own_bisection>);
    static_assert(!runs_on_pair<hides_accepts> && !runs_on_pair<hides_prepare> && !runs_on_pair<hides_id> && !runs_on_pair<hides_options>);
    static_assert(!runs_on_pair<hides_init> && !runs_on_pair<hides_step> && !runs_on_pair<hides_view>);
    static_assert(!runs_on_pair<hides_estimate> && !runs_on_pair<hides_best> && !runs_on_pair<hides_intrinsic>);
    // A braced {lo, hi} is converted to a pair and never asks accepts_v, as before; the rest of the protocol is checked.
    static_assert(!runs_on_braced<hides_prepare> && !runs_on_braced<hides_id> && !runs_on_braced<hides_options>);
    static_assert(!runs_on_braced<hides_init> && !runs_on_braced<hides_step> && !runs_on_braced<hides_view>);
    static_assert(!runs_on_braced<hides_estimate> && !runs_on_braced<hides_best> && !runs_on_braced<hides_intrinsic>);
    // prepare must return a std::expected: one that returns a std::optional problem has a value_type but no error, which
    // run needs on its failure path, so the call was a hard error there.
    struct optional_prepare : nr::bisection<>
    {
        template<class F, class In>
        constexpr auto prepare(const F& fn, const In& in) const
        {
            auto p  = nr::bisection<>::prepare(fn, in);
            using P = std::remove_cvref_t<decltype(*p)>;
            return p ? std::optional<P> { *p } : std::optional<P> {};
        }
    };
    static_assert(!runs_on_pair<optional_prepare> && !runs_on_braced<optional_prepare>);
    // The builders read views through detail::stop_allowed_v: a facade-derived type without it makes with_stop false
    // too, not a hard error inside that variable template.
    struct no_views : nxx::bracketing_facade
    {
    };
    template<class S, class C>
    constexpr bool takes_stop = requires(const S& s, C c) { s.with_stop(c); };
    static_assert(!takes_stop<no_views, nxx::width_tol<double>> && !takes_stop<no_views, nxx::f_tol<double>>);
    static_assert(takes_stop<nr::bisection<>, nxx::width_tol<double>> && takes_stop<nr::brent<>, nxx::f_tol<double>>);
#if !defined(_MSC_VER) || defined(__clang__)
    // A views of the wrong kind is false too. cl rejects these types at with_stop's constraint, with or without the
    // check, so they are left out there.
    struct int_views : nxx::bracketing_facade
    {
        static constexpr int views = 0;
    };
    struct member_views : nxx::bracketing_facade
    {
        nxx::view_kind views = nxx::view_kind::enclosure;
    };
    static_assert(!takes_stop<int_views, nxx::f_tol<double>> && !takes_stop<member_views, nxx::f_tol<double>>);
#endif

    TEST_CASE("solvers: a solver that implements the protocol runs through the facade under its own id")
    {
        const auto res = own_bisection {}(sq2, { 1.0, 2.0 });
        const auto ref = nr::bisection {}(sq2, { 1.0, 2.0 });
        if (res && ref) {
            CHECK(res->by == nxx::algo::user_first);
            CHECK(res->x == ref->x);
            CHECK(res->used.iterations == ref->used.iterations);
            CHECK(res->used.evaluations == ref->used.evaluations);
        }
        else
            FAIL_CHECK("bisection's protocol under its own id failed on x^2 - 2");
    }

    // best_x needs a solution with an x; a search result is a sign_bracket, with two ends.
    template<class R>
    concept has_best_x = requires(const R& r) { nxx::best_x(r); };
    static_assert(has_best_x<nxx::result<nr::root_estimate<double>>>);
    static_assert(!has_best_x<std::expected<nxx::solution<nr::sign_bracket<double>>, nxx::failure<nr::root_estimate<double>>>>);

    TEST_CASE("solvers: the facades accept exactly their inputs")
    {
        using quad_t     = decltype(sq2);
        using pair_t     = std::pair<double, double>;
        using newton_df  = decltype(nr::newton {}.with_derivative(dsq2));
        using bracket_t  = nxx::bracket<double>;
        using made_t     = std::expected<bracket_t, nxx::errc>;
        using estimate_t = nr::root_estimate<double>;

        // Bracketing solvers take bracket-like inputs, not a guess.
        constexpr bool bisection_guess   = std::is_invocable_v<const nr::bisection<>&, const quad_t&, const double&>;
        constexpr bool bisection_pair    = std::is_invocable_v<const nr::bisection<>&, const quad_t&, const pair_t&>;
        constexpr bool bisection_bracket = std::is_invocable_v<const nr::bisection<>&, const quad_t&, const bracket_t&>;
        constexpr bool bisection_made    = std::is_invocable_v<const nr::bisection<>&, const quad_t&, const made_t&>;
        constexpr bool bisection_est     = std::is_invocable_v<const nr::bisection<>&, const quad_t&, const estimate_t&>;
        constexpr bool brent_guess       = std::is_invocable_v<const nr::brent<>&, const quad_t&, const double&>;
        static_assert(!bisection_guess);
        static_assert(bisection_pair && bisection_bracket && bisection_made);
        static_assert(!bisection_est);    // a root_estimate only through .from_enclosure() (phase 3)
        static_assert(!brent_guess);

        // Newton needs a derivative source.
        constexpr bool newton_bare       = std::is_invocable_v<const nr::newton<>&, const quad_t&, const double&>;
        constexpr bool newton_with_df    = std::is_invocable_v<const newton_df&, const quad_t&, const double&>;
        constexpr bool newton_structural = std::is_invocable_v<const nr::newton<>&, const quadratic_poly&, const double&>;
        constexpr bool newton_int        = std::is_invocable_v<const newton_df&, const quad_t&, const int&>;
        static_assert(!newton_bare);
        static_assert(newton_with_df);
        static_assert(newton_structural);
        static_assert(!newton_int);

        // Open methods take a guess of a real type or a root estimate, not an int and not a bracket.
        constexpr bool secant_int     = std::is_invocable_v<const nr::secant<>&, const quad_t&, const int&>;
        constexpr bool secant_double  = std::is_invocable_v<const nr::secant<>&, const quad_t&, const double&>;
        constexpr bool secant_pair    = std::is_invocable_v<const nr::secant<>&, const quad_t&, const pair_t&>;
        constexpr bool secant_est     = std::is_invocable_v<const nr::secant<>&, const quad_t&, const estimate_t&>;
        constexpr bool expand_bracket = std::is_invocable_v<const nr::expand<>&, const quad_t&, const bracket_t&>;
        static_assert(!secant_int);
        static_assert(secant_double && secant_est);
        static_assert(!secant_pair);
        static_assert(expand_bracket);

        CHECK_FALSE(bisection_guess);
        CHECK_FALSE(newton_bare);
        CHECK(newton_with_df);
        CHECK_FALSE(secant_int);
    }
}
