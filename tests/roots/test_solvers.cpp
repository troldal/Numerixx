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
            const auto fixed    = solver(sq2, nxx::bracket { 1.0, 2.0 });    // the literal, checked at compile time
            const auto braced   = solver(sq2, { lo, hi });
            const auto paired   = solver(sq2, std::pair { lo, hi });
            const auto made     = solver(sq2, nxx::bracket<double>::make(lo, hi));
            const auto reversed = solver(sq2, std::pair { hi, lo });    // make() re-orders
            const auto made_rev = solver(sq2, nxx::bracket<double>::make(hi, lo));
            const auto curried  = solver.on({ lo, hi })(sq2);
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
                CHECK(err.best->x > 0.0);    // lo moved towards 0 by lo / 1.6, keeping the smaller |f|
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
