// The multiprecision leg. It checks the premise that later phases build on (DESIGN D16): cpp_bin_float_50 is described
// by std::numeric_limits as an inexact, non-integer type, which is what Numerixx's default scalar_traits accept without
// an adapter. This test therefore uses Boost.Multiprecision directly, not <numerixx/adapters/multiprecision.hpp>.
// Phases 1-8 add their multiprecision instantiations here; the spike adds the root solvers (DESIGN §3.5, §7.2).
#include <numerixx/core.hpp>
#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <doctest/doctest.h>

#include <expected>
#include <limits>
#include <optional>
#include <type_traits>

using mp50 = boost::multiprecision::cpp_bin_float_50;

namespace
{
    namespace nr = nxx::roots;

    const auto mp_quad  = [](const mp50& x) -> mp50 { return x * x - 2; };
    const auto mp_dquad = [](const mp50& x) -> mp50 { return 2 * x; };

    const mp50& mp_root2()
    {
        static const mp50 value = sqrt(mp50 { 2 });
        return value;
    }

    const mp50& mp_tolerance()
    {
        static const mp50 value { "1e-45" };
        return value;
    }

    // A validated tolerance and its parts are not criteria (DESIGN §6.6, §7.2), in any scalar type.
    template<class T>
    concept brent_from = requires(T t) { nr::brent { t }; };
    template<class T>
    concept bisection_from = requires(T t) { nr::bisection { t }; };
    template<class T>
    concept secant_from = requires(T t) { nr::secant { t }; };
    template<class T>
    concept newton_from = requires(T t) { nr::newton { t }; };
    template<class S, class C>
    concept with_stop_accepts = requires(const S& s, const C& c) { s.with_stop(c); };
    static_assert(!brent_from<nxx::tolerance<mp50>> && !brent_from<nxx::abs_tolerance<mp50>> && !brent_from<nxx::rel_tolerance<mp50>>);
    static_assert(!bisection_from<nxx::tolerance<mp50>> && !bisection_from<nxx::abs_tolerance<mp50>> &&
                  !bisection_from<nxx::rel_tolerance<mp50>>);
    static_assert(!secant_from<nxx::tolerance<mp50>> && !secant_from<nxx::abs_tolerance<mp50>> && !secant_from<nxx::rel_tolerance<mp50>>);
    static_assert(!newton_from<nxx::tolerance<mp50>> && !newton_from<nxx::abs_tolerance<mp50>> && !newton_from<nxx::rel_tolerance<mp50>>);
    static_assert(std::is_constructible_v<nr::brent<nxx::width_tol<mp50>>, nxx::tolerance<mp50>>);
    static_assert(std::is_constructible_v<nr::bisection<nxx::options<nxx::width_tol<mp50>>>, nxx::tolerance<mp50>>);
    static_assert(std::is_constructible_v<nr::secant<nxx::options<nxx::x_tol<mp50>>>, nxx::tolerance<mp50>>);
    static_assert(std::is_constructible_v<nr::newton<nxx::options<nxx::x_tol<mp50>>>, nxx::tolerance<mp50>>);
    static_assert(!std::is_constructible_v<nr::bisection<nxx::options<nxx::width_tol<mp50>>>, nxx::rel_tolerance<mp50>>);
    static_assert(!std::is_constructible_v<nr::secant<nxx::options<nxx::x_tol<mp50>>>, nxx::abs_tolerance<mp50>>);
    static_assert(!with_stop_accepts<nr::bisection<>, nxx::tolerance<mp50>> && !with_stop_accepts<nr::secant<>, mp50>);
    static_assert(with_stop_accepts<nr::bisection<>, nxx::width_tol<mp50>> && with_stop_accepts<nr::brent<>, nxx::f_tol<mp50>>);
}    // namespace

TEST_SUITE("multiprecision")
{
    TEST_CASE("cpp_bin_float_50 has the numeric_limits the scalar traits rely on")
    {
        static_assert(std::numeric_limits<mp50>::is_specialized);
        static_assert(!std::numeric_limits<mp50>::is_integer);
        static_assert(!std::numeric_limits<mp50>::is_exact);
        static_assert(std::numeric_limits<mp50>::digits10 >= 50);
    }

    TEST_CASE("cpp_bin_float_50 arithmetic is more precise than double")
    {
        const mp50 third  = mp50 { 1 } / 3;
        const mp50 recomb = third * 3 - 1;
        CHECK(abs(recomb) < mp50 { 1e-45 });

        const mp50 root2 = sqrt(mp50 { 2 });
        CHECK(abs(root2 * root2 - 2) < mp50 { 1e-45 });
    }

    TEST_CASE("cpp_bin_float_50 is a real scalar for the root solvers")
    {
        static_assert(nxx::real<mp50>);
        static_assert(nxx::is_real_v<mp50>);
    }

    // DESIGN §3.5, §10.3 phase 1 (A1): every module default is achievable, >= 4 eps. cpp_bin_float_50 is not a literal
    // type, so the static_asserts that tests/roots/test_solvers.cpp makes for float, double and long double cannot be
    // made here: this run-time check reads each solver's own default criterion, calls its factor and thresholds, and
    // computes 4 eps at run time (the variable template four_eps<T> is constexpr and cannot be instantiated for mp50).
    TEST_CASE("cpp_bin_float_50 the defaults are achievable: the solvers' default thresholds >= 4 eps")
    {
        const mp50 four_eps = mp50 { 4 } * std::numeric_limits<mp50>::epsilon();
        const mp50 one { 1 };
        const mp50 zero { 0 };
        const auto bisection_stop = nr::bisection<> {}.options().stop;    // floored_width
        const auto brent_tol      = nr::brent<> {}.tolerance();           // floored_width
        const auto newton_stop    = nr::newton<> {}.options().stop;       // step_tol<3, 5>
        const auto secant_stop    = nr::secant<> {}.options().stop;       // step_tol<7, 10>
        CHECK(bisection_stop.factor<mp50>() >= four_eps);
        CHECK(brent_tol.factor<mp50>() >= four_eps);
        CHECK(bisection_stop.threshold(one) >= four_eps);
        CHECK(bisection_stop.threshold(zero) >= four_eps);    // the absolute floor at 0
        CHECK(brent_tol.threshold(zero) >= four_eps);
        CHECK(newton_stop.threshold(one) >= four_eps);
        CHECK(secant_stop.threshold(one) >= four_eps);
        CHECK(newton_stop.threshold(zero) >= four_eps);
        CHECK(secant_stop.threshold(zero) >= four_eps);

        // Relative above 1: still >= 4 eps |x|.
        const mp50 big { 1e6 };
        CHECK(bisection_stop.threshold(big) >= four_eps * big);
        CHECK(newton_stop.threshold(big) >= four_eps * big);
        CHECK(secant_stop.threshold(big) >= four_eps * big);

        // The comparison can fail: at p = 168 digits, step_tol<99, 100> asks for 2^-ceil(0.99 p) = 2^-167, below
        // 4 eps = 2^(3 - p) = 2^-165, and so does a threshold of 2 eps, which the integer check on digits approved on
        // 2026-10-04 accepted (DESIGN §12.21 item 11).
        static_assert(std::numeric_limits<mp50>::digits == 168);
        CHECK_FALSE(nxx::step_tol<99, 100>::threshold(one) >= four_eps);
        CHECK_FALSE(mp50 { 2 } * std::numeric_limits<mp50>::epsilon() >= four_eps);
    }

    TEST_CASE("cpp_bin_float_50 bisection converges to the square root of 2")
    {
        // floored_width{} at 168 bits needs 165 halvings of [1, 2]: the default budget (200) is sized for it, so the
        // defaults are achievable in this type too (DESIGN §3.5).
        const auto res = nr::bisection {}(mp_quad, { mp50 { 1 }, mp50 { 2 } });
        if (res) {
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(res->x)>, mp50>);
            CHECK(res->by == nr::algos::bisection);
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(abs(res->x - mp_root2()) < mp_tolerance());
            CHECK(res->used.evaluations == res->used.iterations + 2u);
        }
        else
            FAIL_CHECK("bisection failed on x^2 - 2 in cpp_bin_float_50");
    }

    TEST_CASE("cpp_bin_float_50 brent converges to the square root of 2")
    {
        const auto res = nr::brent {}(mp_quad, { mp50 { 1 }, mp50 { 2 } });
        if (res) {
            CHECK(res->by == nr::algos::brent);
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(abs(res->x - mp_root2()) < mp_tolerance());
        }
        else
            FAIL_CHECK("brent failed on x^2 - 2 in cpp_bin_float_50");
    }

    // Role-typed mixed criteria at run time (DESIGN §6.2): cpp_bin_float_50 is not a literal type, so its paths are
    // make() and the constexpr constructor from a validated tolerance; the relative part is a rel_tolerance<mp50>, and
    // make(T, T) is deleted.
    TEST_CASE("cpp_bin_float_50 mixed width tolerance through make(abs, rel_tolerance)")
    {
        const mp50 abs_part { "1e-40" };
        const mp50 rel_part { "1e-40" };
        const auto R = nxx::rel_tolerance<mp50>::make(rel_part);
        CHECK(R.has_value());
        if (R) {
            const auto tol = nxx::width_tol<mp50>::make(abs_part, *R);
            CHECK(tol.has_value());
            if (tol) {
                CHECK(tol->abs() == abs_part);
                CHECK(tol->rel() == rel_part);
                const auto res = nr::brent { *tol }(mp_quad, { mp50 { 1 }, mp50 { 2 } });
                CHECK(res.has_value());
                if (res) {
                    CHECK(res->how == nxx::stop_reason::criterion);
                    CHECK(abs(res->x - mp_root2()) <= abs_part + rel_part * mp_root2());
                }
            }
            CHECK(nxx::width_tol<mp50>::make(-abs_part, *R) == std::unexpected(nxx::errc::invalid_input));
            CHECK(nxx::x_tol<mp50>::make(mp50 { 0 }, *R).has_value());    // purely relative
        }
        CHECK(nxx::width_tol<mp50>::make(abs_part).has_value());
        CHECK(nxx::x_tol<mp50>::make(mp50 { 0 }) == std::unexpected(nxx::errc::invalid_input));
    }

    // The remedy of brent's reason for a validated tolerance, brent{nxx::width_tol{*tol}}, at run time.
    TEST_CASE("cpp_bin_float_50 brent with a wrapped validated tolerance")
    {
        const auto tol = nxx::tolerance<mp50>::make(mp50 { "1e-40" });
        CHECK(tol.has_value());
        if (tol) {
            const auto res = nr::brent { nxx::width_tol { *tol } }(mp_quad, { mp50 { 1 }, mp50 { 2 } });
            CHECK(res.has_value());
            if (res) {
                CHECK(res->how == nxx::stop_reason::criterion);
                CHECK(abs(res->x - mp_root2()) <= tol->value());
            }
        }
    }

    // A validated tolerance with a relative part (DESIGN §6.2, §12.24): the constexpr constructor takes run-time
    // values, and make(tolerance, rel_tolerance) cannot fail.
    TEST_CASE("cpp_bin_float_50 a validated tolerance with a relative part, through the constructor and make")
    {
        static_assert(std::is_constructible_v<nxx::width_tol<mp50>, nxx::tolerance<mp50>, nxx::rel_tolerance<mp50>>);
        static_assert(std::is_constructible_v<nxx::x_tol<mp50>, nxx::tolerance<mp50>, nxx::rel_tolerance<mp50>>);
        static_assert(!std::is_constructible_v<nxx::width_tol<mp50>, nxx::tolerance<mp50>, mp50>);
        const auto T = nxx::tolerance<mp50>::make(mp50 { "1e-40" });
        const auto R = nxx::rel_tolerance<mp50>::make(mp50 { "1e-40" });
        CHECK((T && R));
        if (T && R) {
            const auto w = nxx::width_tol { *T, *R };
            static_assert(std::is_same_v<decltype(w), const nxx::width_tol<mp50>>);
            CHECK((w.abs() == T->value() && w.rel() == R->value()));
            const auto m = nxx::x_tol<mp50>::make(*T, *R);
            CHECK(m.has_value());
            if (m) { CHECK((m->abs() == T->value() && m->rel() == R->value())); }
            const auto res = nr::brent { w }(mp_quad, { mp50 { 1 }, mp50 { 2 } });
            CHECK(res.has_value());
            if (res) {
                CHECK(res->how == nxx::stop_reason::criterion);
                CHECK(abs(res->x - mp_root2()) <= T->value() + R->value() * mp_root2());
            }
        }
    }

    TEST_CASE("cpp_bin_float_50 newton converges to the square root of 2")
    {
        const auto res = nr::newton {}.with_derivative(mp_dquad)(mp_quad, mp50 { 1 });
        if (res) {
            CHECK(res->by == nr::algos::newton);
            CHECK((res->how == nxx::stop_reason::criterion || res->how == nxx::stop_reason::exact_zero));
            CHECK(abs(res->x - mp_root2()) < mp_tolerance());
            CHECK(res->used.evaluations == 2u * res->used.iterations + 1u);
        }
        else
            FAIL_CHECK("newton failed on x^2 - 2 in cpp_bin_float_50");
    }

    TEST_CASE("cpp_bin_float_50 secant converges to the square root of 2")
    {
        const auto res = nr::secant {}(mp_quad, mp50 { 1 });
        if (res) {
            CHECK(res->by == nr::algos::secant);
            CHECK((res->how == nxx::stop_reason::criterion || res->how == nxx::stop_reason::exact_zero));
            CHECK(abs(res->x - mp_root2()) < mp_tolerance());
        }
        else
            FAIL_CHECK("secant failed on x^2 - 2 in cpp_bin_float_50");
    }

    TEST_CASE("cpp_bin_float_50 non-finite inputs are non_finite_input, equal ends invalid_input")
    {
        // DESIGN §6.3: the input codes do not depend on the scalar type.
        const mp50 nan = std::numeric_limits<mp50>::quiet_NaN();
        const mp50 inf = std::numeric_limits<mp50>::infinity();
        CHECK(nxx::bracket<mp50>::make(nan, mp50 { 1 }) == std::unexpected(nxx::errc::non_finite_input));
        CHECK(nxx::bracket<mp50>::make(mp50 { 1 }, -inf) == std::unexpected(nxx::errc::non_finite_input));
        CHECK(nxx::bracket<mp50>::make(mp50 { 1 }, mp50 { 1 }) == std::unexpected(nxx::errc::invalid_input));

        const auto res = nr::bisection {}(mp_quad, { nan, mp50 { 2 } });
        CHECK_FALSE(res.has_value());
        if (!res) {
            CHECK(res.error().code == nxx::errc::non_finite_input);
            CHECK(res.error().used == nxx::counters {});
            CHECK_FALSE(res.error().best.has_value());
        }

        const auto dx = nxx::deriv::diff(mp_quad, inf);
        CHECK_FALSE(dx.has_value());
        if (!dx) {
            CHECK(dx.error().code == nxx::errc::non_finite_input);
            CHECK(dx.error().evaluations == 0u);
        }
    }

    TEST_CASE("cpp_bin_float_50 failure estimates: the R3 order, and a pole failure carries no enclosure")
    {
        // DESIGN §6.7, §7.2: an enclosure first; the smaller width, the smaller hi/2 - lo/2 when both widths overflow;
        // then the smaller |f(x)|, with a NaN |f(x)| last.
        using est      = nr::root_estimate<mp50>;
        const mp50 nan = std::numeric_limits<mp50>::quiet_NaN();
        const mp50 inf = std::numeric_limits<mp50>::infinity();
        const mp50 m   = (std::numeric_limits<mp50>::max)();
        const auto enc = [](const mp50& lo, const mp50& hi, const mp50& fx) {
            return est { lo, fx, mp50 { hi - lo }, nr::sign_bracket<mp50> { nxx::detail::trust_me {}, lo, mp50 { -1 }, hi, mp50 { 1 } } };
        };
        const est open_one { mp50 { 0 }, mp50 { 1 } };
        const est open_nan { mp50 { 0 }, nan };
        CHECK(nxx::better_than(open_one, open_nan));
        CHECK_FALSE(nxx::better_than(open_nan, open_one));
        CHECK_FALSE(nxx::better_than(open_nan, open_nan));

        const est whole = enc(-m, m, mp50 { 0 });
        const est most  = enc(mp50 { -m / 2 }, m, mp50 { 1 });
        const est half  = enc(mp50 { 0 }, m, mp50 { 1 });
        CHECK(whole.enclosure->width() == inf);
        CHECK(most.enclosure->width() == inf);
        CHECK(nxx::better_than(most, whole));
        CHECK_FALSE(nxx::better_than(whole, most));
        CHECK(nxx::better_than(half, most));
        CHECK_FALSE(nxx::better_than(most, half));
        CHECK(nxx::better_than(enc(mp50 { 1 }, mp50 { 2 }, nan), open_one));

        const auto hyperbola = [](const mp50& x) -> mp50 { return 1 / (x - mp50 { 1 } / 3); };
        const auto res       = nr::bisection {}(hyperbola, { mp50 { 0 }, mp50 { 1 } });
        CHECK_FALSE(res.has_value());
        if (!res && res.error().best) {
            CHECK(res.error().code == nxx::errc::sign_change_not_root);
            CHECK_FALSE(res.error().best->enclosure.has_value());
            CHECK(abs(res.error().best->x - mp50 { 1 } / 3) <= res.error().best->uncertainty);
        }
        else
            FAIL_CHECK("bisection on a pole in cpp_bin_float_50 fails with a best estimate");
    }
}
