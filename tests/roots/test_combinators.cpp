// Combinators (DESIGN §6.10, §9.3): first_of, first_of_with and is_fatal, then, warm_fallback over the spike's root
// solvers. Laziness (a later alternative or stage never runs when it is not needed), the merge of failures (the last
// code and cause, the better best estimate, the total cost), the cost of every stage, seeding without re-evaluation,
// and constexpr evaluation of the headline chain.

#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <expected>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

// Some checks below are bit for bit. Clang would contract this file's x * x - 2.0 into a fused multiply-add where it
// constant-folds and not elsewhere (the library's own code has contraction off), so it is switched off here too.
#if defined(__clang__)
#    pragma clang fp contract(off)
#endif

namespace
{
    namespace r = nxx::roots;
    using est_t = r::root_estimate<double>;

    constexpr auto sq2  = [](double x) { return x * x - 2.0; };    // root sqrt(2)
    constexpr auto dsq2 = [](double x) { return 2.0 * x; };

    namespace comb_user
    {
        enum class eval_error { domain, fatal };

        // Found by ADL from nxx::is_fatal: only `fatal` stops a first_of chain.
        constexpr bool is_fatal(eval_error err) noexcept { return err == eval_error::fatal; }

        enum class plain_error { oops };    // no is_fatal: never fatal
    }    // namespace comb_user

    // x^2 - 2 for x >= 0; the given user error for x < 0.
    constexpr auto g_domain = [](double x) -> std::expected<double, comb_user::eval_error> {
        if (x < 0.0) return std::unexpected(comb_user::eval_error::domain);
        return x * x - 2.0;
    };
    constexpr auto g_fatal = [](double x) -> std::expected<double, comb_user::eval_error> {
        if (x < 0.0) return std::unexpected(comb_user::eval_error::fatal);
        return x * x - 2.0;
    };

    // ---- Alternatives with known failure modes on x^2 - 2 ------------------------------------------------------------
    // f'(0) = 0: zero_derivative after one iteration, best {0, -2}, no enclosure.
    constexpr auto zero_deriv = r::newton {}.with_derivative(dsq2).on(0.0);
    // Two secant steps from 1: budget_exhausted near the root (|f| about 0.04), no enclosure.
    constexpr auto starved_secant = r::secant {}.with_budget(2).on(1.0);
    // Three secant steps from 1: budget_exhausted closer to the root (|f| about 1e-3), no enclosure.
    constexpr auto starved_secant3 = r::secant {}.with_budget(3).on(1.0);
    // Six halvings of [0, 2]: budget_exhausted, enclosure width 2/64, |f| about 0.02.
    constexpr auto narrow = r::bisection { nxx::never {} }.with_budget(6).on(nxx::bracket { 0.0, 2.0 });
    // One halving of [1.4142135, 3]: budget_exhausted, enclosure width about 0.79, |f| about 2e-7.
    constexpr auto wide = r::bisection { nxx::never {} }.with_budget(1).on(nxx::bracket { 1.4142135, 3.0 });
    // Equal endpoints: invalid_input before anything is evaluated, no best estimate.
    constexpr auto invalid = r::bisection {}.on(std::pair { 1.0, 1.0 });

    // The headline chain (DESIGN §6.11): Newton fails (f'(0) = 0), the secant fails within 5 iterations, bisection
    // succeeds. It runs at compile time.
    constexpr auto headline = nxx::first_of(
        r::newton {}.with_derivative(dsq2).on(0.0), r::secant {}.with_budget(5).on(0.0), r::bisection {}.on(nxx::bracket { 0.0, 2.0 }));
    static_assert(headline(sq2).has_value() && headline(sq2)->by == r::algos::bisection);

    // Staging and warm restarts at compile time too.
    constexpr auto pipeline =
        nxx::then(r::expand {}.on(nxx::bracket { 2.0, 2.5 }), r::bisection { nxx::width_tol { 1e-4 } }, r::newton {}.with_derivative(dsq2));
    static_assert(pipeline(sq2).has_value() && pipeline(sq2)->by == r::algos::newton);

    constexpr auto warm =
        nxx::warm_fallback(r::bisection { nxx::never {} }.with_budget(3).on(nxx::bracket { 0.0, 2.0 }), r::newton {}.with_derivative(dsq2));
    static_assert(warm(sq2).has_value() && warm(sq2)->by == r::algos::newton);

    // ---- Helpers -----------------------------------------------------------------------------------------------------
    bool same_bits(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }

    bool same_estimate(const est_t& a, const est_t& b)
    { return same_bits(a.x, b.x) && same_bits(a.fx, b.fx) && same_bits(a.uncertainty, b.uncertainty) && a.enclosure == b.enclosure; }

    bool same_best(const std::optional<est_t>& a, const std::optional<est_t>& b)
    {
        if (a.has_value() != b.has_value()) return false;
        return !a || same_estimate(*a, *b);
    }

    // Whole results, bit for bit: the estimate, the cost, the algorithm and the stop reason, or the whole failure.
    template<class R>
    bool same_result(const R& a, const R& b)
    {
        if (a.has_value() != b.has_value()) return false;
        if (a) return same_estimate(*a, *b) && a->used == b->used && a->by == b->by && a->how == b->how;
        const auto& fa = a.error();
        const auto& fb = b.error();
        return fa.code == fb.code && fa.by == fb.by && fa.used == fb.used && same_best(fa.best, fb.best) && fa.cause == fb.cause;
    }

    // A curried solver or a stage that counts how often it runs.
    template<class S>
    class spy
    {
        S              solver_;
        std::uint32_t* runs_;

    public:
        spy(S solver, std::uint32_t& runs) : solver_(std::move(solver)), runs_(&runs) {}

        template<class... A>
            requires std::is_invocable_v<const S&, const A&...>
        auto operator()(const A&... args) const
        {
            ++*runs_;
            return std::invoke(solver_, args...);
        }
    };

    // A curried solver that counts its own calls of f, in its own counter.
    template<class S>
    class counting
    {
        S              solver_;
        std::uint32_t* calls_;

    public:
        counting(S solver, std::uint32_t& calls) : solver_(std::move(solver)), calls_(&calls) {}

        template<class F>
        auto operator()(const F& fn) const
        { return solver_(nxx::fn::counted(std::cref(fn), *calls_)); }
    };

    struct stop_on_any_failure
    {
        template<class Est, class UE>
        constexpr bool operator()(const nxx::failure<Est, UE>&) const noexcept
        { return false; }
    };

    struct stop_on_input_error
    {
        template<class Est, class UE>
        constexpr bool operator()(const nxx::failure<Est, UE>& err) const noexcept
        { return !nxx::is_input_error(err.code); }
    };
}    // namespace

TEST_SUITE("roots")
{
    // ---- first_of: laziness ------------------------------------------------------------------------------------------
    TEST_CASE("combinators: first_of returns the first success and never runs a later alternative")
    {
        std::uint32_t runs1 = 0;
        std::uint32_t runs2 = 0;
        const auto    s1    = r::brent {}.on(nxx::bracket { 1.0, 2.0 });
        const auto    chain = nxx::first_of(spy { s1, runs1 }, spy { r::bisection {}.on(nxx::bracket { 1.0, 2.0 }), runs2 });
        const auto    got   = chain(sq2);
        CHECK(same_result(got, s1(sq2)));
        CHECK(runs1 == 1u);
        CHECK(runs2 == 0u);
        if (got) { CHECK(got->by == r::algos::brent); }
        else {
            FAIL_CHECK("brent failed on x2 - 2 over [1, 2]");
        }
    }

    TEST_CASE("combinators: first_of equals its first alternative whenever that one succeeds (property)")
    {
        std::mt19937                           gen(20260928u);
        std::uniform_real_distribution<double> coef(0.5, 3.5);    // roots sqrt(c) in (0.7, 1.9), inside [0, 2]
        std::uint32_t                          runs2 = 0;
        for (int i = 0; i < 100; ++i) {
            const double c     = coef(gen);
            const auto   fc    = [c](double x) { return x * x - c; };
            const auto   s1    = r::brent {}.on(nxx::bracket { 0.0, 2.0 });
            const auto   chain = nxx::first_of(s1, spy { r::bisection {}.on(nxx::bracket { 0.0, 2.0 }), runs2 });
            const auto   alone = s1(fc);
            CHECK(alone.has_value());
            CHECK(same_result(chain(fc), alone));
        }
        CHECK(runs2 == 0u);
    }

    // ---- first_of: the merge of failures -----------------------------------------------------------------------------
    TEST_CASE("combinators: first_of on total failure reports the last code and the better best estimate")
    {
        const auto e_sec = starved_secant(sq2);
        const auto e_zd  = zero_deriv(sq2);
        if (e_sec || e_zd || !e_sec.error().best || !e_zd.error().best) {
            FAIL_CHECK("precondition: both alternatives fail with a best estimate");
            return;
        }
        CHECK(e_sec.error().code == nxx::errc::budget_exhausted);
        CHECK(e_zd.error().code == nxx::errc::zero_derivative);
        CHECK(e_zd.error().used == nxx::counters { 1, 2 });    // f(0), then f'(0) = 0
        CHECK(std::abs(e_sec.error().best->fx) < std::abs(e_zd.error().best->fx));

        SUBCASE("the earlier attempt has the better estimate")
        {
            const auto m = nxx::first_of(starved_secant, zero_deriv)(sq2);
            if (m) {
                FAIL_CHECK("the chain cannot succeed");
                return;
            }
            CHECK(m.error().code == nxx::errc::zero_derivative);    // the last code
            CHECK(m.error().by == r::algos::newton);
            CHECK(m.error().used == e_sec.error().used + e_zd.error().used);    // the total cost
            CHECK(same_best(m.error().best, e_sec.error().best));               // the better estimate, not the last
        }
        SUBCASE("the later attempt has the better estimate")
        {
            const auto m = nxx::first_of(zero_deriv, starved_secant)(sq2);
            if (m) {
                FAIL_CHECK("the chain cannot succeed");
                return;
            }
            CHECK(m.error().code == nxx::errc::budget_exhausted);
            CHECK(m.error().by == r::algos::secant);
            CHECK(m.error().used == e_zd.error().used + e_sec.error().used);
            CHECK(same_best(m.error().best, e_sec.error().best));
        }
    }

    TEST_CASE("combinators: first_of merges enclosure-aware - the narrower enclosure wins over the smaller residual")
    {
        const auto e_narrow = narrow(sq2);
        const auto e_wide   = wide(sq2);
        if (e_narrow || e_wide || !e_narrow.error().best || !e_wide.error().best) {
            FAIL_CHECK("precondition: both bisections are starved and fail with a best estimate");
            return;
        }
        const est_t& b_narrow = *e_narrow.error().best;
        const est_t& b_wide   = *e_wide.error().best;
        if (!b_narrow.enclosure || !b_wide.enclosure) {
            FAIL_CHECK("precondition: a starved bisection's best estimate carries its enclosure");
            return;
        }
        CHECK(e_narrow.error().code == nxx::errc::budget_exhausted);
        CHECK(e_wide.error().code == nxx::errc::budget_exhausted);
        CHECK(b_narrow.enclosure->width() == 2.0 / 64.0);
        CHECK(b_narrow.enclosure->width() < b_wide.enclosure->width());
        CHECK(std::abs(b_wide.fx) < std::abs(b_narrow.fx));    // so the residual rule alone would pick the wide one
        CHECK(nxx::better_than(b_narrow, b_wide));
        CHECK_FALSE(nxx::better_than(b_wide, b_narrow));

        for (const auto& m : { nxx::first_of(narrow, wide)(sq2), nxx::first_of(wide, narrow)(sq2) }) {
            if (m) {
                FAIL_CHECK("the chain cannot succeed");
                continue;
            }
            CHECK(m.error().code == nxx::errc::budget_exhausted);
            CHECK(m.error().by == r::algos::bisection);
            CHECK(m.error().used == e_narrow.error().used + e_wide.error().used);
            CHECK(same_best(m.error().best, e_narrow.error().best));
        }
    }

    // better_than is a strict weak order (roots/bracket.hpp, DESIGN §6.7): an estimate with a sign-changing enclosure beats one without,
    // whatever the residuals, so the best estimate of a chain does not depend on grouping or fold order.
    TEST_CASE("combinators: first_of prefers the estimate with an enclosure when only one has one")
    {
        const auto e_narrow = narrow(sq2);
        const auto e_wide   = wide(sq2);
        const auto e_zd     = zero_deriv(sq2);
        const auto e_sec3   = starved_secant3(sq2);
        if (e_narrow || e_wide || e_zd || e_sec3 || !e_narrow.error().best || !e_wide.error().best || !e_zd.error().best ||
            !e_sec3.error().best)
        {
            FAIL_CHECK("precondition: every alternative fails with a best estimate");
            return;
        }
        CHECK_FALSE(e_zd.error().best->enclosure.has_value());
        CHECK_FALSE(e_sec3.error().best->enclosure.has_value());

        SUBCASE("the estimate with the enclosure has the smaller residual")
        {
            CHECK(std::abs(e_wide.error().best->fx) < std::abs(e_zd.error().best->fx));
            for (const auto& m : { nxx::first_of(wide, zero_deriv)(sq2), nxx::first_of(zero_deriv, wide)(sq2) }) {
                if (m) {
                    FAIL_CHECK("the chain cannot succeed");
                    continue;
                }
                CHECK(same_best(m.error().best, e_wide.error().best));
                CHECK(m.error().used == e_wide.error().used + e_zd.error().used);
            }
        }
        SUBCASE("the estimate without an enclosure has the smaller residual: the enclosure still wins")
        {
            CHECK(std::abs(e_sec3.error().best->fx) < std::abs(e_narrow.error().best->fx));
            for (const auto& m : { nxx::first_of(narrow, starved_secant3)(sq2), nxx::first_of(starved_secant3, narrow)(sq2) }) {
                if (m) {
                    FAIL_CHECK("the chain cannot succeed");
                    continue;
                }
                CHECK(same_best(m.error().best, e_narrow.error().best));
                CHECK(m.error().used == e_narrow.error().used + e_sec3.error().used);
            }
        }
    }

    // The enclosure rule's premise is that an enclosure holds a sign change of a root (DESIGN §6.7, §7.2): a pole failure
    // carries its estimate without the enclosure, so a chain does not rank the pole above an open method's estimate.
    // Before, brent's final enclosure around pi/2 outranked the secant's estimate, and the chain's best was the pole
    // (x = 1.5707963267948974, |f| about 1.2e15).
    TEST_CASE("combinators: first_of does not rank a pole failure's estimate as the best")
    {
        const auto tangent  = [](double x) { return std::tan(x); };
        const auto at_pole  = r::brent {}.on({ 1.0, 2.0 });
        const auto starved  = r::secant {}.with_budget(2).on(3.0);
        const auto e_pole   = at_pole(tangent);
        const auto e_secant = starved(tangent);
        if (e_pole || e_secant || !e_pole.error().best || !e_secant.error().best) {
            FAIL_CHECK("precondition: brent finds the pole and the starved secant fails, both with a best estimate");
            return;
        }
        CHECK(e_pole.error().code == nxx::errc::sign_change_not_root);
        CHECK_FALSE(e_pole.error().best->enclosure.has_value());
        CHECK(std::abs(e_pole.error().best->x - std::numbers::pi / 2.0) < 1e-12);
        CHECK(std::abs(e_pole.error().best->fx) > 1e10);
        CHECK(e_secant.error().code == nxx::errc::budget_exhausted);
        CHECK(std::abs(e_secant.error().best->x - std::numbers::pi) < 1e-4);
        CHECK(std::abs(e_secant.error().best->fx) < std::abs(e_pole.error().best->fx));

        for (const auto& m : { nxx::first_of(at_pole, starved)(tangent), nxx::first_of(starved, at_pole)(tangent) }) {
            if (m) {
                FAIL_CHECK("the chain cannot succeed");
                continue;
            }
            CHECK(same_best(m.error().best, e_secant.error().best));
            CHECK(m.error().used == e_pole.error().used + e_secant.error().used);
        }
    }

    // The same rule through then (DESIGN §6.7, §6.10): when stage 2 finds a pole in the search's bracket, the chain keeps
    // stage 2's pole estimate. Before, then merged the search's estimate, a window end whose enclosure [1, 2] holds the
    // pole, and it outranked the secant's estimate in first_of (best x = 1, |f| about 1.56, in both orders).
    TEST_CASE("combinators: a then whose stage 2 finds a pole keeps the pole estimate without an enclosure (regression)")
    {
        const auto tangent = [](double x) { return std::tan(x); };
        const auto search  = r::expand {}.on(nxx::bracket { 1.0, 2.0 });

        const auto check_chain = [&](const auto& alone, const auto& chain) {
            const auto e_alone = alone(tangent);
            const auto got     = chain(tangent);
            if (e_alone || got || !e_alone.error().best) {
                FAIL_CHECK("precondition: the solver alone and the chain fail, the solver with a best estimate");
                return;
            }
            CHECK(e_alone.error().code == nxx::errc::sign_change_not_root);
            CHECK(got.error().code == nxx::errc::sign_change_not_root);
            CHECK(got.error().used == e_alone.error().used);
            CHECK(same_best(got.error().best, e_alone.error().best));
        };
        check_chain(r::brent {}.on({ 1.0, 2.0 }), nxx::then(search, r::brent {}));
        check_chain(r::bisection {}.on({ 1.0, 2.0 }), nxx::then(search, r::bisection {}));

        // The symptom as first seen: a first_of of the chain and a starved secant reported the search's window end.
        const auto starved  = r::secant {}.with_budget(2).on(3.0);
        const auto e_secant = starved(tangent);
        const auto m        = nxx::first_of(nxx::then(search, r::brent {}), starved)(tangent);
        CHECK((!m && !e_secant && same_best(m.error().best, e_secant.error().best)));
    }

    TEST_CASE("combinators: first_of keeps a best estimate when only one attempt has one")
    {
        const auto e_inv = invalid(sq2);
        const auto e_zd  = zero_deriv(sq2);
        CHECK(!e_inv.has_value());
        if (!e_inv) {
            CHECK(e_inv.error().code == nxx::errc::invalid_input);
            CHECK(e_inv.error().used == nxx::counters {});
            CHECK_FALSE(e_inv.error().best.has_value());
        }
        if (e_zd) {
            FAIL_CHECK("newton from 0 cannot succeed");
            return;
        }

        const auto m1 = nxx::first_of(invalid, zero_deriv)(sq2);
        if (m1) { FAIL_CHECK("the chain cannot succeed"); }
        else {
            CHECK(m1.error().code == nxx::errc::zero_derivative);
            CHECK(same_best(m1.error().best, e_zd.error().best));
            CHECK(m1.error().used == e_zd.error().used);
        }

        const auto m2 = nxx::first_of(zero_deriv, invalid)(sq2);
        if (m2) { FAIL_CHECK("the chain cannot succeed"); }
        else {
            CHECK(m2.error().code == nxx::errc::invalid_input);    // the last code
            CHECK(m2.error().by == r::algos::bisection);
            CHECK(same_best(m2.error().best, e_zd.error().best));    // the earlier estimate survives
            CHECK(m2.error().used == e_zd.error().used);
        }
    }

    TEST_CASE("combinators: first_of reports the last attempt's cause")
    {
        const auto bis   = r::bisection {}.on({ -1.0, 2.0 });      // g fails at -1: callback_failed, cause domain
        const auto sec   = r::secant {}.with_budget(1).on(1.0);    // budget_exhausted, no cause
        const auto e_bis = bis(g_domain);
        const auto e_sec = sec(g_domain);
        if (e_bis || e_sec) {
            FAIL_CHECK("precondition: both alternatives fail");
            return;
        }
        CHECK(e_bis.error().code == nxx::errc::callback_failed);
        CHECK(e_bis.error().cause == std::optional { comb_user::eval_error::domain });
        CHECK_FALSE(e_bis.error().best.has_value());
        CHECK(e_sec.error().code == nxx::errc::budget_exhausted);
        CHECK_FALSE(e_sec.error().cause.has_value());

        const auto m1 = nxx::first_of(bis, sec)(g_domain);    // domain is not fatal: the secant runs
        if (m1) { FAIL_CHECK("the chain cannot succeed"); }
        else {
            CHECK(m1.error().code == nxx::errc::budget_exhausted);
            CHECK_FALSE(m1.error().cause.has_value());    // the last cause, not the earlier one
            CHECK(same_best(m1.error().best, e_sec.error().best));
            CHECK(m1.error().used == e_bis.error().used + e_sec.error().used);
        }

        const auto m2 = nxx::first_of(sec, bis)(g_domain);
        if (m2) { FAIL_CHECK("the chain cannot succeed"); }
        else {
            CHECK(m2.error().code == nxx::errc::callback_failed);
            CHECK(m2.error().cause == std::optional { comb_user::eval_error::domain });
            CHECK(same_best(m2.error().best, e_sec.error().best));
            CHECK(m2.error().used == e_sec.error().used + e_bis.error().used);
        }
    }

    TEST_CASE("combinators: first_of success pays for the failed attempts and counts every call of f")
    {
        std::uint32_t c1    = 0;
        std::uint32_t c2    = 0;
        std::uint32_t c3    = 0;
        const auto    good  = r::brent {}.on(nxx::bracket { 1.0, 2.0 });
        const auto    chain = nxx::first_of(counting { starved_secant, c1 }, counting { invalid, c2 }, counting { good, c3 });
        const auto    got   = chain(sq2);
        const auto    e1    = starved_secant(sq2);
        const auto    e2    = invalid(sq2);
        const auto    s3    = good(sq2);
        if (!got || e1 || e2 || !s3) {
            FAIL_CHECK("precondition: the secant and the invalid bisection fail, brent succeeds");
            return;
        }
        CHECK(got->by == r::algos::brent);
        CHECK(same_bits(got->x, s3->x));
        CHECK(got->used == e1.error().used + e2.error().used + s3->used);
        CHECK(c1 == e1.error().used.evaluations);
        CHECK(c2 == 0u);
        CHECK(c3 == s3->used.evaluations);
        CHECK(got->used.evaluations == c1 + c2 + c3);
    }

    TEST_CASE("combinators: the headline chain succeeds by bisection at run time and at compile time")
    {
        const auto got = headline(sq2);
        const auto p1  = r::newton {}.with_derivative(dsq2).on(0.0)(sq2);
        const auto p2  = r::secant {}.with_budget(5).on(0.0)(sq2);
        const auto p3  = r::bisection {}.on(nxx::bracket { 0.0, 2.0 })(sq2);
        if (!got || p1 || p2 || !p3) {
            FAIL_CHECK("precondition: newton and the secant fail, bisection succeeds");
            return;
        }
        CHECK(p1.error().code == nxx::errc::zero_derivative);
        CHECK(got->by == r::algos::bisection);
        CHECK(same_estimate(*got, *p3));
        CHECK(got->how == p3->how);
        CHECK(got->used == p1.error().used + p2.error().used + p3->used);    // the cost includes every attempt
        CHECK(got->used == nxx::counters { 57, 62 });                        // DESIGN §6.11
        CHECK(got->x == doctest::Approx(std::sqrt(2.0)).epsilon(1e-15));
        if (got->enclosure) {
            CHECK(got->enclosure->lo() <= std::sqrt(2.0));
            CHECK(std::sqrt(2.0) <= got->enclosure->hi());
        }
        else {
            FAIL_CHECK("a bisection result carries its enclosure");
        }

        static_assert(headline(sq2).has_value());
        static_assert(headline(sq2)->by == r::algos::bisection);
        static_assert(headline(sq2)->x > 1.4142135 && headline(sq2)->x < 1.4142136);
    }

    // ---- is_fatal and custom policies --------------------------------------------------------------------------------
    TEST_CASE("combinators: is_fatal is found by ADL and defaults to false")
    {
        static_assert(nxx::is_fatal(comb_user::eval_error::fatal));
        static_assert(!nxx::is_fatal(comb_user::eval_error::domain));
        static_assert(!nxx::is_fatal(comb_user::plain_error::oops));
        static_assert(!nxx::is_fatal(42));
        static_assert(!nxx::is_fatal(nxx::errc::stalled));
        CHECK(nxx::is_fatal(comb_user::eval_error::fatal));
    }

    TEST_CASE("combinators: a fatal user error stops first_of and a non-fatal one falls through")
    {
        std::uint32_t runs  = 0;
        const auto    bis   = r::bisection {}.on({ -1.0, 2.0 });
        const auto    rest  = r::brent {}.on({ 0.0, 2.0 });
        const auto    chain = nxx::first_of(bis, spy { rest, runs });

        const auto stopped = chain(g_fatal);
        CHECK(runs == 0u);    // the later alternative never ran
        if (stopped) { FAIL_CHECK("a fatal error must stop the chain"); }
        else {
            CHECK(stopped.error().code == nxx::errc::callback_failed);
            CHECK(stopped.error().by == r::algos::bisection);
            CHECK(stopped.error().cause == std::optional { comb_user::eval_error::fatal });
            CHECK(same_result(stopped, bis(g_fatal)));
        }

        const auto went_on = chain(g_domain);
        CHECK(runs == 1u);
        const auto e_bis  = bis(g_domain);
        const auto s_rest = rest(g_domain);
        if (!went_on || e_bis || !s_rest) { FAIL_CHECK("a non-fatal error falls through to brent, which succeeds"); }
        else {
            CHECK(went_on->by == r::algos::brent);
            CHECK(same_estimate(*went_on, *s_rest));
            CHECK(went_on->used == e_bis.error().used + s_rest->used);
        }

        // A three-alternative chain stops at the fatal error wherever it occurs.
        std::uint32_t runs3  = 0;
        const auto    chain3 = nxx::first_of(r::secant {}.with_budget(1).on(1.0), bis, spy { rest, runs3 });
        const auto    res3   = chain3(g_fatal);
        CHECK(runs3 == 0u);
        if (res3) { FAIL_CHECK("a fatal error must stop the chain"); }
        else {
            CHECK(res3.error().code == nxx::errc::callback_failed);
            CHECK(res3.error().cause == std::optional { comb_user::eval_error::fatal });
        }
    }

    TEST_CASE("combinators: first_of_with applies a custom policy at every alternative")
    {
        std::uint32_t runs   = 0;
        const auto    good   = r::brent {}.on(nxx::bracket { 1.0, 2.0 });
        const auto    strict = nxx::first_of_with(stop_on_any_failure {}, starved_secant, spy { good, runs });
        const auto    res    = strict(sq2);
        CHECK(runs == 0u);
        CHECK(same_result(res, starved_secant(sq2)));    // the first failure, unchanged

        // Continue on numerical errors, stop on input errors: the policy is also applied inside the nested chain.
        std::uint32_t runs3  = 0;
        const auto    chain3 = nxx::first_of_with(stop_on_input_error {}, zero_deriv, invalid, spy { good, runs3 });
        const auto    res3   = chain3(sq2);
        const auto    e_zd   = zero_deriv(sq2);
        CHECK(runs3 == 0u);
        if (res3 || e_zd) { FAIL_CHECK("the chain stops at the invalid input"); }
        else {
            CHECK(res3.error().code == nxx::errc::invalid_input);
            CHECK(res3.error().used == e_zd.error().used);
            CHECK(same_best(res3.error().best, e_zd.error().best));
        }

        // The same policy lets a numerical failure fall through to a success.
        std::uint32_t runs2  = 0;
        const auto    chain2 = nxx::first_of_with(stop_on_input_error {}, zero_deriv, spy { good, runs2 });
        const auto    res2   = chain2(sq2);
        CHECK(runs2 == 1u);
        CHECK(res2.has_value());
    }

    TEST_CASE("combinators: first_of with one alternative returns it unchanged and nests three or more")
    {
        const auto only = r::brent {}.on(nxx::bracket { 1.0, 2.0 });
        using only_t    = std::remove_const_t<decltype(only)>;
        static_assert(std::is_same_v<decltype(nxx::first_of(only)), only_t>);
        static_assert(std::is_same_v<decltype(nxx::first_of_with(stop_on_any_failure {}, only)), only_t>);
        CHECK(same_result(nxx::first_of(only)(sq2), only(sq2)));
        CHECK(same_result(nxx::first_of_with(stop_on_any_failure {}, only)(sq2), only(sq2)));

        // Four alternatives: the first three fail, the last succeeds.
        const auto four = nxx::first_of(zero_deriv, starved_secant, invalid, only);
        using zd_t      = std::remove_const_t<decltype(zero_deriv)>;
        using ss_t      = std::remove_const_t<decltype(starved_secant)>;
        using inv_t     = std::remove_const_t<decltype(invalid)>;
        using P         = nxx::continue_unless_fatal;
        static_assert(std::is_same_v<std::remove_const_t<decltype(four)>,
                                     nxx::first_of_t<P, zd_t, nxx::first_of_t<P, ss_t, nxx::first_of_t<P, inv_t, only_t>>>>);
        const auto got = four(sq2);
        const auto e1  = zero_deriv(sq2);
        const auto e2  = starved_secant(sq2);
        const auto s4  = only(sq2);
        if (!got || e1 || e2 || !s4) {
            FAIL_CHECK("precondition: the first three fail, brent succeeds");
            return;
        }
        CHECK(got->by == r::algos::brent);
        CHECK(same_estimate(*got, *s4));
        CHECK(got->used == e1.error().used + e2.error().used + s4->used);    // invalid costs nothing

        // Grouping does not matter.
        CHECK(same_result(nxx::first_of(nxx::first_of(zero_deriv, starved_secant), nxx::first_of(invalid, only))(sq2), got));
        CHECK(same_result(nxx::first_of(zero_deriv, nxx::first_of(starved_secant, nxx::first_of(invalid, only)))(sq2), got));

        // Four failures: the last code, the total cost, and the best estimate (here the narrow enclosure, which wins
        // every pairwise comparison: it is narrower than the wide one and has a smaller residual than the others).
        const auto all_fail = nxx::first_of(zero_deriv, starved_secant, wide, narrow)(sq2);
        const auto e_w      = wide(sq2);
        const auto e_n      = narrow(sq2);
        if (all_fail || e_w || e_n) {
            FAIL_CHECK("precondition: every alternative fails");
            return;
        }
        CHECK(all_fail.error().code == nxx::errc::budget_exhausted);
        CHECK(all_fail.error().by == r::algos::bisection);
        CHECK(all_fail.error().used == e1.error().used + e2.error().used + e_w.error().used + e_n.error().used);
        CHECK(same_best(all_fail.error().best, e_n.error().best));
    }

    TEST_CASE("combinators: a first_of nested as the first alternative keeps its own storage (cl layout regression)")
    {
        // On cl, [[msvc::no_unique_address]] on first_of_t's leading empty policy overlapped the inner chain with the
        // outer chain's second alternative, so first_of(first_of(a, b), c) ran a corrupted b (DESIGN §5.3).
        const auto a     = r::newton {}.with_derivative(dsq2).on(0.0);
        const auto b     = r::secant {}.with_budget(5).on(0.0);
        const auto c     = r::bisection {}.on(nxx::bracket { 0.0, 2.0 });
        const auto inner = nxx::first_of(a, b);
        const auto left  = nxx::first_of(inner, c);
        static_assert(sizeof(left) >= sizeof(inner) + sizeof(c));

        // x^2 + 1 has no root: every alternative runs, and the merged failure depends on each one's state.
        const auto no_root = [](double x) { return x * x + 1.0; };
        const auto flat    = nxx::first_of(a, b, c)(no_root);
        const auto nested  = left(no_root);
        if (flat || nested) {
            FAIL_CHECK("precondition: every alternative fails on x^2 + 1");
            return;
        }
        CHECK(nested.error().code == flat.error().code);
        CHECK(nested.error().used == flat.error().used);
        CHECK(same_best(nested.error().best, flat.error().best));
        CHECK(same_result(nested, flat));
        CHECK(same_result(left(sq2), nxx::first_of(a, b, c)(sq2)));
    }

    // ---- then --------------------------------------------------------------------------------------------------------
    TEST_CASE("combinators: then runs expand - bisection - newton and adds the cost of every stage")
    {
        const auto s1 = r::expand {}.on(nxx::bracket { 2.0, 2.5 });
        const auto s2 = r::bisection { nxx::width_tol { 1e-4 } };
        const auto s3 = r::newton {}.with_derivative(dsq2);

        const auto got = pipeline(sq2);
        const auto r1  = s1(sq2);
        if (!got || !r1) {
            FAIL_CHECK("precondition: expand finds a sign change and the pipeline succeeds");
            return;
        }
        const auto r2 = s2(sq2, *r1);    // a search result is a bracketing solver's input
        if (!r2) {
            FAIL_CHECK("precondition: bisection succeeds from the search result");
            return;
        }
        const auto r3 = s3(sq2, *r2);    // a root estimate seeds an open method
        if (!r3) {
            FAIL_CHECK("precondition: newton succeeds from the bisection result");
            return;
        }
        CHECK(r1->by == r::algos::expand);
        CHECK(r2->by == r::algos::bisection);
        CHECK(got->by == r::algos::newton);
        CHECK(same_estimate(*got, *r3));
        CHECK(got->how == r3->how);
        CHECK(got->used == r1->used + r2->used + r3->used);
        CHECK(got->x == doctest::Approx(std::sqrt(2.0)).epsilon(1e-15));

        // The sampled search result is not re-evaluated by bisection, the bisection result not by Newton: every call of
        // f (and of f') is accounted for.
        std::uint32_t calls   = 0;
        const auto    counted = nxx::then(s1, s2, r::newton {}.with_derivative(nxx::fn::counted(dsq2, calls)));
        const auto    got_c   = counted(nxx::fn::counted(sq2, calls));
        if (got_c) {
            CHECK(got_c->used == got->used);
            CHECK(calls == got_c->used.evaluations);
        }
        else {
            FAIL_CHECK("the counted pipeline must succeed as well");
        }
    }

    TEST_CASE("combinators: then propagates a stage failure unchanged and skips the later stages")
    {
        std::uint32_t runs   = 0;
        const auto    stage1 = r::bisection {}.on({ 3.0, 4.0 });    // no sign change for x2 - 2
        const auto    chain  = nxx::then(stage1, spy { r::secant {}, runs });
        const auto    got    = chain(sq2);
        const auto    alone  = stage1(sq2);
        static_assert(std::is_same_v<decltype(got), decltype(alone)>);    // one result type through the chain
        CHECK(runs == 0u);
        CHECK(same_result(got, alone));
        if (!got) {
            CHECK(got.error().code == nxx::errc::no_sign_change);
            CHECK(got.error().by == r::algos::bisection);
            CHECK(got.error().best.has_value());
        }
        else {
            FAIL_CHECK("the chain cannot succeed");
        }

        // A failing search: the searcher's failure type is the solvers' failure type.
        std::uint32_t runs_s = 0;
        const auto    search = r::expand {}.with_budget(3).on(nxx::bracket { 1.0, 2.0 });
        const auto    via_s  = nxx::then(search, spy { r::bisection {}, runs_s })([](double x) { return x * x + 1.0; });
        CHECK(runs_s == 0u);
        if (!via_s) {
            CHECK(via_s.error().code == nxx::errc::budget_exhausted);
            CHECK(via_s.error().by == r::algos::expand);
            CHECK(via_s.error().used.iterations == 3u);
        }
        else {
            FAIL_CHECK("x2 + 1 has no sign change");
        }

        // A failing middle stage carries the cost of the stages before it; the last stage never runs.
        std::uint32_t runs3 = 0;
        const auto    s1    = r::expand {}.on(nxx::bracket { 2.0, 2.5 });
        const auto    s2    = r::bisection { nxx::never {} }.with_budget(2);
        const auto    three = nxx::then(s1, s2, spy { r::newton {}.with_derivative(dsq2), runs3 });
        const auto    got3  = three(sq2);
        const auto    r1    = s1(sq2);
        CHECK(runs3 == 0u);
        if (got3 || !r1) {
            FAIL_CHECK("precondition: expand succeeds and the starved bisection fails");
            return;
        }
        const auto r2 = s2(sq2, *r1);
        if (r2) {
            FAIL_CHECK("precondition: the starved bisection fails");
            return;
        }
        CHECK(got3.error().code == nxx::errc::budget_exhausted);
        CHECK(got3.error().by == r::algos::bisection);
        CHECK(got3.error().used == r1->used + r2.error().used);
        CHECK(same_best(got3.error().best, r2.error().best));
    }

    TEST_CASE("combinators: then seeds the secant with the bisection's x and f(x) without re-evaluating it")
    {
        std::vector<double> xs;    // every point at which f is called, in order
        const auto          logged = [&xs](double x) {
            xs.push_back(x);
            return x * x - 2.0;
        };
        const auto stage1 = r::bisection { nxx::width_tol { 1e-3 } }.with_budget(100).on({ 1.0, 2.0 });
        const auto chain  = nxx::then(stage1, r::secant {});
        const auto got    = chain(logged);
        const auto rb     = stage1(sq2);
        if (!got || !rb) {
            FAIL_CHECK("precondition: the coarse bisection and the chain succeed");
            return;
        }
        const auto rs = r::secant {}(sq2, *rb);      // seeded with x and f(x)
        const auto rp = r::secant {}(sq2, rb->x);    // from the bare guess: f(x) is evaluated once more
        if (!rs || !rp) {
            FAIL_CHECK("precondition: the secant polish succeeds");
            return;
        }
        CHECK(got->by == r::algos::secant);
        CHECK(got->used == rb->used + rs->used);
        CHECK(same_estimate(*got, *rs));
        CHECK(xs.size() == got->used.evaluations);    // every call of f is counted, none is hidden

        const std::size_t nb = rb->used.evaluations;
        if (xs.size() > nb) {
            CHECK_FALSE(same_bits(xs[nb], rb->x));    // the secant's first call is its second point, not the seed
        }
        else {
            FAIL_CHECK("the secant stage evaluated nothing");
        }

        CHECK(rp->used.evaluations == rs->used.evaluations + 1u);
        CHECK(rp->used.iterations == rs->used.iterations);
        CHECK(same_bits(rp->x, rs->x));    // same trajectory: the seed's f(x) is the value f returns there
    }

    TEST_CASE("combinators: stage 2 of then and warm_fallback reports an input code as non_finite_value (regression)")
    {
        // Stage 2 starts from stage 1's value, not the caller's input (DESIGN §6.10, §12 item 23). A projection that
        // sends that value off the reals is non_finite_input for Newton alone, but not for the chain.
        const auto off = [](double x) { return x > 1.2 ? std::numeric_limits<double>::quiet_NaN() : x; };
        const auto nt  = r::newton {}.with_derivative(dsq2).with_projection(off);
        const auto s1  = r::secant {}.on(3.0);
        const auto r1  = s1(sq2);
        const auto own = nt(sq2, 3.0);    // the caller's own guess: a genuine input error
        if (!r1 || own) {
            FAIL_CHECK("precondition: the secant succeeds and the projected Newton rejects 3");
            return;
        }
        CHECK(own.error().code == nxx::errc::non_finite_input);

        const auto staged = nxx::then(s1, nt);
        const auto got    = staged(sq2);
        if (got) {
            FAIL_CHECK("stage 2 rejects its start");
            return;
        }
        CHECK(got.error().code == nxx::errc::non_finite_value);
        CHECK(got.error().by == r::algos::newton);
        CHECK(got.error().used == r1->used);
        CHECK(got.error().used == nxx::counters { 8, 10 });
        CHECK(same_best(got.error().best, static_cast<const est_t&>(*r1)));

        const auto s1w = r::secant {}.with_budget(3).on(3.0);
        const auto e1w = s1w(sq2);
        const auto wf  = nxx::warm_fallback(s1w, nt)(sq2);
        if (wf || e1w) {
            FAIL_CHECK("precondition: the starved secant and the fallback fail");
            return;
        }
        CHECK(wf.error().code == nxx::errc::non_finite_value);
        CHECK(wf.error().by == r::algos::newton);
        CHECK(wf.error().used == e1w.error().used);
        CHECK(wf.error().used == nxx::counters { 3, 5 });
        CHECK(same_best(wf.error().best, e1w.error().best));

        // A policy that stops on input errors gives the same outcome in either order: the bisection's root.
        const auto bis    = r::bisection {}.on(nxx::bracket { 1.0, 2.0 });
        const auto rb     = bis(sq2);
        const auto first  = nxx::first_of_with(stop_on_input_error {}, staged, bis)(sq2);
        const auto second = nxx::first_of_with(stop_on_input_error {}, bis, staged)(sq2);
        const auto warm1  = nxx::first_of_with(stop_on_input_error {}, nxx::warm_fallback(s1w, nt), bis)(sq2);
        if (!rb || !first || !second || !warm1) {
            FAIL_CHECK("every first_of_with chain succeeds by bisection");
            return;
        }
        CHECK(first->by == r::algos::bisection);
        CHECK(same_estimate(*first, *rb));
        CHECK(same_estimate(*second, *rb));
        CHECK(same_estimate(*warm1, *rb));
        CHECK(first->used == got.error().used + rb->used);
        CHECK(second->used == rb->used);
        CHECK(warm1->used == wf.error().used + rb->used);
    }

    TEST_CASE("combinators: a then whose stage 2 fails carries stage 1's estimate (regression)")
    {
        // f fails (NaN) for x <= 1. Stage 1 succeeds; stage 2 fails at its own start, before any running best exists.
        // The failure keeps stage 1's estimate (DESIGN §3.4, §6.10), and every call of f is counted.
        const auto    f_nan = [](double x) { return x <= 1.0 ? std::numeric_limits<double>::quiet_NaN() : x * x - 2.0; };
        const auto    s1    = r::secant {}.on(3.0);
        const auto    r1    = s1(f_nan);
        std::uint32_t calls = 0;
        if (!r1) {
            FAIL_CHECK("precondition: the secant from 3 succeeds");
            return;
        }
        const est_t s1_est = *r1;
        CHECK(r1->used == nxx::counters { 8, 10 });

        // Clamped to [0, 1], stage 2's start is 1, where f is NaN: its first evaluation fails.
        const auto nt_clamp = nxx::then(s1, r::newton {}.with_derivative(dsq2).with_projection(r::clamp_to { 0.0, 1.0 }));
        const auto got      = nt_clamp(nxx::fn::counted(f_nan, calls));
        if (got) {
            FAIL_CHECK("stage 2 fails at its clamped start");
            return;
        }
        CHECK(got.error().code == nxx::errc::non_finite_value);
        CHECK(got.error().used == nxx::counters { 8, 11 });
        CHECK(calls == got.error().used.evaluations);
        CHECK(same_best(got.error().best, s1_est));
        CHECK(same_best(nxx::best(got), s1_est));

        calls               = 0;
        const auto sc_clamp = nxx::then(s1, r::secant {}.with_projection(r::clamp_to { 0.0, 1.0 }))(nxx::fn::counted(f_nan, calls));
        if (sc_clamp) {
            FAIL_CHECK("stage 2 fails at its clamped start");
            return;
        }
        CHECK(sc_clamp.error().code == nxx::errc::non_finite_value);
        CHECK(calls == sc_clamp.error().used.evaluations);
        CHECK(same_best(sc_clamp.error().best, s1_est));

        // A search stage 1: its sign_bracket gives the root_estimate that the failure carries.
        const auto s_search = r::expand {}.on(nxx::bracket { 2.0, 2.5 });
        const auto reject   = [](const auto&, const auto&) {
            return nxx::result<est_t> { std::unexpect, nxx::failure<est_t> { nxx::errc::invalid_input, r::algos::bisection, {}, {}, {} } };
        };
        const auto rs = s_search(sq2);
        const auto gs = nxx::then(s_search, reject)(sq2);
        if (!rs || gs) {
            FAIL_CHECK("precondition: expand succeeds and the rejecting stage fails");
            return;
        }
        CHECK(gs.error().code == nxx::errc::non_finite_value);
        CHECK(gs.error().used == rs->used);
        CHECK(same_best(gs.error().best, rs->best()));
    }

    // ---- warm_fallback -----------------------------------------------------------------------------------------------
    TEST_CASE("combinators: warm_fallback restarts newton from a starved bisection's best estimate")
    {
        const auto starved = r::bisection { nxx::never {} }.with_budget(3).on(nxx::bracket { 0.0, 2.0 });
        const auto nt      = r::newton {}.with_derivative(dsq2);
        const auto got     = warm(sq2);
        const auto e1      = starved(sq2);
        if (!got || e1 || !e1.error().best) {
            FAIL_CHECK("precondition: the bisection is starved, with a best estimate, and the fallback succeeds");
            return;
        }
        CHECK(e1.error().code == nxx::errc::budget_exhausted);
        CHECK(e1.error().used == nxx::counters { 3, 5 });
        CHECK(e1.error().best->enclosure.has_value());
        const auto r2 = nt(sq2, *e1.error().best);
        if (!r2) {
            FAIL_CHECK("precondition: newton succeeds from the best estimate");
            return;
        }
        CHECK(got->by == r::algos::newton);
        CHECK(same_estimate(*got, *r2));
        CHECK(got->used == e1.error().used + r2->used);

        // Newton starts from the best estimate's x and f(x): every call of f and f' is counted, none repeated.
        std::uint32_t calls = 0;
        const auto    wf_c  = nxx::warm_fallback(starved, r::newton {}.with_derivative(nxx::fn::counted(dsq2, calls)));
        const auto    got_c = wf_c(nxx::fn::counted(sq2, calls));
        if (got_c) {
            CHECK(calls == got_c->used.evaluations);
            CHECK(got_c->used == got->used);
        }
        else {
            FAIL_CHECK("the counted fallback must succeed as well");
        }

        // Stage 1 succeeds: stage 2 never runs.
        std::uint32_t runs = 0;
        const auto    good = r::bisection {}.on(nxx::bracket { 0.0, 2.0 });
        const auto    ok   = nxx::warm_fallback(good, spy { nt, runs })(sq2);
        CHECK(runs == 0u);
        CHECK(same_result(ok, good(sq2)));
    }

    // A pole is no start for an open method (DESIGN §3.4, §6.10): from the pole failure's x, Newton's step and the
    // secant's are tiny, and their step criteria accepted the pole. Before, warm_fallback(brent, newton) on tan
    // succeeded at |f| about 5.8e14 and (bisection, secant) at |f| about 652, both with stop_reason::criterion.
    TEST_CASE("combinators: warm_fallback returns a stage-1 sign_change_not_root failure as is, without a restart (regression)")
    {
        const auto tangent = [](double x) { return std::tan(x); };
        const auto dtan    = [](double x) { return 1.0 / (std::cos(x) * std::cos(x)); };

        const auto check_pole = [&](const auto& s1, const auto& s2) {
            std::uint32_t runs = 0;
            const auto    e1   = s1(tangent);
            const auto    got  = nxx::warm_fallback(s1, spy { s2, runs })(tangent);
            CHECK(runs == 0u);
            CHECK(same_result(got, e1));
            if (got || !got.error().best) {
                FAIL_CHECK("the chain fails with the pole estimate");
                return;
            }
            CHECK(got.error().code == nxx::errc::sign_change_not_root);
            CHECK(got.error().used == e1.error().used);
            CHECK_FALSE(got.error().best->enclosure.has_value());
            CHECK(std::abs(got.error().best->x - std::numbers::pi / 2.0) < 1e-12);
        };
        check_pole(r::brent {}.on(nxx::bracket { 1.0, 2.0 }), r::newton {}.with_derivative(dtan));
        check_pole(r::bisection {}.on(nxx::bracket { 1.0, 2.0 }), r::secant {});
        check_pole(nxx::then(r::expand {}.on(nxx::bracket { 1.0, 2.0 }), r::brent {}), r::newton {}.with_derivative(dtan));
    }

    TEST_CASE("combinators: warm_fallback returns a stage-1 failure without a best estimate as is")
    {
        std::uint32_t runs = 0;
        const auto    wf   = nxx::warm_fallback(invalid, spy { r::newton {}.with_derivative(dsq2), runs });
        const auto    got  = wf(sq2);
        CHECK(runs == 0u);
        CHECK(same_result(got, invalid(sq2)));
        if (!got) {
            CHECK(got.error().code == nxx::errc::invalid_input);
            CHECK(got.error().used == nxx::counters {});
            CHECK_FALSE(got.error().best.has_value());
        }
        else {
            FAIL_CHECK("equal endpoints are invalid");
        }
    }

    TEST_CASE("combinators: warm_fallback merges two failures")
    {
        // x2 + 1: bisection finds no sign change (best: x = 0, |f| = 1); Newton from there hits f'(0) = 0.
        const auto sq_plus1  = [](double x) { return x * x + 1.0; };
        const auto dsq_plus1 = [](double x) { return 2.0 * x; };
        const auto s1        = r::bisection {}.on(nxx::bracket { 0.0, 2.0 });
        const auto s2        = r::newton {}.with_derivative(dsq_plus1);
        const auto got       = nxx::warm_fallback(s1, s2)(sq_plus1);
        const auto e1        = s1(sq_plus1);
        if (got || e1 || !e1.error().best) {
            FAIL_CHECK("precondition: no sign change with a best estimate");
            return;
        }
        const auto e2 = s2(sq_plus1, *e1.error().best);
        if (e2) {
            FAIL_CHECK("precondition: newton fails from x = 0");
            return;
        }
        CHECK(e1.error().code == nxx::errc::no_sign_change);
        CHECK(e2.error().code == nxx::errc::zero_derivative);
        CHECK(e2.error().used == nxx::counters { 1, 1 });    // seeded: only f'(0) is evaluated
        CHECK(got.error().code == nxx::errc::zero_derivative);
        CHECK(got.error().by == r::algos::newton);
        CHECK(got.error().used == e1.error().used + e2.error().used);
        if (got.error().best) {
            CHECK(got.error().best->x == 0.0);
            CHECK(got.error().best->fx == 1.0);
        }
        else {
            FAIL_CHECK("the merged failure keeps a best estimate");
        }
    }
}
