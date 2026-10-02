// Run-time solver chains (DESIGN §6.10, §9.3): any_solver erases the type of a curried solver; first_of over a run-time
// range of them has the static chain's laziness, merge and cost accounting, and returns exactly what the equivalent
// static chain returns, bit for bit, on success and on total failure.

#include <numerixx/core/any_solver.hpp>
#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <array>
#include <bit>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <expected>
#include <functional>
#include <optional>
#include <random>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

// Results are compared bit for bit: no fused multiply-adds in this file's functions (Clang contracts them otherwise).
#if defined(__clang__)
#    pragma clang fp contract(off)
#endif

namespace
{
    namespace r    = nxx::roots;
    using est_t    = r::root_estimate<double>;
    using fn_t     = std::function<double(double)>;
    using solver_t = nxx::any_solver<fn_t, est_t>;    // plain callbacks: UE = none
    using result_t = solver_t::result_type;

    constexpr auto sq2      = [](double x) { return x * x - 2.0; };    // root sqrt(2)
    constexpr auto sq_plus1 = [](double x) { return x * x + 1.0; };    // no real root
    constexpr auto twice    = [](double x) { return 2.0 * x; };        // the derivative of both

    namespace rt_user
    {
        enum class eval_error { domain, fatal };

        constexpr bool is_fatal(eval_error err) noexcept { return err == eval_error::fatal; }
    }    // namespace rt_user

    using gfn_t     = std::function<std::expected<double, rt_user::eval_error>(double)>;
    using gsolver_t = nxx::any_solver<gfn_t, est_t, rt_user::eval_error>;

    constexpr auto g_domain = [](double x) -> std::expected<double, rt_user::eval_error> {
        if (x < 0.0) return std::unexpected(rt_user::eval_error::domain);
        return x * x - 2.0;
    };
    constexpr auto g_fatal = [](double x) -> std::expected<double, rt_user::eval_error> {
        if (x < 0.0) return std::unexpected(rt_user::eval_error::fatal);
        return x * x - 2.0;
    };

    // The alternatives of the chains below.
    constexpr auto newton0  = r::newton {}.with_derivative(twice).on(0.0);      // f'(0) = 0: zero_derivative
    constexpr auto secant5  = r::secant {}.with_budget(5).on(0.0);              // fails within 5 iterations
    constexpr auto bisect02 = r::bisection {}.on(nxx::bracket { 0.0, 2.0 });    // succeeds on x^2 - 2
    constexpr auto invalid  = r::bisection {}.on(std::pair { 1.0, 1.0 });       // invalid_input, no cost, no best

    auto     static_chain() { return nxx::first_of(newton0, secant5, bisect02); }
    solver_t runtime_chain() { return nxx::first_of(std::vector<solver_t> { newton0, secant5, bisect02 }); }

    // ---- Type properties -------------------------------------------------------------------------------------------
    static_assert(std::is_copy_constructible_v<solver_t>);
    static_assert(std::is_copy_assignable_v<solver_t>);
    static_assert(std::copyable<solver_t>);
    static_assert(!std::is_default_constructible_v<solver_t>);    // never empty
    static_assert(std::is_invocable_r_v<result_t, const solver_t&, const fn_t&>);
    static_assert(!std::is_nothrow_invocable_v<const solver_t&, const fn_t&>);    // f may throw
    static_assert(std::is_same_v<result_t, nxx::result<est_t>>);

    using brent_on_t = decltype(r::brent {}.on(nxx::bracket { 0.0, 2.0 }));
    using searcher_t = decltype(r::expand {}.on(nxx::bracket { 1.0, 2.0 }));
    using static3_t  = decltype(static_chain());
    static_assert(std::is_constructible_v<solver_t, brent_on_t>);
    static_assert(std::is_convertible_v<brent_on_t, solver_t>);       // implicit: every matching solver value is one
    static_assert(std::is_constructible_v<solver_t, static3_t>);      // a static chain converts too
    static_assert(!std::is_constructible_v<solver_t, searcher_t>);    // a search result is another success type
    static_assert(!std::is_convertible_v<searcher_t, solver_t>);
    static_assert(!std::is_constructible_v<nxx::any_solver<fn_t, r::root_estimate<float>>, brent_on_t>);       // another estimate
    static_assert(!std::is_constructible_v<nxx::any_solver<fn_t, est_t, rt_user::eval_error>, brent_on_t>);    // another cause
    static_assert(!std::is_constructible_v<solver_t, int>);    // the deleted constructor, not a hard error
    static_assert(!std::is_constructible_v<solver_t, gsolver_t>);
    static_assert(std::is_constructible_v<gsolver_t, brent_on_t>);    // the same solver value, called with a fallible f

    // first_of over a range is itself an any_solver, whatever the range.
    static_assert(std::is_same_v<decltype(nxx::first_of(std::declval<std::vector<solver_t>>())), solver_t>);
    static_assert(std::is_same_v<decltype(nxx::first_of(std::declval<std::span<const solver_t>>())), solver_t>);
    static_assert(std::is_same_v<decltype(nxx::first_of(std::declval<std::array<solver_t, 2>>())), solver_t>);
    static_assert(
        std::is_same_v<decltype(nxx::first_of_with(nxx::continue_unless_fatal {}, std::declval<std::vector<solver_t>>())), solver_t>);
    static_assert(nxx::is_any_solver_v<solver_t> && !nxx::is_any_solver_v<brent_on_t>);

    // With this header included, the one-argument first_of and first_of_with still return their solver: a non-range
    // argument makes the range overload's constraint false instead of a hard error (a single any_solver included).
    using bisect02_t = std::remove_cvref_t<decltype(bisect02)>;
    static_assert(std::is_same_v<decltype(nxx::first_of(bisect02)), bisect02_t>);
    static_assert(std::is_same_v<decltype(nxx::first_of(std::declval<solver_t>())), solver_t>);
    static_assert(std::is_same_v<decltype(nxx::first_of_with(nxx::continue_unless_fatal {}, bisect02)), bisect02_t>);

    // ---- Helpers ---------------------------------------------------------------------------------------------------
    bool same_bits(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }

    bool same_estimate(const est_t& a, const est_t& b)
    { return same_bits(a.x, b.x) && same_bits(a.fx, b.fx) && same_bits(a.uncertainty, b.uncertainty) && a.enclosure == b.enclosure; }

    bool same_best(const std::optional<est_t>& a, const std::optional<est_t>& b)
    {
        if (a.has_value() != b.has_value()) return false;
        return !a || same_estimate(*a, *b);
    }

    template<class R>
    bool same_result(const R& a, const R& b)
    {
        if (a.has_value() != b.has_value()) return false;
        if (a) return same_estimate(*a, *b) && a->used == b->used && a->by == b->by && a->how == b->how;
        const auto& fa = a.error();
        const auto& fb = b.error();
        return fa.code == fb.code && fa.where == fb.where && fa.used == fb.used && same_best(fa.best, fb.best) && fa.cause == fb.cause;
    }

    // A curried solver that counts how often it runs.
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
    TEST_CASE("any_solver: a run-time chain equals the static chain bit for bit on success")
    {
        const fn_t     f   = sq2;
        const solver_t rt  = runtime_chain();
        const auto     got = rt(f);
        const auto     st  = static_chain()(f);      // the static chain given the same fn_t
        const auto     sl  = static_chain()(sq2);    // ... and given the plain lambda
        CHECK(same_result(got, st));
        CHECK(same_result(got, sl));
        if (!got || !st) {
            FAIL_CHECK("the chain succeeds by bisection");
            return;
        }
        CHECK(same_bits(got->x, st->x));
        CHECK(same_bits(got->fx, st->fx));
        CHECK(got->used.iterations == st->used.iterations);
        CHECK(got->used.evaluations == st->used.evaluations);
        CHECK(got->by == st->by);
        CHECK(got->by == r::algos::bisection);
        CHECK(got->x == doctest::Approx(std::sqrt(2.0)).epsilon(1e-15));

        // The cost includes every attempt, and every call of f and f' is counted.
        const auto e1 = newton0(f);
        const auto e2 = secant5(f);
        const auto s3 = bisect02(f);
        if (e1 || e2 || !s3) {
            FAIL_CHECK("precondition: newton and the secant fail, bisection succeeds");
            return;
        }
        CHECK(got->used == e1.error().used + e2.error().used + s3->used);

        std::uint32_t  calls   = 0;
        const solver_t counted = nxx::first_of(
            std::vector<solver_t> { r::newton {}.with_derivative(nxx::fn::counted(twice, calls)).on(0.0), secant5, bisect02 });
        const fn_t counted_f = nxx::fn::counted(sq2, calls);
        const auto got_c     = counted(counted_f);
        CHECK(same_result(got_c, got));
        if (got_c) { CHECK(calls == got_c->used.evaluations); }
    }

    TEST_CASE("any_solver: a run-time chain equals the static chain bit for bit on total failure")
    {
        const fn_t g   = sq_plus1;
        const auto got = runtime_chain()(g);
        const auto st  = static_chain()(g);
        const auto sl  = static_chain()(sq_plus1);
        CHECK(same_result(got, st));
        CHECK(same_result(got, sl));
        if (got || st) {
            FAIL_CHECK("x2 + 1 has no real root");
            return;
        }
        CHECK(got.error().code == st.error().code);
        CHECK(got.error().code == nxx::errc::no_sign_change);    // the last code: bisection's
        CHECK(got.error().where == r::algos::bisection);
        CHECK(got.error().used == st.error().used);
        CHECK(same_best(got.error().best, st.error().best));
        CHECK(got.error().best.has_value());

        const auto e1 = newton0(g);
        const auto e2 = secant5(g);
        const auto e3 = bisect02(g);
        if (e1 || e2 || e3) {
            FAIL_CHECK("every alternative fails on x2 + 1");
            return;
        }
        CHECK(got.error().used == e1.error().used + e2.error().used + e3.error().used);
    }

    TEST_CASE("any_solver: a run-time chain equals the static chain when the best estimates mix enclosures")
    {
        // Three failures: a wide enclosure with a tiny residual, no enclosure with a small residual, a narrow enclosure
        // with a larger residual. The static and the run-time chain must agree on the best estimate over all attempts.
        constexpr auto wide    = r::bisection { nxx::never {} }.with_budget(1).on(nxx::bracket { 1.4142135, 3.0 });
        constexpr auto secant3 = r::secant {}.with_budget(3).on(1.0);
        constexpr auto narrow  = r::bisection { nxx::never {} }.with_budget(6).on(nxx::bracket { 0.0, 2.0 });

        const fn_t f   = sq2;
        const auto st  = nxx::first_of(wide, secant3, narrow)(f);
        const auto got = nxx::first_of(std::vector<solver_t> { wide, secant3, narrow })(f);
        if (got || st) {
            FAIL_CHECK("every alternative is starved");
            return;
        }
        CHECK(got.error().code == st.error().code);
        CHECK(got.error().used == st.error().used);
        CHECK(same_best(got.error().best, st.error().best));
        CHECK(same_result(got, st));
    }

    TEST_CASE("any_solver: an empty run-time chain fails with invalid_input at zero cost and without a best estimate")
    {
        const solver_t empty = nxx::first_of(std::vector<solver_t> {});
        const fn_t     f     = sq2;
        const auto     got   = empty(f);
        if (got) {
            FAIL_CHECK("an empty chain cannot succeed");
            return;
        }
        CHECK(got.error().code == nxx::errc::invalid_input);
        CHECK(got.error().where == nxx::algo::none);
        CHECK(got.error().used == nxx::counters {});
        CHECK_FALSE(got.error().best.has_value());

        const solver_t empty_with = nxx::first_of_with(stop_on_any_failure {}, std::vector<solver_t> {});
        CHECK(same_result(empty_with(f), got));
    }

    TEST_CASE("any_solver: an alternative after a success never runs (property)")
    {
        std::mt19937                           gen(20260928u);
        std::uniform_real_distribution<double> coef(0.5, 3.5);    // roots sqrt(c) in (0.7, 1.9), inside [0, 2]
        std::uint32_t                          runs1 = 0;
        std::uint32_t                          runs2 = 0;
        const auto                             s1    = r::brent {}.on(nxx::bracket { 0.0, 2.0 });
        const solver_t                         chain = nxx::first_of(std::vector<solver_t> { spy { s1, runs1 }, spy { bisect02, runs2 } });
        for (int i = 0; i < 100; ++i) {
            const double c     = coef(gen);
            const fn_t   fc    = [c](double x) { return x * x - c; };
            const auto   alone = s1(fc);
            CHECK(alone.has_value());
            CHECK(same_result(chain(fc), alone));
        }
        CHECK(runs1 == 100u);
        CHECK(runs2 == 0u);
    }

    TEST_CASE("any_solver: is_fatal stops a run-time chain and a non-fatal error falls through")
    {
        std::uint32_t   runs  = 0;
        const auto      bis   = r::bisection {}.on({ -1.0, 2.0 });    // g fails at -1
        const auto      rest  = r::brent {}.on({ 0.0, 2.0 });
        const gsolver_t chain = nxx::first_of(std::vector<gsolver_t> { bis, spy { rest, runs } });
        const gfn_t     gf    = g_fatal;
        const gfn_t     gd    = g_domain;

        const auto stopped = chain(gf);
        CHECK(runs == 0u);
        CHECK(same_result(stopped, nxx::first_of(bis, rest)(gf)));
        if (stopped) { FAIL_CHECK("a fatal error must stop the chain"); }
        else {
            CHECK(stopped.error().code == nxx::errc::callback_failed);
            CHECK(stopped.error().cause == std::optional { rt_user::eval_error::fatal });
        }

        const auto went_on = chain(gd);
        CHECK(runs == 1u);
        CHECK(same_result(went_on, nxx::first_of(bis, rest)(gd)));
        const auto e_bis  = bis(gd);
        const auto s_rest = rest(gd);
        if (!went_on || e_bis || !s_rest) { FAIL_CHECK("a non-fatal error falls through to brent, which succeeds"); }
        else {
            CHECK(went_on->by == r::algos::brent);
            CHECK(went_on->used == e_bis.error().used + s_rest->used);
        }
    }

    TEST_CASE("any_solver: first_of_with applies its policy to a run-time range")
    {
        const fn_t f = sq2;

        std::uint32_t  runs   = 0;
        const solver_t strict = nxx::first_of_with(stop_on_any_failure {}, std::vector<solver_t> { newton0, spy { bisect02, runs } });
        const auto     res    = strict(f);
        CHECK(runs == 0u);
        CHECK(same_result(res, newton0(f)));
        CHECK(same_result(res, nxx::first_of_with(stop_on_any_failure {}, newton0, bisect02)(f)));

        // Continue on numerical errors, stop on input errors.
        std::uint32_t  runs3 = 0;
        const solver_t lenient =
            nxx::first_of_with(stop_on_input_error {}, std::vector<solver_t> { newton0, invalid, spy { bisect02, runs3 } });
        const auto res3 = lenient(f);
        CHECK(runs3 == 0u);
        CHECK(same_result(res3, nxx::first_of_with(stop_on_input_error {}, newton0, invalid, bisect02)(f)));
        const auto e_zd = newton0(f);
        if (res3 || e_zd) { FAIL_CHECK("the chain stops at the invalid input"); }
        else {
            CHECK(res3.error().code == nxx::errc::invalid_input);
            CHECK(res3.error().used == e_zd.error().used);
            CHECK(same_best(res3.error().best, e_zd.error().best));
        }

        std::uint32_t  runs2 = 0;
        const solver_t goes  = nxx::first_of_with(stop_on_input_error {}, std::vector<solver_t> { newton0, spy { bisect02, runs2 } });
        const auto     res2  = goes(f);
        CHECK(runs2 == 1u);
        CHECK(res2.has_value());
    }

    TEST_CASE("any_solver: a static first_of over a run-time chain and a run-time chain inside another")
    {
        const solver_t tail   = nxx::first_of(std::vector<solver_t> { secant5, bisect02 });
        const auto     mixed  = nxx::first_of(newton0, tail);                              // static over run-time
        const solver_t nested = nxx::first_of(std::vector<solver_t> { newton0, tail });    // run-time over run-time
        const solver_t erased = static_chain();                                            // a static chain, erased

        for (const fn_t& fx : { fn_t { sq2 }, fn_t { sq_plus1 } }) {    // success and total failure
            const auto want = static_chain()(fx);
            CHECK(same_result(mixed(fx), want));
            CHECK(same_result(nested(fx), want));
            CHECK(same_result(erased(fx), want));
        }
    }

    TEST_CASE("any_solver: first_of over std::span, std::array and an lvalue vector owns its alternatives")
    {
        const fn_t f    = sq2;
        const auto want = static_chain()(f);

        const std::array<solver_t, 3> arr { solver_t { newton0 }, solver_t { secant5 }, solver_t { bisect02 } };
        const solver_t                from_array = nxx::first_of(arr);
        CHECK(same_result(from_array(f), want));

        const std::vector<solver_t> vec { newton0, secant5, bisect02 };
        const solver_t              from_span = nxx::first_of(std::span<const solver_t> { vec });
        const solver_t              from_vec  = nxx::first_of(vec);
        CHECK(same_result(from_span(f), want));
        CHECK(same_result(from_vec(f), want));

        // The chain copies the range: it outlives the storage it was built from.
        std::optional<solver_t> survivor;
        {
            std::vector<solver_t> scratch { newton0, secant5, bisect02 };
            survivor.emplace(nxx::first_of(std::span<solver_t> { scratch }));
        }
        if (survivor) { CHECK(same_result((*survivor)(f), want)); }
        else {
            FAIL_CHECK("emplace failed");
        }

        // Also with a policy.
        const solver_t with_policy = nxx::first_of_with(nxx::continue_unless_fatal {}, std::span<const solver_t> { vec });
        CHECK(same_result(with_policy(f), want));
    }

    TEST_CASE("any_solver: the one-argument first_of still returns its solver with this header included")
    {
        const fn_t     f   = sq2;
        const solver_t one = bisect02;
        CHECK(same_result(nxx::first_of(bisect02)(f), bisect02(f)));
        CHECK(same_result(nxx::first_of(one)(f), bisect02(f)));
        CHECK(same_result(nxx::first_of_with(stop_on_any_failure {}, one)(f), bisect02(f)));
    }

    TEST_CASE("any_solver: copies, assignment and moves keep a working chain")
    {
        const fn_t     f    = sq2;
        const solver_t b    = runtime_chain();
        const auto     want = b(f);

        solver_t a = nxx::first_of(std::vector<solver_t> { invalid });
        CHECK_FALSE(a(f).has_value());
        a = b;    // replaces the whole value
        CHECK(same_result(a(f), want));

        solver_t moved = std::move(a);    // a move is a copy: a is never empty
        CHECK(same_result(moved(f), want));
        CHECK(same_result(a(f), want));    // NOLINT(bugprone-use-after-move): the invariant under test

        std::vector<solver_t> v { newton0, bisect02 };
        std::swap(v[0], v[1]);
        CHECK(same_result(v[0](f), bisect02(f)));
        CHECK(same_result(v[1](f), newton0(f)));
    }

    TEST_CASE("any_solver: a fallible run-time chain keeps the user's error when everything fails")
    {
        const auto      sec1  = r::secant {}.with_budget(1).on(1.0);    // budget_exhausted, a best estimate, no cause
        const auto      bis   = r::bisection {}.on({ -1.0, 2.0 });      // callback_failed with the user's error
        const gsolver_t chain = nxx::first_of(std::vector<gsolver_t> { sec1, bis });
        const gfn_t     gd    = g_domain;
        const auto      got   = chain(gd);
        CHECK(same_result(got, nxx::first_of(sec1, bis)(gd)));
        CHECK(same_result(got, nxx::first_of(sec1, bis)(g_domain)));
        const auto e1 = sec1(gd);
        const auto e2 = bis(gd);
        if (got || e1 || e2) {
            FAIL_CHECK("every alternative fails");
            return;
        }
        CHECK(got.error().code == nxx::errc::callback_failed);
        CHECK(got.error().cause == std::optional { rt_user::eval_error::domain });
        CHECK(got.error().where == r::algos::bisection);
        CHECK(same_best(got.error().best, e1.error().best));
        CHECK(got.error().used == e1.error().used + e2.error().used);

        // Every alternative fails with the user's error.
        const gsolver_t all_bad = nxx::first_of(std::vector<gsolver_t> { bis, r::brent {}.on({ -2.0, 1.0 }) });
        const auto      res     = all_bad(gd);
        if (res) { FAIL_CHECK("both alternatives evaluate g below 0"); }
        else {
            CHECK(res.error().code == nxx::errc::callback_failed);
            CHECK(res.error().cause == std::optional { rt_user::eval_error::domain });
            CHECK(res.error().where == r::algos::brent);
            CHECK_FALSE(res.error().best.has_value());
        }
    }
}
