// Results (DESIGN §6.3): nxx::best gives the solution's estimate or the failure's best estimate as one std::optional,
// for every result whose success and failure carry the same estimate, through first_of and any_solver too. A search
// result (a sign_bracket on success, a root_estimate on failure) reaches best's deleted sibling, and a one-shot
// derivative's std::expected<T, fault<>> is not a result. best_x needs a solution with an x. detail::is_result_v
// recognises the results the combinators' classifiers will check (§6.10).

#include <numerixx/core/any_solver.hpp>
#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <expected>
#include <functional>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{
    namespace r = nxx::roots;

    enum class user_error { domain };

    using est_t           = r::root_estimate<double>;
    using root_result_t   = nxx::result<est_t>;
    using cause_result_t  = nxx::result<est_t, user_error>;
    using search_result_t = std::expected<nxx::solution<r::sign_bracket<double>>, nxx::failure<est_t>>;
    using deriv_result_t  = std::expected<double, nxx::fault<>>;    // one-shot differentiation (§6.12)
    using fn_t            = std::function<double(double)>;
    using solver_t        = nxx::any_solver<fn_t, est_t>;

    template<class R>
    concept best_ok = requires(const R& res) { nxx::best(res); };

    template<class R>
    concept best_x_ok = requires(const R& res) { nxx::best_x(res); };

    // ---- nxx::best: accepted exactly for the results whose success and failure carry the same estimate -------------
    static_assert(best_ok<root_result_t>);
    static_assert(best_ok<cause_result_t>);
    static_assert(best_ok<solver_t::result_type>);
    static_assert(!best_ok<search_result_t>);    // the deleted sibling, with its reason (cf.best_search_result)
    static_assert(!best_ok<deriv_result_t>);     // not a result: constrained away, with no reason
    static_assert(!best_ok<std::expected<double, int>>);
    static_assert(std::is_same_v<decltype(nxx::best(std::declval<const root_result_t&>())), std::optional<est_t>>);
    static_assert(std::is_same_v<decltype(nxx::best(std::declval<const cause_result_t&>())), std::optional<est_t>>);

    // ---- best_x: only for results whose solution has an x ----------------------------------------------------------
    static_assert(best_x_ok<root_result_t>);
    static_assert(!best_x_ok<search_result_t>);    // a search result has no x: read r->lo() and r->hi()
    static_assert(!best_x_ok<deriv_result_t>);

    // ---- detail::is_result_v: any std::expected<solution<S>, failure<F, UE>>, a search result included -------------
    static_assert(nxx::detail::is_result_v<root_result_t>);
    static_assert(nxx::detail::is_result_v<search_result_t>);
    static_assert(!nxx::detail::is_result_v<deriv_result_t>);
    static_assert(!nxx::detail::is_result_v<std::expected<double, int>>);

    // ---- best in a constant expression -----------------------------------------------------------------------------
    constexpr est_t         one_est { 1.0, 0.0, 0.0, std::nullopt };
    constexpr root_result_t ce_solved { nxx::solution<est_t> {
        one_est, nxx::counters { 1, 1 }, nxx::algo::user_first, nxx::stop_reason::exact_zero } };
    constexpr root_result_t ce_starved { std::unexpected(
        nxx::failure<est_t> { nxx::errc::budget_exhausted, nxx::algo::user_first, nxx::counters { 1, 1 }, one_est, {} }) };
    constexpr root_result_t ce_rejected { std::unexpected(nxx::failure<est_t> { nxx::errc::invalid_input, nxx::algo::none, {}, {}, {} }) };
    static_assert(nxx::best(ce_solved) == std::optional<est_t>(one_est));
    static_assert(nxx::best(ce_starved) == std::optional<est_t>(one_est));
    static_assert(!nxx::best(ce_rejected).has_value());

    constexpr auto sq2      = [](double x) { return x * x - 2.0; };    // root sqrt(2)
    constexpr auto sq_plus1 = [](double x) { return x * x + 1.0; };    // no real root
    constexpr auto twice    = [](double x) { return 2.0 * x; };        // the derivative of both

    // The alternatives of the chains below: newton fails at f'(0) = 0, secant within 5 iterations, bisection succeeds on
    // x^2 - 2 and fails on x^2 + 1.
    constexpr auto newton0  = r::newton {}.with_derivative(twice).on(0.0);
    constexpr auto secant5  = r::secant {}.with_budget(5).on(0.0);
    constexpr auto bisect02 = r::bisection {}.on(nxx::bracket { 0.0, 2.0 });

    // best is the solution's estimate on success, and the failure's best estimate (possibly none) on failure.
    template<class R>
    void check_best(const R& res)
    {
        const auto b = nxx::best(res);
        if (res) {
            CHECK(b == std::optional<est_t>(static_cast<const est_t&>(*res)));
            CHECK(nxx::best_x(res) == std::optional<double>(res->x));
        }
        else {
            CHECK(b == res.error().best);
            CHECK(nxx::best_x(res) == (b ? std::optional<double>(b->x) : std::nullopt));
        }
    }
}    // namespace

TEST_SUITE("roots")
{
    TEST_CASE("results: best on success is the solution's estimate")
    {
        const auto res = r::brent {}(sq2, { 1.0, 2.0 });
        static_assert(std::is_same_v<decltype(nxx::best(res)), std::optional<est_t>>);
        const auto b = nxx::best(res);
        CHECK(res.has_value());
        if (res && b) {
            CHECK(b->x == res->x);
            CHECK(b->fx == res->fx);
            CHECK(b->uncertainty == res->uncertainty);
            CHECK(b->enclosure == res->enclosure);
        }
        else
            FAIL_CHECK("brent failed on x^2 - 2, or best gave nothing for a success");
    }

    TEST_CASE("results: best on a failure is the failure's best estimate, or nothing if no evaluation succeeded")
    {
        // Starved: a failure with a best estimate.
        const auto starved = r::bisection { nxx::never {} }.with_budget(3)(sq2, nxx::bracket { 0.0, 2.0 });
        CHECK_FALSE(starved.has_value());
        if (!starved) {
            CHECK(starved.error().code == nxx::errc::budget_exhausted);
            CHECK(starved.error().best.has_value());
        }
        CHECK(nxx::best(starved).has_value());
        check_best(starved);

        // Equal ends: rejected before any evaluation, so no best estimate.
        const auto rejected = r::bisection {}.on(std::pair { 1.0, 1.0 })(sq2);
        CHECK_FALSE(rejected.has_value());
        if (!rejected) {
            CHECK(rejected.error().code == nxx::errc::invalid_input);
            CHECK(rejected.error().used == nxx::counters {});
        }
        CHECK(nxx::best(rejected) == std::nullopt);
        check_best(rejected);

        // A fallible callback: the cause type changes nothing.
        const auto g = [](double x) -> std::expected<double, user_error> {
            if (x < 0.0) return std::unexpected(user_error::domain);
            return x * x - 2.0;
        };
        const auto refused = r::bisection {}(g, nxx::bracket { -1.0, 2.0 });
        CHECK_FALSE(refused.has_value());
        if (!refused) { CHECK(refused.error().code == nxx::errc::callback_failed); }
        check_best(refused);
        check_best(r::bisection {}(g, nxx::bracket { 1.0, 2.0 }));
    }

    TEST_CASE("results: best through first_of, on success and on total failure")
    {
        const auto chain = nxx::first_of(newton0, secant5, bisect02);

        const auto ok = chain(sq2);
        CHECK(ok.has_value());
        check_best(ok);

        const auto failed = chain(sq_plus1);
        CHECK_FALSE(failed.has_value());
        CHECK(nxx::best(failed).has_value());    // merged over the alternatives
        check_best(failed);
    }

    TEST_CASE("results: best through any_solver, on success, on total failure and on an empty chain")
    {
        const solver_t chain = nxx::first_of(std::vector<solver_t> { newton0, secant5, bisect02 });
        const fn_t     f     = sq2;
        const fn_t     g     = sq_plus1;

        const auto ok = chain(f);
        static_assert(std::is_same_v<decltype(nxx::best(ok)), std::optional<est_t>>);
        CHECK(ok.has_value());
        check_best(ok);

        const auto failed = chain(g);
        CHECK_FALSE(failed.has_value());
        CHECK(nxx::best(failed).has_value());
        check_best(failed);

        const auto empty = nxx::first_of(std::vector<solver_t> {})(f);
        CHECK_FALSE(empty.has_value());
        CHECK(nxx::best(empty) == std::nullopt);
    }
}
