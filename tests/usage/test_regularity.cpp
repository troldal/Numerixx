// Regularity (DESIGN §3.2, §9.1; spike exit criterion 10): every solver, curried solver, chain, returned function,
// criterion, refined type and any_solver is std::copyable, also when it holds lambdas with captures (whose own copy
// assignment is deleted); std::semiregular where every part is default-constructible. At run time, a chain holding
// capturing lambdas is copy-ASSIGNED into another object of the same type, which then gives the source's result.

#include <numerixx/core/any_solver.hpp>
#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <functional>
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

    template<class T>
    inline constexpr bool copyable_v = std::copyable<std::remove_cvref_t<T>>;
    template<class T>
    inline constexpr bool semiregular_v = std::semiregular<std::remove_cvref_t<T>>;

    constexpr auto f_plain = [](double x) { return x * x - 2.0; };

    template<class R>
    bool same_result(const R& a, const R& b)
    {
        if (a.has_value() != b.has_value()) return false;
        if (a)
            return a->x == b->x && a->fx == b->fx && a->uncertainty == b->uncertainty && a->enclosure == b->enclosure &&
                   a->used == b->used && a->by == b->by && a->how == b->how;
        const auto& ea = a.error();
        const auto& eb = b.error();
        return ea.code == eb.code && ea.where == eb.where && ea.used == eb.used && ea.best.has_value() == eb.best.has_value() &&
               (!ea.best || ea.best->x == eb.best->x);
    }

    // Factories: every call returns the same closure types, so two results with different captures have one type.
    auto make_chain(double k, double x0, double lo, double hi)
    {
        const auto df   = [k](double x) { return k * x; };
        const auto proj = [lo, hi](double x) { return std::clamp(x, lo, hi); };
        return nxx::first_of(r::newton {}.with_derivative(df).on(x0),
                             r::secant {}.with_projection(proj).on(x0),
                             r::brent {}.on({ lo, hi }));
    }

    auto make_pipeline(double k, std::pair<double, double> window, nxx::tolerance<double> tol)    // a refined tolerance (DESIGN §6.2)
    {
        const auto df = [k](double x) { return k * x; };
        return nxx::then(r::bisection { nxx::width_tol { tol } }.on(window), r::newton {}.with_derivative(df));
    }

    auto make_warm(double k, std::pair<double, double> window, nxx::max_iterations budget)
    {
        const auto df = [k](double x) { return k * x; };
        return nxx::warm_fallback(r::bisection {}.with_budget(budget).on(window), r::newton {}.with_derivative(df));
    }

    auto make_derivative(double k)
    {
        return d::derivative_of([k](double x) { return k * x * x; });
    }

    auto make_observed(double lo, double hi, int& seen)
    {
        return r::secant {}.with_projection([lo, hi](double x) { return std::clamp(x, lo, hi); }).with_observer([&seen](const auto&) {
            ++seen;
        });
    }

    // A capture whose copy may throw but whose move cannot: the copyable box copies first, then moves in place (DESIGN §3.2).
    auto make_table_projection(std::vector<double> limits)
    {
        return r::secant {}.with_projection([limits](double x) { return std::clamp(x, limits.front(), limits.back()); });
    }

    // A callable whose copy throws on demand, for the copyable box's exception guarantee.
    struct flaky_copy
    {
        std::vector<double> offset;
        bool*               fail;

        flaky_copy(std::vector<double> o, bool* f) : offset(std::move(o)), fail(f) {}
        flaky_copy(const flaky_copy& other) : offset(other.offset), fail(other.fail)
        {
#if defined(__cpp_exceptions)
            if (*fail) throw 1;
#endif
        }
        flaky_copy(flaky_copy&&) noexcept        = default;
        flaky_copy& operator=(const flaky_copy&) = delete;
        double      operator()(double x) const { return x + offset.front(); }
    };

    using fn_t     = std::function<double(double)>;
    using solver_t = nxx::any_solver<fn_t, r::root_estimate<double>>;
}    // namespace

// ---- Compile-time: types without lambdas -------------------------------------------------------------------------
static_assert(semiregular_v<r::bisection<>>);
static_assert(semiregular_v<r::brent<>>);
static_assert(semiregular_v<r::secant<>>);
static_assert(semiregular_v<r::newton<>>);
static_assert(semiregular_v<r::expand<>>);
static_assert(copyable_v<decltype(r::bisection { nxx::width_tol { 1e-6 } })>);
static_assert(copyable_v<decltype(r::brent { nxx::width_tol { 1e-6 } })>);
static_assert(copyable_v<decltype(r::secant { nxx::x_tol { 1e-6 } })>);
static_assert(copyable_v<decltype(r::newton { nxx::x_tol { 1e-6 } })>);
static_assert(copyable_v<decltype(r::bisection {}.with_budget(10))>);
static_assert(semiregular_v<decltype(r::brent {}.with_stop(nxx::never {}))>);

// Criteria.
static_assert(copyable_v<nxx::x_tol<double>>);
static_assert(copyable_v<nxx::width_tol<double>>);
static_assert(copyable_v<nxx::f_tol<double>>);
static_assert(copyable_v<nxx::max_evaluations>);
static_assert(copyable_v<nxx::min_iterations>);
static_assert(semiregular_v<nxx::floored_width>);
static_assert(semiregular_v<nxx::step_tol<3, 5>>);
static_assert(semiregular_v<nxx::step_tol<7, 10>>);
static_assert(semiregular_v<nxx::never>);
static_assert(copyable_v<decltype(nxx::x_tol { 1e-6 } || nxx::max_evaluations { 10 })>);
static_assert(copyable_v<decltype(nxx::floored_width {} && nxx::f_tol { 1e-6 })>);
static_assert(semiregular_v<decltype(nxx::never {} || nxx::floored_width {})>);

// Refined types: copyable, and without a default value (an illegal state is not representable).
static_assert(copyable_v<nxx::tolerance<double>>);
static_assert(copyable_v<nxx::abs_tolerance<double>>);
static_assert(copyable_v<nxx::rel_tolerance<float>>);
static_assert(copyable_v<nxx::evaluation_budget>);
static_assert(copyable_v<nxx::max_iterations>);
static_assert(copyable_v<nxx::bracket<double>>);
static_assert(!std::default_initializable<nxx::tolerance<double>>);
static_assert(!std::default_initializable<nxx::max_iterations>);
static_assert(!std::default_initializable<nxx::bracket<double>>);

// Records and the rest of the vocabulary.
static_assert(copyable_v<r::root_estimate<double>>);
static_assert(copyable_v<r::sign_bracket<double>>);
static_assert(copyable_v<r::clamp_to<double>>);
static_assert(copyable_v<nxx::solution<r::root_estimate<double>>>);
static_assert(copyable_v<nxx::failure<r::root_estimate<double>, int>>);
static_assert(copyable_v<nxx::result<r::root_estimate<double>>>);
static_assert(semiregular_v<nxx::counters>);
static_assert(semiregular_v<nxx::none>);
static_assert(semiregular_v<nxx::no_derivative>);
static_assert(semiregular_v<nxx::no_projection>);
static_assert(semiregular_v<nxx::no_observer>);

// Derivative policies and functions.
static_assert(semiregular_v<d::numeric<>>);
static_assert(copyable_v<decltype(d::numeric { d::central_1_2, d::relative { 1e-6 } })>);
static_assert(copyable_v<decltype(d::derivative_of(f_plain))>);

// any_solver: copyable, deliberately not default-constructible (never empty, DESIGN §6.10).
static_assert(copyable_v<solver_t>);
static_assert(!std::default_initializable<solver_t>);
static_assert(!semiregular_v<solver_t>);
static_assert(copyable_v<nxx::any_solver<fn_t, r::root_estimate<double>, int>>);

TEST_SUITE("usage")
{
    TEST_CASE("copyable: solvers, curried solvers and combinators holding capturing lambdas")
    {
        const double k    = rt(2.0);
        const double lo   = rt(0.0);
        const double hi   = rt(2.0);
        int          seen = 0;

        const auto          df   = [k](double x) { return k * x; };
        const auto          proj = [lo, hi](double x) { return std::clamp(x, lo, hi); };
        const auto          obs  = [&seen](const auto&) { ++seen; };
        std::vector<double> limits { lo, hi };
        const auto          table = [limits](double x) { return std::clamp(x, limits.front(), limits.back()); };

        // The premise: these lambdas cannot be copy-assigned themselves (on the unqualified closure types: decltype of
        // a const variable is const, and no const type is assignable).
        static_assert(!std::is_copy_assignable_v<std::remove_const_t<decltype(df)>>);
        static_assert(!std::is_copy_assignable_v<std::remove_const_t<decltype(proj)>>);
        static_assert(!std::is_copy_assignable_v<std::remove_const_t<decltype(obs)>>);
        static_assert(!std::is_copy_assignable_v<std::remove_const_t<decltype(table)>>);

        // Solvers.
        const auto nt = r::newton {}.with_derivative(df);
        static_assert(copyable_v<decltype(nt)>);
        const auto observed = r::secant {}.with_projection(proj).with_observer(obs);
        static_assert(copyable_v<decltype(observed)>);
        static_assert(copyable_v<decltype(r::newton {}.with_derivative(df).with_projection(proj).with_observer(obs))>);
        static_assert(copyable_v<decltype(r::bisection {}.with_observer(obs))>);
        static_assert(copyable_v<decltype(r::brent {}.with_observer(obs))>);
        static_assert(copyable_v<decltype(r::secant {}.with_projection(table))>);
        static_assert(copyable_v<decltype(r::newton {}.with_derivative(fn_t { df }))>);

        // Curried solvers.
        static_assert(copyable_v<decltype(nt.on(1.0))>);
        static_assert(copyable_v<decltype(r::brent {}.on({ lo, hi }))>);
        static_assert(copyable_v<decltype(r::bisection {}.on(nxx::bracket { 0.0, 2.0 }))>);
        static_assert(copyable_v<decltype(r::brent {}.on(nxx::bracket<double>::make(lo, hi)))>);
        static_assert(copyable_v<decltype(r::expand {}.on(nxx::bracket { 2.0, 2.5 }))>);

        // Combinators.
        const auto chain = nxx::first_of(nt.on(1.0), observed.on(1.0), r::brent {}.on({ lo, hi }));
        static_assert(copyable_v<decltype(chain)>);
        const auto pipeline = nxx::then(r::expand {}.on(nxx::bracket { 2.0, 2.5 }), r::bisection { nxx::width_tol { 1e-4 } }, nt);
        static_assert(copyable_v<decltype(pipeline)>);
        const auto warm = nxx::warm_fallback(r::bisection {}.with_budget(3).on({ lo, hi }), nt);
        static_assert(copyable_v<decltype(warm)>);
        const std::uint32_t limit  = 100;
        const auto          policy = [limit](const auto& e) { return e.used.evaluations < limit; };
        const auto          with   = nxx::first_of_with(policy, nt.on(1.0), r::brent {}.on({ lo, hi }));
        static_assert(copyable_v<decltype(with)>);

        // Returned functions and instrumentation.
        const auto dfn = d::derivative_of([k](double x) { return k * x * x; });
        static_assert(copyable_v<decltype(dfn)>);
        static_assert(copyable_v<decltype(d::derivative_of(table))>);
        std::uint32_t calls   = 0;
        const auto    counted = nxx::fn::counted(proj, calls);
        static_assert(copyable_v<decltype(counted)>);

        // Run-time chains built from these.
        const solver_t from_chain = chain;
        static_assert(copyable_v<decltype(from_chain)>);
        const auto runtime = nxx::first_of(std::vector<solver_t> { chain, r::bisection {}.on({ lo, hi }) });
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(runtime)>, solver_t>);

        CHECK(chain(f_plain).has_value());
    }

    TEST_CASE("semiregular solvers: default-construct, then assign a configured value")
    {
        r::bisection<> b0;
        b0                 = r::bisection {}.with_budget(nxx::max_iterations { 3 });
        const auto starved = b0(f_plain, { rt(0.0), rt(2.0) });
        CHECK_FALSE(starved.has_value());
        if (!starved) {
            CHECK(starved.error().code == nxx::errc::budget_exhausted);
            CHECK(starved.error().used.iterations == 3);
        }

        r::brent<> br;
        const auto fine = br(f_plain, { rt(0.0), rt(2.0) });
        CHECK(fine.has_value());
    }

    TEST_CASE("copy-assigning a first_of chain that holds capturing lambdas")
    {
        auto a = make_chain(2.0, 1.0, 0.0, 2.0);    // correct f': newton succeeds
        auto b = make_chain(0.0, 1.0, 0.0, 2.0);    // f' == 0: newton fails, the projected secant takes over
        static_assert(std::is_same_v<decltype(a), decltype(b)>);
        const auto ra = a(f_plain);
        const auto rb = b(f_plain);
        CHECK(ra.has_value());
        CHECK(rb.has_value());
        if (ra && rb) {
            CHECK(ra->by == r::algos::newton);
            CHECK(rb->by != r::algos::newton);
        }
        CHECK_FALSE(same_result(ra, rb));

        b = a;    // copy assignment: the captures are destroyed and reconstructed
        CHECK(same_result(b(f_plain), ra));
        CHECK(same_result(a(f_plain), ra));    // the source is unchanged

        auto c = make_chain(0.0, 1.0, 0.0, 2.0);
        c      = std::as_const(b);
        CHECK(same_result(c(f_plain), ra));
    }

    TEST_CASE("copy-assigning then_t and warm_fallback_t chains that hold capturing lambdas")
    {
        const auto tol = nxx::tolerance<double>::make(rt(1e-3));
        CHECK(tol.has_value());
        if (tol) {
            auto       pa = make_pipeline(2.0, { 1.0, 2.0 }, *tol);    // bisection, then newton polish
            auto       pb = make_pipeline(0.0, { 1.0, 2.0 }, *tol);    // newton fails with a zero derivative
            const auto ra = pa(f_plain);
            CHECK(ra.has_value());
            CHECK_FALSE(pb(f_plain).has_value());
            pb = pa;
            CHECK(same_result(pb(f_plain), ra));
        }

        const auto three = nxx::max_iterations::make(3);
        CHECK(three.has_value());
        if (three) {
            auto       wa = make_warm(2.0, { 1.0, 2.0 }, *three);    // starved bisection, then newton from its best
            auto       wb = make_warm(0.0, { 1.0, 2.0 }, *three);
            const auto ra = wa(f_plain);
            CHECK(ra.has_value());
            if (ra) { CHECK(ra->by == r::algos::newton); }
            CHECK_FALSE(wb(f_plain).has_value());
            wb = wa;
            CHECK(same_result(wb(f_plain), ra));
        }
    }

    TEST_CASE("copy-assigning a derivative_fn that holds a capturing lambda")
    {
        auto       da = make_derivative(1.0);    // d/dx x^2 = 2x
        auto       db = make_derivative(3.0);    // d/dx 3x^2 = 6x
        const auto va = da(rt(2.0));
        CHECK(std::abs(va.value_or(0.0) - 4.0) < 1e-8);
        CHECK(std::abs(db(rt(2.0)).value_or(0.0) - 12.0) < 1e-8);
        db = da;
        CHECK(db(rt(2.0)).value_or(0.0) == va.value_or(1.0));
    }

    TEST_CASE("copy-assigning a solver with a capturing projection and a capturing observer")
    {
        int        seen_a = 0;
        int        seen_b = 0;
        auto       sa     = make_observed(0.0, 3.0, seen_a);
        auto       sb     = make_observed(0.0, 1.0, seen_b);    // the root sqrt(2) is outside [0, 1]: stalled
        const auto ra     = sa(f_plain, rt(0.1));
        CHECK(ra.has_value());
        const auto rb = sb(f_plain, rt(0.5));
        CHECK_FALSE(rb.has_value());
        CHECK(seen_a > 0);
        CHECK(seen_b > 0);

        sb                  = sa;
        const int  before_a = seen_a;
        const int  before_b = seen_b;
        const auto rb2      = sb(f_plain, rt(0.1));
        CHECK(same_result(rb2, ra));
        if (ra) {
            CHECK(seen_a - before_a == static_cast<int>(ra->used.iterations));    // the copied observer: one call per iteration
        }
        CHECK(seen_b == before_b);
    }

    TEST_CASE("copy-assigning a solver whose capture is not nothrow-copyable")
    {
        auto       ta = make_table_projection({ 0.0, 3.0 });
        auto       tb = make_table_projection({ 0.0, 1.0 });
        const auto ra = ta(f_plain, rt(0.1));
        CHECK(ra.has_value());
        CHECK_FALSE(tb(f_plain, rt(0.1)).has_value());
        tb = ta;
        CHECK(same_result(tb(f_plain, rt(0.1)), ra));
    }

    TEST_CASE("copyable_box: the assignment strategy follows the capture")
    {
        const std::vector<double> xs { 1.0, 2.0 };
        const auto                by_double = [k = 2.0](double x) { return k * x; };
        const auto                by_vector = [v = xs](double x) { return x + v.front(); };    // a non-const member
        const auto                by_const  = [xs](double x) { return x + xs.front(); };       // a const member
        using by_double_t                   = std::remove_cvref_t<decltype(by_double)>;
        using by_vector_t                   = std::remove_cvref_t<decltype(by_vector)>;
        using by_const_t                    = std::remove_cvref_t<decltype(by_const)>;
        static_assert(nxx::detail::box_kind_v<by_double_t> == 2);    // nothrow copy: in place
        static_assert(nxx::detail::box_kind_v<by_vector_t> == 3);    // throwing copy, nothrow move: copy first
        static_assert(nxx::detail::box_kind_v<flaky_copy> == 3);
        // Capturing a const variable by copy makes a const member, so the closure's move copies it and may throw: the
        // box falls back to std::optional.
        static_assert(nxx::detail::box_kind_v<by_const_t> == 4);

        bool                                  fail = false;
        nxx::detail::copyable_box<flaky_copy> a { flaky_copy { { 1.0 }, &fail } };
        nxx::detail::copyable_box<flaky_copy> b { flaky_copy { { 2.0 }, &fail } };
        a = b;
        CHECK((*a)(0.0) == 2.0);
#if defined(__cpp_exceptions)
        nxx::detail::copyable_box<flaky_copy> c { flaky_copy { { 3.0 }, &fail } };
        fail = true;
        CHECK_THROWS(a = c);
        fail = false;
        CHECK((*a)(0.0) == 2.0);    // a throwing copy leaves the box unchanged
#endif
    }

    TEST_CASE("copy-assigning any_solver values")
    {
        const fn_t fn = f_plain;    // convert f once (DESIGN §6.10)
        solver_t   s1 = make_chain(2.0, 1.0, 0.0, 2.0);
        solver_t   s2 = make_chain(0.0, 1.0, 0.0, 2.0);
        const auto r1 = s1(fn);
        CHECK(r1.has_value());
        CHECK_FALSE(same_result(s2(fn), r1));
        s2 = s1;
        CHECK(same_result(s2(fn), r1));
        CHECK(same_result(r1, make_chain(2.0, 1.0, 0.0, 2.0)(fn)));    // the same as the static chain

        solver_t s3 = nxx::first_of(std::vector<solver_t> { s1, r::bisection {}.on({ 0.0, 2.0 }) });
        s3          = s2;
        CHECK(same_result(s3(fn), r1));
    }
}
