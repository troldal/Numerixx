// Determinism (DESIGN §9.1, §6.1): given the same values of f, a solver takes the same path on every platform, because
// the library's stop tests, step rules and midpoints use only the four basic operations and exact or correctly
// rounded helpers, with floating-point contraction off. Checked two ways:
//   - repeated calls of one solver value, of a copy and of a copy-assigned value give bit-identical results;
//   - golden values: a table of solves whose results are hard-coded as hexadecimal literals with exact counters.
//
// The golden functions use only correctly rounded operations (+, -, *, /, sqrt and an explicit std::fma) and never a
// product followed by an addition, so Clang's contraction of the test's own code cannot change their values: f is
// the same on every compiler, and so must be the solvers' results.
//
// Regenerating the table: compile this file with -DNXX_PRINT_GOLDEN (GCC 16 was used) and run the test case
// "determinism: golden values"; it prints the rows instead of checking them. Then verify them on every toolchain.
#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iterator>
#include <limits>
#include <ranges>
#include <utility>
#include <vector>

#if defined(NXX_PRINT_GOLDEN)
#    include <charconv>
#    include <cstdio>
#    include <string>
#endif

namespace
{
    namespace r     = nxx::roots;
    namespace algos = nxx::roots::algos;
    using nxx::errc;
    using nxx::stop_reason;

    constexpr double inf = std::numeric_limits<double>::infinity();

    // Contraction-free functions of correctly rounded operations.
    constexpr auto fn_sqrt2  = [](double x) { return x - 2.0 / x; };    // root sqrt(2)
    constexpr auto dfn_sqrt2 = [](double x) { return 1.0 + 2.0 / (x * x); };
    constexpr auto fn_cubic  = [](double x) { return (x - 0.3) * (x + 2.1) * (x - 7.0); };    // roots 0.3, -2.1, 7
    constexpr auto dfn_cubic = [](double x) {
        const double a = x - 0.3;
        const double b = x + 2.1;
        const double c = x - 7.0;
        return std::fma(a, b, std::fma(a, c, b * c));
    };
    constexpr auto fn_root    = [](double x) { return std::sqrt(x) - 1.7; };    // root 2.89
    constexpr auto dfn_root   = [](double x) { return 0.5 / std::sqrt(x); };
    constexpr auto fn_horner  = [](double x) { return std::fma(x, std::fma(x, x, -3.0), 1.0); };    // x^3 - 3x + 1
    constexpr auto dfn_horner = [](double x) { return std::fma(3.0 * x, x, -3.0); };

    std::uint64_t bits(double v) { return std::bit_cast<std::uint64_t>(v); }

    // A result, flattened: the solution (or the failure's best estimate) and everything else a solve reports.
    struct estimate_record
    {
        double x           = 0.0;
        double fx          = 0.0;
        double uncertainty = 0.0;
        bool   enclosed    = false;
        double lo          = 0.0;
        double flo         = 0.0;
        double hi          = 0.0;
        double fhi         = 0.0;
    };

    struct record
    {
        const char*     name = "";
        bool            ok   = false;
        nxx::algo       by {};      // solution::by, or failure::where
        stop_reason     how {};     // success only
        errc            code {};    // failure only
        std::uint32_t   iterations   = 0;
        std::uint32_t   evaluations  = 0;
        bool            has_estimate = false;    // the solution, or failure::best
        estimate_record est {};
    };

    estimate_record flatten(const r::root_estimate<double>& e)
    {
        estimate_record out { e.x, e.fx, e.uncertainty, false, 0.0, 0.0, 0.0, 0.0 };
        if (e.enclosure) {
            out.enclosed = true;
            out.lo       = e.enclosure->lo();
            out.flo      = e.enclosure->flo();
            out.hi       = e.enclosure->hi();
            out.fhi      = e.enclosure->fhi();
        }
        return out;
    }

    template<class R>
    record observe(const char* name, const R& res)
    {
        record out;
        out.name = name;
        out.ok   = res.has_value();
        if (res) {
            out.by           = res->by;
            out.how          = res->how;
            out.iterations   = res->used.iterations;
            out.evaluations  = res->used.evaluations;
            out.has_estimate = true;
            out.est          = flatten(*res);
        }
        else {
            const auto& e   = res.error();
            out.by          = e.where;
            out.code        = e.code;
            out.iterations  = e.used.iterations;
            out.evaluations = e.used.evaluations;
            if (e.best) {
                out.has_estimate = true;
                out.est          = flatten(*e.best);
            }
        }
        return out;
    }

    // Bit-for-bit equality of two records (+0.0 and -0.0 differ; so would NaN payloads).
    void check_identical(const record& got, const record& want)
    {
        INFO(want.name);
        CHECK(got.ok == want.ok);
        CHECK(got.by == want.by);
        CHECK(got.how == want.how);
        CHECK(got.code == want.code);
        CHECK(got.iterations == want.iterations);
        CHECK(got.evaluations == want.evaluations);
        CHECK(got.has_estimate == want.has_estimate);
        CHECK(bits(got.est.x) == bits(want.est.x));
        CHECK(bits(got.est.fx) == bits(want.est.fx));
        CHECK(bits(got.est.uncertainty) == bits(want.est.uncertainty));
        CHECK(got.est.enclosed == want.est.enclosed);
        CHECK(bits(got.est.lo) == bits(want.est.lo));
        CHECK(bits(got.est.flo) == bits(want.est.flo));
        CHECK(bits(got.est.hi) == bits(want.est.hi));
        CHECK(bits(got.est.fhi) == bits(want.est.fhi));
    }

    // Calls s(args...) on s itself (three times), on a copy and on a copy-assigned value; every result must be
    // bit-identical to the first, which must succeed or fail as expected (both paths are covered).
    template<class S, class Other, class... A>
    void check_repeatable(const char* name, bool expect_ok, const S& s, const Other& other, const A&... args)
    {
        const record first = observe(name, s(args...));
        INFO(name);
        CHECK(first.ok == expect_ok);
        check_identical(observe(name, s(args...)), first);
        const S copy = s;
        check_identical(observe(name, copy(args...)), first);
        S assigned = other;
        assigned   = s;
        check_identical(observe(name, assigned(args...)), first);
        check_identical(observe(name, s(args...)), first);    // the copies did not disturb the original
    }

    // Every golden solve, in table order.
    std::vector<record> golden_runs()
    {
        std::vector<record> v;
        v.push_back(observe("bisection sqrt2 [1, 2]", r::bisection {}(fn_sqrt2, { 1.0, 2.0 })));
        v.push_back(observe("bisection width_tol{1e-6} cubic [0, 1]", r::bisection { nxx::width_tol { 1e-6 } }(fn_cubic, { 0.0, 1.0 })));
        v.push_back(observe("bisection budget 10 root [0, 10]", r::bisection {}.with_budget(10)(fn_root, { 0.0, 10.0 })));
        v.push_back(observe("bisection no sign change cubic [3, 4]", r::bisection {}(fn_cubic, { 3.0, 4.0 })));
        v.push_back(observe("brent sqrt2 [1, 2]", r::brent {}(fn_sqrt2, { 1.0, 2.0 })));
        v.push_back(observe("brent width_tol{1e-10} cubic [-3, 0]", r::brent { nxx::width_tol { 1e-10 } }(fn_cubic, { -3.0, 0.0 })));
        v.push_back(observe("brent horner [0, 1]", r::brent {}(fn_horner, { 0.0, 1.0 })));
        v.push_back(observe("brent root [0, 10]", r::brent {}(fn_root, { 0.0, 10.0 })));
        v.push_back(
            observe("brent width_tol{0, 1e-12} horner [1, 2]", r::brent { nxx::width_tol { 0.0, 1e-12 } }(fn_horner, { 1.0, 2.0 })));
        v.push_back(observe("secant sqrt2 from 1", r::secant {}(fn_sqrt2, 1.0)));
        v.push_back(observe("secant x_tol{1e-12, 1e-10} horner from 2", r::secant { nxx::x_tol { 1e-12, 1e-10 } }(fn_horner, 2.0)));
        v.push_back(observe("secant cubic from 6", r::secant {}(fn_cubic, 6.0)));
        v.push_back(observe("secant f_tol{1e-9} root from 1", r::secant { nxx::f_tol { 1e-9 } }(fn_root, 1.0)));
        v.push_back(observe("newton sqrt2 from 1", r::newton {}.with_derivative(dfn_sqrt2)(fn_sqrt2, 1.0)));
        v.push_back(observe("newton cubic from 6", r::newton {}.with_derivative(dfn_cubic)(fn_cubic, 6.0)));
        v.push_back(observe("newton horner from -3", r::newton {}.with_derivative(dfn_horner)(fn_horner, -3.0)));
        v.push_back(observe("newton horner from 1 (zero derivative)", r::newton {}.with_derivative(dfn_horner)(fn_horner, 1.0)));
        v.push_back(observe("newton x_tol{1e-14} root from 1", r::newton { nxx::x_tol { 1e-14 } }.with_derivative(dfn_root)(fn_root, 1.0)));
        v.push_back(observe("then(bisection width_tol{1e-3}, secant) sqrt2 [1, 2]",
                            nxx::then(r::bisection { nxx::width_tol { 1e-3 } }.with_budget(100).on({ 1.0, 2.0 }), r::secant {})(fn_sqrt2)));
        v.push_back(observe("first_of(newton, secant, brent) horner",
                            nxx::first_of(r::newton {}.with_derivative(dfn_horner).on(1.0),
                                          r::secant {}.on(1.0),
                                          r::brent {}.on({ 0.0, 1.0 }))(fn_horner)));
        v.push_back(observe("then(expand, brent) root [1, 2]", nxx::then(r::expand {}.on({ 1.0, 2.0 }), r::brent {})(fn_root)));
        v.push_back(observe("first_of(newton, brent) horner, total failure",
                            nxx::first_of(r::newton {}.with_derivative(dfn_horner).on(1.0), r::brent {}.on({ 3.0, 4.0 }))(fn_horner)));
        return v;
    }

#if defined(NXX_PRINT_GOLDEN)
    std::string hex(double v)
    {
        if (v == inf) return "inf";
        if (v == -inf) return "-inf";
        char       buf[64];
        const bool neg = std::signbit(v);
        const auto res = std::to_chars(buf, buf + sizeof buf, neg ? -v : v, std::chars_format::hex);
        return std::string(neg ? "-0x" : "0x") + std::string(buf, res.ptr);
    }

    const char* name_of(nxx::algo a)
    {
        if (a == algos::bisection) return "algos::bisection";
        if (a == algos::brent) return "algos::brent";
        if (a == algos::secant) return "algos::secant";
        if (a == algos::newton) return "algos::newton";
        if (a == algos::expand) return "algos::expand";
        return "nxx::algo::none";
    }

    const char* name_of(stop_reason s)
    {
        switch (s) {
            case stop_reason::exact_zero:
                return "stop_reason::exact_zero";
            case stop_reason::criterion:
                return "stop_reason::criterion";
            case stop_reason::resolution_limit:
                return "stop_reason::resolution_limit";
            case stop_reason::algorithm:
                return "stop_reason::algorithm";
        }
        return "stop_reason{}";
    }

    std::string name_of(errc e)
    {
        switch (e) {
            case errc::no_sign_change:
                return "errc::no_sign_change";
            case errc::budget_exhausted:
                return "errc::budget_exhausted";
            case errc::stalled:
                return "errc::stalled";
            case errc::non_finite_value:
                return "errc::non_finite_value";
            case errc::zero_derivative:
                return "errc::zero_derivative";
            case errc::diverged:
                return "errc::diverged";
            default:
                return e == errc {} ? "errc{}" : "errc{" + std::to_string(static_cast<int>(e)) + "}";
        }
    }

    void print_row(const record& g)
    {
        const auto& e = g.est;
        std::printf("            { \"%s\", %s, %s, %s, %s, %u, %u, %s,\n              { %s, %s, %s, %s, %s, %s, %s, %s } },\n",
                    g.name,
                    g.ok ? "true" : "false",
                    name_of(g.by),
                    g.ok ? name_of(g.how) : "stop_reason{}",
                    name_of(g.code).c_str(),
                    unsigned { g.iterations },
                    unsigned { g.evaluations },
                    g.has_estimate ? "true" : "false",
                    hex(e.x).c_str(),
                    hex(e.fx).c_str(),
                    hex(e.uncertainty).c_str(),
                    e.enclosed ? "true" : "false",
                    hex(e.lo).c_str(),
                    hex(e.flo).c_str(),
                    hex(e.hi).c_str(),
                    hex(e.fhi).c_str());
    }
#endif
}    // namespace

TEST_SUITE("roots")
{
    TEST_CASE("determinism: repeated calls and copies give bit-identical results")
    {
        const auto fn_exp  = [](double x) { return std::exp(x) - 3.0; };    // not correctly rounded: fine on one platform
        const auto dfn_exp = [](double x) { return std::exp(x); };

        SUBCASE("bisection")
        {
            check_repeatable("bisection sqrt2", true, r::bisection {}, r::bisection {}, fn_sqrt2, std::pair { 1.0, 2.0 });
            const auto wt = r::bisection { nxx::width_tol { 1e-9, 1e-12 } };
            check_repeatable("bisection width_tol exp", true, wt, wt.with_budget(3), fn_exp, std::pair { 0.0, 3.0 });
            check_repeatable("bisection budget failure",
                             false,
                             r::bisection {}.with_budget(5),
                             r::bisection {}.with_budget(6),
                             fn_exp,
                             std::pair { 0.0, 3.0 });
            check_repeatable("bisection no sign change", false, r::bisection {}, r::bisection {}, fn_exp, std::pair { 2.0, 3.0 });
        }
        SUBCASE("brent")
        {
            check_repeatable("brent horner", true, r::brent {}, r::brent {}, fn_horner, std::pair { 0.0, 1.0 });
            const auto wt = r::brent { nxx::width_tol { 1e-8 } };
            check_repeatable("brent width_tol exp", true, wt, wt.with_budget(2), fn_exp, std::pair { -1.0, 4.0 });
            check_repeatable("brent budget failure",
                             false,
                             r::brent {}.with_budget(2),
                             r::brent {}.with_budget(3),
                             fn_exp,
                             std::pair { -1.0, 4.0 });
        }
        SUBCASE("secant")
        {
            check_repeatable("secant sqrt2", true, r::secant {}, r::secant {}, fn_sqrt2, 1.0);
            const auto xt = r::secant { nxx::x_tol { 1e-12, 1e-9 } };
            check_repeatable("secant x_tol exp", true, xt, xt.with_budget(1), fn_exp, 0.5);
            check_repeatable("secant clamped",
                             false,
                             r::secant {}.with_projection(r::clamp_to { 0.0, 1.0 }),
                             r::secant {}.with_projection(r::clamp_to { 0.0, 2.0 }),
                             fn_exp,
                             0.5);    // pinned at the edge: stalled
        }
        SUBCASE("newton")
        {
            const auto nt = r::newton {}.with_derivative(dfn_exp);
            check_repeatable("newton exp", true, nt, nt.with_budget(1), fn_exp, 3.0);
            const auto nh = r::newton {}.with_derivative(dfn_horner);
            check_repeatable("newton horner", true, nh, nh.with_budget(2), fn_horner, -3.0);
            check_repeatable("newton zero derivative", false, nh, nh.with_budget(2), fn_horner, 1.0);
        }
        SUBCASE("first_of chain")
        {
            const auto chain =
                nxx::first_of(r::newton {}.with_derivative(dfn_horner).on(1.0), r::secant {}.on(1.0), r::brent {}.on({ 0.0, 1.0 }));
            const auto other =
                nxx::first_of(r::newton {}.with_derivative(dfn_horner).on(2.0), r::secant {}.on(2.0), r::brent {}.on({ 1.0, 2.0 }));
            check_repeatable("first_of horner", true, chain, other, fn_horner);
            const auto failing  = nxx::first_of(r::newton {}.with_derivative(dfn_horner).on(1.0), r::brent {}.on({ 3.0, 4.0 }));
            const auto failing2 = nxx::first_of(r::newton {}.with_derivative(dfn_horner).on(1.0), r::brent {}.on({ 5.0, 6.0 }));
            check_repeatable("first_of total failure", false, failing, failing2, fn_horner);
        }
        SUBCASE("then pipeline")
        {
            const auto pipe  = nxx::then(r::bisection { nxx::width_tol { 1e-3 } }.with_budget(100).on({ 0.0, 3.0 }), r::secant {});
            const auto other = nxx::then(r::bisection { nxx::width_tol { 1e-3 } }.with_budget(100).on({ 0.5, 2.5 }), r::secant {});
            check_repeatable("then bisection-secant exp", true, pipe, other, fn_exp);
            const auto search  = nxx::then(r::expand {}.on({ 1.0, 2.0 }), r::brent {});
            const auto search2 = nxx::then(r::expand {}.on({ 1.0, 3.0 }), r::brent {});
            check_repeatable("then expand-brent root", true, search, search2, fn_root);
        }
        SUBCASE("steps_view: iterating the same view twice gives the same states")
        {
            const auto solver = r::brent {};
            auto       prob   = solver.prepare(std::cref(fn_horner), nxx::bracket<double>::make(0.0, 1.0));
            if (prob) {
                const nxx::steps_view      view { solver, *prob };
                std::vector<std::uint64_t> first;
                std::vector<std::uint64_t> second;
                for (const auto& st : view | std::views::take(40))
                    if (st) first.push_back(bits(solver.estimate(*st).x));
                for (const auto& st : view | std::views::take(40))
                    if (st) second.push_back(bits(solver.estimate(*st).x));
                CHECK(first.size() >= 3);
                CHECK(first == second);
            }
            else
                FAIL_CHECK("brent::prepare failed");
        }
    }

    TEST_CASE("determinism: golden values")
    {
        const std::vector<record> got = golden_runs();
#if defined(NXX_PRINT_GOLDEN)
        for (const record& g : got) print_row(g);
        CHECK(!got.empty());
#else
        // clang-format off
        // { name, ok, by/where, how, code, iterations, evaluations, has_estimate,
        //   { x, fx, uncertainty, enclosed, lo, flo, hi, fhi } }
        const record golden[] = {
            { "bisection sqrt2 [1, 2]", true, algos::bisection, stop_reason::criterion, errc{}, 50, 52, true,
              { 0x1.6a09e667f3bccp+0, -0x1p-52, 0x1p-50, true, 0x1.6a09e667f3bccp+0, -0x1p-52, 0x1.6a09e667f3bdp+0, 0x1.cp-50 } },
            { "bisection width_tol{1e-6} cubic [0, 1]", true, algos::bisection, stop_reason::criterion, errc{}, 20, 22, true,
              { 0x1.33334p-2, -0x1.9ba5e4b4a0409p-19, 0x1p-20, true, 0x1.3333p-2, 0x1.9ba5ddd2d7df4p-17, 0x1.33334p-2, -0x1.9ba5e4b4a0409p-19 } },
            { "bisection budget 10 root [0, 10]", false, algos::bisection, stop_reason{}, errc::budget_exhausted, 10, 12, true,
              { 0x1.72p+1, 0x1.817c2bb578p-13, 0x1.4p-7, true, 0x1.70cp+1, -0x1.60a7d17bdc4p-9, 0x1.72p+1, 0x1.817c2bb578p-13 } },
            { "bisection no sign change cubic [3, 4]", false, algos::bisection, stop_reason{}, errc::no_sign_change, 0, 2, true,
              { 0x1.8p+1, -0x1.b8a3d70a3d70ap+5, 0x1p+0, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
            { "brent sqrt2 [1, 2]", true, algos::brent, stop_reason::criterion, errc{}, 6, 8, true,
              { 0x1.6a09e667f3bcdp+0, 0x1p-52, 0x1.8p-51, true, 0x1.6a09e667f3bcap+0, -0x1.4p-50, 0x1.6a09e667f3bcdp+0, 0x1p-52 } },
            { "brent width_tol{1e-10} cubic [-3, 0]", true, algos::brent, stop_reason::criterion, errc{}, 9, 11, true,
              { -0x1.0ccccccccae76p+1, 0x1.4b4fa3d707a7ep-34, 0x1.b7cep-35, true, -0x1.0ccccccce6644p+1, -0x1.1775b28f7994ap-30, -0x1.0ccccccccae76p+1, 0x1.4b4fa3d707a7ep-34 } },
            { "brent horner [0, 1]", true, algos::brent, stop_reason::criterion, errc{}, 8, 10, true,
              { 0x1.63a1a7e0b7389p-2, 0x1.48ad7df8e4804p-54, 0x1p-51, true, 0x1.63a1a7e0b7389p-2, 0x1.48ad7df8e4804p-54, 0x1.63a1a7e0b7391p-2, -0x1.2f90a536ec39cp-50 } },
            { "brent root [0, 10]", true, algos::brent, stop_reason::exact_zero, errc{}, 3, 5, true,
              { 0x1.71eb851eb851ep+1, 0x0p+0, 0x1.71eb851eb851ep+1, true, 0x0p+0, -0x1.b333333333333p+0, 0x1.71eb851eb851ep+1, 0x0p+0 } },
            { "brent width_tol{0, 1e-12} horner [1, 2]", true, algos::brent, stop_reason::criterion, errc{}, 9, 11, true,
              { 0x1.8836fa2cf5039p+0, -0x1.fee5fe060c6b8p-54, 0x1.af4p-41, true, 0x1.8836fa2cf5039p+0, -0x1.fee5fe060c6b8p-54, 0x1.8836fa2cf5db3p+0, 0x1.b3bb14b8cd1f5p-39 } },
            { "secant sqrt2 from 1", true, algos::secant, stop_reason::criterion, errc{}, 7, 9, true,
              { 0x1.6a09e667f3bccp+0, -0x1p-52, 0x1p-52, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
            { "secant x_tol{1e-12, 1e-10} horner from 2", true, algos::secant, stop_reason::criterion, errc{}, 8, 10, true,
              { 0x1.8836fa2cf5039p+0, -0x1.fee5fe060c6b8p-54, 0x1p-49, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
            { "secant cubic from 6", true, algos::secant, stop_reason::exact_zero, errc{}, 7, 9, true,
              { 0x1.cp+2, 0x0p+0, 0x1.33558p-33, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
            { "secant f_tol{1e-9} root from 1", true, algos::secant, stop_reason::criterion, errc{}, 6, 8, true,
              { 0x1.71eb851eb8489p+1, -0x1.6p-46, 0x1.3066cb8p-26, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
            { "newton sqrt2 from 1", true, algos::newton, stop_reason::criterion, errc{}, 5, 11, true,
              { 0x1.6a09e667f3bccp+0, -0x1p-52, 0x1.c0ep-40, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
            { "newton cubic from 6", true, algos::newton, stop_reason::exact_zero, errc{}, 6, 13, true,
              { 0x1.cp+2, 0x0p+0, 0x1p-50, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
            { "newton horner from -3", true, algos::newton, stop_reason::criterion, errc{}, 6, 13, true,
              { -0x1.e11f642522d1cp+0, -0x1.549d103c20773p-51, 0x1.19aa8p-32, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
            { "newton horner from 1 (zero derivative)", false, algos::newton, stop_reason{}, errc::zero_derivative, 1, 2, true,
              { 0x1p+0, -0x1p+0, inf, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
            { "newton x_tol{1e-14} root from 1", true, algos::newton, stop_reason::exact_zero, errc{}, 5, 11, true,
              { 0x1.71eb851eb851fp+1, 0x0p+0, 0x1.7fd78p-33, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
            { "then(bisection width_tol{1e-3}, secant) sqrt2 [1, 2]", true, algos::secant, stop_reason::criterion, errc{}, 14, 17, true,
              { 0x1.6a09e667f3bccp+0, -0x1p-52, 0x1p-52, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
            { "first_of(newton, secant, brent) horner", true, algos::brent, stop_reason::criterion, errc{}, 59, 64, true,
              { 0x1.63a1a7e0b7389p-2, 0x1.48ad7df8e4804p-54, 0x1p-51, true, 0x1.63a1a7e0b7389p-2, 0x1.48ad7df8e4804p-54, 0x1.63a1a7e0b7391p-2, -0x1.2f90a536ec39cp-50 } },
            { "then(expand, brent) root [1, 2]", true, algos::brent, stop_reason::exact_zero, errc{}, 3, 5, true,
              { 0x1.71eb851eb851ep+1, 0x0p+0, 0x1.c7ae147ae1478p-1, true, 0x1p+1, -0x1.24a5332cfdd98p-2, 0x1.71eb851eb851ep+1, 0x0p+0 } },
            { "first_of(newton, brent) horner, total failure", false, algos::brent, stop_reason{}, errc::no_sign_change, 1, 4, true,
              { 0x1p+0, -0x1p+0, inf, false, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0 } },
        };
        // clang-format on
        CHECK(got.size() == std::size(golden));
        for (std::size_t i = 0; i < got.size() && i < std::size(golden); ++i) {
            CHECK(std::strcmp(got[i].name, golden[i].name) == 0);
            check_identical(got[i], golden[i]);
        }
#endif
    }
}
