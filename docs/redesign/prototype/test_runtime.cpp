// Feasibility prototype: run-time solver chains (nxx/runtime.hpp). A chain is built from configuration
// strings at run time and must give exactly the static chain's results. Does NOT include nxx/pipes.hpp,
// so it needs no FXT (also not for -fno-exceptions).
#ifndef NXX_FIX_UNWRAP_FAILURE
#  define NXX_FIX_UNWRAP_FAILURE   // newton + deriv::numeric needs it (see test_core.cpp)
#endif
#include "nxx/roots.hpp"
#include "nxx/deriv.hpp"
#include "nxx/runtime.hpp"
#include <array>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// ------------------------------------------------ allocation counter (replaces global operator new)
#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic ignored "-Wmismatched-new-delete"   // false positive: our operator new IS malloc
#endif
static std::size_t g_allocs = 0;
void* operator new(std::size_t n) {
    ++g_allocs;
    if (void* p = std::malloc(n ? n : 1)) return p;
    std::abort();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

namespace r = nxx::roots;
namespace d = nxx::deriv;

using fn_t     = std::function<double(double)>;
using est_t    = r::root_estimate<double>;
using solver_t = nxx::any_solver<fn_t, est_t>;
using result_t = solver_t::result_type;

constexpr auto df = [](double x) { return 2.0 * x; };

// the static chains of test_core.cpp [1] and [3]
constexpr auto static_chain = nxx::first_of(r::newton{}.with_derivative(df).on(0.0),        // zero_derivative
                                            r::secant{nxx::default_step{}, 5}.on(0.0),      // budget too small
                                            r::bisection{}.on(nxx::bracket{0.0, 2.0}));     // succeeds
constexpr auto static_chain_fd = nxx::first_of(r::newton{}.with_derivative(d::numeric{}).on(1.0),
                                               r::brent{}.on(nxx::bracket{0.0, 2.0}));

// ------------------------------------------------ compile-time properties
using bisect_t = decltype(r::bisection{}.on(nxx::bracket{0.0, 2.0}));
static_assert(std::is_copy_constructible_v<solver_t> && std::is_copy_assignable_v<solver_t>);
static_assert(!std::is_default_constructible_v<solver_t>);                       // never empty
static_assert(std::is_convertible_v<bisect_t, solver_t>);                        // a solver value IS an any_solver
static_assert(std::is_convertible_v<decltype(static_chain), solver_t>);          // so is a static chain
static_assert(!std::is_copy_assignable_v<decltype(static_chain)>);               // ...which is not assignable (closure)
static_assert(!std::is_constructible_v<solver_t, decltype(r::expand_out{}.on(nxx::bracket{0.0, 1.0}))>);   // Est differs
static_assert(!std::is_constructible_v<solver_t, int>);
static_assert(!std::is_nothrow_invocable_v<const solver_t&, const fn_t&>);        // std::function call is not noexcept
static_assert(std::is_same_v<decltype(nxx::first_of(std::declval<std::vector<solver_t>>())), solver_t>);
static_assert(std::is_same_v<decltype(nxx::first_of(std::declval<std::span<const solver_t>>())), solver_t>);
static_assert(std::is_same_v<decltype(nxx::first_of(std::declval<std::array<solver_t, 2>>())), solver_t>);
static_assert(std::is_same_v<decltype(nxx::first_of(std::declval<solver_t>())), solver_t>);   // core's identity overload

// ------------------------------------------------ run-time configuration -> solvers (no exceptions)
std::optional<solver_t> make_solver(std::string_view name) {
    if (name == "newton")    return r::newton{}.with_derivative(df).on(0.0);
    if (name == "newton_fd") return r::newton{}.with_derivative(d::numeric{}).on(1.0);
    if (name == "secant")    return r::secant{nxx::default_step{}, 5}.on(0.0);
    if (name == "bisection") return r::bisection{}.on(nxx::bracket{0.0, 2.0});
    if (name == "brent")     return r::brent{}.on(nxx::bracket{0.0, 2.0});
    return std::nullopt;
}
std::expected<solver_t, std::string> make_chain(const std::vector<std::string>& names) {
    std::vector<solver_t> alts;
    alts.reserve(names.size());
    for (const auto& n : names) {
        auto s = make_solver(n);
        if (!s) return std::unexpected("unknown solver '" + n + "'");
        alts.push_back(*s);
    }
    return nxx::first_of(std::move(alts));
}

// ------------------------------------------------ helpers
static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { std::printf("  CHECK FAILED line %d: %s\n", __LINE__, #cond); ++failures; } } while (0)

bool same_est(const est_t& a, const est_t& b) { return a.x == b.x && a.fx == b.fx && a.enclosure == b.enclosure; }
template<class R> bool same(const R& a, const R& b) {   // field-by-field, bit-exact doubles
    if (a.has_value() != b.has_value()) return false;
    if (a) return same_est(*a, *b) && a->used == b->used && a->by == b->by && a->how == b->how;
    const auto &ea = a.error(), &eb = b.error();
    return ea.code == eb.code && ea.where == eb.where && ea.used == eb.used && ea.cause == eb.cause &&
           ea.best.has_value() == eb.best.has_value() && (!ea.best || same_est(*ea.best, *eb.best));
}
const char* name(nxx::errc e) {
    switch (e) {
        case nxx::errc::invalid_input: return "invalid_input";
        case nxx::errc::no_sign_change: return "no_sign_change";
        case nxx::errc::budget_exhausted: return "budget_exhausted";
        case nxx::errc::stalled: return "stalled";
        case nxx::errc::zero_derivative: return "zero_derivative";
        case nxx::errc::callback_failed: return "callback_failed";
        default: return "other";
    }
}
const char* stdlib() {
#if defined(_LIBCPP_VERSION)
    return "libc++";
#elif defined(__GLIBCXX__)
    return "libstdc++";
#elif defined(_MSVC_STL_VERSION)
    return "MSVC STL";
#else
    return "?";
#endif
}

int main() {
    const fn_t f  = [](double x) { return x * x - 2.0; };
    const fn_t f2 = [](double x) { return x * x + 1.0; };   // no real root: every alternative fails
    const std::vector<std::string> config{"newton", "secant", "bisection"};

    std::printf("[R1] config {newton, secant, bisection} on x^2-2: ");
    auto chain = make_chain(config);
    CHECK(chain.has_value());
    const result_t rr = (*chain)(f);
    const result_t rs = static_chain(f);
    const result_t rl = static_chain([](double x) { return x * x - 2.0; });   // static chain on a plain lambda
    nxx::counters sum{};
    for (const auto& n : config) { auto one = (*make_solver(n))(f); sum = sum + (one ? one->used : one.error().used); }
    CHECK(rr && rr->by == nxx::algo::bisection && same(rr, rs) && same(rr, rl) && rr->used == sum);
    if (rr) std::printf("x=%.17g iters=%u evals=%u == static: %s, cost == sum of attempts: %s\n", rr->x,
                        rr->used.iterations, rr->used.evaluations, same(rr, rs) && same(rr, rl) ? "yes" : "NO",
                        rr->used == sum ? "yes" : "NO");

    std::printf("[R2] same chain on x^2+1 (all fail): ");
    const result_t fr = (*chain)(f2);
    const result_t fs = static_chain(f2);
    nxx::counters sum2{};
    std::optional<est_t> best;   // expected best: minimum |f|, ties -> the later attempt
    for (const auto& n : config) {
        auto one = (*make_solver(n))(f2);
        sum2 = sum2 + one.error().used;
        if (one.error().best && (!best || !(nxx::math::abs(best->fx) < nxx::math::abs(one.error().best->fx)))) best = one.error().best;
    }
    CHECK(!fr && same(fr, fs) && fr.error().code == nxx::errc::no_sign_change && fr.error().where == nxx::algo::bisection &&
          fr.error().used == sum2 && fr.error().best && best && same_est(*fr.error().best, *best));
    if (!fr) std::printf("%s where=bisection used=%u/%u best x=%g |f|=%g enclosure=%s == static: %s\n", name(fr.error().code),
                         fr.error().used.iterations, fr.error().used.evaluations, fr.error().best->x, fr.error().best->fx,
                         fr.error().best->enclosure ? "yes" : "no", same(fr, fs) ? "yes" : "NO");

    std::printf("[R3] config {newton_fd, brent} == static chain_fd: ");
    auto chain_fd = make_chain({"newton_fd", "brent"});
    const result_t rfd = (*chain_fd)(f);
    CHECK(rfd && rfd->by == nxx::algo::newton && same(rfd, static_chain_fd(f)));
    if (rfd) std::printf("x=%.17g evals=%u equal=%s\n", rfd->x, rfd->used.evaluations, same(rfd, static_chain_fd(f)) ? "yes" : "NO");

    std::printf("[R4] unknown name and empty range: ");
    auto bad = make_chain({"newton", "halley"});
    const result_t re = nxx::first_of(std::vector<solver_t>{})(f);
    CHECK(!bad && !re && re.error().code == nxx::errc::invalid_input && re.error().where == nxx::algo::none &&
          re.error().used == nxx::counters{} && !re.error().best);
    std::printf("%s; empty -> %s used=%u/%u best=%s\n", bad ? "?" : bad.error().c_str(), name(re.error().code),
                re.error().used.iterations, re.error().used.evaluations, re.error().best ? "yes" : "none");

    std::printf("[R5] lazy: alternatives after a success never run: ");
    int probe_calls = 0;
    const solver_t probe = [&probe_calls](const fn_t& g) { ++probe_calls; return r::brent{}(g, nxx::bracket{0.0, 2.0}); };
    const result_t lz = nxx::first_of(std::vector<solver_t>{*make_solver("bisection"), probe})(f);
    const int calls_after_success = probe_calls;
    CHECK(lz && calls_after_success == 0);
    const result_t lz2 = nxx::first_of(std::vector<solver_t>{*make_solver("newton"), probe})(f);
    CHECK(lz2 && lz2->by == nxx::algo::brent && probe_calls == 1);
    std::printf("probe calls after a success=%d, after a failure=%d\n", calls_after_success, probe_calls);

    std::printf("[R6] value semantics: ");
    solver_t a = *make_solver("newton");
    const solver_t b = *make_solver("bisection");
    a = b;                                    // copy-assignment
    CHECK(a(f) && a(f)->by == nxx::algo::bisection);
    solver_t c = std::move(a);                // a move is a copy: a is still usable
    CHECK(a(f) && c(f) && same(a(f), c(f)));
    std::vector<solver_t> v{*make_solver("secant"), static_chain, *chain};   // raw, static chain, runtime chain
    std::swap(v[0], v[2]);
    const solver_t nested = nxx::first_of(std::span<const solver_t>(v));
    const auto mixed = nxx::first_of(r::newton{}.with_derivative(df).on(0.0), nested);   // static over run-time
    const result_t rm = mixed(f);
    const nxx::counters newton_cost = (*make_solver("newton"))(f).error().used;
    CHECK(same(nested(f), rr) && rm && same_est(*rm, *rr) && rm->used == rr->used + newton_cost);
    std::printf("assign/move/swap ok; nested run-time chain == R1: %s; static first_of over it adds the failed newton's cost: %s\n",
                same(nested(f), rr) ? "yes" : "NO", rm && rm->used == rr->used + newton_cost ? "yes" : "NO");

    std::printf("[R7] fallible callback keeps the user's error type: ");
    enum class callback_error { diverged };
    using gfn_t = std::function<std::expected<double, callback_error>(double)>;
    using gsolver_t = nxx::any_solver<gfn_t, est_t, callback_error>;
    const gfn_t g = [](double x) -> std::expected<double, callback_error> {
        if (x < 0) return std::unexpected(callback_error::diverged);
        return x * x - 2.0;
    };
    const std::vector<gsolver_t> galts{r::bisection{}.on(nxx::bracket{-1.0, 2.0}), r::brent{}.on(nxx::bracket{0.0, 2.0})};
    const auto gr = nxx::first_of(galts)(g);
    const auto gs = nxx::first_of(r::bisection{}.on(nxx::bracket{-1.0, 2.0}), r::brent{}.on(nxx::bracket{0.0, 2.0}))(g);
    const auto gf = nxx::first_of(std::span(galts).first(1))(g);
    CHECK(gr && gr->by == nxx::algo::brent && same(gr, gs) && !gf && gf.error().code == nxx::errc::callback_failed &&
          gf.error().cause && *gf.error().cause == callback_error::diverged);
    std::printf("success by brent == static: %s; all-fail cause kept: %s\n", same(gr, gs) ? "yes" : "NO",
                !gf && gf.error().cause ? "yes" : "NO");

    std::printf("[R8] allocations: ");
    std::size_t n0 = g_allocs;
    const solver_t small = r::bisection{}.on(nxx::bracket{0.0, 2.0});
    const std::size_t wrap_small = g_allocs - n0;
    n0 = g_allocs;
    const solver_t big = static_chain;
    const std::size_t wrap_chain = g_allocs - n0;
    n0 = g_allocs;
    const solver_t small2 = small, big2 = big;
    const std::size_t copies = g_allocs - n0;
    n0 = g_allocs;
    const result_t call1 = (*chain)(f);
    const std::size_t per_call = g_allocs - n0;
    n0 = g_allocs;
    const result_t call2 = (*chain)([](double x) { return x * x - 2.0; });   // converts to fn_t per call
    const std::size_t per_call_lambda = g_allocs - n0;
    n0 = g_allocs;
    const std::array<double, 8> pad{};
    const result_t call3 = (*chain)([pad](double x) { return x * x - 2.0 + pad[0]; });   // 64-byte capture
    const std::size_t per_call_fat = g_allocs - n0;
    CHECK(per_call == 0 && same(call1, rr) && same(call2, rr) && same(call3, rr) && same(small2(f), small(f)) && same(big2(f), rs));
    std::printf("wrap bisection=%zu, wrap static chain=%zu, 2 copies=%zu, call(fn_t)=%zu, call(empty lambda)=%zu, call(64B lambda)=%zu\n",
                wrap_small, wrap_chain, copies, per_call, per_call_lambda, per_call_fat);

    std::printf("[R9] sizes (%s, %zu-bit): any_solver=%zu std::function<double(double)>=%zu bisection.on=%zu newton.on=%zu static chain=%zu result=%zu\n",
                stdlib(), sizeof(void*) * 8, sizeof(solver_t), sizeof(fn_t), sizeof(bisect_t),
                sizeof(decltype(r::newton{}.with_derivative(df).on(0.0))), sizeof(static_chain), sizeof(result_t));
    std::printf("%s (%d check failures)\n", failures ? "FAIL" : "PASS", failures);
    return failures ? 1 : 0;
}
