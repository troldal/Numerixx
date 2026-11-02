// TanhSinh of nxx::quadrature (Phase 3b, docs/quadrature/PHASE3B_API.md, with the numerics of PHASE3_API.md): exp()
// and the node tables against mpmath, the rule's accuracy, its steps and evaluation counts, termination and values
// that are not finite, constant evaluation, the substitution for a singularity at an end other than 0, and the
// fallback to the default. The edge cases it shares with Romberg and GaussKronrod<15> are in
// test_quadrature_kronrod.cpp's template, and the protocol's properties in test_quadrature_protocol.cpp.
//
// exp() is checked against mpmath at 60 digits within the derived 1.48 ulp (PHASE3_API.md B: the final rounding of
// 1 + p, 0.5 ulp; Horner's scheme, 0.73; the rounding of r, 0.25); the test allows 1.5. The tables' complements and
// weights are checked per node within the a-priori bound of B. The references are in
// quadrature_tanh_sinh_reference.hpp, from the script in docs/quadrature/design.md.
// Integrals are checked against closed forms or mpmath with the honesty check |value - ref| <= abs_error and within
// the tolerance. Evaluation counts are asserted exactly in double, in total and per level through the observer, for fn
// of + - * / and sqrt, against the prototype of the specification (docs/quadrature/tanh-sinh-prototype/ts.py, Python
// floats, the C++ double table loaded); no count is asserted where a libm function's last bits could move a decision.

#include "quadrature_check.hpp"
#include "quadrature_tanh_sinh_reference.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <span>
#include <utility>
#include <vector>

using namespace quadrature_check;

namespace {

    namespace ref = tanh_sinh_reference;

    // x = k / 2^s, exact in every type tested (the script keeps k within the type's digits).
    template<typename T>
    T dyadic(long long k, int s)
    { return T(T(k) / T(1LL << s)); }

    // |v / radix^e - m| for a reference value m radix^e, 1 <= m < radix, measured in bin100: in units of eps it is the
    // error in ulps of the reference's binade. A radix-2 value is split exactly by frexp first, since
    // cpp_bin_float_50's smallest values lie below cpp_bin_float_100's range.
    template<typename T>
    double scaled_error(T const &v, char const *mantissa, long long exponent)
    {
        bin100 const m(mantissa);
        using std::frexp, std::ldexp;
        if constexpr (std::numeric_limits<T>::radix == 2) {
            int e = 0;
            T const f = frexp(v, &e); // v = f 2^e, 1/2 <= f < 1; std::frexp or Boost's, found by lookup
            return static_cast<double>(bin100(abs(bin100(ldexp(wide(f), static_cast<int>(e - exponent))) - m)));
        }
        else {
            bin100 const scale = pow(bin100(10), static_cast<int>(exponent));
            return static_cast<double>(bin100(abs(bin100(wide(v) / scale) - m)));
        }
    }

    // exp() at the script's points of its domain, in ulps, within the derived 1.48 (the test allows 1.5); 0 at and
    // below x_min, and normal just above it.
    template<typename T>
    double check_exp(std::span<ref::ExpPoint const> points)
    {
        qd::detail::ExpConstants<T> const constants = qd::detail::exp_constants<T>();
        double const eps = eps_of<T>();
        double worst = 0.0;
        for (ref::ExpPoint const &point : points) {
            T const x = dyadic<T>(point.k, point.s);
            double const ulps = scaled_error(qd::detail::exp(x, constants), point.mantissa, point.exponent) / eps;
            INFO("x = " << point.k << " / 2^" << point.s);
            CHECK(ulps <= 1.5);
            worst = std::max(worst, ulps);
        }
        CHECK(qd::detail::exp(constants.x_min, constants) == T(0));
        CHECK(qd::detail::exp(T(constants.x_min - T(1)), constants) == T(0));
        T const inside = T(constants.x_min * T(T(1) - T(T(64) * T(std::numeric_limits<T>::epsilon()))));
        CHECK(qd::detail::exp(inside, constants) >= T((std::numeric_limits<T>::min)()));
        return worst;
    }

    constexpr qd::detail::ExpConstants<double> exp_double = qd::detail::exp_constants<double>();

    double node_t(int level, int i)
    { return level == 0 ? i : (2.0 * i + 1.0) / static_cast<double>(1 << level); }

    constexpr auto level_size = [](auto const &level) { return level.size; };

    // The complement and weight of node i of a level.
    template<typename T>
    std::pair<T, T> node_of(int level, int i)
    {
        return qd::detail::with_tanh_sinh_level<T>(level, [i](qd::detail::TanhSinhLevel<T> const &nodes) {
            REQUIRE(i < nodes.size);
            return std::pair<T, T>{ nodes.complement(i), nodes.weight(i) };
        });
    }

    // Every sampled node of a table within its bound; returns the largest ratio of error to bound.
    template<typename T>
    double check_table(std::span<ref::Node const> nodes, int series_additions)
    {
        double const eps = eps_of<T>();
        double worst = 0.0;
        for (ref::Node const &node : nodes) {
            INFO("level " << node.level << " node " << node.i);
            auto const [comp, w] = node_of<T>(node.level, node.i);
            NodeBound const bound = node_bound(node_t(node.level, node.i), series_additions);
            double const comp_error =
                scaled_error(comp, node.comp, node.comp_exponent) / static_cast<double>(bin100(node.comp)) / eps;
            double const w_error = scaled_error(w, node.w, node.w_exponent) / static_cast<double>(bin100(node.w)) / eps;
            CHECK(comp_error <= bound.comp);
            CHECK(w_error <= bound.w);
            worst = std::max({ worst, comp_error / bound.comp, w_error / bound.w });
        }
        return worst;
    }

    // A built-in floating type B inside a class: a type without a compile-time table, so its pieces are generated
    // (cached at run time, temporaries in constant evaluation) with B's arithmetic, and must equal B's compile-time
    // tables bit for bit.
    template<typename B>
    struct Boxed
    {
        B v;
        constexpr Boxed(int i) : v(static_cast<B>(i)) {}
        constexpr Boxed(B x) : v(x) {}
        friend constexpr Boxed operator+(Boxed x, Boxed y) { return Boxed(x.v + y.v); }
        friend constexpr Boxed operator-(Boxed x, Boxed y) { return Boxed(x.v - y.v); }
        friend constexpr Boxed operator*(Boxed x, Boxed y) { return Boxed(x.v * y.v); }
        friend constexpr Boxed operator/(Boxed x, Boxed y) { return Boxed(x.v / y.v); }
        friend constexpr Boxed operator-(Boxed x) { return Boxed(-x.v); }
        friend constexpr bool operator==(Boxed, Boxed) = default;
        friend constexpr auto operator<=>(Boxed x, Boxed y) { return x.v <=> y.v; }
    };

} // namespace

template<typename B>
struct std::numeric_limits<Boxed<B>> : std::numeric_limits<B>
{
};

namespace {

    // The compile-time pieces of a tabled type equal those generated at run time, bit for bit.
    template<typename T, int L, int... P>
    bool same_level(std::integer_sequence<int, P...>)
    {
        static_assert(nxx::real<Boxed<T>> && !qd::detail::tabled<Boxed<T>>);
        auto const same = []<int Q>(std::integral_constant<int, Q>) {
            auto const &table = qd::detail::tanh_sinh_table<T, L, Q>;
            auto const &generated = qd::detail::cached<&qd::detail::tanh_sinh_piece<Boxed<T>, L, Q>>();
            if (table.size != generated.size) return false;
            for (std::size_t i = 0; i < table.size; ++i)
                if (!(table.comp[i] == generated.comp[i].v) || !(table.w[i] == generated.w[i].v)) return false;
            return true;
        };
        return (same(std::integral_constant<int, P>{}) && ...);
    }

    template<typename T, int... L>
    bool same_tables(std::integer_sequence<int, L...>)
    { return (same_level<T, L>(std::make_integer_sequence<int, qd::detail::level_pieces<T>(L)>{}) && ...); }

} // namespace

// exp() at both ends of its domain in double, within its 1.5 ulp of mpmath's value plus the half ulp of rounding that
// to a double: e^-708 = 3.307553003638408e-308 in [2^-1022, 2^-1021), an ulp 2^-1074, just above min(); e^6.5 =
// 665.1416330443618 in [2^9, 2^10), an ulp 2^-43 (the double table needs t up to 6.56); and 0 at x_min.
static_assert(nxx::math::abs(qd::detail::exp(-708.0, exp_double) - 3.307553003638408e-308) <= 2.0 * 0x1p-1074);
static_assert(nxx::math::abs(qd::detail::exp(6.5, exp_double) - 665.1416330443618) <= 2.0 * 0x1p-43);
static_assert(qd::detail::exp(exp_double.x_min, exp_double) == 0.0);

TEST_CASE("quadrature: tanh_sinh: exp against mpmath")
{
    double const worst_float = check_exp<float>(ref::exp_float);
    double const worst_double = check_exp<double>(ref::exp_double);
    double worst_long = 0.0;
    if constexpr (std::numeric_limits<long double>::min_exponent == -16381)
        worst_long = check_exp<long double>(ref::exp_long_double); // x87 and binary128 share the domain
    else worst_long = check_exp<long double>(ref::exp_double);
    double const worst_bin50 = check_exp<bin50>(ref::exp_bin50);
    double const worst_bin50_et = check_exp<bin50_et>(ref::exp_bin50);
    double const worst_dec50 = check_exp<dec50>(ref::exp_dec50);
    MESSAGE("exp, largest error in ulps: float " << worst_float << ", double " << worst_double << ", long double "
                                                 << worst_long << ", cpp_bin_float_50 " << worst_bin50 << " ("
                                                 << worst_bin50_et << " with expression templates), cpp_dec_float_50 "
                                                 << worst_dec50);
    // The compile-time constants equal the run-time ones (the static_asserts above use them).
    qd::detail::ExpConstants<double> const runtime = qd::detail::exp_constants<double>();
    CHECK(runtime.c1 == exp_double.c1);
    CHECK(runtime.c2 == exp_double.c2);
    CHECK(runtime.x_min == exp_double.x_min);
}

// The double table in constant evaluation: the nodes per level (the script's counts), and level 9 across its two
// pieces (1,024 and 656 nodes), its first node t = 2^-9 and its last t = 3359 / 512 = 6.56, both normal.
static_assert(qd::detail::with_tanh_sinh_level<double>(0, level_size) == 7);
static_assert(qd::detail::with_tanh_sinh_level<double>(4, level_size) == 53);
static_assert(qd::detail::with_tanh_sinh_level<double>(9, level_size) == 1680);
static_assert(qd::detail::with_tanh_sinh_level<double>(9, [](auto const &level) {
    return level.complement(1679) >= std::numeric_limits<double>::min() && level.weight(0) < 1.0;
}));
// A type without a table in constant evaluation: its level 4, generated as a temporary, equals the double table's.
static_assert(qd::detail::with_tanh_sinh_level<Boxed<double>>(4, [](auto const &level) {
    return qd::detail::with_tanh_sinh_level<double>(4, [&level](auto const &table) {
        bool same = level.size == table.size;
        for (int i = 0; same && i < table.size; ++i)
            same = level.complement(i).v == table.complement(i) && level.weight(i).v == table.weight(i);
        return same;
    });
}));

TEST_CASE("quadrature: tanh_sinh: the tables against mpmath")
{
    // The sampled nodes within the bound of B: double (8 series additions below t = 1/2) and cpp_bin_float_50 (18).
    double const worst_double = check_table<double>(ref::table_double, 8);
    double const worst_bin50 = check_table<bin50>(ref::table_bin50, 18);
    MESSAGE("tables, largest ratio of error to bound: double " << worst_double << ", cpp_bin_float_50 "
                                                               << worst_bin50);
    // The nodes per level: those with q = e^(-2 sinh t) >= 2^-m, 2^-m the smallest power of 2 at or above min(), that
    // is 2 sinh t < m ln 2 (the script's counts; no node lies within 4e-5 of that boundary, relative, against x_min's
    // 4 eps). In total float 2,289, double 3,361, long double (80-bit and binary128) 4,781, cpp_bin_float_50 10,814,
    // cpp_dec_float_50 9,655.
    auto const counts = []<typename T>(T, std::span<int const> expected) {
        for (int level = 0; level < qd::detail::tanh_sinh_levels; ++level) {
            INFO("level " << level);
            int const n = qd::detail::with_tanh_sinh_level<T>(level, level_size);
            CHECK(n == expected[static_cast<std::size_t>(level)]);
            CHECK(n <= qd::detail::level_capacity<T>(level)); // the capacity bound holds
        }
    };
    counts(0.0f, ref::nodes_float);
    counts(0.0, ref::nodes_double);
    if constexpr (std::numeric_limits<long double>::min_exponent == -16381) counts(0.0L, ref::nodes_long_double);
    else counts(0.0L, ref::nodes_double);
    counts(bin50{}, ref::nodes_bin50);
    counts(bin50_et{}, ref::nodes_bin50);
    counts(dec50{}, ref::nodes_dec50);
    // Compile-time and run-time tables are bit-identical.
    CHECK(same_tables<float>(std::make_integer_sequence<int, qd::detail::tanh_sinh_levels>{}));
    CHECK(same_tables<double>(std::make_integer_sequence<int, qd::detail::tanh_sinh_levels>{}));
    CHECK(same_tables<long double>(std::make_integer_sequence<int, qd::detail::tanh_sinh_levels>{}));
}

// The method (step 3).
namespace {

    using tanh_sinh = qd::TanhSinh;

    // A run with its evaluations counted.
    template<typename T, typename F>
    std::pair<std::expected<qd::Integral<T>, qd::Failure<T>>, int> run_counted(F fn, T a, T b, qd::Criteria<T> c = {})
    {
        int n = 0;
        auto const r = qd::integrate<tanh_sinh>(Counted<F>{ fn, &n }, a, b, c);
        return { r, n };
    }

    // Each step an observer receives, with the evaluations counted so far.
    struct Recorder
    {
        int const *evaluations;
        std::vector<qd::Step<double>> steps{};
        std::vector<int> counts{};

        void operator()(qd::Step<double> const &s)
        {
            steps.push_back(s);
            counts.push_back(*evaluations);
        }
    };

    // An observed run in double: its result and what the observer saw.
    using Observed = std::pair<std::expected<qd::Integral<double>, qd::Failure<double>>, Recorder>;

    template<typename F>
    Observed run_observed(F fn, double a, double b, qd::Criteria<double> c = {})
    {
        int n = 0;
        Recorder recorder{ &n };
        auto const r = qd::integrate<tanh_sinh>(qd::observer(recorder), Counted<F>{ fn, &n }, a, b, c);
        return { r, recorder };
    }

    double root_of(double x)
    { return std::sqrt(x); }

    // e^x by its Taylor polynomial of degree 20 in Horner's form: + - * / only (its error on [0, 1] is below 1e-19, so
    // e - 1 remains the reference), for exact counts.
    double taylor_exp(double x)
    {
        double s = 1.0;
        for (int n = 20; n > 0; --n) s = 1.0 + x * s / n;
        return s;
    }

    // x^-0.875 = 1 / (x^(1/2) x^(1/4) x^(1/8)) and x^-(127/128) by 7 nested square roots: powers by sqrt alone.
    double power_875(double x)
    { return 1.0 / (root_of(x) * root_of(root_of(x)) * root_of(root_of(root_of(x)))); }

    double power_127(double x)
    {
        double r = x;
        double product = 1.0;
        for (int k = 0; k < 7; ++k) {
            r = root_of(r);
            product *= r;
        }
        return 1.0 / product;
    }

    // A hash of x's bits in [-1, 1]: noise that never resolves (as test_quadrature_kronrod.cpp's).
    double noise_of(double x)
    {
        std::uint64_t h = std::bit_cast<std::uint64_t>(x) * 0x9E3779B97F4A7C15ull;
        h ^= h >> 29;
        h *= 0xBF58476D1CE4E5B9ull;
        h ^= h >> 32;
        return static_cast<double>(h >> 11) / 4503599627370496.0 - 1.0;
    }

    enum class status { ok, not_converged, non_finite };

    template<typename T>
    status status_of(std::expected<qd::Integral<T>, qd::Failure<T>> const &r)
    {
        if (r) return status::ok;
        return r.error().reason == qd::error::not_converged ? status::not_converged : status::non_finite;
    }

    // The fallback of integrate()'s comment (PHASE3B_API.md Q4): the default only where TanhSinh is not_converged with
    // a best, and its success only where it agrees with that best within the two estimates.
    template<typename G>
    std::expected<qd::Integral<double>, qd::Failure<double>> with_fallback(G const &g, double a, double b,
                                                                           qd::Criteria<double> c = {})
    {
        auto const t = qd::integrate<qd::TanhSinh>(g, a, b, c);
        return t.or_else([&](qd::Failure<double> const &f) -> decltype(t) {
            if (f.reason != qd::error::not_converged || !f.best) return t;
            auto const k = qd::integrate(g, a, b, c);
            if (k && std::abs(k->value - f.best->value) <= f.best->abs_error + k->abs_error) return k;
            return t;
        });
    }

    // pi to 60 digits (mpmath).
    constexpr char const *pi = "3.14159265358979323846264338327950288419716939937510582097494";

    constexpr auto ts_square = [](double x) { return x * x; };
    constexpr auto ts_runge = [](double x) { return 1.0 / (1.0 + x * x); };
    constexpr auto ts_reciprocal = [](double x) { return 1.0 / x; };

} // namespace

static_assert(qd::iterative_method<tanh_sinh, double>);

TEST_CASE("quadrature: tanh_sinh: accuracy on smooth and endpoint cases")
{
    // Kahaner's 1, 3, 6, 7 and 19 (e^x, sqrt x, x^1.5, 1 / sqrt x, log x), x^-0.9 and sqrt(x) log x on [0, 1]: ok, the
    // honesty check, and the error within the request rel_tol |ref|. References are closed forms.
    struct Case
    {
        char const *name;
        int fn;
        char const *ref;
    };
    Case const cases[] = {
        { "e^x", 0, e_minus_1 },
        { "sqrt x", 1, "0.6666666666666666666666666666666666666666666666666666666666666666666667" },
        { "x^1.5", 2, "0.4" },
        { "1 / sqrt x", 3, "2" },
        { "log x", 4, "-1" },
        { "x^-0.9", 5, "10" },
        { "sqrt(x) log x", 6, "-0.4444444444444444444444444444444444444444444444444444444444444444444444" },
    };
    auto const check = []<typename T>(T, Case const &c, T const &rel_tol) {
        auto const fn = [&c](T const &x) -> T {
            using std::exp, std::log, std::pow, std::sqrt;
            switch (c.fn) {
            case 0: return T(exp(x));
            case 1: return T(sqrt(x));
            case 2: return T(x * sqrt(x));
            case 3: return T(T(1) / sqrt(x));
            case 4: return T(log(x));
            case 5: return T(pow(x, T(T(-9) / T(10))));
            default: return T(sqrt(x) * log(x));
            }
        };
        INFO(c.name << ", rel_tol " << static_cast<double>(wide(rel_tol)));
        auto const r = qd::integrate<tanh_sinh>(fn, T(0), T(1), { .rel_tol = rel_tol });
        REQUIRE(r.has_value());
        double const error = error_of(r->value, c.ref);
        CHECK(error <= static_cast<double>(wide(r->abs_error))); // the honesty check
        CHECK(error <= static_cast<double>(wide(rel_tol)) * std::abs(static_cast<double>(bin100(c.ref))));
    };
    for (Case const &c : cases) {
        for (double const rel : { 1.5e-8, 1e-12 }) check(0.0, c, rel);
        for (char const *rel : { "1e-25", "1e-45" }) {
            check(bin50{}, c, bin50(rel));
            check(bin50_et{}, c, bin50_et(rel));
            check(dec50{}, c, dec50(rel));
        }
    }
}

TEST_CASE("quadrature: tanh_sinh: steps")
{
    // Report i is level i (PHASE3B_API.md B and E): level 0 has no current and is not ready; ready from level 4; the
    // evaluations after each level as the prototype's.
    auto const [r, seen] = run_observed(taylor_exp, 0.0, 1.0);
    REQUIRE(r.has_value());
    CHECK(seen.counts == std::vector<int>{ 10, 18, 33, 62, 121 });
    REQUIRE(seen.steps.size() == 5);
    CHECK(!seen.steps[0].current.has_value());
    CHECK(!seen.steps[0].ready);
    for (std::size_t i = 1; i < seen.steps.size(); ++i) {
        INFO("step " << i);
        CHECK(seen.steps[i].current.has_value());
        CHECK(seen.steps[i].ready == (i >= 4));
        CHECK(!seen.steps[i].exhausted);
        CHECK(seen.steps[i].trusted);
    }
    CHECK(seen.steps.back().current == *r); // the last report is the result (a < b)
    CHECK(r == qd::integrate<tanh_sinh>(taylor_exp, 0.0, 1.0)); // observed equals unobserved

    // A run to the cap: 10 reports, only the last exhausted; 1/x's estimate is unbounded, so no report has a current.
    auto const [divergent, cap] = run_observed(ts_reciprocal, 0.0, 1.0);
    CHECK(divergent == std::unexpected(qd::Failure<double>{ qd::error::not_converged, std::nullopt }));
    CHECK(cap.counts == std::vector<int>{ 10, 21, 41, 81, 163, 325, 648, 1293, 2582, 5159 });
    REQUIRE(cap.steps.size() == 10);
    for (std::size_t i = 0; i < cap.steps.size(); ++i) {
        INFO("step " << i);
        CHECK(cap.steps[i].exhausted == (i == 9));
        CHECK(!cap.steps[i].current.has_value());
    }

    // The other levels the prototype records: 1 / sqrt x, a rounding cut on [1, 1 + 1e-9], the table's end on
    // x^-(127/128) and the step at 0.3.
    CHECK(run_observed([](double x) { return 1.0 / root_of(x); }, 0.0, 1.0).second.counts
          == std::vector<int>{ 10, 19, 35, 67, 131 });
    CHECK(run_observed([](double) { return 1.0; }, 1.0, 1.0 + 1e-9).second.counts
          == std::vector<int>{ 5, 11, 23, 45, 89 });
    CHECK(run_observed(power_127, 0.0, 1.0, { .rel_tol = 1e-4 }).second.counts
          == std::vector<int>{ 10, 21, 41, 81, 163, 325, 648, 1293, 2582, 5159 });
    CHECK(run_observed([](double x) { return x >= 0.3 ? 1.0 : 0.0; }, 0.0, 1.0).second.counts
          == std::vector<int>{ 10, 15, 23, 39, 72, 137, 267, 525, 1039, 2066 });

    // No point strictly inside: one report, exhausted, with no current, and fn never called. A failure at level 0
    // gives no report.
    auto const [none, single] = run_observed(ts_square, 1.0, std::nextafter(1.0, 2.0));
    CHECK(none == std::unexpected(qd::Failure<double>{ qd::error::not_converged, std::nullopt }));
    REQUIRE(single.steps.size() == 1);
    CHECK(single.steps[0].exhausted);
    CHECK(!single.steps[0].current.has_value());
    CHECK(single.counts == std::vector<int>{ 0 });
    auto const [shifted, silent] = run_observed([](double x) { return x >= 0.5 ? root_of(x - 0.5) : dnan; }, 0.0, 1.0);
    CHECK(status_of(shifted) == status::non_finite);
    CHECK(silent.steps.empty());
}

TEST_CASE("quadrature: tanh_sinh: evaluation counts in double")
{
    // The prototype's counts with the C++ table (docs/quadrature/design.md), with the status and whether a failure
    // carries a best estimate: absent only where the estimate is unbounded or nothing can be sampled.
    double const nan = dnan;
    struct Case
    {
        char const *name;
        std::pair<std::expected<qd::Integral<double>, qd::Failure<double>>, int> run;
        status expected;
        int evaluations;
        bool best;
    };
    Case const cases[] = {
        // The minimum cost: level 4, 121 evaluations, at every tolerance down to rel 0 (the 2F test).
        { "e^x by Taylor [0,1]", run_counted(taylor_exp, 0.0, 1.0), status::ok, 121, false },
        { "e^x by Taylor [0,1] rel 0 (the 2F test)", run_counted(taylor_exp, 0.0, 1.0, { .rel_tol = 0.0 }), status::ok,
          121, false },
        { "1/(1+x^2) [0,1]", run_counted(+ts_runge, 0.0, 1.0), status::ok, 121, false },
        { "x^2 [0,3]", run_counted(+ts_square, 0.0, 3.0), status::ok, 105, false },
        { "1/sqrt x [0,1]", run_counted([](double x) { return 1.0 / root_of(x); }, 0.0, 1.0), status::ok, 131, false },
        { "x^-0.875 [0,1]", run_counted(power_875, 0.0, 1.0), status::ok, 151, false },
        { "x^-0.875 [0,1] 1e-12", run_counted(power_875, 0.0, 1.0, { .rel_tol = 1e-12 }), status::ok, 151, false },
        { "a peak of width 1e-2 at 0.42",
          run_counted([](double x) { return 1.0 / (1e-4 + (x - 0.42) * (x - 0.42)); }, 0.0, 1.0), status::ok, 3602,
          false },
        // The cap, 10 levels, with no stagnation rule: unbounded estimates, no best.
        { "1/x [0,1] (the cap, the table's end)", run_counted(+ts_reciprocal, 0.0, 1.0), status::not_converged, 5159,
          false },
        { "1/x^2 [0,1] (the cap, a value not finite)", run_counted([](double x) { return 1.0 / (x * x); }, 0.0, 1.0),
          status::not_converged, 3022, false },
        // The table's end with the geometric end term (D1) refuses x^-(127/128) at 1e-4, with an honest best (below).
        { "x^-(127/128) [0,1] 1e-4", run_counted(power_127, 0.0, 1.0, { .rel_tol = 1e-4 }), status::not_converged,
          5159, true },
        // The rounding cut and its sliver (D3).
        { "1 on [1, 1 + 1e-9]", run_counted([](double) { return 1.0; }, 1.0, 1.0 + 1e-9), status::ok, 89, false },
        { "1 on [1000, 1001] rel 0", run_counted([](double) { return 1.0; }, 1000.0, 1001.0, { .rel_tol = 0.0 }),
          status::ok, 109, false },
        { "1e300 on [0, 3 denorm_min]", run_counted([](double) { return 1e300; }, 0.0, 3.0 * dmin), status::ok, 27,
          false },
        { "(1 - x)^-0.875 [0,1] (the power-law end term)",
          run_counted([](double x) { return power_875(1.0 - x); }, 0.0, 1.0), status::not_converged, 3690, true },
        // Values that are not finite (PHASE3_API.md D).
        { "sqrt(x - 0.5) [0,1]", run_counted([nan](double x) { return x >= 0.5 ? root_of(x - 0.5) : nan; }, 0.0, 1.0),
          status::non_finite, 5, false },
        { "sqrt(x - 0.01) [0,1]",
          run_counted([nan](double x) { return x >= 0.01 ? root_of(x - 0.01) : nan; }, 0.0, 1.0), status::non_finite,
          6, false },
        { "NaN on |x - 0.3| < 0.05",
          run_counted([nan](double x) { return std::abs(x - 0.3) < 0.05 ? nan : 1.0; }, 0.0, 1.0), status::non_finite,
          15, false },
        // A NaN region of width 0.02 between the nodes: no sampling sees it.
        { "NaN on |x - 0.3| < 0.01 (missed)",
          run_counted([nan](double x) { return std::abs(x - 0.3) < 0.01 ? nan : 1.0; }, 0.0, 1.0), status::ok, 121,
          false },
        { "1/sqrt x, NaN below 1e-300 (a tail)",
          run_counted([nan](double x) { return x < 1e-300 ? nan : 1.0 / root_of(x); }, 0.0, 1.0), status::ok, 131,
          false },
        // Noise that never resolves, at rel 0: the cap, with an honest best.
        { "e^x (1 + 1e-10 h(x)) rel 0",
          run_counted([](double x) { return taylor_exp(x) * (1.0 + 1e-10 * noise_of(x)); }, 0.0, 1.0,
                      { .rel_tol = 0.0 }),
          status::not_converged, 3711, true },
        // Interior features: honest failures (use the default).
        { "the step at 0.3 (Kahaner 2)", run_counted([](double x) { return x >= 0.3 ? 1.0 : 0.0; }, 0.0, 1.0),
          status::not_converged, 2066, true },
        { "|x - 0.37|", run_counted([](double x) { return std::abs(x - 0.37); }, 0.0, 1.0), status::not_converged, 3718,
          true },
        // No point of double strictly inside.
        { "1 on [1, nextafter(1)]", run_counted([](double) { return 1.0; }, 1.0, std::nextafter(1.0, 2.0)),
          status::not_converged, 0, false },
        // A term that overflows ends the run at once: the centre's, (10 / 2) 1e308.
        { "1e308 on [0, 10]", run_counted([](double) { return 1e308; }, 0.0, 10.0), status::non_finite, 1, false },
    };
    for (Case const &c : cases) {
        INFO(c.name);
        auto const &r = c.run.first;
        CHECK(status_of(r) == c.expected);
        CHECK(c.run.second == c.evaluations);
        if (!r) CHECK(r.error().best.has_value() == c.best);
    }
    // The bests are honest: x^-(127/128) = 128, (1 - x)^-0.875 = 8, the noise within its amplitude 1e-10 (e - 1), the
    // step 0.7 and the kink (0.37^2 + 0.63^2) / 2 = 0.2669.
    auto const honest = [](Case const &c, char const *value, double margin) {
        INFO(c.name);
        REQUIRE(c.run.first.error().best.has_value());
        qd::Integral<double> const &best = *c.run.first.error().best;
        CHECK(error_of(best.value, value) <= best.abs_error + margin);
    };
    honest(cases[10], "128", 0.0);
    honest(cases[14], "8", 0.0);
    honest(cases[20], e_minus_1, 1e-10 * 1.7183);
    honest(cases[21], "0.7", 0.0);
    honest(cases[22], "0.2669", 0.0);
    // The rounding cut's results (D3: the raw sum and the slivers kept apart), against 1e-9 (1 + 1e-9 - 1 as doubles)
    // and 1, within the rule's rounding for a constant: the table's weights err by B_w eps each (the a-priori bound),
    // so the value by (b - a) / 2 eps sum h w B_w; the products, the compensated sums and the slivers add 6 eps
    // (b - a); the trapezoidal sum of the weights over the line is 2 up to e^(-pi^2 / h), negligible from h = 1/16.
    // Draft 1's reading, the slivers in the running sum, erred by 6.8e-8 and 5.3e-14 relative.
    double const rounding = deps * (tanh_sinh_sums().w_bound / 2.0 + 6.0);
    auto const &narrow = cases[11].run.first;
    double const width = (1.0 + 1e-9) - 1.0;
    CHECK(error_of(narrow->value, bin100(width)) <= rounding * width);
    CHECK(error_of(cases[12].run.first->value, "1") <= rounding);
    // The first ends ok with an estimate above the request: twice its floor, which holds the two slivers, h W at the
    // last kept value times the half-width per side; at level 4 the last kept node is t = 2.75 (its point 1.7e-16 from
    // the end, the next rounding onto it) and h W = 2.0e-7 beyond it, so each sliver is 1.0e-16 and the floor 2.0e-16,
    // 2.0e-7 of the width (0.9 half-ulps of 1 a side; up to 1.6 on other widths next to 1).
    CHECK(narrow->abs_error > 1.5e-8 * width);
    CHECK(narrow->abs_error <= 2.5e-7 * width);
    CHECK(error_of(cases[13].run.first->value, bin100(3.0 * dmin) * bin100(1e300)) <= cases[13].run.first->abs_error);
}

TEST_CASE("quadrature: tanh_sinh: termination")
{
    // x^-0.99 by pow at 1e-4: a table's end (fn is finite down to 4.3e-312, below the table), whose geometric end term
    // refuses it with a best (PHASE3B_API.md D1); an end term of 0 there, or h Te, returned ok with an error of 0.08.
    auto const strong = qd::integrate<tanh_sinh>([](double x) { return std::pow(x, -0.99); }, 0.0, 1.0,
                                                 { .rel_tol = 1e-4 });
    REQUIRE(!strong.has_value());
    CHECK(strong.error().reason == qd::error::not_converged);
    REQUIRE(strong.error().best.has_value());
    CHECK(error_of(strong.error().best->value, "100") <= strong.error().best->abs_error);
    // x^-0.99 on [0, 1e-40]: the side at 0 ends at the table, then by rounding, then at a tail; the best is honest.
    // The integral is 100 (1e-40)^0.01 = 100 10^-0.4.
    auto const tiny = qd::integrate<tanh_sinh>([](double x) { return std::pow(x, -0.99); }, 0.0, 1e-40,
                                               { .rel_tol = 1e-4 });
    REQUIRE(!tiny.has_value());
    REQUIRE(tiny.error().best.has_value());
    CHECK(error_of(tiny.error().best->value, "39.8107170553497250770252305087752043487677037297380446865284")
          <= tiny.error().best->abs_error);
    // (1 - x)^-0.9 by pow: limited by the spacing of double near 1, not_converged with an honest best.
    auto const far_end = qd::integrate<tanh_sinh>([](double x) { return std::pow(1.0 - x, -0.9); }, 0.0, 1.0);
    REQUIRE(!far_end.has_value());
    CHECK(far_end.error().reason == qd::error::not_converged);
    REQUIRE(far_end.error().best.has_value());
    CHECK(error_of(far_end.error().best->value, "10") <= far_end.error().best->abs_error);
    // Level 0 drops each side's run of outermost small terms, from the outside in (D2): with fn exactly 0 at the left
    // point of t = 2 and large beyond it, the exact value 2/3 - 2 x2; cutting at that first small term instead missed
    // 2.5e-5.
    double const x2 = 0.5 * qd::detail::with_tanh_sinh_level<double>(0, [](auto const &nodes) {
        return nodes.complement(2);
    });
    auto const zero_inside = qd::integrate<tanh_sinh>([x2](double x) { return (x - x2) / root_of(x); }, 0.0, 1.0);
    REQUIRE(zero_inside.has_value());
    CHECK(error_of(zero_inside->value, bin100(2) / 3 - 2 * bin100(x2)) <= zero_inside->abs_error);
    // 1/x in cpp_bin_float_50: the cap, no best, as in double.
    auto const divergent = qd::integrate<tanh_sinh>([](bin50 const &x) { return bin50(1 / x); }, bin50(0), bin50(1));
    CHECK(divergent == std::unexpected(qd::Failure<bin50>{ qd::error::not_converged, std::nullopt }));
}

TEST_CASE("quadrature: tanh_sinh: values not finite")
{
    // In the tail (comp < eps beyond every kept point) a value that is not finite cuts the side, as for 1/x^2 at the
    // points nearest 0; elsewhere it is non_finite, with no best, whatever the type.
    qd::Failure<double> const invalid{ qd::error::non_finite, std::nullopt };
    CHECK(qd::integrate<tanh_sinh>([](double x) { return std::sqrt(x - 0.5); }, 0.0, 1.0) == std::unexpected(invalid));
    CHECK(qd::integrate<tanh_sinh>([](double x) { return std::log(x - 0.25); }, 0.0, 1.0) == std::unexpected(invalid));
    auto const tail = qd::integrate<tanh_sinh>([](double x) { return x < 1e-300 ? dinf : std::log(x); }, 0.0, 1.0);
    REQUIRE(tail.has_value());
    CHECK(error_of(tail->value, "-1") <= tail->abs_error);
    auto const shifted = qd::integrate<tanh_sinh>([](bin50 const &x) { return bin50(sqrt(bin50(x - bin50(0.5)))); },
                                                  bin50(0), bin50(1));
    CHECK(shifted == std::unexpected(qd::Failure<bin50>{ qd::error::non_finite, std::nullopt }));
}

// Constant evaluation: polynomial and rational cases with the honesty check and within the request (the run-time
// twins in "quadrature: tanh_sinh: constant evaluation" check that run time gives the same results, bit for bit); 1/x
// on [0, 1] runs to the cap with an unbounded estimate and ends not_converged with no best, a flag rather than inf
// arithmetic, which would not be a constant expression; a NaN fn is non_finite. x^2 on [0, 3] is 9 (105 evaluations),
// 1 / (1 + x^2) on [0, 1] pi / 4 (121), against the double nearest pi / 4 (its half ulp far inside the estimate).
namespace {
    constexpr auto ts_square_result = qd::integrate<tanh_sinh>(ts_square, 0.0, 3.0);
    constexpr auto ts_runge_result = qd::integrate<tanh_sinh>(ts_runge, 0.0, 1.0);
    constexpr auto ts_reciprocal_result = qd::integrate<tanh_sinh>(ts_reciprocal, 0.0, 1.0);
    constexpr auto ts_square_f = [](float x) { return x * x; };
    constexpr auto ts_square_l = [](long double x) { return x * x; };
    constexpr auto ts_square_float = qd::integrate<tanh_sinh>(ts_square_f, 0.0f, 3.0f);
    constexpr auto ts_square_long = qd::integrate<tanh_sinh>(ts_square_l, 0.0L, 3.0L);
} // namespace
static_assert(ts_square_result && nxx::math::abs(ts_square_result->value - 9.0) <= ts_square_result->abs_error
              && ts_square_result->abs_error <= 1.5e-8 * 9.0);
static_assert(ts_runge_result
              && nxx::math::abs(ts_runge_result->value - 0.78539816339744831) <= ts_runge_result->abs_error
              && ts_runge_result->abs_error <= 1.5e-8 * 0.7854);
static_assert(ts_reciprocal_result == std::unexpected(qd::Failure<double>{ qd::error::not_converged, std::nullopt }));
static_assert(qd::integrate<tanh_sinh>([](double) { return std::numeric_limits<double>::quiet_NaN(); }, 0.0, 1.0)
              == std::unexpected(qd::Failure<double>{ qd::error::non_finite, std::nullopt }));
static_assert(ts_square_float && nxx::math::abs(ts_square_float->value - 9.0f) <= ts_square_float->abs_error);
static_assert(ts_square_long && nxx::math::abs(ts_square_long->value - 9.0L) <= ts_square_long->abs_error);

TEST_CASE("quadrature: tanh_sinh: constant evaluation")
{
    CHECK(qd::integrate<tanh_sinh>(ts_square, 0.0, 3.0) == ts_square_result);
    CHECK(qd::integrate<tanh_sinh>(ts_runge, 0.0, 1.0) == ts_runge_result);
    CHECK(qd::integrate<tanh_sinh>(ts_reciprocal, 0.0, 1.0) == ts_reciprocal_result);
    CHECK(qd::integrate<tanh_sinh>(ts_square_f, 0.0f, 3.0f) == ts_square_float);
    CHECK(qd::integrate<tanh_sinh>(ts_square_l, 0.0L, 3.0L) == ts_square_long);
    // The constant's error against the closed forms, from bin100: the static_asserts' references rounded.
    CHECK(error_of(ts_square_result->value, "9") <= ts_square_result->abs_error);
    CHECK(error_of(ts_runge_result->value, pi_over_4) <= ts_runge_result->abs_error);
}

TEST_CASE("quadrature: tanh_sinh: the substitution for a singularity at an end other than 0")
{
    // cos(x) (1 - x)^-0.9 on [0, 1] = sum_n (-1)^n (cos 1 / ((2n)! (2n + 0.1)) + sin 1 / ((2n + 1)! (2n + 1.1)))
    // (mpmath; the same as 10 times the integral of cos(1 - v^10) on [0, 1]). Directly the singularity at 1 is limited
    // by the spacing of double there; with u = 1 - x it moves to 0, where u^-0.9 is formed from u itself.
    char const *exact = "6.00083408694045601486578515775831407028906508236";
    auto const direct =
        qd::integrate<tanh_sinh>([](double x) { return std::cos(x) * std::pow(1.0 - x, -0.9); }, 0.0, 1.0);
    REQUIRE(!direct.has_value());
    CHECK(direct.error().reason == qd::error::not_converged);
    auto const substituted =
        qd::integrate<tanh_sinh>([](double u) { return std::cos(1.0 - u) * std::pow(u, -0.9); }, 0.0, 1.0);
    REQUIRE(substituted.has_value());
    CHECK(error_of(substituted->value, exact) <= substituted->abs_error);
    CHECK(error_of(substituted->value, exact) <= 1.5e-8 * 6.0008);
    // Singular at both ends: 1 / sqrt(1 - x^2) on [-1, 1] = pi, split at 0, each half with u = 1 -+ x as the integral
    // of 1 / sqrt(u (2 - u)) over [0, 1], pi / 2, with fn written in u (1 - x^2 = u (2 - u)).
    auto const half = qd::integrate<tanh_sinh>([](double u) { return 1.0 / root_of(u * (2.0 - u)); }, 0.0, 1.0);
    REQUIRE(half.has_value());
    CHECK(error_of(2.0 * half->value, pi) <= 2.0 * half->abs_error);
    CHECK(error_of(2.0 * half->value, pi) <= 1.5e-8 * 3.1416);
}

TEST_CASE("quadrature: tanh_sinh: the fallback to the default")
{
    // integrate()'s pattern (PHASE3B_API.md Q4). A step: TanhSinh fails honestly, the default ends ok, and it agrees
    // with TanhSinh's best, so it is returned.
    auto const step = [](double x) { return x >= 0.3 ? 1.0 : 0.0; };
    auto const stepped = with_fallback(step, 0.0, 1.0);
    REQUIRE(stepped.has_value());
    CHECK(error_of(stepped->value, "0.700000000000000011102230246251565404236316681") <= stepped->abs_error);
    CHECK(stepped == qd::integrate(step, 0.0, 1.0));
    // An endpoint singularity: TanhSinh itself ends ok, and the default is not called (as many calls as TanhSinh
    // alone makes).
    auto const log_of = [](double x) { return std::log(x); };
    int calls = 0;
    int alone = 0;
    auto const logarithm = with_fallback(Counted{ log_of, &calls }, 0.0, 1.0);
    (void)qd::integrate<tanh_sinh>(Counted{ log_of, &alone }, 0.0, 1.0);
    REQUIRE(logarithm.has_value());
    CHECK(calls == alone);
    CHECK(error_of(logarithm->value, "-1") <= logarithm->abs_error);
    // x^-0.99 at 1e-4: TanhSinh fails honestly, and the default's ok (99.47 +- 0.0099, an error of 0.53) contradicts
    // TanhSinh's best (99.92 +- 0.099), so TanhSinh's failure stands; or_else alone returned that false success.
    auto const power = [](double x) { return std::pow(x, -0.99); };
    auto const strong = with_fallback(power, 0.0, 1.0, { .rel_tol = 1e-4 });
    CHECK(strong == qd::integrate<tanh_sinh>(power, 0.0, 1.0, { .rel_tol = 1e-4 }));
    auto const plain = qd::integrate(power, 0.0, 1.0, { .rel_tol = 1e-4 });
    REQUIRE(plain.has_value());
    CHECK(error_of(plain->value, "100") > plain->abs_error); // the false success the pattern rejects
    // 1/x: no best, so the default is not called.
    CHECK(with_fallback(ts_reciprocal, 0.0, 1.0)
          == std::unexpected(qd::Failure<double>{ qd::error::not_converged, std::nullopt }));
}

// The independent review's cases (Codex, PHASE3B_API.md): two valid types it found the method failing in, exponents
// that end before eps^2 and exponents that reach further below 0 than above.
namespace {
    using narrow = mp::number<mp::cpp_bin_float<100, mp::digit_base_2, void, int, -100, 100>>;
    using lopsided = mp::number<mp::cpp_bin_float<53, mp::digit_base_2, void, int, -10000, 1000>>;
} // namespace

TEST_CASE("quadrature: tanh_sinh: the independent review's cases")
{
    // The power law's products d |y| below min(), D = denorm_min: (D / (1 - x)) / 4 diverges, and both products at the
    // last two distinct points near 1 were 0, which made the excess 0 and the run ok {3 D, 50 D}. From ratios, the
    // estimate is unbounded: no best.
    auto const divergent =
        qd::integrate<tanh_sinh>([](double x) { return (dmin / (1.0 - x)) / 4.0; }, 0.0, 1.0, { .rel_tol = 0.0 });
    CHECK(divergent == std::unexpected(qd::Failure<double>{ qd::error::not_converged, std::nullopt }));
    // (D (1 - x)^(-1023/1024)) / 4 = 256 D: no longer ok, and its best honest.
    auto const tiny = qd::integrate<tanh_sinh>(
        [](double x) { return (dmin * std::pow(1.0 - x, -1023.0 / 1024.0)) / 4.0; }, 0.0, 1.0, { .rel_tol = 0.0 });
    REQUIRE(!tiny.has_value());
    REQUIRE(tiny.error().best.has_value());
    CHECK(error_of(tiny.error().best->value, bin100(256) * bin100(dmin)) <= tiny.error().best->abs_error);
    // eps^2 / 16 and eps / 8 below min() (min() 2^-100, eps 2^-99): ln 2's series and sinh's ran forever; now the
    // constants and the nodes form, and a run ends.
    static_assert(nxx::real<narrow>);
    auto const one = qd::integrate<tanh_sinh>([](narrow const &) { return narrow(1); }, narrow(0), narrow(1));
    REQUIRE(one.has_value());
    CHECK(error_of(one->value, "1") <= one->abs_error);
    // Exponents to -10000 but only to 1000: 2^1024 overflowed in the search for x_min, which came out as -1418.9 and
    // flushed e^-2000 to 0. Now x_min lies at ln min() = -6931.47 (4 eps inward), e^-2000 = 2.5765e-869 (mpmath) is
    // within exp()'s 1.5 ulp, and log x integrates.
    static_assert(nxx::real<lopsided>);
    qd::detail::ExpConstants<lopsided> const constants = qd::detail::exp_constants<lopsided>();
    CHECK(constants.x_min < lopsided(-6931));
    char const *e_2000 = "1.52617269672275131712597356791827058566759031120209"; // times 2^-2886
    CHECK(scaled_error(qd::detail::exp(lopsided(-2000), constants), e_2000, -2886) / eps_of<lopsided>() <= 1.5);
    auto const log_x = qd::integrate<tanh_sinh>([](lopsided const &x) { return lopsided(log(x)); }, lopsided(0),
                                                lopsided(1));
    REQUIRE(log_x.has_value());
    CHECK(error_of(log_x->value, "-1") <= log_x->abs_error);
}
