// TanhSinh of nxx::quadrature (Phase 3b, docs/quadrature/PHASE3B_API.md, with the numerics of PHASE3_API.md). Step 1:
// detail::exp, the exponential in T from + - * / on the domain the node tables need. Step 2: the node tables in
// pieces.
//
// exp() is checked against mpmath at 60 digits within the derived 1.48 ulp (PHASE3_API.md B: the final rounding of
// 1 + p, 0.5 ulp; Horner's scheme, 0.73; the rounding of r, 0.25); the test allows 1.5. The tables' complements and
// weights are checked per node within the a-priori bound of B. The references are in
// quadrature_tanh_sinh_reference.hpp, from the script in docs/quadrature/design.md.

#include "quadrature_check.hpp"
#include "quadrature_tanh_sinh_reference.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>

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

    // The a-priori bound of PHASE3_API.md B on the relative error of a node's complement and weight, in eps (exp's
    // error taken as 1.5 ulp): for t >= 1/2, 2.5 + 4 sinh t + 3.5 e^-t and 6.2 + 4 sinh t + 3.5 e^-t; for t < 1/2,
    // where sinh and cosh come from their series of N additions, 2.5 + (N + 0.5) sinh t and 6.2 + (N + 0.5) sinh t.
    // The error of s = sinh t moves q = e^(-2s) by 2s times it, which dominates the outer nodes.
    struct NodeBound
    {
        double comp;
        double w;
    };

    NodeBound node_bound(double t, int series_additions)
    {
        if (t >= 0.5) {
            double const common = 4.0 * std::sinh(t) + 3.5 * std::exp(-t);
            return { 2.5 + common, 6.2 + common };
        }
        double const series = (series_additions + 0.5) * std::sinh(t);
        return { 2.5 + series, 6.2 + series };
    }

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
