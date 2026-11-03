// The exp-sinh table of nxx::quadrature (Phase 4 step 1, docs/quadrature/PHASE4_API.md B): one table for exp-sinh on a
// half line and sinh-sinh on the line, with u = c sinh t, c = 3/2. Its sampled nodes against mpmath within a derived
// a-priori bound, the nodes per level and side of every type within the capacity bound, compile-time and run-time
// tables bit-identical, and the latent fix of Phase 3b in exp_constants() (no upward factor past max()). The
// references are in quadrature_exp_sinh_reference.hpp, from the script in docs/quadrature/design.md ("The references
// of the exp-sinh tests").
//
// The bound, in eps, to first order (exp() within 1.5 ulp, PHASE3_API.md B; each other operation within 0.5 ulp):
// - sinh t and cosh t: below t = 1/2 by their series of N additions, each within (N + 1/2) / 2 (Phase 3's model);
//   above, from e = exp(t) and 1/e, ds = 1/2 + (3/2 e + 2/e) / (e - 1/e) and dc = 1/2 + (3/2 e + 2/e) / (e + 1/e).
// - u = c s, c exact: an absolute error a = u (ds + 1/2), which moves q, E and x by a times their logarithmic
//   derivative in u (1, 1, coth u): the same shift of the node in t for each value and its weight, as Phase 3's.
// - q = exp(-u): 3/2 + a. E = 1/q: 2 + a. f = c cosh t: dc + 1/2.
// - x = sinh u: below u = 1/2 its series at the computed u, (N + 1/2) / 2 + a coth u; above, (E - q) / 2,
//   1/2 + ((2 + a) E + (3/2 + a) q) / (E - q).
// - The weights as formed at use: f q, dc + 5/2 + a; f E, dc + 3 + a; (f E + f q) / 2, 1/2 + ((dc + 3 + a) E +
//   (dc + 5/2 + a) q) / (E + q).

#include "quadrature_check.hpp"
#include "quadrature_exp_sinh_reference.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <compare>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <span>
#include <utility>

using namespace quadrature_check;

namespace {

    namespace ref = exp_sinh_reference;

    double node_t(int level, int i)
    { return level == 0 ? i : (2.0 * i + 1.0) / static_cast<double>(1 << level); }

    constexpr auto inward_size = [](auto const &level) { return level.inward; };
    constexpr auto outward_size = [](auto const &level) { return level.outward; };

    // The bound of the file's comment on each value of a node.
    struct Bound
    {
        double q;
        double e;
        double f;
        double x;
        double w_end;
        double w_infinity;
        double w_line;
    };

    Bound bound_of(double t, int series_additions)
    {
        double const series = (series_additions + 0.5) / 2.0;
        double ds = series;
        double dc = series;
        if (t >= 0.5) {
            double const e = std::exp(t);
            double const rounding = 1.5 * e + 2.0 / e;
            ds = 0.5 + rounding / (e - 1.0 / e);
            dc = 0.5 + rounding / (e + 1.0 / e);
        }
        double const u = 1.5 * std::sinh(t);
        double const a = u * (ds + 0.5);
        double const r = std::exp(-2.0 * u); // q / E, so that E itself (up to 2^(2^31)) is never formed
        double const x = u < 0.5 ? series + a / std::tanh(u) : 0.5 + ((2.0 + a) + (1.5 + a) * r) / (1.0 - r);
        double const w_infinity = dc + 3.0 + a;
        double const w_end = dc + 2.5 + a;
        return { 1.5 + a, 2.0 + a, dc + 0.5, x, w_end, w_infinity, 0.5 + (w_infinity + w_end * r) / (1.0 + r) };
    }

    // The relative error of v against a reference value, in eps.
    template<typename T>
    double relative(T const &v, ref::Value const &reference)
    {
        return scaled_error(v, reference.mantissa, reference.exponent) / static_cast<double>(bin100(reference.mantissa))
               / eps_of<T>();
    }

    // Every sampled node of a table within its bound; returns the largest ratio of error to bound.
    template<typename T>
    double check_table(std::span<ref::Node const> nodes, int series_additions)
    {
        double worst = 0.0;
        for (ref::Node const &node : nodes) {
            INFO("level " << node.level << " node " << node.i);
            Bound const b = bound_of(node_t(node.level, node.i), series_additions);
            double const ratio =
                qd::detail::with_exp_sinh_level<T>(node.level, [&](qd::detail::ExpSinhLevel<T> const &level) {
                    REQUIRE(node.i < level.inward);
                    CHECK((node.i < level.outward) == node.outward);
                    T const q = level.q(node.i);
                    T const f = level.f(node.i);
                    T const w_end = T(f * q);
                    std::pair<double, double> const errors[] = { { relative(q, node.q), b.q },
                        { relative(f, node.f), b.f }, { relative(w_end, node.w_end), b.w_end } };
                    double largest = 0.0;
                    for (auto const &[error, bound] : errors) {
                        CHECK(error <= bound);
                        largest = std::max(largest, error / bound);
                    }
                    if (!node.outward || node.i >= level.outward) return largest;
                    T const e = level.e(node.i);
                    T const w_infinity = T(f * e);
                    T const w_line = T(T(w_infinity + w_end) / T(2));
                    std::pair<double, double> const far[] = { { relative(e, node.e), b.e },
                        { relative(level.x(node.i), node.x), b.x },
                        { relative(w_infinity, node.w_infinity), b.w_infinity },
                        { relative(w_line, node.w_line), b.w_line } };
                    for (auto const &[error, bound] : far) {
                        CHECK(error <= bound);
                        largest = std::max(largest, error / bound);
                    }
                    return largest;
                });
            worst = std::max(worst, ratio);
        }
        return worst;
    }

    // t = 0, the centre of both maps, exact in every type: q = E = 1, f = c = 3/2, x = 0.
    template<typename T>
    bool exact_centre()
    {
        return qd::detail::with_exp_sinh_level<T>(0, [](qd::detail::ExpSinhLevel<T> const &level) {
            return level.q(0) == T(1) && level.e(0) == T(1) && level.f(0) == T(T(3) / T(2)) && level.x(0) == T(0);
        });
    }

    // The nodes per level of each side, the script's counts, within the capacity bound.
    template<typename T>
    void check_counts(std::span<int const> inward, std::span<int const> outward)
    {
        for (int level = 0; level < qd::detail::tanh_sinh_levels; ++level) {
            INFO("level " << level);
            std::size_t const k = static_cast<std::size_t>(level);
            int const n_in = qd::detail::with_exp_sinh_level<T>(level, inward_size);
            int const n_out = qd::detail::with_exp_sinh_level<T>(level, outward_size);
            CHECK(n_in == inward[k]);
            CHECK(n_out == outward[k]);
            CHECK(n_out <= n_in);
            CHECK(n_in <= qd::detail::level_capacity<T>(level));
        }
    }

    // The compile-time pieces of a tabled type equal those generated at run time for its boxed twin (std::vectors on
    // the heap: a double piece's four columns, 32 KB, and its twin's would crowd Emscripten's 64 KB stack), bit for
    // bit.
    template<typename T, int L, int... P>
    bool same_level(std::integer_sequence<int, P...>)
    {
        static_assert(nxx::real<Boxed<T>> && !qd::detail::tabled<Boxed<T>>);
        auto const same = []<int Q>(std::integral_constant<int, Q>) {
            auto const &table = qd::detail::exp_sinh_table<T, L, Q>;
            auto const &generated = qd::detail::cached<&qd::detail::exp_sinh_piece<Boxed<T>, L, Q>>();
            if (table.inward != generated.inward || table.outward != generated.outward) return false;
            for (std::size_t i = 0; i < table.inward; ++i)
                if (!(table.q[i] == generated.q[i].v) || !(table.f[i] == generated.f[i].v)) return false;
            for (std::size_t i = 0; i < table.outward; ++i)
                if (!(table.e[i] == generated.e[i].v) || !(table.x[i] == generated.x[i].v)) return false;
            return true;
        };
        return (same(std::integral_constant<int, P>{}) && ...);
    }

    template<typename T, int... L>
    bool same_tables(std::integer_sequence<int, L...>)
    { return (same_level<T, L>(std::make_integer_sequence<int, qd::detail::level_pieces<T>(L)>{}) && ...); }

    // A type whose exponents reach further below 0 than above (max() = 2^100 with double's min()), for the latent fix
    // of exp_constants(): its reach, from min_exponent, needs factors up to 2^512, of which 2^128 and above pass max().
    // Its arithmetic is double's, so where the factors are formed they equal double's, and so does exp().
    struct Narrow
    {
        double v;
        constexpr Narrow(int i) : v(i) {}
        constexpr Narrow(double x) : v(x) {}
        friend constexpr Narrow operator+(Narrow x, Narrow y) { return Narrow(x.v + y.v); }
        friend constexpr Narrow operator-(Narrow x, Narrow y) { return Narrow(x.v - y.v); }
        friend constexpr Narrow operator*(Narrow x, Narrow y) { return Narrow(x.v * y.v); }
        friend constexpr Narrow operator/(Narrow x, Narrow y) { return Narrow(x.v / y.v); }
        friend constexpr Narrow operator-(Narrow x) { return Narrow(-x.v); }
        friend constexpr bool operator==(Narrow, Narrow) = default;
        friend constexpr auto operator<=>(Narrow x, Narrow y) { return x.v <=> y.v; }
    };

} // namespace

template<>
struct std::numeric_limits<Narrow> : std::numeric_limits<double>
{
    static constexpr Narrow max() { return 0x1p100; }
};

// The capacity bound's bits of a radix-10 digit are rounded up: 4, so that n ln 2 >= (1 - min_exponent) ln 10.
static_assert(qd::detail::significand_bits<dec50>() / std::numeric_limits<dec50>::digits == 4);

// The upward factors of exp() in Narrow stop at 2^64, where the next would pass max(); its exp() is double's.
static_assert([] {
    qd::detail::ExpConstants<Narrow> const narrow = qd::detail::exp_constants<Narrow>();
    qd::detail::ExpConstants<double> const wide_range = qd::detail::exp_constants<double>();
    bool same = narrow.up.size() == 11 && narrow.up[7] == Narrow(0x1p64) && narrow.up[8] == Narrow(0)
                && narrow.up[10] == Narrow(0);
    for (double const x : { -700.0, -1.5, 0.0, 0.75, 6.5, 21.5 })
        same = same && qd::detail::exp(Narrow(x), narrow).v == qd::detail::exp(x, wide_range);
    return same;
}());

// The double table in constant evaluation: the nodes per level and side (the script's counts), and level 9 across its
// two pieces (1,024 and 730 nodes): its last node toward the end normal, and toward infinity the last node with
// f <= (max / 2) q, the next refused.
static_assert(qd::detail::with_exp_sinh_level<double>(0, inward_size) == 7);
static_assert(qd::detail::with_exp_sinh_level<double>(0, outward_size) == 7);
static_assert(qd::detail::with_exp_sinh_level<double>(5, inward_size) == 110);
static_assert(qd::detail::with_exp_sinh_level<double>(5, outward_size) == 109);
static_assert(qd::detail::with_exp_sinh_level<double>(9, inward_size) == 1754);
static_assert(qd::detail::with_exp_sinh_level<double>(9, outward_size) == 1752);
static_assert(qd::detail::with_exp_sinh_level<double>(9, [](auto const &level) {
    double const half_max = std::numeric_limits<double>::max() / 2.0;
    return level.q(1753) >= std::numeric_limits<double>::min() && level.f(1751) <= half_max * level.q(1751)
           && level.f(1752) > half_max * level.q(1752);
}));
// A type without a table in constant evaluation: its level 4, generated as a temporary, equals the double table's.
static_assert(qd::detail::with_exp_sinh_level<Boxed<double>>(4, [](auto const &level) {
    return qd::detail::with_exp_sinh_level<double>(4, [&level](auto const &table) {
        bool same = level.inward == table.inward && level.outward == table.outward;
        for (int i = 0; same && i < table.inward; ++i) same = level.q(i).v == table.q(i) && level.f(i).v == table.f(i);
        for (int i = 0; same && i < table.outward; ++i) same = level.e(i).v == table.e(i) && level.x(i).v == table.x(i);
        return same;
    });
}));

TEST_CASE("quadrature: exp_sinh: the table against mpmath")
{
    // The sampled nodes within the bound: double (8 series additions below 1/2) and cpp_bin_float_50 (18).
    double const worst_double = check_table<double>(ref::table_double, 8);
    double const worst_bin50 = check_table<bin50>(ref::table_bin50, 18);
    MESSAGE("exp-sinh table, largest ratio of error to bound: double " << worst_double << ", cpp_bin_float_50 "
                                                                       << worst_bin50);
    CHECK(exact_centre<float>());
    CHECK(exact_centre<double>());
    CHECK(exact_centre<long double>());
    CHECK(exact_centre<bin50>());
    CHECK(exact_centre<bin50_et>());
    CHECK(exact_centre<dec50>());
}

TEST_CASE("quadrature: exp_sinh: the nodes per level")
{
    // Toward the finite end the nodes with q >= 2^-m, 2^-m the smallest power of 2 at or above min(), that is
    // u < m ln 2; toward infinity those with f <= (max / 2) q (the script's counts; no node lies within 1e-4 of
    // either boundary, relative, against x_min's 4 eps and the error of u). With t > 0: float 2,435 and 2,413,
    // double 3,507 and 3,503, long double (80-bit and binary128) 4,928 and 4,927, cpp_bin_float_50 10,961 on both
    // sides, cpp_dec_float_50 9,801 on both.
    check_counts<float>(ref::inward_float, ref::outward_float);
    check_counts<double>(ref::inward_double, ref::outward_double);
    if constexpr (std::numeric_limits<long double>::min_exponent == -16381)
        check_counts<long double>(ref::inward_long_double, ref::outward_long_double);
    else check_counts<long double>(ref::inward_double, ref::outward_double);
    check_counts<bin50>(ref::inward_bin50, ref::outward_bin50);
    check_counts<bin50_et>(ref::inward_bin50, ref::outward_bin50);
    check_counts<dec50>(ref::inward_dec50, ref::outward_dec50);
}

TEST_CASE("quadrature: exp_sinh: compile-time and run-time tables")
{
    CHECK(same_tables<float>(std::make_integer_sequence<int, qd::detail::tanh_sinh_levels>{}));
    CHECK(same_tables<double>(std::make_integer_sequence<int, qd::detail::tanh_sinh_levels>{}));
    CHECK(same_tables<long double>(std::make_integer_sequence<int, qd::detail::tanh_sinh_levels>{}));
}
