// The order of roots' failure estimates (DESIGN §6.6, §6.7, §7.2): nxx::better_than calls root_estimate's hidden friend,
// a strict weak order by the key (e, o, k, n, a): an enclosure first; the smaller width, and the smaller hi/2 - lo/2
// only when both widths overflow; then the smaller |f(x)|, with a NaN |f(x)| last. Being a strict weak order, it makes
// the best estimate of a chain independent of how its alternatives are grouped or folded.

#include <numerixx/core/any_solver.hpp>
#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <functional>
#include <limits>
#include <optional>
#include <random>
#include <tuple>
#include <type_traits>
#include <vector>

namespace
{
    namespace r = nxx::roots;

    template<class T>
    using est = r::root_estimate<T>;

    template<class T>
    constexpr T k_max = (std::numeric_limits<T>::max)();
    template<class T>
    constexpr T k_inf = std::numeric_limits<T>::infinity();
    template<class T>
    constexpr T k_nan = std::numeric_limits<T>::quiet_NaN();
    template<class T>
    constexpr T k_d = std::numeric_limits<T>::denorm_min();

    // An estimate with the enclosure [lo, hi] (finite ends, lo < hi), at lo, with the given f(x).
    template<class T>
    constexpr est<T> enclosed(T lo, T hi, T fx = T(1))
    { return est<T> { lo, fx, hi - lo, r::sign_bracket<T> { nxx::detail::trust_me {}, lo, T(-1), hi, T(1) } }; }

    // An estimate without an enclosure.
    template<class T>
    constexpr est<T> open_at(T fx)
    { return est<T> { T(0), fx }; }

    // ---- The CPO and the hidden friend ----------------------------------------------------------------------------
    static_assert(std::is_invocable_r_v<bool, decltype(nxx::better_than), const est<float>&, const est<float>&>);
    static_assert(std::is_invocable_r_v<bool, decltype(nxx::better_than), const est<double>&, const est<double>&>);
    static_assert(std::is_invocable_r_v<bool, decltype(nxx::better_than), const est<long double>&, const est<long double>&>);
    static_assert(nxx::detail::has_better_than_v<est<double>> && nxx::detail::has_better_than_v<const est<double>&>);
    static_assert(noexcept(nxx::better_than(std::declval<const est<double>&>(), std::declval<const est<double>&>())));
    // Not invocable on a type without an ADL better_than, and not across estimate types.
    static_assert(!std::is_invocable_v<decltype(nxx::better_than), const double&, const double&>);
    static_assert(!std::is_invocable_v<decltype(nxx::better_than), const r::sign_bracket<double>&, const r::sign_bracket<double>&>);
    static_assert(!std::is_invocable_v<decltype(nxx::better_than), const est<double>&, const est<float>&>);

    // A hidden friend: with both namespaces brought in by using-directives, the unqualified name finds only the object
    // nxx::better_than (a namespace-scope roots::better_than would make it ambiguous).
    constexpr bool unqualified_under_using_directives()
    {
        using namespace nxx;
        using namespace nxx::roots;
        return better_than(enclosed(0.0, 1.0), open_at(0.0)) && !better_than(open_at(0.0), enclosed(0.0, 1.0));
    }
    static_assert(unqualified_under_using_directives());

    // Rows below in constant expressions (a width that overflows is not a constant expression, so those rows run only
    // at run time).
    static_assert(nxx::better_than(open_at(1.0), open_at(k_nan<double>)) && !nxx::better_than(open_at(k_nan<double>), open_at(1.0)));
    static_assert(nxx::better_than(enclosed(0.0, k_max<double> / 2), enclosed(0.0, k_max<double>, 0.0)));
    static_assert(nxx::better_than(enclosed(k_d<double>, 3 * k_d<double>), enclosed(2 * k_d<double>, 5 * k_d<double>)));

    // ---- The key (e, o, k, n, a), written independently of the order ------------------------------------------------
    template<class T>
    struct key
    {
        int  e;    // 0 with an enclosure, 1 without
        int  o;    // 0 for a finite width, 1 for a width that overflows to inf
        T    k;    // the width when finite, hi/2 - lo/2 when it overflows
        int  n;    // 1 for a NaN |f(x)|
        T    a;    // |f(x)|, or 0 when NaN
        auto tie() const { return std::tie(e, o, k, n, a); }
    };

    template<class T>
    key<T> key_of(const est<T>& x)
    {
        key<T> out { 1, 0, T(0), 0, T(0) };
        if (x.enclosure) {
            const T w = x.enclosure->width();
            out.e     = 0;
            out.o     = nxx::math::isfinite(w) ? 0 : 1;
            out.k     = out.o == 0 ? w : x.enclosure->hi() / T(2) - x.enclosure->lo() / T(2);
        }
        const T fa = nxx::math::abs(x.fx);
        out.n      = nxx::math::isnan(fa) ? 1 : 0;
        out.a      = out.n == 1 ? T(0) : fa;
        return out;
    }

    // Draws straight from std::mt19937, whose output sequence the standard specifies (the distributions' is not), so
    // every standard library draws the same pool: an index below n, and a real in [lo, hi).
    std::size_t draw_index(std::mt19937& gen, std::size_t n) { return gen() % n; }
    double      draw_real(std::mt19937& gen, double lo, double hi) { return lo + (hi - lo) * (static_cast<double>(gen()) / 4294967296.0); }

    // Estimates that stress the order: f(x) NaN, +-inf, -0, 0, subnormal, ordinary or huge; no enclosure, or ends drawn
    // from extreme (+-max, +-max/2, +-max/4), ordinary, zero and subnormal values. Small pools make ties frequent.
    template<class T>
    std::vector<est<T>> stress_estimates(std::size_t count, std::uint32_t seed)
    {
        const T                 d = k_d<T>;
        const T                 m = k_max<T>;
        const std::array<T, 14> fxs { k_nan<T>, -k_nan<T>, k_inf<T>, -k_inf<T>, T(-0.0), T(0), d,
                                      -3 * d,   T(1e-3),   T(-1e-3), T(1),      T(2),    m,    -m };
        const std::array<T, 20> ends { -m,    -m / 2,  -m / 4, T(-1e30), T(-1), T(-0.5), T(-0.0), d,     2 * d, 3 * d,
                                       5 * d, T(0.25), T(0.5), T(1),     T(2),  T(1e30), m / 4,   m / 2, m,     T(7) };
        std::mt19937            gen(seed);
        const auto              fx_of = [&] {
            return draw_index(gen, 10) == 0 ? static_cast<T>(draw_real(gen, -100.0, 100.0)) : fxs[draw_index(gen, fxs.size())];
        };
        const auto end_of = [&] {
            return draw_index(gen, 10) == 0 ? static_cast<T>(draw_real(gen, -100.0, 100.0)) : ends[draw_index(gen, ends.size())];
        };
        std::vector<est<T>> out;
        out.reserve(count);
        while (out.size() < count) {
            const T fx = fx_of();
            if (draw_index(gen, 10) < 3) {
                out.push_back(open_at(fx));
                continue;
            }
            T lo = end_of();
            T hi = end_of();
            if (!(lo < hi) && !(hi < lo)) continue;    // equal ends (0 and -0 too): not a bracket
            if (hi < lo) std::swap(lo, hi);
            out.push_back(enclosed(lo, hi, fx));
        }
        return out;
    }

    // A failure that a curried solver returns as is, for chains over fixed estimates.
    template<class T>
    struct fixed_failure
    {
        est<T> best;

        template<class F>
        constexpr auto operator()(const F&) const -> nxx::result<est<T>>
        {
            return nxx::result<est<T>> { std::unexpect,
                                         nxx::failure<est<T>> { nxx::errc::budget_exhausted, nxx::algo::user_first, { 1, 1 }, best, {} } };
        }
    };

    template<class T>
    bool same_estimate(const est<T>& a, const est<T>& b)
    {
        const auto same = [](T u, T v) {
            return (nxx::math::isnan(u) && nxx::math::isnan(v)) || (u == v && std::signbit(u) == std::signbit(v));
        };
        return same(a.x, b.x) && same(a.fx, b.fx) && same(a.uncertainty, b.uncertainty) && a.enclosure == b.enclosure;
    }

    // first_of over the three estimates, in every order and grouping, statically and as a run-time chain: every fold
    // must select `want`.
    template<class T>
    void check_every_grouping(const std::array<est<T>, 3>& three, const est<T>& want)
    {
        const auto         f = [](T x) { return x; };
        std::array<int, 3> idx { 0, 1, 2 };
        using solver_t               = nxx::any_solver<std::function<T(T)>, est<T>>;
        const std::function<T(T)> sf = f;
        do {
            const fixed_failure<T>                   a { three[static_cast<std::size_t>(idx[0])] };
            const fixed_failure<T>                   b { three[static_cast<std::size_t>(idx[1])] };
            const fixed_failure<T>                   c { three[static_cast<std::size_t>(idx[2])] };
            const std::array<nxx::result<est<T>>, 4> results { nxx::first_of(a, b, c)(f),
                                                               nxx::first_of(nxx::first_of(a, b), c)(f),
                                                               nxx::first_of(a, nxx::first_of(b, c))(f),
                                                               nxx::first_of(std::vector<solver_t> { a, b, c })(sf) };
            for (const auto& res : results) {
                if (res || !res.error().best) {
                    FAIL_CHECK("a chain of failures fails with a best estimate");
                    continue;
                }
                CHECK(same_estimate(*res.error().best, want));
                CHECK(res.error().used == nxx::counters { 3, 3 });
            }
        }
        while (std::next_permutation(idx.begin(), idx.end()));
    }
}    // namespace

TEST_SUITE("roots")
{
    TEST_CASE_TEMPLATE("order: a NaN |f(x)| ranks last, in both argument orders", T, float, double, long double)
    {
        const T nan = k_nan<T>;
        CHECK(nxx::better_than(open_at(T(1)), open_at(nan)));
        CHECK_FALSE(nxx::better_than(open_at(nan), open_at(T(1))));
        CHECK(nxx::better_than(open_at(k_inf<T>), open_at(nan)));    // even an infinite |f(x)| ranks before NaN
        CHECK_FALSE(nxx::better_than(open_at(nan), open_at(k_inf<T>)));
        CHECK(nxx::better_than(enclosed(T(0), T(1), T(5)), enclosed(T(0), T(1), nan)));    // the same enclosure: |f(x)| decides
        CHECK_FALSE(nxx::better_than(enclosed(T(0), T(1), nan), enclosed(T(0), T(1), T(5))));
        CHECK(nxx::better_than(enclosed(T(0), T(1), nan), open_at(T(0))));    // the enclosure still comes first
        CHECK_FALSE(nxx::better_than(open_at(T(0)), enclosed(T(0), T(1), nan)));
    }

    TEST_CASE_TEMPLATE("order: two NaN |f(x)| are equivalent, and so are -0 and 0", T, float, double, long double)
    {
        const T nan = k_nan<T>;
        CHECK_FALSE(nxx::better_than(open_at(nan), open_at(-nan)));
        CHECK_FALSE(nxx::better_than(open_at(-nan), open_at(nan)));
        CHECK_FALSE(nxx::better_than(open_at(nan), open_at(nan)));
        CHECK_FALSE(nxx::better_than(enclosed(T(0), T(1), nan), enclosed(T(0), T(1), nan)));
        CHECK_FALSE(nxx::better_than(open_at(T(-0.0)), open_at(T(0))));
        CHECK_FALSE(nxx::better_than(open_at(T(0)), open_at(T(-0.0))));
    }

    TEST_CASE_TEMPLATE("order: when both widths overflow, the smaller hi/2 - lo/2 wins", T, float, double, long double)
    {
        const T    m     = k_max<T>;
        const auto whole = enclosed(-m, m);        // width inf, halves max
        const auto most  = enclosed(-m / 2, m);    // width inf, halves 0.75 max
        CHECK(whole.enclosure->width() == k_inf<T>);
        CHECK(most.enclosure->width() == k_inf<T>);
        CHECK(nxx::better_than(most, whole));
        CHECK_FALSE(nxx::better_than(whole, most));
        // The halves decide before |f(x)|: the wider enclosure loses even with the smaller |f(x)|.
        CHECK(nxx::better_than(enclosed(-m / 2, m, T(1)), enclosed(-m, m, T(0))));
    }

    TEST_CASE_TEMPLATE("order: a finite width ranks before a width that overflows", T, float, double, long double)
    {
        const T    m        = k_max<T>;
        const auto finite   = enclosed(T(0), m);      // width max
        const auto overflow = enclosed(-m / 2, m);    // width inf
        CHECK(finite.enclosure->width() == m);
        CHECK(nxx::better_than(finite, overflow));
        CHECK_FALSE(nxx::better_than(overflow, finite));
    }

    // The half-width key approved on 2026-10-04 inverted these two under round-to-nearest, ties-to-even ([d, 3d] has
    // half-width 2d, [2d, 5d] has d); the width key does not (2d < 3d). DESIGN §6.7, §12.21.
    TEST_CASE_TEMPLATE("order: [d, 3d] ranks before [2d, 5d] near the subnormal range", T, float, double, long double)
    {
        const T    d      = k_d<T>;
        const auto narrow = enclosed(d, 3 * d);
        const auto wide   = enclosed(2 * d, 5 * d);
        CHECK(narrow.enclosure->width() == 2 * d);
        CHECK(wide.enclosure->width() == 3 * d);
        CHECK(nxx::better_than(narrow, wide));
        CHECK_FALSE(nxx::better_than(wide, narrow));
    }

    // A fixed-seed property test of the four axioms of a strict weak order, of agreement with the key (e, o, k, n, a),
    // and of the width rule: a strictly nested (hence strictly narrower in exact arithmetic) enclosure never ranks after
    // the one around it unless both have the same (o, k), which the half-width key would also pass the axioms without.
    TEST_CASE_TEMPLATE("order: better_than is a strict weak order that never prefers a wider enclosure (property)",
                       T,
                       float,
                       double,
                       long double)
    {
        constexpr std::uint32_t pool_seed   = 20261006u;
        constexpr std::uint32_t nested_seed = 20261007u;
        CAPTURE(pool_seed);
        CAPTURE(nested_seed);
        const std::vector<est<T>> v  = stress_estimates<T>(160, pool_seed);
        const std::size_t         n  = v.size();
        const auto                lt = [](const est<T>& a, const est<T>& b) { return nxx::better_than(a, b); };

        std::vector<std::uint8_t> less(n * n);
        std::size_t               key_mismatch = 0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                less[i * n + j] = lt(v[i], v[j]) ? 1 : 0;
                if ((less[i * n + j] == 1) != (key_of(v[i]).tie() < key_of(v[j]).tie())) ++key_mismatch;
            }
        CHECK(key_mismatch == 0u);

        std::size_t irreflexive  = 0;
        std::size_t asymmetric   = 0;
        std::size_t transitive   = 0;
        std::size_t incomparable = 0;
        const auto  at           = [&](std::size_t i, std::size_t j) { return less[i * n + j] == 1; };
        const auto  equiv        = [&](std::size_t i, std::size_t j) { return !at(i, j) && !at(j, i); };
        for (std::size_t i = 0; i < n; ++i) {
            if (at(i, i)) ++irreflexive;
            for (std::size_t j = 0; j < n; ++j)
                if (at(i, j) && at(j, i)) ++asymmetric;
        }
        // The two transitivity axioms are O(n^3): they run over every third estimate, 54 of the 160. The key check above
        // covers all 160, and the pool holds an enclosure whose width overflows, so the key's o = 1 is exercised.
        CHECK(std::ranges::any_of(v, [](const est<T>& e) { return e.enclosure && !nxx::math::isfinite(e.enclosure->width()); }));
        for (std::size_t i = 0; i < n; i += 3)
            for (std::size_t j = 0; j < n; j += 3)
                for (std::size_t k = 0; k < n; k += 3) {
                    if (at(i, j) && at(j, k) && !at(i, k)) ++transitive;
                    if (equiv(i, j) && equiv(j, k) && !equiv(i, k)) ++incomparable;
                }
        CHECK(irreflexive == 0u);
        CHECK(asymmetric == 0u);
        CHECK(transitive == 0u);
        CHECK(incomparable == 0u);

        // Nested pairs: inner = [lo2, hi2] inside outer = [lo, hi], with lo <= lo2 < hi2 <= hi and not both equal.
        std::mt19937 gen(nested_seed);
        std::size_t  nested   = 0;
        std::size_t  inverted = 0;
        for (const est<T>& outer : v) {
            if (!outer.enclosure) continue;
            const T                lo  = outer.enclosure->lo();
            const T                hi  = outer.enclosure->hi();
            const T                mid = nxx::math::midpoint(lo, hi);
            const std::array<T, 5> cut { lo, nxx::math::midpoint(lo, mid), mid, nxx::math::midpoint(mid, hi), hi };
            for (int trial = 0; trial < 8; ++trial) {
                T a = cut[draw_index(gen, cut.size())];
                T b = cut[draw_index(gen, cut.size())];
                if (b < a) std::swap(a, b);
                if (!(a < b) || (a == lo && b == hi)) continue;
                for (const est<T>& other : v) {    // every f(x) for the inner one, NaN included
                    const est<T> inner = enclosed(a, b, other.fx);
                    ++nested;
                    const key<T> ko = key_of(outer);
                    const key<T> ki = key_of(inner);
                    if (lt(outer, inner) && !(ko.o == ki.o && ko.k == ki.k)) ++inverted;
                }
            }
        }
        CHECK(nested > 1000u);
        CHECK(inverted == 0u);
    }

    TEST_CASE_TEMPLATE("order: first_of over three extreme estimates selects the same best in every grouping",
                       T,
                       float,
                       double,
                       long double)
    {
        const T m = k_max<T>;
        const T d = k_d<T>;
        // Widths max (finite), inf with halves 0.75 max, inf with halves max.
        check_every_grouping<T>({ enclosed(-m, m), enclosed(-m / 2, m), enclosed(T(0), m) }, enclosed(T(0), m));
        // Both widths overflow, against an estimate without an enclosure and an exact zero.
        check_every_grouping<T>({ enclosed(-m, m), open_at(T(0)), enclosed(-m / 2, m) }, enclosed(-m / 2, m));
        // No enclosures: NaN, inf and a subnormal |f(x)|.
        check_every_grouping<T>({ open_at(k_nan<T>), open_at(k_inf<T>), open_at(d) }, open_at(d));
        // Subnormal enclosures, with a NaN |f(x)| on the narrower width.
        check_every_grouping<T>({ enclosed(d, 3 * d, k_nan<T>), enclosed(2 * d, 5 * d, T(0)), enclosed(d, 3 * d, T(1)) },
                                enclosed(d, 3 * d, T(1)));
    }
}
