/**
 * @file
 * @brief Numerical derivatives of a scalar function by finite differences.
 *
 * derivative() estimates f'(x) from values of fn near x, by a stencil: `Central<P>`, or `Forward<P>` and `Backward<P>`
 * where fn is defined on one side of x only. Results are std::expected; Numerixx itself never throws.
 *
 * @code
 * auto const d = nxx::deriv::derivative<nxx::deriv::Central<4>>([](double x) { return std::exp(x); }, 1.0);
 * if (d) use(*d); // e to about 1e-13
 * @endcode
 */

#pragma once

#include <numerixx/math.hpp>

#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <expected>
#include <functional>
#include <limits>
#include <numeric>
#include <type_traits>

/** @brief Numerical derivatives of scalar functions. */
namespace nxx::deriv {

    /** @brief Why derivative() returned no estimate. */
    enum class error {
        /**
         * No finite estimate: x is not finite, fn returned an infinite or NaN value, the sum or the quotient
         * overflowed, a point x + k h would overflow, or h is below the resolution of x.
         */
        non_finite
    };

    /** @cond INTERNAL */
    namespace detail {

        // The orders the stencil families accept. Beyond the caps the gain is under half a digit even for 50-digit
        // types (central), or the rounding amplification is too large (one-sided: 27.7 for P = 6).
        constexpr bool valid_central(int p) { return p >= 2 && p <= 10 && p % 2 == 0; }
        constexpr bool valid_one_sided(int p) { return p >= 1 && p <= 6; }

    } // namespace detail
    /** @endcond */

    /**
     * @brief Central stencil of accuracy order @p P, the default stencil of derivative().
     *
     * Nodes -P/2 ... P/2 without 0: P evaluations of fn, and an error of about eps^(P/(P+1)) with the default step.
     * `Central<2>` is the central difference (fn(x + h) - fn(x - h)) / 2h.
     *
     * @tparam P The accuracy order: even, 2 <= P <= 10.
     */
    template<int P>
    struct Central
    {
        static_assert(detail::valid_central(P), "Central<P> needs an even order P from 2 to 10");
    };

    /**
     * @brief One-sided stencil of accuracy order @p P that never evaluates fn below x.
     *
     * Nodes 0 ... P: P + 1 evaluations of fn, including fn(x). For fn defined on one side of x only: the error is
     * larger than for a Central stencil of the same order, and order 1 reaches only about eps^(1/2).
     *
     * @tparam P The accuracy order: 1 <= P <= 6.
     */
    template<int P>
    struct Forward
    {
        static_assert(detail::valid_one_sided(P), "Forward<P> needs an order P from 1 to 6");
    };

    /**
     * @brief `Forward<P>` mirrored: nodes -P ... 0, so fn is never evaluated above x.
     *
     * @tparam P The accuracy order: 1 <= P <= 6.
     */
    template<int P>
    struct Backward
    {
        static_assert(detail::valid_one_sided(P), "Backward<P> needs an order P from 1 to 6");
    };

    /** @cond INTERNAL */
    namespace detail {

        template<typename S>
        inline constexpr bool is_stencil = false;
        template<int P>
        inline constexpr bool is_stencil<Central<P>> = valid_central(P);
        template<int P>
        inline constexpr bool is_stencil<Forward<P>> = valid_one_sided(P);
        template<int P>
        inline constexpr bool is_stencil<Backward<P>> = valid_one_sided(P);

        template<typename S>
        concept stencil = is_stencil<S>;

        // f'(x) ~= sum_k numerators[k] * fn(x + nodes[k] * h) / (denominator * h), with error O(h^order). The weights
        // are integers over one denominator, so they are exact in every real type and the evaluator divides once.
        template<std::size_t N>
        struct Table
        {
            std::array<int, N> nodes; // strictly increasing
            std::array<int, N> numerators;
            int denominator;
            int order;
        };

        constexpr long long factorial(std::size_t n)
        {
            long long f = 1;
            for (std::size_t i = 2; i <= n; ++i) f *= static_cast<long long>(i);
            return f;
        }

        // make_table builds the Table of a stencil at compile time. The weights of a first-derivative rule are fixed by
        // its nodes: they are the unique solution of the moment conditions M_0 = 0, M_1 = denominator and M_j = 0 for
        // 2 <= j <= order (see derivative() and docs/DERIV_TUTORIAL.md, section 8). Equivalently, they come from
        // differentiating the polynomial that interpolates fn at the nodes. For the three families that solution has a
        // closed form, so each overload evaluates the closed form in exact rational arithmetic (long long numerators
        // and denominators, reduced with gcd), then puts all weights over one common denominator (their lcm), so that
        // the table holds integers only. The largest intermediate is 10! times 5 (Central<10>), far inside long long,
        // and the largest numerator is 2100, exact in every real type.

        // Central, P = 2m, on the nodes -m ... -1, 1 ... m. A first derivative is odd, so the weights are
        // antisymmetric, w(-k) = -w(k), and the rule reads f'(x) ~= sum_{k=1..m} a_k (f(x + k h) - f(x - k h)) / h with
        //
        //     a_k = (-1)^(k+1) (m!)^2 / (k (m-k)! (m+k)!).
        //
        // Worked example, Central<4> (m = 2): a_1 = (2!)^2 / (1 * 1! * 3!) = 4/6 = 2/3 and a_2 = -(2!)^2 / (2 * 0! *
        // 4!) = -4/48 = -1/12. The common denominator is lcm(3, 12) = 12, so a_1 = 8/12 and a_2 = -1/12, and the table
        // is nodes {-2, -1, 1, 2}, numerators {1, -8, 8, -1}, denominator 12.
        template<int P>
        constexpr Table<std::size_t{ P }> make_table(Central<P>)
        {
            constexpr std::size_t m = P / 2;

            // Step 1: a_k = num[k] / den[k] in lowest terms, for k = 1 ... m (index 0 is unused), and the running lcm
            // of the denominators.
            std::array<long long, m + 1> num{};
            std::array<long long, m + 1> den{};
            long long common = 1;
            for (std::size_t k = 1; k <= m; ++k) {
                long long n = factorial(m) * factorial(m);
                long long d = static_cast<long long>(k) * factorial(m - k) * factorial(m + k);
                long long const g = std::gcd(n, d);
                n /= g;
                d /= g;
                num[k] = k % 2 == 1 ? n : -n; // the sign (-1)^(k+1)
                den[k] = d;
                common = std::lcm(common, d);
            }

            // Step 2: scale every a_k to the common denominator, w = a_k * common, and store the pair of nodes -k and
            // +k with the weights -w and +w. The nodes increase through the table: -m ... -1 at indices 0 ... m - 1,
            // then 1 ... m at indices m ... 2m - 1, so node -k goes to index m - k and node +k to index m + k - 1.
            Table<std::size_t{ P }> t{};
            for (std::size_t k = 1; k <= m; ++k) {
                int const node = static_cast<int>(k);
                int const w = static_cast<int>(num[k] * (common / den[k])); // exact: den[k] divides common
                t.nodes[m - k] = -node;
                t.numerators[m - k] = -w;
                t.nodes[m + k - 1] = node;
                t.numerators[m + k - 1] = w;
            }
            t.denominator = static_cast<int>(common);
            t.order = P;
            return t;
        }

        // Forward, P, on the nodes 0, 1, ..., P. The rule is f'(x) ~= sum_{k=0..P} w_k f(x + k h) / h with
        //
        //     w_0 = -(1 + 1/2 + ... + 1/P)   and   w_k = (-1)^(k+1) C(P, k) / k,  k = 1 ... P.
        //
        // The denominator of C(P, k) / k divides k, and that of the harmonic sum in w_0 divides lcm(1, ..., P), so
        // lcm(1, ..., P) is a common denominator of all the weights and no gcd reduction is needed.
        //
        // Worked example, Forward<2>: common = lcm(1, 2) = 2. k = 1: C(2, 1) = 2, so w_1 = +2/1 = 4/2; k = 2: C(2, 2) =
        // 1, so w_2 = -1/2; w_0 = -(1 + 1/2) = -(2/2 + 1/2) = -3/2. The table is nodes {0, 1, 2}, numerators
        // {-3, 4, -1}, denominator 2: f'(x) ~= (-3 f(x) + 4 f(x + h) - f(x + 2h)) / 2h.
        template<int P>
        constexpr Table<std::size_t{ P + 1 }> make_table(Forward<P>)
        {
            constexpr std::size_t p = P;
            long long common = 1;
            for (std::size_t k = 1; k <= p; ++k) common = std::lcm(common, static_cast<long long>(k));

            Table<std::size_t{ P + 1 }> t{};
            long long w0 = 0; // w_0 * common, accumulated term by term
            long long binomial = 1; // C(P, k), updated from C(P, k - 1)
            for (std::size_t k = 1; k <= p; ++k) {
                long long const kk = static_cast<long long>(k);
                // C(P, k) = C(P, k - 1) (P - k + 1) / k; the division is exact, since C(P, k - 1) (P - k + 1) = k C(P,
                // k).
                binomial = binomial * (P - kk + 1) / kk;
                // w_k * common = (-1)^(k+1) C(P, k) * (common / k); common / k is exact because k divides common.
                long long const w = (k % 2 == 1 ? 1 : -1) * binomial * (common / kk);
                t.nodes[k] = static_cast<int>(k);
                t.numerators[k] = static_cast<int>(w);
                w0 -= common / kk; // the term -1/k of w_0, times common
            }
            t.nodes[0] = 0;
            t.numerators[0] = static_cast<int>(w0);
            t.denominator = static_cast<int>(common);
            t.order = P;
            return t;
        }

        // Backward, P, on the nodes -P, ..., -1, 0: Forward<P> with h replaced by -h. Substituting -h into the forward
        // rule gives f'(x) ~= sum_k w_k f(x - k h) / (-h) = sum_k (-w_k) f(x + (-k) h) / h, so each node k becomes -k
        // and each weight w_k becomes -w_k. The forward table is read from its end so that the nodes still increase,
        // which the evaluator relies on (nodes[0] and nodes[size - 1] are the outermost nodes).
        //
        // Worked example, Backward<2>: from Forward<2>, nodes {0, 1, 2} with numerators {-3, 4, -1}, to nodes
        // {-2, -1, 0} with numerators {1, -4, 3}, denominator 2: f'(x) ~= (f(x - 2h) - 4 f(x - h) + 3 f(x)) / 2h.
        template<int P>
        constexpr Table<std::size_t{ P + 1 }> make_table(Backward<P>)
        {
            constexpr std::size_t p = P;
            Table<std::size_t{ P + 1 }> const f = detail::make_table(Forward<P>{});
            Table<std::size_t{ P + 1 }> t{};
            for (std::size_t i = 0; i <= p; ++i) {
                t.nodes[i] = -f.nodes[p - i]; // index i of the result is index p - i of the forward table, mirrored
                t.numerators[i] = -f.numerators[p - i];
            }
            t.denominator = f.denominator;
            t.order = f.order;
            return t;
        }

        // A rule of order p has truncation error O(h^p) and rounding error O(eps / h). A step of
        // eps^(1/(p+1)) * max(|x|, 1) balances the two, scaled so that the step is relative for large |x|.
        //
        // In detail, the error is about C h^p |f^(p+1)| + S eps |f| / h, where C is the stencil's error constant and S
        // the sum of its absolute weights. It is smallest at h = (S eps |f| / (p C |f^(p+1)|))^(1/(p+1)), with an
        // error of order eps^(p/(p+1)). f and its derivatives are unknown, so the step assumes that f varies on the
        // scale max(|x|, 1) and drops the constant (S / (p C))^(1/(p+1)), which lies between 1.4 and 2.0 for these
        // stencils and costs at most a factor 1.64 in error. For Central<2>: h = eps^(1/3) * max(|x|, 1).
        template<int Order, real T>
        constexpr T default_step(T const &x)
        {
            T const scale = math::abs(x) < T(1) ? T(1) : math::abs(x);
            return math::root_epsilon<Order + 1, T>() * scale;
        }

    } // namespace detail
    /** @endcond */

    /**
     * @brief The derivative of @p fn at @p x by stencil @p S with step @p h.
     *
     * fn is called once per node of S, and only at finite points. Failures are detected before the arithmetic that
     * would produce them, so they are error values in constant expressions too, except an overflow of the final sum or
     * quotient, which GCC rejects at compile time. Exceptions thrown by fn or by T pass through unchanged.
     *
     * @tparam S The stencil: `Central<P>` (by default `Central<2>`), `Forward<P>` or `Backward<P>`.
     * @tparam T The real type, deduced from @p x.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @param fn The function to differentiate.
     * @param x The point.
     * @param h The step; it converts to T, so a double step works for a multiprecision T.
     * @pre h > 0, and h well above ulp(x).
     * @return The estimate of f'(x), or error::non_finite.
     */
    // By Taylor's theorem, sum_k w_k f(x + n_k h) = sum_j M_j h^j f^(j)(x) / j!, with the moments M_j = sum_k w_k
    // n_k^j. The weights make M_0 = 0, M_1 = denominator and M_j = 0 for 2 <= j <= order, so the sum divided by
    // denominator * h is f'(x) plus a truncation error of order h^order. For Central<2>, the even terms cancel:
    // (f(x + h) - f(x - h)) / 2h = f'(x) + h^2/6 f'''(x) + O(h^4).
    //
    // F & in the constraint, because fn is called as a named reference, an lvalue. h is a non-deduced context
    // (type_identity_t), so T comes from x alone and a double step converts to a multiprecision T instead of making
    // deduction fail.
    template<typename S = Central<2>, real T, typename F>
    requires detail::stencil<S> && std::invocable<F &, T const &>
             && std::convertible_to<std::invoke_result_t<F &, T const &>, T>
    [[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
    {
        assert(h > T(0)); // a non-positive step is a bug in the caller, not a run-time failure

        constexpr auto table = detail::make_table(S{});
        constexpr std::size_t size = table.nodes.size();
        constexpr int lo = table.nodes[0];
        constexpr int hi = table.nodes[size - 1];

        // With K = max(|lo|, |hi|) and Q = max(2K, denominator), the check h <= max / 2Q and |x| <= max - Q h keeps
        // every point and the divisor finite, so nothing overflows before fn is called. A point lies within K h of x,
        // so it stays at least about K h below max: max - Q h is rounded by at most half an ulp of max, and K h of
        // margin keeps a rounded point from tying up to infinity, as it could with a margin of 0 (Q = K). The divisor
        // is at most Q h <= max / 2. The check cannot overflow itself: Q h <= max / 2, and max - Q h only subtracts.
        // For Central<2>, Q = 2.
        constexpr int k_max = -lo > hi ? -lo : hi;
        constexpr int q = 2 * k_max > table.denominator ? 2 * k_max : table.denominator;
        T const largest = (std::numeric_limits<T>::max)();
        if (!(h <= largest / T(2 * q) && math::abs(x) <= largest - T(q) * h)) return std::unexpected(error::non_finite);

        // x + k h; node 0 is x itself.
        auto const point = [&x, &h](int k) -> T {
            if (k == 0) return x;
            return T(x + T(k) * h);
        };

        if constexpr (size == 2) {
            // Divide by the distance between the points actually evaluated, not by (hi - lo) h: the points are rounded,
            // and at x = 1 with h = 1e-5 their spacing can differ from 2h by about 1e-11 relative, as much as the error
            // aimed for. With two nodes, the weights are -1 and +1 over hi - lo, so this is the whole formula.
            T const x_lo = point(lo);
            T const x_hi = point(hi);
            T const dx = x_hi - x_lo;
            if (dx == T(0)) return std::unexpected(error::non_finite); // h is below the resolution of x
            // fn may return an expression template, and every intermediate is a T const, never auto: an unevaluated
            // expression would hold references to its operands, which can dangle. fn is called through std::invoke
            // with a T const lvalue, exactly as the constraint std::invocable<F &, T const &> promises, so member
            // pointers work and so does a callable that rejects rvalues.
            T const f_lo = static_cast<T>(std::invoke(fn, x_lo));
            T const f_hi = static_cast<T>(std::invoke(fn, x_hi));
            // Checked before subtracting, so inf - inf is never formed: constant evaluation rejects arithmetic that
            // produces a NaN (and on GCC an overflow), so the failure would be a compile error instead of an error
            // value.
            if (!math::isfinite(f_lo) || !math::isfinite(f_hi)) return std::unexpected(error::non_finite);
            T const d = (f_hi - f_lo) / dx;
            if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the difference or quotient overflowed
            return d;
        } else {
            // With more nodes, dividing by the outer-node distance would cancel the outer-node rounding only for
            // two-node rules and enlarge it for the others (4x for Central<4>), so the divisor is denominator * h, in a
            // single division: with dyadic x and a power-of-two h the result is then exact for polynomials of degree
            // up to the order. Node rounding adds about S / 2 * eps * |x f'| / h, the same order as fn's rounding.
            // h is below the resolution of x. Both sides are checked: at a power of the radix the spacing differs on
            // the two sides of x, so h can move x one way and not the other (at x = -1, h = 0x1.8p-54 moves x up but
            // not down), which would leave Backward and Central rules with coinciding points and a wrong result.
            if (T(x + h) == x || T(x - h) == x) return std::unexpected(error::non_finite);
            T sum = T(0);
            for (std::size_t k = 0; k < size; ++k) {
                T const xk = point(table.nodes[k]); // a T const lvalue, as for two nodes
                T const fk = static_cast<T>(std::invoke(fn, xk));
                // Checked before it is added, so inf - inf is never formed (see above).
                if (!math::isfinite(fk)) return std::unexpected(error::non_finite);
                sum = T(sum + T(table.numerators[k]) * fk);
            }
            T const d = sum / (T(table.denominator) * h);
            if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the sum or the quotient overflowed
            return d;
        }
    }

    /**
     * @brief The derivative of @p fn at @p x by stencil @p S with the default step.
     *
     * The step is eps^(1/(P+1)) * max(|x|, 1) for a stencil of order P. It balances truncation against rounding when
     * fn varies on the scale max(|x|, 1), and the error is then about eps^(P/(P+1)); pass a step when fn varies
     * faster, or when the derivative is much smaller than that error.
     *
     * @tparam S The stencil: `Central<P>` (by default `Central<2>`), `Forward<P>` or `Backward<P>`.
     * @tparam T The real type, deduced from @p x.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @param fn The function to differentiate.
     * @param x The point.
     * @return The estimate of f'(x), or error::non_finite.
     */
    template<typename S = Central<2>, real T, typename F>
    requires detail::stencil<S> && std::invocable<F &, T const &>
             && std::convertible_to<std::invoke_result_t<F &, T const &>, T>
    [[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x)
    {
        // Checked first: a NaN x would give a NaN step and trip the h > 0 precondition of the overload below.
        if (!math::isfinite(x)) return std::unexpected(error::non_finite);
        // Qualified, so that argument-dependent lookup cannot pick a derivative() from the namespace of T or F.
        return nxx::deriv::derivative<S>(fn, x, detail::default_step<detail::make_table(S{}).order>(x));
    }

} // namespace nxx::deriv
