/**
 * @file
 * @brief Numerical derivatives of a scalar function by finite differences.
 *
 * derivative() estimates f'(x) from values of fn near x, by a method: a stencil, `Central<P>`, or `Forward<P>` and
 * `Backward<P>` where fn is defined on one side of x only; or `Ridders`, which adapts the step and also estimates the
 * error. Every method returns a std::expected<Estimate<T>, error>; Numerixx itself never throws.
 *
 * @code
 * auto const d = nxx::deriv::derivative<nxx::deriv::Ridders>([](double x) { return std::exp(x); }, 1.0);
 * if (d) use(d->value, *d->abs_error); // e to about 1e-15, with an error estimate of about 6e-14
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
#include <optional>
#include <type_traits>
#include <utility>

/** @brief Numerical derivatives of scalar functions. */
namespace nxx::deriv {

    /**
     * @brief Why derivative() returned no estimate.
     *
     * More values may be added; a switch on it should have a default.
     */
    enum class error {
        /**
         * No finite estimate: x is not finite, fn returned an infinite or NaN value, the sum or the quotient
         * overflowed, a point x + k h would overflow, or h is below the resolution of x.
         */
        non_finite
    };

    /**
     * @brief A derivative estimate: the value, and an estimate of its absolute error when the method provides one.
     *
     * Converts implicitly to T. Where a conversion is not considered (auto, std::format, and the operator templates of
     * Boost.Multiprecision), write `e.value` or `T(e)`. The conversion also applies where it may surprise: `if (*r)`
     * tests the value against 0, not whether there is a result, and `int n = *r` truncates the value.
     *
     * @code
     * auto const r = nxx::deriv::derivative<nxx::deriv::Ridders>(fn, 1.0);
     * if (r) {
     *     double const slope = *r;                          // the conversion
     *     double const bound = r->abs_error.value_or(0.0); // present for Ridders, empty for the stencils
     * }
     * @endcode
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Estimate
    {
        /** @brief The estimate of f'(x). */
        T value;

        /** @brief An estimate of |value - f'(x)|, not a guaranteed bound; empty for the stencils. */
        std::optional<T> abs_error;

        /**
         * @brief The value, so that an Estimate can be used where a T is expected.
         * @return value.
         */
        constexpr operator T() const { return value; }

        /** @brief Compares both members; without it, == would compare only the values, through the conversion. */
        friend constexpr bool operator==(Estimate const &, Estimate const &) = default;
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

    /**
     * @brief Ridders' extrapolation, with an error estimate.
     *
     * Central differences at h, h/2, h/4, ... are extrapolated in Richardson's tableau until rounding stops the
     * improvement, so the step need not be tuned: the initial step should be larger than the best step of the central
     * difference, and by default it is max(|x|, 1) / 32. The method stops when an extrapolated entry equals its parent,
     * after two levels that improve nothing, or at the resolution of x. It restarts below a step at which fn is not
     * finite, and below one whose result has at most two correct bits by its own estimate while that estimate is far
     * above rounding (a pole between the samples, aliasing).
     *
     * It costs 2 evaluations of fn per level, typically 10 to 16 for double and at most 64. The result's abs_error is
     * the change of the best entry from its parent, but at least the rounding error of the difference it came from and
     * at least its distance to the latest diagonal entry. It assumes that fn is accurate to about epsilon, and is
     * optimistic for a noisier fn. Like the stencils, it runs in constant expressions.
     *
     * @code
     * // The default initial step 1/32 reaches below 0, where log is NaN, so the method restarts below it.
     * auto const r = nxx::deriv::derivative<nxx::deriv::Ridders>([](double x) { return std::log(x); }, 0.01);
     * if (r) use(r->value, *r->abs_error); // 100 to about 1e-13 relative
     * @endcode
     */
    struct Ridders
    {
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

        /**
         * @brief The methods derivative() accepts: the stencils and Ridders.
         *
         * Closed on purpose: users cannot add a method.
         *
         * @tparam M The method type.
         */
        template<typename M>
        concept method = stencil<M> || std::same_as<M, Ridders>;

        /**
         * @brief A function derivative() can differentiate at a T.
         *
         * F & in the constraint, because fn is called as a named reference, an lvalue, through std::invoke with a T
         * const lvalue: member pointers work, and so does a callable that rejects rvalues.
         *
         * @tparam F The callable type.
         * @tparam T The real type.
         */
        template<typename F, typename T>
        concept callable =
            std::invocable<F &, T const &> && std::convertible_to<std::invoke_result_t<F &, T const &>, T>;

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

        /**
         * @brief max(|x|, 1): the scale on which the default steps assume fn varies.
         *
         * So that a step is relative for large |x| and absolute near 0.
         *
         * @tparam T The real type.
         * @param x The point.
         * @return max(|x|, 1).
         */
        template<real T>
        constexpr T scale(T const &x)
        { return math::abs(x) < T(1) ? T(1) : math::abs(x); }

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
        { return math::root_epsilon<Order + 1, T>() * detail::scale(x); }

        /**
         * @brief Whether every point x + k h and the divisor of a rule are finite: h <= max / 2q and |x| <= max - q h.
         *
         * With K = max(|lo|, |hi|) the largest node and q = max(2K, denominator), nothing overflows before fn is
         * called. A point lies within K h of x, so it stays at least about K h below max: max - q h is rounded by at
         * most half an ulp of max, and K h of margin keeps a rounded point from tying up to infinity, as it could with
         * a margin of 0 (q = K). The divisor is at most q h <= max / 2. The check cannot overflow itself: q h <= max /
         * 2, and max - q h only subtracts. For Central<2>, and for every level of Ridders, q = 2.
         *
         * @tparam T The real type.
         * @param x The point, finite.
         * @param h The step, > 0 and finite.
         * @param q max(2K, denominator) of the rule.
         * @return true if the rule can be evaluated without overflow.
         */
        template<real T>
        constexpr bool fits(T const &x, T const &h, int q)
        {
            T const largest = (std::numeric_limits<T>::max)();
            return h <= largest / T(2 * q) && math::abs(x) <= largest - T(q) * h;
        }

        /**
         * @brief A two-point difference, with the values and the spacing it came from.
         *
         * Ridders needs the values and the spacing for its rounding floor; the two-node stencils use the value only.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct Difference
        {
            /** @brief (f_hi - f_lo) / dx. */
            T value;

            /** @brief fn at the lower point. */
            T f_lo;

            /** @brief fn at the upper point. */
            T f_hi;

            /** @brief The distance between the points as evaluated, > 0. */
            T dx;
        };

        /**
         * @brief The rounding floor of a difference: its rounding error when fn is accurate to about eps.
         *
         * (eps |f_lo| + eps |f_hi|) / dx, eps applied term by term so that the sum cannot overflow, and saturated at
         * max where the division would overflow (a large fn over a tiny step), which constant evaluation would reject.
         *
         * @tparam T The real type.
         * @param d The difference.
         * @return The floor, at most max.
         */
        template<real T>
        constexpr T rounding_floor(Difference<T> const &d)
        {
            T const eps = std::numeric_limits<T>::epsilon();
            T const largest = (std::numeric_limits<T>::max)();
            T const sum = eps * math::abs(d.f_lo) + eps * math::abs(d.f_hi);
            if (T(sum / largest) >= T(d.dx / T(2))) return largest;
            return sum / d.dx;
        }

        /**
         * @brief The difference (fn(x_hi) - fn(x_lo)) / (x_hi - x_lo), the central difference and both two-node rules.
         *
         * Divides by the distance between the points actually evaluated, not by the nominal one: the points are
         * rounded, and at x = 1 with h = 1e-5 their spacing can differ from 2h by about 1e-11 relative, as much as the
         * error aimed for.
         *
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x_lo The lower point, finite.
         * @param x_hi The upper point, finite.
         * @return The difference with fn's values and the spacing, or error::non_finite if the points coincide, fn is
         *         not finite at one of them, or the quotient overflows.
         */
        template<real T, typename F>
        constexpr std::expected<Difference<T>, error> difference(F &fn, T const &x_lo, T const &x_hi)
        {
            T const dx = x_hi - x_lo;
            if (dx == T(0)) return std::unexpected(error::non_finite); // h is below the resolution of x
            // fn may return an expression template, and every intermediate is a T const, never auto: an unevaluated
            // expression would hold references to its operands, which can dangle.
            T const f_lo = static_cast<T>(std::invoke(fn, x_lo));
            T const f_hi = static_cast<T>(std::invoke(fn, x_hi));
            // Checked before subtracting, so inf - inf is never formed: constant evaluation rejects arithmetic that
            // produces a NaN (and on GCC an overflow), so the failure would be a compile error instead of an error
            // value.
            if (!math::isfinite(f_lo) || !math::isfinite(f_hi)) return std::unexpected(error::non_finite);
            T const d = (f_hi - f_lo) / dx;
            if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the difference or quotient overflowed
            return Difference<T>{ d, f_lo, f_hi, dx };
        }

        /**
         * @brief The derivative of fn at x by stencil S with step h.
         *
         * By Taylor's theorem, sum_k w_k f(x + n_k h) = sum_j M_j h^j f^(j)(x) / j!, with the moments M_j =
         * sum_k w_k n_k^j. The weights make M_0 = 0, M_1 = denominator and M_j = 0 for 2 <= j <= order, so the sum
         * divided by denominator * h is f'(x) plus a truncation error of order h^order. For Central<2>, the even
         * terms cancel: (f(x + h) - f(x - h)) / 2h = f'(x) + h^2/6 f'''(x) + O(h^4).
         *
         * @tparam S The stencil.
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point; a NaN x fails the overflow check.
         * @param h The step, > 0 and finite.
         * @return The estimate, or error::non_finite.
         */
        template<typename S, real T, typename F>
        constexpr std::expected<T, error> stencil_derivative(F &fn, T const &x, T const &h)
        {
            constexpr auto table = detail::make_table(S{});
            constexpr std::size_t size = table.nodes.size();
            constexpr int lo = table.nodes[0];
            constexpr int hi = table.nodes[size - 1];
            constexpr int k_max = -lo > hi ? -lo : hi;
            constexpr int q = 2 * k_max > table.denominator ? 2 * k_max : table.denominator;
            if (!detail::fits(x, h, q)) return std::unexpected(error::non_finite);

            // x + k h; node 0 is x itself.
            auto const point = [&x, &h](int k) -> T {
                if (k == 0) return x;
                return T(x + T(k) * h);
            };

            if constexpr (size == 2) {
                // With two nodes, the weights are -1 and +1 over hi - lo, so the difference is the whole formula.
                auto const d = detail::difference(fn, point(lo), point(hi));
                if (!d) return std::unexpected(d.error());
                return d->value;
            } else {
                // With more nodes, dividing by the outer-node distance would cancel the outer-node rounding only for
                // two-node rules and enlarge it for the others (4x for Central<4>), so the divisor is denominator * h,
                // in a single division: with dyadic x and a power-of-two h the result is then exact for polynomials of
                // degree up to the order. Node rounding adds about S / 2 * eps * |x f'| / h, the same order as fn's
                // rounding. h is below the resolution of x. Both sides are checked: at a power of the radix the
                // spacing differs on the two sides of x, so h can move x one way and not the other (at x = -1,
                // h = 0x1.8p-54 moves x up but not down), which would leave Backward and Central rules with coinciding
                // points and a wrong result.
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
         * @brief An array of N copies of v, so that T need not be default-constructible.
         *
         * @tparam N The size.
         * @tparam T The real type.
         * @param v The value.
         * @return The array.
         */
        template<std::size_t N, real T>
        constexpr std::array<T, N> filled(T const &v)
        {
            return [&v]<std::size_t... I>(std::index_sequence<I...>) { return std::array<T, N>{ ((void)I, v)... }; }(
                       std::make_index_sequence<N>{});
        }

        /**
         * @brief The levels Ridders may attempt in all, restarts included: at most 2 * 32 = 64 evaluations of fn.
         *
         * Only a guard: a result takes 5 to 11 levels for double and float, 10 to 19 for 50 and 100 digits, up to 18
         * near an edge.
         */
        inline constexpr int ridders_levels = 32;

        /**
         * @brief Ridders' extrapolation of the central difference at x, from the initial step h.
         *
         * From Ridders (1982) and Richardson extrapolation with its a posteriori error estimate (Burden et al.,
         * Numerical Analysis; Sidi, Practical Extrapolation Methods, 2003). The constants are our own, chosen with a
         * prototype (docs/PHASE3_API.md): ratio 2, two levels of stall, the rounding floor, the untrusted-result test.
         *
         * Level i is the central difference D_i at h_i = h_start / 2^i. Its error expansion has even powers of h only,
         * so column j of the tableau removes the h^(2j) term, in Burden's difference form
         *
         *     A[i][j] = A[i][j-1] + (A[i][j-1] - A[i-1][j-1]) / (4^j - 1),
         *
         * which cannot overflow below max (the product form (4^j A[i][j-1] - A[i-1][j-1]) / (4^j - 1) can) and is
         * exact when the two parents are equal. Only the previous row is kept.
         *
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point, finite.
         * @param h The initial step, > 0 and finite.
         * @return The best entry with its error estimate; if the levels run out during restarts, the last untrusted
         *         result; error::non_finite if no tableau reached two finite levels.
         */
        template<real T, typename F>
        constexpr std::expected<Estimate<T>, error> ridders(F &fn, T const &x, T const &h)
        {
            std::array<T, ridders_levels> prev = detail::filled<ridders_levels>(T(0));
            std::array<T, ridders_levels> row = prev;
            std::optional<Estimate<T>> fallback; // the last untrusted result, if the levels run out
            T h_start = h;
            int budget = ridders_levels;

            while (budget > 0) {
                // One tableau, from h_start. The best entry and the stall counter belong to it.
                bool have_best = false;
                bool restart = false;
                bool resolution = false;
                T best = T(0);
                T best_error = T(0);
                T best_floor = T(0);
                T diagonal = T(0);
                int stall = 0;
                T h_i = h_start;
                for (std::size_t i = 0; budget > 0; ++i, h_i = T(h_i / T(2))) {
                    // Checked before evaluating, and after the overflow check, so that x +- h_i is only formed when it
                    // is finite: below the resolution of x the points coincide, and a tiny h could underflow to 0.
                    bool const fit = detail::fits(x, h_i, 2);
                    if (!(h_i > T(0)) || (fit && (T(x + h_i) == x || T(x - h_i) == x))) {
                        resolution = true;
                        break;
                    }
                    --budget;
                    std::expected<Difference<T>, error> const level =
                        fit ? detail::difference(fn, T(x - h_i), T(x + h_i))
                            : std::expected<Difference<T>, error>(std::unexpect, error::non_finite);
                    bool finite = level.has_value();
                    if (finite) {
                        row[0] = level->value;
                        T factor = T(1);
                        for (std::size_t j = 1; j <= i; ++j) {
                            factor = T(factor * T(4));
                            row[j] = T(row[j - 1] + (row[j - 1] - prev[j - 1]) / (factor - T(1)));
                            finite = finite && math::isfinite(row[j]);
                        }
                    }
                    // fn is not finite within h_i of x (or finite values overflowed): every earlier level's interval
                    // contained the bad point, so the whole tableau goes, and a new one starts below the failing step.
                    if (!finite) {
                        restart = true;
                        h_start = T(h_i / T(2));
                        break;
                    }

                    // The estimate of an extrapolated entry is its change from its diagonal parent (the textbook
                    // successive-approximation estimate), but never less than the rounding floor: the minimum of noisy
                    // differences is biased low. Without the floor, the estimate was below the actual error in 58% of
                    // the prototype's random cases. An entry equal to its parent ends the tableau: nothing is left to
                    // learn (a derivative that is exactly 0 then costs 4 to 6 evaluations instead of 64).
                    T const level_floor = detail::rounding_floor(*level);
                    bool improved = false;
                    bool exact = false;
                    for (std::size_t j = 1; j <= i; ++j) {
                        T const change = math::abs(T(row[j] - prev[j - 1]));
                        if (change == T(0)) {
                            best = row[j];
                            best_error = level_floor;
                            best_floor = level_floor;
                            have_best = exact = true;
                            break;
                        }
                        T const estimate = change < level_floor ? level_floor : change;
                        if (!have_best || estimate < best_error) {
                            best = row[j];
                            best_error = estimate;
                            best_floor = level_floor;
                            have_best = improved = true;
                        }
                    }
                    diagonal = row[i];
                    std::swap(prev, row);
                    if (exact) break;
                    // Each halving doubles the rounding error of D_i, so once two levels improve nothing, more levels
                    // only add rounding. One level is not enough: in the prototype a Runge-function case then stopped
                    // early with an error of 3.7e-5.
                    if (i > 0) {
                        stall = improved ? 0 : stall + 1;
                        if (stall == 2) break;
                    }
                }
                if (restart) continue;
                if (!have_best) break; // the resolution of x or the level cap came before two successful levels

                // The best entry must also agree with the most recent diagonal entry, which catches one that only
                // looked good (aliasing). This cut the prototype's uncovered cases from 3.9% to under 1%.
                T const disagreement = math::abs(T(best - diagonal));
                T const abs_error = best_error < disagreement ? disagreement : best_error;
                Estimate<T> const result{ best, abs_error };

                // At most two correct bits, and more error than rounding explains: fn misbehaves within h_start without
                // producing a non-finite value (a pole between the samples, aliasing), so start again below it. A
                // derivative that is 0 has an estimate near the floor and is not affected. In the prototype, the factor
                // 1 instead of 4 left 11 of 400 poles 1/(x - a) uncovered, 4 and 16 none; 16 to 256 for the floor
                // gave identical results.
                // Written with divisions, which cannot overflow: the floor, and so the estimate, saturates at max.
                bool const untrusted = abs_error >= T(math::abs(best) / T(4)) && T(abs_error / T(64)) > best_floor;
                if (!untrusted || resolution) return result;
                fallback = result;
                h_start = T(h_start / T(2));
            }
            if (fallback) return *fallback;
            return std::unexpected(error::non_finite);
        }

    } // namespace detail
    /** @endcond */

    /**
     * @brief The derivative of @p fn at @p x by method @p M with step @p h.
     *
     * fn is called once per node of a stencil, and only at finite points; Ridders calls it 2 times per level, at most
     * 64 times. Failures are detected before the arithmetic that would produce them, so they are error values in
     * constant expressions too, except an overflow of a final sum or quotient, which GCC rejects at compile time.
     * Exceptions thrown by fn or by T pass through unchanged.
     *
     * @tparam M The method: `Central<P>` (by default `Central<2>`), `Forward<P>`, `Backward<P>` or `Ridders`.
     * @tparam T The real type, deduced from @p x.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @param fn The function to differentiate.
     * @param x The point.
     * @param h The step of a stencil, or the initial step of Ridders; it converts to T, so a double step works for a
     *          multiprecision T.
     * @pre h > 0 and h finite; for a stencil, h well above ulp(x).
     * @return The estimate of f'(x), with an error estimate for Ridders, or error::non_finite.
     */
    // h is a non-deduced context (type_identity_t), so T comes from x alone and a double step converts to a
    // multiprecision T instead of making deduction fail.
    template<typename M = Central<2>, real T, typename F>
    requires detail::method<M> && detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Estimate<T>, error>
        derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
    {
        assert(h > T(0) && math::isfinite(h)); // a bad step is a bug in the caller, not a run-time failure
        if constexpr (std::same_as<M, Ridders>) {
            if (!math::isfinite(x)) return std::unexpected(error::non_finite);
            return detail::ridders(fn, x, h);
        } else {
            return detail::stencil_derivative<M>(fn, x, h).transform(
                [](T const &d) { return Estimate<T>{ d, std::nullopt }; });
        }
    }

    /**
     * @brief The derivative of @p fn at @p x by method @p M with its default step.
     *
     * For a stencil of order P the step is eps^(1/(P+1)) * max(|x|, 1). It balances truncation against rounding when
     * fn varies on the scale max(|x|, 1), and the error is then about eps^(P/(P+1)); pass a step when fn varies
     * faster, or when the derivative is much smaller than that error. For Ridders the initial step is
     * max(|x|, 1) / 32, and the method adapts it.
     *
     * @tparam M The method: `Central<P>` (by default `Central<2>`), `Forward<P>`, `Backward<P>` or `Ridders`.
     * @tparam T The real type, deduced from @p x.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @param fn The function to differentiate.
     * @param x The point.
     * @return The estimate of f'(x), with an error estimate for Ridders, or error::non_finite.
     */
    template<typename M = Central<2>, real T, typename F>
    requires detail::method<M> && detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x)
    {
        // Checked first: a NaN x would give a NaN step and trip the precondition of the overload below.
        if (!math::isfinite(x)) return std::unexpected(error::non_finite);
        // Qualified, so that argument-dependent lookup cannot pick a derivative() from the namespace of T or F.
        if constexpr (std::same_as<M, Ridders>) {
            // max(|x|, 1) / 32: in the prototype s/8 and s/32 were equally accurate and s/32 needed fewer
            // evaluations; s/2 was fooled by sin(20x) and s/128 lost a digit. A power of two, so the steps are exact.
            return nxx::deriv::derivative<M>(fn, x, T(detail::scale(x) / T(32)));
        } else {
            return nxx::deriv::derivative<M>(fn, x, detail::default_step<detail::make_table(M{}).order>(x));
        }
    }

} // namespace nxx::deriv
