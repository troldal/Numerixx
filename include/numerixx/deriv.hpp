/**
 * @file
 * @brief Numerical derivatives of a scalar function, of order 1 to 4, by finite differences.
 *
 * derivative<M, N>() estimates the N-th derivative of fn at x from values of fn near x, by a method M: a stencil,
 * `Central<P>`, or `Forward<P>` and `Backward<P>` where fn is defined on one side of x only; or `Ridders`, which adapts
 * the step and also estimates the error. N is 1 by default, an integer or an `order`. Every method returns a
 * std::expected<Estimate<T>, error>; Numerixx itself never throws.
 *
 * @code
 * auto const d = nxx::deriv::derivative<nxx::deriv::Ridders>([](double x) { return std::exp(x); }, 1.0);
 * if (d) use(d->value, *d->abs_error); // e to about 1e-15, with an error estimate of about 6e-14
 * auto const d2 = nxx::deriv::derivative<nxx::deriv::Ridders, 2>([](double x) { return std::exp(x); }, 1.0);
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
         * No finite estimate: x is not finite, fn returned an infinite or NaN value, a difference, term, sum or
         * quotient would overflow, a point x + k h would overflow, or h is below the resolution of x.
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
        /** @brief The estimate of the derivative. */
        T value;

        /** @brief An estimate of |value - exact|, not a guaranteed bound; empty for the stencils. */
        std::optional<T> abs_error;

        /**
         * @brief The value, so that an Estimate can be used where a T is expected.
         * @return value.
         */
        constexpr operator T() const { return value; }

        /** @brief Compares both members; without it, == would compare only the values, through the conversion. */
        friend constexpr bool operator==(Estimate const &, Estimate const &) = default;
    };

    /**
     * @brief The order of a derivative, for readability: `derivative<Ridders, order::second>` is
     * `derivative<Ridders, 2>`.
     *
     * Scoped, so its names stay out of nxx::deriv and it does not convert to int elsewhere; derivative() takes it or
     * an integer.
     */
    enum class order : int {
        first = 1, /**< f' */
        second = 2, /**< f'' */
        third = 3, /**< f''' */
        fourth = 4 /**< f'''' */
    };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief Whether `Central<p>` exists: p even, 2 <= p <= 10.
         *
         * Beyond the cap the gain is under half a digit even for 50-digit types.
         *
         * @param p The accuracy order.
         * @return true if the order is accepted.
         */
        constexpr bool valid_central(int p) { return p >= 2 && p <= 10 && p % 2 == 0; }

        /**
         * @brief Whether `Forward<p>` and `Backward<p>` exist: 1 <= p <= 6.
         *
         * Beyond the cap the rounding amplification, the sum of the absolute weights, is too large (27.7 for p = 6
         * and the first derivative, 4019 for the fourth).
         *
         * @param p The accuracy order.
         * @return true if the order is accepted.
         */
        constexpr bool valid_one_sided(int p) { return p >= 1 && p <= 6; }

        /**
         * @brief The derivative orders derivative() accepts: an integer (not bool) or an `order`, from 1 to 4.
         *
         * Any integer type, so that generic code can pass an int or a std::size_t from an index sequence. An integer
         * is compared in its own type, so a value beyond int's range is not wrapped into 1 to 4.
         *
         * @tparam N The order.
         */
        template<auto N>
        concept valid_order =
            (std::integral<decltype(N)> && !std::same_as<decltype(N), bool> && N >= 1 && N <= 4)
            || (std::same_as<decltype(N), order> && static_cast<int>(N) >= 1 && static_cast<int>(N) <= 4);

        /**
         * @brief Whether T holds every integer up to v exactly: radix^digits > v, decided with integers.
         *
         * The tables are integers, exact in every type with at least 24 significant bits; a narrower type (half
         * precision) would round the larger weights, and the rule would lose its order.
         *
         * @tparam T The real type.
         * @param v The largest integer, >= 0.
         * @return true if T holds 0 ... v exactly.
         */
        template<typename T>
        constexpr bool holds_integers(long long v)
        {
            long long power = 1;
            for (int i = 0; i < std::numeric_limits<T>::digits; ++i) {
                power *= std::numeric_limits<T>::radix;
                if (power > v) return true;
            }
            return false;
        }

        /**
         * @brief Whether T's epsilon is coarser than 1e-7, decided with integers: radix^(digits - 1) < 10^7.
         *
         * float (2^23), decimal32 (10^6) and bfloat16 are coarse, double is not. Integers, because epsilon() need
         * not be constexpr for a user type, and `digits10` would not count binary and decimal types alike.
         *
         * @tparam T The real type.
         * @return true for a coarse type.
         */
        template<typename T>
        constexpr bool coarse()
        {
            long long power = 1;
            for (int i = 1; i < std::numeric_limits<T>::digits; ++i) {
                power *= std::numeric_limits<T>::radix;
                if (power >= 10'000'000) return false;
            }
            return true;
        }

#ifdef NUMERIXX_NO_PRECISION_WARNING
        inline constexpr bool precision_warning = false;
#else
        /**
         * @brief Whether the low-precision warning is on; defining NUMERIXX_NO_PRECISION_WARNING turns it off.
         *
         * The macro must be defined the same way in every translation unit, best on the command line: overload
         * resolution of derivative() would otherwise differ between them.
         */
        inline constexpr bool precision_warning = true;
#endif

        /**
         * @brief A 4th derivative in a coarse type: mostly rounding noise (`Central<2>` about 8e-2 for float, and
         * Ridders a median of 2e-3 with a tenth of the cases more than 100% off).
         *
         * Selects the deprecated overloads of derivative(), whose warning is reported at the caller's line, so a
         * pragma there silences it, and a system-header include does not hide it.
         *
         * @tparam N The order.
         * @tparam T The real type.
         */
        template<auto N, typename T>
        concept noisy = precision_warning && static_cast<int>(N) == 4 && detail::coarse<T>();

    } // namespace detail
    /** @endcond */

    /**
     * @brief Central stencil of accuracy order @p P, the default stencil of derivative().
     *
     * P is the accuracy order, the O(h^P) of the error, not the number of points. For the N-th derivative the nodes
     * are the P + N - 1 integers around 0, without 0 for odd N and with 0 for even N: P + N - 1 evaluations of fn,
     * and an error of about eps^(P/(P+N)) with the default step. `Central<2>` is (fn(x + h) - fn(x - h)) / 2h for f'
     * and (fn(x + h) - 2 fn(x) + fn(x - h)) / h^2 for f''; for f''' it uses x +- h and x +- 2h, `Central<4>` adds
     * x +- 3h (Chapra's O(h^2) and O(h^4) formulas).
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
     * For the N-th derivative, nodes 0 ... P + N - 1: P + N evaluations of fn, including fn(x). For fn defined on one
     * side of x only: the error is larger than for a Central stencil of the same order, and grows quickly with N
     * (`Forward<1>` reaches about eps^(1/2) for f', eps^(1/5) for f'''').
     *
     * @tparam P The accuracy order: 1 <= P <= 6.
     */
    template<int P>
    struct Forward
    {
        static_assert(detail::valid_one_sided(P), "Forward<P> needs an order P from 1 to 6");
    };

    /**
     * @brief `Forward<P>` mirrored: nodes -(P + N - 1) ... 0, so fn is never evaluated above x.
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
     * The `Central<2>` rule of the requested order, at h, h/2, h/4, ..., is extrapolated in Richardson's tableau until
     * rounding stops the improvement, so the step need not be tuned: the initial step should be larger than the best
     * step of that rule, and by default it is max(|x|, 1) / 32 for f' and max(|x|, 1) / 4 for f'' to f''''. The
     * method stops when an extrapolated entry equals its parent, after two levels that improve nothing, or at the
     * resolution of x. It restarts below a step at which fn is not finite, and below one whose result has at most two
     * correct bits by its own estimate while that estimate is far above rounding (a pole between the samples,
     * aliasing).
     *
     * It costs 2 evaluations of fn per level for f' and f'', 4 for f''' and f'''', plus fn(x) once for f'' and f'''':
     * typically 10 to 16 for f' in double, 13 to 23 for f'', 25 to 45 for f''' and f'''', and at most 64, 65, 128 and
     * 129. The result's abs_error is the change of the best entry from its parent, but at least the rounding error of
     * the difference it came from and at least its distance to the latest diagonal entry. It assumes that fn is
     * accurate to about epsilon, and is optimistic for a noisier fn. Like the stencils, it runs in constant
     * expressions.
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

    namespace detail {

        /**
         * @brief Whether S is a stencil: false except for the three families with an accepted order.
         *
         * @tparam S The type.
         */
        template<typename S>
        inline constexpr bool is_stencil = false;

        /** @brief `Central<P>` is a stencil if P is accepted. @tparam P The accuracy order. */
        template<int P>
        inline constexpr bool is_stencil<Central<P>> = valid_central(P);

        /** @brief `Forward<P>` is a stencil if P is accepted. @tparam P The accuracy order. */
        template<int P>
        inline constexpr bool is_stencil<Forward<P>> = valid_one_sided(P);

        /** @brief `Backward<P>` is a stencil if P is accepted. @tparam P The accuracy order. */
        template<int P>
        inline constexpr bool is_stencil<Backward<P>> = valid_one_sided(P);

        /**
         * @brief The stencils: `Central<P>`, `Forward<P>` and `Backward<P>` with an accepted order.
         *
         * @tparam S The type.
         */
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

        /**
         * @brief The weights of a stencil for the N-th derivative:
         * f^(N)(x) ~= sum_k numerators[k] * fn(x + nodes[k] * h) / (denominator * h^N).
         *
         * The error is O(h^order). The weights are integers over one denominator, so they are exact in every real
         * type and the evaluator divides once by denominator * h, then N - 1 times by h.
         *
         * @tparam Size The number of nodes.
         */
        template<std::size_t Size>
        struct Table
        {
            /** @brief The nodes, in units of h, strictly increasing. */
            std::array<int, Size> nodes;

            /** @brief The weight of each node, times denominator. */
            std::array<int, Size> numerators;

            /** @brief The common denominator of the weights, > 0. */
            int denominator;

            /** @brief The accuracy order. */
            int order;
        };

        /**
         * @brief n!, exact for n <= 20; the tables need at most 4!.
         *
         * @param n The argument.
         * @return n!.
         */
        constexpr long long factorial(std::size_t n)
        {
            long long f = 1;
            for (std::size_t i = 2; i <= n; ++i) f *= static_cast<long long>(i);
            return f;
        }

        /**
         * @brief Reached only when a table value needs more than 24 bits: as a call to a function that is not
         * constexpr, it makes the generator's constant evaluation fail, so a weight is never narrowed silently.
         */
        inline void table_value_too_large() {}

        /**
         * @brief The table of the rule of accuracy order p for the n-th derivative on the given nodes.
         *
         * The weights of a rule are fixed by its nodes: by the Lagrange formula (Burden et al., Numerical Analysis,
         * numerical differentiation),
         *
         *     w_k = n! [t^n] prod_{j != k} (t - n_j) / (n_k - n_j),
         *
         * the n-th derivative at 0 of the polynomial that interpolates fn at the nodes. Equivalently, they solve the
         * moment conditions M_j = sum_k w_k n_k^j = 0 for j != n up to n + p - 1, and M_n = n!; Chapra derives the
         * same rules from Taylor series. Each weight is computed in exact rational arithmetic (long long, reduced
         * with gcd), then all are put over their lcm, so that the table holds integers only. The largest
         * intermediate is 12! times 4! (`Central<10>`, n = 4), far inside long long, and every value must fit 24
         * bits, which keeps it exact in float (the largest numerator is 6,222,216).
         *
         * Worked example, `Central<4>` for f' on the nodes -2, -1, 1, 2: for the node -1, prod_{j != k} (t - n_j) =
         * (t + 2)(t - 1)(t - 2) = t^3 - t^2 - 4t + 4, whose t^1 coefficient is -4, and prod (n_k - n_j) =
         * (-1 + 2)(-1 - 1)(-1 - 2) = 6, so w = 1! * -4 / 6 = -2/3. The four weights 1/12, -2/3, 2/3, -1/12 over their
         * lcm 12 give numerators {1, -8, 8, -1}, denominator 12. Chapra's backward f'''' of O(h^2) on -5 ... 0 comes
         * out as {-2, 11, -24, 26, -14, 3}, denominator 1.
         *
         * @tparam Size The number of nodes, > n.
         * @param nodes The nodes, distinct and increasing.
         * @param n The derivative order.
         * @param p The accuracy order, stored in the table.
         * @return The table.
         */
        template<std::size_t Size>
        consteval Table<Size> lagrange_table(std::array<int, Size> const &nodes, int n, int p)
        {
            std::array<long long, Size> num{};
            std::array<long long, Size> den{};
            long long common = 1;
            for (std::size_t k = 0; k < Size; ++k) {
                // poly = prod_{j != k} (t - n_j), coefficients from t^0 up, built one factor at a time;
                // q = prod_{j != k} (n_k - n_j), nonzero for distinct nodes.
                std::array<long long, Size> poly{};
                poly[0] = 1;
                std::size_t degree = 0;
                long long q = 1;
                for (std::size_t j = 0; j < Size; ++j) {
                    if (j == k) continue;
                    ++degree;
                    for (std::size_t i = degree; i > 0; --i) poly[i] = poly[i - 1] - nodes[j] * poly[i];
                    poly[0] = -nodes[j] * poly[0];
                    q *= nodes[k] - nodes[j];
                }
                long long w = detail::factorial(static_cast<std::size_t>(n)) * poly[static_cast<std::size_t>(n)];
                long long const g = std::gcd(w, q);
                w /= g;
                q /= g;
                if (q < 0) {
                    w = -w;
                    q = -q;
                }
                num[k] = w;
                den[k] = q;
                common = std::lcm(common, q);
            }
            constexpr long long limit = 1LL << 24;
            if (common >= limit) detail::table_value_too_large();
            Table<Size> t{};
            for (std::size_t k = 0; k < Size; ++k) {
                long long const v = num[k] * (common / den[k]); // exact: den[k] divides common
                if (v >= limit || v <= -limit) detail::table_value_too_large();
                t.nodes[k] = nodes[k];
                t.numerators[k] = static_cast<int>(v);
            }
            t.denominator = static_cast<int>(common);
            t.order = p;
            return t;
        }

        /**
         * @brief The table of `Central<P>` for the N-th derivative: the P + N - 1 nodes around 0, without 0 for odd
         * N (an odd derivative gives node 0 no weight) and with 0 for even N.
         *
         * @tparam N The derivative order.
         * @tparam P The accuracy order.
         * @return The table.
         */
        template<int N, int P>
        consteval Table<std::size_t{ P + N - 1 }> make_table(Central<P>)
        {
            constexpr int count = P + N - 1;
            constexpr int m = count / 2; // count is even for odd N (P is even), odd for even N
            std::array<int, std::size_t{ count }> nodes{};
            std::size_t i = 0;
            for (int k = -m; k <= m; ++k)
                if (k != 0 || N % 2 == 0) nodes[i++] = k;
            return detail::lagrange_table(nodes, N, P);
        }

        /**
         * @brief The table of `Forward<P>` for the N-th derivative, on the nodes 0 ... P + N - 1.
         *
         * @tparam N The derivative order.
         * @tparam P The accuracy order.
         * @return The table.
         */
        template<int N, int P>
        consteval Table<std::size_t{ P + N }> make_table(Forward<P>)
        {
            std::array<int, std::size_t{ P + N }> nodes{};
            for (std::size_t i = 0; i < nodes.size(); ++i) nodes[i] = static_cast<int>(i);
            return detail::lagrange_table(nodes, N, P);
        }

        /**
         * @brief The table of `Backward<P>` for the N-th derivative, on the nodes -(P + N - 1) ... 0.
         *
         * The generator handles the mirroring: the weights are (-1)^N times those of `Forward<P>` read backwards.
         *
         * @tparam N The derivative order.
         * @tparam P The accuracy order.
         * @return The table.
         */
        template<int N, int P>
        consteval Table<std::size_t{ P + N }> make_table(Backward<P>)
        {
            std::array<int, std::size_t{ P + N }> nodes{};
            for (std::size_t i = 0; i < nodes.size(); ++i) nodes[i] = static_cast<int>(i) - (P + N - 1);
            return detail::lagrange_table(nodes, N, P);
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

        /**
         * @brief The default step of a stencil of order p for the N-th derivative: eps^(1/(p+N)) * max(|x|, 1).
         *
         * A rule of order p has truncation error O(h^p) and rounding error O(eps / h^N), and this step balances the
         * two, scaled so that the step is relative for large |x|.
         *
         * In detail, the error is about C h^p |f^(p+N)| + S eps |f| / h^N, where C is the stencil's error constant
         * and S the sum of its absolute weights. It is smallest at h = (N S eps |f| / (p C |f^(p+N)|))^(1/(p+N)),
         * with an error of order eps^(p/(p+N)). f and its derivatives are unknown, so the step assumes that f varies
         * on the scale max(|x|, 1) and drops the constant (N S / (p C))^(1/(p+N)). For f' that constant lies between
         * 1.4 and 2.0 and costs at most a factor 1.64 in error; for higher N it costs more by the model (3.5 for
         * `Central<2>` and f'', 11 for f''''), but test functions that vary faster than max(|x|, 1) did better with
         * smaller steps for every order, so the convention is kept. For `Central<2>`: h = eps^(1/3) * max(|x|, 1)
         * for f', eps^(1/4) * max(|x|, 1) for f''.
         *
         * @tparam Exponent p + N.
         * @tparam T The real type.
         * @param x The point, finite.
         * @return The step.
         */
        template<int Exponent, real T>
        constexpr T default_step(T const &x)
        { return math::root_epsilon<Exponent, T>() * detail::scale(x); }

        /**
         * @brief Whether every point x + k h and the divisor of a rule are finite: h <= max / 2q and |x| <= max - q h.
         *
         * With K = max(|lo|, |hi|) the largest node and q = max(2K, denominator), nothing overflows before fn is
         * called. A point lies within K h of x, so it stays at least about K h below max: max - q h is rounded by at
         * most half an ulp of max, and K h of margin keeps a rounded point from tying up to infinity, as it could with
         * a margin of 0 (q = K). The divisor is at most q h <= max / 2. The check cannot overflow itself: q h <= max /
         * 2, and max - q h only subtracts. For Central<2> and f', and for Ridders' levels of f' and f'', q = 2; for
         * Ridders' levels of f''' and f'''', which reach x +- 2h, q = 4.
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
         * @brief Whether a + b is finite, decided without forming it.
         *
         * Dividing by the radix is exact, so a + b exceeds max exactly when a / radix + b / radix exceeds max / radix,
         * and the scaled sum cannot overflow. Finite values of opposite sign can differ by more than max, and
         * constant evaluation would reject the overflow (or the NaN of inf - inf that follows).
         *
         * @tparam T The real type.
         * @param a A finite value.
         * @param b A finite value.
         * @return true if a + b is finite.
         */
        template<real T>
        constexpr bool sum_fits(T const &a, T const &b)
        {
            T const radix = static_cast<T>(std::numeric_limits<T>::radix);
            T const largest = (std::numeric_limits<T>::max)();
            return math::abs(T(a / radix + b / radix)) <= T(largest / radix);
        }

        /**
         * @brief |a - b|, saturated at max where it would overflow, like the rounding floors.
         *
         * @tparam T The real type.
         * @param a A finite value.
         * @param b A finite value.
         * @return |a - b|, at most max.
         */
        template<real T>
        constexpr T distance(T const &a, T const &b)
        {
            if (!detail::sum_fits(a, T(T(0) - b))) return (std::numeric_limits<T>::max)();
            return math::abs(T(a - b));
        }

        /**
         * @brief a / b for a >= 0 and b > 0, saturated at max where the division would overflow (a large value over
         * a tiny step), which constant evaluation would reject.
         *
         * @tparam T The real type.
         * @param a The dividend, >= 0 and finite.
         * @param b The divisor, > 0 and finite.
         * @return a / b, at most max.
         */
        template<real T>
        constexpr T saturating_quotient(T const &a, T const &b)
        {
            T const largest = (std::numeric_limits<T>::max)();
            if (T(a / largest) >= T(b / T(2))) return largest;
            return a / b;
        }

        /**
         * @brief Whether a / b is finite, for b > 0, decided without forming it: |a| / max < b / 2.
         *
         * Conservative: it also rejects quotients from max / 2 to max. For the intermediate quotients of the second
         * difference, which would otherwise overflow in constant evaluation (GCC rejects it) before their check.
         *
         * @tparam T The real type.
         * @param a The dividend, finite.
         * @param b The divisor, > 0 and finite.
         * @return true if a / b is below max / 2 in magnitude.
         */
        template<real T>
        constexpr bool quotient_fits(T const &a, T const &b)
        {
            T const largest = (std::numeric_limits<T>::max)();
            return T(math::abs(a) / largest) < T(b / T(2));
        }

        /**
         * @brief a + b for a, b >= 0, saturated at max.
         *
         * @tparam T The real type.
         * @param a A value, >= 0 and finite.
         * @param b A value, >= 0 and finite.
         * @return a + b, at most max.
         */
        template<real T>
        constexpr T saturating_sum(T const &a, T const &b)
        {
            if (!detail::sum_fits(a, b)) return (std::numeric_limits<T>::max)();
            return T(a + b);
        }

        /**
         * @brief a * b for a, b >= 0, saturated at max.
         *
         * @tparam T The real type.
         * @param a A value, >= 0 and finite.
         * @param b A value, >= 0 and finite.
         * @return a * b, at most max.
         */
        template<real T>
        constexpr T saturating_product(T const &a, T const &b)
        {
            T const largest = (std::numeric_limits<T>::max)();
            // >=, since max / b can round up (see weighted_sum).
            if (b > T(1) && a >= T(largest / b)) return largest;
            return T(a * b);
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
         * max.
         *
         * @tparam T The real type.
         * @param d The difference.
         * @return The floor, at most max.
         */
        template<real T>
        constexpr T rounding_floor(Difference<T> const &d)
        {
            T const eps = std::numeric_limits<T>::epsilon();
            return detail::saturating_quotient(T(eps * math::abs(d.f_lo) + eps * math::abs(d.f_hi)), d.dx);
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
         *         not finite at one of them, or the difference or the quotient overflows.
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
            // Checked before subtracting, so neither inf - inf nor an overflowing difference is formed: constant
            // evaluation rejects arithmetic that produces a NaN (and on GCC an overflow), so the failure would be a
            // compile error instead of an error value.
            if (!math::isfinite(f_lo) || !math::isfinite(f_hi)) return std::unexpected(error::non_finite);
            if (!detail::sum_fits(f_hi, T(T(0) - f_lo))) return std::unexpected(error::non_finite);
            T const d = (f_hi - f_lo) / dx;
            if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the quotient overflowed
            return Difference<T>{ d, f_lo, f_hi, dx };
        }

        /**
         * @brief A second difference, with the values and the offsets it came from.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct SecondDifference
        {
            /** @brief The estimate of f''. */
            T value;

            /** @brief fn at the lower point. */
            T f_lo;

            /** @brief fn at the upper point. */
            T f_hi;

            /** @brief x - x_lo as evaluated, > 0. */
            T h_lo;

            /** @brief x_hi - x as evaluated, > 0. */
            T h_hi;
        };

        /**
         * @brief The second difference ((f_hi - f0) / h_hi - (f0 - f_lo) / h_lo) / ((h_hi + h_lo) / 2), from the
         * offsets of the points actually evaluated.
         *
         * Exact for quadratics whatever the rounding of the points, as long as the offsets are exact, which by
         * Sterbenz's lemma holds for 0 < h <= |x| in any radix; the rounding of the points, about eps |x f'| / h^2,
         * then drops out. In the Ridders prototype the nominal step instead gave a worst case of 8.7 times the
         * estimate for float, against 1.4. h^2 is never formed: it underflows for h below about 1.5e-154 in double.
         *
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point, finite.
         * @param f0 fn(x), finite.
         * @param x_lo The lower point, finite and below x.
         * @param x_hi The upper point, finite and above x.
         * @return The second difference with fn's values and the offsets, or error::non_finite if fn is not finite
         *         at a point, or a difference or quotient overflows.
         */
        template<real T, typename F>
        constexpr std::expected<SecondDifference<T>, error>
            second_difference(F &fn, T const &x, T const &f0, T const &x_lo, T const &x_hi)
        {
            T const h_lo = x - x_lo;
            T const h_hi = x_hi - x;
            T const f_lo = static_cast<T>(std::invoke(fn, x_lo));
            T const f_hi = static_cast<T>(std::invoke(fn, x_hi));
            // Every difference is checked before it is formed, as in difference(): two first differences that
            // overflow with the same sign would otherwise give inf - inf.
            if (!math::isfinite(f_lo) || !math::isfinite(f_hi)) return std::unexpected(error::non_finite);
            if (!detail::sum_fits(f_hi, T(T(0) - f0)) || !detail::sum_fits(f0, T(T(0) - f_lo)))
                return std::unexpected(error::non_finite);
            T const up = T(f_hi - f0);
            T const down = T(f0 - f_lo);
            if (!detail::quotient_fits(up, h_hi) || !detail::quotient_fits(down, h_lo))
                return std::unexpected(error::non_finite);
            T const d_hi = up / h_hi;
            T const d_lo = down / h_lo;
            if (!detail::sum_fits(d_hi, T(T(0) - d_lo))) return std::unexpected(error::non_finite);
            T const d = T(d_hi - d_lo) / T(T(h_hi + h_lo) / T(2));
            if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the quotient overflowed
            return SecondDifference<T>{ d, f_lo, f_hi, h_lo, h_hi };
        }

        /**
         * @brief The rounding floor of a second difference: ((eps |f_hi| + eps |f0|) / h_hi + (eps |f0| +
         * eps |f_lo|) / h_lo) / ((h_hi + h_lo) / 2), its first-order sensitivity to fn's rounding.
         *
         * Saturated at every step: one saturation is not enough, since fn = 1e300 over h = 1e-300 overflows the inner
         * division, and fn = 1 over h = 1e-170 only the last.
         *
         * @tparam T The real type.
         * @param d The second difference.
         * @param f0 fn(x).
         * @return The floor, at most max.
         */
        template<real T>
        constexpr T second_floor(SecondDifference<T> const &d, T const &f0)
        {
            T const eps = std::numeric_limits<T>::epsilon();
            T const e0 = eps * math::abs(f0);
            T const hi = detail::saturating_quotient(T(eps * math::abs(d.f_hi) + e0), d.h_hi);
            T const lo = detail::saturating_quotient(T(e0 + eps * math::abs(d.f_lo)), d.h_lo);
            return detail::saturating_quotient(detail::saturating_sum(hi, lo), T(T(d.h_hi + d.h_lo) / T(2)));
        }

        /**
         * @brief A weighted sum of fn's values, and the sum of their rounding errors for Ridders' floor.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct WeightedSum
        {
            /** @brief sum_k numerators[k] * f_k. */
            T value;

            /** @brief sum_k |numerators[k]| eps |f_k|, when requested; 0 otherwise. */
            T magnitude;

            /** @brief fn at node -1, if the rule has one; 0 otherwise. */
            T f_minus;

            /** @brief fn at node +1, if the rule has one; 0 otherwise. */
            T f_plus;
        };

        /**
         * @brief sum_k numerators[k] * fn(x + nodes[k] h), with nothing overflowing before it is checked.
         *
         * Before each term is formed, |f_k| <= max / |numerators[k]| is checked, then the partial sum with sum_fits:
         * finite values can give an infinite term or sum, and the next term of opposite sign then a NaN, which
         * constant evaluation rejects. A constant 1e308 under `Central<4>` used to do exactly that (-8e308, then
         * inf - inf). A sum that would not fit is non_finite; every other sum is formed exactly as before, term by
         * term in node order.
         *
         * @tparam Magnitude Whether to accumulate the magnitude, which only Ridders needs; its bases have numerators
         *                   of at most 6, so for any eps below 1/32 the magnitude cannot overflow.
         * @tparam Size The number of nodes.
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param t The table.
         * @param x The point, finite.
         * @param h The step, > 0, with every point finite.
         * @param f0 fn(x), if the caller already has it; otherwise node 0 evaluates it.
         * @return The sum, or error::non_finite.
         */
        template<bool Magnitude, std::size_t Size, real T, typename F>
        constexpr std::expected<WeightedSum<T>, error>
            weighted_sum(F &fn, Table<Size> const &t, T const &x, T const &h, std::optional<T> const &f0)
        {
            T const largest = (std::numeric_limits<T>::max)();
            T sum = T(0);
            T magnitude = T(0);
            T f_minus = T(0);
            T f_plus = T(0);
            for (std::size_t k = 0; k < Size; ++k) {
                int const node = t.nodes[k];
                int const c = t.numerators[k];
                T fk = T(0);
                if (node == 0 && f0) {
                    fk = *f0;
                } else {
                    T const xk = node == 0 ? x : T(x + T(node) * h); // a T const lvalue, as for two nodes
                    fk = static_cast<T>(std::invoke(fn, xk));
                }
                if (!math::isfinite(fk)) return std::unexpected(error::non_finite);
                if (node == -1) f_minus = fk;
                if (node == 1) f_plus = fk;
                T const size = T(c < 0 ? -c : c);
                // >=, not >: max / |c| can round up, and then |c| times it overflows (6 fl(max / 6) ties to infinity
                // in double).
                if (c != 0 && math::abs(fk) >= T(largest / size)) return std::unexpected(error::non_finite);
                T const term = T(T(c) * fk);
                if (!detail::sum_fits(sum, term)) return std::unexpected(error::non_finite);
                sum = T(sum + term);
                if constexpr (Magnitude) {
                    T const eps = std::numeric_limits<T>::epsilon();
                    magnitude = T(magnitude + size * T(eps * math::abs(fk)));
                }
            }
            return WeightedSum<T>{ sum, magnitude, f_minus, f_plus };
        }

        /**
         * @brief The N-th derivative of fn at x by stencil S with step h.
         *
         * By Taylor's theorem, sum_k w_k f(x + n_k h) = sum_j M_j h^j f^(j)(x) / j!, with the moments M_j =
         * sum_k w_k n_k^j. The weights make M_N = N! denominator and M_j = 0 for every other j up to N + order - 1, so
         * the sum divided by denominator * h^N is f^(N)(x) plus a truncation error of order h^order. For Central<2>
         * and f', the even terms cancel: (f(x + h) - f(x - h)) / 2h = f'(x) + h^2/6 f'''(x) + O(h^4).
         *
         * @tparam S The stencil.
         * @tparam N The derivative order.
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point; a NaN x fails the overflow check.
         * @param h The step, > 0 and finite.
         * @return The estimate, or error::non_finite.
         */
        template<typename S, int N, real T, typename F>
        constexpr std::expected<T, error> stencil_derivative(F &fn, T const &x, T const &h)
        {
            constexpr auto table = detail::make_table<N>(S{});
            constexpr std::size_t size = table.nodes.size();
            constexpr int lo = table.nodes[0];
            constexpr int hi = table.nodes[size - 1];
            constexpr int k_max = -lo > hi ? -lo : hi;
            constexpr int q = 2 * k_max > table.denominator ? 2 * k_max : table.denominator;
            constexpr long long largest_value = [](auto const &t) {
                long long v = t.denominator;
                for (int const c : t.numerators) v = c < -v ? -c : (c > v ? c : v);
                return v;
            }(table);
            static_assert(detail::holds_integers<T>(largest_value),
                "this rule's integer weights are not exact in T: use a type with more significant digits");
            if (!detail::fits(x, h, q)) return std::unexpected(error::non_finite);

            if constexpr (size == 2) {
                // With two nodes (f' only), the weights are -1 and +1 over hi - lo, so the difference is the whole
                // formula.
                auto const point = [&x, &h](int k) -> T {
                    if (k == 0) return x;
                    return T(x + T(k) * h);
                };
                auto const d = detail::difference(fn, point(lo), point(hi));
                if (!d) return std::unexpected(d.error());
                return d->value;
            } else {
                // h is below the resolution of x. Both sides are checked: at a power of the radix the spacing differs
                // on the two sides of x, so h can move x one way and not the other (at x = -1, h = 0x1.8p-54 moves x
                // up but not down), which would leave Backward and Central rules with coinciding points and a wrong
                // result.
                if (T(x + h) == x || T(x - h) == x) return std::unexpected(error::non_finite);
                if constexpr (N == 2 && std::same_as<S, Central<2>>) {
                    // The three-point f'' from the actual offsets, as the two-node f' rules use the actual distance.
                    T const f0 = static_cast<T>(std::invoke(fn, x));
                    if (!math::isfinite(f0)) return std::unexpected(error::non_finite);
                    auto const d = detail::second_difference(fn, x, f0, T(x - h), T(x + h));
                    if (!d) return std::unexpected(d.error());
                    return d->value;
                } else {
                    // With more nodes, dividing by the outer-node distance would cancel the outer-node rounding only
                    // for two-node rules and enlarge it for the others (4x for Central<4>), so the divisor is
                    // denominator * h, then h N - 1 more times, never h^N (h^4 underflows below 1.2e-77 in double):
                    // with dyadic x and a power-of-two h the result is then exact for polynomials of degree up to
                    // N + order - 1. The divisions cannot overflow spuriously: for h >= 1 each quotient is at most
                    // the one before, otherwise the last is the largest. Node rounding adds about
                    // S / 2 * eps * |x f'| / h^N, the same order as fn's rounding. Neighbouring points must differ:
                    // at h near 0.6 ulp(x), x + h and x + 2h round to the same value, and the rule would return a
                    // finite but wrong value.
                    for (std::size_t k = 1; k < size; ++k) {
                        int const a = table.nodes[k - 1];
                        int const b = table.nodes[k];
                        T const x_a = a == 0 ? x : T(x + T(a) * h);
                        T const x_b = b == 0 ? x : T(x + T(b) * h);
                        if (!(x_a < x_b)) return std::unexpected(error::non_finite);
                    }
                    std::optional<T> const none;
                    auto const s = detail::weighted_sum<false>(fn, table, x, h, none);
                    if (!s) return std::unexpected(s.error());
                    T d = s->value / (T(table.denominator) * h);
                    for (int i = 1; i < N; ++i) d = T(d / h);
                    if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the quotient overflowed
                    return d;
                }
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
         * @brief The levels Ridders may attempt in all, restarts included: at most 64 evaluations of fn for f', 65
         * for f'', 128 for f''' and 129 for f''''.
         *
         * Only a guard: a result takes 5 to 11 levels for f' in double and float, 10 to 19 for 50 and 100 digits, up
         * to 18 near an edge; f''' and f'''' take up to 20 levels for 100 digits, and 28 after restarts at a pole.
         */
        inline constexpr int ridders_levels = 32;

        /**
         * @brief The divisor of Ridders' default initial step max(|x|, 1) / divisor, for orders 1 to 4.
         *
         * For f' s/8 and s/32 were equally accurate in the prototype and s/32 needed fewer evaluations; s/2 was
         * fooled by sin(20x) and s/128 lost a digit. The base's rounding error grows as 1/h^N, so a larger start pays
         * for higher orders: for f'' s/4 was 1.6 times as accurate as s/8 for one more evaluation (7 times s/32), for
         * f''' and f'''' 1.6 and 2.6 times. Powers of two, so the steps are exact.
         */
        inline constexpr std::array<int, 4> ridders_divisors{ 32, 4, 4, 4 };

        /**
         * @brief One level of Ridders' tableau: the base difference and its rounding floor.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct Level
        {
            /** @brief The `Central<2>` estimate of the derivative at this step. */
            T value;

            /** @brief Its rounding error when fn is accurate to about eps, at most max. */
            T floor;
        };

        /**
         * @brief One level of Ridders' tableau for the N-th derivative at step h: the `Central<2>` rule of order N.
         *
         * f' and f'' use the actual distance and offsets (difference(), second_difference()): with the nominal step,
         * the prototype's worst case was 22 against 3.6 times the estimate for f', and 8.7 against 1.4 for float f''.
         * f''' and f'''' use the nominal step, (f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)) / 2h^3 and
         * (f(x+2h) - 4f(x+h) + 6f(x) - 4f(x-h) + f(x-2h)) / h^4, with errors h^2 f^(5)/4 and h^2 f^(6)/6 at the lowest
         * order: their expansions are even, as the tableau needs. Their floor has two parts, each saturated at every
         * step: fn's rounding, sum |w_k| eps |f_k| / (denominator h^N), and the rounding of the points, which the
         * nominal step leaves in: each point x + k h is off by up to eps/2 (|x| + 2h), which moves its value by about
         * that times f', so eps/2 (|x| + 2h) |f'| sum |w_k| / (denominator h^N), with f' from the values at x +- h.
         * Without the second part, the prototype's f''' of exp(10x) on [1, 3] was covered in 96.3% of the cases, with
         * it in 99.0%.
         *
         * @tparam N The derivative order.
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point, finite.
         * @param f0 fn(x) for even N; unused for odd N.
         * @param h The step, with every point finite and distinct.
         * @return The level, or error::non_finite.
         */
        template<int N, real T, typename F>
        constexpr std::expected<Level<T>, error> ridders_level(F &fn, T const &x, T const &f0, T const &h)
        {
            if constexpr (N == 1) {
                auto const d = detail::difference(fn, T(x - h), T(x + h));
                if (!d) return std::unexpected(d.error());
                return Level<T>{ d->value, detail::rounding_floor(*d) };
            } else if constexpr (N == 2) {
                auto const d = detail::second_difference(fn, x, f0, T(x - h), T(x + h));
                if (!d) return std::unexpected(d.error());
                return Level<T>{ d->value, detail::second_floor(*d, f0) };
            } else {
                constexpr auto table = detail::make_table<N>(Central<2>{});
                std::optional<T> const cached = f0;
                auto const s = detail::weighted_sum<true>(fn, table, x, h, cached);
                if (!s) return std::unexpected(s.error());
                constexpr int weights = [](auto const &t) {
                    int w = 0;
                    for (int const c : t.numerators) w += c < 0 ? -c : c;
                    return w;
                }(table);
                T const eps = std::numeric_limits<T>::epsilon();
                T const slope = detail::saturating_quotient(detail::distance(s->f_plus, s->f_minus), T(T(2) * h));
                T const spread = T(T(eps / T(2)) * T(math::abs(x) + T(2) * h)) * T(weights);
                T const divisor = T(table.denominator) * h;
                T value = s->value / divisor;
                T rounding = detail::saturating_quotient(s->magnitude, divisor);
                T points = detail::saturating_quotient(detail::saturating_product(spread, slope), divisor);
                for (int i = 1; i < N; ++i) {
                    value = T(value / h);
                    rounding = detail::saturating_quotient(rounding, h);
                    points = detail::saturating_quotient(points, h);
                }
                if (!math::isfinite(value)) return std::unexpected(error::non_finite); // the quotient overflowed
                return Level<T>{ value, detail::saturating_sum(rounding, points) };
            }
        }

        /**
         * @brief Ridders' extrapolation of the `Central<2>` rule of order N at x, from the initial step h.
         *
         * From Ridders (1982) and Richardson extrapolation with its a posteriori error estimate (Burden et al.,
         * Numerical Analysis; Sidi, Practical Extrapolation Methods, 2003). The constants are our own, chosen with a
         * prototype (docs/PHASE3_API.md, docs/PHASE4_API.md): ratio 2, two levels of stall, the rounding floor, the
         * untrusted-result test; the same for every order.
         *
         * Level i is the base rule D_i at h_i = h_start / 2^i. Its error expansion has even powers of h only, so
         * column j of the tableau removes the h^(2j) term, in Burden's difference form
         *
         *     A[i][j] = A[i][j-1] + (A[i][j-1] - A[i-1][j-1]) / (4^j - 1),
         *
         * which is exact when the two parents are equal. Only the previous row is kept.
         *
         * @tparam N The derivative order.
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point, finite.
         * @param h The initial step, > 0 and finite.
         * @return The best entry with its error estimate; if the levels run out during restarts, the last untrusted
         *         result; error::non_finite if fn(x) is not finite (even N) or no tableau reached two finite levels.
         */
        template<int N, real T, typename F>
        constexpr std::expected<Estimate<T>, error> ridders(F &fn, T const &x, T const &h)
        {
            // The levels of f' and f'' reach x +- h, those of f''' and f'''' x +- 2h.
            constexpr int q = N <= 2 ? 2 : 4;
            // fn(x) is part of every level for even N: evaluated once, and kept across restarts. If it is not finite,
            // no step avoids it.
            T const f0 = N % 2 == 0 ? static_cast<T>(std::invoke(fn, x)) : T(0);
            if (!math::isfinite(f0)) return std::unexpected(error::non_finite);

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
                    // Checked before evaluating, and after the overflow check, so that x +- h_i (and x +- 2 h_i) is
                    // only formed when it is finite: below the resolution of x the points coincide, and a tiny h could
                    // underflow to 0. For f''' and f'''' the inner and outer points must differ too (at h near 0.6 ulp,
                    // x + h and x + 2h round to the same value).
                    bool const fit = detail::fits(x, h_i, q);
                    bool const coincide =
                        fit
                        && (T(x + h_i) == x || T(x - h_i) == x
                            || (N >= 3 && (T(x + T(2) * h_i) == T(x + h_i) || T(x - T(2) * h_i) == T(x - h_i))));
                    if (!(h_i > T(0)) || coincide) {
                        resolution = true;
                        break;
                    }
                    --budget;
                    std::expected<Level<T>, error> const level =
                        fit ? detail::ridders_level<N>(fn, x, f0, h_i)
                            : std::expected<Level<T>, error>(std::unexpect, error::non_finite);
                    bool finite = level.has_value();
                    if (finite) {
                        row[0] = level->value;
                        T factor = T(1);
                        // Both the difference and the sum are checked before they are formed: finite entries of
                        // opposite sign can differ by more than max.
                        for (std::size_t j = 1; j <= i && finite; ++j) {
                            factor = T(factor * T(4));
                            finite = detail::sum_fits(row[j - 1], T(T(0) - prev[j - 1]));
                            if (!finite) break;
                            T const correction = T((row[j - 1] - prev[j - 1]) / (factor - T(1)));
                            finite = detail::sum_fits(row[j - 1], correction);
                            if (!finite) break;
                            row[j] = T(row[j - 1] + correction);
                        }
                    }
                    // fn is not finite within h_i of x (or an entry would exceed max): every earlier level's interval
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
                    T const level_floor = level->floor;
                    bool improved = false;
                    bool exact = false;
                    for (std::size_t j = 1; j <= i; ++j) {
                        T const change = detail::distance(row[j], prev[j - 1]);
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
                    // Each halving multiplies the rounding error of D_i by 2^N, so once two levels improve nothing,
                    // more levels only add rounding. One level is not enough: in the prototype a Runge-function case
                    // then stopped early with an error of 3.7e-5.
                    if (i > 0) {
                        stall = improved ? 0 : stall + 1;
                        if (stall == 2) break;
                    }
                }
                if (restart) continue;
                if (!have_best) break; // the resolution of x or the level cap came before two successful levels

                // The best entry must also agree with the most recent diagonal entry, which catches one that only
                // looked good (aliasing). This cut the prototype's uncovered cases from 3.9% to under 1%.
                T const disagreement = detail::distance(best, diagonal);
                T const abs_error = best_error < disagreement ? disagreement : best_error;
                Estimate<T> const result{ best, abs_error };

                // At most two correct bits, and more error than rounding explains: fn misbehaves within h_start without
                // producing a non-finite value (a pole between the samples, aliasing), so start again below it. A
                // derivative that is 0 has an estimate near the floor and is not affected. In the prototype, the factor
                // 1 instead of 4 left 11 of 400 poles 1/(x - a) uncovered for f' (32 for f''), 4 and 16 none; 16 to
                // 256 for the floor gave identical results.
                // Written with divisions, which cannot overflow: the floor, and so the estimate, saturates at max.
                bool const untrusted = abs_error >= T(math::abs(best) / T(4)) && T(abs_error / T(64)) > best_floor;
                if (!untrusted || resolution) return result;
                fallback = result;
                h_start = T(h_start / T(2));
            }
            if (fallback) return *fallback;
            return std::unexpected(error::non_finite);
        }

        /**
         * @brief The N-th derivative of fn at x by method M with step h, behind the public overloads.
         *
         * @tparam M The method.
         * @tparam N The derivative order, 1 to 4.
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point.
         * @param h The step, > 0 and finite.
         * @return The estimate, or error::non_finite.
         */
        template<typename M, int N, real T, typename F>
        constexpr std::expected<Estimate<T>, error> derivative_with_step(F &fn, T const &x, T const &h)
        {
            assert(h > T(0) && math::isfinite(h)); // a bad step is a bug in the caller, not a run-time failure
            if constexpr (std::same_as<M, Ridders>) {
                if (!math::isfinite(x)) return std::unexpected(error::non_finite);
                return detail::ridders<N>(fn, x, h);
            } else {
                return detail::stencil_derivative<M, N>(fn, x, h).transform(
                    [](T const &d) { return Estimate<T>{ d, std::nullopt }; });
            }
        }

        /**
         * @brief The N-th derivative of fn at x by method M with its default step, behind the public overloads.
         *
         * @tparam M The method.
         * @tparam N The derivative order, 1 to 4.
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point.
         * @return The estimate, or error::non_finite.
         */
        template<typename M, int N, real T, typename F>
        constexpr std::expected<Estimate<T>, error> derivative_default(F &fn, T const &x)
        {
            // Checked first: a NaN x would give a NaN step and trip the precondition.
            if (!math::isfinite(x)) return std::unexpected(error::non_finite);
            if constexpr (std::same_as<M, Ridders>) {
                constexpr int divisor = ridders_divisors[static_cast<std::size_t>(N - 1)];
                return detail::derivative_with_step<M, N>(fn, x, T(detail::scale(x) / T(divisor)));
            } else {
                constexpr int exponent = detail::make_table<N>(M{}).order + N;
                return detail::derivative_with_step<M, N>(fn, x, detail::default_step<exponent>(x));
            }
        }

    } // namespace detail

    /**
     * @brief The @p N-th derivative of @p fn at @p x by method @p M with step @p h.
     *
     * fn is called once per node of a stencil, and only at finite points; Ridders calls it 2 or 4 times per level, at
     * most 64 to 129 times. Failures are detected before the arithmetic that would produce them, so they are error
     * values in constant expressions too, except an overflow of a final quotient, which GCC rejects at compile time.
     * Exceptions thrown by fn or by T pass through unchanged.
     *
     * @code
     * auto const d3 = derivative<Central<2>, order::third>(fn, x, h); // f''' on x +- h, x +- 2h: Chapra's O(h^2)
     * auto const d4 = derivative<Backward<2>, 4>(fn, x, h);           // f'''' on x - 5h ... x
     * @endcode
     *
     * @tparam M The method: `Central<P>` (by default `Central<2>`), `Forward<P>`, `Backward<P>` or `Ridders`.
     * @tparam N The derivative order: 1 (the default) to 4, an integer or an `order`.
     * @tparam T The real type, deduced from @p x.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @param fn The function to differentiate.
     * @param x The point.
     * @param h The step of a stencil, or the initial step of Ridders; it converts to T, so a double step works for a
     *          multiprecision T.
     * @pre h > 0 and h finite; for a stencil, h well above ulp(x).
     * @return The estimate of the derivative, with an error estimate for Ridders, or error::non_finite.
     */
    // h is a non-deduced context (type_identity_t), so T comes from x alone and a double step converts to a
    // multiprecision T instead of making deduction fail. N is auto, because the scoped order does not convert to int.
    template<typename M = Central<2>, auto N = 1, real T, typename F>
    requires detail::method<M> && detail::valid_order<N> && detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Estimate<T>, error>
        derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
    {
        // Qualified, so that argument-dependent lookup cannot pick a function from the namespace of T or F.
        return detail::derivative_with_step<M, static_cast<int>(N)>(fn, x, h);
    }

    /**
     * @brief The @p N-th derivative of @p fn at @p x by method @p M with its default step.
     *
     * For a stencil of order P the step is eps^(1/(P+N)) * max(|x|, 1). It balances truncation against rounding when
     * fn varies on the scale max(|x|, 1), and the error is then about eps^(P/(P+N)); pass a step when fn varies
     * faster, or when the derivative is much smaller than that error. The sensitivity grows with N: for f'''',
     * misjudging fn's scale by a factor 2 costs up to 16 in error. For Ridders the initial step is max(|x|, 1) / 32
     * for f' and max(|x|, 1) / 4 for f'' to f'''', and the method adapts it. The default steps differ between methods
     * and orders, so generic code that varies N should use this overload rather than pass one h for every order.
     *
     * @tparam M The method: `Central<P>` (by default `Central<2>`), `Forward<P>`, `Backward<P>` or `Ridders`.
     * @tparam N The derivative order: 1 (the default) to 4, an integer or an `order`.
     * @tparam T The real type, deduced from @p x.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @param fn The function to differentiate.
     * @param x The point.
     * @return The estimate of the derivative, with an error estimate for Ridders, or error::non_finite.
     */
    template<typename M = Central<2>, auto N = 1, real T, typename F>
    requires detail::method<M> && detail::valid_order<N> && detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x)
    { return detail::derivative_default<M, static_cast<int>(N)>(fn, x); }

    /**
     * @brief The overload with a step, for a 4th derivative in a type with an epsilon above 1e-7 (float): the same
     * result, with a deprecation warning, since that result is mostly rounding noise.
     *
     * Chosen over the general overload by its additional constraint. The warning is reported at the caller's line,
     * so a pragma there silences it (`-Wdeprecated-declarations`, C4996 on MSVC); NUMERIXX_NO_PRECISION_WARNING turns
     * it off everywhere.
     *
     * @tparam M The method.
     * @tparam N The derivative order, 4.
     * @tparam T The real type.
     * @tparam F The callable.
     * @param fn The function to differentiate.
     * @param x The point.
     * @param h The step, or the initial step of Ridders.
     * @return As for the general overload.
     */
    template<typename M = Central<2>, auto N = 1, real T, typename F>
    requires detail::method<M> && detail::valid_order<N> && detail::callable<F, T> && detail::noisy<N, T>
    [[deprecated(
        "the 4th derivative in a type with 7 or fewer significant digits is mostly rounding noise: use "
        "double, or silence this warning at the call or with NUMERIXX_NO_PRECISION_WARNING")]] [[nodiscard]]
    constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
    { return detail::derivative_with_step<M, static_cast<int>(N)>(fn, x, h); }

    /**
     * @brief The overload with the default step, for a 4th derivative in a type with an epsilon above 1e-7 (float):
     * the same result, with a deprecation warning, since that result is mostly rounding noise.
     *
     * @tparam M The method.
     * @tparam N The derivative order, 4.
     * @tparam T The real type.
     * @tparam F The callable.
     * @param fn The function to differentiate.
     * @param x The point.
     * @return As for the general overload.
     */
    template<typename M = Central<2>, auto N = 1, real T, typename F>
    requires detail::method<M> && detail::valid_order<N> && detail::callable<F, T> && detail::noisy<N, T>
    [[deprecated(
        "the 4th derivative in a type with 7 or fewer significant digits is mostly rounding noise: use "
        "double, or silence this warning at the call or with NUMERIXX_NO_PRECISION_WARNING")]] [[nodiscard]]
    constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x)
    { return detail::derivative_default<M, static_cast<int>(N)>(fn, x); }

} // namespace nxx::deriv
