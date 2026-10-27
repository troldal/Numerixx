/**
 * @file
 * @brief Numerical derivatives of a scalar function, of order 1 to 4, by finite differences.
 *
 * derivative<M, N>() estimates the N-th derivative of fn at x from values of fn near x, by a method M: a stencil,
 * `Central<P>`, or `Forward<P>` and `Backward<P>` where fn is defined on one side of x only; or `Ridders`, which adapts
 * the step and also estimates the error, and its one-sided forms `ForwardRidders` and `BackwardRidders`. N is 1 by
 * default, an integer or an `order`. Every method returns a
 * std::expected<Estimate<T>, error>; Numerixx itself never throws, but the result's own value() throws
 * std::bad_expected_access on a failure, so test the result and read it with `->` or `*`.
 *
 * @code
 * auto const d = nxx::deriv::derivative<nxx::deriv::Ridders>([](double x) { return std::exp(x); }, 1.0);
 * if (d) use(d->value, *d->abs_error); // e to about 1e-15, with an error estimate of about 6e-14
 * auto const d2 = nxx::deriv::derivative<nxx::deriv::Ridders, 2>([](double x) { return std::exp(x); }, 1.0);
 * @endcode
 *
 * derivative_of<M, N>(fn) gives the same derivative as a callable, for algorithms that take one, such as Newton's
 * method in nxx::roots.
 */

#pragma once

#include <numerixx/math.hpp>

#include <array>
#include <concepts>
#include <cstddef>
#include <expected>
#include <format>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>

/** @brief Numerical derivatives of scalar functions. */
namespace nxx::deriv {

    /**
     * @brief Why derivative() returned no estimate.
     *
     * More values may be added, and a listed cause may move to a new value: a switch on it should have a default, and a
     * caller should test for success rather than for a particular error. It is formattable: std::format("{}", e) writes
     * the enumerator's name, such as `non_finite`.
     */
    enum class error {
        /**
         * No finite estimate: x is not finite, the step is not positive and finite, fn returned an infinite or NaN
         * value, a difference, term, sum or quotient would overflow, a point x + k h would overflow, or h is below the
         * resolution of x: for the two-node rules of f' (`Central<2>`, `Forward<1>`, `Backward<1>`) when their two
         * points coincide, for every other rule when x + h or x - h rounds to x.
         */
        non_finite
    };

    /**
     * @brief A derivative estimate: the value, and an estimate of its absolute error when the method provides one.
     *
     * Converts implicitly to T. Where a conversion is not considered (auto, std::format, and the arithmetic operator
     * templates of Boost.Multiprecision, whose comparisons do accept it), write `e.value` or `T(e)`. The conversion
     * also applies where it may surprise: `if (*r)` tests the value against 0, not whether there is a result (and Clang
     * rejects it under -Wfloat-conversion with -Werror; write `r->value != 0`), and `int n = *r` truncates the value.
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
         * The cap is a design decision, not a measured limit: by the model's exponent alone, `Central<12>` would cut
         * the error of f' by a factor of about 1.7 in double and about 5 in a 50-digit type.
         *
         * @param p The accuracy order.
         * @return true if the order is accepted.
         */
        constexpr bool valid_central(int p) { return p >= 2 && p <= 10 && p % 2 == 0; }

        /**
         * @brief Whether `Forward<p>` and `Backward<p>` exist: 1 <= p <= 6.
         *
         * The cap is a design decision: the rounding amplification, the sum of the absolute weights, is already 27.7
         * for p = 6 and the first derivative (4019 for the fourth) and grows about as fast as the order pays; by the
         * model, `Forward<7>` would lower the best error for f' by about 14 % and `Forward<8>` not at all.
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
         * @brief A derivative above the first in a coarse type: few correct digits, whatever the method.
         *
         * Median relative errors in float, over nine functions at 300 points each: f' 3e-5 with `Central<2>` and 2e-6
         * with `Ridders`; f'' 4e-4 and 1e-5; f''' 5e-3 and 1e-4; f'''' 2e-2 and 1e-3, with a tenth of Ridders' f''''
         * more than 100% off. `ForwardRidders` gives 8e-4 for f'', 2e-2 for f''' and noise for f''''.
         *
         * Selects the deprecated overloads of derivative(), whose warning is reported at the caller's line, so a
         * pragma there silences it, and a system-header include does not hide it.
         *
         * @tparam N The order.
         * @tparam T The real type.
         */
        template<auto N, typename T>
        concept noisy = precision_warning && static_cast<int>(N) >= 2 && detail::coarse<T>();

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
     * resolution of x. It restarts below a step at which fn is not finite, and, after a result that has at most two
     * correct bits by its own estimate while that estimate is far above rounding (a pole near x, aliasing), below the
     * smallest step it tried.
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

    /**
     * @brief Ridders' extrapolation of the `Forward<1>` rule of the requested order: fn is evaluated at x and above
     * only, with an error estimate.
     *
     * For fn that is undefined or different below x: x at a domain edge, where `Ridders` fails (its samples reach
     * below x), the lower end of a correlation's range, or the right-hand derivative at a kink. The rule's error has
     * every power of h, so each column of the tableau gains one order rather than two: it is less accurate than
     * `Ridders`, typically 2e-13 relative for f' in double, 1e-10 for f'', 4e-8 for f''' and 2e-6 for f'''', against
     * 1e-14 to 7e-10. Where fn is defined on both sides, use `Ridders`, also near a domain edge.
     *
     * The tableau, its stopping rules and restarts are those of `Ridders`, with the initial step max(|x|, 1) / 32 for
     * every order; the error estimate also covers the rounding that the tableau amplifies. From the default step the
     * method needs fn to vary on a scale of at least about 2^-30 max(|x|, 1) (sqrt or log at 1e-9 from 0); closer, its
     * result is wrong but flagged by an estimate of at least |value| / 4: pass a step on fn's own scale. It costs fn(x)
     * once and N evaluations per level: typically 10 for f' and 17, 22 and 29 for f'' to f'''' in double, at most 33,
     * 65, 97 and 129.
     *
     * @code
     * // log is NaN below 0.5 here, and x is the edge: Ridders returns non_finite, ForwardRidders the derivative.
     * auto const fn = [](double x) { return x < 0.5 ? std::nan("") : std::log(x); };
     * auto const r = nxx::deriv::derivative<nxx::deriv::ForwardRidders>(fn, 0.5);
     * if (r) use(r->value, *r->abs_error); // 2 to about 1e-13 relative
     * @endcode
     */
    struct ForwardRidders
    {
    };

    /**
     * @brief `ForwardRidders` mirrored: fn is evaluated at x and below only, with an error estimate.
     *
     * For fn that is undefined or different above x, such as the upper end of a correlation's range. It is
     * `ForwardRidders` applied to fn(-t) at -x, which gives the same points x - k h and the same accuracy, cost and
     * limits.
     */
    struct BackwardRidders
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
         * @brief The extrapolation methods: `Ridders`, `ForwardRidders` and `BackwardRidders`.
         *
         * @tparam M The method type.
         */
        template<typename M>
        concept extrapolation =
            std::same_as<M, Ridders> || std::same_as<M, ForwardRidders> || std::same_as<M, BackwardRidders>;

        /**
         * @brief The methods derivative() accepts: the stencils and the extrapolation methods.
         *
         * Closed on purpose: users cannot add a method.
         *
         * @tparam M The method type.
         */
        template<typename M>
        concept method = stencil<M> || extrapolation<M>;

        /**
         * @brief The stencil whose table method M evaluates: M itself for a stencil, `Central<2>` for `Ridders`,
         * `Forward<1>` for the one-sided extrapolations (`BackwardRidders` runs it on the mirrored function).
         *
         * @tparam M The method.
         */
        template<typename M>
        using base = std::
            conditional_t<std::same_as<M, Ridders>, Central<2>, std::conditional_t<extrapolation<M>, Forward<1>, M>>;

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
         * moment conditions M_j = sum_k w_k n_k^j = 0 for j != n up to n + p - 1, and M_n = n!; Chapra derives the same
         * rules from Taylor series. Each weight is computed in exact rational arithmetic (long long, reduced with gcd),
         * then all are put over their lcm, so that the table holds integers only. The largest intermediate is 12!
         * (`Central<10>`, n = 4; 12! times 4! bounds it), far inside long long, and every value must fit 24 bits, which
         * keeps it exact in float (the largest numerator is 6,222,216).
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
         * In detail, the error is about C h^p |f^(p+N)| + S eps |f| / h^N, where C is the stencil's error constant and
         * S the sum of its absolute weights. It is smallest at h = (N S eps |f| / (p C |f^(p+N)|))^(1/(p+N)), with an
         * error of order eps^(p/(p+N)). f and its derivatives are unknown, so the step assumes that f varies on the
         * scale max(|x|, 1) and drops the constant (N S / (p C))^(1/(p+N)). For f' that constant lies between 1.4 and
         * 2.0 and costs at most a factor of about 1.6 in error (1.63 for `Central<10>`); for higher N it costs more by
         * the model (3.5 for `Central<2>` and f'', 11 for f''''), but test functions that vary faster than max(|x|, 1)
         * did better with smaller steps for every order, so the convention is kept. For `Central<2>`: h = eps^(1/3) *
         * max(|x|, 1) for f', eps^(1/4) * max(|x|, 1) for f''.
         *
         * @tparam Exponent p + N.
         * @tparam T The real type.
         * @param x The point, finite.
         * @return The step.
         */
        template<int Exponent, real T>
        constexpr T stencil_step(T const &x)
        { return math::root_epsilon<Exponent, T>() * detail::scale(x); }

        /**
         * @brief The largest integer in the weight table that method M uses for the N-th derivative: the largest
         * numerator in magnitude, or the denominator.
         *
         * An extrapolation method uses the table of its base stencil, as ridders_level does for f''' and f'''' (the
         * f' and f'' levels need no larger integers than that table).
         *
         * @tparam M The method.
         * @tparam N The derivative order, 1 to 4.
         */
        template<typename M, int N>
        inline constexpr long long largest_weight_integer = [] {
            constexpr auto table = detail::make_table<N>(detail::base<M>{});
            long long v = table.denominator;
            for (int const c : table.numerators) v = c < -v ? -c : (c > v ? c : v);
            return v;
        }();

        /**
         * @brief Whether T holds every integer of the weight table of method M for the N-th derivative exactly, so
         * that the rule keeps its order.
         *
         * Every type with 24 significant bits or more holds every table (the largest integer is 6,222,216, of
         * `Central<10>` for f''''); a narrower type, such as half precision, needs a rule with smaller weights: a lower
         * order, or Ridders. A constraint of derivative(), so that a requires-expression sees it; it must follow
         * method<M> and valid_order<N> there, because it builds the table.
         *
         * @tparam M The method.
         * @tparam N The derivative order, 1 to 4, an integer or an `order`.
         * @tparam T The real type.
         */
        template<typename M, auto N, typename T>
        concept exact_weights = detail::holds_integers<T>(detail::largest_weight_integer<M, static_cast<int>(N)>);

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
         * @brief Whether a * b is finite, for b >= 1, decided without forming it.
         *
         * With q = max / b, rounded: a < q gives a b <= max (if q rounded up, the value below it is at most the exact
         * max / b), and a > q a product beyond max. Only a = q is undecided by the comparison, and there the scaled
         * product (q / radix) * b, exact up to its own rounding and near max / radix, cannot overflow: dividing by the
         * radix is exact, as in sum_fits. So nothing is rejected that fits, unlike a plain a >= q, which rejected
         * max * 1, and nothing that overflows is accepted, unlike a > q, which accepted 6 fl(max / 6) (a tie to
         * infinity in double).
         *
         * @tparam T The real type.
         * @param a A value, >= 0 and finite.
         * @param b A value, >= 1 and finite.
         * @return true if a * b is finite.
         */
        template<real T>
        constexpr bool product_fits(T const &a, T const &b)
        {
            T const radix = static_cast<T>(std::numeric_limits<T>::radix);
            T const largest = (std::numeric_limits<T>::max)();
            T const q = largest / b;
            if (a < q) return true;
            if (a > q) return false;
            return T(T(q / radix) * b) <= T(largest / radix);
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
            if (b > T(1) && !detail::product_fits(a, b)) return (std::numeric_limits<T>::max)();
            return T(a * b);
        }

        /**
         * @brief The rounding error of a value of fn that is accurate to about eps: eps |f|, but at least the spacing
         * of the numbers below the normal range, where rounding is absolute; 0 for f = 0.
         *
         * Without the absolute part, a fn whose values are subnormal had a rounding floor of 0, and Ridders reported
         * a result 192 denorm_min off with an estimate of 0. For normal values this is exactly eps |f|.
         *
         * @tparam T The real type.
         * @param f A value of fn, finite.
         * @return The rounding error, >= 0.
         */
        template<real T>
        constexpr T value_rounding(T const &f)
        {
            T const e = std::numeric_limits<T>::epsilon() * math::abs(f);
            T const tiny = T(std::numeric_limits<T>::denorm_min());
            if (f != T(0) && e < tiny) return tiny;
            return e;
        }

        /**
         * @brief a / (dh h^(N-1)) with dh = denominator * h, without forming h^(N-1) and without a spurious
         * underflow or overflow where the order of the divisions can avoid it.
         *
         * Dividing by dh first, then by h, keeps every intermediate between a and the result when dh <= 1 or h >= 1.
         * When dh > 1 > h, that order shrinks the value before the divisions by h grow it back, and the first quotient
         * can underflow to 0 although the result is normal (`Central<10>` for f'''' with h = 19 * 2^-22 returned 0
         * for a result of 2.6e-308). Then the divisions by h come first and dh last: the intermediates exceed the
         * result by at most the factor dh, so only a result within that factor of max can overflow spuriously, which
         * is reported as non_finite rather than as a silent 0. For N = 1 there is a single division either way.
         *
         * @tparam N The derivative order.
         * @tparam Saturating Whether each division saturates at max (for the rounding floors).
         * @tparam T The real type.
         * @param a The dividend.
         * @param dh denominator * h, > 0.
         * @param h The step, > 0.
         * @return The quotient.
         */
        template<int N, bool Saturating, real T>
        constexpr T divide_by_powers(T const &a, T const &dh, T const &h)
        {
            auto const divide = [](T const &p, T const &q) -> T {
                if constexpr (Saturating)
                    return detail::saturating_quotient(p, q);
                else
                    return T(p / q);
            };
            T v = a;
            if (dh > T(1) && h < T(1)) {
                for (int i = 1; i < N; ++i) v = divide(v, h);
                return divide(v, dh);
            }
            v = divide(v, dh);
            for (int i = 1; i < N; ++i) v = divide(v, h);
            return v;
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
            return detail::saturating_quotient(
                T(detail::value_rounding(d.f_lo) + detail::value_rounding(d.f_hi)), d.dx);
        }

        /**
         * @brief The difference (f_hi - f_lo) / dx of two values of fn, checked: the arithmetic of difference(), and
         * the f' level of `ForwardRidders`, which has fn(x) already.
         *
         * @tparam T The real type.
         * @param f_lo fn at the lower point.
         * @param f_hi fn at the upper point.
         * @param dx The distance between the points as evaluated, > 0.
         * @return The difference with the values and the spacing, or error::non_finite if a value is not finite, or
         *         the difference or the quotient overflows.
         */
        template<real T>
        constexpr std::expected<Difference<T>, error> difference_of(T const &f_lo, T const &f_hi, T const &dx)
        {
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
            return detail::difference_of(f_lo, f_hi, dx);
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

            /** @brief The middle point minus the lower point, as evaluated, > 0. */
            T h_lo;

            /** @brief The upper point minus the middle point, as evaluated, > 0. */
            T h_hi;
        };

        /**
         * @brief The second difference ((f_hi - f_mid) / h_hi - (f_mid - f_lo) / h_lo) / ((h_hi + h_lo) / 2) of three
         * values of fn, checked: the arithmetic of second_difference(), and the f'' level of `ForwardRidders`, whose
         * middle point is x + h.
         *
         * Exact for quadratics whatever the rounding of the points, as long as the offsets are exact, which by
         * Sterbenz's lemma holds for 0 < h <= |x| (0 < 2h <= |x| for the forward points) in any radix; the rounding of
         * the points, about eps |x f'| / h^2, then drops out. In the Ridders prototype the nominal step instead gave a
         * worst case of 8.7 times the estimate for float, against 1.4. h^2 is never formed: it underflows for h below
         * about 1.5e-154 in double.
         *
         * @tparam T The real type.
         * @param f_lo fn at the lower point.
         * @param f_mid fn at the middle point.
         * @param f_hi fn at the upper point.
         * @param h_lo The middle point minus the lower point, as evaluated, > 0.
         * @param h_hi The upper point minus the middle point, as evaluated, > 0.
         * @return The second difference with the outer values and the offsets, or error::non_finite if a value is not
         *         finite, or a difference or quotient overflows.
         */
        template<real T>
        constexpr std::expected<SecondDifference<T>, error>
            second_difference_of(T const &f_lo, T const &f_mid, T const &f_hi, T const &h_lo, T const &h_hi)
        {
            // Every difference is checked before it is formed, as in difference(): two first differences that
            // overflow with the same sign would otherwise give inf - inf.
            if (!math::isfinite(f_lo) || !math::isfinite(f_mid) || !math::isfinite(f_hi))
                return std::unexpected(error::non_finite);
            if (!detail::sum_fits(f_hi, T(T(0) - f_mid)) || !detail::sum_fits(f_mid, T(T(0) - f_lo)))
                return std::unexpected(error::non_finite);
            T const up = T(f_hi - f_mid);
            T const down = T(f_mid - f_lo);
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
         * @brief The second difference of fn on x_lo < x < x_hi, from the offsets of the points actually evaluated
         * (second_difference_of()).
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
            return detail::second_difference_of(f_lo, f0, f_hi, h_lo, h_hi);
        }

        /**
         * @brief The rounding floor of a second difference: ((eps |f_hi| + eps |f_mid|) / h_hi + (eps |f_mid| +
         * eps |f_lo|) / h_lo) / ((h_hi + h_lo) / 2), its first-order sensitivity to fn's rounding (the Lagrange weights
         * of the three points times eps |f|).
         *
         * Saturated at every step: one saturation is not enough, since fn = 1e300 over h = 1e-300 overflows the inner
         * division, and fn = 1 over h = 1e-170 only the last.
         *
         * @tparam T The real type.
         * @param d The second difference.
         * @param f_mid fn at the middle point.
         * @return The floor, at most max.
         */
        template<real T>
        constexpr T second_floor(SecondDifference<T> const &d, T const &f_mid)
        {
            T const e0 = detail::value_rounding(f_mid);
            T const hi = detail::saturating_quotient(T(detail::value_rounding(d.f_hi) + e0), d.h_hi);
            T const lo = detail::saturating_quotient(T(e0 + detail::value_rounding(d.f_lo)), d.h_lo);
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
                if (c != 0 && !detail::product_fits(math::abs(fk), size)) return std::unexpected(error::non_finite);
                T const term = T(T(c) * fk);
                if (!detail::sum_fits(sum, term)) return std::unexpected(error::non_finite);
                sum = T(sum + term);
                if constexpr (Magnitude) magnitude = T(magnitude + size * detail::value_rounding(fk));
            }
            return WeightedSum<T>{ sum, magnitude, f_minus, f_plus };
        }

        /**
         * @brief The largest node of a table in magnitude, K: its points lie within K h of x.
         *
         * @tparam Size The number of nodes.
         * @param t The table.
         * @return K.
         */
        template<std::size_t Size>
        constexpr int largest_node(Table<Size> const &t)
        {
            int const lo = t.nodes[0];
            int const hi = t.nodes[Size - 1];
            return -lo > hi ? -lo : hi;
        }

        /**
         * @brief max(2K, denominator) of a table, the q of fits(): 2, 2, 4, 4 for `Central<2>` and f' to f'''', 2N
         * for `Forward<1>`.
         *
         * @tparam Size The number of nodes.
         * @param t The table.
         * @return q.
         */
        template<std::size_t Size>
        constexpr int reach(Table<Size> const &t)
        {
            int const k = 2 * detail::largest_node(t);
            return k > t.denominator ? k : t.denominator;
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
            if (!detail::fits(x, h, detail::reach(table))) return std::unexpected(error::non_finite);

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
                    // N + order - 1. The divisions are ordered so that they neither underflow nor overflow
                    // spuriously (divide_by_powers). Node rounding adds about
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
                    T const d = detail::divide_by_powers<N, false>(s->value, T(T(table.denominator) * h), h);
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
         * @brief The levels an extrapolation method may attempt in all, restarts included: at most 64 evaluations of
         * fn for f', 65 for f'', 128 for f''' and 129 for f'''' by `Ridders`, and 33, 65, 97 and 129 by the one-sided
         * methods (fn(x) once, then N per level).
         *
         * Only a guard: a result takes 5 to 11 levels for f' in double and float, 10 to 19 for 50 and 100 digits, up
         * to 18 near an edge; f''' and f'''' take up to 20 levels for 100 digits, and 28 after restarts at a pole.
         * Near a domain edge the one-sided methods may use all of them.
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
         * @brief The divisor of the one-sided methods' default initial step max(|x|, 1) / divisor, for every order.
         *
         * Their rules reach N h. In the prototype (docs/deriv/PHASE5_API.md) s/32 was, for f' to f''', the more
         * accurate of the two starts without a gross error (s/64 lost 1.4 to 2 times; s/16 had one for f'' and f''',
         * s/4 five), and for f'''' it had the fewest gross errors (5 of 2,700, against 9 at s/16 and 34 at s/64) and a
         * worst case of 1.8 times the estimate, against 12 at s/16.
         */
        inline constexpr int one_sided_divisor = 32;

        /**
         * @brief One level of an extrapolation tableau: the base rule's estimate and its rounding floor.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct Level
        {
            /** @brief The base rule's estimate of the derivative at this step. */
            T value;

            /** @brief Its rounding error when fn is accurate to about eps, at most max. */
            T floor;
        };

        /**
         * @brief One level of an extrapolation tableau for the N-th derivative at step h: the rule of order N of the
         * base S, `Central<2>` (`Ridders`) or `Forward<1>` (the one-sided methods).
         *
         * f' and f'' use the actual distance and offsets (difference_of(), second_difference_of()): with the nominal
         * step, the prototype's worst case was 22 against 3.6 times the estimate for central f', and 8.7 against 1.4
         * for float f''. The forward f' is (f(x+h) - f(x)) / h, its f'' the second difference on x, x + h, x + 2h.
         *
         * f''' and f'''' use the nominal step with S's table: centrally (f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)) / 2h^3
         * and (f(x+2h) - 4f(x+h) + 6f(x) - 4f(x-h) + f(x-2h)) / h^4, with errors h^2 f^(5)/4 and h^2 f^(6)/6 at the
         * lowest order, whose expansions are even as the central tableau needs; forward the N-th forward difference
         * over h^N, with an error N/2 h f^(N+1) + O(h^2). Their floor has two parts, each saturated at every step:
         * fn's rounding, sum |w_k| eps |f_k| / (denominator h^N), and the rounding of the points, which the nominal
         * step leaves in: each point x + k h is off by up to eps/2 (|x| + K h), K the largest node, which moves its
         * value by about that times f', so eps/2 (|x| + K h) |f'| sum |w_k| / (denominator h^N), with f' from the
         * values at x +- h (central) or at x and x + h (forward). Without the second part, the prototype's central f'''
         * of exp(10x) on [1, 3] was covered in 96.3% of the cases, with it in 99.0%.
         *
         * @tparam S The base stencil.
         * @tparam N The derivative order.
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point, finite.
         * @param f0 fn(x) where the table has node 0 (even N centrally, always forward); unused otherwise.
         * @param h The step, with every point finite and distinct.
         * @return The level, or error::non_finite.
         */
        template<typename S, int N, real T, typename F>
        constexpr std::expected<Level<T>, error> ridders_level(F &fn, T const &x, T const &f0, T const &h)
        {
            constexpr bool central = std::same_as<S, Central<2>>;
            if constexpr (N == 1 && central) {
                auto const d = detail::difference(fn, T(x - h), T(x + h));
                if (!d) return std::unexpected(d.error());
                return Level<T>{ d->value, detail::rounding_floor(*d) };
            } else if constexpr (N == 1) {
                T const x_1 = T(x + h);
                T const dx = T(x_1 - x);
                T const f_1 = static_cast<T>(std::invoke(fn, x_1));
                auto const d = detail::difference_of(f0, f_1, dx);
                if (!d) return std::unexpected(d.error());
                return Level<T>{ d->value, detail::rounding_floor(*d) };
            } else if constexpr (N == 2 && central) {
                auto const d = detail::second_difference(fn, x, f0, T(x - h), T(x + h));
                if (!d) return std::unexpected(d.error());
                return Level<T>{ d->value, detail::second_floor(*d, f0) };
            } else if constexpr (N == 2) {
                T const x_1 = T(x + h);
                T const x_2 = T(x + T(2) * h);
                T const f_1 = static_cast<T>(std::invoke(fn, x_1));
                T const f_2 = static_cast<T>(std::invoke(fn, x_2));
                auto const d = detail::second_difference_of(f0, f_1, f_2, T(x_1 - x), T(x_2 - x_1));
                if (!d) return std::unexpected(d.error());
                return Level<T>{ d->value, detail::second_floor(*d, f_1) };
            } else {
                constexpr auto table = detail::make_table<N>(S{});
                constexpr int reach = detail::largest_node(table);
                std::optional<T> const cached = f0;
                auto const s = detail::weighted_sum<true>(fn, table, x, h, cached);
                if (!s) return std::unexpected(s.error());
                constexpr int weights = [](auto const &t) {
                    int w = 0;
                    for (int const c : t.numerators) w += c < 0 ? -c : c;
                    return w;
                }(table);
                T const eps = std::numeric_limits<T>::epsilon();
                T const slope = central
                                    ? detail::saturating_quotient(detail::distance(s->f_plus, s->f_minus), T(T(2) * h))
                                    : detail::saturating_quotient(detail::distance(s->f_plus, f0), h);
                T const spread = T(T(eps / T(2)) * T(math::abs(x) + T(reach) * h)) * T(weights);
                T const divisor = T(table.denominator) * h;
                T const value = detail::divide_by_powers<N, false>(s->value, divisor, h);
                T const rounding = detail::divide_by_powers<N, true>(s->magnitude, divisor, h);
                T const points =
                    detail::divide_by_powers<N, true>(detail::saturating_product(spread, slope), divisor, h);
                if (!math::isfinite(value)) return std::unexpected(error::non_finite); // the quotient overflowed
                return Level<T>{ value, detail::saturating_sum(rounding, points) };
            }
        }

        /**
         * @brief The rounding floor of an extrapolated entry, from those of its parents: (r f_new + f_old) / (r - 1)
         * for the column factor r, saturated at every step, and max (unknown) stays max.
         *
         * A[i][j] = (r A[i][j-1] - A[i-1][j-1]) / (r - 1), and the coefficients of the levels in an entry alternate in
         * sign, so this is exactly sum |c_k| floor_k, the first-order bound on the entry's rounding.
         *
         * @tparam T The real type.
         * @param f_new The floor of A[i][j-1].
         * @param f_old The floor of A[i-1][j-1].
         * @param factor r, 2^j for the one-sided methods.
         * @return The floor of A[i][j], at most max.
         */
        template<real T>
        constexpr T propagated_floor(T const &f_new, T const &f_old, T const &factor)
        {
            T const largest = (std::numeric_limits<T>::max)();
            T const sum = detail::saturating_sum(detail::saturating_product(f_new, factor), f_old);
            if (sum == largest) return largest;
            return T(sum / T(factor - T(1)));
        }

        /**
         * @brief Ridders' extrapolation of the rule of order N of the base S at x, from the initial step h:
         * `Central<2>` for `Ridders`, `Forward<1>` for the one-sided methods.
         *
         * From Ridders (1982) and Richardson extrapolation with its a posteriori error estimate (Burden et al.,
         * Numerical Analysis; Sidi, Practical Extrapolation Methods, 2003). The constants are our own, chosen with a
         * prototype (docs/deriv/PHASE3_API.md, docs/deriv/PHASE4_API.md, docs/deriv/PHASE5_API.md): ratio 2, two
         * levels of stall, the rounding floor, the untrusted-result test; the same for every order and base.
         *
         * Level i is the base rule D_i at h_i = h_start / 2^i. The central rule's error expansion has even powers of h
         * only, so column j of the tableau removes the h^(2j) term, with r = 4^j; the forward rule's has every power,
         * so column j removes h^j, with r = 2^j. In Burden's difference form
         *
         *     A[i][j] = A[i][j-1] + (A[i][j-1] - A[i-1][j-1]) / (r - 1),
         *
         * which is exact when the two parents are equal. Only the previous row is kept.
         *
         * What differs between the bases follows from S's table (the reach q, whether fn(x) is needed, the points
         * that must differ), except two choices: the factor r, and the floor of an entry. Centrally it is the floor
         * of the newest level, as tuned in Phase 3; one-sided it is carried through the tableau (propagated_floor()),
         * since columns divided by 2^j - 1 amplify the newest level's rounding about 3.5 times rather than 1.4: with
         * the level's floor the prototype's forward f' was covered in 97.5% of the cases (worst 8.9 times the
         * estimate), with the propagated one in 99.4% (worst 2.6).
         *
         * @tparam S The base stencil, `Central<2>` or `Forward<1>`.
         * @tparam N The derivative order.
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point, finite.
         * @param h The initial step, > 0 and finite.
         * @return The best entry with its error estimate; if the levels run out during restarts, the last untrusted
         *         result; error::non_finite if fn(x) is not finite (where the table needs it) or no tableau reached
         *         two finite levels.
         */
        template<typename S, int N, real T, typename F>
        constexpr std::expected<Estimate<T>, error> ridders(F &fn, T const &x, T const &h)
        {
            constexpr auto table = detail::make_table<N>(S{});
            constexpr std::size_t size = table.nodes.size();
            // The points of a level reach x +- K h (K = 1 or 2 centrally, N forward).
            constexpr int q = detail::reach(table);
            constexpr bool one_sided = table.nodes[0] == 0;
            constexpr int ratio = one_sided ? 2 : 4;
            // The nodes and 0 are the consecutive integers lo ... hi for both bases: these points must differ.
            constexpr int lo = table.nodes[0] < 0 ? table.nodes[0] : 0;
            constexpr int hi = table.nodes[size - 1] > 0 ? table.nodes[size - 1] : 0;
            constexpr bool needs_f0 = [](auto const &t) {
                for (int const k : t.nodes)
                    if (k == 0) return true;
                return false;
            }(table);
            // fn(x) is part of every level where the table has node 0 (even N centrally, always forward): evaluated
            // once, and kept across restarts. If it is not finite, no step avoids it.
            T const f0 = needs_f0 ? static_cast<T>(std::invoke(fn, x)) : T(0);
            if (!math::isfinite(f0)) return std::unexpected(error::non_finite);
            auto const point = [&x](T const &step, int k) -> T {
                if (k == 0) return x;
                return T(x + T(k) * step);
            };

            std::array<T, ridders_levels> prev = detail::filled<ridders_levels>(T(0));
            std::array<T, ridders_levels> row = prev;
            // The floors of the entries, carried through the tableau; used by the one-sided methods only.
            std::array<T, ridders_levels> floor_prev = prev;
            std::array<T, ridders_levels> floor_row = prev;
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
                    // Checked before evaluating, and after the overflow check, so that the points are only formed when
                    // they are finite: below the resolution of x the points coincide, and a tiny h could underflow to
                    // 0. Neighbouring points must differ too (at h near 0.6 ulp, x + h and x + 2h round to the same
                    // value).
                    bool const fit = detail::fits(x, h_i, q);
                    bool coincide = false;
                    if (fit)
                        for (int k = lo; k < hi && !coincide; ++k) coincide = point(h_i, k) == point(h_i, k + 1);
                    if (!(h_i > T(0)) || coincide) {
                        resolution = true;
                        break;
                    }
                    --budget;
                    std::expected<Level<T>, error> const level =
                        fit ? detail::ridders_level<S, N>(fn, x, f0, h_i)
                            : std::expected<Level<T>, error>(std::unexpect, error::non_finite);
                    bool finite = level.has_value();
                    if (finite) {
                        row[0] = level->value;
                        floor_row[0] = level->floor;
                        T factor = T(1);
                        // Both the difference and the sum are checked before they are formed: finite entries of
                        // opposite sign can differ by more than max.
                        for (std::size_t j = 1; j <= i && finite; ++j) {
                            factor = T(factor * T(ratio));
                            finite = detail::sum_fits(row[j - 1], T(T(0) - prev[j - 1]));
                            if (!finite) break;
                            T const correction = T((row[j - 1] - prev[j - 1]) / (factor - T(1)));
                            finite = detail::sum_fits(row[j - 1], correction);
                            if (!finite) break;
                            row[j] = T(row[j - 1] + correction);
                            floor_row[j] = detail::propagated_floor(floor_row[j - 1], floor_prev[j - 1], factor);
                        }
                    }
                    // fn is not finite within reach of x (or an entry would exceed max): every earlier level's
                    // interval contained the bad point, so the whole tableau goes, and a new one starts below the
                    // failing step.
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
                    bool improved = false;
                    bool exact = false;
                    for (std::size_t j = 1; j <= i; ++j) {
                        T const entry_floor = one_sided ? floor_row[j] : level->floor;
                        T const change = detail::distance(row[j], prev[j - 1]);
                        if (change == T(0)) {
                            best = row[j];
                            best_error = entry_floor;
                            best_floor = entry_floor;
                            have_best = exact = true;
                            break;
                        }
                        T const estimate = change < entry_floor ? entry_floor : change;
                        if (!have_best || estimate < best_error) {
                            best = row[j];
                            best_error = estimate;
                            best_floor = entry_floor;
                            have_best = improved = true;
                        }
                    }
                    diagonal = row[i];
                    std::swap(prev, row);
                    std::swap(floor_prev, floor_row);
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
                // producing a non-finite value (a pole near x, aliasing, or, for a one-sided base, fn varying on a
                // scale far below h_start), so start again. A derivative that is 0 has an estimate near the floor and
                // is not affected. In the prototype, the factor 1 instead of 4 left 11 of 400 poles 1/(x - a)
                // uncovered for f' (32 for f''), 4 and 16 none; 16 to 256 for the floor gave identical results.
                // Written with divisions, which cannot overflow: the floor, and so the estimate, saturates at max.
                bool const untrusted = abs_error >= T(math::abs(best) / T(4)) && T(abs_error / T(64)) > best_floor;
                if (!untrusted || resolution) return result;
                fallback = result;
                // The new tableau starts below the smallest step this one tried, h_i (the loop ended on a stall or
                // exact agreement, before h_i was halved), as after a non-finite level: restarting at h_start / 2
                // re-evaluated the same large steps, and the levels ran out before the step reached fn's scale. With
                // it, central Ridders missed poles 1e-3 and 1e-4 from x (f'' to f'''' off by 1 to 11, uncovered) and
                // the one-sided methods failed at sqrt and log at 1e-6 and 1e-8. A case the slower descent handled can
                // lose up to 8 times in accuracy; restarts overlapping one or two levels did worse on every measure.
                h_start = T(h_i / T(2));
            }
            if (fallback) return *fallback;
            return std::unexpected(error::non_finite);
        }

        /**
         * @brief fn(-t), for `BackwardRidders`: fn is called with the T const lvalue -t, as everywhere else.
         *
         * A struct with a pointer rather than a lambda that captures fn by reference: MSVC 19.51 rejects such a capture
         * of a function reference in constant evaluation ("a read of an uninitialized symbol"), where GCC and Clang
         * accept it.
         *
         * @tparam F The callable type, possibly a function type.
         * @tparam T The real type.
         */
        template<typename F, real T>
        struct Mirrored
        {
            /** @brief The function. */
            F *fn;

            /**
             * @brief fn at -t.
             * @param t The point.
             * @return fn(-t), converted to T.
             */
            constexpr T operator()(T const &t) const
            {
                T const point = -t;
                return static_cast<T>(std::invoke(*fn, point));
            }
        };

        /**
         * @brief The N-th derivative of fn at x by method M with step h, behind the public overloads.
         *
         * @tparam M The method.
         * @tparam N The derivative order, 1 to 4.
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point.
         * @param h The step.
         * @return The estimate, or error::non_finite.
         */
        template<typename M, int N, real T, typename F>
        constexpr std::expected<Estimate<T>, error> derivative_with_step(F &fn, T const &x, T const &h)
        {
            // The step is the caller's input, often computed: validated here, before fn is called, the same way in
            // debug and release builds and in constant expressions (comparisons only, so a NaN is no arithmetic).
            if (!(h > T(0)) || !math::isfinite(h)) return std::unexpected(error::non_finite);
            if constexpr (std::same_as<M, BackwardRidders>) {
                if (!math::isfinite(x)) return std::unexpected(error::non_finite);
                // ForwardRidders on fn(-t) at -x. Negation is exact and rounding symmetric, so fn receives exactly the
                // points x - k h, and the result is that of a backward base by construction; a Backward<1> table would
                // sum its terms from x - N h and round differently. The N-th derivative changes sign for odd N.
                Mirrored<F, T> mirrored{ std::addressof(fn) };
                return detail::ridders<Forward<1>, N>(mirrored, T(-x), h).transform([](Estimate<T> e) {
                    if constexpr (N % 2 == 1) e.value = T(-e.value);
                    return e;
                });
            } else if constexpr (detail::extrapolation<M>) {
                if (!math::isfinite(x)) return std::unexpected(error::non_finite);
                return detail::ridders<detail::base<M>, N>(fn, x, h);
            } else {
                return detail::stencil_derivative<M, N>(fn, x, h).transform(
                    [](T const &d) { return Estimate<T>{ d, std::nullopt }; });
            }
        }

    } // namespace detail

    /**
     * @brief The step that `derivative<M, N>(fn, x)` uses when no step is passed.
     *
     * For a stencil of accuracy order P it is eps^(1/(P+N)) * max(|x|, 1); for Ridders it is the initial step,
     * max(|x|, 1) / 32 for f' and max(|x|, 1) / 4 for f'' to f'''', and for `ForwardRidders` and `BackwardRidders`
     * max(|x|, 1) / 32 for every order. For every x, `derivative<M, N>(fn, x)` and
     * `derivative<M, N>(fn, x, default_step<M, N>(x))` give the same result (for a non-finite x a non-finite step,
     * and error::non_finite from both), so a caller can repeat a rule at twice or
     * half the default step, or apply it by hand, without rebuilding the formula.
     *
     * @code
     * double const h = default_step<Central<2>>(x);            // what derivative(fn, x) uses
     * auto const coarse = derivative<Central<2>>(fn, x, 2 * h); // the same rule at twice the step
     * @endcode
     *
     * @tparam M The method: `Central<P>` (by default `Central<2>`), `Forward<P>`, `Backward<P>`, `Ridders`,
     *           `ForwardRidders` or `BackwardRidders`.
     * @tparam N The derivative order: 1 (the default) to 4, an integer or an `order`.
     * @tparam T The real type, deduced from @p x.
     * @param x The point.
     * @return The step: positive and finite for a finite x, otherwise not finite.
     */
    template<typename M = Central<2>, auto N = 1, real T>
    requires detail::method<M> && detail::valid_order<N>
    [[nodiscard]] constexpr T default_step(T const &x)
    {
        constexpr int n = static_cast<int>(N);
        if constexpr (std::same_as<M, Ridders>) {
            constexpr int divisor = detail::ridders_divisors[static_cast<std::size_t>(n - 1)];
            return T(detail::scale(x) / T(divisor));
        } else if constexpr (detail::extrapolation<M>) {
            return T(detail::scale(x) / T(detail::one_sided_divisor));
        } else {
            constexpr int exponent = detail::make_table<n>(M{}).order + n;
            return detail::stencil_step<exponent>(x);
        }
    }

    namespace detail {

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
            // Checked first, so that a non-finite x is an error value in constant expressions too: its step would be
            // NaN or infinite, and Clang rejects arithmetic that produces a NaN at compile time.
            if (!math::isfinite(x)) return std::unexpected(error::non_finite);
            // The public default_step, so that derivative(fn, x) and derivative(fn, x, default_step(x)) agree.
            return detail::derivative_with_step<M, N>(fn, x, nxx::deriv::default_step<M, N>(x));
        }

    } // namespace detail

    /**
     * @brief The @p N-th derivative of @p fn at @p x by method @p M with step @p h.
     *
     * fn is called once per node of a stencil, and only at finite points; Ridders calls it 2 or 4 times per level, at
     * most 64 to 129 times, and the one-sided extrapolations fn(x) once and N times per level, at most 33 to 129
     * times, never on the other side of x. Failures are detected before the arithmetic that would produce them, so they
     * are error values in constant expressions too, except an overflow of a final quotient, which GCC rejects at
     * compile time. Exceptions thrown by fn or by T pass through unchanged, and the result's value() throws on a
     * failure.
     *
     * fn is assumed smooth on the points the method samples: a kink or a jump among them gives a wrong value without
     * an error. A kink at x itself is invisible even to Ridders' estimate: a central rule returns the mean of the two
     * one-sided slopes there, and Ridders converges to it; `ForwardRidders` and `BackwardRidders` return the one-sided
     * derivatives.
     *
     * @code
     * auto const d3 = derivative<Central<2>, order::third>(fn, x, h); // f''' on x +- h, x +- 2h: Chapra's O(h^2)
     * auto const d4 = derivative<Backward<2>, 4>(fn, x, h);           // f'''' on x - 5h ... x
     * @endcode
     *
     * @tparam M The method: `Central<P>` (by default `Central<2>`), `Forward<P>`, `Backward<P>`, `Ridders`,
     *           `ForwardRidders` or `BackwardRidders`.
     * @tparam N The derivative order: 1 (the default) to 4, an integer or an `order`.
     * @tparam T The real type, deduced from @p x.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @param fn The function to differentiate.
     * @param x The point.
     * @param h The step of a stencil, or the initial step of an extrapolation method; it converts to T, so a
     *          double step works for a multiprecision T. A step that is not positive and finite gives
     *          error::non_finite before fn is called; a stencil needs h well above ulp(x).
     * @return The estimate of the derivative, with an error estimate for the extrapolation methods, or
     *         error::non_finite.
     */
    // h is a non-deduced context (type_identity_t), so T comes from x alone and a double step converts to a
    // multiprecision T instead of making deduction fail. N is auto, because the scoped order does not convert to int.
    template<typename M = Central<2>, auto N = 1, real T, typename F>
    requires detail::method<M> && detail::valid_order<N> && nxx::detail::callable<F, T>
             && detail::exact_weights<M, N, T>
    [[nodiscard]] constexpr std::expected<Estimate<T>, error>
        derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
    {
        // Qualified, so that argument-dependent lookup cannot pick a function from the namespace of T or F.
        return detail::derivative_with_step<M, static_cast<int>(N)>(fn, x, h);
    }

    /**
     * @brief The @p N-th derivative of @p fn at @p x by method @p M with its default step.
     *
     * For a stencil of order P the step is eps^(1/(P+N)) * max(|x|, 1). It balances truncation against rounding when fn
     * varies on the scale max(|x|, 1), and the error is then about eps^(P/(P+N)); pass a step when fn varies faster, or
     * when the derivative is much smaller than that error. The sensitivity grows with N: for f'''', misjudging fn's
     * scale by a factor 2 costs up to 16 in error. For Ridders the initial step is max(|x|, 1) / 32 for f' and
     * max(|x|, 1) / 4 for f'' to f'''', for its one-sided forms max(|x|, 1) / 32, and the method adapts it;
     * default_step<M, N>(x) returns the step. The default steps differ between methods and orders, so generic code that
     * varies N should use this overload rather than pass one h for every order.
     *
     * @tparam M The method: `Central<P>` (by default `Central<2>`), `Forward<P>`, `Backward<P>`, `Ridders`,
     *           `ForwardRidders` or `BackwardRidders`.
     * @tparam N The derivative order: 1 (the default) to 4, an integer or an `order`.
     * @tparam T The real type, deduced from @p x.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @param fn The function to differentiate.
     * @param x The point.
     * @return The estimate of the derivative, with an error estimate for the extrapolation methods, or
     *         error::non_finite.
     */
    template<typename M = Central<2>, auto N = 1, real T, typename F>
    requires detail::method<M> && detail::valid_order<N> && nxx::detail::callable<F, T>
             && detail::exact_weights<M, N, T>
    [[nodiscard]] constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x)
    { return detail::derivative_default<M, static_cast<int>(N)>(fn, x); }

    /**
     * @brief The overload with a step, for a derivative above the first in a type with an epsilon above 1e-7 (float):
     * the same result, with a deprecation warning, since that result has few correct digits (detail::noisy).
     *
     * Chosen over the general overload by its additional constraint. The warning is reported at the caller's line,
     * so a pragma there silences it (`-Wdeprecated-declarations`, C4996 on MSVC); NUMERIXX_NO_PRECISION_WARNING turns
     * it off everywhere.
     *
     * @tparam M The method.
     * @tparam N The derivative order, 2 to 4.
     * @tparam T The real type.
     * @tparam F The callable.
     * @param fn The function to differentiate.
     * @param x The point.
     * @param h The step, or the initial step of an extrapolation method.
     * @return As for the general overload.
     */
    template<typename M = Central<2>, auto N = 1, real T, typename F>
    requires detail::method<M> && detail::valid_order<N> && nxx::detail::callable<F, T>
             && detail::exact_weights<M, N, T> && detail::noisy<N, T>
    [[deprecated(
        "derivatives above the first in a type with 7 or fewer significant digits keep few correct digits: "
        "use double, or silence this warning at the call or with NUMERIXX_NO_PRECISION_WARNING")]] [[nodiscard]]
    constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
    { return detail::derivative_with_step<M, static_cast<int>(N)>(fn, x, h); }

    /**
     * @brief The overload with the default step, for a derivative above the first in a type with an epsilon above
     * 1e-7 (float): the same result, with a deprecation warning, since that result has few correct digits.
     *
     * @tparam M The method.
     * @tparam N The derivative order, 2 to 4.
     * @tparam T The real type.
     * @tparam F The callable.
     * @param fn The function to differentiate.
     * @param x The point.
     * @return As for the general overload.
     */
    template<typename M = Central<2>, auto N = 1, real T, typename F>
    requires detail::method<M> && detail::valid_order<N> && nxx::detail::callable<F, T>
             && detail::exact_weights<M, N, T> && detail::noisy<N, T>
    [[deprecated(
        "derivatives above the first in a type with 7 or fewer significant digits keep few correct digits: "
        "use double, or silence this warning at the call or with NUMERIXX_NO_PRECISION_WARNING")]] [[nodiscard]]
    constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x)
    { return detail::derivative_default<M, static_cast<int>(N)>(fn, x); }

    namespace detail {

        /**
         * @brief The callable that derivative_of() returns: a copy of fn and the method, called as derivative<M, N>.
         *
         * The call operator has an explicit object parameter, so fn is called with the callable's own constness, as
         * the result of std::bind_front calls its target: a `mutable` fn can be called through a non-const callable,
         * not through a const one. Its constraint states derivative()'s conditions directly rather than naming
         * derivative(), whose overload resolution would select the deprecated twin for a coarse type above f'.
         *
         * @tparam M The method.
         * @tparam N The derivative order, 1 to 4.
         * @tparam G The stored callable: fn decayed.
         */
        template<typename M, int N, typename G>
        struct DerivativeOf
        {
            /** @brief The function, a copy (or a std::reference_wrapper of it). */
            G fn;

            /**
             * @brief Direct-initialises fn, so that an explicit copy or move constructor of G is used, as
             * derivative_of's constraint (std::constructible_from) promises; an aggregate's member would be
             * copy-initialised (independent review of Phase 3).
             * @tparam H The argument's type.
             * @param f The function.
             */
            template<typename H>
            constexpr DerivativeOf(std::in_place_t, H &&f) : fn(std::forward<H>(f))
            {}

            /**
             * @brief derivative<M, N>(fn, x), with the default step.
             * @tparam Self The callable's own type, with its constness and value category.
             * @tparam T The real type, deduced from @p x.
             * @param self The callable.
             * @param x The point.
             * @return The estimate of the derivative, or error::non_finite.
             */
            template<typename Self, real T>
            requires nxx::detail::callable<decltype((std::declval<std::remove_reference_t<Self> &>().fn)), T>
                     && detail::exact_weights<M, N, T>
            constexpr std::expected<Estimate<T>, error> operator()(this Self &&self, T const &x)
            { return detail::derivative_default<M, N>(self.fn, x); }
        };

    } // namespace detail

    /**
     * @brief A derivative of @p fn as a callable: d(x) is `derivative<M, N>(fn, x)`, with the default step.
     *
     * For an algorithm that takes a derivative as a callable, such as Newton's method in nxx::roots:
     * `solve<Newton>(fn, derivative_of(fn), x0)`. d(x) is valid for every real T for which `derivative<M, N>(fn, x)` is,
     * and gives the same result, bit for bit, a std::expected<Estimate<T>, error> whose error estimate (Ridders') is
     * passed through.
     *
     * The callable holds a copy of fn, decayed (a free function becomes a pointer, a lambda is copied or moved), so it
     * stays valid when fn's scope ends. A large fn, such as one holding a table, is copied too: pass `std::ref(fn)` to
     * share it, and keep fn alive while the callable is used. A stateful fn (a `mutable` lambda, a cache) is copied, so
     * its state diverges from the original's; `std::ref(fn)` shares it. A copy of a `mutable` fn is called through a
     * non-const callable (`auto d = derivative_of(g)`), not through a const one. `derivative_of(std::move(fn))` moves
     * from fn.
     *
     * The step is derivative()'s default, which assumes that fn varies on the scale max(|x|, 1); for another step,
     * write the lambda `[&](double x) { return derivative<M, N>(fn, x, h); }`. Unlike a direct call of derivative(), a
     * call of the callable does not warn for a derivative above the first in a coarse type such as float, where it
     * keeps few correct digits: the type is known only where an algorithm calls it, inside a header, where the warning
     * could not be silenced at the caller's line.
     *
     * @code
     * auto const d = derivative_of(fn);                          // f' by Central<2>
     * auto const d2 = derivative_of<Ridders, order::second>(fn); // f'' by Ridders
     * if (auto const r = d(1.0)) use(r->value);
     * @endcode
     *
     * @tparam M The method: `Central<P>` (by default `Central<2>`), `Forward<P>`, `Backward<P>`, `Ridders`,
     *           `ForwardRidders` or `BackwardRidders`.
     * @tparam N The derivative order: 1 (the default) to 4, an integer or an `order`.
     * @tparam F The callable type; its decayed type must be constructible from it.
     * @param fn The function to differentiate.
     * @return The callable, of an unspecified type.
     */
    template<typename M = Central<2>, auto N = 1, typename F>
    requires detail::method<M> && detail::valid_order<N> && std::constructible_from<std::decay_t<F>, F>
    [[nodiscard]] constexpr auto derivative_of(F &&fn)
    { return detail::DerivativeOf<M, static_cast<int>(N), std::decay_t<F>>{ std::in_place, std::forward<F>(fn) }; }

} // namespace nxx::deriv

/** @brief The standard library's formatter, specialised for nxx::deriv::error. */
namespace std {

    /**
     * @brief Formats an nxx::deriv::error as the name of its enumerator, such as `non_finite`, with the format
     * specification of a string.
     *
     * @code
     * if (!d) std::println("no estimate: {}", d.error());   // no estimate: non_finite
     * @endcode
     */
    template<>
    struct formatter<nxx::deriv::error> : formatter<string_view>
    {
        /**
         * @brief Writes the name of @p e.
         *
         * @tparam FormatContext The format context.
         * @param e The error.
         * @param ctx The format context.
         * @return The end of the output.
         */
        template<typename FormatContext>
        auto format(nxx::deriv::error e, FormatContext &ctx) const
        {
            // No default: a new enumerator makes the compilers warn here (-Wswitch) until its name is added.
            string_view name = "unknown";
            switch (e) {
            case nxx::deriv::error::non_finite:
                name = "non_finite";
                break;
            }
            return formatter<string_view>::format(name, ctx);
        }
    };

} // namespace std
