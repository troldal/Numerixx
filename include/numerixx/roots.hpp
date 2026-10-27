/**
 * @file
 * @brief Roots of a real function of one variable: bracketing methods and open methods.
 *
 * solve<M>(fn, start) finds a point where fn changes sign, by a method M. From a bracket, an interval at whose ends fn
 * has opposite signs: `Bisection`, which halves it, `Brent`, which interpolates and falls back on bisection, and
 * `Newton` and `Halley`, whose derivative steps the bracket safeguards. From a guess: `Newton`, `Halley` and `Secant`,
 * which iterate until a sign change appears between two iterates and then finish inside it. Every method returns a
 * std::expected<Root<T>, error>: the root and a bound on its distance to a sign change of fn. Numerixx itself never
 * throws, but the result's own value() throws std::bad_expected_access on a failure, so test the result and read it
 * with `->` or `*`.
 *
 * Newton and Halley take fn's derivatives either from fn, which then returns a std::pair, std::tuple or std::array of
 * exactly T (the value, f', f''), or from separate callables passed after fn.
 *
 * @code
 * namespace rt = nxx::roots;
 * auto const r = rt::solve<rt::Brent>([](double x) { return x * x - 2.0; }, rt::Bracket{1.0, 2.0});
 * if (r) use(r->value, r->abs_error); // sqrt(2) to about 1e-16, after 9 evaluations of fn
 * auto const fdf = [](double x) { return std::pair{x * x - 2.0, 2.0 * x}; };
 * auto const n = rt::solve<rt::Newton>(fdf, 1.0); // the same root from a guess, after 7 evaluations
 * @endcode
 */

#pragma once

#include <numerixx/math.hpp>

#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <expected>
#include <format>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

/** @brief Roots of real functions of one variable. */
namespace nxx::roots {

    /**
     * @brief Why solve() returned no root.
     *
     * More values will be added, and a listed cause may move to a new value: a switch on it should have a default, and
     * a caller should test for success rather than for a particular error. It is formattable: std::format("{}", e)
     * writes the enumerator's name, such as `no_sign_change`.
     */
    enum class error {
        /**
         * A bracket end, the guess, a value of fn or a tolerance is not finite, or a tolerance is negative. A value of
         * fn that is infinite at a bracket end is reported here; a later version may accept it for its sign. The
         * derivatives are not checked: one that is 0 or not finite makes a derivative step unusable (see Newton).
         */
        non_finite,

        /** fn has the same sign at both ends of the bracket, and is 0 at neither. */
        no_sign_change,

        /**
         * An open method found no root from its guess: no sign change between its iterates within its iteration cap,
         * or no usable step (f' = 0 or not finite, Secant's two values equal, a step that would overflow).
         */
        not_converged
    };

    /**
     * @brief An interval [lower, upper] in which fn changes sign. The ends may be given in either order.
     *
     * @code
     * auto const b = nxx::roots::Bracket{300.0, 400.0}; // deduces Bracket<double>
     * @endcode
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Bracket
    {
        /** @brief One end of the interval, usually the lower. */
        T lower;

        /** @brief The other end, usually the upper. */
        T upper;

        /** @brief Compares both ends. */
        friend constexpr bool operator==(Bracket const &, Bracket const &) = default;
    };

    /**
     * @brief A root, and a bound on its distance to a sign change of fn.
     *
     * abs_error bounds the distance from value to a point where fn, as evaluated, changes sign, and is at most
     * 2 tol(value) (see Criteria), for every method. It saturates at max where that distance exceeds what T holds
     * (only with tolerances near max). The root of the exact function may differ by fn's own rounding error divided by
     * |f'|, which no method can see. Where fn evaluates to exactly 0, abs_error is the requested tolerance tol(value),
     * not 0: such a point is only as good as fn's evaluation (`(x + 1e10) - (1e10 + 1.0 / 3)` is 0 on an interval
     * about 2e-6 wide). Read the members by name: members may be added in a later version, which would break a
     * structured binding.
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Root
    {
        /** @brief The root. */
        T value;

        /** @brief A bound on the distance from value to a sign change of fn, or the tolerance at an exact zero. */
        T abs_error;

        /** @brief Compares both members. */
        friend constexpr bool operator==(Root const &, Root const &) = default;
    };

    /**
     * @brief When to stop: when the bracket is at most 2 tol(x) wide, with tol(x) = abs_tol + rel_tol |x| at the
     * current estimate x, or when fn is exactly 0.
     *
     * The defaults ask for the full precision of T at every magnitude and assume no scale, so a root at 1e-12 comes
     * out to about 2 eps relative. They are also the minimum: smaller tolerances are raised to them, and larger ones
     * are capped at max / 4 and 1. Their cost falls on a sign change near 0 that fn does not evaluate to exactly 0: a
     * root that fn misses at 0, a tiny root or a jump at 0 is approached down to T's smallest numbers, about
     * log2(width / denorm_min) halvings (about 1100 in double, but about 2e9 in a Boost.Multiprecision type, whose
     * exponent range is much wider). Set abs_tol to the distance from 0 that is close enough there.
     *
     * @code
     * auto const r = nxx::roots::solve<nxx::roots::Brent>(fn, nxx::roots::Bracket{0.0, 1.0}, {.abs_tol = 1e-12});
     * @endcode
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Criteria
    {
        /** @brief The absolute tolerance, for roots near 0; at least denorm_min. */
        T abs_tol = T(std::numeric_limits<T>::denorm_min());

        /** @brief The relative tolerance; at least 2 epsilon. */
        T rel_tol = T(T(2) * T(std::numeric_limits<T>::epsilon()));

        /** @brief Compares both members. */
        friend constexpr bool operator==(Criteria const &, Criteria const &) = default;
    };

    /**
     * @brief Bisection: halves the bracket at every step, one bit per evaluation of fn, guaranteed.
     *
     * About log2(width / (2 tol)) evaluations, 52 for [1, 2] in double, whatever fn is.
     */
    struct Bisection
    {
    };

    /**
     * @brief Brent's method: inverse quadratic and linear interpolation, safeguarded by bisection (Brent 1973).
     *
     * Superlinear on a smooth fn with a simple root, typically 8 to 12 evaluations in double, and never much slower
     * than Bisection: at most (k + 1)^2 - 2 evaluations, k being Bisection's count. Slower on a multiple root
     * ((x - 1)^3: 156 evaluations, against Bisection's 54).
     */
    struct Brent
    {
    };

    /**
     * @brief Newton's method: the step -f / f'. Needs fn and f'.
     *
     * On a bracket it is Brent's method with Newton's step in place of the interpolation, safeguarded the same way:
     * x^2 - 2 on [1, 2] in 8 evaluations. From a guess it iterates until a sign change appears between two iterates,
     * then continues on that bracket; every result is then a verified sign change, at about one evaluation more than
     * the textbook stop (7 from 1 for x^2 - 2). From a guess it returns a root, the first sign change its iterates
     * cross, not necessarily the one nearest the guess, and error::not_converged where it finds none (a cycle, a
     * divergence, f' = 0). Its iteration cap is 4 times the bits of T (212 in double): beyond it lie a multiple root
     * at 0 with the default abs_tol (x^3 from 1 shrinks by 2/3 per step until x^3 underflows), a guess very far from
     * the root, and roots of multiplicity above about 5; a bracket, or abs_tol, solves them.
     */
    struct Newton
    {
    };

    /**
     * @brief Halley's method: the step -t / (1 - t f'' / (2 f')) with t = f / f'. Needs fn, f' and f''.
     *
     * Cubic near a simple root: x^2 - 2 from 1 in 5 evaluations. Where its correction cannot be formed or would point
     * against Newton's step (far from a root), it takes Newton's step, as Boost.Math's Halley does. On a bracket and
     * from a guess as Newton.
     */
    struct Halley
    {
    };

    /**
     * @brief The secant method: the step through the latest two iterates. Needs fn only; from a guess.
     *
     * Its second point is x0 (1 + eps^(1/4)) (eps^(1/4) itself at x0 = 0). At the first sign change it continues as
     * Brent's method on that bracket, so its results are verified sign changes, as Newton's are: x^2 - 2 from 1 in 9
     * evaluations. Superlinear (order 1.62) on a simple root.
     */
    struct Secant
    {
    };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief The methods solve() accepts on a bracket.
         * @tparam M The method.
         */
        template<typename M>
        concept bracketing =
            std::same_as<M, Bisection> || std::same_as<M, Brent> || std::same_as<M, Newton> || std::same_as<M, Halley>;

        /**
         * @brief The methods solve() accepts from a guess.
         * @tparam M The method.
         */
        template<typename M>
        concept open = std::same_as<M, Newton> || std::same_as<M, Halley> || std::same_as<M, Secant>;

        /**
         * @brief The number of values a method needs from fn at a point: the value, and the derivatives it uses.
         * @tparam M The method.
         */
        template<typename M>
        inline constexpr std::size_t values_needed = std::same_as<M, Newton>   ? 2
                                                     : std::same_as<M, Halley> ? 3
                                                                               : 1;

        /**
         * @brief Whether R is a std::pair, std::tuple or std::array: the tuple-likes fn may return.
         * @tparam R The type.
         */
        template<typename R>
        inline constexpr bool is_standard_tuple = false;

        /**
         * @brief A std::pair is a tuple-like fn may return.
         * @tparam A The first element type.
         * @tparam B The second element type.
         */
        template<typename A, typename B>
        inline constexpr bool is_standard_tuple<std::pair<A, B>> = true;

        /**
         * @brief A std::tuple is a tuple-like fn may return.
         * @tparam Ts The element types.
         */
        template<typename... Ts>
        inline constexpr bool is_standard_tuple<std::tuple<Ts...>> = true;

        /**
         * @brief A std::array is a tuple-like fn may return.
         * @tparam U The element type.
         * @tparam N The size.
         */
        template<typename U, std::size_t N>
        inline constexpr bool is_standard_tuple<std::array<U, N>> = true;

        /**
         * @brief Whether every element of the tuple-like R is exactly T.
         * @tparam R The tuple-like type.
         * @tparam T The real type.
         * @tparam I The indices of the elements.
         * @return true if each `std::tuple_element_t<I, R>` is T.
         */
        template<typename R, typename T, std::size_t... I>
        consteval bool elements_are(std::index_sequence<I...>)
        { return (std::same_as<std::tuple_element_t<I, R>, T> && ...); }

        /**
         * @brief A function returning, by value, a std::pair, std::tuple or std::array whose elements are all exactly
         * T: the value of fn and its derivatives, in order.
         *
         * Exactly T, not convertible to T: an element of a reference type (std::tie of fn's locals) or of an expression
         * type (Boost.Multiprecision's expression templates) could refer to fn's locals and dangle once fn returns, and
         * converting it afterwards is too late.
         *
         * @tparam F The callable type.
         * @tparam T The real type.
         */
        template<typename F, typename T>
        concept tuple_valued =
            std::invocable<F &, T const &> && !std::is_reference_v<std::invoke_result_t<F &, T const &>>
            && is_standard_tuple<std::remove_cv_t<std::invoke_result_t<F &, T const &>>>
            && detail::elements_are<std::remove_cv_t<std::invoke_result_t<F &, T const &>>, T>(
                std::make_index_sequence<std::tuple_size_v<std::remove_cv_t<std::invoke_result_t<F &, T const &>>>>{});

        /**
         * @brief fn gives at least K values at a point: its value alone (K = 1, a result convertible to T), or a
         * tuple-like of at least K elements, each exactly T. A tuple-like result is never read as a value.
         *
         * @tparam F The callable type.
         * @tparam T The real type.
         * @tparam K The number of values needed.
         */
        template<typename F, typename T, std::size_t K>
        concept gives =
            (K == 1 && nxx::detail::callable<F, T>)
            || (tuple_valued<F, T> && std::tuple_size_v<std::remove_cv_t<std::invoke_result_t<F &, T const &>>> >= K);

        /**
         * @brief The largest finite value of T.
         * @tparam T The real type.
         * @return std::numeric_limits<T>::max().
         */
        template<real T>
        constexpr T largest()
        {
            // Parenthesised, so that a max macro, as defined by <windows.h>, cannot expand.
            return T((std::numeric_limits<T>::max)());
        }

        /**
         * @brief The open methods' iteration cap: 4 times the bits of T, digits times the smallest b with
         * 2^b >= radix (212 for double, 672 for cpp_bin_float_50).
         * @tparam T The real type.
         * @return The cap.
         */
        template<real T>
        consteval int iteration_cap()
        {
            int bits_per_digit = 0;
            while ((1 << bits_per_digit) < std::numeric_limits<T>::radix) ++bits_per_digit;
            return 4 * std::numeric_limits<T>::digits * bits_per_digit;
        }

        /**
         * @brief The validated and clamped tolerances, with tol(x) = absolute + relative |x|.
         * @tparam T The real type.
         */
        template<real T>
        struct Tolerance
        {
            /** @brief The absolute tolerance, in [denorm_min, max / 4]. */
            T absolute;

            /** @brief The relative tolerance, in [2 eps, 1]. */
            T relative;

            /**
             * @brief tol(x), saturated at max.
             * @param x The point.
             * @return absolute + relative |x|, or max if that would overflow.
             */
            constexpr T at(T const &x) const
            {
                // scaled is at most |x|, since relative <= 1. Above max / 2 it saturates: max - absolute would round,
                // and could admit a sum that overflows; below it the sum is at most 3/4 max, since absolute <= max / 4.
                T const scaled = T(relative * math::abs(x));
                if (scaled > T(detail::largest<T>() / T(2))) return detail::largest<T>();
                return T(absolute + scaled);
            }
        };

        /**
         * @brief Twice @p t, saturated at max.
         * @param t A tolerance, non-negative.
         * @return 2 t, or max if that would overflow.
         */
        template<real T>
        constexpr T twice(T const &t)
        {
            if (t > T(detail::largest<T>() / T(2))) return detail::largest<T>();
            return T(T(2) * t);
        }

        /**
         * @brief Validates and clamps the caller's criteria.
         *
         * The lower clamps make every method end: Brent's minimum step, tol(b) >= 2 eps |b|, is at least two units of
         * T's spacing at b in radix 2 and 10 (at least denorm_min at 0), so it moves b. Brent requires at least
         * eps |b| (his 2 macheps for rounded arithmetic) and advises against less. The upper clamps keep tol and 2 tol
         * meaningful; the saturating helpers keep them finite.
         *
         * @param criteria The caller's criteria.
         * @return The tolerances, or error::non_finite for a tolerance that is negative or not finite.
         */
        template<real T>
        constexpr std::expected<Tolerance<T>, error> checked_tolerance(Criteria<T> const &criteria)
        {
            // An assert would vanish under NDEBUG, and a NaN tolerance would then never stop Brent: validated here.
            if (!math::isfinite(criteria.abs_tol) || !math::isfinite(criteria.rel_tol) || criteria.abs_tol < T(0)
                || criteria.rel_tol < T(0))
                return std::unexpected(error::non_finite);
            T const tiny = T(std::numeric_limits<T>::denorm_min());
            assert(tiny > T(0) && "the real type's denorm_min() must be positive");
            T const absolute_cap = T(detail::largest<T>() / T(4));
            T const relative_floor = T(T(2) * T(std::numeric_limits<T>::epsilon()));
            T absolute = criteria.abs_tol;
            if (absolute < tiny) absolute = tiny;
            if (absolute > absolute_cap) absolute = absolute_cap;
            T relative = criteria.rel_tol;
            if (relative < relative_floor) relative = relative_floor;
            if (relative > T(1)) relative = T(1);
            return Tolerance<T>{ absolute, relative };
        }

        /**
         * @brief Whether a and b have the same sign; neither may be 0.
         * @param a A nonzero value.
         * @param b A nonzero value.
         * @return true if both are negative or both positive.
         */
        template<real T>
        constexpr bool same_sign(T const &a, T const &b)
        { return (a < T(0)) == (b < T(0)); }

        /**
         * @brief An upper bound on |u - v|, never below it, formed without overflow: every abs_error, stop test and
         * guard uses it.
         *
         * - Exact where the difference is: ends of one sign (0 counts as positive) whose computed difference is below
         *   the smaller magnitude, so that the true one is too and Sterbenz's lemma applies (in radix 2 also at its
         *   boundary, one end twice the other), and with one end 0. This is every stop with the default tolerances.
         * - Otherwise the computed distance d raised by max(d eps, denorm_min), at least one unit of T's spacing at d,
         *   which covers its rounding (loose tolerances stop with ends far apart, where the difference rounds).
         * - denorm_min where ends that differ give a difference of 0: a type without subnormals, such as
         *   cpp_bin_float, flushes a difference below its smallest number to 0.
         * - max where the distance may exceed what T holds: ends of opposite signs, one of them above max / 4.
         *
         * @param u A finite point.
         * @param v A finite point.
         * @return A bound on the distance, saturated at max.
         */
        template<real T>
        constexpr T distance_bound(T const &u, T const &v)
        {
            if (u == v) return T(0);
            T const max = detail::largest<T>();
            T const tiny = T(std::numeric_limits<T>::denorm_min());
            T const au = math::abs(u);
            T const av = math::abs(v);
            if (au == T(0)) return av; // a distance from 0 is exact
            if (av == T(0)) return au;
            T distance = T(0);
            if ((u < T(0)) == (v < T(0))) {
                // Of one sign, |u - v| <= max(|u|, |v|): finite.
                distance = math::abs(T(u - v));
                if (distance == T(0)) return tiny;
                T const smaller = au < av ? au : av;
                if (distance < smaller) return distance;
                // The boundary case of Sterbenz's lemma, the larger end exactly twice the smaller, is exact too. It is
                // recognised in radix 2 only, where doubling is exact; in radix 10 it is raised, conservatively.
                if constexpr (std::numeric_limits<T>::radix == 2) {
                    T const larger = au < av ? av : au;
                    if (distance == smaller && (smaller > T(max / T(2)) || larger <= T(T(2) * smaller)))
                        return distance;
                }
            } else {
                // Of opposite signs the distance is |u| + |v|, finite (at most max / 2) when both are at most max / 4.
                T const quarter = T(max / T(4));
                if (au > quarter || av > quarter) return max;
                distance = T(au + av);
            }
            T raise = T(distance * T(std::numeric_limits<T>::epsilon()));
            if (raise < tiny) raise = tiny;
            if (raise >= T(max - distance)) return max; // max - distance is exact above max / 2 (Sterbenz)
            return T(distance + raise);
        }

        /**
         * @brief u - v, or max where |u - v| exceeds @p limit (decided without forming the difference).
         *
         * Brent's initial and reset step b - a; max then forces bisection steps until a step can be formed.
         *
         * @param u A finite point.
         * @param v A finite point.
         * @param limit The largest distance to form.
         * @return u - v, or max.
         */
        template<real T>
        constexpr T bounded_difference(T const &u, T const &v, T const &limit)
        {
            // A bound at most limit (at most max / 16) means the true distance is too: u - v is finite.
            if (detail::distance_bound(u, v) > limit) return detail::largest<T>();
            return T(u - v);
        }

        /**
         * @brief The midpoint of [a, b], without overflow.
         * @param a A finite point.
         * @param b A finite point.
         * @return (a + b) / 2 for ends of opposite signs, a + (b - a) / 2 for ends of one sign.
         */
        template<real T>
        constexpr T midpoint(T const &a, T const &b)
        {
            if ((a < T(0)) != (b < T(0))) return T(T(a + b) / T(2));
            return T(a + T(T(b - a) / T(2)));
        }

        /**
         * @brief a / b, where it can be formed: b is not 0 and the quotient is at most about max / 2.
         *
         * The test |a| / max < |b| / 2 (nxx::deriv's) is conservative by a factor 2 and forms nothing that overflows.
         *
         * @param a A finite value.
         * @param b A finite value.
         * @return a / b, or nothing.
         */
        template<real T>
        constexpr std::optional<T> quotient(T const &a, T const &b)
        {
            if (b == T(0)) return std::nullopt;
            if (!(T(math::abs(a) / detail::largest<T>()) < T(math::abs(b) / T(2)))) return std::nullopt;
            return T(a / b);
        }

        /**
         * @brief a * b, where it cannot overflow: |a| <= 1, or |b| <= (max / 2) / |a|.
         * @param a A finite value.
         * @param b A finite value.
         * @return a * b, or nothing.
         */
        template<real T>
        constexpr std::optional<T> product(T const &a, T const &b)
        {
            T const aa = math::abs(a);
            if (aa > T(1) && math::abs(b) > T(T(detail::largest<T>() / T(2)) / aa)) return std::nullopt;
            return T(a * b);
        }

        /**
         * @brief Whether x + s cannot overflow: opposite signs; or both at most max / 2; or |x| above max / 2 and |s|
         * at most max - |x| (exact there, by Sterbenz's lemma).
         * @param x A finite point.
         * @param s A finite step.
         * @return true if x + s may be formed.
         */
        template<real T>
        constexpr bool sum_fits(T const &x, T const &s)
        {
            if ((x < T(0)) != (s < T(0))) return true;
            T const half = T(detail::largest<T>() / T(2));
            T const ax = math::abs(x);
            T const as = math::abs(s);
            if (ax <= half) return as <= half;
            return as <= T(detail::largest<T>() - ax);
        }

        /**
         * @brief The values fn gives at a point, as a std::array of K values of T, converted where fn is evaluated.
         * @tparam K The number of values: 1 (the value), 2 (with f') or 3 (with f'').
         * @tparam T The real type.
         * @tparam F The callable type.
         * @tparam I The indices of the values.
         * @param fn The function.
         * @param x The point, a T const lvalue as the constraints promise.
         * @return fn's first K values at x.
         */
        template<std::size_t K, real T, typename F, std::size_t... I>
        constexpr std::array<T, K> evaluate(F &fn, T const &x, std::index_sequence<I...>)
        {
            if constexpr (tuple_valued<F, T>) {
                std::invoke_result_t<F &, T const &> const values = std::invoke(fn, x);
                return std::array<T, K>{ std::get<I>(values)... };
            } else {
                // A value: expression templates are converted to T before any arithmetic.
                return std::array<T, K>{ static_cast<T>(std::invoke(fn, x)) };
            }
        }

        /**
         * @brief The values fn gives at a point (see the overload above).
         * @tparam K The number of values.
         * @tparam T The real type.
         * @tparam F The callable type.
         * @param fn The function.
         * @param x The point.
         * @return fn's first K values at x.
         */
        template<std::size_t K, real T, typename F>
        constexpr std::array<T, K> evaluate(F &fn, T const &x)
        { return detail::evaluate<K>(fn, x, std::make_index_sequence<K>{}); }

        /**
         * @brief fn and its separate derivative callables, joined into one function that returns their values as an
         * array: the separate form becomes the tuple form, so each method has one implementation.
         *
         * Pointers rather than a lambda that captures the callables by reference: MSVC rejects such a capture of a
         * function reference in constant evaluation (as at nxx::deriv's Mirrored).
         *
         * @tparam T The real type.
         * @tparam F The function's type, possibly a function type.
         * @tparam Ds The derivatives' types.
         */
        template<real T, typename F, typename... Ds>
        struct Joined
        {
            /** @brief The function. */
            F *fn;

            /** @brief The derivative callables, f' then f''. */
            std::tuple<Ds *...> derivatives;

            /**
             * @brief The values at x: fn's, then each derivative's, converted to T.
             * @tparam I The indices of the derivatives.
             * @param x The point.
             * @return The values.
             */
            template<std::size_t... I>
            constexpr std::array<T, 1 + sizeof...(Ds)> values(T const &x, std::index_sequence<I...>) const
            {
                return std::array<T, 1 + sizeof...(Ds)>{ static_cast<T>(std::invoke(*fn, x)),
                    static_cast<T>(std::invoke(*std::get<I>(derivatives), x))... };
            }

            /**
             * @brief The values at x: fn's, then each derivative's.
             * @param x The point.
             * @return The values.
             */
            constexpr std::array<T, 1 + sizeof...(Ds)> operator()(T const &x) const
            { return values(x, std::index_sequence_for<Ds...>{}); }
        };

        /**
         * @brief A validated bracket, in order, with fn's values at its ends; or a root found during validation.
         * @tparam T The real type.
         * @tparam K The number of values at each end.
         */
        template<real T, std::size_t K>
        struct Start
        {
            /** @brief The root, if validation found one (fn exactly 0 at an end or at 0). */
            std::optional<Root<T>> root;

            /** @brief The lower end. */
            T lower;

            /** @brief fn's values at the lower end; the value finite and nonzero. */
            std::array<T, K> at_lower;

            /** @brief The upper end. */
            T upper;

            /** @brief fn's values at the upper end; the value finite, nonzero and of the other sign. */
            std::array<T, K> at_upper;
        };

        /**
         * @brief Validates the bracket, puts it in order, evaluates fn at its ends and, for a bracket around 0, at 0.
         *
         * Zero first: a root at 0 costs one evaluation, where the default tolerance would otherwise approach it down
         * to T's smallest numbers, and the rest of the method sees a bracket of one sign. Skipped where fn(0) is not
         * finite, since removable singularities at 0 (sin(x)/x, x/(e^x - 1)) are common and an evaluation the caller
         * did not ask for must not fail the solve. Only fn's value is checked; the derivatives at the ends are not.
         *
         * @tparam K The number of values the method needs.
         * @param fn The function.
         * @param bracket The caller's bracket.
         * @param tolerance The tolerances.
         * @return The start, or error::non_finite or error::no_sign_change.
         */
        template<std::size_t K, real T, typename F>
        constexpr std::expected<Start<T, K>, error>
            prepare(F &fn, Bracket<T> const &bracket, Tolerance<T> const &tolerance)
        {
            if (!math::isfinite(bracket.lower) || !math::isfinite(bracket.upper))
                return std::unexpected(error::non_finite);
            T lower = bracket.lower;
            T upper = bracket.upper;
            if (upper < lower) {
                lower = bracket.upper;
                upper = bracket.lower;
            }
            std::array<T, K> at_lower = detail::evaluate<K>(fn, lower);
            std::array<T, K> at_upper = detail::evaluate<K>(fn, upper);
            if (!math::isfinite(at_lower[0]) || !math::isfinite(at_upper[0])) return std::unexpected(error::non_finite);
            if (at_lower[0] == T(0))
                return Start<T, K>{ Root<T>{ lower, tolerance.at(lower) }, lower, at_lower, upper, at_upper };
            if (at_upper[0] == T(0))
                return Start<T, K>{ Root<T>{ upper, tolerance.at(upper) }, lower, at_lower, upper, at_upper };
            if (detail::same_sign(at_lower[0], at_upper[0])) return std::unexpected(error::no_sign_change);
            if (lower < T(0) && T(0) < upper) {
                T const zero = T(0);
                std::array<T, K> const at_zero = detail::evaluate<K>(fn, zero);
                if (at_zero[0] == T(0))
                    return Start<T, K>{ Root<T>{ zero, tolerance.at(zero) }, lower, at_lower, upper, at_upper };
                if (math::isfinite(at_zero[0])) {
                    if (detail::same_sign(at_zero[0], at_lower[0])) {
                        lower = zero;
                        at_lower = at_zero;
                    } else {
                        upper = zero;
                        at_upper = at_zero;
                    }
                }
            }
            return Start<T, K>{ std::nullopt, lower, at_lower, upper, at_upper };
        }

        /**
         * @brief Bisection on a validated bracket.
         * @param fn The function.
         * @param s The start: a bracket in order with a sign change and no root found.
         * @param tolerance The tolerances.
         * @return The root, or error::non_finite if fn is not finite at a midpoint.
         */
        template<std::size_t K, real T, typename F>
        constexpr std::expected<Root<T>, error> bisection(F &fn, Start<T, K> const &s, Tolerance<T> const &tolerance)
        {
            T a = s.lower;
            T fa = s.at_lower[0];
            T b = s.upper;
            while (true) {
                T const m = detail::midpoint(a, b);
                // The larger of the distances to the ends is a true bound even where half the width rounds, as among
                // the subnormals or with ends far apart. A midpoint not strictly inside means T's resolution.
                T const left = detail::distance_bound(m, a);
                T const right = detail::distance_bound(b, m);
                T const half = right < left ? left : right;
                if (!(a < m && m < b) || half <= tolerance.at(m)) return Root<T>{ m, half };
                T const fm = detail::evaluate<K>(fn, m)[0];
                if (!math::isfinite(fm)) return std::unexpected(error::non_finite);
                if (fm == T(0)) return Root<T>{ m, tolerance.at(m) };
                if (detail::same_sign(fm, fa)) {
                    a = m;
                    fa = fm;
                } else {
                    b = m;
                }
            }
        }

        /**
         * @brief Newton's step -f / f' from the values (f, f', ...), where it can be formed: f' nonzero and finite,
         * and the quotient representable.
         * @param v The values at the point.
         * @return The step, or nothing.
         */
        template<real T, std::size_t K>
        constexpr std::optional<T> newton_step(std::array<T, K> const &v)
        {
            static_assert(K >= 2);
            if (v[1] == T(0) || !math::isfinite(v[1])) return std::nullopt;
            std::optional<T> const ratio = detail::quotient(v[0], v[1]);
            if (!ratio) return std::nullopt;
            return T(-*ratio);
        }

        /**
         * @brief Halley's step -t / (1 - u), with t = f / f' and u = t (f'' / f') / 2, from the values (f, f', f''):
         * the form without f'^2, which overflows sooner. Where u cannot be formed, 1 - u is 0, or the step would point
         * against Newton's (u > 1, far from a root), Newton's step, as Boost.Math's Halley does.
         * @param v The values at the point.
         * @return The step, or nothing where Newton's cannot be formed either.
         */
        template<real T>
        constexpr std::optional<T> halley_step(std::array<T, 3> const &v)
        {
            std::optional<T> const newton = detail::newton_step(v);
            if (!newton) return std::nullopt;
            if (!math::isfinite(v[2])) return newton;
            std::optional<T> const curvature = detail::quotient(v[2], v[1]); // f'' / f'
            if (!curvature) return newton;
            std::optional<T> const u = detail::product(T(-*newton), T(*curvature / T(2))); // t f'' / (2 f')
            if (!u) return newton;
            T const denominator = T(T(1) - *u);
            if (denominator == T(0)) return newton;
            std::optional<T> const step = detail::quotient(*newton, denominator);
            if (!step || (*step > T(0)) != (*newton > T(0))) return newton;
            return step;
        }

        /**
         * @brief The derivative step of method M from the values at a point (Newton's or Halley's).
         * @tparam M The method, Newton or Halley.
         * @param v The values at the point.
         * @return The step, or nothing where it cannot be formed.
         */
        template<typename M, real T, std::size_t K>
        constexpr std::optional<T> derivative_step(std::array<T, K> const &v)
        {
            if constexpr (std::same_as<M, Halley>)
                return detail::halley_step(v);
            else
                return detail::newton_step(v);
        }

        /**
         * @brief Brent's procedure zero (Brent 1973, ch. 4, section 6) on a validated bracket, with the step chosen by
         * the method: Brent's interpolation, or Newton's or Halley's step from b, the end with the smaller |f|.
         *
         * Implemented from the book, with three changes: the tolerance is the caller's tol(b) rather than
         * 2 macheps |b| + t; the stop test is distance_bound(c, b) <= 2 tol(b) rather than |c - b| / 2 <= tol(b),
         * since the bound is exact there while the half may round in radix 10; and a step is a bisection step whenever
         * |c - b|, |b - a| or |e| exceeds max / 16, so that no product of the interpolation overflows (3 m q reaches 12
         * |c - b|, 2 p 4 (|c - b| + |b - a|)). The quotients cannot overflow: in the interpolation |fb| < |fa| <= |fc|.
         * A derivative step is accepted by Brent's test written for a step s: towards c, |s| < 3 |m| / 2 - tol / 2,
         * and |s| < |e| / 2; one that cannot be formed is a bisection step. The bookkeeping, and so the termination
         * argument, is Brent's for every method.
         *
         * @tparam M The method: Brent, Newton or Halley.
         * @param fn The function.
         * @param s The start: a bracket in order with a sign change and no root found.
         * @param tolerance The tolerances.
         * @return The root, or error::non_finite if fn is not finite at a point.
         */
        template<typename M, std::size_t K, real T, typename F>
        constexpr std::expected<Root<T>, error> brent(F &fn, Start<T, K> const &s, Tolerance<T> const &tolerance)
        {
            T const big = T(detail::largest<T>() / T(16));

            T a = s.lower;
            std::array<T, K> va = s.at_lower;
            T b = s.upper;
            std::array<T, K> vb = s.at_upper;
            T c = a;
            std::array<T, K> vc = va;
            T d = detail::bounded_difference(b, a, big);
            T e = d;
            while (true) {
                // [b, c] is the bracket, with fn(b) and fn(c) of opposite signs; a is the previous b.
                if (detail::same_sign(vb[0], vc[0])) {
                    c = a;
                    vc = va;
                    d = detail::bounded_difference(b, a, big);
                    e = d;
                }
                if (math::abs(vc[0]) < math::abs(vb[0])) {
                    a = b;
                    b = c;
                    c = a;
                    va = vb;
                    vb = vc;
                    vc = va;
                }
                T const fa = va[0];
                T const fb = vb[0];
                T const fc = vc[0];
                T const tol = tolerance.at(b);
                T const width = detail::distance_bound(c, b);
                if (width <= detail::twice(tol)) return Root<T>{ b, width };

                bool const wide = width > big || detail::distance_bound(b, a) > big || math::abs(e) > big;
                // Half the bracket, towards c: exact Brent below big, the overflow-free midpoint above it.
                T const m = wide ? T(detail::midpoint(b, c) - b) : T(T(c - b) / T(2));
                if (wide || math::abs(e) < tol || math::abs(fa) <= math::abs(fb)) {
                    d = m;
                    e = m;
                } else if constexpr (std::same_as<M, Brent>) {
                    T const ratio = T(fb / fa);
                    T p = T(0);
                    T q = T(0);
                    if (a == c) {
                        // Linear interpolation (the secant).
                        p = T(T(2) * m * ratio);
                        q = T(T(1) - ratio);
                    } else {
                        // Inverse quadratic interpolation.
                        T const qa = T(fa / fc);
                        T const r = T(fb / fc);
                        p = T(ratio * T(T(T(2) * m * qa * T(qa - r)) - T(T(b - a) * T(r - T(1)))));
                        q = T(T(qa - T(1)) * T(r - T(1)) * T(ratio - T(1)));
                    }
                    if (p > T(0))
                        q = T(-q);
                    else
                        p = T(-p);
                    T const before_last = e; // Brent's s := e
                    e = d;
                    // Accept the interpolated point only if it falls well inside the bracket and the steps shrink. A
                    // NaN or q = 0 fails both tests, and falls back on bisection.
                    if (T(T(2) * p) < T(T(T(3) * m * q) - math::abs(T(tol * q)))
                        && p < math::abs(T(before_last * q / T(2)))) {
                        d = T(p / q);
                    } else {
                        d = m;
                        e = m;
                    }
                } else {
                    std::optional<T> const step = detail::derivative_step<M>(vb);
                    T const before_last = e;
                    e = d;
                    // Brent's acceptance for a step: towards c, well inside the bracket (so that b + step cannot round
                    // onto c), and shorter than half the step before last.
                    T const inside = T(T(T(3) * math::abs(m)) / T(2) - T(tol / T(2)));
                    if (step && (*step > T(0)) == (m > T(0)) && math::abs(*step) < inside
                        && math::abs(*step) < T(math::abs(before_last) / T(2))) {
                        d = *step;
                    } else {
                        d = m;
                        e = m;
                    }
                }
                a = b;
                va = vb;
                if (math::abs(d) > tol)
                    b = T(b + d);
                else if (m > T(0))
                    b = T(b + tol);
                else
                    b = T(b - tol);
                vb = detail::evaluate<K>(fn, b);
                if (!math::isfinite(vb[0])) return std::unexpected(error::non_finite);
                if (vb[0] == T(0)) return Root<T>{ b, tolerance.at(b) };
            }
        }

        /**
         * @brief The bracket loop of method M on two consecutive iterates of an open method, whose values have
         * opposite signs: the hand-over. No zero first: they are adjacent iterates, not the caller's bracket.
         * @tparam M The bracketing method: Brent (for Secant), Newton or Halley.
         * @param fn The function.
         * @param x One iterate.
         * @param vx fn's values there.
         * @param y The other iterate.
         * @param vy fn's values there.
         * @param tolerance The tolerances.
         * @return The root, or error::non_finite.
         */
        template<typename M, std::size_t K, real T, typename F>
        constexpr std::expected<Root<T>, error> hand_over(F &fn,
            T const &x,
            std::array<T, K> const &vx,
            T const &y,
            std::array<T, K> const &vy,
            Tolerance<T> const &tolerance)
        {
            if (x < y) return detail::brent<M>(fn, Start<T, K>{ std::nullopt, x, vx, y, vy }, tolerance);
            return detail::brent<M>(fn, Start<T, K>{ std::nullopt, y, vy, x, vx }, tolerance);
        }

        /**
         * @brief Newton or Halley from a guess: the method's step, at least tol(x) long, until fn's values at two
         * consecutive iterates have opposite signs; then the bracket loop on them (the hand-over).
         *
         * The minimum step takes its direction from the signs of f and f', not from a step that may have underflowed to
         * 0. An x + step that would overflow, or a step that cannot be formed, is error::not_converged.
         *
         * @tparam M The method, Newton or Halley.
         * @param fn The function, giving the values M needs.
         * @param guess The starting point.
         * @param tolerance The tolerances.
         * @return The root, error::non_finite or error::not_converged.
         */
        template<typename M, real T, typename F>
        constexpr std::expected<Root<T>, error> open_derivative(F &fn, T const &guess, Tolerance<T> const &tolerance)
        {
            constexpr std::size_t K = values_needed<M>;
            if (!math::isfinite(guess)) return std::unexpected(error::non_finite);
            T x = guess;
            std::array<T, K> vx = detail::evaluate<K>(fn, x);
            if (!math::isfinite(vx[0])) return std::unexpected(error::non_finite);
            if (vx[0] == T(0)) return Root<T>{ x, tolerance.at(x) };
            for (int i = 0; i < detail::iteration_cap<T>(); ++i) {
                std::optional<T> const step = detail::derivative_step<M>(vx);
                if (!step) return std::unexpected(error::not_converged);
                T const tol = tolerance.at(x);
                bool const upward = !detail::same_sign(vx[0], vx[1]); // Newton's direction, -f / f'
                T const move = math::abs(*step) > tol ? *step : (upward ? tol : T(-tol));
                if (!detail::sum_fits(x, move)) return std::unexpected(error::not_converged);
                T const next = T(x + move);
                std::array<T, K> const vn = detail::evaluate<K>(fn, next);
                if (!math::isfinite(vn[0])) return std::unexpected(error::non_finite);
                if (vn[0] == T(0)) return Root<T>{ next, tolerance.at(next) };
                if (!detail::same_sign(vn[0], vx[0])) return detail::hand_over<M>(fn, x, vx, next, vn, tolerance);
                x = next;
                vx = vn;
            }
            return std::unexpected(error::not_converged);
        }

        /**
         * @brief Secant's second point: x0 (1 + eta), eta = eps^(1/4); eta at x0 = 0; x0 (1 - eta) where x0 (1 + eta)
         * would overflow. Relative, so no scale is assumed.
         * @param x0 The guess, finite.
         * @return The second point.
         */
        template<real T>
        constexpr T secant_second_point(T const &x0)
        {
            T const eta = math::root_epsilon<4, T>();
            if (x0 == T(0)) return eta;
            if (math::abs(x0) > T(detail::largest<T>() / T(T(1) + eta))) return T(x0 * T(T(1) - eta));
            return T(x0 * T(T(1) + eta));
        }

        /**
         * @brief Secant from a guess: the step (x - xp) (f / (fp - f)) through the latest two iterates, at least tol(x)
         * long, until fn's values at two consecutive iterates have opposite signs; then Brent on them.
         *
         * fp and f have one sign before the hand-over, so fp - f cannot overflow; the ratio first keeps a large f from
         * overflowing the product (1e300 (x - 1) from 1e7). The minimum step takes its direction from the signs of f,
         * fp - f and x - xp.
         *
         * @param fn The function.
         * @param guess The starting point.
         * @param tolerance The tolerances.
         * @return The root, error::non_finite or error::not_converged.
         */
        template<real T, typename F>
        constexpr std::expected<Root<T>, error> open_secant(F &fn, T const &guess, Tolerance<T> const &tolerance)
        {
            if (!math::isfinite(guess)) return std::unexpected(error::non_finite);
            T xp = guess;
            std::array<T, 1> vp = detail::evaluate<1>(fn, xp);
            if (!math::isfinite(vp[0])) return std::unexpected(error::non_finite);
            if (vp[0] == T(0)) return Root<T>{ xp, tolerance.at(xp) };
            T x = detail::secant_second_point(xp);
            T dx = T(x - xp); // small: x0 eta, or eta
            std::array<T, 1> vx = detail::evaluate<1>(fn, x);
            for (int i = 0; i < detail::iteration_cap<T>(); ++i) {
                if (!math::isfinite(vx[0])) return std::unexpected(error::non_finite);
                if (vx[0] == T(0)) return Root<T>{ x, tolerance.at(x) };
                if (!detail::same_sign(vx[0], vp[0])) return detail::hand_over<Brent>(fn, xp, vp, x, vx, tolerance);
                T const difference = T(vp[0] - vx[0]);
                if (difference == T(0)) return std::unexpected(error::not_converged);
                std::optional<T> const ratio = detail::quotient(vx[0], difference);
                if (!ratio) return std::unexpected(error::not_converged);
                std::optional<T> const step = detail::product(dx, *ratio);
                if (!step) return std::unexpected(error::not_converged);
                T const tol = tolerance.at(x);
                // The step's sign is that of dx * f / (fp - f), taken from the signs in case the step underflowed.
                bool const upward = ((dx > T(0)) == (vx[0] > T(0))) == (difference > T(0));
                T const move = math::abs(*step) > tol ? *step : (upward ? tol : T(-tol));
                if (!detail::sum_fits(x, move)) return std::unexpected(error::not_converged);
                xp = x;
                vp = vx;
                x = T(x + move);
                dx = move;
                vx = detail::evaluate<1>(fn, x);
            }
            return std::unexpected(error::not_converged);
        }

    } // namespace detail
    /** @endcond */

    /**
     * @brief A root of @p fn in @p bracket by the method @p M: Bisection, Brent, or Newton or Halley safeguarded by the
     * bracket.
     *
     * The bracket's ends may come in either order; fn must have opposite signs at them, or be 0 at one. When the
     * bracket contains 0 in its interior, fn is also evaluated at 0 first: a root at 0 is then found at once, and the
     * bracket is halved at 0 (skipped where fn(0) is not finite). The method then narrows the bracket until it is at
     * most 2 tol(x) wide (see Criteria), or fn is exactly 0 at a point. Newton and Halley take their derivatives from
     * fn, which then returns a std::pair, std::tuple or std::array of exactly T (the value, f', f''); the other methods
     * accept such a fn too and use its first element.
     *
     * A pole inside the bracket changes sign like a root, and is returned as one: where fn may have poles, check
     * |fn(value)| (a later version may report a pole as an error). fn's exceptions pass through unchanged.
     *
     * @code
     * namespace rt = nxx::roots;
     * auto const r = rt::solve<rt::Brent>([](double t) { return psat(t) - 101325.0; }, rt::Bracket{300.0, 400.0});
     * auto const s = rt::solve<rt::Bisection>(fn, rt::Bracket{0.0, 1.0}, {.rel_tol = 1e-8});
     * auto const n = rt::solve<rt::Newton>([](double x) { return std::pair{x * x - 2.0, 2.0 * x}; },
     *                                      rt::Bracket{1.0, 2.0});
     * @endcode
     *
     * @tparam M The method: `Bisection`, `Brent`, `Newton` or `Halley`.
     * @tparam T The real type, deduced from @p bracket.
     * @tparam F A callable taking T const & and returning something convertible to T, or (for Newton and Halley, which
     *           need it) a std::pair, std::tuple or std::array of exactly T: the value and the derivatives. For Newton
     *           on a fn without a derivative, pass the derivative as a callable after fn, use
     *           nxx::deriv::derivative_of, or use Secant or Brent.
     * @param fn The function.
     * @param bracket An interval in which fn changes sign.
     * @param criteria When to stop; it converts to T, so `{.rel_tol = 1e-10}` works for any T.
     * @return The root and a bound on its distance to a sign change of fn; error::non_finite if an end, a tolerance or
     *         a value of fn is not finite or a tolerance is negative; error::no_sign_change if fn has the same sign at
     *         both ends.
     */
    // criteria is a non-deduced context (type_identity_t), so T comes from the bracket alone and a braced list of
    // doubles converts to a multiprecision T.
    template<typename M, real T, typename F>
    requires detail::bracketing<M> && detail::gives<F, T, detail::values_needed<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        constexpr std::size_t K = detail::values_needed<M>;
        // Qualified, so that argument-dependent lookup cannot pick a function from the namespace of T or F.
        std::expected<detail::Tolerance<T>, error> const tolerance = detail::checked_tolerance(criteria);
        if (!tolerance) return std::unexpected(tolerance.error());
        std::expected<detail::Start<T, K>, error> const start = detail::prepare<K>(fn, bracket, *tolerance);
        if (!start) return std::unexpected(start.error());
        if (start->root) return *start->root;
        if constexpr (std::same_as<M, Bisection>)
            return detail::bisection(fn, *start, *tolerance);
        else
            return detail::brent<M>(fn, *start, *tolerance);
    }

    /**
     * @brief A root of @p fn from @p guess by the open method @p M: Newton, Halley or Secant.
     *
     * The method iterates from the guess until fn's values at two consecutive iterates have opposite signs, then
     * continues inside that bracket (Newton and Halley safeguarded, Secant as Brent), so the result is a verified sign
     * change with a bound, as for the bracketing methods. It returns a root, the first sign change its iterates cross,
     * not necessarily the one nearest the guess; from a guess already within tol of a root it takes one step past it (2
     * evaluations). error::not_converged where it finds no sign change within its cap (4 times the bits of T: 212
     * iterations in double) or cannot take a step: a cycle, a divergence, f' = 0, a multiple root at 0 with the default
     * abs_tol. A bracket solves those.
     *
     * @code
     * namespace rt = nxx::roots;
     * auto const fdf = [](double x) { return std::pair{x * x - 2.0, 2.0 * x}; };
     * auto const n = rt::solve<rt::Newton>(fdf, 1.0);                              // 7 evaluations
     * auto const s = rt::solve<rt::Secant>([](double x) { return x * x - 2.0; }, 1.0); // 9
     * @endcode
     *
     * @tparam M The method: `Newton`, `Halley` or `Secant`.
     * @tparam T The real type, deduced from @p guess.
     * @tparam F A callable taking T const &: for Secant returning something convertible to T, for Newton and Halley a
     *           std::pair, std::tuple or std::array of exactly T (the value, f', f''). For Newton on a fn without a
     *           derivative, pass the derivative as a callable after fn, use nxx::deriv::derivative_of, or use Secant.
     * @param fn The function.
     * @param guess The starting point.
     * @param criteria When to stop; it converts to T.
     * @return The root and a bound on its distance to a sign change of fn; error::non_finite if the guess, a tolerance
     *         or a value of fn is not finite or a tolerance is negative; error::not_converged if the method found no
     *         root.
     */
    template<typename M, real T, typename F>
    requires detail::open<M> && detail::gives<F, T, detail::values_needed<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        std::expected<detail::Tolerance<T>, error> const tolerance = detail::checked_tolerance(criteria);
        if (!tolerance) return std::unexpected(tolerance.error());
        if constexpr (std::same_as<M, Secant>)
            return detail::open_secant(fn, guess, *tolerance);
        else
            return detail::open_derivative<M>(fn, guess, *tolerance);
    }

    /**
     * @brief Newton on @p bracket with the derivative as a separate callable: the same as Newton on a fn returning
     * (f, f'), with each callable evaluated once per point.
     * @tparam M Newton.
     * @tparam T The real type, deduced from @p bracket.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @tparam D The derivative: a callable as F.
     * @param fn The function.
     * @param dfn Its derivative.
     * @param bracket An interval in which fn changes sign.
     * @param criteria When to stop.
     * @return As for a bracket.
     */
    template<typename M, real T, typename F, typename D>
    requires std::same_as<M, Newton> && nxx::detail::callable<F, T> && nxx::detail::callable<D, T>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, D &&dfn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const joined{ std::addressof(fn),
            { std::addressof(dfn) } };
        return nxx::roots::solve<M>(joined, bracket, criteria);
    }

    /**
     * @brief Newton from @p guess with the derivative as a separate callable.
     * @tparam M Newton.
     * @tparam T The real type, deduced from @p guess.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @tparam D The derivative: a callable as F.
     * @param fn The function.
     * @param dfn Its derivative.
     * @param guess The starting point.
     * @param criteria When to stop.
     * @return As from a guess.
     */
    template<typename M, real T, typename F, typename D>
    requires std::same_as<M, Newton> && nxx::detail::callable<F, T> && nxx::detail::callable<D, T>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, D &&dfn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const joined{ std::addressof(fn),
            { std::addressof(dfn) } };
        return nxx::roots::solve<M>(joined, guess, criteria);
    }

    /**
     * @brief Halley on @p bracket with the first and second derivatives as separate callables.
     * @tparam M Halley.
     * @tparam T The real type, deduced from @p bracket.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @tparam D The first derivative: a callable as F.
     * @tparam D2 The second derivative: a callable as F.
     * @param fn The function.
     * @param dfn Its first derivative.
     * @param d2fn Its second derivative.
     * @param bracket An interval in which fn changes sign.
     * @param criteria When to stop.
     * @return As for a bracket.
     */
    template<typename M, real T, typename F, typename D, typename D2>
    requires std::same_as<M, Halley> && nxx::detail::callable<F, T> && nxx::detail::callable<D, T>
             && nxx::detail::callable<D2, T>
    [[nodiscard]] constexpr std::expected<Root<T>, error> solve(F &&fn,
        D &&dfn,
        D2 &&d2fn,
        Bracket<T> const &bracket,
        Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>, std::remove_reference_t<D2>> const
            joined{ std::addressof(fn), { std::addressof(dfn), std::addressof(d2fn) } };
        return nxx::roots::solve<M>(joined, bracket, criteria);
    }

    /**
     * @brief Halley from @p guess with the first and second derivatives as separate callables.
     * @tparam M Halley.
     * @tparam T The real type, deduced from @p guess.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @tparam D The first derivative: a callable as F.
     * @tparam D2 The second derivative: a callable as F.
     * @param fn The function.
     * @param dfn Its first derivative.
     * @param d2fn Its second derivative.
     * @param guess The starting point.
     * @param criteria When to stop.
     * @return As from a guess.
     */
    template<typename M, real T, typename F, typename D, typename D2>
    requires std::same_as<M, Halley> && nxx::detail::callable<F, T> && nxx::detail::callable<D, T>
             && nxx::detail::callable<D2, T>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, D &&dfn, D2 &&d2fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>, std::remove_reference_t<D2>> const
            joined{ std::addressof(fn), { std::addressof(dfn), std::addressof(d2fn) } };
        return nxx::roots::solve<M>(joined, guess, criteria);
    }

} // namespace nxx::roots

/** @brief The standard library's formatter, specialised for nxx::roots::error. */
namespace std {

    /**
     * @brief Formats an nxx::roots::error as the name of its enumerator, such as `no_sign_change`, with the format
     * specification of a string.
     *
     * @code
     * if (!r) std::println("no root: {}", r.error());   // no root: no_sign_change
     * @endcode
     */
    template<>
    struct formatter<nxx::roots::error> : formatter<string_view>
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
        auto format(nxx::roots::error e, FormatContext &ctx) const
        {
            // No default: a new enumerator makes the compilers warn here (-Wswitch) until its name is added.
            string_view name = "unknown";
            switch (e) {
            case nxx::roots::error::non_finite:
                name = "non_finite";
                break;
            case nxx::roots::error::no_sign_change:
                name = "no_sign_change";
                break;
            case nxx::roots::error::not_converged:
                name = "not_converged";
                break;
            }
            return formatter<string_view>::format(name, ctx);
        }
    };

} // namespace std
