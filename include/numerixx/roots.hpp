/**
 * @file
 * @brief Roots of a real function of one variable: bracketing methods and open methods.
 *
 * solve<M>(fn, start) finds a point where fn changes sign, by a method M. From a bracket, an interval at whose ends fn
 * has opposite signs: `Bisection`, which halves it, `Brent`, which interpolates and falls back on bisection,
 * `Toms748` and `AndersonBjorck`, two more interpolating methods, `Ridders`, `RegulaFalsi` and `Illinois`, kept as
 * teaching examples, and `Newton` and `Halley`, whose derivative steps the bracket safeguards. From a guess:
 * `Newton`, `Halley`, `Secant` and `Steffensen`,
 * which iterate until a sign change appears between two iterates and then finish inside it. Every method returns a
 * std::expected<Root<T>, error>: the root and a bound on its distance to a sign change of fn. Numerixx itself never
 * throws, but the result's own value() throws std::bad_expected_access on a failure, so test the result and read it
 * with `->` or `*`.
 *
 * Newton and Halley take fn's derivatives either from fn, which then returns a std::pair, std::tuple or std::array of
 * exactly T (the value, f', f''), or from separate callables passed after fn, which may return a std::expected:
 * nxx::deriv::derivative_of(fn) gives a numerical derivative as such a callable.
 *
 * search(fn, x0) finds a bracket from a guess, by an outward search, for a bracketing method to finish:
 * `search(fn, x0).and_then([&](Bracket<double> b) { return solve<Brent>(fn, b); })`.
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

        /**
         * fn has the same sign at both ends of the bracket, and is 0 at neither; or search() found no sign change
         * between the points it evaluated, within its cap or before both its sides ended.
         */
        no_sign_change,

        /**
         * An open method found no root from its guess: no sign change between its iterates within its iteration cap,
         * or no usable step (f' = 0, not finite or a failed std::expected, Secant's two values equal, a step that would
         * overflow). With a numerical derivative also where its stencil leaves fn's domain at an iterate, or where |f'|
         * is as small as the derivative's own error (a multiple root): see Newton with separate derivatives.
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
     * @brief Alefeld, Potra and Shi's Algorithm 748 (1995, Algorithm 4.2): inverse cubic interpolation with quadratic
     * and secant fallbacks, a double-length secant step, and a bisection whenever an iteration has not halved the
     * bracket. A bracket only; no derivative form.
     *
     * Costs about what Brent costs on smooth fn (the paper reports 4% fewer evaluations; on 27 of its test cases at
     * full precision Brent took 309, Toms748 325), with a better worst case: at most about 4 evaluations per halving
     * of the bracket, where Brent's bound is quadratic in the number of halvings. Like every bracketing method it
     * returns a point and a bound, not the final bracket (Boost.Math's toms748_solve returns the bracket).
     */
    struct Toms748
    {
    };

    /**
     * @brief Regula falsi with Anderson and Bjorck's rescaling (1973): the secant through the latest point and the
     * retained end, whose value is scaled down by 1 - f(c) / f(b) (or 1/2) whenever the same end is retained again. A
     * bracket only.
     *
     * Guarded by a bisection after three steps that have not halved the bracket, so at most about 4 evaluations per
     * halving (the published method took 65 million evaluations on x^20 - 1 over [0, 5]; guarded, 30). Far cheaper than
     * Brent where fn is strongly curved near one end of the bracket (ln(x / 1e-12) on [1e-20, 1]: 18 evaluations against
     * 45) and on a walk to a sign change at 0; dearer on high powers and very wide brackets.
     */
    struct AndersonBjorck
    {
    };

    /**
     * @brief Ridders' method (1979): the midpoint, then the zero of an exponential fit through the three values. A
     * bracket. Not nxx::deriv::Ridders, the derivative method of the same author.
     *
     * Each iteration halves the bracket at a cost of 2 evaluations (so at most twice Bisection's count) and converges
     * quadratically per iteration on a simple root. A teaching example; the cheapest method here on roots of odd
     * multiplicity and on flat fn ((x - 1)^3: 86 evaluations against Brent's 156), dearer than Brent otherwise (416
     * against 309 on 27 of Alefeld, Potra and Shi's test cases). The fit's zero is formed without the cancellation of
     * the textbook formula, which loses roots near an end of the bracket.
     */
    struct Ridders
    {
    };

    /**
     * @brief False position: the zero of the secant through the bracket's ends, guarded. A bracket.
     *
     * A teaching example. As published, one end of the bracket stalls while the other converges: x^20 - 1 over [0, 5]
     * took over 200,000 evaluations. Guarded here, as AndersonBjorck, by a bisection after three steps that have not
     * halved the bracket (45 evaluations there), so the stall shows only in its counts: 645 on 27 of Alefeld, Potra and
     * Shi's test cases, against Illinois' 382 and AndersonBjorck's 327.
     */
    struct RegulaFalsi
    {
    };

    /**
     * @brief The Illinois method (Dowell and Jarratt 1971): false position with the retained end's value halved
     * whenever that end is retained again, the simplest cure for its stall; guarded, as AndersonBjorck. A bracket.
     */
    struct Illinois
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

    /**
     * @brief Steffensen's method: Newton's with the slope (f(z) - f(x)) / f(x) at z = x + f(x). Needs fn only; from a
     * guess, 2 evaluations per iteration, at the first sign change continuing as Brent's method, as Secant.
     *
     * A teaching example: quadratic without a derivative, at the price of scale dependence. f(x) is added to x, so fn's
     * scale sets the cost and even which root is found: 1000 (x^2 - 2) from 1 returns -sqrt 2 after 27 evaluations,
     * where x^2 - 2 returns sqrt 2 after 10 (z = 1 - 1000 brackets the other root). From a far guess it is slow (x^2 - 2
     * from 10: 39 evaluations, Secant 13). Prefer Secant.
     */
    struct Steffensen
    {
    };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief The methods solve() accepts on a bracket.
         * @tparam M The method.
         */
        template<typename M>
        concept bracketing = std::same_as<M, Bisection> || std::same_as<M, Brent> || std::same_as<M, Toms748>
                             || std::same_as<M, AndersonBjorck> || std::same_as<M, Ridders>
                             || std::same_as<M, RegulaFalsi> || std::same_as<M, Illinois> || std::same_as<M, Newton>
                             || std::same_as<M, Halley>;

        /**
         * @brief The methods solve() accepts from a guess.
         * @tparam M The method.
         */
        template<typename M>
        concept open = std::same_as<M, Newton> || std::same_as<M, Halley> || std::same_as<M, Secant>
                       || std::same_as<M, Steffensen>;

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
            && !std::is_volatile_v<std::invoke_result_t<F &, T const &>>
            && is_standard_tuple<std::remove_cv_t<std::invoke_result_t<F &, T const &>>>
            && detail::elements_are<std::remove_cv_t<std::invoke_result_t<F &, T const &>>, T>(
                std::make_index_sequence<std::tuple_size_v<std::remove_cv_t<std::invoke_result_t<F &, T const &>>>>{});

        /**
         * @brief A function returning one value: its result converts to T and is not a std::pair, std::tuple or
         * std::array (a tuple-like result is never read as a value, even where T converts from it).
         *
         * @tparam F The callable type.
         * @tparam T The real type.
         */
        template<typename F, typename T>
        concept value_valued = nxx::detail::callable<F, T>
                               && !is_standard_tuple<std::remove_cvref_t<std::invoke_result_t<F &, T const &>>>;

        /**
         * @brief Whether R is a std::expected.
         * @tparam R The type.
         */
        template<typename R>
        inline constexpr bool is_expected = false;

        /**
         * @brief A std::expected is one.
         * @tparam U The value type.
         * @tparam E The error type.
         */
        template<typename U, typename E>
        inline constexpr bool is_expected<std::expected<U, E>> = true;

        /**
         * @brief A separate derivative callable: one returning a value (value_valued), or a std::expected whose value
         * converts to T, such as nxx::deriv::derivative_of's callable. The expected's error is never read.
         *
         * The value is the expression derivative_value converts, `*r` on the result as returned (by value: an rvalue;
         * by reference: with its reference and constness), so that a conversion the constraint admits is the one the
         * body performs (independent review of Phase 3).
         *
         * @tparam D The callable type.
         * @tparam T The real type.
         */
        // is_expected comes first, so that *r is formed only for an expected.
        template<typename D, typename T>
        concept derivative_valued =
            value_valued<D, T>
            || (std::invocable<D &, T const &> && is_expected<std::remove_cvref_t<std::invoke_result_t<D &, T const &>>>
                && std::convertible_to<decltype(*std::declval<std::invoke_result_t<D &, T const &>>()), T>);

        /**
         * @brief fn gives at least K values at a point: its value alone (K = 1, value_valued), or a tuple-like of at
         * least K elements, each exactly T.
         *
         * @tparam F The callable type.
         * @tparam T The real type.
         * @tparam K The number of values needed.
         */
        template<typename F, typename T, std::size_t K>
        concept gives =
            (K == 1 && value_valued<F, T>)
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
            // In 64 bits: an admitted radix may exceed 2^30, where an int shift would overflow.
            long long bits_per_digit = 0;
            long long reach = 1;
            while (reach < std::numeric_limits<T>::radix) {
                reach *= 2;
                ++bits_per_digit;
            }
            long long const cap = 4LL * std::numeric_limits<T>::digits * bits_per_digit;
            return cap > (std::numeric_limits<int>::max)() ? (std::numeric_limits<int>::max)() : static_cast<int>(cap);
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
         * - max where the distance may exceed what T holds: ends of opposite signs whose halves sum to more than
         *   max / 2. Ends of opposite signs above max / 4 are measured by halves, raised by two units (in radix 10 the
         *   halvings and the doubling may round).
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
            T units = T(1);
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
                // Of opposite signs the distance is |u| + |v|: finite (at most max / 2) when both are at most max / 4;
                // above that, by halves. Saturating there instead left a bracket around 0 under huge tolerances wider
                // than its true width forever (independent review of Phase 2).
                T const quarter = T(max / T(4));
                if (au <= quarter && av <= quarter) {
                    distance = T(au + av);
                } else {
                    T const half_sum = T(T(au / T(2)) + T(av / T(2)));
                    if (half_sum > T(max / T(2))) return max;
                    distance = T(T(2) * half_sum);
                    units = T(2);
                }
            }
            T raise = T(T(distance * T(std::numeric_limits<T>::epsilon())) * units);
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
         * Formed directly where |b| >= 1 (then |a / b| <= |a|); otherwise only if |a| < (max / 2) |b|, which neither
         * overflows nor underflows to 0 for |b| < 1. The one-line test |a| / max < |b| / 2 rejected representable
         * quotients where both sides underflow to 0 (d (x - 1) / d with d = denorm_min; independent review).
         *
         * @param a A finite value.
         * @param b A finite value.
         * @return a / b, or nothing.
         */
        template<real T>
        constexpr std::optional<T> quotient(T const &a, T const &b)
        {
            if (b == T(0)) return std::nullopt;
            T const ab = math::abs(b);
            if (ab >= T(1)) return T(a / b);
            if (!(math::abs(a) < T(T(detail::largest<T>() / T(2)) * ab))) return std::nullopt;
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
         * @brief Whether x + s cannot overflow: opposite signs; or both at most max / 2; or the smaller magnitude at
         * most max less the larger (exact there, by Sterbenz's lemma, since the larger exceeds max / 2). The test is on
         * the magnitudes, not on x's: 0 + max fits (a first form tested |s| <= max / 2 wherever |x| was, and rejected
         * Steffensen's first point on max - x from 0; independent review of Phase 6).
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
            T const larger = ax < as ? as : ax;
            T const smaller = ax < as ? ax : as;
            if (larger <= half) return true;
            return smaller <= T(detail::largest<T>() - larger);
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
         * @brief A separate derivative's value at x: its value converted to T, or 0 for a std::expected without one.
         *
         * A failed derivative counts as one that is not finite, and 0 has the same effect in both steps: f' = 0 makes
         * Newton's step unusable, as f' = inf does, and f'' = 0 gives exactly Newton's step in Halley's, as f'' = inf
         * does. 0, because every nxx::real has it, where a NaN or an infinity need not exist.
         *
         * @tparam T The real type.
         * @tparam D The callable type.
         * @param dfn The derivative.
         * @param x The point.
         * @return The value.
         */
        template<real T, typename D>
        constexpr T derivative_value(D &dfn, T const &x)
        {
            if constexpr (value_valued<D, T>) {
                return static_cast<T>(std::invoke(dfn, x));
            } else {
                // Called once, and the result kept as returned: a copy would require a copyable expected, and a const
                // copy would refuse a value that converts only as an rvalue (independent review of Phase 3).
                decltype(auto) r = std::invoke(dfn, x);
                if (!r) return T(0);
                return static_cast<T>(*std::forward<decltype(r)>(r));
            }
        }

        /**
         * @brief fn and its separate derivative callables, joined into one function that returns their values as an
         * array: the separate form becomes the tuple form, so each method has one implementation.
         *
         * Pointers rather than a lambda that captures the callables by reference: MSVC rejects such a capture of a
         * function reference in constant evaluation (as at nxx::deriv's Mirrored).
         *
         * @tparam T The real type.
         * @tparam F The function's type, possibly a function type.
         * @tparam Ds The derivatives' types: each returns a value or a std::expected (derivative_valued).
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
                    detail::derivative_value(*std::get<I>(derivatives), x)... };
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
                T next = b;
                if (math::abs(d) > tol)
                    next = T(b + d);
                else if (m > T(0))
                    next = T(b + tol);
                else
                    next = T(b - tol);
                // A point not strictly inside the bracket means it is already at the tolerance's resolution: no step
                // can make progress, and evaluating there again would repeat the iteration forever.
                if (!((b < next && next < c) || (c < next && next < b))) return Root<T>{ b, width };
                a = b;
                va = vb;
                b = next;
                vb = detail::evaluate<K>(fn, b);
                if (!math::isfinite(vb[0])) return std::unexpected(error::non_finite);
                if (vb[0] == T(0)) return Root<T>{ b, tolerance.at(b) };
            }
        }

        /**
         * @brief x - y, where it cannot overflow (sum_fits).
         * @param x A finite value.
         * @param y A finite value.
         * @return x - y, or nothing.
         */
        template<real T>
        constexpr std::optional<T> difference(T const &x, T const &y)
        {
            if (!detail::sum_fits(x, T(-y))) return std::nullopt;
            return T(x - y);
        }

        /**
         * @brief x + y, where both exist and the sum cannot overflow.
         * @param x A value, or nothing.
         * @param y A value, or nothing.
         * @return x + y, or nothing.
         */
        template<real T>
        constexpr std::optional<T> sum(std::optional<T> const &x, std::optional<T> const &y)
        {
            if (!x || !y || !detail::sum_fits(*x, *y)) return std::nullopt;
            return T(*x + *y);
        }

        /**
         * @brief The result where the bracket [a, b] meets the stop test: the end with the smaller |f|, with the bound
         * on the bracket's width.
         * @param a The lower end.
         * @param fa fn there.
         * @param b The upper end.
         * @param fb fn there.
         * @param tolerance The tolerances.
         * @return The root, or nothing to go on.
         */
        template<real T>
        constexpr std::optional<Root<T>>
            stop_at(T const &a, T const &fa, T const &b, T const &fb, Tolerance<T> const &tolerance)
        {
            T const u = math::abs(fa) < math::abs(fb) ? a : b;
            T const width = detail::distance_bound(a, b);
            if (width <= detail::twice(tolerance.at(u))) return Root<T>{ u, width };
            return std::nullopt;
        }

        /**
         * @brief The state of Toms748 and AndersonBjorck: the bracket a < b with fn's values, and the ends dropped by
         * the last two steps, d and then e (Toms748's interpolation points).
         * @tparam T The real type.
         */
        template<real T>
        struct Enclosure
        {
            /** @brief The lower end. */
            T a;

            /** @brief fn at a. */
            T fa;

            /** @brief The upper end. */
            T b;

            /** @brief fn at b; of the other sign than fa. */
            T fb;

            /** @brief The end the last step dropped. */
            T d;

            /** @brief fn at d. */
            T fd;

            /** @brief The d before it. */
            T e;

            /** @brief fn at e. */
            T fe;
        };

        /**
         * @brief One step of enclose(): the root where the search ends, the new state, and the point evaluated.
         * @tparam T The real type.
         */
        template<real T>
        struct Enclosed
        {
            /** @brief The root, where the step found fn exactly 0 or the bracket met the stop test. */
            std::optional<Root<T>> root;

            /** @brief The state after the step. */
            Enclosure<T> next;

            /** @brief The point evaluated. */
            T c;

            /** @brief fn at c. */
            T fc;
        };

        /**
         * @brief Whether the bracket is wider than max / 16, where every step is a bisection step and no other point
         * is formed (Brent's guard): every difference of points is then at most max / 4.
         * @param s The state.
         * @return true while the bracket is that wide.
         */
        template<real T>
        constexpr bool wide(Enclosure<T> const &s)
        { return detail::distance_bound(s.a, s.b) > T(detail::largest<T>() / T(16)); }

        /**
         * @brief Where a proposed point is evaluated: the paper's modified bracket (Alefeld, Potra and Shi 1995, p.
         * 341) with our tolerance.
         *
         * delta = 0.7 tol(u), u the end with the smaller |f|; at least denorm_min, since in a type whose denorm_min()
         * is min() (Boost.Multiprecision's) 0.7 min flushes to 0; max / 4 where tol(u) exceeds max / 2. The midpoint
         * where the bracket is at most 4 delta wide, where there is no point, and where the point kept 2 delta inside
         * the ends is not strictly inside (2 delta rounds away at an end of much larger magnitude than u). A point at
         * least 2 delta inside either ends the search or removes at least 1.4 tol.
         *
         * @param s The state.
         * @param proposed The method's point, or nothing.
         * @param tolerance The tolerances.
         * @return The point.
         */
        template<real T>
        constexpr T place(Enclosure<T> const &s, std::optional<T> const &proposed, Tolerance<T> const &tolerance)
        {
            T const u = math::abs(s.fa) < math::abs(s.fb) ? s.a : s.b;
            T const t = tolerance.at(u);
            T delta = T(detail::largest<T>() / T(4));
            if (t <= T(detail::largest<T>() / T(2))) {
                delta = T(T(T(7) / T(10)) * t);
                T const tiny = T(std::numeric_limits<T>::denorm_min());
                if (delta < tiny) delta = tiny;
            }
            T const middle = detail::midpoint(s.a, s.b);
            if (!proposed || T(detail::distance_bound(s.a, s.b) / T(4)) <= delta) return middle;
            T c = *proposed;
            T const margin = T(T(2) * delta);
            if (c <= T(s.a + margin))
                c = T(s.a + margin);
            else if (c >= T(s.b - margin))
                c = T(s.b - margin);
            if (s.a < c && c < s.b) return c;
            return middle;
        }

        /**
         * @brief One step of Toms748 and AndersonBjorck, the paper's bracket subroutine: the point (the midpoint while
         * the bracket is wide, else the proposed point placed as place() says), e := d where asked, evaluate fn, keep
         * the half with the sign change and record the dropped end as d, and stop on an exact zero or the stop test.
         * @param fn The function.
         * @param s The state.
         * @param proposed The method's point, computed only where the bracket is not wide; or nothing.
         * @param shift Whether e takes d's value before d changes.
         * @param tolerance The tolerances.
         * @return The step, or error::non_finite if fn is not finite at the point.
         */
        template<real T, typename F>
        constexpr std::expected<Enclosed<T>, error> enclose(F &fn,
            Enclosure<T> const &s,
            std::optional<T> const &proposed,
            bool shift,
            Tolerance<T> const &tolerance)
        {
            T const c = detail::wide(s) ? detail::midpoint(s.a, s.b) : detail::place(s, proposed, tolerance);
            // Unreachable: a midpoint not strictly inside means adjacent ends, which met the stop test before.
            assert(s.a < c && c < s.b);
            if (!(s.a < c && c < s.b))
                return Enclosed<T>{ Root<T>{ c, detail::distance_bound(s.a, s.b) }, s, c, s.fa };
            Enclosure<T> next = s;
            if (shift) {
                next.e = s.d;
                next.fe = s.fd;
            }
            T const fc = detail::evaluate<1>(fn, c)[0];
            if (!math::isfinite(fc)) return std::unexpected(error::non_finite);
            if (fc == T(0)) return Enclosed<T>{ Root<T>{ c, tolerance.at(c) }, next, c, fc };
            if (detail::same_sign(s.fa, fc)) {
                next.d = s.a;
                next.fd = s.fa;
                next.a = c;
                next.fa = fc;
            } else {
                next.d = s.b;
                next.fd = s.fb;
                next.b = c;
                next.fb = fc;
            }
            return Enclosed<T>{ detail::stop_at(next.a, next.fa, next.b, next.fb, tolerance), next, c, fc };
        }

        /**
         * @brief The secant's zero through (a, fa) and (b, fb), fa and fb of opposite signs.
         * @param a The lower end.
         * @param fa fn at a.
         * @param b The upper end.
         * @param fb fn at b.
         * @return The point, or nothing where fb - fa would overflow.
         */
        template<real T>
        constexpr std::optional<T> secant_point(T const &a, T const &fa, T const &b, T const &fb)
        {
            std::optional<T> const spread = detail::difference(fb, fa);
            if (!spread) return std::nullopt;
            std::optional<T> const r = detail::quotient(fa, *spread);
            if (!r) return std::nullopt;
            return T(a - T(*r * T(b - a))); // |r| < 1, and b - a at most max / 16
        }

        /**
         * @brief Newton-Quadratic (Alefeld, Potra and Shi 1995, p. 330): k Newton steps on the quadratic through (a,
         * fa), (b, fb), (d, fd), from the end where it is convex towards the root; the secant where its leading
         * coefficient is 0, where a quantity cannot be formed, or as soon as an iterate leaves (a, b) (in exact
         * arithmetic they never do).
         * @param s The state; d is outside [a, b].
         * @param k The number of Newton steps, 2 or 3.
         * @return The point, or nothing.
         */
        template<real T>
        constexpr std::optional<T> newton_quadratic(Enclosure<T> const &s, int k)
        {
            std::optional<T> const secant = detail::secant_point(s.a, s.fa, s.b, s.fb);
            std::optional<T> const rise = detail::difference(s.fb, s.fa);
            std::optional<T> const rise_d = detail::difference(s.fd, s.fb);
            if (!rise || !rise_d) return secant;
            std::optional<T> const b1 = detail::quotient(*rise, T(s.b - s.a)); // f[a, b]
            std::optional<T> const bd = detail::quotient(*rise_d, T(s.d - s.b)); // f[b, d]
            if (!b1 || !bd) return secant;
            std::optional<T> const curve = detail::difference(*bd, *b1);
            if (!curve) return secant;
            std::optional<T> const a2 = detail::quotient(*curve, T(s.d - s.a)); // f[a, b, d]
            if (!a2 || *a2 == T(0)) return secant;
            T c = detail::same_sign(*a2, s.fa) ? s.a : s.b;
            for (int i = 0; i < k; ++i) {
                // P(c) = fa + (B + A (c - b)) (c - a), P'(c) = B + A ((c - a) + (c - b)): 2c - a - b would overflow
                // above max / 2.
                std::optional<T> const ab = detail::product(*a2, T(c - s.b));
                std::optional<T> const slope = ab ? detail::sum(b1, ab) : std::nullopt;
                std::optional<T> const term = slope ? detail::product(*slope, T(c - s.a)) : std::nullopt;
                std::optional<T> const value = term ? detail::sum(std::optional<T>(s.fa), term) : std::nullopt;
                std::optional<T> const along = detail::product(*a2, T(T(c - s.a) + T(c - s.b)));
                std::optional<T> const derivative = along ? detail::sum(b1, along) : std::nullopt;
                std::optional<T> const step =
                    value && derivative ? detail::quotient(*value, *derivative) : std::nullopt;
                if (!step || !detail::sum_fits(c, T(-*step))) return secant;
                T const next = T(c - *step);
                if (!(s.a < next && next < s.b)) return secant;
                c = next;
            }
            return c;
        }

        /**
         * @brief ipzero (Alefeld, Potra and Shi 1995, pp. 333-334): the zero of the inverse cubic through fn's values
         * at a, b, d and e, by Aitken-Neville; each term a distance times a ratio of values, the ratio first.
         * @param s The state; fn's four values distinct.
         * @return The point, or nothing where a quantity cannot be formed or the point is not inside (a, b).
         */
        template<real T>
        constexpr std::optional<T> inverse_cubic(Enclosure<T> const &s)
        {
            auto const term = [](std::optional<T> const &dx, T const &f, std::optional<T> const &df) -> std::optional<T> {
                if (!dx || !df) return std::nullopt;
                std::optional<T> const r = detail::quotient(f, *df);
                if (!r) return std::nullopt;
                return detail::product(*dx, *r);
            };
            auto const minus = [](std::optional<T> const &x, std::optional<T> const &y) -> std::optional<T> {
                if (!x || !y) return std::nullopt;
                return detail::difference(*x, *y);
            };
            std::optional<T> const q11 = term(T(s.d - s.e), s.fd, detail::difference(s.fe, s.fd));
            std::optional<T> const q21 = term(T(s.b - s.d), s.fb, detail::difference(s.fd, s.fb));
            std::optional<T> const q31 = term(T(s.a - s.b), s.fa, detail::difference(s.fb, s.fa));
            std::optional<T> const d21 = term(T(s.b - s.d), s.fd, detail::difference(s.fd, s.fb));
            std::optional<T> const d31 = term(T(s.a - s.b), s.fb, detail::difference(s.fb, s.fa));
            std::optional<T> const q22 = term(minus(d21, q11), s.fb, detail::difference(s.fe, s.fb));
            std::optional<T> const q32 = term(minus(d31, q21), s.fa, detail::difference(s.fd, s.fa));
            std::optional<T> const d32 = term(minus(d31, q21), s.fd, detail::difference(s.fd, s.fa));
            std::optional<T> const q33 = term(minus(d32, q22), s.fa, detail::difference(s.fe, s.fa));
            std::optional<T> const c = detail::sum(std::optional<T>(s.a), detail::sum(detail::sum(q31, q32), q33));
            if (!c || !(s.a < *c && *c < s.b)) return std::nullopt;
            return c;
        }

        /**
         * @brief Steps 4.2.3 (k = 2) and 4.2.5 (k = 3): the inverse cubic where fn's values at a, b, d and e are
         * distinct and its point is inside, otherwise Newton-Quadratic with k steps.
         * @param s The state.
         * @param k 2 or 3.
         * @param quadratic Whether only Newton-Quadratic may be used (the first iteration, n = 2, before e exists).
         * @return The point, or nothing.
         */
        template<real T>
        constexpr std::optional<T> interpolate(Enclosure<T> const &s, int k, bool quadratic)
        {
            bool const distinct = s.fa != s.fb && s.fa != s.fd && s.fa != s.fe && s.fb != s.fd && s.fb != s.fe
                                  && s.fd != s.fe;
            if (!quadratic && distinct) {
                std::optional<T> const c = detail::inverse_cubic(s);
                if (c) return c;
            }
            return detail::newton_quadratic(s, k);
        }

        /**
         * @brief Steps 4.2.7-4.2.9: twice the secant step from the end u with the smaller |f|, or nothing (the
         * midpoint) where that moves more than half the bracket.
         * @param s The state.
         * @return The point, or nothing.
         */
        template<real T>
        constexpr std::optional<T> double_secant(Enclosure<T> const &s)
        {
            bool const at_a = math::abs(s.fa) < math::abs(s.fb);
            T const u = at_a ? s.a : s.b;
            T const fu = at_a ? s.fa : s.fb;
            std::optional<T> const spread = detail::difference(s.fb, s.fa);
            if (!spread) return std::nullopt;
            std::optional<T> const r = detail::quotient(fu, *spread); // |r| <= 1/2: fa and fb have opposite signs
            if (!r) return std::nullopt;
            T const width = T(s.b - s.a);
            T const c = T(u - T(T(T(2) * *r) * width));
            if (math::abs(T(c - u)) > T(width / T(2))) return std::nullopt;
            return c;
        }

        /**
         * @brief A step of Toms748 with its point, computed only where the bracket is not wide.
         *
         * A function rather than a lambda capturing fn: MSVC rejects a lambda's capture of a function reference in
         * constant evaluation (as at Joined).
         *
         * @param fn The function.
         * @param s The state.
         * @param point The method's point from the state, or nothing.
         * @param shift Whether e takes d's value.
         * @param tolerance The tolerances.
         * @return The step, or error::non_finite.
         */
        template<real T, typename F, typename P>
        constexpr std::expected<Enclosed<T>, error>
            advance(F &fn, Enclosure<T> const &s, P const &point, bool shift, Tolerance<T> const &tolerance)
        {
            std::optional<T> const proposed = detail::wide(s) ? std::optional<T>() : std::optional<T>(point(s));
            return detail::enclose(fn, s, proposed, shift, tolerance);
        }

        /**
         * @brief Alefeld, Potra and Shi's Algorithm 4.2 (1995, pp. 334-335) on a validated bracket, each point taken
         * through enclose(): a secant step; then per iteration two interpolation steps (inverse cubic, or
         * Newton-Quadratic with 2 and then 3 Newton steps), a double-length secant step, and a bisection step where the
         * iteration has not brought the bracket below half its width at its start. At most 4 evaluations per iteration,
         * each iteration at least halving the bracket.
         * @param fn The function.
         * @param start The start: a bracket in order with a sign change and no root found.
         * @param tolerance The tolerances.
         * @return The root, or error::non_finite if fn is not finite at a point.
         */
        template<real T, typename F>
        constexpr std::expected<Root<T>, error> toms748(F &fn, Start<T, 1> const &start, Tolerance<T> const &tolerance)
        {
            T const lower = start.lower;
            T const upper = start.upper;
            T const f_lower = start.at_lower[0];
            T const f_upper = start.at_upper[0];
            if (std::optional<Root<T>> const r = detail::stop_at(lower, f_lower, upper, f_upper, tolerance)) return *r;
            // d and e are read only once steps have set them (the first iteration uses Newton-Quadratic, without e).
            Enclosure<T> s{ lower, f_lower, upper, f_upper, upper, f_upper, upper, f_upper };
            std::expected<Enclosed<T>, error> r = detail::advance(
                fn, s, [](Enclosure<T> const &e) { return detail::secant_point(e.a, e.fa, e.b, e.fb); }, false,
                tolerance); // 4.2.1-4.2.2
            for (bool first = true;; first = false) {
                if (!r) return std::unexpected(r.error());
                if (r->root) return *r->root;
                s = r->next;
                T const a0 = s.a;
                T const b0 = s.b;
                // 4.2.3-4.2.4: e := d before the bracket.
                r = detail::advance(
                    fn, s, [first](Enclosure<T> const &e) { return detail::interpolate(e, 2, first); }, true, tolerance);
                if (!r) return std::unexpected(r.error());
                if (r->root) return *r->root;
                // 4.2.5-4.2.6.
                r = detail::advance(
                    fn, r->next, [](Enclosure<T> const &e) { return detail::interpolate(e, 3, false); }, false,
                    tolerance);
                if (!r) return std::unexpected(r.error());
                if (r->root) return *r->root;
                // 4.2.7-4.2.10, e := d-bar (4.1.9).
                r = detail::advance(
                    fn, r->next, [](Enclosure<T> const &e) { return detail::double_secant(e); }, true, tolerance);
                if (!r) return std::unexpected(r.error());
                if (r->root) return *r->root;
                // 4.2.11: a bisection where the iteration has not halved the bracket, e := d-hat.
                if (detail::distance_bound(r->next.a, r->next.b) < T(detail::distance_bound(a0, b0) / T(2))) continue;
                r = detail::advance(
                    fn, r->next, [](Enclosure<T> const &e) { return std::optional<T>(detail::midpoint(e.a, e.b)); },
                    true, tolerance);
            }
        }

        /**
         * @brief False position on a validated bracket, the retained end's value scaled by M's factor, guarded by a
         * bisection after three consecutive steps that have not brought the bracket below half its width at the last
         * halving (a guard bisection counts as one); each point taken through enclose().
         *
         * b is the latest point and a the retained end, whose value is replaced by g_a. The point is the secant's
         * zero through (a, g_a) and (b, f_b). Where f_c has f_b's sign, a is retained and g_a scaled by M's factor:
         * AndersonBjorck's m = 1 - f_c / f_b, or 1/2 where m is not positive (formed only where |f_c| < |f_b|, so
         * nothing overflows; Anderson and Bjorck 1973); Illinois' 1/2; RegulaFalsi's 1. Otherwise b becomes the
         * retained end with its true value. c becomes b. Nothing depends on the factor being below 1: g_a keeps its
         * sign, so the step's ratio stays in (0, 1).
         *
         * @tparam M AndersonBjorck, Illinois or RegulaFalsi.
         * @param fn The function.
         * @param start The start: a bracket in order with a sign change and no root found.
         * @param tolerance The tolerances.
         * @return The root, or error::non_finite if fn is not finite at a point.
         */
        template<typename M, real T, typename F>
        constexpr std::expected<Root<T>, error>
            false_position(F &fn, Start<T, 1> const &start, Tolerance<T> const &tolerance)
        {
            T retained = start.lower;
            T f_retained = start.at_lower[0];
            T scaled = f_retained;
            T latest = start.upper;
            T f_latest = start.at_upper[0];
            if (std::optional<Root<T>> const r = detail::stop_at(retained, f_retained, latest, f_latest, tolerance))
                return *r;
            T reference = detail::distance_bound(retained, latest);
            int unhalved = 0;
            while (true) {
                bool const ordered = retained < latest;
                Enclosure<T> const s{ ordered ? retained : latest, ordered ? f_retained : f_latest,
                    ordered ? latest : retained, ordered ? f_latest : f_retained, retained, f_retained, retained,
                    f_retained };
                bool const bisect = unhalved >= 3;
                std::optional<T> proposed;
                if (!detail::wide(s)) {
                    if (bisect) {
                        proposed = detail::midpoint(s.a, s.b);
                    } else if (std::optional<T> const spread = detail::difference(f_latest, scaled)) {
                        // f_latest and scaled have opposite signs, so 0 < r < 1.
                        if (std::optional<T> const r = detail::quotient(f_latest, *spread))
                            proposed = T(latest - T(*r * T(latest - retained)));
                    }
                }
                std::expected<Enclosed<T>, error> const step = detail::enclose(fn, s, proposed, false, tolerance);
                if (!step) return std::unexpected(step.error());
                if (step->root) return *step->root;
                if (!detail::same_sign(step->fc, f_latest)) {
                    retained = latest;
                    f_retained = f_latest;
                    scaled = f_latest;
                } else if constexpr (std::same_as<M, AndersonBjorck>) {
                    T const m = math::abs(step->fc) < math::abs(f_latest) ? T(T(1) - T(step->fc / f_latest))
                                                                          : T(T(1) / T(2));
                    scaled = T(scaled * m);
                } else if constexpr (std::same_as<M, Illinois>) {
                    scaled = T(scaled * T(T(1) / T(2)));
                } // RegulaFalsi: the factor 1
                latest = step->c;
                f_latest = step->fc;
                T const width = detail::distance_bound(retained < latest ? retained : latest,
                    retained < latest ? latest : retained);
                if (bisect || width <= T(reference / T(2))) {
                    reference = width;
                    unhalved = 0;
                } else {
                    ++unhalved;
                }
            }
        }

        /**
         * @brief The square root of x >= 0, faithful (within one ulp), from + - * /, the same at compile time and run time.
         *
         * x is scaled into [1/radix^2, 1] by even powers of the radix, squared up so that the scaling takes a number of
         * steps logarithmic in the exponent (single factors would take about 10^9 in cpp_bin_float_50), then Newton's
         * iteration runs from 1, above the root, in the correction form y + (s/y - y)/2, as root_epsilon, until rounding
         * stops it, and the result is unscaled exactly. The textbook form (y + s/y)/2 is up to 3 ulps off in radix 10,
         * where y + s/y falls in the decade above; the correction form stays within one (exhaustive in decimal at 4 to 6
         * digits). Correct rounding would need the exact residual x - y^2, which Dekker's splitting gives in radix 2 but
         * not for every nxx::real.
         *
         * @param x The value, not negative.
         * @return sqrt x; x itself for 0, +inf and NaN, and for a negative x (a precondition, asserted), so that the
         *         loops end.
         */
        template<real T>
        constexpr T sqrt(T const &x)
        {
            assert(!(x < T(0)));
            if (!(x > T(0)) || !math::isfinite(x)) return x;
            T const r = static_cast<T>(std::numeric_limits<T>::radix);
            T const top = detail::largest<T>();
            T s = x;
            T scale = T(1);
            // A type too narrow to hold radix^2 (a Boost type whose max is 3.98) is not scaled, and Newton starts at
            // max(1, s), above the root (independent review of Phase 6: radix^2 overflowed, and sqrt 2 came out 1.9e-6).
            bool const scalable = r <= T(top / r);
            T const r2 = scalable ? T(r * r) : r;
            // p = q^2 is squared up while it fits and keeps s on its side of the interval; s is scaled by p, scale by q.
            while (scalable && s > T(1)) {
                T p = r2;
                T q = r;
                while (p <= T(top / p) && T(s / p) > p) {
                    p = T(p * p);
                    q = T(q * q);
                }
                s = T(s / p);
                scale = T(scale * q);
            }
            while (scalable && s < T(T(1) / r2)) {
                T p = r2;
                T q = r;
                while (p <= T(top / p) && T(s * p) < T(T(1) / p)) {
                    p = T(p * p);
                    q = T(q * q);
                }
                s = T(s * p);
                scale = T(scale / q);
            }
            T y = s > T(1) ? s : T(1);
            while (true) {
                T const next = T(y + T(T(T(s / y) - y) / T(2)));
                if (!(next < y)) break;
                y = next;
            }
            return T(y * scale);
        }

        /**
         * @brief Ridders' point: the zero of the exponential fit through the old bracket's ends and its midpoint m,
         * from the state after the midpoint step, whose bracket has m at one end and e at the other (e's value has the
         * other sign than m's), and the dropped end d.
         *
         * x4 = m + (m - a) sign(f(a) - f(b)) f(m) / sqrt(f(m)^2 - f(a) f(b)) (Ridders 1979), written as
         * x4 = e + w (m - e), w = 1 - |f(m)| / root = -g_e g_d / (root (root + |g_m|)), with g = f / s for
         * s = max(|f(e)|, |f(d)|, |f(m)|), so that no square overflows (root is in [|g_m|, sqrt 2]). The textbook form
         * rounds |f(m)| / root to 1 where f(e) f(d) is below eps of f(m)^2, which put the point on the end: x - 1e-300
         * over [0, 1] bisected through 1922 evaluations, against 4.
         *
         * @param s The state after the midpoint step.
         * @param m The midpoint.
         * @param fm fn there.
         * @return The point, or nothing where w cannot be formed (root underflows).
         */
        template<real T>
        constexpr std::optional<T> ridders_point(Enclosure<T> const &s, T const &m, T const &fm)
        {
            bool const m_is_a = s.a == m;
            T const e = m_is_a ? s.b : s.a;
            T const fe = m_is_a ? s.fb : s.fa;
            T scale = math::abs(fe);
            if (math::abs(s.fd) > scale) scale = math::abs(s.fd);
            if (math::abs(fm) > scale) scale = math::abs(fm);
            T const ge = T(fe / scale);
            T const gd = T(s.fd / scale);
            T const gm = T(fm / scale);
            T const root = detail::sqrt(T(T(gm * gm) - T(ge * gd)));
            std::optional<T> const w = detail::quotient(T(-T(ge * gd)), T(root * T(root + math::abs(gm))));
            if (!w) return std::nullopt;
            return T(e + T(*w * T(m - e)));
        }

        /**
         * @brief Ridders' method on a validated bracket: per iteration enclose() at the midpoint, then at Ridders' point
         * from the state that step leaves. Each iteration halves the bracket, at most 2 evaluations.
         * @param fn The function.
         * @param start The start: a bracket in order with a sign change and no root found.
         * @param tolerance The tolerances.
         * @return The root, or error::non_finite if fn is not finite at a point.
         */
        template<real T, typename F>
        constexpr std::expected<Root<T>, error> ridders(F &fn, Start<T, 1> const &start, Tolerance<T> const &tolerance)
        {
            T const lower = start.lower;
            T const upper = start.upper;
            T const f_lower = start.at_lower[0];
            T const f_upper = start.at_upper[0];
            if (std::optional<Root<T>> const r = detail::stop_at(lower, f_lower, upper, f_upper, tolerance)) return *r;
            Enclosure<T> s{ lower, f_lower, upper, f_upper, upper, f_upper, upper, f_upper };
            while (true) {
                // The midpoint: enclose's point where none is proposed.
                std::expected<Enclosed<T>, error> const half = detail::enclose(fn, s, std::optional<T>(), false, tolerance);
                if (!half) return std::unexpected(half.error());
                if (half->root) return *half->root;
                std::optional<T> const proposed =
                    detail::wide(half->next) ? std::optional<T>() : detail::ridders_point(half->next, half->c, half->fc);
                std::expected<Enclosed<T>, error> const fit =
                    detail::enclose(fn, half->next, proposed, false, tolerance);
                if (!fit) return std::unexpected(fit.error());
                if (fit->root) return *fit->root;
                s = fit->next;
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
         * would overflow; x0 + tol(x0) where x0 (1 + eta) rounds back to x0 (among the subnormals). Relative, so no
         * scale is assumed.
         * @param x0 The guess, finite.
         * @param tolerance The tolerances.
         * @return The second point, distinct from x0.
         */
        template<real T>
        constexpr T secant_second_point(T const &x0, Tolerance<T> const &tolerance)
        {
            T const eta = math::root_epsilon<4, T>();
            if (x0 == T(0)) return eta;
            if (math::abs(x0) > T(detail::largest<T>() / T(T(1) + eta))) return T(x0 * T(T(1) - eta));
            T const second = T(x0 * T(T(1) + eta));
            if (second == x0) return T(x0 + tolerance.at(x0)); // tiny x0: tol(x0) is about denorm_min
            return second;
        }

        /**
         * @brief The secant's step from x: dx (f / (fp - f)), dx = x - xp, through (x, f) and (xp, fp) of one sign, at
         * least tol(x) long, its direction from the signs of f, fp - f and dx in case the step underflowed. Secant's and
         * Steffensen's.
         * @param x The latest point.
         * @param f fn there.
         * @param dx x - xp.
         * @param difference fp - f, of one sign's values (it cannot overflow).
         * @param tolerance The tolerances.
         * @return The move from x, or nothing where it cannot be formed (fp = f, a quotient or product out of range, or
         *         x + move overflowing).
         */
        template<real T>
        constexpr std::optional<T>
            secant_move(T const &x, T const &f, T const &dx, T const &difference, Tolerance<T> const &tolerance)
        {
            if (difference == T(0)) return std::nullopt;
            std::optional<T> const ratio = detail::quotient(f, difference);
            if (!ratio) return std::nullopt;
            std::optional<T> const step = detail::product(dx, *ratio);
            if (!step) return std::nullopt;
            T const tol = tolerance.at(x);
            bool const upward = ((dx > T(0)) == (f > T(0))) == (difference > T(0));
            T const move = math::abs(*step) > tol ? *step : (upward ? tol : T(-tol));
            if (!detail::sum_fits(x, move)) return std::nullopt;
            return move;
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
            T x = detail::secant_second_point(xp, tolerance);
            T dx = T(x - xp); // small: x0 eta, or eta
            std::array<T, 1> vx = detail::evaluate<1>(fn, x);
            for (int i = 0;; ++i) {
                // The checks come first, so that the last point evaluated before the cap is checked too.
                if (!math::isfinite(vx[0])) return std::unexpected(error::non_finite);
                if (vx[0] == T(0)) return Root<T>{ x, tolerance.at(x) };
                if (!detail::same_sign(vx[0], vp[0])) return detail::hand_over<Brent>(fn, xp, vp, x, vx, tolerance);
                if (i == detail::iteration_cap<T>()) return std::unexpected(error::not_converged);
                std::optional<T> const step = detail::secant_move(x, vx[0], dx, T(vp[0] - vx[0]), tolerance);
                if (!step) return std::unexpected(error::not_converged);
                T const move = *step;
                xp = x;
                vp = vx;
                x = T(x + move);
                dx = move;
                vx = detail::evaluate<1>(fn, x);
            }
        }

        /**
         * @brief Steffensen from a guess: z = x + f(x), at least tol(x) away from x in f(x)'s direction; the secant's
         * step through (x, f(x)) and (z, f(z)) (Steffensen's -f(x)^2 / (f(z) - f(x)) where z - x = f(x)); Brent at the
         * first sign change, between x and z or x and the next iterate.
         *
         * z - x is taken as formed where it is exact (|z - x| <= |x|, by Sterbenz), else as the intended offset, where
         * the difference could overflow.
         *
         * @param fn The function.
         * @param guess The starting point.
         * @param tolerance The tolerances.
         * @return The root, error::non_finite or error::not_converged.
         */
        template<real T, typename F>
        constexpr std::expected<Root<T>, error> open_steffensen(F &fn, T const &guess, Tolerance<T> const &tolerance)
        {
            if (!math::isfinite(guess)) return std::unexpected(error::non_finite);
            T x = guess;
            std::array<T, 1> vx = detail::evaluate<1>(fn, x);
            if (!math::isfinite(vx[0])) return std::unexpected(error::non_finite);
            if (vx[0] == T(0)) return Root<T>{ x, tolerance.at(x) };
            for (int i = 0; i < detail::iteration_cap<T>(); ++i) {
                T const tol = tolerance.at(x);
                T const offset = math::abs(vx[0]) > tol ? vx[0] : (vx[0] > T(0) ? tol : T(-tol));
                if (!detail::sum_fits(x, offset)) return std::unexpected(error::not_converged);
                T const z = T(x + offset);
                std::array<T, 1> const vz = detail::evaluate<1>(fn, z);
                if (!math::isfinite(vz[0])) return std::unexpected(error::non_finite);
                if (vz[0] == T(0)) return Root<T>{ z, tolerance.at(z) };
                if (!detail::same_sign(vz[0], vx[0])) return detail::hand_over<Brent>(fn, x, vx, z, vz, tolerance);
                T const width = math::abs(offset) <= math::abs(x) ? T(z - x) : offset;
                std::optional<T> const move = detail::secant_move(x, vx[0], T(-width), T(vz[0] - vx[0]), tolerance);
                if (!move) return std::unexpected(error::not_converged);
                T const next = T(x + *move);
                std::array<T, 1> const vn = detail::evaluate<1>(fn, next);
                if (!math::isfinite(vn[0])) return std::unexpected(error::non_finite);
                if (vn[0] == T(0)) return Root<T>{ next, tolerance.at(next) };
                if (!detail::same_sign(vn[0], vx[0])) return detail::hand_over<Brent>(fn, x, vx, next, vn, tolerance);
                x = next;
                vx = vn;
            }
            return std::unexpected(error::not_converged);
        }

    } // namespace detail
    /** @endcond */

    /**
     * @brief A root of @p fn in @p bracket by the method @p M: a bracketing method (Bisection, Brent, Toms748,
     * AndersonBjorck, Ridders, RegulaFalsi or Illinois), or Newton or Halley safeguarded by the bracket.
     *
     * The bracket's ends may come in either order; fn must have opposite signs at them, or be 0 at one. When the
     * bracket contains 0 in its interior, fn is also evaluated at 0 first: a root at 0 is then found at once, and the
     * bracket is halved at 0 (skipped where fn(0) is not finite). The method then narrows the bracket until it is at
     * most 2 tol(x) wide (see Criteria), or fn is exactly 0 at a point. Newton and Halley take their derivatives from
     * fn, which then returns a std::pair, std::tuple or std::array of exactly T (the value, f', f''); the other methods
     * accept such a fn too and use its first element.
     *
     * A pole inside the bracket changes sign like a root, and is returned as one: where fn may have poles, check
     * |fn(value)| (a later version may report a pole as an error); a point that lands exactly on a pole, where fn is
     * infinite, is error::non_finite. fn's exceptions pass through unchanged.
     *
     * Choosing a method: Brent by default. Toms748 costs about the same, with a worst case of about 4 evaluations per
     * halving of the bracket where Brent's grows with the square of the halvings. AndersonBjorck where fn is strongly
     * curved near one end of the bracket (a logarithm, 1/x) or the sign change is at 0. Bisection where fn is rough
     * and its count must be predictable. Ridders where the root has odd multiplicity or fn is flat ((x - 1)^3: 86
     * evaluations against Brent's 156). RegulaFalsi, Illinois and Steffensen are teaching examples. For a value-only
     * fn, not a derivative method with a numerical derivative.
     * Each returns a point and a bound on its distance to the sign change, not the final bracket.
     *
     * @code
     * namespace rt = nxx::roots;
     * auto const r = rt::solve<rt::Brent>([](double t) { return psat(t) - 101325.0; }, rt::Bracket{300.0, 400.0});
     * auto const s = rt::solve<rt::Bisection>(fn, rt::Bracket{0.0, 1.0}, {.rel_tol = 1e-8});
     * auto const n = rt::solve<rt::Newton>([](double x) { return std::pair{x * x - 2.0, 2.0 * x}; },
     *                                      rt::Bracket{1.0, 2.0});
     * @endcode
     *
     * @tparam M The method: `Bisection`, `Brent`, `Toms748`, `AndersonBjorck`, `Ridders`, `RegulaFalsi`, `Illinois`,
     *           `Newton` or `Halley`.
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
        else if constexpr (std::same_as<M, Toms748>)
            return detail::toms748(fn, *start, *tolerance);
        else if constexpr (std::same_as<M, AndersonBjorck> || std::same_as<M, Illinois> || std::same_as<M, RegulaFalsi>)
            return detail::false_position<M>(fn, *start, *tolerance);
        else if constexpr (std::same_as<M, Ridders>)
            return detail::ridders(fn, *start, *tolerance);
        else
            return detail::brent<M>(fn, *start, *tolerance);
    }

    /**
     * @brief A root of @p fn from @p guess by the open method @p M: Newton, Halley, Secant or Steffensen.
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
     * @tparam M The method: `Newton`, `Halley`, `Secant` or `Steffensen`.
     * @tparam T The real type, deduced from @p guess.
     * @tparam F A callable taking T const &: for Secant and Steffensen returning something convertible to T, for Newton
     *           and Halley a std::pair, std::tuple or std::array of exactly T (the value, f', f''). For Newton on a fn
     *           without a derivative, pass the derivative as a callable after fn, use nxx::deriv::derivative_of, or use
     *           Secant.
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
        else if constexpr (std::same_as<M, Steffensen>)
            return detail::open_steffensen(fn, guess, *tolerance);
        else
            return detail::open_derivative<M>(fn, guess, *tolerance);
    }

    /**
     * @brief Newton on @p bracket with the derivative as a separate callable: the same as Newton on a fn returning
     * (f, f'), with each callable evaluated once per point.
     *
     * The derivative may return a value, or a std::expected whose value converts to T (its error is never read), such
     * as the callable nxx::deriv::derivative_of(fn) or a lambda around nxx::deriv::derivative(fn, x, h) with a step of
     * one's own. A std::expected without a value counts as a derivative that is not finite: on a bracket the step is a
     * bisection step, from a guess the method ends with error::not_converged, and Halley takes Newton's step where only
     * f'' failed. A numerical derivative fails where fn is not finite at a node or the stencil cannot be formed (a
     * point near max, a step below the resolution of x); a fn that stays finite outside its domain gives a finite,
     * wrong derivative instead, so return NaN there. Only the value is read: nxx::deriv::Ridders' error estimate does
     * not steer the search.
     *
     * A numerical derivative is a convenience. With nxx::deriv's default `Central<2>` a point costs 3 calls of fn for
     * Newton and 6 for Halley, about twice Secant from a guess and two to three times Brent on a bracket for a
     * value-only fn; with nxx::deriv::Ridders, 6 to 11 times Secant. The bound never depends on the derivative, but from a guess
     * the iterates can: where |f'| is as small as the derivative's own error, at a multiple root or a stationary point,
     * Newton can stall into error::not_converged ((x - 1)^3 from 2) or jump far. A bracket always converges.
     *
     * @code
     * namespace rt = nxx::roots;
     * auto const r = rt::solve<rt::Newton>(fn, nxx::deriv::derivative_of(fn), rt::Bracket{0.0, 1.0});
     * auto const s = rt::solve<rt::Newton>(fn, [&](double x) { return nxx::deriv::derivative(fn, x, 1e-4); }, 1.0);
     * // f with its numerical derivative from one callable, evaluating f(x) once:
     * auto const d = nxx::deriv::derivative_of(fn);
     * auto const fd = [&](double x) {
     *     auto const df = d(x);
     *     return std::pair{fn(x), df ? df->value : std::numeric_limits<double>::quiet_NaN()};
     * };
     * @endcode
     *
     * A misuse fails one of the constraints, named in the compiler's message: detail::value_valued (fn returns a
     * value, not a tuple or a std::expected), detail::derivative_valued (the derivatives), detail::bracketing and
     * detail::open (the method and the start).
     * @tparam M Newton.
     * @tparam T The real type, deduced from @p bracket.
     * @tparam F A callable taking T const & and returning something convertible to T.
     * @tparam D The derivative: a callable as F, or one returning a std::expected whose value converts to T.
     * @param fn The function.
     * @param dfn Its derivative.
     * @param bracket An interval in which fn changes sign.
     * @param criteria When to stop.
     * @return As for a bracket.
     */
    template<typename M, real T, typename F, typename D>
    requires std::same_as<M, Newton> && detail::value_valued<F, T> && detail::derivative_valued<D, T>
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
     * @tparam D The derivative: a callable as F, or one returning a std::expected whose value converts to T.
     * @param fn The function.
     * @param dfn Its derivative.
     * @param guess The starting point.
     * @param criteria When to stop.
     * @return As from a guess.
     */
    template<typename M, real T, typename F, typename D>
    requires std::same_as<M, Newton> && detail::value_valued<F, T> && detail::derivative_valued<D, T>
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
     * @tparam D The first derivative: a callable as F, or one returning a std::expected whose value converts to T.
     * @tparam D2 The second derivative: as D.
     * @param fn The function.
     * @param dfn Its first derivative.
     * @param d2fn Its second derivative.
     * @param bracket An interval in which fn changes sign.
     * @param criteria When to stop.
     * @return As for a bracket.
     */
    template<typename M, real T, typename F, typename D, typename D2>
    requires std::same_as<M, Halley> && detail::value_valued<F, T> && detail::derivative_valued<D, T>
             && detail::derivative_valued<D2, T>
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
     * @tparam D The first derivative: a callable as F, or one returning a std::expected whose value converts to T.
     * @tparam D2 The second derivative: as D.
     * @param fn The function.
     * @param dfn Its first derivative.
     * @param d2fn Its second derivative.
     * @param guess The starting point.
     * @param criteria When to stop.
     * @return As from a guess.
     */
    template<typename M, real T, typename F, typename D, typename D2>
    requires std::same_as<M, Halley> && detail::value_valued<F, T> && detail::derivative_valued<D, T>
             && detail::derivative_valued<D2, T>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, D &&dfn, D2 &&d2fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>, std::remove_reference_t<D2>> const
            joined{ std::addressof(fn), { std::addressof(dfn), std::addressof(d2fn) } };
        return nxx::roots::solve<M>(joined, guess, criteria);
    }

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief One side of search(): growing away from x0, or probing towards fn's domain edge once a value was not
         * finite.
         * @tparam T The real type.
         */
        template<real T>
        struct SearchSide
        {
            /** @brief -1 below x0, +1 above. */
            int direction;

            /** @brief Whether the side still searches. */
            bool alive;

            /** @brief The side's last point where fn was finite, with fn(x0)'s sign. */
            T last;

            /** @brief The nearest point where fn was not finite, once one was met: the side probes towards it. */
            std::optional<T> edge;

            /** @brief The probes taken towards the edge. */
            int probes;
        };

        /**
         * @brief The outward search behind search() (`PHASE5_API.md`).
         * @param fn The function.
         * @param x0 The start.
         * @return The bracket, or error::non_finite or error::no_sign_change.
         */
        template<real T, typename F>
        constexpr std::expected<Bracket<T>, error> outward(F &fn, T const &x0)
        {
            if (!math::isfinite(x0)) return std::unexpected(error::non_finite);
            T const f0 = detail::evaluate<1>(fn, x0)[0];
            if (!math::isfinite(f0)) return std::unexpected(error::non_finite);
            if (f0 == T(0)) return Bracket<T>{ x0, x0 };
            // |x0| / 50 rather than |x0| (1/50), whose rounding differs; at least denorm_min, since |x0| / 50
            // underflows to 0 for |x0| below 25 denorm_min.
            T h = x0 == T(0) ? T(T(1) / T(50)) : T(math::abs(x0) / T(50));
            if (h == T(0)) h = T(std::numeric_limits<T>::denorm_min());
            // A type whose denorm_min() is 0, where the standard says min() (independent reviews of Phase 5).
            if (h == T(0)) h = math::abs(x0);
            T const growth = T(T(3) / T(2));
            T const half_max = T(detail::largest<T>() / T(2));
            // Probes towards a domain edge: the bits of T, enough to approach an edge at 0 by as many halvings.
            int const probe_limit = detail::iteration_cap<T>() / 4;
            std::array<SearchSide<T>, 2> sides{ SearchSide<T>{ -1, true, x0, std::nullopt, 0 },
                SearchSide<T>{ 1, true, x0, std::nullopt, 0 } };
            for (int step = 0; step < detail::iteration_cap<T>(); ++step) {
                for (SearchSide<T> &side : sides) {
                    if (!side.alive) continue;
                    T x = x0;
                    if (side.edge) {
                        if (side.probes >= probe_limit) {
                            side.alive = false;
                            continue;
                        }
                        x = detail::midpoint(side.last, *side.edge);
                        bool const between =
                            side.last < *side.edge ? side.last < x && x < *side.edge : *side.edge < x && x < side.last;
                        if (!between) {
                            side.alive = false;
                            continue;
                        }
                        ++side.probes;
                    } else {
                        T const s = side.direction < 0 ? T(-h) : h;
                        if (!detail::sum_fits(x0, s)) {
                            side.alive = false;
                            continue;
                        }
                        x = T(x0 + s);
                    }
                    T const fx = detail::evaluate<1>(fn, x)[0];
                    if (!math::isfinite(fx)) {
                        side.edge = x;
                        continue;
                    }
                    if (fx == T(0)) return Bracket<T>{ x, x };
                    if (!detail::same_sign(fx, f0)) {
                        if (side.last < x) return Bracket<T>{ side.last, x };
                        return Bracket<T>{ x, side.last };
                    }
                    side.last = x;
                }
                if (!sides[0].alive && !sides[1].alive) return std::unexpected(error::no_sign_change);
                // h grows by 3/2 where that fits: h (3/2) below max / 2, and above it h + h / 2, checked (the same
                // value in binary, where h / 2 is exact). A guard h <= max / (3/2) let 1.5 h overflow where the quotient
                // rounds up, as in double and std::float16_t (independent reviews of Phase 5). Else no larger step.
                if (h <= half_max)
                    h = T(h * growth);
                else if (detail::sum_fits(h, T(h / T(2))))
                    h = T(h + T(h / T(2)));
                else
                    return std::unexpected(error::no_sign_change);
            }
            return std::unexpected(error::no_sign_change);
        }

    } // namespace detail
    /** @endcond */

    /**
     * @brief A bracket around a sign change of @p fn, searched outward from @p x0.
     *
     * fn is evaluated at x0, then at x0 - h and x0 + h for h = |x0| / 50 (1/50 at x0 = 0) growing by 3/2 per step,
     * until fn's sign differs from its sign at x0. The result is the narrowest bracket that shows the sign change: the
     * last point on that side where fn had x0's sign, and the new point, in order; {x, x} where fn is exactly 0 at x,
     * which solve() returns at once as a root. x0 sets the scale of the search, so a quantity searches at its own
     * scale; from a tiny x0 a root at scale 1 can lie beyond the cap. At most 4 times the bits of T steps (212 in
     * double, so at most 425 evaluations), and no step that would pass max: from x0 = max, whose largest step is
     * 0.77 max, the search reaches down to 0.23 max.
     *
     * Where fn is not finite at a point (outside its domain), that side bisects towards its last finite point, at most
     * as many times as T has bits, to find a root near the domain's edge. Where fn is finite outside its domain,
     * return NaN there, so that the search does not cross the edge: a pole at the edge, as van der Waals' at its
     * co-volume, would otherwise be bracketed.
     *
     * What it cannot find: a root of even multiplicity, two roots within one step (x^2 - 2 from 100), a root beyond a
     * region where fn is not finite on the same side, and a sign change beyond the cap. A pole brackets like a root (tan
     * x from 1.4), and solve() then returns the pole: check |fn(value)|. The bracket need not hold the root nearest x0,
     * and where several roots lie in it, solve() returns one of them.
     *
     * @code
     * namespace rt = nxx::roots;
     * auto const f = [](double t) { return std::exp(20.0 - 5000.0 / t) - 1.0; };
     * auto const r = rt::search(f, 300.0) // [231.7, 254.4], 14 evaluations
     *                    .and_then([&](rt::Bracket<double> b) { return rt::solve<rt::Brent>(f, b); }); // 250
     * @endcode
     *
     * @tparam T The real type, deduced from @p x0 (search(fn, 0) deduces int: write 0.0).
     * @tparam F A callable taking T const & and returning something convertible to T, or a std::pair, std::tuple or
     *           std::array of exactly T, whose first element is used.
     * @param fn The function.
     * @param x0 The start.
     * @return The bracket, its ends in order; error::non_finite if x0 or fn(x0) is not finite; error::no_sign_change
     *         if no sign change was found.
     */
    template<real T, typename F>
    requires detail::gives<F, T, 1>
    [[nodiscard]] constexpr std::expected<Bracket<T>, error> search(F &&fn, T const &x0)
    { return detail::outward(fn, x0); }

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
