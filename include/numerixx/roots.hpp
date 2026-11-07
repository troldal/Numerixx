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
 * A method computes and solve() decides: the method reports two points, fn's values there and its point after each
 * step (a Step), and solve() checks the report, applies the stop test and forms the bound. So a method of one's own,
 * a type that is a bracketing_method or an open_method, gets the same checks and the same guarantee as the
 * library's: `solve<M>(fn, start)`, or `solve(m, fn, start)` for a method with parameters. `observer(o)` before fn
 * reports each step to o.
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
         * A method found no root. From a guess: no sign change between its iterates within its iteration cap, or no
         * usable step (f' = 0, not finite or a failed std::expected, Secant's two values equal, a step that would
         * overflow). With a numerical derivative also where its stencil leaves fn's domain at an iterate, or where |f'|
         * is as small as the derivative's own error (a multiple root): see Newton with separate derivatives. On a
         * bracket: a method that could place no point strictly inside it, which the library's methods do not reach.
         */
        not_converged,

        /**
         * A method of one's own broke the contract of bracketing_method: a step with its point outside its points,
         * outside the caller's bracket, without a sign change, or not nested in the previous step. The
         * library's methods never return it.
         */
        invalid_step
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
     * log2(width / denorm_min) halvings: about 1100 in double, but far more in a Boost.Multiprecision type, whose
     * exponent range is much wider (about 2e9 in cpp_bin_float_50, 2e8 in cpp_dec_float_50). Set abs_tol to the
     * distance from 0 that is close enough there.
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
     * @brief The caller's tolerances after validation, as solve() passes them to a method: tol(x) = absolute +
     * relative |x|, so that a method can scale its own steps to the request.
     *
     * solve() accepts a bracket at most 2 tol(x) wide around the method's point x. A method uses at(x) for its minimum
     * step and its placement (Brent steps at least tol(b); Toms748 keeps its points 0.7 tol inside the ends), so that
     * it never proposes points the request cannot tell apart. solve() builds it from Criteria: absolute in
     * [denorm_min, max / 4], relative in [2 eps, 1].
     *
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
         * @return absolute + relative |x|, or max where relative |x| exceeds max / 2, before the sum could
         *         overflow (independent review of Phase 7, finding 3).
         */
        constexpr T at(T const &x) const
        {
            // scaled is at most |x|, since relative <= 1. Above max / 2 it saturates: max - absolute would round,
            // and could admit a sum that overflows; below it the sum is at most 3/4 max, since absolute <= max / 4.
            // Parenthesised, so that a max macro, as defined by <windows.h>, cannot expand.
            T const max = T((std::numeric_limits<T>::max)());
            T const scaled = T(relative * math::abs(x));
            if (scaled > T(max / T(2))) return max;
            return T(absolute + scaled);
        }

        /** @brief Compares both members. */
        friend constexpr bool operator==(Tolerance const &, Tolerance const &) = default;
    };

    /**
     * @brief A bracket as solve() hands it to a bracketing method, so that the method starts from a checked interval:
     * lower < upper, both finite, fn's values there finite, nonzero and of opposite signs.
     *
     * solve() has already validated the caller's bracket, put it in order, and evaluated fn at its ends and, for a
     * bracket around 0, at 0 (zero first), so the bracket may be half the caller's. A method that needs fn's
     * derivatives (values_of<M> > 1) receives them too: at_lower[0] is fn's value, at_lower[1] f', at_lower[2] f''.
     * The derivatives are not checked.
     *
     * @tparam T The real type.
     * @tparam K The number of values at each end: values_of<M>.
     */
    template<real T, std::size_t K = 1>
    struct Start
    {
        /** @brief The lower end. */
        T lower;

        /** @brief fn's values at the lower end. */
        std::array<T, K> at_lower;

        /** @brief The upper end. */
        T upper;

        /** @brief fn's values at the upper end; the value of the other sign. */
        std::array<T, K> at_upper;

        /** @brief Compares every member. */
        friend constexpr bool operator==(Start const &, Start const &) = default;
    };

    /**
     * @brief What a stepped method reports after its start and after each step, so that solve() can check it and
     * decide: the two points it holds, fn's values there, and its point.
     *
     * On a bracket: lower < upper with fn's values nonzero and of opposite signs, inside the caller's bracket and
     * nested in the previous step's (neither end outside it, and not both ends the same); or point an end at which fn
     * is exactly 0. solve() returns point once the bracket is at most 2 tol(point) wide or point is within tol(point)
     * of both ends, with abs_error the smaller of the bound on the width and the larger distance from point to an
     * end, and point with abs_error tol(point) where fn is exactly 0 there. From a guess,
     * before a sign change: the last two iterates in order (one point after a start that evaluates only the guess;
     * Secant's start holds its first two points), point the latest; such a step
     * counts against the iteration cap. A method that cannot go on returns error::not_converged from advance.
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Step
    {
        /** @brief The lower point. */
        T lower;

        /** @brief fn's value at lower: its true value, not one the method has scaled. */
        T f_lower;

        /** @brief The upper point. */
        T upper;

        /** @brief fn's value at upper. */
        T f_upper;

        /** @brief The method's estimate: lower, upper, or a point between them (Bisection's midpoint). */
        T point;

        /** @brief Compares every member. */
        friend constexpr bool operator==(Step const &, Step const &) = default;
    };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief M::values where M declares it as a static integral constant of at least 1, else 1: the value of
         * values_of, in a function so that the documentation's C++ parser reads values_of's initializer.
         * @tparam M The method.
         * @return The number of values.
         */
        // Only a static integral constant counts: a method with a data member that happens to be called values (a
        // table of coefficients) needs one value (numerics review of Phase 7's draft 1, finding 7). At least 1: fn's
        // value is always read, so values = 0 would index an empty array (review of step 1, finding 4).
        template<typename M>
        consteval std::size_t declared_values()
        {
            if constexpr (requires {
                              requires std::integral<std::remove_cvref_t<decltype(M::values)>>;
                              typename std::integral_constant<std::size_t, M::values>;
                              requires M::values >= 1;
                          })
                return std::size_t{ M::values };
            else
                return std::size_t{ 1 };
        }

    } // namespace detail
    /** @endcond */

    /**
     * @brief How many values a method needs from fn at a point: M::values where M declares it as a static integral
     * constant of at least 1 (Newton 2, Halley 3), else 1, fn's value alone.
     *
     * solve() evaluates that many at the bracket's ends for the method's Start, and requires fn to give them; a method
     * with `values = 2` gets the form of solve() with the derivative as a separate callable.
     *
     * @tparam M The method.
     */
    template<typename M>
    inline constexpr std::size_t values_of = detail::declared_values<M>();

    /** @cond INTERNAL */
    namespace detail {

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
         * @brief A function returning one value: nxx::detail::callable (for a class type T, the result is T, not an
         * expression), and not a std::pair, std::tuple or std::array (a tuple-like result is never read as a value,
         * even where T converts from it).
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
         * @brief A separate derivative callable: one returning a value (value_valued, so T and not an expression for a
         * class type T), or a std::expected whose value converts to T, such as nxx::deriv::derivative_of's callable,
         * whose value is an nxx::deriv::Estimate<T>. The expected's error is never read. The expected holds its value,
         * so the value need not be T: it outlives the callable, unless it is itself an expression type, which only a
         * deliberately spelled std::expected of one can be.
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
            T const m = T(a + T(T(b - a) / T(2)));
            if ((a < m && m < b) || (b < m && m < a)) return m;
            // Ends of one sign with no point formed between them: adjacent (then none exists), or, in a type without
            // subnormals, closer than 2 min, where (b - a) / 2 flushes to 0 (numerics review of Phase 7's draft 1,
            // finding 1). Then both are tiny and their sum fits; the guard keeps the sum from overflowing near max,
            // which constant evaluation would reject.
            T const quarter = T(detail::largest<T>() / T(4));
            if (!(math::abs(a) <= quarter && math::abs(b) <= quarter)) return m;
            T const sum = T(T(a + b) / T(2));
            if ((a < sum && sum < b) || (b < sum && sum < a)) return sum;
            return m;
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
         * @brief A validated bracket for the method, or a root found during validation, so that solve() runs the
         * method only where validation left a sign change.
         * @tparam T The real type.
         * @tparam K The number of values at each end.
         */
        template<real T, std::size_t K>
        struct Prepared
        {
            /** @brief The root, if validation found one (fn exactly 0 at an end or at 0). */
            std::optional<Root<T>> root;

            /** @brief The bracket, in order, with fn's values at its ends; meaningful only without a root. */
            Start<T, K> start;
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
         * @return The start or the root found, or error::non_finite or error::no_sign_change.
         */
        template<std::size_t K, real T, typename F>
        constexpr std::expected<Prepared<T, K>, error>
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
            Start<T, K> const ends{ lower, at_lower, upper, at_upper };
            if (at_lower[0] == T(0)) return Prepared<T, K>{ Root<T>{ lower, tolerance.at(lower) }, ends };
            if (at_upper[0] == T(0)) return Prepared<T, K>{ Root<T>{ upper, tolerance.at(upper) }, ends };
            if (detail::same_sign(at_lower[0], at_upper[0])) return std::unexpected(error::no_sign_change);
            if (lower < T(0) && T(0) < upper) {
                T const zero = T(0);
                std::array<T, K> const at_zero = detail::evaluate<K>(fn, zero);
                if (at_zero[0] == T(0)) return Prepared<T, K>{ Root<T>{ zero, tolerance.at(zero) }, ends };
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
            return Prepared<T, K>{ std::nullopt, Start<T, K>{ lower, at_lower, upper, at_upper } };
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
         * @brief The step of a bracket [a, b] in order whose point is the end with the smaller |f| (b on a tie): the
         * report of Toms748, AndersonBjorck, RegulaFalsi, Illinois and Ridders.
         * @param a The lower end.
         * @param fa fn there.
         * @param b The upper end.
         * @param fb fn there.
         * @return The step.
         */
        template<real T>
        constexpr Step<T> bracket_step(T const &a, T const &fa, T const &b, T const &fb)
        { return Step<T>{ a, fa, b, fb, math::abs(fa) < math::abs(fb) ? a : b }; }


        // ---- Bisection --------------------------------------------------------------------------------------------

        /**
         * @brief Bisection's state: the bracket and fn's values at its ends.
         * @tparam T The real type.
         */
        template<real T>
        struct BisectionState
        {
            /** @brief The lower end. */
            T a;
            /** @brief fn at a. */
            T fa;
            /** @brief The upper end. */
            T b;
            /** @brief fn at b. */
            T fb;
            /** @brief The report: the bracket, its midpoint the point. */
            Step<T> step;
        };

        /**
         * @brief Bisection's report on [a, b]: the midpoint as its point.
         * @param a The lower end.
         * @param fa fn there.
         * @param b The upper end.
         * @param fb fn there.
         * @return The step.
         */
        template<real T>
        constexpr Step<T> bisection_step(T const &a, T const &fa, T const &b, T const &fb)
        {
            return Step<T>{ a, fa, b, fb, detail::midpoint(a, b) };
        }

        // ---- Brent's procedure zero, with Brent's, Newton's or Halley's step ---------------------------------------

        /** @brief The step inside Brent's safeguard: his interpolation, or Newton's or Halley's derivative step. */
        enum class Safeguard { interpolation, newton, halley };

        /**
         * @brief The derivative step of the safeguard G from the values at a point.
         * @tparam G Safeguard::newton or Safeguard::halley.
         * @param v The values at the point.
         * @return The step, or nothing where it cannot be formed.
         */
        template<Safeguard G, real T, std::size_t K>
        constexpr std::optional<T> derivative_step(std::array<T, K> const &v)
        {
            if constexpr (G == Safeguard::halley)
                return detail::halley_step(v);
            else
                return detail::newton_step(v);
        }

        /**
         * @brief The state of Brent's procedure zero (Brent 1973, ch. 4, section 6): [b, c] the bracket with |fn(b)| <=
         * |fn(c)|, a the previous b, d the last step and e the one before.
         * @tparam T The real type.
         * @tparam K The number of values at each point.
         */
        template<real T, std::size_t K>
        struct SafeguardedState
        {
            /** @brief The previous b. */
            T a;
            /** @brief fn's values at a. */
            std::array<T, K> va;
            /** @brief The best point. */
            T b;
            /** @brief fn's values at b. */
            std::array<T, K> vb;
            /** @brief The other end of the bracket. */
            T c;
            /** @brief fn's values at c. */
            std::array<T, K> vc;
            /** @brief The last step. */
            T d;
            /** @brief The step before it. */
            T e;
            /** @brief The report: the bracket [b, c] in order, b the point. */
            Step<T> step;
        };

        /**
         * @brief The top of Brent's loop, then the report: c reset to a where fn(b) and fn(c) have one sign, b and c
         * swapped where |fn(c)| < |fn(b)|.
         * @param s The state, its b just evaluated (fn nonzero there).
         */
        template<real T, std::size_t K>
        constexpr void normalize(SafeguardedState<T, K> &s)
        {
            T const big = T(detail::largest<T>() / T(16));
            if (detail::same_sign(s.vb[0], s.vc[0])) {
                s.c = s.a;
                s.vc = s.va;
                s.d = detail::bounded_difference(s.b, s.a, big);
                s.e = s.d;
            }
            if (math::abs(s.vc[0]) < math::abs(s.vb[0])) {
                s.a = s.b;
                s.b = s.c;
                s.c = s.a;
                s.va = s.vb;
                s.vb = s.vc;
                s.vc = s.va;
            }
            s.step = s.b < s.c ? Step<T>{ s.b, s.vb[0], s.c, s.vc[0], s.b }
                               : Step<T>{ s.c, s.vc[0], s.b, s.vb[0], s.b };
        }

        /**
         * @brief Brent's state on a validated bracket.
         * @param s The bracket.
         * @return The state, reported.
         */
        template<real T, std::size_t K>
        constexpr SafeguardedState<T, K> safeguarded_start(Start<T, K> const &s)
        {
            T const big = T(detail::largest<T>() / T(16));
            T const d = detail::bounded_difference(s.upper, s.lower, big);
            SafeguardedState<T, K> state{ s.lower, s.at_lower, s.upper, s.at_upper, s.lower, s.at_lower, d, d,
                Step<T>{ s.lower, s.at_lower[0], s.upper, s.at_upper[0], s.upper } };
            detail::normalize(state);
            return state;
        }

        /**
         * @brief One iteration of Brent's procedure zero from the step's choice to the next point's evaluation, with
         * the step chosen by G (detail::brent of Phase 1 to 6, its stop test moved to solve()).
         *
         * Implemented from the book, with three changes: the tolerance is the caller's tol(b) rather than
         * 2 macheps |b| + t; the stop test (in solve()) is distance_bound(c, b) <= 2 tol(b) rather than |c - b| / 2 <=
         * tol(b); and a step is a bisection step whenever |c - b|, |b - a| or |e| exceeds max / 16, so that no product
         * of the interpolation overflows. A derivative step is accepted by Brent's test written for a step s: towards
         * c, |s| < 3 |m| / 2 - tol / 2, and |s| < |e| / 2; one that cannot be formed is a bisection step.
         *
         * @tparam G The step.
         * @param fn The function.
         * @param s The state, normalized and not yet meeting the stop test.
         * @param tolerance The tolerances.
         * @return Nothing, or error::non_finite if fn is not finite at the new point.
         */
        template<Safeguard G, std::size_t K, real T, typename F>
        constexpr std::expected<void, error>
            safeguarded_advance(F &fn, SafeguardedState<T, K> &s, Tolerance<T> const &tolerance)
        {
            T const big = T(detail::largest<T>() / T(16));
            T const fa = s.va[0];
            T const fb = s.vb[0];
            T const fc = s.vc[0];
            T const tol = tolerance.at(s.b);
            T const width = detail::distance_bound(s.c, s.b);
            bool const wide = width > big || detail::distance_bound(s.b, s.a) > big || math::abs(s.e) > big;
            // Half the bracket, towards c: exact Brent below big, the overflow-free midpoint above it.
            T const m = wide ? T(detail::midpoint(s.b, s.c) - s.b) : T(T(s.c - s.b) / T(2));
            if (wide || math::abs(s.e) < tol || math::abs(fa) <= math::abs(fb)) {
                s.d = m;
                s.e = m;
            } else if constexpr (G == Safeguard::interpolation) {
                T const ratio = T(fb / fa);
                T p = T(0);
                T q = T(0);
                if (s.a == s.c) {
                    // Linear interpolation (the secant).
                    p = T(T(2) * m * ratio);
                    q = T(T(1) - ratio);
                } else {
                    // Inverse quadratic interpolation.
                    T const qa = T(fa / fc);
                    T const r = T(fb / fc);
                    p = T(ratio * T(T(T(2) * m * qa * T(qa - r)) - T(T(s.b - s.a) * T(r - T(1)))));
                    q = T(T(qa - T(1)) * T(r - T(1)) * T(ratio - T(1)));
                }
                if (p > T(0))
                    q = T(-q);
                else
                    p = T(-p);
                T const before_last = s.e; // Brent's s := e
                s.e = s.d;
                // Accept the interpolated point only if it falls well inside the bracket and the steps shrink. A NaN
                // or q = 0 fails both tests, and falls back on bisection.
                if (T(T(2) * p) < T(T(T(3) * m * q) - math::abs(T(tol * q)))
                    && p < math::abs(T(before_last * q / T(2)))) {
                    s.d = T(p / q);
                } else {
                    s.d = m;
                    s.e = m;
                }
            } else {
                std::optional<T> const step = detail::derivative_step<G>(s.vb);
                T const before_last = s.e;
                s.e = s.d;
                // Brent's acceptance for a step: towards c, well inside the bracket (so that b + step cannot round
                // onto c), and shorter than half the step before last.
                T const inside = T(T(T(3) * math::abs(m)) / T(2) - T(tol / T(2)));
                if (step && (*step > T(0)) == (m > T(0)) && math::abs(*step) < inside
                    && math::abs(*step) < T(math::abs(before_last) / T(2))) {
                    s.d = *step;
                } else {
                    s.d = m;
                    s.e = m;
                }
            }
            T next = s.b;
            if (math::abs(s.d) > tol)
                next = T(s.b + s.d);
            else if (s.b < s.c) // towards c; m > 0 says the same wherever m is not 0, and m flushes to 0 in a type
                                // without subnormals where c - b is below 2 min (numerics review of draft 1, finding 1)
                next = T(s.b + tol);
            else
                next = T(s.b - tol);
            // A point not strictly inside the bracket: no step can make progress. Not reached by the library's
            // arithmetic where the stop test failed (PHASE7_API.md B); not_converged rather than a root outside the
            // bound.
            if (!((s.b < next && next < s.c) || (s.c < next && next < s.b)))
                return std::unexpected(error::not_converged);
            s.a = s.b;
            s.va = s.vb;
            s.b = next;
            s.vb = detail::evaluate<K>(fn, s.b);
            if (!math::isfinite(s.vb[0])) return std::unexpected(error::non_finite);
            if (s.vb[0] == T(0)) {
                s.step = s.b < s.c ? Step<T>{ s.b, s.vb[0], s.c, s.vc[0], s.b }
                                   : Step<T>{ s.c, s.vc[0], s.b, s.vb[0], s.b };
                return {};
            }
            detail::normalize(s);
            return {};
        }

        // ---- The bracket subroutine of Toms748, AndersonBjorck and Ridders -----------------------------------------

        /**
         * @brief One step of the paper's bracket subroutine, without its stop test (solve()'s): the new state, the
         * point evaluated, and whether fn is exactly 0 there or no point was strictly inside.
         * @tparam T The real type.
         */
        template<real T>
        struct Narrowed
        {
            /** @brief The state after the step. */
            Enclosure<T> next;
            /** @brief The point evaluated. */
            T c;
            /** @brief fn at c. */
            T fc;
            /** @brief fn is exactly 0 at c. */
            bool zero;
            /** @brief No point strictly inside (not reached: the stop test meets such a bracket first). */
            bool stuck;
        };

        /**
         * @brief One step of Toms748, AndersonBjorck, the false-position methods and Ridders, the paper's bracket
         * subroutine: the point (the midpoint while the bracket is wide, else the proposed point placed as place()
         * says), e := d where asked, evaluate fn, keep the half with the sign change and record the dropped end as d.
         * @param fn The function.
         * @param s The state.
         * @param proposed The method's point, computed only where the bracket is not wide; or nothing.
         * @param shift Whether e takes d's value before d changes.
         * @param tolerance The tolerances.
         * @return The step, or error::non_finite if fn is not finite at the point.
         */
        template<real T, typename F>
        constexpr std::expected<Narrowed<T>, error> narrow(F &fn,
            Enclosure<T> const &s,
            std::optional<T> const &proposed,
            bool shift,
            Tolerance<T> const &tolerance)
        {
            T const c = detail::wide(s) ? detail::midpoint(s.a, s.b) : detail::place(s, proposed, tolerance);
            // Not reached: a midpoint not strictly inside means adjacent ends, which met the stop test before
            // (midpoint() falls back on (a + b) / 2 where a + (b - a) / 2 flushes, numerics review of draft 1, finding
            // 1).
            assert(s.a < c && c < s.b);
            if (!(s.a < c && c < s.b)) return Narrowed<T>{ s, c, s.fa, false, true };
            Enclosure<T> next = s;
            if (shift) {
                next.e = s.d;
                next.fe = s.fd;
            }
            T const fc = detail::evaluate<1>(fn, c)[0];
            if (!math::isfinite(fc)) return std::unexpected(error::non_finite);
            if (fc == T(0)) return Narrowed<T>{ next, c, fc, true, false };
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
            return Narrowed<T>{ next, c, fc, false, false };
        }

        /**
         * @brief The report after a step of narrow() that found a point inside: the new bracket, or [a, c] with c the
         * point where fn is 0 there.
         * @param before The state before the step.
         * @param r The step.
         * @return The step.
         */
        template<real T>
        constexpr Step<T> narrowed_step(Enclosure<T> const &before, Narrowed<T> const &r)
        {
            if (r.zero) return Step<T>{ before.a, before.fa, r.c, r.fc, r.c };
            return detail::bracket_step(r.next.a, r.next.fa, r.next.b, r.next.fb);
        }

        /**
         * @brief Toms748's state: the bracket subroutine's state, where Algorithm 4.2 stands (stage 0 the first
         * secant step; 1 and 2 the interpolation steps; 3 the double-length secant step; 4 the bisection step), and
         * the bracket at the start of the iteration.
         * @tparam T The real type.
         */
        template<real T>
        struct Toms748State
        {
            /** @brief The bracket and the dropped ends. */
            Enclosure<T> s;
            /** @brief The lower end at the start of the iteration. */
            T a0;
            /** @brief The upper end at the start of the iteration. */
            T b0;
            /** @brief The next step of the algorithm. */
            int stage;
            /** @brief The first iteration, which uses Newton-Quadratic without e. */
            bool first;
            /** @brief The report. */
            Step<T> step;
        };

        /**
         * @brief One step of Alefeld, Potra and Shi's Algorithm 4.2 (1995, pp. 334-335): one evaluation, through
         * narrow(). Stage 0: a secant step (4.2.1-4.2.2); per iteration, two interpolation steps (inverse cubic, or
         * Newton-Quadratic with 2 and then 3 Newton steps; 4.2.3-4.2.6), a double-length secant step (4.2.7-4.2.10)
         * and a bisection step where the iteration has not brought the bracket below half its width at its start
         * (4.2.11).
         * @param fn The function.
         * @param st The state.
         * @param tolerance The tolerances.
         * @return Nothing, error::non_finite where fn is not finite at a point, or error::not_converged where the
         *         method cannot go on.
         */
        template<real T, typename F>
        constexpr std::expected<void, error> toms748_advance(F &fn, Toms748State<T> &st, Tolerance<T> const &tolerance)
        {
            Enclosure<T> const before = st.s;
            bool const w = detail::wide(before);
            std::optional<T> proposed;
            bool shift = false;
            switch (st.stage) {
            case 0:
                if (!w) proposed = detail::secant_point(before.a, before.fa, before.b, before.fb);
                break;
            case 1: // e := d before the bracket
                if (!w) proposed = detail::interpolate(before, 2, st.first);
                shift = true;
                break;
            case 2:
                if (!w) proposed = detail::interpolate(before, 3, false);
                break;
            case 3: // e := d-bar (4.1.9)
                if (!w) proposed = detail::double_secant(before);
                shift = true;
                break;
            default: // e := d-hat
                if (!w) proposed = detail::midpoint(before.a, before.b);
                shift = true;
                break;
            }
            std::expected<Narrowed<T>, error> const r = detail::narrow(fn, before, proposed, shift, tolerance);
            if (!r) return std::unexpected(r.error());
            if (r->stuck) return std::unexpected(error::not_converged);
            st.step = detail::narrowed_step(before, *r);
            if (r->zero) return {};
            st.s = r->next;
            switch (st.stage) {
            case 0:
                st.stage = 1;
                st.a0 = st.s.a;
                st.b0 = st.s.b;
                break;
            case 1:
                st.stage = 2;
                break;
            case 2:
                st.stage = 3;
                break;
            case 3:
                if (detail::distance_bound(st.s.a, st.s.b) < T(detail::distance_bound(st.a0, st.b0) / T(2))) {
                    st.first = false;
                    st.a0 = st.s.a;
                    st.b0 = st.s.b;
                    st.stage = 1;
                } else {
                    st.stage = 4;
                }
                break;
            default:
                st.first = false;
                st.a0 = st.s.a;
                st.b0 = st.s.b;
                st.stage = 1;
                break;
            }
            return {};
        }

        // ---- False position --------------------------------------------------------------------------------------

        /** @brief The factor by which false position scales the retained end's value when it is retained again. */
        enum class Scaling { none, illinois, anderson_bjorck };

        /**
         * @brief False position's state: b the latest point, a the retained end, whose value is replaced by g_a.
         * @tparam T The real type.
         */
        template<real T>
        struct FalsePositionState
        {
            /** @brief The retained end. */
            T retained;
            /** @brief fn there. */
            T f_retained;
            /** @brief Its scaled value, g_a. */
            T scaled;
            /** @brief The latest point. */
            T latest;
            /** @brief fn there. */
            T f_latest;
            /** @brief The bracket's width at the last halving. */
            T reference;
            /** @brief Steps since then. */
            int unhalved;
            /** @brief The report. */
            Step<T> step;
        };

        /**
         * @brief One step of false position, the retained end's value scaled by S's factor, guarded by a bisection
         * after three consecutive steps that have not brought the bracket below half its width at the last halving (a
         * guard bisection counts as one); the point taken through narrow().
         *
         * The point is the secant's zero through (a, g_a) and (b, f_b). Where f_c has f_b's sign, a is retained and g_a
         * scaled: AndersonBjorck's m = 1 - f_c / f_b, or 1/2 where m is not positive (formed only where |f_c| < |f_b|;
         * Anderson and Bjorck 1973); Illinois' 1/2; RegulaFalsi's 1. Otherwise b becomes the retained end with its true
         * value. c becomes b.
         *
         * @tparam S The scaling.
         * @param fn The function.
         * @param st The state.
         * @param tolerance The tolerances.
         * @return Nothing, error::non_finite where fn is not finite at a point, or error::not_converged where the
         *         method cannot go on.
         */
        template<Scaling S, real T, typename F>
        constexpr std::expected<void, error>
            false_position_advance(F &fn, FalsePositionState<T> &st, Tolerance<T> const &tolerance)
        {
            bool const ordered = st.retained < st.latest;
            Enclosure<T> const s{ ordered ? st.retained : st.latest, ordered ? st.f_retained : st.f_latest,
                ordered ? st.latest : st.retained, ordered ? st.f_latest : st.f_retained, st.retained, st.f_retained,
                st.retained, st.f_retained };
            bool const bisect = st.unhalved >= 3;
            std::optional<T> proposed;
            if (!detail::wide(s)) {
                if (bisect) {
                    proposed = detail::midpoint(s.a, s.b);
                } else if (std::optional<T> const spread = detail::difference(st.f_latest, st.scaled)) {
                    // f_latest and scaled have opposite signs, so 0 < r < 1.
                    if (std::optional<T> const r = detail::quotient(st.f_latest, *spread))
                        proposed = T(st.latest - T(*r * T(st.latest - st.retained)));
                }
            }
            std::expected<Narrowed<T>, error> const r = detail::narrow(fn, s, proposed, false, tolerance);
            if (!r) return std::unexpected(r.error());
            if (r->stuck) return std::unexpected(error::not_converged);
            st.step = detail::narrowed_step(s, *r);
            if (r->zero) return {};
            if (!detail::same_sign(r->fc, st.f_latest)) {
                st.retained = st.latest;
                st.f_retained = st.f_latest;
                st.scaled = st.f_latest;
            } else if constexpr (S == Scaling::anderson_bjorck) {
                T const m = math::abs(r->fc) < math::abs(st.f_latest) ? T(T(1) - T(r->fc / st.f_latest))
                                                                        : T(T(1) / T(2));
                st.scaled = T(st.scaled * m);
            } else if constexpr (S == Scaling::illinois) {
                st.scaled = T(st.scaled * T(T(1) / T(2)));
            } // RegulaFalsi: the factor 1
            st.latest = r->c;
            st.f_latest = r->fc;
            T const width = detail::distance_bound(st.retained < st.latest ? st.retained : st.latest,
                st.retained < st.latest ? st.latest : st.retained);
            if (bisect || width <= T(st.reference / T(2))) {
                st.reference = width;
                st.unhalved = 0;
            } else {
                ++st.unhalved;
            }
            return {};
        }

        // ---- Ridders ---------------------------------------------------------------------------------------------

        /**
         * @brief Ridders' state: the bracket subroutine's state, and the midpoint just evaluated where the fit's point
         * comes next.
         * @tparam T The real type.
         */
        template<real T>
        struct RiddersState
        {
            /** @brief The bracket. */
            Enclosure<T> s;
            /** @brief The midpoint evaluated by the last step (stage 1). */
            T m;
            /** @brief fn there. */
            T fm;
            /** @brief 0: the midpoint next; 1: the fit's point next. */
            int stage;
            /** @brief The report. */
            Step<T> step;
        };

        /**
         * @brief One step of Ridders' method: alternately narrow() at the midpoint, and at Ridders' point from the
         * state that step left. Each pair halves the bracket.
         * @param fn The function.
         * @param st The state.
         * @param tolerance The tolerances.
         * @return Nothing, error::non_finite where fn is not finite at a point, or error::not_converged where the
         *         method cannot go on.
         */
        template<real T, typename F>
        constexpr std::expected<void, error> ridders_advance(F &fn, RiddersState<T> &st, Tolerance<T> const &tolerance)
        {
            Enclosure<T> const before = st.s;
            std::optional<T> proposed;
            if (st.stage == 1 && !detail::wide(before)) proposed = detail::ridders_point(before, st.m, st.fm);
            std::expected<Narrowed<T>, error> const r = detail::narrow(fn, before, proposed, false, tolerance);
            if (!r) return std::unexpected(r.error());
            if (r->stuck) return std::unexpected(error::not_converged);
            st.step = detail::narrowed_step(before, *r);
            if (r->zero) return {};
            st.s = r->next;
            if (st.stage == 0) {
                st.m = r->c;
                st.fm = r->fc;
                st.stage = 1;
            } else {
                st.stage = 0;
            }
            return {};
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
         * @brief The step of two points x and y, y the latest and the point, in order: an open method's report before
         * its sign change, or where fn is exactly 0 at y.
         * @param x The earlier point.
         * @param fx fn there.
         * @param y The latest point.
         * @param fy fn there.
         * @return The step.
         */
        template<real T>
        constexpr Step<T> pair_step(T const &x, T const &fx, T const &y, T const &fy)
        {
            if (x < y) return Step<T>{ x, fx, y, fy, y };
            return Step<T>{ y, fy, x, fx, y };
        }

        // ---- The open methods ------------------------------------------------------------------------------------

        /**
         * @brief Brent's procedure on two iterates of an open method, whose values have opposite signs: the
         * hand-over. No zero first: they are adjacent iterates, not the caller's bracket.
         * @param x One iterate.
         * @param vx fn's values there.
         * @param y The other.
         * @param vy fn's values there.
         * @return The state, reported.
         */
        template<real T, std::size_t K>
        constexpr SafeguardedState<T, K>
            hand_over(T const &x, std::array<T, K> const &vx, T const &y, std::array<T, K> const &vy)
        {
            if (x < y) return detail::safeguarded_start(Start<T, K>{ x, vx, y, vy });
            return detail::safeguarded_start(Start<T, K>{ y, vy, x, vx });
        }

        /**
         * @brief Newton's or Halley's state from a guess: the latest iterate and, once two iterates have shown a sign
         * change, Brent's procedure with the same step on them.
         * @tparam T The real type.
         * @tparam K The number of values: 2 or 3.
         */
        template<real T, std::size_t K>
        struct OpenDerivativeState
        {
            /** @brief The latest iterate. */
            T x;
            /** @brief fn's values there. */
            std::array<T, K> vx;
            /** @brief The bracket loop after the hand-over. */
            std::optional<SafeguardedState<T, K>> finish;
            /** @brief The report. */
            Step<T> step;
        };

        /**
         * @brief One step of Newton or Halley from a guess: the method's step, at least tol(x) long, its direction from
         * the signs of f and f'; at the first sign change, the hand-over; then the bracket loop.
         * @tparam G Safeguard::newton or Safeguard::halley.
         * @param fn The function.
         * @param st The state.
         * @param tolerance The tolerances.
         * @return Nothing, error::non_finite where fn is not finite at a point, or error::not_converged where the
         *         method cannot go on.
         */
        template<Safeguard G, real T, std::size_t K, typename F>
        constexpr std::expected<void, error>
            open_derivative_advance(F &fn, OpenDerivativeState<T, K> &st, Tolerance<T> const &tolerance)
        {
            if (st.finish) {
                std::expected<void, error> const r = detail::safeguarded_advance<G>(fn, *st.finish, tolerance);
                st.step = st.finish->step;
                return r;
            }
            std::optional<T> const step = detail::derivative_step<G>(st.vx);
            if (!step) return std::unexpected(error::not_converged);
            T const tol = tolerance.at(st.x);
            bool const upward = !detail::same_sign(st.vx[0], st.vx[1]); // Newton's direction, -f / f'
            T const move = math::abs(*step) > tol ? *step : (upward ? tol : T(-tol));
            if (!detail::sum_fits(st.x, move)) return std::unexpected(error::not_converged);
            T const next = T(st.x + move);
            std::array<T, K> const vn = detail::evaluate<K>(fn, next);
            if (!math::isfinite(vn[0])) return std::unexpected(error::non_finite);
            if (vn[0] != T(0) && !detail::same_sign(vn[0], st.vx[0])) {
                st.finish = detail::hand_over(st.x, st.vx, next, vn);
                st.step = st.finish->step;
                return {};
            }
            st.step = detail::pair_step(st.x, st.vx[0], next, vn[0]);
            st.x = next;
            st.vx = vn;
            return {};
        }

        /**
         * @brief The secant method's state from a guess: the latest two iterates and the last move, then Brent.
         * @tparam T The real type.
         */
        template<real T>
        struct SecantState
        {
            /** @brief The previous iterate. */
            T xp;
            /** @brief fn there. */
            T fp;
            /** @brief The latest iterate. */
            T x;
            /** @brief fn there. */
            T fx;
            /** @brief x - xp. */
            T dx;
            /** @brief Brent's procedure after the hand-over. */
            std::optional<SafeguardedState<T, 1>> finish;
            /** @brief The report. */
            Step<T> step;
        };

        /**
         * @brief The secant method's report on a new iterate x after xp: the hand-over at a sign change, else the
         * two iterates (the zero where fn is 0 at x).
         * @param st The state, x just evaluated and finite.
         */
        template<real T>
        constexpr void secant_report(SecantState<T> &st)
        {
            if (st.fx != T(0) && !detail::same_sign(st.fx, st.fp)) {
                st.finish = detail::hand_over(st.xp, std::array<T, 1>{ st.fp }, st.x, std::array<T, 1>{ st.fx });
                st.step = st.finish->step;
                return;
            }
            st.step = detail::pair_step(st.xp, st.fp, st.x, st.fx);
        }

        /**
         * @brief One step of the secant method from a guess: the step (x - xp) (f / (fp - f)) through the latest two
         * iterates, at least tol(x) long; at the first sign change Brent on them.
         * @param fn The function.
         * @param st The state.
         * @param tolerance The tolerances.
         * @return Nothing, error::non_finite where fn is not finite at a point, or error::not_converged where the
         *         method cannot go on.
         */
        template<real T, typename F>
        constexpr std::expected<void, error> secant_advance(F &fn, SecantState<T> &st, Tolerance<T> const &tolerance)
        {
            if (st.finish) {
                std::expected<void, error> const r =
                    detail::safeguarded_advance<Safeguard::interpolation>(fn, *st.finish, tolerance);
                st.step = st.finish->step;
                return r;
            }
            std::optional<T> const move = detail::secant_move(st.x, st.fx, st.dx, T(st.fp - st.fx), tolerance);
            if (!move) return std::unexpected(error::not_converged);
            st.xp = st.x;
            st.fp = st.fx;
            st.x = T(st.x + *move);
            st.dx = *move;
            st.fx = detail::evaluate<1>(fn, st.x)[0];
            if (!math::isfinite(st.fx)) return std::unexpected(error::non_finite);
            detail::secant_report(st);
            return {};
        }

        /**
         * @brief Steffensen's state from a guess: the latest iterate, then Brent.
         * @tparam T The real type.
         */
        template<real T>
        struct SteffensenState
        {
            /** @brief The latest iterate. */
            T x;
            /** @brief fn there. */
            T fx;
            /** @brief Brent's procedure after the hand-over. */
            std::optional<SafeguardedState<T, 1>> finish;
            /** @brief The report. */
            Step<T> step;
        };

        /**
         * @brief One iteration of Steffensen's method from a guess: z = x + f(x), at least tol(x) away from x in f(x)'s
         * direction; the secant's step through (x, f(x)) and (z, f(z)); Brent at the first sign change, between x and z
         * or x and the next iterate. Up to two evaluations.
         * @param fn The function.
         * @param st The state.
         * @param tolerance The tolerances.
         * @return Nothing, error::non_finite where fn is not finite at a point, or error::not_converged where the
         *         method cannot go on.
         */
        template<real T, typename F>
        constexpr std::expected<void, error>
            steffensen_advance(F &fn, SteffensenState<T> &st, Tolerance<T> const &tolerance)
        {
            if (st.finish) {
                std::expected<void, error> const r =
                    detail::safeguarded_advance<Safeguard::interpolation>(fn, *st.finish, tolerance);
                st.step = st.finish->step;
                return r;
            }
            T const x = st.x;
            T const fx = st.fx;
            T const tol = tolerance.at(x);
            T const offset = math::abs(fx) > tol ? fx : (fx > T(0) ? tol : T(-tol));
            if (!detail::sum_fits(x, offset)) return std::unexpected(error::not_converged);
            T const z = T(x + offset);
            T const fz = detail::evaluate<1>(fn, z)[0];
            if (!math::isfinite(fz)) return std::unexpected(error::non_finite);
            if (fz == T(0)) {
                st.step = detail::pair_step(x, fx, z, fz);
                return {};
            }
            if (!detail::same_sign(fz, fx)) {
                st.finish = detail::hand_over(x, std::array<T, 1>{ fx }, z, std::array<T, 1>{ fz });
                st.step = st.finish->step;
                return {};
            }
            T const width = math::abs(offset) <= math::abs(x) ? T(z - x) : offset;
            std::optional<T> const move = detail::secant_move(x, fx, T(-width), T(fz - fx), tolerance);
            if (!move) return std::unexpected(error::not_converged);
            T const next = T(x + *move);
            T const fn_next = detail::evaluate<1>(fn, next)[0];
            if (!math::isfinite(fn_next)) return std::unexpected(error::non_finite);
            if (fn_next != T(0) && !detail::same_sign(fn_next, fx)) {
                st.finish = detail::hand_over(x, std::array<T, 1>{ fx }, next, std::array<T, 1>{ fn_next });
                st.step = st.finish->step;
                return {};
            }
            st.step = detail::pair_step(x, fx, next, fn_next);
            st.x = next;
            st.fx = fn_next;
            return {};
        }

    } // namespace detail
    /** @endcond */

    /**
     * @brief Bisection: halves the bracket at every step, one bit per evaluation of fn, guaranteed.
     *
     * About log2(width / (2 tol)) evaluations, 52 for [1, 2] in double, whatever fn is. Its point is the midpoint of
     * the bracket.
     */
    struct Bisection
    {
        /**
         * @brief The state on a validated bracket, reporting it with its midpoint.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param s The bracket.
         * @return The state.
         */
        template<real T, typename F>
        constexpr std::expected<detail::BisectionState<T>, error>
            start(F &, Start<T> const &s, Tolerance<T> const &) const
        {
            T const fa = s.at_lower[0];
            T const fb = s.at_upper[0];
            return detail::BisectionState<T>{ s.lower, fa, s.upper, fb,
                detail::bisection_step(s.lower, fa, s.upper, fb) };
        }

        /**
         * @brief One halving: fn at the midpoint, and the half with the sign change kept.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param s The state.
         * @return Nothing, or error::non_finite where fn is not finite at the midpoint.
         */
        template<real T, typename F>
        constexpr std::expected<void, error> advance(F &fn, detail::BisectionState<T> &s, Tolerance<T> const &) const
        {
            T const m = s.step.point;
            // Not reached where the stop test failed: the midpoint is then strictly inside (PHASE7_API.md B).
            if (!(s.a < m && m < s.b)) return std::unexpected(error::not_converged);
            T const fm = detail::evaluate<1>(fn, m)[0];
            if (!math::isfinite(fm)) return std::unexpected(error::non_finite);
            if (fm == T(0)) {
                s.step = Step<T>{ s.a, s.fa, m, fm, m };
                return {};
            }
            if (detail::same_sign(fm, s.fa)) {
                s.a = m;
                s.fa = fm;
            } else {
                s.b = m;
                s.fb = fm;
            }
            s.step = detail::bisection_step(s.a, s.fa, s.b, s.fb);
            return {};
        }
    };

    /**
     * @brief Brent's method: inverse quadratic and linear interpolation, safeguarded by bisection (Brent 1973).
     *
     * Superlinear on a smooth fn with a simple root, typically 8 to 12 evaluations in double, and never much slower
     * than Bisection: at most (k + 1)^2 - 2 evaluations, k being Bisection's count. Slower on a multiple root
     * ((x - 1)^3: 156 evaluations, against Bisection's 54). Its point is b, the end with the smaller |f|.
     */
    struct Brent
    {
        /**
         * @brief The state on a validated bracket: Brent's a, b, c, d and e.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param s The bracket.
         * @return The state.
         */
        template<real T, typename F>
        constexpr std::expected<detail::SafeguardedState<T, 1>, error>
            start(F &, Start<T> const &s, Tolerance<T> const &) const
        { return detail::safeguarded_start(s); }

        /**
         * @brief One iteration of Brent's procedure zero: an interpolated or bisection point, at least tol(b) from b,
         * evaluated.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param s The state.
         * @param t The tolerances.
         * @return Nothing, error::non_finite where fn is not finite at the point, or error::not_converged where no
         *         point lies strictly inside (not reached in the arithmetic solve() assumes).
         */
        template<real T, typename F>
        constexpr std::expected<void, error>
            advance(F &fn, detail::SafeguardedState<T, 1> &s, Tolerance<T> const &t) const
        { return detail::safeguarded_advance<detail::Safeguard::interpolation>(fn, s, t); }
    };

    /**
     * @brief Alefeld, Potra and Shi's Algorithm 748 (1995, Algorithm 4.2): inverse cubic interpolation with quadratic
     * and secant fallbacks, a double-length secant step, and a bisection whenever an iteration has not halved the
     * bracket. A bracket only; no derivative form.
     *
     * Costs about what Brent costs on smooth fn (the paper reports 4% fewer evaluations; on 27 of its test cases at
     * full precision Brent took 309, Toms748 325), with a better worst case: at most about 4 evaluations per halving
     * of the bracket, where Brent's bound is quadratic in the number of halvings. Like every bracketing method it
     * returns a point and a bound, not the final bracket (Boost.Math's toms748_solve returns the bracket). Its point is
     * the end with the smaller |f|; each step is one evaluation, so an iteration of the algorithm is up to four steps.
     */
    struct Toms748
    {
        /**
         * @brief The state on a validated bracket, before the algorithm's first secant step.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param s The bracket.
         * @return The state.
         */
        template<real T, typename F>
        constexpr std::expected<detail::Toms748State<T>, error>
            start(F &, Start<T> const &s, Tolerance<T> const &) const
        {
            T const fa = s.at_lower[0];
            T const fb = s.at_upper[0];
            // d and e are read only once steps have set them (the first iteration uses Newton-Quadratic, without e).
            return detail::Toms748State<T>{ detail::Enclosure<T>{ s.lower, fa, s.upper, fb, s.upper, fb, s.upper, fb },
                s.lower, s.upper, 0, true, detail::bracket_step(s.lower, fa, s.upper, fb) };
        }

        /**
         * @brief One evaluation: the algorithm's next secant, interpolation, double-length secant or bisection point.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param s The state.
         * @param t The tolerances.
         * @return Nothing, error::non_finite where fn is not finite at the point, or error::not_converged where no
         *         point lies strictly inside (not reached in the arithmetic solve() assumes).
         */
        template<real T, typename F>
        constexpr std::expected<void, error> advance(F &fn, detail::Toms748State<T> &s, Tolerance<T> const &t) const
        { return detail::toms748_advance(fn, s, t); }
    };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief False position with a scaling of the retained end's value, so that AndersonBjorck, Illinois and
         * RegulaFalsi share one implementation and differ only in the factor.
         * @tparam S The scaling.
         */
        template<Scaling S>
        struct FalsePositionMethod
        {
            /**
             * @brief The state on a validated bracket: the lower end retained, the upper the latest point.
             * @tparam T The real type.
             * @tparam F fn's type.
             * @param s The bracket.
             * @return The state.
             */
            template<real T, typename F>
            constexpr std::expected<FalsePositionState<T>, error>
                start(F &, Start<T> const &s, Tolerance<T> const &) const
            {
                T const fa = s.at_lower[0];
                T const fb = s.at_upper[0];
                return FalsePositionState<T>{ s.lower, fa, fa, s.upper, fb, detail::distance_bound(s.lower, s.upper), 0,
                    detail::bracket_step(s.lower, fa, s.upper, fb) };
            }

            /**
             * @brief One evaluation: the secant's zero through the latest point and the retained end's scaled value,
             * or a guard bisection.
             * @tparam T The real type.
             * @tparam F fn's type.
             * @param fn The function.
             * @param s The state.
             * @param t The tolerances.
             * @return Nothing, error::non_finite where fn is not finite at the point, or error::not_converged where no
             *         point lies strictly inside (not reached in the arithmetic solve() assumes).
             */
            template<real T, typename F>
            constexpr std::expected<void, error> advance(F &fn, FalsePositionState<T> &s, Tolerance<T> const &t) const
            { return detail::false_position_advance<S>(fn, s, t); }
        };

        /**
         * @brief A method with a derivative step (Newton's or Halley's): inside Brent's safeguard on a bracket, and
         * from a guess until two iterates show a sign change, then inside it, so that Newton and Halley share Brent's
         * bookkeeping and differ only in the step.
         * @tparam G The step.
         * @tparam K The values it needs: the value and K - 1 derivatives.
         */
        template<Safeguard G, std::size_t K>
        struct DerivativeMethod
        {
            /** @brief fn's values the method needs: the value and K - 1 derivatives (values_of). */
            static constexpr std::size_t values = K;

            /**
             * @brief The state on a validated bracket: Brent's a, b, c, d and e, with K values at each point.
             * @tparam T The real type.
             * @tparam F fn's type.
             * @param s The bracket.
             * @return The state.
             */
            template<real T, typename F>
            constexpr std::expected<SafeguardedState<T, K>, error>
                start(F &, Start<T, K> const &s, Tolerance<T> const &) const
            { return detail::safeguarded_start(s); }

            /**
             * @brief One iteration of Brent's procedure with the derivative step from b, where Brent's test accepts
             * it, else a bisection step.
             * @tparam T The real type.
             * @tparam F fn's type.
             * @param fn The function.
             * @param s The state.
             * @param t The tolerances.
             * @return Nothing, error::non_finite where fn is not finite at the point, or error::not_converged where no
             *         point lies strictly inside (not reached in the arithmetic solve() assumes).
             */
            template<real T, typename F>
            constexpr std::expected<void, error> advance(F &fn, SafeguardedState<T, K> &s, Tolerance<T> const &t) const
            { return detail::safeguarded_advance<G>(fn, s, t); }

            /**
             * @brief The state from a guess: fn's values there, the guess its point.
             * @tparam T The real type.
             * @tparam F fn's type.
             * @param fn The function.
             * @param guess The guess, finite.
             * @return The state, or error::non_finite where fn is not finite at the guess.
             */
            template<real T, typename F>
            constexpr std::expected<OpenDerivativeState<T, K>, error>
                start(F &fn, T const &guess, Tolerance<T> const &) const
            {
                std::array<T, K> const v = detail::evaluate<K>(fn, guess);
                if (!math::isfinite(v[0])) return std::unexpected(error::non_finite);
                return OpenDerivativeState<T, K>{ guess, v, std::nullopt, Step<T>{ guess, v[0], guess, v[0], guess } };
            }

            /**
             * @brief One iteration from a guess: the derivative step, at least tol(x) long, to the next iterate; at the
             * first sign change between two iterates, Brent's procedure with the same step on them; then its steps.
             * @tparam T The real type.
             * @tparam F fn's type.
             * @param fn The function.
             * @param s The state.
             * @param t The tolerances.
             * @return Nothing, error::non_finite where fn is not finite at the iterate, or error::not_converged where
             *         the step cannot be formed (f' = 0 or not finite, or x + step would overflow).
             */
            template<real T, typename F>
            constexpr std::expected<void, error>
                advance(F &fn, OpenDerivativeState<T, K> &s, Tolerance<T> const &t) const
            { return detail::open_derivative_advance<G>(fn, s, t); }
        };

    } // namespace detail
    /** @endcond */

    /**
     * @brief Regula falsi with Anderson and Bjorck's rescaling (1973): the secant through the latest point and the
     * retained end, whose value is scaled down by 1 - f(c) / f(b) (or 1/2) whenever the same end is retained again. A
     * bracket only.
     *
     * Guarded by a bisection after three steps that have not halved the bracket, so at most about 4 evaluations per
     * halving (the published method took 65 million evaluations on x^20 - 1 over [0, 5]; guarded, 30). Far cheaper than
     * Brent where fn is strongly curved near one end of the bracket (ln(x / 1e-12) on [1e-20, 1]: 18 evaluations against
     * 45) and on a walk to a sign change at 0; dearer on high powers and very wide brackets. Its point is the end with
     * the smaller |f|; each step is one evaluation.
     */
    struct AndersonBjorck : detail::FalsePositionMethod<detail::Scaling::anderson_bjorck>
    {
    };

    /**
     * @brief Ridders' method (1979): the midpoint, then the zero of an exponential fit through the three values. A
     * bracket. Not nxx::deriv::Ridders, the derivative method of the same author.
     *
     * Each iteration halves the bracket at a cost of 2 evaluations (so at most twice Bisection's count) and converges
     * quadratically per iteration on a simple root. A teaching example; of the interpolating methods the cheapest on
     * roots of odd multiplicity ((x - 1)^3: 86 evaluations against Brent's 156, though Bisection takes 54), dearer than
     * Brent otherwise (416 against 309 on 27 of Alefeld, Potra and Shi's test cases). The fit's zero is formed without
     * the cancellation of the textbook formula, which loses roots near an end of the bracket. Its point is the end with
     * the smaller |f|; each step is one evaluation, the midpoint or the fit's zero.
     */
    struct Ridders
    {
        /**
         * @brief The state on a validated bracket, the midpoint next.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param s The bracket.
         * @return The state.
         */
        template<real T, typename F>
        constexpr std::expected<detail::RiddersState<T>, error>
            start(F &, Start<T> const &s, Tolerance<T> const &) const
        {
            T const fa = s.at_lower[0];
            T const fb = s.at_upper[0];
            return detail::RiddersState<T>{ detail::Enclosure<T>{ s.lower, fa, s.upper, fb, s.upper, fb, s.upper, fb },
                s.lower, fa, 0, detail::bracket_step(s.lower, fa, s.upper, fb) };
        }

        /**
         * @brief One evaluation: the midpoint, or the zero of the fit through the bracket's ends and the midpoint.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param s The state.
         * @param t The tolerances.
         * @return Nothing, error::non_finite where fn is not finite at the point, or error::not_converged where no
         *         point lies strictly inside (not reached in the arithmetic solve() assumes).
         */
        template<real T, typename F>
        constexpr std::expected<void, error> advance(F &fn, detail::RiddersState<T> &s, Tolerance<T> const &t) const
        { return detail::ridders_advance(fn, s, t); }
    };

    /**
     * @brief False position: the zero of the secant through the bracket's ends, guarded. A bracket.
     *
     * A teaching example. As published, one end of the bracket stalls while the other converges: x^20 - 1 over [0, 5]
     * took over 200,000 evaluations. Guarded here, as AndersonBjorck, by a bisection after three steps that have not
     * halved the bracket (45 evaluations there), so the stall shows only in its counts: 645 on 27 of Alefeld, Potra and
     * Shi's test cases, against Illinois' 382 and AndersonBjorck's 327.
     */
    struct RegulaFalsi : detail::FalsePositionMethod<detail::Scaling::none>
    {
    };

    /**
     * @brief The Illinois method (Dowell and Jarratt 1971): false position with the retained end's value halved
     * whenever that end is retained again, the simplest cure for its stall; guarded, as AndersonBjorck. A bracket.
     */
    struct Illinois : detail::FalsePositionMethod<detail::Scaling::illinois>
    {
    };

    /**
     * @brief Newton's method: the step -f / f'. Needs fn and f' (values = 2).
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
    struct Newton : detail::DerivativeMethod<detail::Safeguard::newton, 2>
    {
    };

    /**
     * @brief Halley's method: the step -t / (1 - t f'' / (2 f')) with t = f / f'. Needs fn, f' and f'' (values = 3).
     *
     * Cubic near a simple root: x^2 - 2 from 1 in 5 evaluations. Where its correction cannot be formed or would point
     * against Newton's step (far from a root), it takes Newton's step, as Boost.Math's Halley does. On a bracket and
     * from a guess as Newton.
     */
    struct Halley : detail::DerivativeMethod<detail::Safeguard::halley, 3>
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
        /**
         * @brief The state from a guess: fn at the guess and at the second point.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param guess The guess, finite.
         * @param tolerance The tolerances.
         * @return The state, or error::non_finite where fn is not finite at either point.
         */
        template<real T, typename F>
        constexpr std::expected<detail::SecantState<T>, error>
            start(F &fn, T const &guess, Tolerance<T> const &tolerance) const
        {
            T const fp = detail::evaluate<1>(fn, guess)[0];
            if (!math::isfinite(fp)) return std::unexpected(error::non_finite);
            Step<T> const at_guess{ guess, fp, guess, fp, guess };
            if (fp == T(0)) return detail::SecantState<T>{ guess, fp, guess, fp, T(0), std::nullopt, at_guess };
            T const x = detail::secant_second_point(guess, tolerance);
            T const fx = detail::evaluate<1>(fn, x)[0];
            if (!math::isfinite(fx)) return std::unexpected(error::non_finite);
            detail::SecantState<T> s{ guess, fp, x, fx, T(x - guess), std::nullopt, at_guess };
            detail::secant_report(s);
            return s;
        }

        /**
         * @brief One iteration: the secant step through the latest two iterates, at least tol(x) long; at the first
         * sign change, Brent's procedure on them; then its steps.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param s The state.
         * @param t The tolerances.
         * @return Nothing, error::non_finite where fn is not finite at the iterate, or error::not_converged where the
         *         step cannot be formed (the two values equal, or a quotient, product or sum out of range).
         */
        template<real T, typename F>
        constexpr std::expected<void, error> advance(F &fn, detail::SecantState<T> &s, Tolerance<T> const &t) const
        { return detail::secant_advance(fn, s, t); }
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
        /**
         * @brief The state from a guess: fn's value there.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param guess The guess, finite.
         * @return The state, or error::non_finite where fn is not finite at the guess.
         */
        template<real T, typename F>
        constexpr std::expected<detail::SteffensenState<T>, error>
            start(F &fn, T const &guess, Tolerance<T> const &) const
        {
            T const fx = detail::evaluate<1>(fn, guess)[0];
            if (!math::isfinite(fx)) return std::unexpected(error::non_finite);
            return detail::SteffensenState<T>{ guess, fx, std::nullopt, Step<T>{ guess, fx, guess, fx, guess } };
        }

        /**
         * @brief One iteration, up to two evaluations: z = x + f(x), at least tol(x) from x, then the secant step
         * through x and z; at the first sign change, Brent's procedure on it; then its steps.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param s The state.
         * @param t The tolerances.
         * @return Nothing, error::non_finite where fn is not finite at a point, or error::not_converged where the
         *         step cannot be formed.
         */
        template<real T, typename F>
        constexpr std::expected<void, error> advance(F &fn, detail::SteffensenState<T> &s, Tolerance<T> const &t) const
        { return detail::steffensen_advance(fn, s, t); }
    };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief The type of fn the concepts assume by default, so that a method can be checked without a particular
         * fn: a plain function of T returning the values a method needs (T, or a std::array of K values of T).
         * @tparam T The real type.
         * @tparam K The number of values.
         */
        template<typename T, std::size_t K>
        using archetype = std::conditional_t<K == 1, T (*)(T const &), std::array<T, K> (*)(T const &)>;

        /**
         * @brief What a method's start returns for fn of type F from a start X, from which the state alias takes the
         * state.
         * @tparam M The method.
         * @tparam T The real type.
         * @tparam F fn's type.
         * @tparam X The start: a Start, or T for a guess.
         */
        template<typename M, typename T, typename F, typename X>
        using start_result = decltype(std::declval<M const &>().start(
            std::declval<F &>(), std::declval<X const &>(), std::declval<Tolerance<T> const &>()));

        /**
         * @brief Whether R is std::expected<S, error> for a state S that is not a reference.
         * @tparam R The type.
         */
        template<typename R>
        inline constexpr bool is_state_result = false;

        /**
         * @brief It is, for a state that is not a reference.
         * @tparam S The state.
         */
        template<typename S>
        inline constexpr bool is_state_result<std::expected<S, error>> = !std::is_reference_v<S>;

        /**
         * @brief A start's result: std::expected<S, error>, so that the concepts can require that shape without
         * naming the method's state type.
         * @tparam R The type.
         */
        template<typename R>
        concept state_result = is_state_result<R>;

    } // namespace detail
    /** @endcond */

    /**
     * @brief The state type of a bracketing method M for fn of type F, for a method of one's own that holds a
     * library method's state in its own: what `m.start(fn, start, tolerance)` returns inside its std::expected.
     * @tparam M The method.
     * @tparam T The real type.
     * @tparam F fn's type, reference removed; by default a plain function of T returning values_of<M> values.
     */
    template<typename M, typename T, typename F = detail::archetype<T, values_of<M>>>
    using bracketing_state_t = typename detail::start_result<M, T, F, Start<T, values_of<M>>>::value_type;

    /**
     * @brief What a type must provide to work with solve() on a bracket, one step at a time: the method computes and
     * reports, solve() checks each report and decides when to stop.
     *
     * The nine bracketing methods of the library are bracketing methods; a user's type can be another, given as a
     * type, `solve<M>(fn, bracket)`, or as a value, `solve(m, fn, bracket)`. Two const members, each returning a
     * std::expected with error:
     * - `m.start(fn, start, tolerance)` returns a state S of the method's own type, whose member `step` is the Step<T>
     *   after the start; `start` is a Start<T, values_of<M>>;
     * - `m.advance(fn, s, tolerance)` takes one more step (usually one evaluation of fn) and updates `s.step`.
     *
     * advance takes the state by non-const reference: one taken by value or by const reference fails the concept,
     * since the loop would then step a copy and never end.
     *
     * solve() runs the loop: after the start and after each advance it checks the step (all five numbers finite, else
     * error::non_finite; lower <= point <= upper, inside the caller's bracket, fn's values nonzero and of opposite
     * signs at lower < upper or exactly 0 at the point, and nested in the previous step's bracket: neither end
     * outside it, and not both ends the same, else error::invalid_step), then returns point where fn is 0 there
     * (abs_error tol(point)) or where the bracket is at most 2 tol(point) wide or point is within tol(point) of both
     * ends (abs_error the smaller of the bound on the width and the larger distance to an end, so at most 2
     * tol(point)). Since every bracket is nested in the last, every method ends.
     *
     * The contract, which solve() relies on:
     * - **Given:** the validated Tolerance; fn as an lvalue, valid for the call only (keep no reference, pointer or
     *   copy of it in the state); a Start with lower < upper, finite, fn's value finite, nonzero and of opposite
     *   signs at the ends, values_of<M> values at each, the derivatives unchecked. advance is called only on a step
     *   that does not meet the request; nothing is called after an error.
     * - **Required:** fn evaluated only inside the start's bracket, as the library does: with std::invoke and an
     *   argument of type T, its result converted to T at once; a value of fn that is not finite returned as
     *   error::non_finite; start and every advance leave `step` describing the state as above, with f_lower and
     *   f_upper fn's true values at lower and upper (not a scaled value, such as Pegasus's), which solve() cannot
     *   check and which make the result a verified sign change; error::not_converged from advance where the method
     *   cannot go on. Exceptions, fn's and the method's own, pass through.
     * - **For convergence, not checked:** keep each new point about tol(point) inside the ends. A method that moves
     *   only one end otherwise stalls at the stop test, the other end never moving.
     * - **Arithmetic:** IEEE round-to-nearest with subnormals, or a type that flushes below min() as
     *   Boost.Multiprecision's do; flush-to-zero modes (-ffast-math) are outside the derivations.
     *
     * Check a method where it is defined: `static_assert(nxx::roots::bracketing_method<MyMethod, double>)`; a failure
     * there names the cause (a member that is not const, a return type, the state's `step`), which solve()'s
     * diagnostic does not.
     *
     * @tparam M The method.
     * @tparam T The real type.
     * @tparam F fn's type, reference removed; by default a plain function of T returning values_of<M> values.
     */
    template<typename M, typename T, typename F = detail::archetype<T, values_of<M>>>
    concept bracketing_method =
        requires(M const &m, F &fn, Start<T, values_of<M>> const &s, Tolerance<T> const &t) {
            { m.start(fn, s, t) } -> detail::state_result;
        } && requires(M const &m, F &fn, bracketing_state_t<M, T, F> &st, Tolerance<T> const &t) {
            { m.advance(fn, st, t) } -> std::same_as<std::expected<void, error>>;
            { std::as_const(st).step } -> std::same_as<Step<T> const &>;
            // advance updates the state: one taken by value or by const reference would loop for ever on a copy
            // (PHASE7_API.md, Q17).
            requires !requires { m.advance(fn, std::move(st), t); };
        };

    /**
     * @brief The state type of an open method M for fn of type F, for a method of one's own that holds a library
     * method's state in its own: what `m.start(fn, guess, tolerance)` returns inside its std::expected.
     * @tparam M The method.
     * @tparam T The real type.
     * @tparam F fn's type, reference removed; by default a plain function of T returning values_of<M> values.
     */
    template<typename M, typename T, typename F = detail::archetype<T, values_of<M>>>
    using open_state_t = typename detail::start_result<M, T, F, T>::value_type;

    /**
     * @brief What a type must provide to work with solve() from a guess, one step at a time: the method iterates,
     * solve() checks each report and returns once two of its points bracket a sign change that meets the request.
     *
     * Newton, Halley, Secant and Steffensen are open methods; a user's type can be another, given as a type,
     * `solve<M>(fn, guess)`, or as a value, `solve(m, fn, guess)`. The members are those of bracketing_method, with
     * `m.start(fn, guess, tolerance)` taking a finite guess of type T in place of a Start.
     *
     * The step reports two points and fn's values there. Before a sign change they are the last two iterates, in order
     * (one point after a start that evaluates only the guess; Secant's holds its first two), and point the latest.
     * solve() checks each step as for bracketing_method, without a
     * caller's bracket: a step whose values have opposite signs at lower < upper is bracketed and is tested as on a
     * bracket, and a bracketed step after a bracketed one must be nested in it; fn exactly 0 at the point ends
     * the solve. A step without a sign change counts against the iteration cap, 4 times the bits of T (212 in double):
     * at the cap solve() returns error::not_converged. So every open method ends, and every success is a bracketed
     * step that met the request, a verified sign change.
     *
     * The contract: as for bracketing_method, and
     * - **Given:** a finite guess; fn may be evaluated anywhere.
     * - **Allowed:** a sign change may be lost again (a plain iteration need not keep it); the steps without one keep
     *   counting against the cap.
     * - **For convergence, not checked:** steps of at least tol(x), or an iterate converging from one side never
     *   brackets; and once a sign change is found, keep each point about tol inside it, as on a bracket. Holding a
     *   bracketing method's state (`std::optional<bracketing_state_t<Brent, T, F>>`) and forwarding to it after the
     *   first sign change is one way, which the library's open methods take.
     *
     * Check a method where it is defined: `static_assert(nxx::roots::open_method<MyMethod, double>)`.
     *
     * @tparam M The method.
     * @tparam T The real type.
     * @tparam F fn's type, reference removed; by default a plain function of T returning values_of<M> values.
     */
    template<typename M, typename T, typename F = detail::archetype<T, values_of<M>>>
    concept open_method = requires(M const &m, F &fn, T const &guess, Tolerance<T> const &t) {
        { m.start(fn, guess, t) } -> detail::state_result;
    } && requires(M const &m, F &fn, open_state_t<M, T, F> &st, Tolerance<T> const &t) {
        { m.advance(fn, st, t) } -> std::same_as<std::expected<void, error>>;
        { std::as_const(st).step } -> std::same_as<Step<T> const &>;
        // advance updates the state: one taken by value or by const reference would loop on a copy (PHASE7_API.md,
        // Q17).
        requires !requires { m.advance(fn, std::move(st), t); };
    };

    /**
     * @brief The observer's marker, nxx::Observer, shared with the other modules: solve() tells an observer from fn
     * by this type.
     */
    using nxx::Observer;

    /**
     * @brief nxx::observer, shared with the other modules: makes a callable the observer of a solve by a stepped
     * method, `o(step)` called with the Step<T> after the method's start and after each step.
     *
     * An observed run computes exactly what an unobserved one does. See nxx::observer for the rest.
     *
     * @param o A callable taking Step<T> const &, returning void.
     * @return The observer, for solve().
     */
    using nxx::observer;

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief Whether O is an Observer, so that the forms of solve() that take a method value do not take an
         * observer for one.
         * @tparam O The type.
         */
        template<typename O>
        inline constexpr bool is_observer = false;

        /**
         * @brief It is.
         * @tparam O The callable's type.
         */
        template<typename O>
        inline constexpr bool is_observer<Observer<O>> = true;

        /**
         * @brief The method solve() runs when it is given an observer: M wrapped so that each step is reported,
         * itself a stepped method with M's state, so that the observed run is the same computation and solve() needs
         * no second loop.
         * @tparam M The method.
         * @tparam O The observer's callable, as given to observer().
         */
        template<typename M, typename O>
        struct Observed
        {
            /**
             * @brief The method, held for the call to solve() that observes it, so that a method that cannot be
             * copied can be observed.
             */
            M const &method;

            /** @brief The observer, called as a non-const callable: a mutable lambda can count its calls. */
            mutable nxx::detail::kept<O> observer;

            /** @brief The method's values, forwarded, so that solve() evaluates what M needs. */
            static constexpr std::size_t values = values_of<M>;

            // Self = M makes each declaration dependent, so that an observer on a method without these members fails
            // the concept softly instead of where the class is instantiated.

            /**
             * @brief The method's start on a bracket, then a report of its first step.
             * @tparam T The real type.
             * @tparam K The number of values at each end.
             * @tparam F fn's type.
             * @tparam Self The method, as a dependent type.
             * @param fn The function.
             * @param s The bracket.
             * @param t The tolerances.
             * @return The method's state, or its error.
             */
            template<real T, std::size_t K, typename F, typename Self = M>
            requires bracketing_method<Self, T, F>
            constexpr auto start(F &fn, Start<T, K> const &s, Tolerance<T> const &t) const
                -> decltype(std::declval<Self const &>().start(fn, s, t))
            {
                auto started = method.start(fn, s, t);
                if (started) report(started->step);
                return started;
            }

            /**
             * @brief The method's start from a guess, then a report of its first step.
             * @tparam T The real type.
             * @tparam F fn's type.
             * @tparam Self The method, as a dependent type.
             * @param fn The function.
             * @param guess The guess.
             * @param t The tolerances.
             * @return The method's state, or its error.
             */
            template<real T, typename F, typename Self = M>
            requires open_method<Self, T, F>
            constexpr auto start(F &fn, T const &guess, Tolerance<T> const &t) const
                -> decltype(std::declval<Self const &>().start(fn, guess, t))
            {
                auto started = method.start(fn, guess, t);
                if (started) report(started->step);
                return started;
            }

            /**
             * @brief The method's step, then a report of it.
             * @tparam F fn's type.
             * @tparam S The method's state.
             * @tparam T The real type.
             * @tparam Self The method, as a dependent type.
             * @param fn The function.
             * @param s The state.
             * @param t The tolerances.
             * @return Nothing, or the method's error.
             */
            template<typename F, typename S, real T, typename Self = M>
            constexpr auto advance(F &fn, S &s, Tolerance<T> const &t) const
                -> decltype(std::declval<Self const &>().advance(fn, s, t))
            {
                auto const r = method.advance(fn, s, t);
                if (r) report(s.step);
                return r;
            }

            /**
             * @brief Calls the observer with a step, with one readable diagnostic where the observer cannot take it.
             * @tparam T The real type.
             * @param step The step.
             */
            template<real T>
            constexpr void report(Step<T> const &step) const
            {
                // Each requirement is checked only where the one before it holds, so a wrong observer gives one
                // diagnostic, not a cascade.
                constexpr bool callable = std::invocable<nxx::detail::kept<O> &, Step<T> const &>;
                static_assert(callable, "nxx::roots::observer(o): o must be callable as o(Step<T> const &)");
                if constexpr (callable) {
                    constexpr bool returns_void =
                        std::is_void_v<std::invoke_result_t<nxx::detail::kept<O> &, Step<T> const &>>;
                    static_assert(returns_void,
                        "nxx::roots::observer(o): o must return void (a result is reserved for stopping the run)");
                    if constexpr (returns_void) std::invoke(observer, step);
                }
            }
        };

        /**
         * @brief The method m observed by o, for the forms of solve() that take an observer.
         * @param m The method, held by reference: it must outlive the result (a temporary does, for the full
         * expression that solves with it).
         * @param o The observer.
         * @return The observed method.
         */
        template<typename M, typename O>
        constexpr Observed<M, O> observed(M const &m, Observer<O> &&o)
        { return Observed<M, O>{ m, std::move(o.callable) }; }

        /**
         * @brief The caller's bracket after validation, inside which every step and result must lie.
         * @tparam T The real type.
         */
        template<real T>
        struct Within
        {
            /** @brief The lower end. */
            T lower;

            /** @brief The upper end. */
            T upper;
        };

        /**
         * @brief The loop for a stepped method: check each step, return a root once it meets the request, else
         * advance (bracketing_method and open_method state the checks).
         *
         * The stop test is Criteria's rule, the bracket at most 2 tol(point) wide, or point within tol(point) of both
         * ends. For a method whose point is an end the two are one test; Bisection's midpoint can meet the second
         * where the unit distance_bound adds to a width that rounds fails the first (numerics review of Phase 7's
         * draft 1, finding 2). Each bracketed step after a bracketed step nested in it (neither end outside, not
         * both the same): finitely many such
         * brackets of T's numbers exist (the same review, finding 3); from a guess the steps without a sign change
         * count against the cap. So every method ends.
         *
         * @param m The method.
         * @param fn The function.
         * @param state The method's state, started.
         * @param tolerance The tolerances.
         * @param within The caller's bracket, for a bracketing method; nothing from a guess.
         * @return The root, or an error.
         */
        template<typename M, real T, typename F, typename S>
        constexpr std::expected<Root<T>, error>
            run(M const &m, F &fn, S &state, Tolerance<T> const &tolerance, std::optional<Within<T>> const &within)
        {
            std::optional<Within<T>> previous; // the last step's bracket, where it was one
            int unbracketed = 0;               // from a guess: the steps without a sign change so far
            while (true) {
                Step<T> const &s = std::as_const(state).step;
                if (!math::isfinite(s.lower) || !math::isfinite(s.upper) || !math::isfinite(s.point)
                    || !math::isfinite(s.f_lower) || !math::isfinite(s.f_upper))
                    return std::unexpected(error::non_finite);
                if (!(s.lower <= s.point && s.point <= s.upper)) return std::unexpected(error::invalid_step);
                if (within && (s.lower < within->lower || within->upper < s.upper))
                    return std::unexpected(error::invalid_step);
                if ((s.point == s.lower && s.f_lower == T(0)) || (s.point == s.upper && s.f_upper == T(0)))
                    return Root<T>{ s.point, tolerance.at(s.point) };
                bool const bracketed =
                    s.f_lower != T(0) && s.f_upper != T(0) && !detail::same_sign(s.f_lower, s.f_upper);
                if (bracketed) {
                    if (!(s.lower < s.upper)) return std::unexpected(error::invalid_step);
                    if (previous
                        && (s.lower < previous->lower || previous->upper < s.upper
                            || (s.lower == previous->lower && s.upper == previous->upper)))
                        return std::unexpected(error::invalid_step);
                    previous = Within<T>{ s.lower, s.upper };
                    T const tol = tolerance.at(s.point);
                    T const below = detail::distance_bound(s.point, s.lower);
                    T const above = detail::distance_bound(s.upper, s.point);
                    T const larger = below < above ? above : below;
                    T const width = detail::distance_bound(s.lower, s.upper);
                    // Both bound the distance from the point to the sign change; the smaller keeps Root's bound,
                    // 2 tol, where the width test passed but a point next to an end rounds its far distance up
                    // (independent review of Phase 7, finding 1).
                    if (width <= detail::twice(tol) || larger <= tol)
                        return Root<T>{ s.point, width < larger ? width : larger };
                } else {
                    // On a bracket the sign change must be kept; from a guess a step without one counts against the
                    // cap, so that a method that finds and loses sign changes still ends (simplicity review of
                    // draft 1, finding 3).
                    if (within) return std::unexpected(error::invalid_step);
                    previous.reset();
                    if (unbracketed == detail::iteration_cap<T>()) return std::unexpected(error::not_converged);
                    ++unbracketed;
                }
                std::expected<void, error> const r = m.advance(fn, state, tolerance);
                if (!r) return std::unexpected(r.error());
            }
        }

        /**
         * @brief solve() on a bracket: the tolerances and the bracket validated, zero first, then the method's start
         * and the loop.
         * @param m The method.
         * @param fn The function.
         * @param bracket The caller's bracket.
         * @param criteria The caller's criteria.
         * @return The root, or an error.
         */
        template<typename M, real T, typename F>
        constexpr std::expected<Root<T>, error>
            on_bracket(M const &m, F &fn, Bracket<T> const &bracket, Criteria<T> const &criteria)
        {
            constexpr std::size_t K = values_of<M>;
            // Qualified, so that argument-dependent lookup cannot pick a function from the namespace of T or F.
            std::expected<Tolerance<T>, error> const tolerance = detail::checked_tolerance(criteria);
            if (!tolerance) return std::unexpected(tolerance.error());
            std::expected<Prepared<T, K>, error> const prepared = detail::prepare<K>(fn, bracket, *tolerance);
            if (!prepared) return std::unexpected(prepared.error());
            if (prepared->root) return *prepared->root;
            Start<T, K> const &start = prepared->start;
            auto state = m.start(fn, start, *tolerance);
            if (!state) return std::unexpected(state.error());
            std::optional<Within<T>> const within = Within<T>{ start.lower, start.upper };
            return detail::run(m, fn, *state, *tolerance, within);
        }

        /**
         * @brief solve() from a guess: the tolerances and the guess validated, then the method's start and the loop.
         * @param m The method.
         * @param fn The function.
         * @param guess The caller's guess.
         * @param criteria The caller's criteria.
         * @return The root, or an error.
         */
        template<typename M, real T, typename F>
        constexpr std::expected<Root<T>, error> on_guess(M const &m, F &fn, T const &guess, Criteria<T> const &criteria)
        {
            std::expected<Tolerance<T>, error> const tolerance = detail::checked_tolerance(criteria);
            if (!tolerance) return std::unexpected(tolerance.error());
            if (!math::isfinite(guess)) return std::unexpected(error::non_finite);
            auto state = m.start(fn, guess, *tolerance);
            if (!state) return std::unexpected(state.error());
            return detail::run(m, fn, *state, *tolerance, std::optional<Within<T>>());
        }

    } // namespace detail
    /** @endcond */

    /**
     * @brief A root of @p fn in @p bracket by the method @p M: a bracketing method (Bisection, Brent, Toms748,
     * AndersonBjorck, Ridders, RegulaFalsi or Illinois), Newton or Halley safeguarded by the bracket, or a
     * bracketing_method of one's own.
     *
     * The bracket's ends may come in either order; fn must have opposite signs at them, or be 0 at one. When the
     * bracket contains 0 in its interior, fn is also evaluated at 0 first: a root at 0 is then found at once, and the
     * bracket is halved at 0 (skipped where fn(0) is not finite). The method then narrows the bracket until it is at
     * most 2 tol(x) wide (see Criteria), or fn is exactly 0 at a point. Newton and Halley take their derivatives from
     * fn, which then returns a std::pair, std::tuple or std::array of exactly T (the value, f', f''); the other methods
     * accept such a fn too and use its first element.
     *
     * A pole inside the bracket changes sign like a root, and is returned as one (a later version may report a pole as
     * an error). |fn(value)| alone does not tell them apart: a well-conditioned root can have a large |fn| nearby, a
     * pole a small one. Where fn may have poles, compare |fn(value)| with |fn| a little to each side, at value +-
     * 1024 abs_error: it grows away from a root and falls away from a pole. A point that lands exactly on a pole, where
     * fn is infinite, is error::non_finite. fn's exceptions pass through unchanged.
     *
     * Choosing a method: Brent by default. Toms748 costs about the same, with a worst case of about 4 evaluations per
     * halving of the bracket where Brent's grows with the square of the halvings. AndersonBjorck where fn is strongly
     * curved near one end of the bracket (a logarithm, 1/x) or the sign change is at 0. Bisection where fn is rough,
     * its count must be predictable, or the root has odd multiplicity ((x - 1)^3: 54 evaluations, against Ridders' 86
     * and Brent's 156). Ridders, RegulaFalsi and Illinois are teaching examples. For a value-only fn, not a derivative
     * method with a numerical derivative.
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
     *           `Newton`, `Halley`, or a default-constructible bracketing_method (a method with parameters is given
     *           as a value).
     * @tparam T The real type, deduced from @p bracket.
     * @tparam F A callable taking T const & and returning T (for a class type T, T itself or a built-in arithmetic
     *           type, not an expression template; for a built-in T, anything convertible to T), or (for Newton and
     *           Halley, which need it) a std::pair, std::tuple or std::array of exactly T: the value and the
     *           derivatives. For Newton on a fn without a derivative, pass the derivative as a callable after fn, use
     *           nxx::deriv::derivative_of, or use Secant or Brent.
     * @param fn The function.
     * @param bracket An interval in which fn changes sign.
     * @param criteria When to stop; it converts to T, so `{.rel_tol = 1e-10}` works for any T.
     * @return The root and a bound on its distance to a sign change of fn; error::non_finite if an end, a tolerance or
     *         a value of fn is not finite or a tolerance is negative; error::no_sign_change if fn has the same sign at
     *         both ends; error::not_converged if the method cannot go on (not reached by the library's methods);
     *         error::invalid_step if a method of one's own broke the contract of bracketing_method. Any other error a
     *         method of one's own returns passes through.
     */
    // criteria is a non-deduced context (type_identity_t), so T comes from the bracket alone and a braced list of
    // doubles converts to a multiprecision T.
    template<typename M, real T, typename F>
    requires std::default_initializable<M> && bracketing_method<M, T, std::remove_reference_t<F>>
             && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_bracket(M{}, fn, bracket, criteria); }

    /**
     * @brief A root of @p fn in @p bracket by the method @p m, given as a value: for a method of one's own with
     * parameters, as solve<M>(fn, bracket) otherwise.
     *
     * @code
     * namespace rt = nxx::roots;
     * auto const r = rt::solve(rt::Toms748{}, fn, rt::Bracket{1.0, 2.0}, {.rel_tol = 1e-12});
     * @endcode
     *
     * @tparam M The method: a bracketing_method.
     * @tparam T The real type, deduced from @p bracket.
     * @tparam F As for solve<M>(fn, bracket).
     * @param m The method.
     * @param fn The function.
     * @param bracket An interval in which fn changes sign.
     * @param criteria When to stop.
     * @return As for solve<M>(fn, bracket).
     */
    template<typename M, real T, typename F>
    requires(!detail::is_observer<M>) && bracketing_method<M, T, std::remove_reference_t<F>>
            && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(M const &m, F &&fn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_bracket(m, fn, bracket, criteria); }

    /**
     * @brief As solve<M>(fn, bracket), each step reported to an observer: `o(step)` with the Step<T> after the
     * method's start and after each step, the ends in order.
     *
     * The observed run computes exactly what an unobserved one does. An lvalue observer is kept by reference, an
     * rvalue by value; it is called as a non-const callable, must return void, and its exceptions pass through.
     *
     * @code
     * namespace rt = nxx::roots;
     * auto const print = [](rt::Step<double> const &s) { std::println("[{}, {}] {}", s.lower, s.upper, s.point); };
     * auto const r = rt::solve<rt::Brent>(rt::observer(print), fn, rt::Bracket{1.0, 2.0});
     * @endcode
     *
     * @tparam M The method, as for solve<M>(fn, bracket).
     * @tparam O The observer's callable.
     * @tparam T The real type, deduced from @p bracket.
     * @tparam F As for solve<M>(fn, bracket).
     * @param o The observer, made by observer().
     * @param fn The function.
     * @param bracket An interval in which fn changes sign.
     * @param criteria When to stop.
     * @return As for solve<M>(fn, bracket).
     */
    template<typename M, typename O, real T, typename F>
    requires std::default_initializable<M> && bracketing_method<M, T, std::remove_reference_t<F>>
             && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(Observer<O> o, F &&fn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_bracket(detail::observed(M{}, std::move(o)), fn, bracket, criteria); }

    /**
     * @brief As solve(m, fn, bracket), each step reported to an observer, as for solve<M>(observer(o), fn, bracket).
     * @tparam M The method: a bracketing_method.
     * @tparam O The observer's callable.
     * @tparam T The real type, deduced from @p bracket.
     * @tparam F As for solve<M>(fn, bracket).
     * @param m The method.
     * @param o The observer, made by observer().
     * @param fn The function.
     * @param bracket An interval in which fn changes sign.
     * @param criteria When to stop.
     * @return As for solve<M>(fn, bracket).
     */
    template<typename M, typename O, real T, typename F>
    requires bracketing_method<M, T, std::remove_reference_t<F>> && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error> solve(M const &m,
        Observer<O> o,
        F &&fn,
        Bracket<T> const &bracket,
        Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_bracket(detail::observed(m, std::move(o)), fn, bracket, criteria); }

    /**
     * @brief A root of @p fn from @p guess by the open method @p M: Newton, Halley, Secant, Steffensen, or an
     * open_method of one's own.
     *
     * The method iterates from the guess until fn's values at two consecutive iterates have opposite signs, then
     * continues inside that bracket (Newton and Halley safeguarded, Secant and Steffensen as Brent), so the result is a
     * verified sign change with a bound, as for the bracketing methods. It returns a root, the first sign change its
     * iterates cross, not necessarily the one nearest the guess; from a guess already within tol of a root it takes one
     * step past it (2 evaluations). error::not_converged where it finds no sign change within the cap (4 times the bits
     * of T: 212 steps in double) or cannot take a step: a cycle, a divergence, f' = 0, a multiple root at 0 with the
     * default abs_tol. A bracket solves those.
     *
     * @code
     * namespace rt = nxx::roots;
     * auto const fdf = [](double x) { return std::pair{x * x - 2.0, 2.0 * x}; };
     * auto const n = rt::solve<rt::Newton>(fdf, 1.0);                              // 7 evaluations
     * auto const s = rt::solve<rt::Secant>([](double x) { return x * x - 2.0; }, 1.0); // 9
     * @endcode
     *
     * @tparam M The method: `Newton`, `Halley`, `Secant`, `Steffensen`, or a default-constructible open_method (a
     *           method with parameters is given as a value).
     * @tparam T The real type, deduced from @p guess.
     * @tparam F A callable taking T const &: for Secant and Steffensen returning T (for a class type T, T itself or a
     *           built-in arithmetic type, not an expression template; for a built-in T, anything convertible to T),
     *           for Newton and Halley a std::pair, std::tuple or std::array of exactly T (the value, f', f''). For
     *           Newton on a fn without a derivative, pass the derivative as a callable after fn, use
     *           nxx::deriv::derivative_of, or use Secant.
     * @param fn The function.
     * @param guess The starting point.
     * @param criteria When to stop; it converts to T.
     * @return The root and a bound on its distance to a sign change of fn; error::non_finite if the guess, a tolerance
     *         or a value of fn is not finite or a tolerance is negative; error::not_converged if the method found no
     *         root; error::invalid_step if a method of one's own broke the contract of open_method. Any other error a
     *         method of one's own returns passes through.
     */
    template<typename M, real T, typename F>
    requires std::default_initializable<M> && open_method<M, T, std::remove_reference_t<F>>
             && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_guess(M{}, fn, guess, criteria); }

    /**
     * @brief A root of @p fn from @p guess by the open method @p m, given as a value: for a method of one's own with
     * parameters, as solve<M>(fn, guess) otherwise.
     * @tparam M The method: an open_method.
     * @tparam T The real type, deduced from @p guess.
     * @tparam F As for solve<M>(fn, guess).
     * @param m The method.
     * @param fn The function.
     * @param guess The starting point.
     * @param criteria When to stop.
     * @return As for solve<M>(fn, guess).
     */
    template<typename M, real T, typename F>
    requires(!detail::is_observer<M>) && open_method<M, T, std::remove_reference_t<F>>
            && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(M const &m, F &&fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_guess(m, fn, guess, criteria); }

    /**
     * @brief As solve<M>(fn, guess), each step reported to an observer: `o(step)` with the Step<T> after the method's
     * start and after each step.
     *
     * Before the first sign change the step holds the last two iterates, in order (one point after a start that
     * evaluates only the guess; Secant's start holds its first two), its point
     * the latest; after it, the bracket. The observed run computes exactly what an unobserved one does.
     *
     * @tparam M The method, as for solve<M>(fn, guess).
     * @tparam O The observer's callable.
     * @tparam T The real type, deduced from @p guess.
     * @tparam F As for solve<M>(fn, guess).
     * @param o The observer, made by observer().
     * @param fn The function.
     * @param guess The starting point.
     * @param criteria When to stop.
     * @return As for solve<M>(fn, guess).
     */
    template<typename M, typename O, real T, typename F>
    requires std::default_initializable<M> && open_method<M, T, std::remove_reference_t<F>>
             && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(Observer<O> o, F &&fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_guess(detail::observed(M{}, std::move(o)), fn, guess, criteria); }

    /**
     * @brief As solve(m, fn, guess), each step reported to an observer, as for solve<M>(observer(o), fn, guess).
     * @tparam M The method: an open_method.
     * @tparam O The observer's callable.
     * @tparam T The real type, deduced from @p guess.
     * @tparam F As for solve<M>(fn, guess).
     * @param m The method.
     * @param o The observer, made by observer().
     * @param fn The function.
     * @param guess The starting point.
     * @param criteria When to stop.
     * @return As for solve<M>(fn, guess).
     */
    template<typename M, typename O, real T, typename F>
    requires open_method<M, T, std::remove_reference_t<F>> && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error> solve(M const &m,
        Observer<O> o,
        F &&fn,
        T const &guess,
        Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_guess(detail::observed(m, std::move(o)), fn, guess, criteria); }

    /**
     * @brief Newton on @p bracket with the derivative as a separate callable: the same as Newton on a fn returning
     * (f, f'), with each callable evaluated once per point. Also for a method of one's own that declares
     * `values = 2`.
     *
     * The derivative may return a value, as fn does, or a std::expected whose value converts to T (its error is never
     * read), such as the callable nxx::deriv::derivative_of(fn) or a lambda around nxx::deriv::derivative(fn, x, h)
     * with a step of one's own. A std::expected without a value counts as a derivative that is not finite: on a
     * bracket the step is a
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
     * The separate-derivative forms take the method by name only: with a method value or an observer, fn is one
     * callable returning the values, as in the last example.
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
     * value, not a tuple or a std::expected), detail::derivative_valued (the derivatives), bracketing_method and
     * open_method (the method and the start).
     * @tparam M Newton, or a default-constructible bracketing_method with `values = 2`.
     * @tparam T The real type, deduced from @p bracket.
     * @tparam F A callable taking T const & and returning T (for a class type T, T itself or a built-in arithmetic
     *           type, not an expression template; for a built-in T, anything convertible to T).
     * @tparam D The derivative: a callable as F, or one returning a std::expected whose value converts to T.
     * @param fn The function.
     * @param dfn Its derivative.
     * @param bracket An interval in which fn changes sign.
     * @param criteria When to stop.
     * @return As for a bracket.
     */
    template<typename M, real T, typename F, typename D>
    requires(values_of<M> == 2) && detail::value_valued<F, T> && detail::derivative_valued<D, T>
            && std::default_initializable<M>
            && bracketing_method<M, T, detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, D &&dfn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const joined{ std::addressof(fn),
            { std::addressof(dfn) } };
        return nxx::roots::solve<M>(joined, bracket, criteria);
    }

    /**
     * @brief Newton from @p guess with the derivative as a separate callable. Also for an open method of one's own
     * that declares `values = 2`.
     * @tparam M Newton, or a default-constructible open_method with `values = 2`.
     * @tparam T The real type, deduced from @p guess.
     * @tparam F A callable taking T const & and returning T (for a class type T, T itself or a built-in arithmetic
     *           type, not an expression template; for a built-in T, anything convertible to T).
     * @tparam D The derivative: a callable as F, or one returning a std::expected whose value converts to T.
     * @param fn The function.
     * @param dfn Its derivative.
     * @param guess The starting point.
     * @param criteria When to stop.
     * @return As from a guess.
     */
    template<typename M, real T, typename F, typename D>
    requires(values_of<M> == 2) && detail::value_valued<F, T> && detail::derivative_valued<D, T>
            && std::default_initializable<M>
            && open_method<M, T, detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, D &&dfn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const joined{ std::addressof(fn),
            { std::addressof(dfn) } };
        return nxx::roots::solve<M>(joined, guess, criteria);
    }

    /**
     * @brief Halley on @p bracket with the first and second derivatives as separate callables. Also for a method of
     * one's own that declares `values = 3`.
     * @tparam M Halley, or a default-constructible bracketing_method with `values = 3`.
     * @tparam T The real type, deduced from @p bracket.
     * @tparam F A callable taking T const & and returning T (for a class type T, T itself or a built-in arithmetic
     *           type, not an expression template; for a built-in T, anything convertible to T).
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
    requires(values_of<M> == 3) && detail::value_valued<F, T> && detail::derivative_valued<D, T>
            && detail::derivative_valued<D2, T> && std::default_initializable<M>
            && bracketing_method<M,
                T,
                detail::Joined<T,
                    std::remove_reference_t<F>,
                    std::remove_reference_t<D>,
                    std::remove_reference_t<D2>> const>
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
     * @brief Halley from @p guess with the first and second derivatives as separate callables. Also for an open
     * method of one's own that declares `values = 3`.
     * @tparam M Halley, or a default-constructible open_method with `values = 3`.
     * @tparam T The real type, deduced from @p guess.
     * @tparam F A callable taking T const & and returning T (for a class type T, T itself or a built-in arithmetic
     *           type, not an expression template; for a built-in T, anything convertible to T).
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
    requires(values_of<M> == 3) && detail::value_valued<F, T> && detail::derivative_valued<D, T>
            && detail::derivative_valued<D2, T> && std::default_initializable<M>
            && open_method<M,
                T,
                detail::Joined<T,
                    std::remove_reference_t<F>,
                    std::remove_reference_t<D>,
                    std::remove_reference_t<D2>> const>
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
     * What it cannot find: a root of even multiplicity (unless fn is exactly 0 at a point it tries), two roots within
     * one step (x^2 - 2 from 100), a root beyond a region where fn is not finite on the same side, and a sign change
     * beyond the cap. A pole brackets like a root (tan x from 1.4), and solve() then returns the pole; solve()'s
     * comment says how to tell a pole from a root. The bracket need not hold the root nearest x0,
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
     * @tparam F A callable taking T const & and returning T (for a class type T, T itself or a built-in arithmetic
     *           type, not an expression template; for a built-in T, anything convertible to T), or a std::pair,
     *           std::tuple or std::array of exactly T, whose first element is used.
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
            case nxx::roots::error::invalid_step:
                name = "invalid_step";
                break;
            }
            return formatter<string_view>::format(name, ctx);
        }
    };

} // namespace std
