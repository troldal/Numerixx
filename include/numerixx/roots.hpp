/**
 * @file
 * @brief Roots of a real function of one variable: bracketing methods.
 *
 * solve<M>(fn, bracket) finds a point where fn changes sign inside a bracket, an interval at whose ends fn has opposite
 * signs, by a method M: `Bisection`, which halves the bracket, or `Brent`, which interpolates and falls back on
 * bisection. It returns a std::expected<Root<T>, error>: the root and a bound on its distance to a sign change of fn.
 * Numerixx itself never throws, but the result's own value() throws std::bad_expected_access on a failure, so test
 * the result and read it with `->` or `*`.
 *
 * @code
 * namespace rt = nxx::roots;
 * auto const r = rt::solve<rt::Brent>([](double x) { return x * x - 2.0; }, rt::Bracket{1.0, 2.0});
 * if (r) use(r->value, r->abs_error); // sqrt(2) to about 1e-16, after 9 evaluations of fn
 * @endcode
 */

#pragma once

#include <numerixx/math.hpp>

#include <cassert>
#include <concepts>
#include <expected>
#include <format>
#include <functional>
#include <limits>
#include <optional>
#include <string_view>
#include <type_traits>

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
         * A bracket end, a value of fn or a tolerance is not finite, or a tolerance is negative. A value of fn that is
         * infinite at a bracket end is reported here; a later version may accept it for its sign.
         */
        non_finite,

        /** fn has the same sign at both ends of the bracket, and is 0 at neither. */
        no_sign_change
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
     * 2 tol(value) (see Criteria). It saturates at max where that distance exceeds what T holds (only with tolerances
     * near max). The root of the exact function may differ by fn's own rounding error divided by
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

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief The bracketing methods solve() accepts: a closed set.
         * @tparam M The method.
         */
        template<typename M>
        concept bracketing = std::same_as<M, Bisection> || std::same_as<M, Brent>;

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
         * @brief fn at x, converted to T before any arithmetic (expression templates do not escape).
         * @param fn The function.
         * @param x The point, a T const lvalue as the constraint promises.
         * @return fn(x) as a T.
         */
        template<real T, typename F>
        constexpr T evaluate(F &fn, T const &x)
        { return static_cast<T>(std::invoke(fn, x)); }

        /**
         * @brief A validated bracket, in order, with fn's values at its ends; or a root found during validation.
         * @tparam T The real type.
         */
        template<real T>
        struct Start
        {
            /** @brief The root, if validation found one (fn exactly 0 at an end or at 0). */
            std::optional<Root<T>> root;

            /** @brief The lower end. */
            T lower;

            /** @brief fn at the lower end, finite and nonzero. */
            T f_lower;

            /** @brief The upper end. */
            T upper;

            /** @brief fn at the upper end, finite, nonzero and of the other sign. */
            T f_upper;
        };

        /**
         * @brief Validates the bracket, puts it in order, evaluates fn at its ends and, for a bracket around 0, at 0.
         *
         * Zero first: a root at 0 costs one evaluation, where the default tolerance would otherwise approach it down
         * to T's smallest numbers, and the rest of the method sees a bracket of one sign. Skipped where fn(0) is not
         * finite, since removable singularities at 0 (sin(x)/x, x/(e^x - 1)) are common and an evaluation the caller
         * did not ask for must not fail the solve.
         *
         * @param fn The function.
         * @param bracket The caller's bracket.
         * @param tolerance The tolerances.
         * @return The start, or error::non_finite or error::no_sign_change.
         */
        template<real T, typename F>
        constexpr std::expected<Start<T>, error>
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
            T f_lower = detail::evaluate(fn, lower);
            T f_upper = detail::evaluate(fn, upper);
            if (!math::isfinite(f_lower) || !math::isfinite(f_upper)) return std::unexpected(error::non_finite);
            if (f_lower == T(0))
                return Start<T>{ Root<T>{ lower, tolerance.at(lower) }, lower, f_lower, upper, f_upper };
            if (f_upper == T(0))
                return Start<T>{ Root<T>{ upper, tolerance.at(upper) }, lower, f_lower, upper, f_upper };
            if (detail::same_sign(f_lower, f_upper)) return std::unexpected(error::no_sign_change);
            if (lower < T(0) && T(0) < upper) {
                T const zero = T(0);
                T const f_zero = detail::evaluate(fn, zero);
                if (f_zero == T(0))
                    return Start<T>{ Root<T>{ zero, tolerance.at(zero) }, lower, f_lower, upper, f_upper };
                if (math::isfinite(f_zero)) {
                    if (detail::same_sign(f_zero, f_lower)) {
                        lower = zero;
                        f_lower = f_zero;
                    } else {
                        upper = zero;
                        f_upper = f_zero;
                    }
                }
            }
            return Start<T>{ std::nullopt, lower, f_lower, upper, f_upper };
        }

        /**
         * @brief Bisection on a validated bracket.
         * @param fn The function.
         * @param s The start: a bracket in order with a sign change and no root found.
         * @param tolerance The tolerances.
         * @return The root, or error::non_finite if fn is not finite at a midpoint.
         */
        template<real T, typename F>
        constexpr std::expected<Root<T>, error> bisection(F &fn, Start<T> const &s, Tolerance<T> const &tolerance)
        {
            T a = s.lower;
            T fa = s.f_lower;
            T b = s.upper;
            while (true) {
                T const m = detail::midpoint(a, b);
                // The larger of the distances to the ends is a true bound even where half the width rounds, as among
                // the subnormals or with ends far apart. A midpoint not strictly inside means T's resolution.
                T const left = detail::distance_bound(m, a);
                T const right = detail::distance_bound(b, m);
                T const half = right < left ? left : right;
                if (!(a < m && m < b) || half <= tolerance.at(m)) return Root<T>{ m, half };
                T const fm = detail::evaluate(fn, m);
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
         * @brief Brent's procedure zero (Brent 1973, ch. 4, section 6) on a validated bracket.
         *
         * Implemented from the book, with three changes: the tolerance is the caller's tol(b) rather than
         * 2 macheps |b| + t; the stop test is distance_bound(c, b) <= 2 tol(b) rather than |c - b| / 2 <= tol(b),
         * since the bound is exact there while the half may round in radix 10; and a step is a bisection step whenever
         * |c - b|, |b - a| or |e| exceeds max / 16, so that no product of the interpolation overflows (3 m q reaches 12
         * |c - b|, 2 p 4 (|c - b| + |b - a|)). The quotients cannot overflow: in the interpolation |fb| < |fa| <= |fc|.
         *
         * @param fn The function.
         * @param s The start: a bracket in order with a sign change and no root found.
         * @param tolerance The tolerances.
         * @return The root, or error::non_finite if fn is not finite at a point.
         */
        template<real T, typename F>
        constexpr std::expected<Root<T>, error> brent(F &fn, Start<T> const &s, Tolerance<T> const &tolerance)
        {
            T const big = T(detail::largest<T>() / T(16));

            T a = s.lower;
            T fa = s.f_lower;
            T b = s.upper;
            T fb = s.f_upper;
            T c = a;
            T fc = fa;
            T d = detail::bounded_difference(b, a, big);
            T e = d;
            while (true) {
                // [b, c] is the bracket, with fb and fc of opposite signs; a is the previous b.
                if (detail::same_sign(fb, fc)) {
                    c = a;
                    fc = fa;
                    d = detail::bounded_difference(b, a, big);
                    e = d;
                }
                if (math::abs(fc) < math::abs(fb)) {
                    a = b;
                    b = c;
                    c = a;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }
                T const tol = tolerance.at(b);
                T const width = detail::distance_bound(c, b);
                if (width <= detail::twice(tol)) return Root<T>{ b, width };

                bool const wide = width > big || detail::distance_bound(b, a) > big || math::abs(e) > big;
                // Half the bracket, towards c: exact Brent below big, the overflow-free midpoint above it.
                T const m = wide ? T(detail::midpoint(b, c) - b) : T(T(c - b) / T(2));
                if (wide || math::abs(e) < tol || math::abs(fa) <= math::abs(fb)) {
                    d = m;
                    e = m;
                } else {
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
                }
                a = b;
                fa = fb;
                if (math::abs(d) > tol)
                    b = T(b + d);
                else if (m > T(0))
                    b = T(b + tol);
                else
                    b = T(b - tol);
                fb = detail::evaluate(fn, b);
                if (!math::isfinite(fb)) return std::unexpected(error::non_finite);
                if (fb == T(0)) return Root<T>{ b, tolerance.at(b) };
            }
        }

    } // namespace detail
    /** @endcond */

    /**
     * @brief A root of @p fn in @p bracket by the bracketing method @p M.
     *
     * The bracket's ends may come in either order; fn must have opposite signs at them, or be 0 at one. When the
     * bracket contains 0 in its interior, fn is also evaluated at 0 first: a root at 0 is then found at once, and the
     * bracket is halved at 0 (skipped where fn(0) is not finite). The method then narrows the bracket until it is at
     * most 2 tol(x) wide (see Criteria), or fn is exactly 0 at a point.
     *
     * A pole inside the bracket changes sign like a root, and is returned as one: where fn may have poles, check
     * |fn(value)| (a later version may report a pole as an error). fn's exceptions pass through unchanged.
     *
     * @code
     * namespace rt = nxx::roots;
     * auto const r = rt::solve<rt::Brent>([](double t) { return psat(t) - 101325.0; }, rt::Bracket{300.0, 400.0});
     * auto const s = rt::solve<rt::Bisection>(fn, rt::Bracket{0.0, 1.0}, {.rel_tol = 1e-8});
     * @endcode
     *
     * @tparam M The method: `Bisection` or `Brent`.
     * @tparam T The real type, deduced from @p bracket.
     * @tparam F A callable taking T const & and returning something convertible to T.
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
    requires detail::bracketing<M> && nxx::detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        // Qualified, so that argument-dependent lookup cannot pick a function from the namespace of T or F.
        std::expected<detail::Tolerance<T>, error> const tolerance = detail::checked_tolerance(criteria);
        if (!tolerance) return std::unexpected(tolerance.error());
        std::expected<detail::Start<T>, error> const start = detail::prepare(fn, bracket, *tolerance);
        if (!start) return std::unexpected(start.error());
        if (start->root) return *start->root;
        if constexpr (std::same_as<M, Bisection>)
            return detail::bisection(fn, *start, *tolerance);
        else
            return detail::brent(fn, *start, *tolerance);
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
            }
            return formatter<string_view>::format(name, ctx);
        }
    };

} // namespace std
