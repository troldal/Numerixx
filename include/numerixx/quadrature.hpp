/**
 * @file
 * @brief Integrals of a real function of one variable, computed from its values: the module nxx::quadrature.
 *
 * Use it when you need the integral of f from a to b and cannot write down an antiderivative, or f exists only as a
 * program: a callable that returns f(x). integrate(fn, a, b) returns the integral and an estimate of its error by the
 * method that suits most integrands on a finite interval; integrate<M>(fn, a, b) chooses another method M:
 * - `GaussKronrod<15>`, the default: adaptive, it subdivides the interval where fn is hard, so it copes with kinks,
 *   peaks and mild endpoint singularities without being told where they are.
 * - `TanhSinh`: for a singularity at an end (log x, 1 / sqrt x), for many digits in a multiprecision type, and the one
 *   method for an infinite interval.
 * - `Romberg`: extrapolated trapezoidal sums, for teaching and comparison.
 * - `GaussLegendre<N>`: a fixed rule of N points, with no error estimate and no tolerance, so that its result is a
 *   smooth function of the parameters inside fn: the rule for an integral inside a root solve, an optimiser or another
 *   integral.
 *
 * integrate(m, fn, a, b) takes the method as a value instead, for a method with parameters. The methods with an
 * estimate return an Integral (the value and its estimated error) or a Failure (why, and the best estimate reached);
 * GaussLegendre returns a plain number or an error. Every form returns a std::expected. Numerixx itself never throws,
 * except std::bad_alloc where GaussKronrod, or TanhSinh in a type without a compile-time table, cannot allocate; the
 * result's own value() throws std::bad_expected_access on a failure, so test the result and read it with `->` or `*`.
 *
 * A method is a type, so a type of one's own works as the library's do: a fixed_rule, an iterative_method or an
 * estimating_method, whose concepts state the contract. An iterative method (GaussKronrod<15>, TanhSinh, Romberg)
 * computes one step at a time and reports a Step; integrate() runs the loop and decides when the criteria are met, and
 * observer(o) shows each Step to o, to watch a method converge.
 *
 * @code
 * namespace qd = nxx::quadrature;
 * auto const r = qd::integrate([](double x) { return std::exp(x); }, 0.0, 1.0); // e - 1, 150 evaluations
 * if (r) use(r->value, r->abs_error);
 * auto const s = qd::integrate<qd::Romberg>(qd::observer([](qd::Step<double> const &step) { log(step); }),
 *                                           [](double x) { return std::exp(x); }, 0.0, 1.0);
 * auto const g = qd::integrate<qd::TanhSinh>([](double x) { return std::exp(-x * x); }, -inf, inf); // sqrt(pi)
 * @endcode
 *
 * The ends are oriented: b < a gives the negated integral, and a == b gives 0 without calling fn. An infinite end is
 * for TanhSinh only (a method opts in, see estimating_method); every other method returns error::non_finite. An
 * interval narrow against its position (|a + b| far above |b - a|) rounds each point at which fn is evaluated by about
 * eps |a + b| / 2, so the result errs by about that times the integral of |f'|: fn's conditioning, which no method
 * sees.
 */

#pragma once

#include <numerixx/math.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <expected>
#include <format>
#include <functional>
#include <limits>
#include <optional>
#include <ranges>
#include <span>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

/**
 * @brief Numerical integration (quadrature): the integral of a real function of one variable over a finite or infinite
 * interval, computed from the function's values at chosen points.
 */
namespace nxx::quadrature {

    /**
     * @brief Why integrate() returned no result, so that the caller can decide what to try next.
     *
     * non_finite means that a number was not finite where the method cannot accept one: in the input or in fn's
     * values. A different tolerance does not help; the interval, the method or the way fn is written must change.
     * not_converged means that the method ran out of work before its error estimate met the request: the integral may
     * need another method or a substitution, or it may not exist.
     *
     * More values may be added, and a listed cause may move to a new value: a switch on it should have a default. It
     * is formattable: std::format("{}", e) writes the enumerator's name, such as `not_converged`.
     */
    enum class error {
        /**
         * A number was not finite where the method cannot accept one, so no estimate could be formed.
         *
         * The causes:
         * - An end is NaN, or infinite for a method that does not take infinite ends (every method but TanhSinh).
         * - A tolerance is negative or not finite.
         * - fn returned a value that is not finite where the method cannot step around it. GaussKronrod: on a panel
         *   of its initial mesh that no short bisection avoids, such as a region where fn is NaN. TanhSinh: anywhere
         *   but within eps times the half-width of an end and beyond the points kept towards it, so that 1 / sqrt x
         *   on [0, 1] is ok, 1/x not_converged and sqrt(x - 0.5) non_finite; towards an infinite end, anywhere its
         *   tail rule does not apply (beyond the kept points, at a distance above 1/eps, or after a term already
         *   negligible), so that x^50 e^-x on [0, inf) is ok and x^51 e^-x, whose last term before the overflow is
         *   not yet negligible, non_finite.
         * - A sum overflowed.
         * - A method of one's own returned a result that is not finite, or reported its own failure so.
         */
        non_finite,

        /**
         * The method ran out of work before its error estimate met the request; Failure::best holds the estimate it
         * reached, where it has one.
         *
         * The ways a method runs out:
         * - Its cap (for TanhSinh, 10 levels).
         * - For GaussKronrod, stagnation: its bisections stop lowering the estimate's excess over the floor (fn
         *   noisier than the request, or a divergent integral).
         * - For GaussKronrod, every panel final with their errors above the tolerance (a strong endpoint singularity,
         *   or a piece where fn is not finite that a short bisection did not avoid).
         *
         * Romberg and GaussKronrod always give a best; a method of one's own may have none. TanhSinh has none where
         * the interval has no point of T strictly inside, or where its bound beyond an end is unbounded: the terms grow
         * towards the end (1/x on [0, 1] or [1, inf), a divergent integral) or converge too slowly for T's range
         * (x^-0.999 on [0, 1]; substitute x = u^k with fn written in u, since fn(u^k) underflows and gains nothing, or
         * use a type with a wider exponent range; the default is no better). A best from a divergent integral is a
         * number with no meaning.
         */
        not_converged
    };

    /**
     * @brief What a method with an error estimate returns on success: the integral, and an estimate of how far it is
     * from the exact value, to compare with the accuracy you need.
     *
     * abs_error estimates |value - exact|; it is not a bound. No method that samples fn at finitely many points can
     * bound the error: a narrow peak between the points is missed, and an fn that vanishes at every point sampled
     * gives 0 with the floor as its estimate, whatever its integral. For the library's methods and a != b, abs_error is
     * never below the rounding floor (see Criteria), so an integral that is exactly 0 reports the floor, not 0; a == b
     * gives {0, 0}. Read the members by name: members may be added in a later version, which would break a structured
     * binding.
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Integral
    {
        /** @brief The integral. */
        T value;

        /** @brief An estimate of |value - exact|, not a bound; for the library's methods never below the floor. */
        T abs_error;

        /** @brief Compares both members. */
        friend constexpr bool operator==(Integral const &, Integral const &) = default;
    };

    /**
     * @brief What a method with an error estimate returns instead of an Integral: why it failed, and, where it ran out
     * of work, the best estimate it reached, which may still serve a caller who can accept its error.
     *
     * Formattable: std::format("{}", f) writes the reason's name, such as `not_converged`.
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Failure
    {
        /** @brief Why there is no result. */
        error reason {};

        /**
         * @brief For error::not_converged the last estimate and its error estimate, if the method has one; empty for
         * non_finite. Meaningless where the integral diverges.
         */
        std::optional<Integral<T>> best {};

        /** @brief Compares both members. */
        friend constexpr bool operator==(Failure const &, Failure const &) = default;
    };

    /**
     * @brief The accuracy you ask of a method with an error estimate: it stops once its estimate is at most
     * max(abs_tol, rel_tol |value|), or at the least error it can certify.
     *
     * The defaults ask for about half the type's digits (rel_tol = sqrt(eps), 1.5e-8 in double), not full precision:
     * no estimate from samples certifies the last digits. Set rel_tol for more digits; a caller who chose a 50-digit
     * type for its digits must. Set abs_tol where the integral may be near 0 and a distance from 0 is close enough: it
     * saves evaluations. nxx::roots::Criteria asks for full precision by default; under using-directives for both
     * namespaces this name and `error` collide with nxx::roots', so qualify them.
     *
     * Below a rounding floor no request can be met, since the rounding of fn's values and of the sum then dominates the
     * method's error. The floor of the library's methods is 50 max(eps L1, denorm_min), L1 the integral of |fn| by the
     * same method on the same evaluations (for GaussKronrod the sum of every panel's): about 1e-14 relative in double
     * for an fn of one sign. GaussKronrod and TanhSinh accept at twice the floor, Romberg at the floor, and each ends
     * an integral that cancels to 0 there, so abs_tol = 0 works. TanhSinh's floor includes the slivers between its last
     * points and the ends, each up to about an ulp of its end, so on an interval a few million ulps wide or less its
     * abs_error can exceed the request. A method of one's own certifies its own minimum (Step::attainable), or none.
     *
     * Every estimate assumes fn accurate to about eps: noise in fn above the floor is invisible to the samples before
     * it shows, so a noisy fn can stop early with an estimate below its true error, or run to the cap.
     *
     * @code
     * auto const r = nxx::quadrature::integrate<nxx::quadrature::Romberg>(fn, 0.0, 1.0, {.rel_tol = 1e-12});
     * @endcode
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Criteria
    {
        /** @brief The absolute tolerance, for integrals near 0; 0 asks for the relative tolerance alone. */
        T abs_tol = T(0);

        /** @brief The relative tolerance: sqrt(eps) by default, 1.5e-8 in double; at most 1. */
        T rel_tol = nxx::math::root_epsilon<2, T>();

        /** @brief Compares both members. */
        friend constexpr bool operator==(Criteria const &, Criteria const &) = default;
    };

    /**
     * @brief What an iterative method reports after its start and after each step, so that integrate() can decide
     * whether to stop and an observer can watch the method converge.
     *
     * A method knows its own limits (the minimum work before an estimate is trusted, the least error it can certify,
     * when it cannot go on); integrate() owns the decision against the caller's criteria. Step carries the one to the
     * other, so that no method sees the criteria and every method is stopped by the same test. The defaults are the
     * neutral values, so a method sets only the limits it has:
     * `Step<T>{.current = Integral<T>{value, estimate}, .attainable = floor, .ready = level >= 4}`.
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Step
    {
        /**
         * @brief The value and its estimate now, for the ends in order (a < b), as integrate() would return them; empty
         * where the method has no usable estimate yet or at all.
         */
        std::optional<Integral<T>> current{};

        /**
         * @brief The least error the method can certify, so that a request below it still ends: integrate() accepts an
         * estimate at most max(abs_tol, rel_tol |value|, attainable).
         *
         * Romberg's is its rounding floor F, Gauss-Kronrod's and TanhSinh's 2 F (TanhSinh's F with its slivers).
         * Finite and not negative; 0, the default, certifies nothing beyond the request.
         */
        T attainable = T(0);

        /**
         * @brief Whether the method has done the minimum work for its estimate to be trusted (Romberg, TanhSinh: level
         * 4), so that an early estimate that is small by chance does not stop it; true by default.
         */
        bool ready = true;

        /**
         * @brief Whether the method cannot go on (its cap, stagnation, every panel final, or for TanhSinh no point
         * strictly inside), so that integrate() returns not_converged instead of asking for another step; false by
         * default.
         */
        bool exhausted = false;

        /**
         * @brief Whether a result that meets the request may be returned as a success; false where the method knows its
         * estimate is unreliable (Gauss-Kronrod after it gave up on a piece of the interval). True by default.
         */
        bool trusted = true;
    };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief The types whose rules and tables are computed at compile time, once per program, so that the common
         * types pay no run-time setup; every other type builds them at first use and keeps them.
         *
         * float, double and long double, named, since std::is_floating_point_v also admits std::float16_t and
         * std::bfloat16_t, whose tables have not been measured.
         *
         * @tparam T The real type.
         */
        template<typename T>
        concept tabled = std::same_as<T, float> || std::same_as<T, double> || std::same_as<T, long double>;

        /**
         * @brief The name of an error's enumerator: the text both formatters write, kept in one place.
         * @param e The error.
         * @return The name, such as `non_finite`.
         */
        constexpr std::string_view name(error e)
        {
            // No default: a new enumerator makes the compilers warn here (-Wswitch) until its name is added.
            switch (e) {
            case error::non_finite:
                return "non_finite";
            case error::not_converged:
                return "not_converged";
            }
            return "unknown";
        }

        /**
         * @brief The larger of @p a and @p b, by value: the maxima of the floor and of the stop test.
         * @param a A value.
         * @param b A value.
         * @return b if a < b, else a.
         */
        template<real T>
        constexpr T larger(T const &a, T const &b)
        { return a < b ? b : a; }

        /**
         * @brief Whether |@p x| > max / 2: the threshold of the constant-evaluation guards (design.md, "Constant
         * evaluation").
         *
         * An overflow is not a constant expression, so in constant evaluation the methods must refuse a sum before it
         * could overflow rather than detect the overflow after it. A sum or difference of two values at most max / 2
         * cannot overflow; a value above max / 2 is refused instead, and the method returns error::non_finite.
         *
         * @param x A finite value.
         * @return |x| > max / 2.
         */
        template<real T>
        constexpr bool large(T const &x)
        { return nxx::math::abs(x) > T(T((std::numeric_limits<T>::max)()) / T(2)); }

        /**
         * @brief Whether the product @p a @p b stays within max / 2, decided without forming an overflow: the
         * constant-evaluation guard for a product, as large() is for a sum.
         *
         * For |a| <= 1 the product, at most |b|, is formed and compared; else |b| <= (max / 2) / |a|, which bounds the
         * product by (max / 2)(1 + eps) < max.
         *
         * @param a A finite value.
         * @param b A finite value.
         * @return Whether a b is safe for the constant-evaluation guard.
         */
        template<real T>
        constexpr bool small_product(T const &a, T const &b)
        {
            T const m = nxx::math::abs(a);
            T const half_max = T(T((std::numeric_limits<T>::max)()) / T(2));
            if (m <= T(1)) return T(m * nxx::math::abs(b)) <= half_max;
            return nxx::math::abs(b) <= T(half_max / m);
        }

        /**
         * @brief Whether @p a + @p b may be formed, so that one call site serves run time and constant evaluation.
         *
         * At run time always: an overflow is found by a later check of the sum, which costs nothing per addition. In
         * constant evaluation, where an overflow is not a constant expression, only if neither exceeds max / 2.
         *
         * @param a A finite value.
         * @param b A finite value.
         * @return Whether the sum is safe to form.
         */
        template<real T>
        constexpr bool fits(T const &a, T const &b)
        {
            if consteval {
                return !detail::large(a) && !detail::large(b);
            }
            return true;
        }

        /**
         * @brief An array of M elements, each T(0), for the rules' and states' arrays: nxx::real does not require T to
         * be default constructible, so std::array<T, M>{} may not compile.
         * @tparam M The size.
         * @return The array.
         */
        template<real T, std::size_t M>
        constexpr std::array<T, M> zeros()
        {
            return []<std::size_t... I>(std::index_sequence<I...>) {
                return std::array<T, M>{ ((void)I, T(0))... };
            }(std::make_index_sequence<M>{});
        }

        /**
         * @brief fn at @p x, converted to T: the one way the methods call fn, so that fn may return any type
         * convertible to T (a double for a float T) and the conversion happens in one place.
         * @param fn The function.
         * @param x The point, a T const lvalue as the constraint promises.
         * @return fn(x) as a T.
         */
        template<real T, typename F>
        constexpr T evaluate(F &fn, T const &x)
        { return static_cast<T>(std::invoke(fn, x)); }

        /**
         * @brief A sum compensated by Neumaier's variant of Kahan's method, so that a sum of many terms errs by about
         * one rounding of the result rather than by the number of terms times eps (design.md, "Sums").
         *
         * The rounding error of each addition is accumulated separately and added at the end. The methods sum hundreds
         * to thousands of terms, and their rounding floor, 50 eps L1, assumes the sum's own error is far below it.
         *
         * A non-finite term makes the sum, and value(), non-finite (the correction becomes NaN or infinite), as does
         * an overflow of the running sum, so the caller checks value() once.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct Sum
        {
            /** @brief The running sum. */
            T sum = T(0);

            /** @brief The accumulated rounding errors of the additions. */
            T correction = T(0);

            /**
             * @brief Adds @p x.
             * @param x The term.
             */
            constexpr void add(T const &x)
            {
                T const t = T(sum + x);
                // The rounding error of sum + x is exact as (sum - t) + x when |sum| >= |x|, else as (x - t) + sum.
                if (nxx::math::abs(sum) >= nxx::math::abs(x)) correction = T(correction + T(T(sum - t) + x));
                else correction = T(correction + T(T(x - t) + sum));
                sum = t;
            }

            /**
             * @brief The compensated sum.
             * @return sum + correction.
             */
            constexpr T value() const
            { return T(sum + correction); }
        };

        /**
         * @brief The compensated value of @p s, or nothing where constant evaluation could not form it, so that a sum
         * near max / 2 becomes error::non_finite there instead of a failed compilation.
         *
         * sum + correction overflows only where |sum| > max / 2 and the correction is nonzero (the correction of terms
         * within max / 2 stays far below max / 2), and an overflow is not a constant expression. At run time always the
         * value.
         *
         * @param s The sum.
         * @return s.value(), or std::nullopt.
         */
        template<real T>
        constexpr std::optional<T> compensated(Sum<T> const &s)
        {
            if consteval {
                if (detail::large(s.sum) && !(s.correction == T(0))) return std::nullopt;
            }
            return s.value();
        }

        /**
         * @brief A finite interval [a, b], a != b, in the form the methods need: where to place a point, and how to
         * scale a term, without ever overflowing or losing a tiny value (design.md, "The interval" and "The scale of
         * each term").
         *
         * The obvious forms fail at the extremes: (a + b) / 2 and b - a overflow on [-max, max], and a scale
         * (b - a) w / 2 formed first can be subnormal on a tiny interval or make a tiny value subnormal on a wide one.
         * Interval keeps the centre and half-width, which never overflow, b - a where it is finite, and clamps every
         * point into [lo, hi]. interval_of() builds it once per integral, and every finite-interval method uses it.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct Interval
        {
            /** @brief The centre, a / 2 + b / 2: the origin of the points. */
            T c;

            /** @brief The half-width, b / 2 - a / 2, signed; neither it nor c can overflow, on [-max, max] included. */
            T hw;

            /** @brief b - a, where finite; the scale of the terms. */
            T width;

            /** @brief Whether b - a is finite, which decides how a term is formed. */
            bool finite;

            /** @brief Whether |b - a| >= 1, which decides the order of a term's products. */
            bool wide;

            /** @brief min(a, b), the lower bound of every point. */
            T lo;

            /** @brief max(a, b), the upper bound of every point. */
            T hi;

            /**
             * @brief The point @p x, clamped into [lo, hi], so that fn is never called outside the interval.
             *
             * c + hw u can leave the interval where the halves round: on [denorm_min, 3 denorm_min], a / 2 = 0 and
             * b / 2 = 2 denorm_min, so c = hw = 2 denorm_min and the points reach 0 and 4 denorm_min. Clamped, a point
             * may coincide with an end, as it may anyway on an interval a few ulps wide.
             *
             * @param x A point.
             * @return x, or the nearer end.
             */
            constexpr T point(T const &x) const
            { return x < lo ? lo : (hi < x ? hi : x); }

            /**
             * @brief The term hw w y, formed before it is summed, in the order that keeps every intermediate normal
             * unless the term itself is subnormal or overflows.
             *
             * Forming the term before the sum matters: a sum of values formed first could overflow. The order of the
             * products is chosen by the interval. Where b - a is finite, the term is (b - a) (w / 2) y: with
             * |b - a| >= 1 the scale (b - a) (w / 2) is formed first, which lies between w / 2 and max (w / 2 <= 1), so
             * that a tiny y on a wide interval is not made subnormal first; with |b - a| < 1 the value is multiplied
             * first, so that a subnormal scale (b - a) w / 2 is never formed. Where b - a overflows, |hw| > max / 2:
             * (hw y) w for |y| <= 1, where hw y cannot overflow and is normal unless y is 0, else hw (w y), where
             * |w y| >= w is normal (hw w would overflow for w > 1). w / 2 is exact in radix 2.
             *
             * @param w The weight relative to the half-width: a Gauss weight on [-1, 1], or Romberg's 2^(1-k).
             * @param y fn's value.
             * @return The term.
             */
            constexpr T term(T const &w, T const &y) const
            {
                if (!finite) {
                    if (nxx::math::abs(y) <= T(1)) return T(T(hw * y) * w);
                    return T(hw * T(w * y));
                }
                T const half = T(w / T(2));
                if (wide) return T(T(width * half) * y);
                return T(width * T(half * y));
            }

            /**
             * @brief The term, or nothing where constant evaluation could not form it within max / 2, so that such a
             * term becomes error::non_finite there instead of a failed compilation.
             *
             * In constant evaluation an overflow, or arithmetic on a NaN, is not a constant expression (GCC and Clang
             * reject them at different operations), so a value of fn that is not finite, or a product of term() that
             * could exceed max / 2, gives nothing, and the method error::non_finite: conservative, where run time
             * would give the value. At run time the term is always formed, and the caller checks its sum once.
             *
             * @param w The weight relative to the half-width.
             * @param y fn's value.
             * @return The term, or std::nullopt.
             */
            constexpr std::optional<T> checked_term(T const &w, T const &y) const
            {
                if consteval {
                    if (!nxx::math::isfinite(y)) return std::nullopt;
                    // The same four orders as term(), each product checked before it is formed.
                    bool fits = false;
                    T const half = T(w / T(2));
                    if (!finite && nxx::math::abs(y) <= T(1))
                        fits = detail::small_product(hw, y) && detail::small_product(T(hw * y), w);
                    else if (!finite)
                        fits = detail::small_product(w, y) && detail::small_product(hw, T(w * y));
                    else if (wide)
                        fits = detail::small_product(width, half) && detail::small_product(T(width * half), y);
                    else
                        fits = detail::small_product(half, y) && detail::small_product(width, T(half * y));
                    if (!fits) return std::nullopt;
                }
                return term(w, y);
            }
        };

        /**
         * @brief The midpoint a / 2 + b / 2, the form that cannot overflow: the centre of an interval, and the point at
         * which Gauss-Kronrod bisects a panel. (a + b) / 2 overflows where a + b does.
         * @param a An end.
         * @param b The other end.
         * @return The midpoint.
         */
        template<real T>
        constexpr T midpoint(T const &a, T const &b)
        { return T(T(a / T(2)) + T(b / T(2))); }

        /**
         * @brief The rounding floor: the least error estimate a method trusts, below which the rounding of fn's values
         * and of the sum, not the method, dominates the error (design.md, "The floor").
         *
         * 50 max(eps L1, denorm_min), with QUADPACK's constant 50. A request below the floor cannot be met, so the
         * methods accept an estimate at the floor (or twice it) instead of running to their cap; an integral that
         * cancels to 0 then ends too, without an abs_tol.
         *
         * @param l1 The integral of |fn| by the method on the same evaluations.
         * @return The floor.
         */
        template<real T>
        constexpr T rounding_floor(T const &l1)
        {
            T const eps = T(std::numeric_limits<T>::epsilon());
            T const tiny = T(std::numeric_limits<T>::denorm_min());
            return T(T(50) * detail::larger(T(eps * l1), tiny));
        }

        /**
         * @brief The largest error estimate at which an integration may stop: the caller's request, or the error the
         * method can certify where the request is below it. integrate()'s stop test compares the estimate with it.
         * @param abs_tol The absolute tolerance.
         * @param rel_tol The relative tolerance, at most 1.
         * @param value The integral.
         * @param floor The least error the method can certify (Step::attainable).
         * @return max(abs_tol, rel_tol |value|, floor).
         */
        template<real T>
        constexpr T wanted(T const &abs_tol, T const &rel_tol, T const &value, T const &floor)
        {
            T const relative = T(rel_tol * nxx::math::abs(value));
            return detail::larger(detail::larger(abs_tol, relative), floor);
        }

        /**
         * @brief The Interval of [a, b], built once per integral by the finite-interval methods.
         * @param a The lower end, finite.
         * @param b The upper end, finite, not a.
         * @return Its centre, half-width and scale.
         */
        template<real T>
        constexpr Interval<T> interval_of(T const &a, T const &b)
        {
            // Whether b - a overflows, decided before it is formed (constant evaluation rejects the overflow): ends of
            // one sign never overflow; ends of opposite signs do exactly when |a| / 2 + |b| / 2 > max / 2. In radix 2
            // the halves of such ends are exact, and the rounded half-sum exceeds max / 2 exactly when the rounded
            // |a| + |b| exceeds max: the same rounding, scaled by 2.
            bool const opposite = (a < T(0)) != (b < T(0));
            T const half_sum = T(T(nxx::math::abs(a) / T(2)) + T(nxx::math::abs(b) / T(2)));
            bool const finite = !opposite || !(half_sum > T(T((std::numeric_limits<T>::max)()) / T(2)));
            T const width = finite ? T(b - a) : T(0);
            return Interval<T>{ detail::midpoint(a, b),
                                T(T(b / T(2)) - T(a / T(2))),
                                width,
                                finite,
                                finite && nxx::math::abs(width) >= T(1),
                                a < b ? a : b,
                                a < b ? b : a };
        }

        /**
         * @brief An upper bound on the bits of T's significand, which sizes the iteration caps that must scale with
         * T's precision: the bisections of the node generators, Gauss-Kronrod's cap (64 per bit) and exp()'s reach.
         *
         * digits times the bits of a digit of T's radix, rounded up: exact in radix 2, 4 bits per digit in radix 10
         * (cpp_dec_float_50). It needs only std::numeric_limits, which every nxx::real provides, and no logarithm.
         *
         * @tparam T The real type.
         * @return The bits.
         */
        template<real T>
        constexpr int significand_bits()
        {
            int bits_per_digit = 0;
            long long reach = 1;
            while (reach < std::numeric_limits<T>::radix) {
                reach *= 2;
                ++bits_per_digit;
            }
            return std::numeric_limits<T>::digits * bits_per_digit;
        }

        /**
         * @brief A polynomial's value and derivative at a point, which Newton's method needs together: what
         * bracketed_zero() asks of P_N for the Gauss nodes and of E_8 for the Kronrod nodes, and what the weights'
         * formulas read.
         * @tparam T The real type.
         */
        template<real T>
        struct Legendre
        {
            /** @brief The value, P_N(x) or E_8(x). */
            T p;

            /** @brief The derivative, P_N'(x) or E_8'(x). */
            T dp;
        };

        /**
         * @brief One step of the Legendre recurrence, shared by every evaluation of a Legendre polynomial or series,
         * so that the Gauss nodes, the Sturm count and the Kronrod polynomial all round the same way.
         *
         * P_(k+1)(x) from P_k(x) and P_(k-1)(x): (k + 1) P_(k+1) = (2k + 1) x P_k - k P_(k-1). The recurrence is the
         * stable way to evaluate P_n on [-1, 1], and it needs only + - * / in T (design.md, "GaussLegendre<N> (D)").
         *
         * @param k The degree of p1, at least 1.
         * @param x The point.
         * @param p1 P_k(x).
         * @param p0 P_(k-1)(x).
         * @return P_(k+1)(x).
         */
        template<real T>
        constexpr T next_legendre(int k, T const &x, T const &p1, T const &p0)
        { return T(T(T(T(T(2 * k + 1) * x) * p1) - T(T(k) * p0)) / T(k + 1)); }

        /**
         * @brief One step of the differentiated recurrence, for the derivative that Newton's method and the weights
         * need: (k + 1) P'_(k+1) = (2k + 1) (P_k + x P'_k) - k P'_(k-1).
         * @param k The degree of p1, at least 1.
         * @param x The point.
         * @param p1 P_k(x).
         * @param d1 P'_k(x).
         * @param d0 P'_(k-1)(x).
         * @return P'_(k+1)(x).
         */
        template<real T>
        constexpr T next_derivative(int k, T const &x, T const &p1, T const &d1, T const &d0)
        { return T(T(T(T(2 * k + 1) * T(p1 + T(x * d1))) - T(T(k) * d0)) / T(k + 1)); }

        /**
         * @brief P_n(x) and P_n'(x), the polynomial whose zeros are the Gauss-Legendre nodes: for Newton's method on
         * them, and for the weights 2 / ((1 - x)(1 + x) P_n'(x)^2).
         *
         * By the three-term recurrence and its derivative, from P_0 = 1 and P_1 = x.
         *
         * @param n The degree, at least 1.
         * @param x The point, in [-1, 1].
         * @return P_n(x) and P_n'(x).
         */
        template<real T>
        constexpr Legendre<T> legendre(int n, T const &x)
        {
            // |P_k| <= 1 on [-1, 1], so nothing overflows.
            T p0 = T(1);
            T p1 = x;
            T d0 = T(0);
            T d1 = T(1);
            for (int k = 1; k < n; ++k) {
                T const p2 = detail::next_legendre(k, x, p1, p0);
                T const d2 = detail::next_derivative(k, x, p1, d1, d0);
                p0 = p1;
                p1 = p2;
                d0 = d1;
                d1 = d2;
            }
            return Legendre<T>{ p1, d1 };
        }

        /**
         * @brief The number of zeros of P_n strictly above @p x, so that the node generator can isolate each zero by
         * bisection, by counting alone, before Newton's method refines it (design.md, "GaussLegendre<N> (D)").
         *
         * The usual starting guesses for the nodes need pi and cos, which an nxx::real need not have; a count of zeros
         * above a point needs only the recurrence, and bisection on it brackets exactly one zero.
         *
         * The count is the number of sign changes of P_0(x), ..., P_n(x). The sequence is a Sturm sequence (the
         * recurrence of orthogonal polynomials): its sign changes fall by one as x passes a zero of P_n upward, from
         * n below -1 to 0 above 1. A zero term is skipped and the last nonzero sign kept: where an intermediate P_k
         * vanishes its neighbours have opposite signs, and where P_n vanishes at x, x itself is not counted.
         *
         * @param n The degree, at least 1.
         * @param x The point.
         * @return The count.
         */
        template<real T>
        constexpr int zeros_above(int n, T const &x)
        {
            int changes = 0;
            bool negative = false; // the sign of P_0 = 1
            auto const count = [&changes, &negative](T const &p) {
                if (p == T(0)) return;
                bool const below = p < T(0);
                if (below != negative) ++changes;
                negative = below;
            };
            T p0 = T(1);
            T p1 = x;
            count(p1);
            for (int k = 1; k < n; ++k) {
                T const p2 = detail::next_legendre(k, x, p1, p0);
                p0 = p1;
                p1 = p2;
                count(p1);
            }
            return changes;
        }

        /**
         * @brief The one zero of a polynomial in (@p lower, @p upper), to T's precision: the step that turns an
         * isolated zero into a node, for P_n (the Gauss nodes) and for the Stieltjes polynomial E_8 (the Kronrod
         * nodes).
         *
         * Newton's method converges fast but can leave the bracket; bisection cannot, so a Newton step is taken only
         * where it lands strictly inside, and the bracket otherwise halves (design.md, "GaussLegendre<N> (D)"). The
         * first iterate is the midpoint. It stops when a Newton step is at most 2 eps |x|, or when a step below
         * sqrt(eps) |x| stops shrinking (rounding); the second stop never fires on the library's rules but guards a
         * type whose rounding stalls Newton above 2 eps |x|.
         *
         * The orientation of the bracket is given, not computed: a computed value of exactly 0 at an end would turn it
         * around. The result lies strictly inside the bracket: every iterate does.
         *
         * @param polynomial Its value and derivative at a point, as a Legendre<T>.
         * @param lower The lower end of a bracket that holds one zero.
         * @param upper The upper end.
         * @param negative_below Whether the polynomial is negative between lower and the zero.
         * @return The zero.
         */
        template<real T, typename P>
        constexpr T bracketed_zero(P const &polynomial, T const &lower, T const &upper, bool negative_below)
        {
            T const eps = T(std::numeric_limits<T>::epsilon());
            T const sqrt_eps = nxx::math::root_epsilon<2, T>();
            // A safeguard only: each bisection halves the bracket and Newton converges quadratically inside it, so the
            // loop ends long before at T's resolution.
            int const cap = 2 * (detail::significand_bits<T>() + 8);
            T lo = lower;
            T hi = upper;
            T x = T(T(lo + hi) / T(2));
            if (!(lo < x && x < hi)) return x; // no point of T inside: the bracket is at T's resolution already
            T previous = T(0);
            bool previous_newton = false;
            for (int i = 0; i < cap; ++i) {
                Legendre<T> const e = polynomial(x);
                if (e.p == T(0)) break;
                if ((e.p < T(0)) == negative_below) lo = x;
                else hi = x;
                T next = x;
                bool newton = false;
                // |p| < |dp| (hi - lo) means the step is shorter than the bracket: it cannot overflow, and dp != 0.
                if (nxx::math::abs(e.p) < T(nxx::math::abs(e.dp) * T(hi - lo))) {
                    next = T(x - T(e.p / e.dp));
                    // A step that rounds back to x, now an end of the bracket, is convergence (a step of 0), not a
                    // reason to bisect away from the zero.
                    newton = next == x || (lo < next && next < hi);
                }
                if (!newton) {
                    next = T(T(lo + hi) / T(2));
                    if (!(lo < next && next < hi)) break; // the bracket is at T's resolution
                }
                T const step = nxx::math::abs(T(next - x));
                x = next;
                if (newton && step <= T(T(T(2) * eps) * nxx::math::abs(x))) break;
                // A step that no longer shrinks is rounding, but only once it is small: two growing steps near an
                // extremum of the polynomial are not convergence.
                if (newton && previous_newton && step >= previous && step <= T(sqrt_eps * nxx::math::abs(x))) break;
                previous = step;
                previous_newton = newton;
            }
            assert(lower < x && x < upper && "a zero must lie strictly inside its bracket");
            return x;
        }

        /**
         * @brief The N-point Gauss-Legendre rule as the methods store and apply it: the non-negative half, since the
         * rule is symmetric about 0, with the nodes ascending and their weights.
         *
         * apply() evaluates fn at c + hw x and c - hw x for each stored x, and once at c for x = 0 (odd N).
         *
         * @tparam T The real type.
         * @tparam N The number of points.
         */
        template<real T, int N>
        struct Rule
        {
            /** @brief The number of non-negative nodes, (N + 1) / 2. */
            static constexpr std::size_t size = static_cast<std::size_t>((N + 1) / 2);

            /** @brief The non-negative nodes, ascending; x[0] = 0 for odd N. */
            std::array<T, size> x;

            /** @brief The weights of the nodes. */
            std::array<T, size> w;
        };

        /**
         * @brief The N-point Gauss-Legendre rule computed in T, so that any nxx::real gets nodes and weights to its own
         * precision, from + - * / alone (no pi, no cos), at compile time or at run time (design.md,
         * "GaussLegendre<N> (D)").
         *
         * The positive zeros of P_N are isolated one by one, from the largest, by bisection on zeros_above(), then
         * found by bracketed_zero(); the negative nodes are their mirror images, and 0 is a node for odd N. The weights
         * are 2 / ((1 - x)(1 + x) P_N'(x)^2); (1 - x)(1 + x) rather than 1 - x^2, whose subtraction cancels near the
         * outer nodes, while 1 - x is exact there. float, double and long double read the result from a compile-time
         * table (gauss_legendre_table); other types generate it at first use and keep it (with_generated()).
         *
         * @tparam N The number of points.
         * @tparam T The real type.
         * @return The non-negative half of the rule.
         */
        template<int N, real T>
        constexpr Rule<T, N> gauss_legendre_rule()
        {
            constexpr int m = N / 2;   // the positive zeros
            constexpr int odd = N % 2; // 0 is a zero of P_N for odd N
            // The isolation's cap: a type with few digits could miscount near a zero forever; after it the bracket is
            // used as it is.
            int const cap = detail::significand_bits<T>() + 8;
            Rule<T, N> rule{ detail::zeros<T, Rule<T, N>::size>(), detail::zeros<T, Rule<T, N>::size>() };
            T upper = T(1); // above the k-th largest zero x_k and below x_(k-1)
            for (int k = 1; k <= m; ++k) {
                T lo = T(0);
                T hi = upper;
                bool isolated = false;
                for (int step = 0; step < cap; ++step) {
                    T const mid = T(T(lo + hi) / T(2));
                    int const above = detail::zeros_above(N, mid);
                    if (above >= k) {
                        lo = mid;
                        isolated = above == k;
                        if (isolated) break; // (lo, hi) holds x_k and no other zero
                    }
                    else hi = mid;
                }
                // P_N > 0 above its largest zero x_1, and its sign between x_(k+1) and x_k is (-1)^k: lo lies there.
                auto const p = [](T const &v) { return detail::legendre(N, v); };
                T const x = detail::bracketed_zero(p, lo, hi, k % 2 == 1);
                rule.x[static_cast<std::size_t>(odd + m - k)] = x;
                // lo's count is exactly k, so x_(k+1) <= lo < x_k: an end found by counting, where the computed zero
                // could lie on either side of x_k. The computed zero only where the cap ended the isolation.
                upper = isolated ? lo : x;
            }
            if constexpr (odd == 1) rule.x[0] = T(0);
            for (std::size_t j = 0; j < Rule<T, N>::size; ++j) {
                T const x = rule.x[j];
                Legendre<T> const e = detail::legendre(N, x);
                rule.w[j] = T(T(2) / T(T(T(T(1) - x) * T(T(1) + x)) * T(e.dp * e.dp)));
            }
            return rule;
        }

        /**
         * @brief The N-point Gauss-Legendre rule for float, double and long double (tabled), computed at compile time
         * once per program, so that integrate<GaussLegendre<N>> costs N evaluations of fn and nothing more at run time.
         * @tparam N The number of points.
         * @tparam T float, double or long double.
         */
        template<int N, typename T>
        inline constexpr Rule<T, N> gauss_legendre_table = detail::gauss_legendre_rule<N, T>();

        /**
         * @brief The rule @p generate makes, built at its first use and kept, so that a type without a compile-time
         * table, such as a Boost.Multiprecision number, pays for its nodes once per program rather than on every call
         * (design.md, "The node cache").
         *
         * One function-local static per generator, so per rule and T, whatever fn is (initialised once,
         * thread-safely, by the language; never changed). Keyed on the generator alone, not on a function templated on
         * fn's type, where each lambda would build its own copy. A type that changed its precision at run time would
         * see the nodes of its first; no nxx::real type does. Not for constant evaluation, where a static is not
         * allowed.
         * @tparam generate The generator, such as gauss_legendre_rule<N, T>.
         * @return The rule.
         */
        template<auto generate>
        auto const &cached()
        {
            static auto const rule = generate();
            return rule;
        }

        /**
         * @brief @p use applied to the rule @p generate makes, for a type without a compile-time table: one call site
         * for both contexts, since constant evaluation cannot read the cache.
         *
         * The generator itself in constant evaluation (each evaluation generates the rule anew), the cached rule at
         * run time.
         * @tparam generate The generator.
         * @param use A callable taking the rule.
         * @return What use returns.
         */
        template<auto generate, typename U>
        constexpr decltype(auto) with_generated(U &&use)
        {
            if consteval {
                return use(generate());
            }
            else {
                return use(detail::cached<generate>());
            }
        }

        /**
         * @brief @p use applied to the N-point Gauss-Legendre rule in T, so that GaussLegendre<N> need not know where
         * its rule comes from: the compile-time table for float, double and long double, with_generated() for every
         * other T.
         * @param use A callable taking the rule.
         * @return What use returns.
         */
        template<int N, real T, typename U>
        constexpr decltype(auto) with_gauss_legendre(U &&use)
        {
            if constexpr (tabled<T>) return use(detail::gauss_legendre_table<N, T>);
            else return detail::with_generated<&detail::gauss_legendre_rule<N, T>>(use);
        }

        /**
         * @brief The coefficients of the Stieltjes polynomial E_8 of the G7-K15 pair in the Legendre basis, whose zeros
         * are the 8 Kronrod nodes added to the 7 Gauss nodes (design.md, "Gauss-Kronrod (Phase 2)").
         *
         * E_8 = sum_k a_k P_k, a_8 = 1, orthogonal to P_7 P_j for every j < 8 (Kronrod 1965; Laurie 1997). Computing it
         * in T, rather than reading published nodes, gives the pair to T's precision for any nxx::real.
         *
         * The integrals of three Legendre polynomials have Adams' closed form,
         * int_{-1}^{1} P_a P_b P_c = 2 A(s - a) A(s - b) A(s - c) / ((2s + 1) A(s)), a + b + c = 2s, with
         * A(m) = (2m - 1)!! / (2m)!!, and 0 where a + b + c is odd or the triangle inequality fails. E_8 is even, so
         * the conditions with even j vanish identically (7 + k + j odd); the condition with odd j involves a_k for
         * k >= 7 - j only (the triangle inequality), and fixes a_(7-j) from the coefficients above it: a triangular
         * system, solved from the top, with no linear solve.
         *
         * @return a_0 ... a_8 (the odd ones 0).
         */
        template<real T>
        constexpr std::array<T, 9> stieltjes()
        {
            constexpr int n = 7;
            constexpr int smax = (3 * n + 1) / 2; // the largest s = (7 + k + j) / 2 asked for
            std::array<T, smax + 1> A = detail::zeros<T, smax + 1>();
            A[0] = T(1);
            for (int i = 1; i <= smax; ++i) {
                std::size_t const u = static_cast<std::size_t>(i);
                A[u] = T(T(A[u - 1] * T(2 * i - 1)) / T(2 * i));
            }
            auto const at = [&A](int i) -> T const & { return A[static_cast<std::size_t>(i)]; };
            auto const triple = [&at](int a, int b, int c) {
                int const half_sum = (a + b + c) / 2;
                return T(T(T(T(T(2) * at(half_sum - a)) * at(half_sum - b)) * at(half_sum - c))
                         / T(T(2 * half_sum + 1) * at(half_sum)));
            };
            std::array<T, 9> a = detail::zeros<T, 9>();
            a[8] = T(1);
            for (int j = 1; j <= n; j += 2) {
                int const k0 = n - j;
                T sum = T(0);
                for (int k = k0 + 2; k <= n + 1; k += 2)
                    sum = T(sum + T(a[static_cast<std::size_t>(k)] * triple(n, k, j)));
                a[static_cast<std::size_t>(k0)] = T(T(-sum) / triple(n, k0, j));
            }
            return a;
        }

        /**
         * @brief E_8(x) and E_8'(x) from E_8's Legendre coefficients, for Newton's method on the Kronrod nodes and for
         * the Kronrod weights.
         *
         * sum_k a_k P_k(x) and its derivative, with each P_k from the three-term recurrence as legendre() forms it.
         *
         * @param a The coefficients a_0 ... a_8.
         * @param x The point, in [-1, 1].
         * @return The value and the derivative.
         */
        template<real T>
        constexpr Legendre<T> legendre_series(std::array<T, 9> const &a, T const &x)
        {
            T p0 = T(1);
            T p1 = x;
            T d0 = T(0);
            T d1 = T(1);
            T sum = T(a[0] + T(a[1] * p1));
            T derivative = a[1];
            for (int k = 1; k < 8; ++k) {
                std::size_t const next = static_cast<std::size_t>(k + 1);
                T const p2 = detail::next_legendre(k, x, p1, p0);
                T const d2 = detail::next_derivative(k, x, p1, d1, d0);
                p0 = p1;
                p1 = p2;
                d0 = d1;
                d1 = d2;
                sum = T(sum + T(a[next] * p1));
                derivative = T(derivative + T(a[next] * d1));
            }
            return Legendre<T>{ sum, derivative };
        }

        /**
         * @brief The G7-K15 pair as Gauss-Kronrod's panels apply it: two rules on the same 15 points, so that their
         * difference estimates the error at no extra evaluation.
         *
         * The non-negative half, since both rules are symmetric about 0: the 8 nodes ascending, the 7-point Gauss
         * nodes at the even indices (x[0] = 0), the Kronrod nodes at the odd ones, with both rules' weights (the Gauss
         * weight 0 where a node is Kronrod's only).
         *
         * @tparam T The real type.
         */
        template<real T>
        struct Kronrod
        {
            /** @brief The number of non-negative nodes. */
            static constexpr std::size_t size = 8;

            /** @brief The non-negative nodes, ascending. */
            std::array<T, size> x;

            /** @brief The 15-point Kronrod weights: the rule whose value a panel returns. */
            std::array<T, size> wk;

            /** @brief The 7-point Gauss weights at the even indices, 0 at the odd ones: the rule compared with. */
            std::array<T, size> wg;
        };

        /**
         * @brief The G7-K15 pair computed in T, from + - * / alone, so that any nxx::real gets the pair to its own
         * precision (design.md, "Gauss-Kronrod (Phase 2)").
         *
         * The Gauss part is GaussLegendre<7>'s. The Kronrod nodes are the zeros of E_8, each between two Gauss nodes or
         * between the largest one and 1, found by bracketed_zero() with the orientation from the node's index (E_8 > 0
         * above its largest zero, sign (-1)^k between its k-th and (k+1)-th largest). The weights are the interpolatory
         * formulas, 2 / (8 P_7 E_8') at a Kronrod node and w_G + 2 / (8 P_7' E_8) at a Gauss node.
         *
         * @return The non-negative half of the pair.
         */
        template<real T>
        constexpr Kronrod<T> gauss_kronrod_rule()
        {
            Rule<T, 7> const gauss = detail::gauss_legendre_rule<7, T>();
            std::array<T, 9> const a = detail::stieltjes<T>();
            auto const e8 = [&a](T const &x) { return detail::legendre_series(a, x); };
            T const c = T(T(1) / T(4)); // 2 / (n + 1), n = 7
            Kronrod<T> rule{ detail::zeros<T, 8>(), detail::zeros<T, 8>(), detail::zeros<T, 8>() };
            for (std::size_t k = 1; k <= 4; ++k) {
                T const lo = gauss.x[4 - k];
                T const hi = k == 1 ? T(1) : gauss.x[5 - k];
                rule.x[2 * (4 - k) + 1] = detail::bracketed_zero(e8, lo, hi, k % 2 == 1);
            }
            for (std::size_t j = 0; j < 4; ++j) {
                rule.x[2 * j] = gauss.x[j];
                rule.wg[2 * j] = gauss.w[j];
            }
            for (std::size_t i = 0; i < Kronrod<T>::size; ++i) {
                Legendre<T> const e = e8(rule.x[i]);
                Legendre<T> const p = detail::legendre(7, rule.x[i]);
                if (i % 2 == 0) rule.wk[i] = T(rule.wg[i] + T(c / T(p.dp * e.p)));
                else rule.wk[i] = T(c / T(p.p * e.dp));
            }
            return rule;
        }

        /**
         * @brief The G7-K15 pair for float, double and long double, computed at compile time once per program, so
         * that the default method pays nothing for its rule at run time.
         * @tparam T float, double or long double.
         */
        template<typename T>
        inline constexpr Kronrod<T> gauss_kronrod_table = detail::gauss_kronrod_rule<T>();

        /**
         * @brief @p use applied to the G7-K15 pair in T, so that Gauss-Kronrod need not know where its rule comes from:
         * the compile-time table for float, double and long double, with_generated() for every other T, as
         * with_gauss_legendre().
         * @param use A callable taking the pair.
         * @return What use returns.
         */
        template<real T, typename U>
        constexpr decltype(auto) with_gauss_kronrod(U &&use)
        {
            if constexpr (tabled<T>) return use(detail::gauss_kronrod_table<T>);
            else return detail::with_generated<&detail::gauss_kronrod_rule<T>>(use);
        }

        /**
         * @brief GaussLegendre<N>'s whole computation: the N-point rule applied to fn on an interval, as one
         * compensated sum of N terms.
         *
         * sum over the nodes of hw w_i fn(c + hw x_i), each term formed by Interval::checked_term() in the order that
         * keeps it normal, each point clamped into the interval. The centre first for odd N, then each pair c + d,
         * c - d with the nodes ascending.
         *
         * @param rule The non-negative half of the rule.
         * @param fn The function.
         * @param interval The interval.
         * @return The integral, or error::non_finite if the sum is not finite.
         */
        template<real T, int N, typename F>
        constexpr std::expected<T, error> apply(Rule<T, N> const &rule, F &fn, Interval<T> const &interval)
        {
            Sum<T> sum;
            // Adds hw w fn(x); false only in constant evaluation, where a term or the sum could leave max / 2
            // (Interval::checked_term).
            auto const add = [&](T const &w, T const &x) {
                std::optional<T> const t = interval.checked_term(w, detail::evaluate(fn, interval.point(x)));
                if (!t) return false;
                if consteval {
                    if (detail::large(sum.sum)) return false;
                }
                sum.add(*t);
                return true;
            };
            std::size_t first = 0;
            if constexpr (N % 2 == 1) {
                if (!add(rule.w[0], interval.c)) return std::unexpected(error::non_finite);
                first = 1;
            }
            for (std::size_t j = first; j < Rule<T, N>::size; ++j) {
                T const d = T(interval.hw * rule.x[j]);
                if (!add(rule.w[j], T(interval.c + d)) || !add(rule.w[j], T(interval.c - d)))
                    return std::unexpected(error::non_finite);
            }
            // At run time every scale is finite and nonzero, so a non-finite value of fn, or an overflow, makes the sum
            // non-finite: one check, outside the loop.
            std::optional<T> const result = detail::compensated(sum);
            if (!result || !nxx::math::isfinite(*result)) return std::unexpected(error::non_finite);
            return *result;
        }

        /**
         * @brief What Romberg's method keeps between its levels, so that it can run as an iterative method whose loop
         * integrate() owns (PHASE3A_API.md): the state after level k, T_k's Richardson row and the row before it.
         *
         * Level k is the trapezoidal sum T_k on 2^k panels, which reuses T_(k-1)'s points and adds the 2^(k-1) new
         * midpoints; Richardson's extrapolation then turns T_0 ... T_k into R(k, k). Only two rows of the table are
         * needed, so the state is a fixed size (design.md, "Romberg (E)").
         *
         * @tparam T The real type.
         */
        template<real T>
        struct RombergState
        {
            /**
             * @brief The last level: 2^16 + 1 = 65,537 evaluations, then not_converged, so that an fn Romberg suits
             * badly (sqrt x, a kink) cannot run on. Fixed, not from T: smooth fn need at most 9 levels in double and 14
             * at 50 digits for 1e-45 (PHASE1_API.md, "The cap").
             */
            static constexpr int cap = 16;

            /**
             * @brief The first level at which the method may stop, 17 evaluations, so that a few samples that happen
             * to agree do not end it: with fewer, the samples of sin^2 on [0, 2 pi] all vanish (Numerixx 1.x returned
             * 0 for it after level 2).
             */
            static constexpr int min_level = 4;

            /** @brief The length of a Richardson row: R(k, 0..k) for k up to the cap. */
            static constexpr std::size_t row_size = cap + 1;

            /** @brief The interval [a, b], a < b: the points and the terms' scale. */
            Interval<T> interval;

            /**
             * @brief The weight of level k's new terms relative to the half-width: 1 for T_0 and level 1, then
             * 2^(1-k), the trapezoidal rule's step halving with each level.
             */
            T weight;

            /** @brief The trapezoidal sum of |terms| at level k, L1, from which the rounding floor is formed. */
            T l1;

            /** @brief The Richardson row of level k - 1, R(k-1, 0..k-1): what level k extrapolates against. */
            std::array<T, row_size> previous;

            /** @brief The Richardson row of level k, R(k, 0..k); R(k, k) is the method's value. */
            std::array<T, row_size> row;

            /** @brief The level. */
            int k;

            /** @brief The step: {R(k, k), max(estimate, floor)}, the floor attainable, ready from min_level. */
            Step<T> step;
        };

        /**
         * @brief Level k's new evaluations: fn at the midpoints of level k - 1's panels, the only points level k adds,
         * so that each level costs 2^(k-1) evaluations and reuses every earlier one.
         *
         * The midpoints are u_j = (2j - 1 - m) / m on [-1, 1], m = 2^(k-1), mapped to c + hw u_j. Each term is added
         * to sum (compensated, for T_k) and its magnitude to abs_sum (for L1), in order; in constant evaluation each
         * term and sum is kept within max / 2 (Interval::checked_term).
         *
         * @param fn The function.
         * @param s The state.
         * @param sum The compensated sum of the terms.
         * @param abs_sum The plain sum of |terms|: its terms are positive, so it errs by about m eps relative.
         * @return false where a term or a sum could not be formed: non_finite.
         */
        template<real T, typename F>
        constexpr bool add_midpoints(F &fn, RombergState<T> const &s, Sum<T> &sum, T &abs_sum)
        {
            int const m = 1 << (s.k - 1);
            return std::ranges::all_of(std::views::iota(1, m + 1), [&](int j) {
                T const u = T(T(2 * j - 1 - m) / T(m)); // exact: at most 16 bits
                T const x = s.interval.point(T(s.interval.c + T(s.interval.hw * u)));
                std::optional<T> const t = s.interval.checked_term(s.weight, detail::evaluate(fn, x));
                if (!t) return false;
                if consteval {
                    if (detail::large(sum.sum) || detail::large(abs_sum)) return false;
                }
                sum.add(*t);
                abs_sum = T(abs_sum + nxx::math::abs(*t));
                return true;
            });
        }

        /**
         * @brief Richardson's extrapolation of level k, which removes the trapezoidal rule's error terms in h^2, h^4,
         * ... one by one, so that R(k, k) converges far faster than T_k on smooth fn.
         *
         * R(k, i) = R(k, i-1) + (R(k, i-1) - R(k-1, i-1)) / (4^i - 1) for i = 1 ... k, from R(k, 0) = T_k. 4^i - 1 is
         * exact in radix 2 up to the cap for double, and in radix 10 (cpp_dec_float_50), and for float up to i = 12,
         * beyond which the correction is at rounding level anyway.
         *
         * @param s The state, with row[0] formed.
         * @return false where constant evaluation could not keep an entry within max / 2: non_finite.
         */
        template<real T>
        constexpr bool extrapolate(RombergState<T> &s)
        {
            T power = T(1);
            return std::ranges::all_of(std::views::iota(std::size_t{ 1 }, static_cast<std::size_t>(s.k) + 1),
                [&](std::size_t i) {
                    if consteval {
                        if (detail::large(s.row[i - 1]) || detail::large(s.previous[i - 1])) return false;
                    }
                    power = T(power * T(4));
                    s.row[i] = T(s.row[i - 1] + T(T(s.row[i - 1] - s.previous[i - 1]) / T(power - T(1))));
                    return true;
                });
        }

        /**
         * @brief One level of Romberg's method, shared by its start and its steps: T_k, its Richardson row, the
         * estimate and the floor, as Phase 1's loop formed them, reported in s.step for integrate()'s stop test.
         *
         * T_k = T_(k-1) / 2 plus the new midpoints' terms; L1 likewise, which is also the level's finiteness check.
         * The estimate is |R(k, k) - R(k-1, k-1)| and abs_error max(estimate, floor) (design.md, "Romberg (E)").
         *
         * @param fn The function.
         * @param s The state, with level k - 1 in previous and k set.
         * @return Nothing, or non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> romberg_level(F &fn, RombergState<T> &s)
        {
            Failure<T> const overflow{ error::non_finite, std::nullopt };
            // Level 1's midpoint has T_0's weight, 1; each later level halves it, as h halves.
            if (s.k > 1) s.weight = T(s.weight / T(2));
            Sum<T> sum;
            T abs_sum = T(0);
            if (!detail::add_midpoints(fn, s, sum, abs_sum)) return std::unexpected(overflow);
            std::optional<T> const level = detail::compensated(sum);
            if (!level) return std::unexpected(overflow);
            if consteval {
                if (detail::large(*level) || detail::large(abs_sum)) return std::unexpected(overflow);
            }
            s.row[0] = T(T(s.previous[0] / T(2)) + *level);
            // A non-finite value of fn makes L1 non-finite, so L1 is the check per level.
            s.l1 = T(T(s.l1 / T(2)) + abs_sum);
            if (!nxx::math::isfinite(s.l1)) return std::unexpected(overflow);
            if (!detail::extrapolate(s)) return std::unexpected(overflow);
            std::size_t const diagonal = static_cast<std::size_t>(s.k);
            if consteval {
                if (detail::large(s.row[diagonal])) return std::unexpected(overflow);
            }
            // The difference of the last two diagonal entries: it estimates the error of R(k-1, k-1), so it is
            // pessimistic for R(k, k) on smooth fn, and stays honest where Richardson's assumptions fail (sqrt x),
            // where the subdiagonal |R(k, k) - R(k, k-1)| would claim 4^k - 1 times too little. An overflow of T_k or
            // of the extrapolation reaches it (R(k-1, k-1) is finite), so it is checked once.
            T const estimate = nxx::math::abs(T(s.row[diagonal] - s.previous[diagonal - 1]));
            if (!nxx::math::isfinite(estimate)) return std::unexpected(overflow);
            T const floor = detail::rounding_floor(s.l1);
            s.step = Step<T>{ .current = Integral<T>{ s.row[diagonal], detail::larger(estimate, floor) },
                .attainable = floor,
                .ready = s.k >= RombergState<T>::min_level,
                .exhausted = s.k == RombergState<T>::cap };
            return {};
        }

        /**
         * @brief Romberg's start (iterative_method's start): T_0 from the ends, fn(a) then fn(b), then level 1, so
         * that the first Step reported has an estimate.
         *
         * The ends are evaluated as they are, not clamped or moved: a singularity at an end is error::non_finite.
         *
         * @param fn The function.
         * @param a The lower end.
         * @param b The upper end, a < b.
         * @return The state after level 1, or non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<RombergState<T>, Failure<T>> romberg_start(F &fn, T const &a, T const &b)
        {
            Failure<T> const overflow{ error::non_finite, std::nullopt };
            RombergState<T> s{ detail::interval_of(a, b), T(1), T(0), detail::zeros<T, RombergState<T>::row_size>(),
                detail::zeros<T, RombergState<T>::row_size>(), 1, Step<T>{} };
            std::optional<T> const ta = s.interval.checked_term(s.weight, detail::evaluate(fn, a));
            std::optional<T> const tb = s.interval.checked_term(s.weight, detail::evaluate(fn, b));
            if (!ta || !tb) return std::unexpected(overflow);
            s.l1 = T(nxx::math::abs(*ta) + nxx::math::abs(*tb));
            if (!nxx::math::isfinite(s.l1)) return std::unexpected(overflow);
            s.previous[0] = T(*ta + *tb);
            std::expected<void, Failure<T>> const first = detail::romberg_level(fn, s);
            if (!first) return std::unexpected(first.error());
            return s;
        }

        /**
         * @brief Romberg's step (iterative_method's advance): the current row becomes the previous one, and the next
         * level is formed.
         * @param fn The function.
         * @param s The state, below the cap.
         * @return Nothing, or non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> romberg_advance(F &fn, RombergState<T> &s)
        {
            s.previous = s.row;
            ++s.k;
            return detail::romberg_level(fn, s);
        }

        /**
         * @brief A panel of the adaptive Gauss-Kronrod loop: a piece of [a, b] with its own 15-point value and error
         * estimate, the unit the loop sums, ranks and bisects (design.md, "Gauss-Kronrod", "A panel").
         *
         * The integral is the sum of the panels' values and the estimate the sum of their estimates, so a panel keeps
         * what each running sum needs (value, floor, excess) and what the heap needs to rank it (difference, order).
         *
         * @tparam T The real type.
         */
        template<real T>
        struct Panel
        {
            /** @brief The lower end. */
            T lo;

            /** @brief The upper end. */
            T hi;

            /** @brief The 15-point Kronrod value: the panel's share of the integral. */
            T value;

            /**
             * @brief |K15 - G7|, the panel's raw estimate: it ranks the panel in the heap, so that the panel with the
             * largest error is bisected first.
             */
            T difference;

            /**
             * @brief The rounding floor, 50 max(eps L1, denorm_min), L1 the panel's integral of |fn|: the least error
             * the panel can certify; summed into F.
             */
            T floor;

            /**
             * @brief The estimate's excess over the floor, max(difference, floor) - floor: what bisection can still
             * reduce; summed into X, which the stagnation rule watches.
             */
            T excess;

            /**
             * @brief The order of creation, the heap's tie-breaker, so that the order of bisection is deterministic:
             * panels of equal difference leave the heap in it.
             */
            std::size_t order;
        };

        /**
         * @brief Why a panel has no estimate, which decides what the loop does next: a descent around a value that is
         * not finite, or the end of the run (design.md, "Values that are not finite").
         */
        enum class fault {
            not_finite, ///< L1 not finite: fn not finite at a point of the panel (or an overflow of its terms); the
                        ///< loop bisects around it (resolve())
            overflow    ///< a finite L1 but a value or difference that overflowed, or (in constant evaluation) a value
                        ///< or a sum the checks refuse: the run is non_finite
        };

        /**
         * @brief The G7-K15 pair on [lo, hi], lo < hi: one panel's value, its estimate and its floor from 15
         * evaluations of fn (PHASE2_API.md B and C).
         *
         * Both rules run on the same 15 evaluations: the 7 Gauss nodes are among the 15 Kronrod nodes (at the even
         * indices), so the estimate |K15 - G7| costs no extra evaluation. Each rule's value is a compensated sum; L1,
         * the sum of |term| that sets the floor, is a plain sum, since its terms are not negative and do not cancel.
         * The points are evaluated in a fixed order, the centre then c + d and c - d for the nodes ascending, which
         * the prototype's evaluation counts reproduce. One finiteness check per panel, on L1, separates a value of fn
         * that is not finite (not_finite: the loop can bisect around it) from an overflow of the sums (overflow).
         *
         * @param rule The pair.
         * @param fn The function.
         * @param lo The lower end.
         * @param hi The upper end.
         * @param order The panel's order of creation.
         * @return The panel, or why it has no estimate.
         */
        template<real T, typename F>
        constexpr std::expected<Panel<T>, fault>
            sample(Kronrod<T> const &rule, F &fn, T const &lo, T const &hi, std::size_t order)
        {
            Interval<T> const interval = detail::interval_of(lo, hi);
            Sum<T> kronrod;
            Sum<T> gauss;
            T l1 = T(0);
            bool refused = false;
            // One node: its term enters K15, L1 and, at a Gauss node (even index), G7. In constant evaluation a term or
            // a sum the guards refuse stops the panel (overflow) instead of failing to compile.
            auto const add = [&](std::size_t i, T const &x) {
                if (refused) return;
                T const y = detail::evaluate(fn, interval.point(x));
                std::optional<T> const tk = interval.checked_term(rule.wk[i], y);
                std::optional<T> const tg = i % 2 == 0 ? interval.checked_term(rule.wg[i], y) : std::optional<T>(T(0));
                if (!tk || !tg || !detail::fits(kronrod.sum, *tk) || !detail::fits(gauss.sum, *tg)
                    || !detail::fits(l1, *tk)) {
                    refused = true;
                    return;
                }
                kronrod.add(*tk);
                l1 = T(l1 + nxx::math::abs(*tk));
                if (i % 2 == 0) gauss.add(*tg);
            };
            add(0, interval.c);
            for (std::size_t i = 1; i < Kronrod<T>::size; ++i) {
                T const d = T(interval.hw * rule.x[i]);
                add(i, T(interval.c + d));
                add(i, T(interval.c - d));
            }
            if (refused) return std::unexpected(fault::overflow);
            // A value of fn that is not finite makes L1 not finite: the one check that tells it from an overflow.
            if (!nxx::math::isfinite(l1)) return std::unexpected(fault::not_finite);
            std::optional<T> const kv = detail::compensated(kronrod);
            std::optional<T> const gv = detail::compensated(gauss);
            if (!kv || !gv || !detail::fits(*kv, *gv)) return std::unexpected(fault::overflow);
            T const k = *kv;
            T const g = *gv;
            T const difference = nxx::math::abs(T(k - g));
            // Beyond the spec's check of L1: a finite L1 whose value or difference overflows (the 7-point weights
            // are larger than the 15-point ones at the Gauss nodes) ends the run non_finite.
            if (!nxx::math::isfinite(k) || !nxx::math::isfinite(difference)) return std::unexpected(fault::overflow);
            T const floor = detail::rounding_floor(l1);
            return Panel<T>{ lo, hi, k, difference, floor, T(detail::larger(difference, floor) - floor), order };
        }

        /**
         * @brief Whether [lo, hi] can be bisected usefully, so that a panel too narrow to split becomes final instead
         * of being bisected forever (PHASE2_API.md C).
         *
         * It can where its midpoint lies strictly inside and the outer nodes of both children lie strictly inside
         * them. On a panel a few ulps wide the children's nodes round onto their ends: a split would then evaluate fn
         * at the ends, where an endpoint singularity may lie, and gain nothing.
         *
         * @param rule The pair.
         * @param lo The lower end.
         * @param hi The upper end.
         * @return Whether the panel can be split.
         */
        template<real T>
        constexpr bool splittable(Kronrod<T> const &rule, T const &lo, T const &hi)
        {
            T const mid = detail::midpoint(lo, hi);
            if (!(lo < mid && mid < hi)) return false;
            T const outer = rule.x[Kronrod<T>::size - 1];
            auto const inside = [&outer](T const &p, T const &q) {
                Interval<T> const child = detail::interval_of(p, q);
                T const d = T(child.hw * outer);
                return p < T(child.c - d) && T(child.c + d) < q;
            };
            return inside(lo, mid) && inside(mid, hi);
        }

        /**
         * @brief A piece waiting to be sampled in a descent, with its depth, so that the descent stops at its limit of
         * bisections below the panel that started it (KronrodState::depth).
         * @tparam T The real type.
         */
        template<real T>
        struct Pending
        {
            /** @brief The lower end. */
            T lo;

            /** @brief The upper end. */
            T hi;

            /** @brief The bisections below the panel that started the descent. */
            int level;
        };

        /**
         * @brief The three sums of the adaptive loop, compensated: the value K, the floor F and the excess
         * X = sum (e_p - floor_p); the estimate is E = F + X.
         *
         * F and X are kept apart, rather than E alone, because the loop needs each: the stop test accepts at 2 F, and
         * stagnation watches X, which falls as bisection works, where E barely moves once it nears the floor
         * (PHASE2_API.md C and E).
         *
         * @tparam T The real type.
         */
        template<real T>
        struct KronrodTotals
        {
            /** @brief The sum of the panels' 15-point values: the integral. */
            Sum<T> value;

            /** @brief The sum of the panels' floors: the least error the run can certify. */
            Sum<T> floor;

            /** @brief The sum of the panels' excesses over their floors: the error bisection can still reduce. */
            Sum<T> excess;
        };

        /**
         * @brief The values of a set of KronrodTotals, checked finite: what a Step and the stagnation rule are formed
         * from.
         * @tparam T The real type.
         */
        template<real T>
        struct KronrodValues
        {
            /** @brief K, the integral. */
            T value;

            /** @brief F, the floor. */
            T floor;

            /** @brief X, the excess over the floor. */
            T excess;

            /** @brief E = F + X, the estimate. */
            T estimate;
        };

        /**
         * @brief Adaptive Gauss-Kronrod as an iterative method (PHASE2_API.md C to E, PHASE3A_API.md): the state that
         * lets GaussKronrod<15> take one bisection per step and report a Step, while integrate() owns the stop test.
         *
         * The state holds every panel, the heap of the active ones (those above their floor that can still be split)
         * and the running sums. A step bisects the active panel with the largest difference. The rule is copied,
         * because for a type without a table, in constant evaluation, with_gauss_kronrod passes a temporary.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct KronrodState
        {
            /**
             * @brief The panels of the initial mesh, so that a feature one panel's 15 points straddle is still seen
             * (PHASE2_API.md C and Q1, as Shampine's quadgk).
             *
             * One panel can accept K = G on a peak it missed. On 100 Gaussian peaks of width 1e-3 the prototype had 36
             * false successes with 1 panel, 2 with 2, 1 with 4 and none with 10. No mesh cures a feature narrower than
             * the spacing of the nodes.
             */
            static constexpr int mesh = 10;

            /**
             * @brief The bisections below a panel whose L1 is not finite, so that a value of fn that is not finite at a
             * node is stepped around, but a region where fn is not finite is not walked into (PHASE2_API.md Q8).
             */
            static constexpr int depth = 2;

            /**
             * @brief The bisections without a new minimum of the excess X after which the run ends not_converged, so
             * that a noisy fn or a divergent integral stops long before the cap (PHASE2_API.md Q2).
             */
            static constexpr int stagnation = 32;

            /** @brief The G7-K15 pair. */
            Kronrod<T> rule;

            /**
             * @brief The cap on bisections and descent splits together, 64 per bit of T, so that every run ends: x^-0.9
             * on [0, 1] at cpp_bin_float_50's default needs 3,550 of its 10,752.
             */
            int cap;

            /**
             * @brief Every current panel, active and final, in the order it was committed. Removing a parent's
             * contribution from a compensated sum is not exact (a small contribution absorbed into the correction
             * beside much larger ones is lost when they are subtracted again), so every decision to stop is made on
             * sums formed afresh from these panels (settle; independent review, 2098416).
             */
            std::vector<Panel<T>> panels;

            /**
             * @brief Whether each panel was bisected: a bisected parent is flagged rather than erased, so that the
             * indices in the heap stay valid and fresh() skips it.
             */
            std::vector<bool> removed;

            /** @brief The running sums, updated per bisection: cheap, and what stagnation watches between settles. */
            KronrodTotals<T> running;

            /** @brief Whether constant evaluation refused a sum outside max / 2: the run then ends non_finite. */
            bool refused;

            /**
             * @brief The active panels, a heap of indices into panels by difference (the estimate's order), so that
             * the panel with the largest error is found in logarithmic time; a final panel never enters it.
             */
            std::vector<std::size_t> heap;

            /** @brief The order of the next panel sampled, the heap's tie-breaker. */
            std::size_t order;

            /** @brief Bisections and descent splits so far, against the cap. */
            int work;

            /**
             * @brief The panels of the current bisection or mesh panel, before they are committed: held apart, so that
             * a bisection whose descent fails can be undone by not committing them.
             */
            std::vector<Panel<T>> pieces;

            /** @brief The work list of a descent below a panel whose L1 is not finite. */
            std::vector<Pending<T>> stack;

            /**
             * @brief Whether a piece no descent avoided was undone: then the run cannot end as a success, since the
             * estimate no longer covers that piece (Step::trusted).
             */
            bool tainted;

            /** @brief The least excess so far, from the completed mesh on: stagnation's reference. */
            T least;

            /** @brief Bisections since the excess last reached a new minimum; at stagnation the run is exhausted. */
            int since;

            /** @brief The step: {K, max(E, F)}, 2F attainable, exhausted by an empty heap, the cap or stagnation. */
            Step<T> step;

            /**
             * @brief An empty state with the rule and the cap of T; kronrod_start() samples the mesh into it.
             * @param r The G7-K15 pair.
             */
            constexpr explicit KronrodState(Kronrod<T> const &r)
                : rule(r), cap(64 * detail::significand_bits<T>()), panels(), removed(), running(), refused(false),
                  heap(), order(0), work(0), pieces(), stack(), tainted(false), least(T(0)), since(0), step{}
            {}

            /**
             * @brief Adds a panel's contribution to a set of sums, or subtracts it (a bisected parent's); refuses
             * (refused) a sum constant evaluation could not form within max / 2.
             * @param t The sums.
             * @param p The panel.
             * @param add Whether to add (else subtract).
             */
            constexpr void contribute(KronrodTotals<T> &t, Panel<T> const &p, bool add)
            {
                T const v = add ? p.value : T(-p.value);
                T const f = add ? p.floor : T(-p.floor);
                T const x = add ? p.excess : T(-p.excess);
                if (!detail::fits(t.value.sum, v) || !detail::fits(t.floor.sum, f) || !detail::fits(t.excess.sum, x)) {
                    refused = true;
                    return;
                }
                t.value.add(v);
                t.floor.add(f);
                t.excess.add(x);
            }

            /**
             * @brief The sums formed afresh from the current panels, in the order they were committed: free of the
             * drift that subtracting parents leaves in the running sums, for every decision to stop.
             * @return The sums.
             */
            constexpr KronrodTotals<T> fresh()
            {
                KronrodTotals<T> t;
                std::ranges::for_each(std::views::iota(std::size_t{ 0 }, panels.size())
                                          | std::views::filter([this](std::size_t i) { return !removed[i]; }),
                    [&](std::size_t i) { contribute(t, panels[i], true); });
                return t;
            }

            /**
             * @brief The heap's order: by difference, the same order as the estimate (only panels above their floor
             * enter the heap), then by creation. With std::push_heap's max-heap, the panel with the largest difference
             * is bisected next, and of panels with equal differences the one created first.
             * @param p A panel's index.
             * @param q Another's.
             * @return Whether p ranks below q.
             */
            constexpr bool precedes(std::size_t p, std::size_t q) const
            {
                Panel<T> const &u = panels[p];
                Panel<T> const &v = panels[q];
                return u.difference < v.difference || (u.difference == v.difference && v.order < u.order);
            }

            /**
             * @brief [p, q] sampled into pieces, stepping around values of fn that are not finite: a singularity at a
             * node is usually avoided by one or two bisections, which make it an end of the children, where no node
             * lies (design.md, "Values that are not finite").
             *
             * The pieces are the panel itself, or, where its L1 is not finite, those of a descent of at most depth
             * bisections below it, each counted against the cap, each child first its lower half. The descent is
             * bounded, since an unbounded one walks towards a region where fn is not finite (sqrt(x - 0.5) on [0, 1]
             * walked to 0.5 through 1,063 levels in double).
             *
             * @param fn The function.
             * @param p The lower end.
             * @param q The upper end.
             * @return The fault that remains, if any: overflow ends the run; not_finite, a piece no descent avoided.
             */
            template<typename F>
            constexpr std::optional<fault> resolve(F &fn, T const &p, T const &q)
            {
                stack.clear();
                stack.push_back(Pending<T>{ p, q, 0 });
                bool missing = false;
                while (!stack.empty()) {
                    Pending<T> const top = stack.back();
                    stack.pop_back();
                    std::expected<Panel<T>, fault> const s = detail::sample(rule, fn, top.lo, top.hi, order++);
                    if (s) {
                        pieces.push_back(*s);
                        continue;
                    }
                    if (s.error() == fault::overflow) return fault::overflow;
                    if (top.level < depth && detail::splittable(rule, top.lo, top.hi) && work < cap) {
                        ++work;
                        T const mid = detail::midpoint(top.lo, top.hi);
                        // The upper half is pushed first, so the lower half is sampled first.
                        stack.push_back(Pending<T>{ mid, top.hi, top.level + 1 });
                        stack.push_back(Pending<T>{ top.lo, mid, top.level + 1 });
                    }
                    else missing = true;
                }
                if (missing) return fault::not_finite;
                return std::nullopt;
            }

            /**
             * @brief The pieces enter the panels and the running sums, and the heap unless final: at its floor
             * (splitting cannot lower it below its floor, since L1 adds up), or not splittable. A final panel still
             * counts in the sums; it is never bisected again.
             */
            constexpr void commit()
            {
                auto const before = [this](std::size_t p, std::size_t q) { return precedes(p, q); };
                std::ranges::for_each(pieces, [&](Panel<T> const &piece) {
                    panels.push_back(piece);
                    removed.push_back(false);
                    contribute(running, piece, true);
                    if (!(piece.difference <= piece.floor) && detail::splittable(rule, piece.lo, piece.hi)) {
                        heap.push_back(panels.size() - 1);
                        std::push_heap(heap.begin(), heap.end(), before);
                    }
                });
            }

            /**
             * @brief The values of a set of sums, checked once: a value of fn or an overflow that made a sum not
             * finite shows here.
             * @param t The sums.
             * @return The values, or nothing where they are not finite or constant evaluation refuses them.
             */
            static constexpr std::optional<KronrodValues<T>> values_of(KronrodTotals<T> const &t)
            {
                std::optional<T> const v = detail::compensated(t.value);
                std::optional<T> const f = detail::compensated(t.floor);
                std::optional<T> const x = detail::compensated(t.excess);
                if (!v || !f || !x || !detail::fits(*f, *x)) return std::nullopt;
                T const e = T(*f + *x);
                if (!nxx::math::isfinite(*v) || !nxx::math::isfinite(e) || !nxx::math::isfinite(*f))
                    return std::nullopt;
                return KronrodValues<T>{ *v, *f, *x, e };
            }

            /**
             * @brief The step from a set of values: {K, max(E, F)} (a Neumaier residue can leave E a hair below F),
             * and twice the floor as the attainable error, since a kink or a jump never brings every panel to its
             * floor (X = 0), and on a smooth fn the test cannot accept an error above max(request, 2 F).
             * @param v The values.
             */
            constexpr void report(KronrodValues<T> const &v)
            {
                step = Step<T>{ .current = Integral<T>{ v.value, detail::larger(v.estimate, v.floor) },
                    .attainable = T(T(2) * v.floor),
                    .exhausted = heap.empty() || work >= cap || since >= stagnation,
                    .trusted = !tainted };
            }
        };

        /**
         * @brief Gauss-Kronrod's start: the initial mesh of 10 equal panels, 150 evaluations, so that a narrow feature
         * a single panel would miss is seen before any estimate is trusted (PHASE2_API.md C, Q1).
         *
         * The mesh points are clamped into [a, b], and the last is b itself, so the panels tile the interval. Equal
         * adjacent points (an interval too narrow for 10) give empty panels, which are skipped. A mesh panel with no
         * finite estimate, after its descent, is non_finite: there is no estimate of the whole to fall back on.
         *
         * @param rule The G7-K15 pair.
         * @param fn The function.
         * @param a The lower end.
         * @param b The upper end, a < b.
         * @return The state after the mesh, or non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<KronrodState<T>, Failure<T>>
            kronrod_start(Kronrod<T> const &rule, F &fn, T const &a, T const &b)
        {
            Failure<T> const invalid{ error::non_finite, std::nullopt };
            constexpr int mesh = KronrodState<T>::mesh;
            KronrodState<T> s(rule);
            Interval<T> const whole = detail::interval_of(a, b);
            T left = a;
            bool const meshed = std::ranges::all_of(std::views::iota(1, mesh + 1), [&](int i) {
                T const u = T(T(-1) + T(T(2 * i) / T(mesh)));
                T const right = i == mesh ? b : whole.point(T(whole.c + T(whole.hw * u)));
                if (left < right) {
                    s.pieces.clear();
                    if (s.resolve(fn, left, right)) return false;
                    s.commit();
                    if (s.refused) return false;
                }
                left = right;
                return true;
            });
            if (!meshed) return std::unexpected(invalid);
            std::optional<KronrodValues<T>> const after_mesh = KronrodState<T>::values_of(s.running);
            if (!after_mesh) return std::unexpected(invalid);
            // Stagnation counts from the completed mesh's excess (independent review of Phase 2, finding 3).
            s.least = after_mesh->excess;
            s.since = 0;
            s.report(*after_mesh);
            return s;
        }

        /**
         * @brief Gauss-Kronrod's step: the panel with the largest estimate bisected, where the error is, and the
         * stagnation count updated (PHASE2_API.md C to E).
         *
         * A bisection whose descent leaves a piece where fn is not finite is undone: the parent stays, final, with its
         * own estimate, and the run is tainted. Stagnation, 32 bisections without a new minimum of the excess X, ends
         * a run that bisection no longer improves: fn noisier than the request, or a divergent integral.
         *
         * @param fn The function.
         * @param s The state, not exhausted (so the heap is not empty).
         * @return Nothing, or non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> kronrod_advance(F &fn, KronrodState<T> &s)
        {
            Failure<T> const invalid{ error::non_finite, std::nullopt };
            auto const before = [&s](std::size_t p, std::size_t q) { return s.precedes(p, q); };
            std::pop_heap(s.heap.begin(), s.heap.end(), before);
            std::size_t const index = s.heap.back();
            s.heap.pop_back();
            Panel<T> const parent = s.panels[index];
            T const mid = detail::midpoint(parent.lo, parent.hi);
            ++s.work;
            s.pieces.clear();
            std::optional<fault> const lower_half = s.resolve(fn, parent.lo, mid);
            if (lower_half == fault::overflow) return std::unexpected(invalid);
            std::optional<fault> const upper_half = s.resolve(fn, mid, parent.hi);
            if (upper_half == fault::overflow) return std::unexpected(invalid);
            if (lower_half || upper_half) {
                // Undone: the pieces are not committed, and the parent, already out of the heap, stays as it was.
                s.tainted = true;
            }
            else {
                s.removed[index] = true;
                s.contribute(s.running, parent, false);
                s.commit();
                if (s.refused) return std::unexpected(invalid);
                // A removed contribution larger than what remains of a running sum leaves that sum mostly rounding:
                // resynchronise, so that stagnation, too, sees the current values.
                std::optional<KronrodValues<T>> const after = KronrodState<T>::values_of(s.running);
                if (!after) return std::unexpected(invalid);
                if (nxx::math::abs(parent.value) > nxx::math::abs(after->value) || parent.floor > after->floor
                    || parent.excess > after->excess) {
                    s.running = s.fresh();
                    if (s.refused) return std::unexpected(invalid);
                }
            }
            std::optional<KronrodValues<T>> const current = KronrodState<T>::values_of(s.running);
            if (!current) return std::unexpected(invalid);
            // X, not E: E barely moves near the floor, and a peak found late raises it (PHASE2_API.md E).
            if (current->excess < s.least) {
                s.least = current->excess;
                s.since = 0;
            }
            else ++s.since;
            s.report(*current);
            return {};
        }

        /**
         * @brief Gauss-Kronrod's settle: the step from the sums formed afresh, which then replace the running ones, so
         * that integrate() decides, and the result is taken, on sums free of the running sums' drift.
         *
         * Where the fresh sums fail the stop test that the running ones passed, the loop simply goes on from them.
         *
         * @param s The state.
         * @return Nothing, or non_finite.
         */
        template<real T>
        constexpr std::expected<void, Failure<T>> kronrod_settle(KronrodState<T> &s)
        {
            Failure<T> const invalid{ error::non_finite, std::nullopt };
            KronrodTotals<T> const exact = s.fresh();
            if (s.refused) return std::unexpected(invalid);
            std::optional<KronrodValues<T>> const checked = KronrodState<T>::values_of(exact);
            if (!checked) return std::unexpected(invalid);
            s.running = exact;
            s.report(*checked);
            return {};
        }

        /**
         * @brief A value carried to about twice T's precision, as the unevaluated sum hi + lo of two T, so that ln 2 is
         * known well enough for exp()'s reduction in any T (design.md, "TanhSinh": exp in T).
         *
         * exp() reduces x by k ln 2 with k up to a few hundred, or millions in a type with a wide exponent range; ln 2
         * rounded to T would make the reduced argument err by k times its rounding. ln 2 cannot come from a literal: a
         * double literal holds no more than a double's digits, which is too few for a 50-digit type. So ln_two() sums
         * its series in Pairs, and exp_constants() splits the Pair into c1 and c2. |lo| is at most about half an ulp
         * of hi.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct Pair
        {
            /** @brief The leading part: the value rounded to T. */
            T hi;

            /** @brief The rest: what rounding to T left out. */
            T lo;
        };

        /**
         * @brief @p a + @p b as a Pair, exactly, for operands of any magnitudes (Knuth's two-sum): the building block
         * of pair_sum().
         * @param a A value.
         * @param b A value.
         * @return The rounded sum and its rounding error.
         */
        template<real T>
        constexpr Pair<T> two_sum(T const &a, T const &b)
        {
            T const s = T(a + b);
            T const b_part = T(s - a);
            return Pair<T>{ s, T(T(a - T(s - b_part)) + T(b - b_part)) };
        }

        /**
         * @brief @p a + @p b as a Pair, exactly, for |a| >= |b| (Dekker's fast two-sum): the cheaper form, used where
         * the order of the magnitudes is known, as when a Pair is renormalised.
         * @param a The larger value.
         * @param b The smaller value.
         * @return The rounded sum and its rounding error.
         */
        template<real T>
        constexpr Pair<T> fast_two_sum(T const &a, T const &b)
        {
            T const s = T(a + b);
            return Pair<T>{ s, T(b - T(s - a)) };
        }

        /**
         * @brief radix^n, exactly: the scale of Veltkamp's splitter in two_product(), in T's own radix, so that the
         * split holds in radix 10 as in radix 2.
         * @param n The power, at least 0.
         * @return radix^n.
         */
        template<real T>
        constexpr T radix_power(int n)
        {
            T const radix = T(std::numeric_limits<T>::radix);
            T power = T(1);
            for (int i = 0; i < n; ++i) power = T(power * radix);
            return power;
        }

        /**
         * @brief @p a @p b as a Pair, exactly, without a fused multiply-add: the product that pair_quotient() and
         * pair_product() need, from + - * / alone as nxx::real allows.
         *
         * Dekker's product: each factor is split by Veltkamp's method into a high part of about half its digits and the
         * rest, with the splitter radix^ceil(digits / 2) + 1, so that the products of the parts are exact, and their
         * sum less the rounded product is its rounding error. It holds in any radix.
         *
         * @param a A value.
         * @param b A value.
         * @return The rounded product and its rounding error.
         */
        template<real T>
        constexpr Pair<T> two_product(T const &a, T const &b)
        {
            T const splitter = T(detail::radix_power<T>((std::numeric_limits<T>::digits + 1) / 2) + T(1));
            T const p = T(a * b);
            T const ca = T(splitter * a);
            T const a_hi = T(ca - T(ca - a));
            T const a_lo = T(a - a_hi);
            T const cb = T(splitter * b);
            T const b_hi = T(cb - T(cb - b));
            T const b_lo = T(b - b_hi);
            T e = T(T(a_hi * b_hi) - p);
            e = T(e + T(a_hi * b_lo));
            e = T(e + T(a_lo * b_hi));
            e = T(e + T(a_lo * b_lo));
            return Pair<T>{ p, e };
        }

        /**
         * @brief @p a + @p b for Pairs (the double-word addition): the running sum of ln_two()'s series.
         * @param a A Pair.
         * @param b A Pair.
         * @return The sum, renormalised.
         */
        template<real T>
        constexpr Pair<T> pair_sum(Pair<T> const &a, Pair<T> const &b)
        {
            // The leading parts and the trailing parts are added exactly, then the errors are folded in and the
            // result renormalised so that |lo| is again at most about half an ulp of hi.
            Pair<T> const s = detail::two_sum(a.hi, b.hi);
            Pair<T> const t = detail::two_sum(a.lo, b.lo);
            T e = T(s.lo + t.hi);
            Pair<T> const u = detail::fast_two_sum(s.hi, e);
            e = T(t.lo + u.lo);
            return detail::fast_two_sum(u.hi, e);
        }

        /**
         * @brief @p a / @p n for a small integer n held exactly in T: the divisions by 3 and by 2i + 1 in ln_two()'s
         * series.
         * @param a A Pair.
         * @param n The divisor.
         * @return The quotient.
         */
        template<real T>
        constexpr Pair<T> pair_quotient(Pair<T> const &a, T const &n)
        {
            // A first quotient, then the exact remainder a - q1 n (two_product), divided again for the low part.
            T const q1 = T(a.hi / n);
            Pair<T> const p = detail::two_product(q1, n);
            T const r = T(T(T(a.hi - p.hi) - p.lo) + a.lo);
            return detail::fast_two_sum(q1, T(r / n));
        }

        /**
         * @brief @p a @p n for a small integer n held exactly in T: the final doubling in ln_two().
         * @param a A Pair.
         * @param n The factor.
         * @return The product.
         */
        template<real T>
        constexpr Pair<T> pair_product(Pair<T> const &a, T const &n)
        {
            Pair<T> const p = detail::two_product(a.hi, n);
            return detail::fast_two_sum(p.hi, T(p.lo + T(a.lo * n)));
        }

        /**
         * @brief ln 2 to about twice T's precision, computed in T itself, so that exp() needs no constant beyond what
         * + - * / can form (BRIEF.md 3.10, D11).
         *
         * ln 2 = 2 atanh(1/3) = 2 sum_i 3^-(2i+1) / (2i + 1): a series of rational terms that shrink by a factor of
         * about 9 each, summed in Pairs to a term below eps^2 / 16. It works in any radix and any precision, where a
         * literal would hold only a double's digits.
         *
         * @return ln 2 as a Pair.
         */
        template<real T>
        constexpr Pair<T> ln_two()
        {
            T const eps = T(std::numeric_limits<T>::epsilon());
            // At least min(): in a type whose exponents end before eps^2 (min() 2^-100, eps 2^-99) eps^2 / 16 is 0,
            // which no term goes below.
            T const stop = detail::larger(T(T(eps * eps) / T(16)), T((std::numeric_limits<T>::min)()));
            T const three = T(3);
            Pair<T> power = detail::pair_quotient(Pair<T>{ T(1), T(0) }, three);
            Pair<T> sum{ T(0), T(0) };
            for (int i = 0;; ++i) {
                Pair<T> const term = detail::pair_quotient(power, T(2 * i + 1));
                if (term.hi < stop) break;
                sum = detail::pair_sum(sum, term);
                power = detail::pair_quotient(detail::pair_quotient(power, three), three);
            }
            return detail::pair_product(sum, T(2));
        }

        /**
         * @brief Everything exp() needs that depends only on T, computed once by exp_constants() and passed in, so that
         * the many calls that build a table do not recompute it.
         *
         * The reduction x = k ln 2 + r by Cody and Waite's split of ln 2 into c1 + c2; tables of powers of 2 that
         * stand in for k, since T has no conversion to an integer (nxx::real does not require one); and the length
         * of the Taylor series of expm1 on r.
         *
         * The domain is the tables' (PHASE3_API.md B): from x_min = ln 2^-m moved a few eps inward, 2^-m the smallest
         * power of 2 at or above min(), up to about 22. k counts powers of 2, at most -log2 min() <= (1 - min_exponent)
         * times the bits of a digit; reach adds a margin.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct ExpConstants
        {
            /** @brief The largest |k| the reduction meets, with a margin: the size of the tables of powers of 2. */
            static constexpr long long reach =
                (1LL - std::numeric_limits<T>::min_exponent)
                    * static_cast<long long>(detail::significand_bits<T>() / std::numeric_limits<T>::digits)
                + 4;

            /** @brief The bits of reach: the length of the tables, one entry per bit of k. */
            static constexpr std::size_t bits =
                static_cast<std::size_t>(std::bit_width(static_cast<unsigned long long>(reach)));

            /** @brief 1 / ln 2, rounded: k is the nearest integer to x / ln 2. */
            T inverse_ln2;

            /**
             * @brief ln 2 cut to digits - m radix digits, m the digits of reach, so that k c1 is exact for every k the
             * domain needs and x - k c1 loses nothing.
             */
            T c1;

            /** @brief The rest of ln 2, from the Pair: ln 2 - c1, whose product with k carries the only rounding. */
            T c2;

            /**
             * @brief Where the domain ends below: exp(x) is 0 for x <= x_min, so that every nonzero result is at least
             * min() and no table holds a subnormal value.
             */
            T x_min;

            /** @brief The terms of the Taylor series of expm1 for |r| <= ln 2 / 2: until one is below eps / 8. */
            int terms;

            /** @brief 2^i: the values the greedy search for k compares against. */
            std::array<T, bits> powers;

            /**
             * @brief The factors of 2^k for k > 0: 2 for bit 0, 2^(2^(i-1)) for bit i applied twice, so that 2^k is
             * formed as two factors that each stay in range; 0 past max(), never taken.
             */
            std::array<T, bits> up;

            /** @brief Their inverses, for k < 0. */
            std::array<T, bits> down;
        };

        /**
         * @brief The constants of exp() in T, from + - * / alone, so that they can be formed in constant evaluation and
         * in any nxx::real.
         * @return The constants.
         */
        template<real T>
        constexpr ExpConstants<T> exp_constants()
        {
            using limits = std::numeric_limits<T>;
            constexpr std::size_t bits = ExpConstants<T>::bits;
            T const radix = T(limits::radix);
            T const eps = T(limits::epsilon());
            Pair<T> const ln2 = detail::ln_two<T>();

            // Cody and Waite: k c1 must be exact for every |k| <= reach, so c1 keeps digits - m radix digits of ln 2,
            // m the radix digits of reach, cut by subtracting units of the leading digit's place downward.
            int m = 0;
            for (long long held = 1; held <= ExpConstants<T>::reach; held *= limits::radix) ++m;
            int const keep = limits::digits - m;
            T leading = T(1);
            while (leading > ln2.hi) leading = T(leading / radix);
            while (T(leading * radix) <= ln2.hi) leading = T(leading * radix);
            T rest = ln2.hi;
            T c1 = T(0);
            T unit = leading;
            for (int i = 0; i < keep; ++i) {
                while (rest >= unit) {
                    rest = T(rest - unit);
                    c1 = T(c1 + unit);
                }
                unit = T(unit / radix);
            }
            T const c2 = T(rest + ln2.lo);

            std::array<T, bits> powers = detail::zeros<T, bits>();
            std::array<T, bits> up = detail::zeros<T, bits>();
            std::array<T, bits> down = detail::zeros<T, bits>();
            powers[0] = T(1);
            for (std::size_t i = 1; i < bits; ++i) powers[i] = T(powers[i - 1] + powers[i - 1]);
            // 2^k as two factors, each about 2^(k/2): bit 0 contributes 2 once, bit i > 0 contributes 2^(2^(i-1))
            // twice, so that neither factor leaves the range where the product is still inside it.
            up[0] = T(2);
            down[0] = T(T(1) / T(2));
            up[1] = up[0];
            down[1] = down[0];
            // An upward factor only where it fits (up[i - 1] <= max / up[i - 1]), else 0: in a type whose exponents
            // reach further below 0 than above, reach outruns the largest power of 2 (the Phase 4 spike's type with
            // max_exponent 100 formed 2^128). exp() never takes such a factor: its arguments, up to about 22, need
            // k < 2^5.
            T const largest = T((limits::max)());
            for (std::size_t i = 2; i < bits; ++i) {
                up[i] = up[i - 1] > T(0) && up[i - 1] <= T(largest / up[i - 1]) ? T(up[i - 1] * up[i - 1]) : T(0);
                down[i] = T(down[i - 1] * down[i - 1]);
            }

            // m, the largest with 2^-m >= min(), by a greedy descent over its bits from min() itself, not min_exponent
            // (Boost's cpp_bin_float_50 has min() = 2^min_exponent, a binade above radix^(min_exponent - 1)): bit i is
            // taken if 2^-(m + 2^i) >= min(), formed from the downward factors. Nothing overflows, in a type whose
            // exponents reach further below 0 than above too (2^(2^i) did there); a product that underflows lies below
            // min() and is rejected, as it should be. x_min = -m ln 2, moved 4 eps inward, so that the rounded
            // reduction cannot give a result below min().
            T const smallest = T((limits::min)());
            T scale = T(1); // 2^-m
            T exponent = T(0);
            for (std::size_t i = bits; i-- > 0;) {
                T const candidate = i > 0 ? T(T(scale * down[i]) * down[i]) : T(scale * down[0]);
                if (candidate >= smallest) {
                    scale = candidate;
                    exponent = T(exponent + powers[i]);
                }
            }
            T const ln_min = T(T(exponent * c1) + T(exponent * c2));
            T const x_min = T(-T(ln_min - T(T(4) * T(eps * ln_min))));

            // |r| <= ln 2 / 2: the Taylor terms of expm1 until the next one is below eps / 8.
            T const half_ln2 = T(ln2.hi / T(2));
            T term = T(1);
            int terms = 0;
            while (term > T(eps / T(8))) {
                ++terms;
                term = T(T(term * half_ln2) / T(terms));
            }
            return ExpConstants<T>{ T(T(1) / ln2.hi), c1, c2, x_min, terms, powers, up, down };
        }

        /**
         * @brief e^x on the tables' domain, the library's own: the tanh-sinh and exp-sinh tables are built from it
         * (design.md, "TanhSinh": exp in T; PHASE3_API.md B and Q10).
         *
         * std::exp will not do. T may be any nxx::real, which offers + - * / and comparisons but no exp; the tables
         * are built in constant evaluation, where std::exp is not constexpr on Clang and MSVC; and a library's exp
         * can be poor where the tables need it (Boost's in cpp_bin_float_50 errs by 1,568 ulp at x = 1096). This one
         * uses + - * / in T alone, so the same code gives every table, at compile time or at run time.
         *
         * x = k ln 2 + r with k the nearest integer to x / ln 2, found by a greedy comparison against powers of 2,
         * r = (x - k c1) - k c2, and e^x = (1 + p(r)) 2^k, p the series by Horner's scheme. The error: the final
         * rounding of 1 + p (0.5 ulp), Horner's (0.73), the rounding of r (0.25): at most 1.48 ulp.
         *
         * Only the tables' domain, as Q10 decided, which keeps it short: no NaN test, no overflow guard, no subnormal
         * results. x <= x_min gives 0, and x up to about 22 (the tables need t_max) is far from overflow.
         *
         * @param x The argument, x_min < x <= 22 (or x <= x_min for 0).
         * @param constants The constants.
         * @return e^x.
         */
        template<real T>
        constexpr T exp(T const &x, ExpConstants<T> const &constants)
        {
            if (!(x > constants.x_min)) return T(0);
            T const y = T(x * constants.inverse_ln2);
            bool const negative = y < T(0);
            // k = round(|y|): add 1/2, then take the powers of 2 it holds, largest first (T has no floor or
            // conversion to an integer). Each power taken also takes its factor of 2^k.
            T rest = T(T(negative ? T(-y) : y) + T(T(1) / T(2)));
            std::array<T, ExpConstants<T>::bits> const &factors = negative ? constants.down : constants.up;
            T k = T(0);
            T scale_1 = T(1);
            T scale_2 = T(1);
            for (std::size_t i = constants.powers.size(); i-- > 0;) {
                if (rest >= constants.powers[i]) {
                    rest = T(rest - constants.powers[i]);
                    k = T(k + constants.powers[i]);
                    scale_1 = T(scale_1 * factors[i]);
                    if (i > 0) scale_2 = T(scale_2 * factors[i]);
                }
            }
            assert(scale_1 > T(0)); // no upward factor left 0 by exp_constants() was taken
            T const signed_k = negative ? T(-k) : k;
            T const r = T(T(x - T(signed_k * constants.c1)) - T(signed_k * constants.c2));
            // expm1(r) by Horner's scheme: p = r (1 + r/2 (1 + r/3 (1 + ...))).
            T p = T(0);
            for (int n = constants.terms; n >= 1; --n) p = T(T(r / T(n)) * T(T(1) + p));
            return T(T(T(T(1) + p) * scale_1) * scale_2);
        }

        /**
         * @brief sinh t and cosh t together, as the tables need both for a node: sinh t places it, cosh t weights it.
         * @tparam T The real type.
         */
        template<real T>
        struct Hyperbolic
        {
            /** @brief sinh t. */
            T s;

            /** @brief cosh t. */
            T c;
        };

        /**
         * @brief sinh t and cosh t for t >= 0, accurate near 0 as well, from one call of exp().
         *
         * Below t = 1/2, by their Taylor series, since there e - 1/e cancels and would lose the digits of a small
         * sinh t; else from e = exp(t) as (e - 1/e) / 2 and (e + 1/e) / 2.
         *
         * @param t The argument, t >= 0.
         * @param constants exp()'s constants.
         * @return sinh t and cosh t.
         */
        template<real T>
        constexpr Hyperbolic<T> sinh_cosh(T const &t, ExpConstants<T> const &constants)
        {
            if (t < T(T(1) / T(2))) {
                // At least min(), as ln_two()'s stop.
                T const eighth = detail::larger(T(T(std::numeric_limits<T>::epsilon()) / T(8)),
                    T((std::numeric_limits<T>::min)()));
                T const tt = T(t * t);
                T odd = t;
                T even = T(1);
                T s = t;
                T c = T(1);
                // Stop when the next terms are below eps / 8 of their sums: cosh t >= 1, so the even term's test is
                // relative too.
                for (int n = 1;; ++n) {
                    even = T(T(even * tt) / T(T(2 * n - 1) * T(2 * n)));
                    odd = T(T(odd * tt) / T(T(2 * n) * T(2 * n + 1)));
                    c = T(c + even);
                    s = T(s + odd);
                    if (even < eighth && odd <= T(s * eighth)) break;
                }
                return Hyperbolic<T>{ s, c };
            }
            T const e = detail::exp(t, constants);
            T const inverse = T(T(1) / e);
            return Hyperbolic<T>{ T(T(e - inverse) / T(2)), T(T(e + inverse) / T(2)) };
        }

        /**
         * @brief The levels of the tanh-sinh rule and of the exp-sinh table: h = 1 down to 2^-9, the same for every
         * type (PHASE3_API.md Q2: nine levels lose nothing in double, but in binary128 at rel_tol 1e-32 they lose 39 of
         * 340 cases, and 168 bits needs the tenth too).
         */
        inline constexpr int tanh_sinh_levels = 10;

        /**
         * @brief The most nodes in a piece, the unit of the tables: each piece is its own constexpr initializer, so
         * that each stays within the compiler's constexpr budget, which a whole table exceeds (PHASE3_API.md Q3: in
         * pieces the largest takes at most 42 % of Clang's or Emscripten's default limit).
         */
        inline constexpr int piece_size = 1024;

        /**
         * @brief Ten times an upper bound on the tables' t_max, computed from T's limits in integers, so that the
         * pieces' arrays can be sized at compile time before any node is generated.
         *
         * A node is kept while 2 sinh t < L in the tanh-sinh table and while (3/2) sinh t < L in the exp-sinh table
         * (both of its sides keep q >= min()), L = -ln min() <= n ln 2, n = (1 - min_exponent) times the bits of a
         * digit rounded up (4 in radix 10, so that n ln 2 >= (1 - min_exponent) ln 10), so t_max <= asinh(2L / 3) <=
         * ln(4L / 3 + 1) <= ln n (n >= 14) < 0.7 bit_width(n).
         *
         * @return 10 times the bound.
         */
        template<real T>
        constexpr int reach_tenths()
        {
            unsigned long long const n = static_cast<unsigned long long>(1LL - std::numeric_limits<T>::min_exponent)
                                         * static_cast<unsigned long long>(detail::significand_bits<T>()
                                                                           / std::numeric_limits<T>::digits);
            return 7 * static_cast<int>(std::bit_width(n));
        }

        /**
         * @brief An upper bound on the nodes of a level, which sizes its pieces: t = 0, 1, ... at level 0, the odd
         * multiples of 2^-level above, up to reach_tenths() / 10.
         * @param level The level.
         * @return The bound.
         */
        template<real T>
        constexpr int level_capacity(int level)
        {
            int const tenths = detail::reach_tenths<T>();
            if (level == 0) return tenths / 10 + 1;
            return ((tenths << level) / 10 + 1) / 2;
        }

        /**
         * @brief The pieces of a level: how many initializers hold it.
         * @param level The level.
         * @return ceil(level_capacity / piece_size).
         */
        template<real T>
        constexpr int level_pieces(int level)
        { return (detail::level_capacity<T>(level) + piece_size - 1) / piece_size; }

        /**
         * @brief The capacity of a piece: piece_size, or, for a level's last piece, what is left of the level's bound,
         * so that no piece's array is larger than it can need.
         * @param level The level.
         * @param piece The piece.
         * @return The capacity.
         */
        template<real T>
        constexpr std::size_t piece_capacity(int level, int piece)
        {
            int const left = detail::level_capacity<T>(level) - piece * piece_size;
            return static_cast<std::size_t>(left < piece_size ? left : piece_size);
        }

        /**
         * @brief The storage of a column of a piece of N nodes, chosen by where the piece is built.
         *
         * std::array for float, double and long double, whose pieces are compile-time tables: a constexpr variable
         * can hold an array, but not a std::vector, whose allocation cannot outlive constant evaluation. std::vector
         * for every other T, whose pieces are built at run time, so that a piece lives on the heap: a column of 1,024
         * cpp_bin_float_50 values takes about 40 KB, and a piece holds two (tanh-sinh) or four (exp-sinh), past
         * Emscripten's default stack of 64 KB.
         *
         * @tparam T The real type.
         * @tparam N The capacity.
         */
        template<real T, std::size_t N>
        using Column = std::conditional_t<tabled<T>, std::array<T, N>, std::vector<T>>;

        /**
         * @brief N nodes of value 0 in either storage, the starting value of a piece's column (T need not be default
         * constructible).
         * @return The column.
         */
        template<real T, std::size_t N>
        constexpr Column<T, N> column()
        {
            if constexpr (tabled<T>) return detail::zeros<T, N>();
            else return std::vector<T>(N, T(0));
        }

        /**
         * @brief What a step of the tanh-sinh rule reads per node, for a piece of a level, t ascending: each node's
         * complement 1 - x and its weight dx/dt (without the step h).
         *
         * The complement is stored instead of x because near the ends x rounds to 1 and loses its distance to the
         * end, which is what an integrand singular at an end needs: the points are formed as b - hw comp and
         * a + hw comp.
         *
         * @tparam T The real type.
         * @tparam N The capacity.
         */
        template<real T, std::size_t N>
        struct Nodes
        {
            /** @brief 1 - x = 2q / (1 + q), q = exp(-2 sinh t), formed directly, never as 1 - x. */
            Column<T, N> comp;

            /** @brief dx/dt = cosh t 4q / (1 + q)^2. */
            Column<T, N> w;

            /** @brief The nodes held: N, or fewer where the table ends. */
            std::size_t size;
        };

        /**
         * @brief Generates a piece of a level of the tanh-sinh table, the nodes x = tanh(sinh t), from exp() alone; the
         * one generator for both the compile-time and the run-time tables (PHASE3_API.md B).
         *
         * The map is x = tanh(c sinh t) with c = 1 rather than the literature's pi/2: no pi is needed in T, and it
         * costs fewer evaluations on smooth integrands (PHASE3_API.md B). Level 0 holds t = i, level L >= 1 the odd
         * multiples t = (2i + 1) 2^-L, i = piece 1024 + index. The table ends at the first node whose
         * q = exp(-2 sinh t) is 0, that is below x_min: q >= min() for the nodes kept, so that no stored value is
         * subnormal (subnormal values would lose their digits for no change in any result). Nothing overflows: q only
         * underflows.
         *
         * @tparam T The real type.
         * @tparam L The level, 0 to 9.
         * @tparam P The piece.
         * @return The piece.
         */
        template<real T, int L, int P>
        constexpr Nodes<T, detail::piece_capacity<T>(L, P)> tanh_sinh_piece()
        {
            constexpr std::size_t capacity = detail::piece_capacity<T>(L, P);
            return detail::with_generated<&detail::exp_constants<T>>([](ExpConstants<T> const &constants) {
                Nodes<T, capacity> nodes{ detail::column<T, capacity>(), detail::column<T, capacity>(), 0 };
                T h = T(1);
                for (int i = 0; i < L; ++i) h = T(h / T(2));
                for (std::size_t index = 0; index < capacity; ++index) {
                    int const i = P * piece_size + static_cast<int>(index);
                    T const t = L == 0 ? T(i) : T(T(2 * i + 1) * h);
                    Hyperbolic<T> const hyperbolic = detail::sinh_cosh(t, constants);
                    T const q = detail::exp(T(-T(hyperbolic.s + hyperbolic.s)), constants);
                    if (!(q > T(0))) break;
                    T const denominator = T(T(1) + q);
                    nodes.comp[index] = T(T(T(2) * q) / denominator);
                    nodes.w[index] = T(T(T(T(4) * hyperbolic.c) * q) / T(denominator * denominator));
                    nodes.size = index + 1;
                }
                return nodes;
            });
        }

        /**
         * @brief A piece of the tanh-sinh table for float, double and long double, computed at compile time once per
         * program, so that these types pay nothing at run time; each piece is its own initializer, within the
         * constexpr budget.
         * @tparam T float, double or long double.
         * @tparam L The level.
         * @tparam P The piece.
         */
        template<typename T, int L, int P>
        inline constexpr Nodes<T, piece_capacity<T>(L, P)> tanh_sinh_table = detail::tanh_sinh_piece<T, L, P>();

        /**
         * @brief @p use applied to a piece of the tanh-sinh table in T, so that the caller need not know where the
         * piece comes from: the compile-time table for float, double and long double, with_generated() (keyed per
         * level, piece and T) for every other T.
         * @param use A callable taking the piece.
         * @return What use returns.
         */
        template<real T, int L, int P, typename U>
        constexpr decltype(auto) with_tanh_sinh_piece(U &&use)
        {
            if constexpr (tabled<T>) return use(detail::tanh_sinh_table<T, L, P>);
            else return detail::with_generated<&detail::tanh_sinh_piece<T, L, P>>(use);
        }

        /**
         * @brief A view of one level of the tanh-sinh table across its pieces, so that a step reads its level once for
         * both sides and the slivers' weights (PHASE3B_API.md C).
         *
         * Node i is t = i at level 0 and t = (2i + 1) 2^-level above, t ascending. The view holds spans into the
         * pieces, valid only inside with_tanh_sinh_level()'s callback: in constant evaluation a piece of a type
         * without a table exists only for the call that generated it.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct TanhSinhLevel
        {
            /** @brief The most pieces of any level: the finest level's, the size of the arrays of spans. */
            static constexpr std::size_t most =
                static_cast<std::size_t>(detail::level_pieces<T>(tanh_sinh_levels - 1));

            /** @brief The complements of each piece, its nodes held. */
            std::array<std::span<T const>, most> comp;

            /** @brief The weights of each piece. */
            std::array<std::span<T const>, most> w;

            /** @brief The nodes of the level. */
            int size;

            /**
             * @brief @p use applied to piece P of level L, for with_level(): how the generic assembly names this
             * table's pieces.
             * @param use A callable taking the piece.
             * @return What use returns.
             */
            template<int L, int P, typename U>
            static constexpr decltype(auto) with_piece(U &&use)
            { return detail::with_tanh_sinh_piece<T, L, P>(use); }

            /**
             * @brief Adds the nodes of piece @p p, for with_level().
             * @param p The piece, after those added.
             * @param nodes Its nodes.
             * @return Whether the piece is full, so that the level may go on: a piece that is not full holds the end of
             *         the table.
             */
            template<std::size_t N>
            constexpr bool add(std::size_t p, Nodes<T, N> const &nodes)
            {
                comp[p] = std::span<T const>(nodes.comp.data(), nodes.size);
                w[p] = std::span<T const>(nodes.w.data(), nodes.size);
                size += static_cast<int>(nodes.size);
                return nodes.size == N;
            }

            /**
             * @brief The complement 1 - x of node @p i, by its index in the level, for the steps that read single
             * nodes.
             * @param i The node, 0 <= i < size.
             * @return Its complement.
             */
            constexpr T const &complement(int i) const
            {
                assert(0 <= i && i < size);
                return comp[static_cast<std::size_t>(i / piece_size)][static_cast<std::size_t>(i % piece_size)];
            }

            /**
             * @brief The weight dx/dt of node @p i, by its index in the level.
             * @param i The node, 0 <= i < size.
             * @return Its weight.
             */
            constexpr T const &weight(int i) const
            {
                assert(0 <= i && i < size);
                return w[static_cast<std::size_t>(i / piece_size)][static_cast<std::size_t>(i % piece_size)];
            }

            /**
             * @brief @p visit(i, comp, w) on the nodes from @p first, t ascending, until it returns false: the sweep of
             * a side, a walk through the pieces in order, cheaper than complement() and weight() per node in constant
             * evaluation.
             * @param first The first node, at least 0.
             * @param visit The visitor.
             */
            template<typename V>
            constexpr void each(int first, V &&visit) const
            {
                // Plain loops, not std::ranges::all_of over iota: this is the hot path of constant evaluation, where
                // the algorithm took 14 % more steps per evaluation and 1/x to the cap 98 % of Clang's default limit.
                std::size_t const from = static_cast<std::size_t>(first);
                for (std::size_t p = from / piece_size; p < most; ++p) {
                    std::size_t const start = p * piece_size;
                    for (std::size_t index = from > start ? from - start : 0; index < comp[p].size(); ++index)
                        if (!visit(static_cast<int>(start + index), comp[p][index], w[p][index])) return;
                }
            }
        };

        /**
         * @brief What the exp-sinh and sinh-sinh maps read per node, for a piece of a level, t ascending: one table for
         * both infinite maps (PHASE4_API.md B, Q5).
         *
         * With c = 3/2 and u = c sinh t: on [a, inf) the nodes are a + E toward infinity and a + q toward a, on the
         * line +-x. The weights are formed at use from these columns: f E, f q, and (f E + f q) / 2 on the line.
         * c = 3/2 rather than 1 or pi/2: it needs no pi, and it measured fewer evaluations and fewer false successes
         * than c = 1 (PHASE4_API.md Q2).
         *
         * The two sides end at different t: the nodes toward infinity are a prefix of those toward the finite end, so
         * their columns hold the first `outward` entries, the others the first `inward`.
         *
         * @tparam T The real type.
         * @tparam N The capacity.
         */
        template<real T, std::size_t N>
        struct ExpSinhNodes
        {
            /** @brief q = exp(-u), the distance toward the finite end. */
            Column<T, N> q;

            /** @brief f = c cosh t, the factor of every weight. */
            Column<T, N> f;

            /** @brief E = 1 / q, the distance toward infinity. */
            Column<T, N> e;

            /**
             * @brief The line's node x = sinh u: by its series where u < 1/2, where E - q cancels, else (E - q) / 2;
             * stored rather than formed at use, so that the line costs one load per node and no branch.
             */
            Column<T, N> x;

            /** @brief The nodes toward the finite end: N, or fewer where the table ends. */
            std::size_t inward;

            /** @brief The nodes toward infinity and on the line, at most inward. */
            std::size_t outward;
        };

        /**
         * @brief Generates a piece of a level of the exp-sinh table, from exp() on its domain alone; the one generator
         * for both the compile-time and the run-time tables (PHASE4_API.md B).
         *
         * c = 3/2 is formed as T(3) / T(2), exact in even radices (2 and 10 among them). Levels and positions are the
         * tanh-sinh table's. The side toward the finite end ends at the first node whose q = exp(-u) is 0, below x_min:
         * q >= min() for the nodes kept. The side toward infinity ends at the first node with f > (max / 2) q, decided
         * before E = 1/q is formed, so that 1/q cannot overflow and E and f E stay at most about max / 2 in every T.
         *
         * @tparam T The real type.
         * @tparam L The level, 0 to 9.
         * @tparam P The piece.
         * @return The piece.
         */
        template<real T, int L, int P>
        constexpr ExpSinhNodes<T, detail::piece_capacity<T>(L, P)> exp_sinh_piece()
        {
            constexpr std::size_t capacity = detail::piece_capacity<T>(L, P);
            return detail::with_generated<&detail::exp_constants<T>>([](ExpConstants<T> const &constants) {
                ExpSinhNodes<T, capacity> nodes{ detail::column<T, capacity>(), detail::column<T, capacity>(),
                    detail::column<T, capacity>(), detail::column<T, capacity>(), 0, 0 };
                T const c = T(T(3) / T(2));
                T const half = T(T(1) / T(2));
                T const smallest = T((std::numeric_limits<T>::min)());
                T const half_max = T(T((std::numeric_limits<T>::max)()) / T(2));
                T h = T(1);
                for (int i = 0; i < L; ++i) h = T(h / T(2));
                bool outward = true;
                for (std::size_t index = 0; index < capacity; ++index) {
                    int const i = P * piece_size + static_cast<int>(index);
                    T const t = L == 0 ? T(i) : T(T(2 * i + 1) * h);
                    Hyperbolic<T> const hyperbolic = detail::sinh_cosh(t, constants);
                    T const u = T(c * hyperbolic.s);
                    T const q = detail::exp(T(-u), constants);
                    if (q < smallest) break;
                    T const f = T(c * hyperbolic.c);
                    nodes.q[index] = q;
                    nodes.f[index] = f;
                    nodes.inward = index + 1;
                    // Once a node is refused toward infinity, so is every later one: the side ends at its first.
                    outward = outward && f <= T(half_max * q);
                    if (!outward) continue;
                    T const e = T(T(1) / q);
                    nodes.e[index] = e;
                    nodes.x[index] = u < half ? detail::sinh_cosh(u, constants).s : T(T(e - q) / T(2));
                    nodes.outward = index + 1;
                }
                return nodes;
            });
        }

        /**
         * @brief A piece of the exp-sinh table for float, double and long double, computed at compile time once per
         * program, each piece its own initializer, as tanh_sinh_table.
         * @tparam T float, double or long double.
         * @tparam L The level.
         * @tparam P The piece.
         */
        template<typename T, int L, int P>
        inline constexpr ExpSinhNodes<T, piece_capacity<T>(L, P)> exp_sinh_table = detail::exp_sinh_piece<T, L, P>();

        /**
         * @brief @p use applied to a piece of the exp-sinh table in T, wherever the piece comes from, as
         * with_tanh_sinh_piece().
         * @param use A callable taking the piece.
         * @return What use returns.
         */
        template<real T, int L, int P, typename U>
        constexpr decltype(auto) with_exp_sinh_piece(U &&use)
        {
            if constexpr (tabled<T>) return use(detail::exp_sinh_table<T, L, P>);
            else return detail::with_generated<&detail::exp_sinh_piece<T, L, P>>(use);
        }

        /**
         * @brief A view of one level of the exp-sinh table across its pieces, so that a step of the infinite maps reads
         * its level once, as TanhSinhLevel does for the finite one.
         *
         * Positions as TanhSinhLevel's: q and f for the nodes toward the finite end, E and x for the first `outward` of
         * them. Valid only inside with_exp_sinh_level()'s callback, for the same reason as TanhSinhLevel.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct ExpSinhLevel
        {
            /** @brief The most pieces of any level: the finest level's, the size of the array of pieces. */
            static constexpr std::size_t most =
                static_cast<std::size_t>(detail::level_pieces<T>(tanh_sinh_levels - 1));

            /** @brief The columns of a piece, each of its nodes held on its side. */
            struct Piece
            {
                /** @brief q, toward the finite end. */
                std::span<T const> q;

                /** @brief f, toward the finite end. */
                std::span<T const> f;

                /** @brief E, toward infinity. */
                std::span<T const> e;

                /** @brief x, toward infinity. */
                std::span<T const> x;
            };

            /** @brief The pieces. */
            std::array<Piece, most> pieces;

            /** @brief The nodes of the level toward the finite end. */
            int inward;

            /** @brief The nodes of the level toward infinity and on the line, at most inward. */
            int outward;

            /**
             * @brief @p use applied to piece P of level L, for with_level(): how the generic assembly names this
             * table's pieces.
             * @param use A callable taking the piece.
             * @return What use returns.
             */
            template<int L, int P, typename U>
            static constexpr decltype(auto) with_piece(U &&use)
            { return detail::with_exp_sinh_piece<T, L, P>(use); }

            /**
             * @brief Adds the nodes of piece @p p, for with_level().
             *
             * The side toward infinity of a piece is taken only while every earlier piece's was full, so that the side
             * ends at its first refused node in t.
             *
             * @param p The piece, after those added.
             * @param nodes Its nodes.
             * @return Whether the piece is full toward the finite end, so that the level may go on.
             */
            template<std::size_t N>
            constexpr bool add(std::size_t p, ExpSinhNodes<T, N> const &nodes)
            {
                std::size_t const held = static_cast<std::size_t>(outward) == p * piece_size ? nodes.outward : 0;
                pieces[p] = Piece{ std::span<T const>(nodes.q.data(), nodes.inward),
                    std::span<T const>(nodes.f.data(), nodes.inward), std::span<T const>(nodes.e.data(), held),
                    std::span<T const>(nodes.x.data(), held) };
                inward += static_cast<int>(nodes.inward);
                outward += static_cast<int>(held);
                return nodes.inward == N;
            }

            /**
             * @brief The piece of node @p i, for the accessors below.
             * @param i The node, at least 0.
             * @return Its piece.
             */
            constexpr Piece const &piece_of(int i) const
            { return pieces[static_cast<std::size_t>(i / piece_size)]; }

            /**
             * @brief The index of node @p i in its piece, for the accessors below.
             * @param i The node, at least 0.
             * @return Its index.
             */
            static constexpr std::size_t index_of(int i)
            { return static_cast<std::size_t>(i % piece_size); }

            // One accessor per column, each indexing its own: a pointer to member shared by the four crashed cl 19.51
            // in constant evaluation (design.md, "The exp-sinh table").

            /**
             * @brief q = exp(-u) of node @p i, by its index in the level.
             * @param i The node, 0 <= i < inward.
             * @return Its q.
             */
            constexpr T const &q(int i) const
            {
                assert(0 <= i && i < inward);
                return piece_of(i).q[index_of(i)];
            }

            /**
             * @brief f = c cosh t of node @p i, by its index in the level.
             * @param i The node, 0 <= i < inward.
             * @return Its f.
             */
            constexpr T const &f(int i) const
            {
                assert(0 <= i && i < inward);
                return piece_of(i).f[index_of(i)];
            }

            /**
             * @brief E = 1/q of node @p i, by its index in the level.
             * @param i The node, 0 <= i < outward.
             * @return Its E.
             */
            constexpr T const &e(int i) const
            {
                assert(0 <= i && i < outward);
                return piece_of(i).e[index_of(i)];
            }

            /**
             * @brief The line's node x = sinh u of node @p i, by its index in the level.
             * @param i The node, 0 <= i < outward.
             * @return Its x.
             */
            constexpr T const &x(int i) const
            {
                assert(0 <= i && i < outward);
                return piece_of(i).x[index_of(i)];
            }

            /**
             * @brief @p visit(i, piece, index) on the nodes of one side from @p first, t ascending, until it returns
             * false: the sweep of a side, a walk through the pieces in order, as TanhSinhLevel::each.
             *
             * It passes the piece's columns and the node's index in the piece, rather than values, so that a sweep
             * reads only the columns its side needs.
             *
             * @param first The first node, at least 0.
             * @param infinite Whether the side is toward infinity (or the line's): the first `outward` nodes; else
             *        the `inward` nodes toward the finite end.
             * @param visit The visitor.
             */
            template<typename V>
            constexpr void each(int first, bool infinite, V &&visit) const
            {
                // Plain loops, as TanhSinhLevel::each's: the hot path of constant evaluation (design.md). A piece
                // after the side's last holds none of its nodes.
                std::size_t const from = static_cast<std::size_t>(first);
                for (std::size_t p = from / piece_size; p < most; ++p) {
                    Piece const &piece = pieces[p];
                    std::size_t const size = infinite ? piece.e.size() : piece.q.size();
                    std::size_t const start = p * piece_size;
                    for (std::size_t index = from > start ? from - start : 0; index < size; ++index)
                        if (!visit(static_cast<int>(start + index), piece, index)) return;
                }
            }
        };

        /**
         * @brief @p use applied to level L of a table once pieces P onward are added to @p level: the assembly of a
         * level from its pieces, shared by both tables (design.md, "The exp-sinh table": a level at once,
         * generalised).
         *
         * Each piece is a distinct template (its level and index are template arguments), and in constant evaluation a
         * generated piece lives only for the callback it is passed to. So the pieces are visited by nested callbacks,
         * one per piece, and use runs innermost, where every piece the view refers to is still alive. A piece that is
         * not full ends the level.
         *
         * @tparam T The real type.
         * @tparam L The level.
         * @tparam P The next piece.
         * @param level The view (TanhSinhLevel<T> or ExpSinhLevel<T>), holding pieces 0 to P - 1.
         * @param use A callable taking the view.
         * @return What use returns.
         */
        template<real T, int L, int P, typename View, typename U>
        constexpr decltype(auto) with_pieces(View &level, U &use)
        {
            if constexpr (P == detail::level_pieces<T>(L)) return use(std::as_const(level));
            else {
                return View::template with_piece<L, P>([&](auto const &nodes) -> decltype(auto) {
                    if (!level.add(P, nodes)) return use(std::as_const(level));
                    return detail::with_pieces<T, L, P + 1>(level, use);
                });
            }
        }

        /**
         * @brief @p use applied to a level of a table in T, every piece of it at once: the entry from a level known at
         * run time to the pieces, which are named at compile time.
         *
         * The level is matched against L = 0, 1, ... by recursion, so that the matching L can name the pieces as
         * template arguments. The pieces are the compile-time tables for float, double and long double,
         * with_generated() for every other T (in constant evaluation, each piece generated as a temporary for the
         * call).
         *
         * @tparam T The real type.
         * @tparam View The view of a level: TanhSinhLevel<T> or ExpSinhLevel<T>.
         * @param level The level, 0 to 9.
         * @param use A callable taking a View const &.
         * @return What use returns.
         */
        template<real T, typename View, int L = 0, typename U>
        constexpr decltype(auto) with_level(int level, U &&use)
        {
            assert(0 <= level && level < tanh_sinh_levels);
            if constexpr (L + 1 < tanh_sinh_levels) {
                if (level != L) return detail::with_level<T, View, L + 1>(level, use);
            }
            View nodes{};
            return detail::with_pieces<T, L, 0>(nodes, use);
        }

        /**
         * @brief @p use applied to a level of the tanh-sinh table in T, as with_level(): what a step of TanhSinh on a
         * finite interval calls.
         * @param level The level, 0 to 9.
         * @param use A callable taking a TanhSinhLevel<T> const &.
         * @return What use returns.
         */
        template<real T, typename U>
        constexpr decltype(auto) with_tanh_sinh_level(int level, U &&use)
        { return detail::with_level<T, TanhSinhLevel<T>>(level, use); }

        /**
         * @brief @p use applied to a level of the exp-sinh table in T, as with_level(): what a step of TanhSinh on a
         * half line or the line calls.
         * @param level The level, 0 to 9.
         * @param use A callable taking an ExpSinhLevel<T> const &.
         * @return What use returns.
         */
        template<real T, typename U>
        constexpr decltype(auto) with_exp_sinh_level(int level, U &&use)
        { return detail::with_level<T, ExpSinhLevel<T>>(level, use); }

        /**
         * @brief The finest level, whose step 2^-9 makes every node's t times 2^9 an integer: positions are kept so,
         * which lets a side compare its kept points across levels exactly, without arithmetic in T.
         */
        inline constexpr int finest_level = tanh_sinh_levels - 1;

        /**
         * @brief The position of node @p i of a level, as t 2^9, the integer form in which a side records its points.
         * @param level The level.
         * @param i The node's index in its level.
         * @return t 2^9.
         */
        constexpr int position(int level, int i)
        { return level == 0 ? i << finest_level : (2 * i + 1) << (finest_level - level); }

        /**
         * @brief The index of the first node of a level beyond the position @p after (t 2^9): where a side's sweep at a
         * new level starts, past the outermost point it kept.
         * @param level The level.
         * @param after The position, at least 0.
         * @return The index.
         */
        constexpr int first_beyond(int level, int after)
        {
            if (level == 0) return (after >> finest_level) + 1;
            return ((after >> (finest_level - level)) + 1) / 2;
        }

        /**
         * @brief How a side of the tanh-sinh sum ended, which decides its end term: the bound on what lies beyond the
         * side's last kept point (PHASE3_API.md E, PHASE3B_API.md C, PHASE4_API.md D).
         *
         * A side may change kind as levels refine it: the table's end becomes any kind; rounding becomes rounding,
         * pruning or a tail; pruning and a tail become each other, never rounding, since their cut node was strictly
         * inside and so is every node before it.
         */
        enum class cut_by {
            table,    ///< at the table's end, no cut: geometric end term (PHASE3B_API.md D1)
            rounding, ///< a point not strictly inside the interval: the sliver and the power-law end term
            value,    ///< a value of fn not finite in the tail: geometric end term
            pruning   ///< a term below eps L1 beyond the kept ones: end term 0
        };

        /**
         * @brief A point of a side whose term entered the sums, kept so that the side's end term can be formed from
         * its last points without storing them all.
         * @tparam T The real type.
         */
        template<real T>
        struct Kept
        {
            /**
             * @brief Its position, t 2^9: an integer on the finest level's grid, so that nodes of every level compare
             * exactly (detail::position).
             */
            int position;

            /** @brief The point, for its distance to the end in the power-law end term. */
            T point;

            /** @brief fn's value there, for the sliver and the power-law end term. */
            T value;

            /**
             * @brief |term| without the step h, for the geometric end term, which applies the current h itself, and
             * for the pruning threshold.
             */
            T size;
        };

        /**
         * @brief A side of the tanh-sinh sum, towards b (right) or a (left), the centre included in both: what the
         * method must remember of a side between levels so that each level evaluates only its new nodes and still
         * knows where the side ends and how to bound what lies beyond (PHASE3B_API.md C).
         *
         * It holds the cut and its kind, the outermost kept point and the two others its end term needs, and W, the
         * slivers' weights. For exp-sinh and sinh-sinh the sides are as mapping says; a side towards infinity is never
         * cut by rounding, so its W and distinct point go unused.
         *
         * The points of a side are monotone in t (adjacent nodes differ in q by far more than a rounding), so the
         * distance to the end does not grow with t. Two points then stand for the whole side: the kept point with the
         * largest t below the outermost (the geometric end term's Tp: at the current spacing no node lies between
         * them), and the kept point with the largest t farther from the end than the outermost (the power law's second
         * point). Each is updated as points are kept, which makes the same choice as a search of every kept point at
         * each level, in constant time and with no storage per point.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct Side
        {
            /**
             * @brief The cut, t 2^9; -1 for none: nodes at or beyond it are not sampled, so later levels stop the
             * side's sweep there.
             */
            int cut;

            /** @brief How the side ended, which decides its end term. */
            cut_by cause;

            /** @brief The outermost kept point (position 0: the centre): Te of the end terms. */
            Kept<T> outer;

            /** @brief The kept point with the largest t below the outermost (Tp); absent while that is the centre. */
            std::optional<Kept<T>> below;

            /**
             * @brief The kept point with the largest t farther from the end than the outermost: the power law's second
             * point, which needs a distinct distance (points near the end can round to the same one).
             */
            std::optional<Kept<T>> distinct;

            /**
             * @brief W: the table's weights beyond the outermost kept node, every level so far, so that a side cut by
             * rounding can restore the sliver between its last point and the end without calling fn.
             *
             * Added level by level in increasing t (the suffix sum of PHASE3_API.md E, term for term), while the side
             * is uncut or cut by rounding; a side cut by pruning or a tail never returns to rounding, since its cut
             * node was inside, so it never needs W again.
             */
            T weights;

            /**
             * @brief Records a kept point, updating the outermost point and the two the end terms need.
             * @param k The point.
             * @param end The end the side approaches.
             */
            constexpr void keep(Kept<T> const &k, T const &end)
            {
                if (k.position > outer.position) {
                    // The old outermost is the largest t below the new one; it is also the farthest-out point at a
                    // distance above the new outermost's if its own distance is larger (distances only shrink outward).
                    below = outer;
                    if (nxx::math::abs(T(end - outer.point)) > nxx::math::abs(T(end - k.point))) distinct = outer;
                    outer = k;
                    return;
                }
                // A point inside the outermost (a later level fills the gaps): it may be nearer the outermost than the
                // current below, or the farthest-out point still at a distance above the outermost's.
                if (!below || k.position > below->position) below = k;
                if (nxx::math::abs(T(end - k.point)) > nxx::math::abs(T(end - outer.point))
                    && (!distinct || k.position > distinct->position))
                    distinct = k;
            }

            /**
             * @brief Adds level @p level's weights beyond the outermost kept node to W, after the side's sweep, so that
             * W covers every node of every level beyond the last kept point.
             *
             * No node of an earlier level lies between the old and the new outermost node, so adding each level's
             * nodes beyond the new outermost gives the suffix sum term for term.
             *
             * @param nodes The level's nodes.
             * @param level The level.
             */
            constexpr void add_weights(TanhSinhLevel<T> const &nodes, int level)
            {
                // A side cut by pruning or a tail never needs W again (see weights).
                if (cause != cut_by::table && cause != cut_by::rounding) return;
                nodes.each(detail::first_beyond(level, outer.position), [this](int, T const &, T const &w) {
                    weights = T(weights + w);
                    return true;
                });
            }

            /**
             * @brief Adds level @p level's weights f q beyond the outermost kept node to W, for exp-sinh's side
             * toward its finite end, which can be cut by rounding as a tanh-sinh side can; the side toward infinity
             * has no W (PHASE4_API.md D).
             * @param nodes The level's nodes.
             * @param level The level.
             */
            constexpr void add_weights(ExpSinhLevel<T> const &nodes, int level)
            {
                if (cause != cut_by::table && cause != cut_by::rounding) return;
                nodes.each(detail::first_beyond(level, outer.position), false,
                    [this](int, typename ExpSinhLevel<T>::Piece const &piece, std::size_t i) {
                        weights = T(weights + T(piece.f[i] * piece.q[i]));
                        return true;
                    });
            }
        };

        /**
         * @brief A side's end term: the bound on the integral beyond the side's last kept point, which the estimate
         * adds so that a cut side cannot hide the mass it did not sample (PHASE3_API.md E).
         *
         * Unbounded is a flag, not inf: constant evaluation rejects a computed overflow, other compilers reject
         * inf - inf and inf times 0, and nxx::real does not promise an infinity. An unbounded end term leaves the
         * step without a current result, and a run that ends so has no best.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct EndTerm
        {
            /** @brief The bound on what lies beyond the side's end. */
            T value;

            /** @brief Whether no bound exists. */
            bool unbounded;
        };

        /**
         * @brief The end term of a side cut by a value of fn that is not finite, or at the table's end: the bound on
         * the terms beyond the cut, from the last two kept ones (PHASE3_API.md E; PHASE3B_API.md D1 for the table's
         * end).
         *
         * The geometric extrapolation h Te / (1 - r), r = Te / Tp, bounds the tail of a sequence of terms that is
         * log-concave, as x^-alpha's are near a singular end and as decaying tails towards infinity are where they
         * were measured (PHASE4_API.md D, "Soundness"). It is unbounded where the terms do not decrease (1/x: the
         * integral diverges), or where the bound would exceed max / 2. Bailey's h Te alone, which this replaced, gave
         * false successes near singular ends, where r is close to 1.
         *
         * @param h The step.
         * @param outer Te, the outermost kept term without h.
         * @param inner Tp, the next one inward.
         * @return The end term.
         */
        template<real T>
        constexpr EndTerm<T> geometric_end(T const &h, T const &outer, T const &inner)
        {
            if (!(outer < inner)) return EndTerm<T>{ T(0), true };
            T const r = T(outer / inner);
            T const room = T(T(1) - r);
            T const scaled = T(h * outer);
            // h Te / (1 - r) <= max / 2 decided before it is formed (1 - r goes down to an ulp).
            if (scaled > T(room * T(T((std::numeric_limits<T>::max)()) / T(2)))) return EndTerm<T>{ T(0), true };
            return EndTerm<T>{ T(scaled * T(T(1) / room)), false };
        }

        /**
         * @brief The end term of a side cut by rounding: what the sliver between the last kept point and the end may
         * hold beyond the share the sliver term already restored, so that an integrable singularity at the end is not
         * reported as converged too early (PHASE3_API.md E).
         *
         * The sliver term restores W y_e, the rule's share of the sliver at the last kept value; for a function that
         * grows towards the end, the true mass of the sliver is larger. The excess is that of a power law through the
         * last two kept points of distinct distances: with d_e < d_p their distances to the end and y_e, y_p their
         * values, R = d |y|, rho = d_e / d_p and q = R_e / R_p, it is R_e max(0, (1 - rho) q / (rho (1 - q)) - 1):
         * 0 for a constant and an upper bound for x^a. Unbounded if q >= 1 (the mass does not shrink towards the end,
         * as for 1/x), or with no two such points. A first form without it ended (1 - x)^-0.9 ok with an error 8 times
         * its estimate.
         *
         * @param side The side.
         * @param end The end it approaches.
         * @return The end term.
         */
        template<real T>
        constexpr EndTerm<T> power_end(Side<T> const &side, T const &end)
        {
            EndTerm<T> const unbounded{ T(0), true };
            T const half_max = T(T((std::numeric_limits<T>::max)()) / T(2));
            if (!side.distinct) return unbounded;
            // The distances are formed as end - p, exact there; every product is checked before it is formed.
            T const d_e = nxx::math::abs(T(end - side.outer.point));
            T const d_p = nxx::math::abs(T(end - side.distinct->point));
            T const y_e = nxx::math::abs(side.outer.value);
            T const y_p = nxx::math::abs(side.distinct->value);
            if (!detail::small_product(d_e, y_e) || !detail::small_product(d_p, y_p)) return unbounded;
            T const r_e = T(d_e * y_e);
            T const r_p = T(d_p * y_p);
            // A product below min() has lost its relative precision: both 0 made 1/x's excess 0 on a scale near
            // denorm_min, a false success. There q and the excess come from ratios, rho |y_e| / |y_p| and
            // d_e (|y_e| factor); every other case keeps the products.
            T const smallest = T((std::numeric_limits<T>::min)());
            bool const lost = (r_e < smallest && y_e > T(0)) || (r_p < smallest && y_p > T(0));
            // No power law through a zero inner point: fn 0 at both points bounds nothing beyond (0), else unbounded.
            if (lost ? y_p == T(0) : r_p == T(0)) return !lost && r_e == T(0) ? EndTerm<T>{ T(0), false } : unbounded;
            T const rho = T(d_e / d_p);
            T const scaled = lost ? T(rho * y_e) : r_e;
            T const base = lost ? y_p : r_p;
            if (!(scaled < base)) return unbounded;
            T const q = T(scaled / base);
            if (q >= T(1)) return unbounded;
            T const numerator = T(T(T(1) - rho) * q);
            T const denominator = T(rho * T(T(1) - q));
            // The quotient within max / 2, decided before it is formed (1 - q goes down to an ulp).
            if (!(denominator > T(0)) || numerator > T(denominator * half_max)) return unbounded;
            T const ratio = T(T(numerator / denominator) - T(1));
            T const factor = ratio < T(0) ? T(0) : ratio;
            if (lost) {
                if (!detail::small_product(y_e, factor)) return unbounded;
                T const excess = T(y_e * factor);
                if (!detail::small_product(d_e, excess)) return unbounded;
                return EndTerm<T>{ T(d_e * excess), false };
            }
            if (!detail::small_product(r_e, factor)) return unbounded;
            return EndTerm<T>{ T(r_e * factor), false };
        }

        /**
         * @brief TanhSinh's map, chosen by which ends are infinite, so that one method and one state serve finite and
         * infinite intervals (PHASE4_API.md B; design.md, "TanhSinh on infinite intervals").
         *
         * It also fixes the roles of the two sides: the right side is swept first except on (-inf, b], where the side
         * toward -inf is (each map sweeps its side toward infinity first, PHASE4_API.md C).
         */
        enum class mapping {
            interval, ///< tanh-sinh on [a, b]: the right side toward b, the left toward a
            upper,    ///< exp-sinh on [a, inf): a + E toward infinity on the right, a + q toward a on the left
            lower,    ///< exp-sinh on (-inf, b]: b - E toward -inf on the left, b - q toward b on the right
            line      ///< sinh-sinh on (-inf, inf): x on the right, -x on the left, both toward infinity
        };

        /**
         * @brief Whether a side of exp-sinh or sinh-sinh is toward infinity, which decides its rules: its nodes
         * (E or q), its tail rule, the underflow rule, and whether it can be cut by rounding.
         *
         * Both sides on the line, else the right on [a, inf) and the left on (-inf, b].
         *
         * @param map The map, not interval.
         * @param right Whether the side is the right one.
         * @return Whether it is toward infinity.
         */
        constexpr bool toward_infinity(mapping map, bool right)
        { return map == mapping::line || (map == mapping::upper) == right; }

        /**
         * @brief A term w y of exp-sinh or sinh-sinh, formed so that an overflow becomes error::non_finite in constant
         * evaluation and at run time alike (PHASE4_API.md C).
         *
         * There is no half-width: the weight carries the scale and reaches max / 2, so every term is checked; the
         * tanh-sinh shortcut past the check (safe_value, which assumes weights at most 1) does not apply.
         *
         * @param w The weight, finite.
         * @param y fn's value, finite.
         * @return The term, or nothing where it is not finite (at run time) or could exceed max / 2 (in constant
         *         evaluation, where an overflow is not a constant expression): non_finite.
         */
        template<real T>
        constexpr std::optional<T> exp_sinh_term(T const &w, T const &y)
        {
            if consteval {
                if (!detail::small_product(w, y)) return std::nullopt;
                return T(w * y);
            }
            else {
                T const t = T(w * y);
                if (!nxx::math::isfinite(t)) return std::nullopt;
                return t;
            }
        }

        /**
         * @brief TanhSinh's state between levels: everything a level needs from the earlier ones, so that each level
         * evaluates fn only at its new nodes (PHASE3B_API.md C).
         *
         * Each level halves the step and adds the nodes between the earlier ones, so the state keeps the sums of every
         * term so far (R_k, L1_k), the previous reported value for the estimate, and per side where it ends and the
         * few points its end term needs. It holds no table and no reference to fn, is independent of fn's type, of
         * fixed size, and allocates nothing.
         *
         * One state for the three maps (start returns one type): the map, and for exp-sinh and sinh-sinh an interval
         * that holds only the finite end, the origin of the points (0 on the line).
         *
         * @tparam T The real type.
         */
        template<real T>
        struct TanhSinhState
        {
            /**
             * @brief The first level at which a result may be returned, as Romberg's (PHASE3_API.md E), so that an
             * estimate from the early levels, which can fall below the error, does not end the run: level 4 rather
             * than Boost's 3 removed the battery's remaining estimates below the error.
             */
            static constexpr int min_level = 4;

            /** @brief The map; tanh-sinh where the ends are finite. */
            mapping map = mapping::interval;

            /**
             * @brief The interval [a, b], a < b, which places the points and scales the terms; for exp-sinh and
             * sinh-sinh the finite end (0 on the line) as c, lo and hi, with hw and width 0: b - a is never formed
             * there.
             */
            Interval<T> interval;

            /**
             * @brief max / 2, formed once for the constant-evaluation guards of the sums: in constant evaluation every
             * step counts against the compiler's limit.
             */
            T half_max;

            /**
             * @brief The largest |fn| for which a term cannot make Interval::checked_term refuse, so that constant
             * evaluation can skip checked_term's guards on the hot path.
             *
             * max / 8 over max(|hw|, 1), every weight being at most 1; below it, term() alone. 0 for exp-sinh and
             * sinh-sinh, whose terms are always checked.
             */
            T safe_value;

            /** @brief The side towards b (for the other maps, as mapping says). */
            Side<T> right;

            /** @brief The side towards a (for the other maps, as mapping says). */
            Side<T> left;

            /** @brief The step, 2^-k, which scales each new term before it is summed. */
            T h;

            /**
             * @brief R_k, the trapezoidal sum of the kept terms, without the slivers (PHASE3B_API.md D3): the slivers
             * change with each level, so they are added to the reported value only, never carried in the running sum.
             */
            T value;

            /** @brief L1_k, the same sum of |terms|, for the rounding floor and the pruning threshold. */
            T l1;

            /**
             * @brief V_(k-1), the previous level's reported value (with its slivers), for the estimate
             * |V_k - V_(k-1)|; empty at level 0.
             */
            std::optional<T> previous;

            /** @brief The level. */
            int k;

            /** @brief The step: {V_k, E_k} from level 1 unless unbounded, attainable 2 F_k, ready from min_level. */
            Step<T> step;
        };

        /**
         * @brief A term of the tanh-sinh sum, formed so that an overflow becomes error::non_finite in constant
         * evaluation and at run time alike, at the least cost on the hot path.
         *
         * In constant evaluation Interval::checked_term, or term() alone where |y| <= safe_value; at run time term(),
         * refused where it is not finite.
         *
         * @param s The state.
         * @param w The weight.
         * @param y fn's value, finite.
         * @return The term, or empty where it overflows (at run time) or constant evaluation could not form it within
         *         max / 2.
         */
        template<real T>
        constexpr std::optional<T> tanh_sinh_term(TanhSinhState<T> const &s, T const &w, T const &y)
        {
            if consteval {
                // checked_term refuses an overflow here, and term() alone below safe_value cannot overflow.
                if (nxx::math::abs(y) <= s.safe_value) return s.interval.term(w, y);
                return s.interval.checked_term(w, y);
            }
            else {
                T const t = s.interval.term(w, y);
                if (!nxx::math::isfinite(t)) return std::nullopt;
                return t;
            }
        }

        /**
         * @brief x added to @p sum, unless constant evaluation could not keep the sum within max / 2: detail::fits
         * for the sweeps' sums, with max / 2 passed in, formed once, since in constant evaluation every step counts.
         * @param sum The sum.
         * @param x The term.
         * @param half_max max / 2.
         * @return false where refused: non_finite.
         */
        template<real T>
        constexpr bool add_within(Sum<T> &sum, T const &x, T const &half_max)
        {
            if consteval {
                if (nxx::math::abs(sum.sum) > half_max || nxx::math::abs(x) > half_max) return false;
            }
            sum.add(x);
            return true;
        }

        /**
         * @brief The point of a node on a side, formed from the node's complement 1 - x, never from x, so that the
         * points near an end stay distinct where x itself rounds to 1.
         *
         * b - hw comp on the right, a + hw comp on the left: c + hw x with x = 1 - comp, rewritten so that the
         * complement, which the table holds to full relative accuracy, sets the distance to the end.
         *
         * @param s The state.
         * @param right Whether the side is the right one.
         * @param comp The node's complement.
         * @return The point.
         */
        template<real T>
        constexpr T tanh_sinh_point(TanhSinhState<T> const &s, bool right, T const &comp)
        {
            return right ? T(s.interval.hi - T(s.interval.hw * comp)) : T(s.interval.lo + T(s.interval.hw * comp));
        }

        /**
         * @brief Level k's step from its sums, for every map: what the level reports to integrate(), so that one
         * function forms the result, the estimate and the limits for every level and map (PHASE3_API.md E,
         * PHASE3B_API.md D3).
         *
         * It forms the slivers of the sides cut by rounding, the end terms, the reported value V_k = R_k + the
         * slivers, the floor F_k with the slivers' sizes, and the estimate E_k = max(|V_k - V_(k-1)|, F_k) + the end
         * terms. The plain level difference serves as the estimate because successive levels roughly square the
         * error, so the difference bounds the previous level's error. Level 0, which has no previous value, and a level
         * whose end terms are unbounded report no current result.
         *
         * @param s The state, with R_k and L1_k formed.
         * @return Nothing, or non_finite.
         */
        template<real T>
        constexpr std::expected<void, Failure<T>> tanh_sinh_report(TanhSinhState<T> &s)
        {
            Failure<T> const invalid{ error::non_finite, std::nullopt };
            if (!nxx::math::isfinite(s.value) || !nxx::math::isfinite(s.l1)) return std::unexpected(invalid);
            T const floor = detail::rounding_floor(s.l1);
            T sliver = T(0);
            T sliver_size = T(0);
            EndTerm<T> ends{ T(0), false };
            bool refused = false;
            auto const end_of = [&](Side<T> const &side, T const &end) {
                if (side.cause == cut_by::rounding) {
                    // The sliver between the last kept point and the end, narrower than an ulp of the end, keeps its
                    // share of the rule: W at the last kept value (no fn calls); on exp-sinh's side toward its finite
                    // end the term has no half-width.
                    std::optional<T> const piece =
                        s.map == mapping::interval
                            ? s.interval.checked_term(T(s.h * side.weights), side.outer.value)
                            : detail::exp_sinh_term(T(s.h * side.weights), side.outer.value);
                    if (!piece || !detail::fits(sliver, *piece) || !detail::fits(sliver_size, *piece)) {
                        refused = true;
                        return;
                    }
                    sliver = T(sliver + *piece);
                    sliver_size = T(sliver_size + nxx::math::abs(*piece));
                }
                // The end term by how the side ended: rounding, the power-law excess; pruning, 0 (the dropped term is
                // below eps L1, F / 50); a tail or the table's end, the geometric bound, or h Te where only the centre
                // was kept; a last term of 0, 0.
                EndTerm<T> end_term{ T(0), false };
                if (side.cause == cut_by::rounding) end_term = detail::power_end(side, end);
                else if (side.cause == cut_by::pruning || side.outer.size == T(0)) end_term = EndTerm<T>{ T(0), false };
                else if (!side.below) end_term = EndTerm<T>{ T(s.h * side.outer.size), false }; // only the centre kept
                else end_term = detail::geometric_end(s.h, side.outer.size, side.below->size);
                if (end_term.unbounded) ends.unbounded = true;
                else if (!detail::fits(ends.value, end_term.value)) refused = true; // each about max / 2 at most
                else ends.value = T(ends.value + end_term.value);
            };
            end_of(s.right, s.interval.hi);
            end_of(s.left, s.interval.lo);
            if (refused || !detail::fits(s.value, sliver) || !detail::fits(floor, sliver_size))
                return std::unexpected(invalid);
            T const total = T(s.value + sliver);
            T const total_floor = T(floor + sliver_size);
            if (!nxx::math::isfinite(total) || !nxx::math::isfinite(total_floor)) return std::unexpected(invalid);
            if consteval {
                // The attainable error below is 2 F, which must not overflow in a constant expression.
                if (detail::large(total_floor)) return std::unexpected(invalid);
            }
            std::optional<Integral<T>> current;
            // Level 0 has no previous value, so no estimate: its step has no current result.
            if (s.previous) {
                if (!detail::fits(total, T(-*s.previous))) return std::unexpected(invalid);
                T const difference = nxx::math::abs(T(total - *s.previous));
                T const base = detail::larger(difference, total_floor);
                if (!ends.unbounded && !detail::fits(base, ends.value)) return std::unexpected(invalid);
                T const estimate = ends.unbounded ? base : T(base + ends.value);
                if (!nxx::math::isfinite(estimate)) return std::unexpected(invalid);
                if (!ends.unbounded) current = Integral<T>{ total, detail::larger(estimate, total_floor) };
            }
            // Twice the floor (Phase 2's test): E >= F plus the end terms, so F alone could not pass. The cap is the
            // table's finest level: exhausted there.
            s.step = Step<T>{ .current = current,
                .attainable = T(T(2) * total_floor),
                .ready = s.k >= TanhSinhState<T>::min_level,
                .exhausted = s.k == tanh_sinh_levels - 1 };
            s.previous = total;
            return {};
        }

        /**
         * @brief TanhSinh's start on a finite interval: level 0, which samples the centre and each side out to where
         * the side ends, and sets the cuts and sums every later level builds on (PHASE3_API.md C and D, PHASE3B_API.md
         * D2).
         *
         * The centre first, then each side outward at t = 1, 2, ..., right side first, until its first point not
         * strictly inside (a rounding cut), its first non-finite value in the tail (a tail cut), or the table's end.
         * The terms are stored, not summed, because the pruning threshold needs the L1 of the whole sweep: each side
         * then drops its run of outermost terms below eps L1, from the outside in, and V_0 and L1_0 are the sums of
         * the kept terms only. With no point strictly inside, a step exhausted with no current, and fn never called.
         *
         * @param fn The function.
         * @param a The lower end.
         * @param b The upper end, a < b.
         * @return The state after level 0, or non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<TanhSinhState<T>, Failure<T>> tanh_sinh_start(F &fn, T const &a, T const &b)
        {
            using result = std::expected<TanhSinhState<T>, Failure<T>>;
            Failure<T> const invalid{ error::non_finite, std::nullopt };
            Interval<T> const interval = detail::interval_of(a, b);
            T const half_max = T(T((std::numeric_limits<T>::max)()) / T(2));
            T const safe_value =
                T(T(T((std::numeric_limits<T>::max)()) / T(8)) / detail::larger(nxx::math::abs(interval.hw), T(1)));
            Kept<T> const unsampled{ 0, interval.c, T(0), T(0) };
            Side<T> const fresh{ -1, cut_by::table, unsampled, std::nullopt, std::nullopt, T(0) };
            TanhSinhState<T> s{ .interval = interval,
                .half_max = half_max,
                .safe_value = safe_value,
                .right = fresh,
                .left = fresh,
                .h = T(1),
                .value = T(0),
                .l1 = T(0),
                .previous = std::nullopt,
                .k = 0,
                .step = Step<T>{} };
            // No point of T strictly inside: nothing can be sampled, and fn is never called at an end.
            if (!(a < interval.c && interval.c < b)) {
                s.step = Step<T>{ .ready = false, .exhausted = true };
                return s;
            }
            return detail::with_tanh_sinh_level<T>(0, [&](TanhSinhLevel<T> const &nodes) -> result {
                T const eps = T(std::numeric_limits<T>::epsilon());
                T const centre_value = detail::evaluate(fn, interval.c);
                if (!nxx::math::isfinite(centre_value)) return std::unexpected(invalid);
                std::optional<T> const centre_term = interval.checked_term(nodes.weight(0), centre_value);
                if (!centre_term || !nxx::math::isfinite(*centre_term)) return std::unexpected(invalid);
                T const centre_size = nxx::math::abs(*centre_term);
                Kept<T> const middle{ 0, interval.c, centre_value, centre_size };
                s.right.outer = middle;
                s.left.outer = middle;

                struct Sample
                {
                    int position;
                    T point;
                    T value;
                    T term;
                };
                // Level 0's samples of a side in a fixed array (at most level_capacity(0) nodes): no allocation.
                constexpr std::size_t most = static_cast<std::size_t>(detail::level_capacity<T>(0));
                struct Swept
                {
                    std::array<std::optional<Sample>, most> samples{};
                    std::size_t count = 0;

                    constexpr std::span<std::optional<Sample> const> held() const
                    { return std::span(samples).first(count); }
                };
                // A side outward from t = 1; false where fn is not finite outside a tail, or a term overflows or is
                // refused: non_finite.
                auto const sweep = [&](Side<T> &side, bool right, Swept &swept) {
                    bool failed = false;
                    nodes.each(1, [&](int i, T const &comp, T const &w) {
                        int const j = detail::position(0, i);
                        T const p = detail::tanh_sinh_point(s, right, comp);
                        // A point that rounds onto the end, or past it, is never evaluated: fn may be singular there.
                        if (!(a < p && p < b)) {
                            side.cut = j;
                            side.cause = cut_by::rounding;
                            return false;
                        }
                        T const y = detail::evaluate(fn, p);
                        if (!nxx::math::isfinite(y)) {
                            // A tail, the endpoint singularity the method exists for: within eps hw of the end and
                            // beyond every kept point of the side (level 0 sweeps outward, t > 0). Elsewhere a value
                            // that is not finite is non_finite (PHASE3_API.md D).
                            if (comp < eps) {
                                side.cut = j;
                                side.cause = cut_by::value;
                            }
                            else failed = true;
                            return false;
                        }
                        std::optional<T> const t = detail::tanh_sinh_term(s, w, y);
                        if (!t) { // a term that overflows, or is refused: non_finite at once
                            failed = true;
                            return false;
                        }
                        assert(swept.count < most);
                        swept.samples[swept.count++] = Sample{ j, p, y, *t };
                        return true;
                    });
                    return !failed;
                };
                Swept right_swept;
                Swept left_swept;
                if (!sweep(s.right, true, right_swept) || !sweep(s.left, false, left_swept))
                    return std::unexpected(invalid);
                // The pruning threshold: eps times the L1 of the whole sweep, the centre included (h = 1), summed in
                // the order centre, right side, left side.
                Sum<T> sizes;
                auto const add_sizes = [&](Swept const &swept) {
                    return std::ranges::all_of(swept.held(), [&](std::optional<Sample> const &x) {
                        return detail::add_within(sizes, nxx::math::abs(x->term), half_max);
                    });
                };
                bool refused = !detail::add_within(sizes, centre_size, half_max) || !add_sizes(right_swept)
                               || !add_sizes(left_swept);
                std::optional<T> const l1_sweep = detail::compensated(sizes);
                if (refused || !l1_sweep || !nxx::math::isfinite(*l1_sweep)) return std::unexpected(invalid);
                // Each side drops its run of outermost small terms, from the outside in, and is cut at the innermost
                // dropped node; a small term with a larger one beyond it is kept (PHASE3B_API.md D2), since cutting at
                // the first small term going outward gave a false success where one term was exactly 0. A pruned
                // side's end term is 0, below the floor's eps L1 by a factor 50. Pruning saves the evaluations a
                // smooth fn does not need: 2.8 to 3.4 times fewer in double (PHASE3_API.md C).
                T const small = T(eps * *l1_sweep);
                auto const prune = [&](Side<T> &side, Swept &swept) {
                    auto const outward = swept.held() | std::views::reverse;
                    auto const large = std::ranges::find_if(outward, [&](std::optional<Sample> const &x) {
                        return !(nxx::math::abs(x->term) < small);
                    });
                    std::size_t const dropped = static_cast<std::size_t>(std::ranges::distance(outward.begin(), large));
                    if (dropped == 0) return;
                    swept.count -= dropped;
                    side.cut = swept.samples[swept.count]->position;
                    side.cause = cut_by::pruning;
                };
                prune(s.right, right_swept);
                prune(s.left, left_swept);
                // V_0 and L1_0 from the kept terms only: the centre, the right side, the left side; each kept point
                // recorded on its side, and W formed beyond the outermost.
                Sum<T> values;
                Sum<T> magnitudes;
                auto const add_kept = [&](Side<T> &side, Swept const &swept, T const &end) {
                    bool const added = std::ranges::all_of(swept.held(), [&](std::optional<Sample> const &x) {
                        T const size = nxx::math::abs(x->term);
                        side.keep(Kept<T>{ x->position, x->point, x->value, size }, end);
                        return detail::add_within(values, x->term, half_max)
                               && detail::add_within(magnitudes, size, half_max);
                    });
                    side.add_weights(nodes, 0);
                    return added;
                };
                refused = !detail::add_within(values, *centre_term, half_max)
                          || !detail::add_within(magnitudes, centre_size, half_max)
                          || !add_kept(s.right, right_swept, b) || !add_kept(s.left, left_swept, a);
                std::optional<T> const v_0 = detail::compensated(values);
                std::optional<T> const a_0 = detail::compensated(magnitudes);
                if (refused || !v_0 || !a_0) return std::unexpected(invalid);
                s.value = *v_0;
                s.l1 = *a_0;
                return detail::tanh_sinh_report(s).transform([&] { return s; });
            });
        }

        /**
         * @brief TanhSinh's next level k on a finite interval: the step halved, fn evaluated only at the new nodes, so
         * that every earlier evaluation is reused (PHASE3_API.md C, PHASE3B_API.md D2).
         *
         * Per side, from the centre outward, the odd multiples of 2^-k up to its cut. At most one new node lies beyond
         * the side's outermost kept one (PHASE3B_API.md C); it cuts the side where it is not strictly inside
         * (rounding), where fn is not finite there and comp < eps (a tail), or where |h term| < eps L1_(k-1) / 2
         * (pruning). A value that is not finite at any other node is non_finite. Then R_k = R_(k-1) / 2 + the new
         * terms, L1 likewise, W, and the step.
         *
         * @param fn The function.
         * @param s The state, not exhausted.
         * @return Nothing, or non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> tanh_sinh_advance(F &fn, TanhSinhState<T> &s)
        {
            Failure<T> const invalid{ error::non_finite, std::nullopt };
            T const eps = T(std::numeric_limits<T>::epsilon());
            s.k += 1;
            s.h = T(s.h / T(2));
            // The pruning threshold in term scale: |h term| < eps L1_(k-1) / 2, eps times the running L1 at the start
            // of the level (PHASE3B_API.md D2). Known before the level starts, it is the same for both sides and
            // independent of the order of evaluation.
            T const threshold = T(T(eps * s.l1) / T(T(2) * s.h));
            using result = std::expected<void, Failure<T>>;
            return detail::with_tanh_sinh_level<T>(s.k, [&](TanhSinhLevel<T> const &nodes) -> result {
                Sum<T> new_values;
                Sum<T> new_magnitudes;
                // A side from the centre outward; false where fn is not finite outside a tail, or a term overflows or
                // is refused: non_finite.
                auto const sweep = [&](Side<T> &side, bool right) {
                    T const &end = right ? s.interval.hi : s.interval.lo;
                    bool failed = false;
                    bool refused = false;
                    nodes.each(0, [&](int i, T const &comp, T const &w) {
                        int const j = detail::position(s.k, i);
                        // A cut side is never sampled at or beyond its cut again.
                        if (side.cut >= 0 && j >= side.cut) return false;
                        T const p = detail::tanh_sinh_point(s, right, comp);
                        // Only a node beyond every kept point can be a tail or be pruned; a node between kept points
                        // refines the sum and is always kept, and a value there that is not finite is non_finite.
                        bool const outer = j > side.outer.position;
                        if (!(s.interval.lo < p && p < s.interval.hi)) {
                            side.cut = j;
                            side.cause = cut_by::rounding;
                            return false;
                        }
                        T const y = detail::evaluate(fn, p);
                        if (!nxx::math::isfinite(y)) {
                            if (outer && comp < eps) {
                                side.cut = j;
                                side.cause = cut_by::value;
                            }
                            else failed = true;
                            return false;
                        }
                        std::optional<T> const t = detail::tanh_sinh_term(s, w, y);
                        if (!t) { // a term that overflows, or is refused: non_finite at once
                            failed = true;
                            return false;
                        }
                        T const size = nxx::math::abs(*t);
                        if (outer && size < threshold) {
                            side.cut = j;
                            side.cause = cut_by::pruning;
                            return false;
                        }
                        refused = refused || !detail::add_within(new_values, T(s.h * *t), s.half_max);
                        refused = refused || !detail::add_within(new_magnitudes, T(s.h * size), s.half_max);
                        side.keep(Kept<T>{ j, p, y, size }, end);
                        return true;
                    });
                    side.add_weights(nodes, s.k);
                    return !failed && !refused;
                };
                if (!sweep(s.right, true) || !sweep(s.left, false)) return std::unexpected(invalid);
                std::optional<T> const added = detail::compensated(new_values);
                std::optional<T> const added_l1 = detail::compensated(new_magnitudes);
                if (!added || !added_l1) return std::unexpected(invalid);
                // Halving the step halves every earlier term: R_k = R_(k-1) / 2 + the new terms, each already
                // multiplied by the new h.
                T const half_value = T(s.value / T(2));
                T const half_l1 = T(s.l1 / T(2));
                if (!detail::fits(half_value, *added) || !detail::fits(half_l1, *added_l1))
                    return std::unexpected(invalid);
                s.value = T(half_value + *added);
                s.l1 = T(half_l1 + *added_l1);
                return detail::tanh_sinh_report(s);
            });
        }

        /**
         * @brief A node of exp-sinh or sinh-sinh on a side, formed at use from the shared table's columns, so that one
         * table serves both maps and both kinds of side (PHASE4_API.md B).
         *
         * Its distance from the origin, as the tail rule reads it, its point and its weight.
         *
         * @tparam T The real type.
         */
        template<real T>
        struct ExpSinhNode
        {
            /** @brief E toward infinity, q toward the finite end, x on the line. */
            T distance;

            /** @brief a + E, a + q; b - E, b - q; x, -x. */
            T point;

            /** @brief f E, f q, or (f E + f q) / 2 on the line. */
            T weight;
        };

        /**
         * @brief Node @p i of @p piece on a side of exp-sinh or sinh-sinh: the one place the maps' points and weights
         * are formed from the table.
         *
         * The line's weight is dx/dt = f cosh u = (f E + f q) / 2, from the columns the half lines use. On (-inf, b]
         * the point is b - d, formed directly: the same value as the prototype's reflection, fn(-x) on [-b, inf),
         * since round to nearest is symmetric (design.md, "TanhSinh on infinite intervals").
         *
         * @param s The state, its map not interval.
         * @param right Whether the side is the right one.
         * @param infinite Whether it is toward infinity.
         * @param piece The node's piece.
         * @param i The node's index in it.
         * @return The node.
         */
        template<real T>
        constexpr ExpSinhNode<T> exp_sinh_node(TanhSinhState<T> const &s, bool right, bool infinite,
            typename ExpSinhLevel<T>::Piece const &piece, std::size_t i)
        {
            T const &f = piece.f[i];
            if (s.map == mapping::line) {
                T const &x = piece.x[i];
                return ExpSinhNode<T>{ x, right ? x : T(-x), T(T(T(f * piece.e[i]) + T(f * piece.q[i])) / T(2)) };
            }
            T const &d = infinite ? piece.e[i] : piece.q[i];
            T const &origin = s.interval.lo;
            return ExpSinhNode<T>{ d, s.map == mapping::upper ? T(origin + d) : T(origin - d), T(f * d) };
        }

        /**
         * @brief Whether a point of exp-sinh or sinh-sinh is strictly inside, so that fn is never called at the finite
         * end or at an infinity.
         *
         * Beyond the finite end, and finite (the line's +-x, at most max / 6, always is). Toward the finite end a
         * point fails by rounding onto the end; toward infinity only by overflowing, where the finite end is large.
         *
         * @param s The state, its map not interval.
         * @param p The point.
         * @return Whether fn may be evaluated there.
         */
        template<real T>
        constexpr bool exp_sinh_inside(TanhSinhState<T> const &s, T const &p)
        {
            if (s.map == mapping::upper) return s.interval.lo < p && nxx::math::isfinite(p);
            if (s.map == mapping::lower) return p < s.interval.lo && nxx::math::isfinite(p);
            return true;
        }

        /**
         * @brief 1/eps, the distance beyond which the tail rule toward infinity accepts a value that is not finite,
         * where T can hold it; empty otherwise, so that a type with a small max() neither overflows nor misreads the
         * rule.
         *
         * Present where eps max() > 1, tested without overflow. In a type whose exponents end below 1/eps no finite
         * node lies beyond it, so that clause of the rule never holds (Phase 4's independent review).
         *
         * @return 1/eps, or empty.
         */
        template<real T>
        constexpr std::optional<T> inverse_epsilon()
        {
            T const eps = T(std::numeric_limits<T>::epsilon());
            if (!(T(eps * T((std::numeric_limits<T>::max)())) > T(1))) return std::nullopt;
            return T(T(1) / eps);
        }

        /**
         * @brief TanhSinh's start on an infinite interval: level 0 of exp-sinh or sinh-sinh, which samples the centre
         * and each side out to where it ends, and sets the cuts and sums every later level builds on (PHASE4_API.md C
         * and D).
         *
         * The maps: exp-sinh on [a, inf), x = a + E toward infinity and a + q toward a, or (-inf, b], x = b - E and
         * b - q; sinh-sinh on the line, x = +-sinh u. Level 0 evaluates the centre (a + 1, b - 1 or 0, weight
         * c = 3/2), then each side outward at t = 1, 2, ..., the side toward infinity first (on the line the positive
         * side), each node's term stored. A point not strictly inside ends the side. A value not finite is a tail
         * where its distance is beyond 1/eps (toward infinity) or below eps (toward the finite end), is pending
         * toward infinity, else non_finite. The rest waits for the threshold, eps L1 of the centre and every stored
         * term: each pending value is a tail where the side's last stored |term| is below it, else non_finite; then
         * the underflow rule toward infinity in increasing t; then Phase 3's pruning from the outside in. With no
         * point strictly inside (a + 1 rounds onto a), a step exhausted with no current, and fn never called.
         *
         * @param fn The function.
         * @param a The lower end, -inf or finite.
         * @param b The upper end, inf or finite, a < b, not both finite.
         * @return The state after level 0, or non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<TanhSinhState<T>, Failure<T>> exp_sinh_start(F &fn, T const &a, T const &b)
        {
            using result = std::expected<TanhSinhState<T>, Failure<T>>;
            Failure<T> const invalid{ error::non_finite, std::nullopt };
            bool const upward = !nxx::math::isfinite(b);
            mapping const map = upward && !nxx::math::isfinite(a) ? mapping::line
                                                                   : (upward ? mapping::upper : mapping::lower);
            T const origin = map == mapping::line ? T(0) : (upward ? a : b);
            T const half_max = T(T((std::numeric_limits<T>::max)()) / T(2));
            Kept<T> const unsampled{ 0, origin, T(0), T(0) };
            Side<T> const fresh{ -1, cut_by::table, unsampled, std::nullopt, std::nullopt, T(0) };
            TanhSinhState<T> s{ .map = map,
                .interval = Interval<T>{ origin, T(0), T(0), true, false, origin, origin },
                .half_max = half_max,
                .safe_value = T(0),
                .right = fresh,
                .left = fresh,
                .h = T(1),
                .value = T(0),
                .l1 = T(0),
                .previous = std::nullopt,
                .k = 0,
                .step = Step<T>{} };
            // The centre a + 1 (b - 1) is formed only where it cannot overflow, at most max - 1 from 0 (in double
            // max - 1 is max itself, so this excludes nothing there); beyond, in a type whose max is small against
            // 1 / eps, a + 1 would be infinite (Phase 4's independent review), and nothing is sampled either.
            T const room = T(T((std::numeric_limits<T>::max)()) - T(1));
            bool const fits =
                map == mapping::upper ? !(room < origin) : (map == mapping::lower ? !(origin < T(-room)) : true);
            T const centre = !fits ? origin
                                   : (map == mapping::upper ? T(origin + T(1))
                                                            : (map == mapping::lower ? T(origin - T(1)) : T(0)));
            // No point strictly inside (from |a| about 2 / eps): nothing can be sampled.
            bool const inside = fits && nxx::math::isfinite(centre)
                                && (map == mapping::upper ? origin < centre
                                                          : (map == mapping::lower ? centre < origin : true));
            if (!inside) {
                s.step = Step<T>{ .ready = false, .exhausted = true };
                return s;
            }
            if consteval {
                // origin +- E stays within max where |origin| <= max / 2, since E <= max / 3 (f E <= max / 2 and
                // f >= 3/2); beyond it a + 1 rounds onto a in any type with max above 4 / eps, so this refuses
                // nothing that run time would sample.
                if (detail::large(origin)) return std::unexpected(invalid);
            }
            // The right side is swept first, except on (-inf, b], where the side toward -inf is the left.
            bool const left_first = map == mapping::lower;
            return detail::with_exp_sinh_level<T>(0, [&](ExpSinhLevel<T> const &nodes) -> result {
                T const eps = T(std::numeric_limits<T>::epsilon());
                std::optional<T> const inverse_eps = detail::inverse_epsilon<T>();
                T const smallest = T((std::numeric_limits<T>::min)());
                T const centre_value = detail::evaluate(fn, centre);
                if (!nxx::math::isfinite(centre_value)) return std::unexpected(invalid);
                std::optional<T> const centre_term = detail::exp_sinh_term(nodes.f(0), centre_value); // f(0) = c
                if (!centre_term) return std::unexpected(invalid);
                T const centre_size = nxx::math::abs(*centre_term);
                Kept<T> const middle{ 0, centre, centre_value, centre_size };
                s.right.outer = middle;
                s.left.outer = middle;

                struct Sample
                {
                    int position;
                    T point;
                    T value;
                    T term;
                    T weight;
                };
                // Level 0's samples of a side in a fixed array, as tanh_sinh_start's, and the position of a value not
                // finite toward infinity whose decision waits for the threshold.
                constexpr std::size_t most = static_cast<std::size_t>(detail::level_capacity<T>(0));
                struct Swept
                {
                    std::array<std::optional<Sample>, most> samples{};
                    std::size_t count = 0;
                    std::optional<int> pending{};

                    constexpr std::span<std::optional<Sample> const> held() const
                    { return std::span(samples).first(count); }
                };
                // A side outward from t = 1; false where fn is not finite outside a tail and not pending, or a term
                // is refused: non_finite.
                auto const sweep = [&](Side<T> &side, bool right, Swept &swept) {
                    bool const infinite = detail::toward_infinity(map, right);
                    bool failed = false;
                    nodes.each(1, infinite,
                        [&](int i, typename ExpSinhLevel<T>::Piece const &piece, std::size_t index) {
                            int const j = detail::position(0, i);
                            ExpSinhNode<T> const node = detail::exp_sinh_node(s, right, infinite, piece, index);
                            // Toward infinity a point fails only by overflowing, which ends the side as the table's
                            // end would; toward the finite end it rounds onto the end, a cut by rounding.
                            if (!detail::exp_sinh_inside(s, node.point)) {
                                side.cut = j;
                                side.cause = infinite ? cut_by::table : cut_by::rounding;
                                return false;
                            }
                            T const y = detail::evaluate(fn, node.point);
                            if (!nxx::math::isfinite(y)) {
                                // The tail rule's first clause, by distance; toward infinity its second clause needs
                                // the threshold, known only after both sweeps, so the value waits as pending.
                                if (infinite ? inverse_eps && node.distance > *inverse_eps : node.distance < eps) {
                                    side.cut = j;
                                    side.cause = cut_by::value;
                                }
                                else if (infinite) swept.pending = j;
                                else failed = true;
                                return false;
                            }
                            std::optional<T> const t = detail::exp_sinh_term(node.weight, y);
                            if (!t) {
                                failed = true;
                                return false;
                            }
                            assert(swept.count < most);
                            swept.samples[swept.count++] = Sample{ j, node.point, y, *t, node.weight };
                            return true;
                        });
                    return !failed;
                };
                Side<T> &first = left_first ? s.left : s.right;
                Side<T> &second = left_first ? s.right : s.left;
                Swept first_swept;
                Swept second_swept;
                if (!sweep(first, !left_first, first_swept) || !sweep(second, left_first, second_swept))
                    return std::unexpected(invalid);
                // The threshold: eps times the L1 of the whole sweep, the centre included (h = 1), summed in the order
                // centre, first side, second side.
                Sum<T> sizes;
                auto const add_sizes = [&](Swept const &swept) {
                    return std::ranges::all_of(swept.held(), [&](std::optional<Sample> const &x) {
                        return detail::add_within(sizes, nxx::math::abs(x->term), half_max);
                    });
                };
                bool refused = !detail::add_within(sizes, centre_size, half_max) || !add_sizes(first_swept)
                               || !add_sizes(second_swept);
                std::optional<T> const l1_sweep = detail::compensated(sizes);
                if (refused || !l1_sweep || !nxx::math::isfinite(*l1_sweep)) return std::unexpected(invalid);
                T const small = T(eps * *l1_sweep);
                // The tail rule's second clause (PHASE4_API.md Q3): a pending value is a tail where the side's last
                // stored term is already below the threshold.
                auto const decide = [&](Side<T> &side, Swept const &swept) {
                    if (!swept.pending) return true;
                    T const last = swept.count > 0 ? nxx::math::abs(swept.samples[swept.count - 1]->term) : centre_size;
                    if (!(last < small)) return false;
                    side.cut = *swept.pending;
                    side.cause = cut_by::value;
                    return true;
                };
                if (!decide(first, first_swept) || !decide(second, second_swept)) return std::unexpected(invalid);
                // The underflow rule toward infinity (PHASE4_API.md Q4), in increasing t: the first stored node with
                // |y| < min(), w min() at least the threshold, after a value of at least min() (the centre's for the
                // first), cuts the side as a table's end, dropped with every node beyond it. Toward infinity the
                // weights reach max / 2, so a value that underflowed is no evidence of a negligible term; where w min()
                // is below the threshold the term is negligible anyway, and pruning below drops it.
                auto const truncate = [&](Side<T> &side, bool right, Swept &swept) {
                    if (!detail::toward_infinity(map, right)) return;
                    auto const indices = std::views::iota(std::size_t{ 0 }, swept.count);
                    auto const at = std::ranges::find_if(indices, [&](std::size_t i) {
                        Sample const &x = *swept.samples[i];
                        T const &before = i > 0 ? swept.samples[i - 1]->value : centre_value;
                        return nxx::math::abs(x.value) < smallest && T(x.weight * smallest) >= small
                               && nxx::math::abs(before) >= smallest;
                    });
                    if (at == indices.end()) return;
                    side.cut = swept.samples[*at]->position;
                    side.cause = cut_by::table;
                    swept.count = *at;
                };
                truncate(first, !left_first, first_swept);
                truncate(second, left_first, second_swept);
                // Each side drops its run of outermost small terms, from the outside in (Phase 3).
                auto const prune = [&](Side<T> &side, Swept &swept) {
                    auto const outward = swept.held() | std::views::reverse;
                    auto const large = std::ranges::find_if(outward, [&](std::optional<Sample> const &x) {
                        return !(nxx::math::abs(x->term) < small);
                    });
                    std::size_t const dropped = static_cast<std::size_t>(std::ranges::distance(outward.begin(), large));
                    if (dropped == 0) return;
                    swept.count -= dropped;
                    side.cut = swept.samples[swept.count]->position;
                    side.cause = cut_by::pruning;
                };
                prune(first, first_swept);
                prune(second, second_swept);
                // V_0 and L1_0 from the kept terms only: the centre, the first side, the second side; W for the side
                // toward the finite end.
                Sum<T> values;
                Sum<T> magnitudes;
                auto const add_kept = [&](Side<T> &side, bool right, Swept const &swept) {
                    bool const added = std::ranges::all_of(swept.held(), [&](std::optional<Sample> const &x) {
                        T const size = nxx::math::abs(x->term);
                        side.keep(Kept<T>{ x->position, x->point, x->value, size }, origin);
                        return detail::add_within(values, x->term, half_max)
                               && detail::add_within(magnitudes, size, half_max);
                    });
                    if (!detail::toward_infinity(map, right)) side.add_weights(nodes, 0);
                    return added;
                };
                refused = !detail::add_within(values, *centre_term, half_max)
                          || !detail::add_within(magnitudes, centre_size, half_max)
                          || !add_kept(first, !left_first, first_swept) || !add_kept(second, left_first, second_swept);
                std::optional<T> const v_0 = detail::compensated(values);
                std::optional<T> const a_0 = detail::compensated(magnitudes);
                if (refused || !v_0 || !a_0) return std::unexpected(invalid);
                s.value = *v_0;
                s.l1 = *a_0;
                return detail::tanh_sinh_report(s).transform([&] { return s; });
            });
        }

        /**
         * @brief TanhSinh's next level k on an infinite interval: the step halved, fn evaluated only at the new nodes,
         * with the tail and underflow rules of the side toward infinity (PHASE4_API.md D).
         *
         * Per side, in the order of the start, the odd multiples of 2^-k up to its cut. A new node beyond the
         * outermost kept one cuts the side, in this order:
         * - where it is not strictly inside;
         * - where fn is not finite there and its distance is beyond 1/eps (toward infinity) or below eps (toward the
         *   finite end), or toward infinity where the outermost kept |term| is below the threshold eps L1_(k-1) / (2h)
         *   (a tail);
         * - toward infinity, where |y| < min() after a kept value of at least min() (pruning if w min() is below the
         *   threshold, else a table's end; the node not kept);
         * - where its |term| is below the threshold (pruning).
         *
         * Elsewhere a value not finite is non_finite. Then the sums and the step, as tanh_sinh_advance.
         *
         * @param fn The function.
         * @param s The state, not exhausted, its map not interval.
         * @return Nothing, or non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> exp_sinh_advance(F &fn, TanhSinhState<T> &s)
        {
            Failure<T> const invalid{ error::non_finite, std::nullopt };
            T const eps = T(std::numeric_limits<T>::epsilon());
            std::optional<T> const inverse_eps = detail::inverse_epsilon<T>();
            T const smallest = T((std::numeric_limits<T>::min)());
            s.k += 1;
            s.h = T(s.h / T(2));
            T const threshold = T(T(eps * s.l1) / T(T(2) * s.h));
            using result = std::expected<void, Failure<T>>;
            return detail::with_exp_sinh_level<T>(s.k, [&](ExpSinhLevel<T> const &nodes) -> result {
                Sum<T> new_values;
                Sum<T> new_magnitudes;
                auto const sweep = [&](Side<T> &side, bool right) {
                    bool const infinite = detail::toward_infinity(s.map, right);
                    bool failed = false;
                    bool refused = false;
                    nodes.each(0, infinite,
                        [&](int i, typename ExpSinhLevel<T>::Piece const &piece, std::size_t index) {
                            int const j = detail::position(s.k, i);
                            if (side.cut >= 0 && j >= side.cut) return false;
                            ExpSinhNode<T> const node = detail::exp_sinh_node(s, right, infinite, piece, index);
                            // As in tanh_sinh_advance, only a node beyond every kept point can cut the side.
                            bool const outer = j > side.outer.position;
                            if (!detail::exp_sinh_inside(s, node.point)) {
                                side.cut = j;
                                side.cause = infinite ? cut_by::table : cut_by::rounding;
                                return false;
                            }
                            T const y = detail::evaluate(fn, node.point);
                            if (!nxx::math::isfinite(y)) {
                                // The tail rule: by distance, or toward infinity after a term already negligible
                                // (PHASE4_API.md Q3); the threshold is known here, so nothing waits.
                                bool const far =
                                    infinite ? inverse_eps && node.distance > *inverse_eps : node.distance < eps;
                                if (outer && (far || (infinite && side.outer.size < threshold))) {
                                    side.cut = j;
                                    side.cause = cut_by::value;
                                }
                                else failed = true;
                                return false;
                            }
                            std::optional<T> const t = detail::exp_sinh_term(node.weight, y);
                            if (!t) {
                                failed = true;
                                return false;
                            }
                            // The underflow rule (PHASE4_API.md Q4), as at level 0: a value below min() is no evidence
                            // of a negligible term where the weight is large.
                            if (outer && infinite && nxx::math::abs(y) < smallest
                                && !(nxx::math::abs(side.outer.value) < smallest)) {
                                side.cut = j;
                                side.cause = T(node.weight * smallest) < threshold ? cut_by::pruning : cut_by::table;
                                return false;
                            }
                            T const size = nxx::math::abs(*t);
                            if (outer && size < threshold) {
                                side.cut = j;
                                side.cause = cut_by::pruning;
                                return false;
                            }
                            refused = refused || !detail::add_within(new_values, T(s.h * *t), s.half_max);
                            refused = refused || !detail::add_within(new_magnitudes, T(s.h * size), s.half_max);
                            // The end is the finite end (0 on the line): distances matter only on the side toward it.
                            side.keep(Kept<T>{ j, node.point, y, size }, s.interval.lo);
                            return true;
                        });
                    // W only toward the finite end: a side toward infinity is never cut by rounding.
                    if (!infinite) side.add_weights(nodes, s.k);
                    return !failed && !refused;
                };
                bool const left_first = s.map == mapping::lower;
                Side<T> &first = left_first ? s.left : s.right;
                Side<T> &second = left_first ? s.right : s.left;
                if (!sweep(first, !left_first) || !sweep(second, left_first)) return std::unexpected(invalid);
                std::optional<T> const added = detail::compensated(new_values);
                std::optional<T> const added_l1 = detail::compensated(new_magnitudes);
                if (!added || !added_l1) return std::unexpected(invalid);
                // R_k = R_(k-1) / 2 + the new terms, as in tanh_sinh_advance.
                T const half_value = T(s.value / T(2));
                T const half_l1 = T(s.l1 / T(2));
                if (!detail::fits(half_value, *added) || !detail::fits(half_l1, *added_l1))
                    return std::unexpected(invalid);
                s.value = T(half_value + *added);
                s.l1 = T(half_l1 + *added_l1);
                return detail::tanh_sinh_report(s);
            });
        }

    } // namespace detail
    /** @endcond */

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief Whether R is std::expected<S, Failure<T>> for some state S, so that iterative_method can require that
         * shape of start()'s result without naming the method's state type.
         * @tparam R The type.
         * @tparam T The real type.
         */
        template<typename R, typename T>
        concept expected_state = requires { typename R::value_type; }
                                 && std::same_as<R, std::expected<typename R::value_type, Failure<T>>>;

        /**
         * @brief How an Observer stores its callable: an lvalue by reference, so that a recorder the caller owns sees
         * every step; an rvalue by value, so that a temporary lambda lives as long as the Observer.
         * @tparam O The callable's type as given to observer().
         */
        template<typename O>
        using kept = std::conditional_t<std::is_lvalue_reference_v<O>,
            std::reference_wrapper<std::remove_reference_t<O>>, std::remove_cvref_t<O>>;

        /**
         * @brief The type an iterative method's start returns for an fn of type F, from which state_t takes the state:
         * std::expected<S, Failure<T>> for its state S.
         * @tparam M The method.
         * @tparam T The real type.
         * @tparam F fn's type as the method receives it.
         */
        template<typename M, typename T, typename F>
        using start_result = decltype(std::declval<M const &>().start(
            std::declval<F &>(), std::declval<T const &>(), std::declval<T const &>()));

        /**
         * @brief fn's type as integrate() passes it to a method, so that the concepts are checked against what a
         * method actually receives: a function as a pointer to it, any other callable as itself (reference removed).
         * @tparam F fn's type as integrate() received it.
         */
        template<typename F>
        using passed = std::conditional_t<std::is_function_v<std::remove_reference_t<F>>,
            std::add_pointer_t<std::remove_reference_t<F>>, std::remove_reference_t<F>>;

        /**
         * @brief fn in the form integrate() passes to a method, so that a plain function works in constant evaluation
         * on every compiler.
         *
         * A function is passed as a pointer to it, so that no method's lambda captures a reference to a function
         * (MSVC's constant evaluator reads such a capture as uninitialized); any other callable as the lvalue it is.
         *
         * @param fn The function.
         * @return A pointer to fn, or fn.
         */
        template<typename F>
        constexpr decltype(auto) as_object(F &fn)
        {
            if constexpr (std::is_function_v<F>) return &fn;
            else return (fn);
        }

        /**
         * @brief Whether M has opted in to infinite ends (its static member `infinite_ends` is true), so that
         * integrate() refuses an infinite end with non_finite for every other method before the method could form
         * b - a (PHASE4_API.md Q1).
         * @tparam M The method.
         */
        template<typename M>
        concept takes_infinite_ends = requires {
            { M::infinite_ends } -> std::convertible_to<bool>;
        } && bool(M::infinite_ends);

        /**
         * @brief Whether @p x is NaN, so that integrate() can refuse a NaN end while letting an infinite one through to
         * a method that takes it: nxx::math has no isnan, and isfinite does not tell the two apart.
         * @param x A value.
         * @return x != x.
         */
        template<real T>
        constexpr bool is_nan(T const &x)
        { return !(x == x); } // NOLINT(misc-redundant-expression)

    } // namespace detail
    /** @endcond */

    /**
     * @brief What a type must provide to work with integrate() as a fixed rule: a rule of your own, such as a midpoint
     * rule with a number of panels, that returns the integral without an error estimate.
     *
     * A fixed rule M has a const member `m.apply(fn, a, b)` returning `std::expected<T, error>`, the integral of fn
     * from a to b. GaussLegendre<N> is one. integrate() validates and orders the ends, so the rule handles only the
     * computation.
     *
     * The contract, which integrate() relies on:
     * - **Given:** a < b, both finite, possibly extreme (b - a may overflow; the ends may be subnormal, adjacent, or of
     *   opposite signs so that b - a rounds). integrate() orders the ends and negates the result for b < a. fn is an
     *   lvalue (a function arrives as a pointer to it), valid for the call only: call it with std::invoke and an
     *   argument of type T, and convert its result with T(...).
     * - **Required:** fn evaluated only at points of [a, b] (clamp a point formed as a + h t; form a midpoint as
     *   a / 2 + b / 2); a finite value or an error (integrate() turns a value that is not finite into
     *   error::non_finite); the method's parameters are its own preconditions; exceptions pass through.
     *
     * A type that is both a fixed rule and a method with an estimate is outside the contract. Check a method where it
     * is defined: `static_assert(nxx::quadrature::fixed_rule<MyRule, double>)`; a failure there names the cause (a
     * member that is not const, a return type that is not exactly std::expected<T, error>).
     *
     * @tparam M The method.
     * @tparam T The real type.
     * @tparam F fn's type as a method receives it (a function as a pointer to it); by default a pointer to a plain
     *           function of T, as integrate() passes one.
     */
    template<typename M, typename T, typename F = T (*)(T const &)>
    concept fixed_rule = requires(M const &m, F &fn, T const &a, T const &b) {
        { m.apply(fn, a, b) } -> std::same_as<std::expected<T, error>>;
    };

    /**
     * @brief What a type must provide to work with integrate() as a method with an error estimate that does its whole
     * computation in one call: for a method that does not step, such as a recursion or a combinator around another
     * method (a substitution, a split of the interval).
     *
     * An estimating method M has a const member `m.estimate(fn, a, b, criteria)` returning
     * `std::expected<Integral<T>, Failure<T>>`. A method that improves its estimate step by step is better written as
     * an iterative_method: integrate() then owns its loop and stop test, and observer() can watch its steps.
     *
     * The contract: as for fixed_rule, and
     * - **Given:** criteria with both members finite and not negative, rel_tol at most 1, abs_tol unclamped (it enters
     *   only a maximum), possibly both 0.
     * - **Required:** it terminates for every criteria it can be given, both tolerances 0 included (a cap, or a floor
     *   of its own); an Integral with a finite value and a finite, non-negative abs_error, or a Failure: non_finite
     *   without a best, or not_converged with or without one (integrate() turns anything else into non_finite);
     *   abs_error estimates |value - exact| and covers the rounding of the method's own sum; success means its own stop
     *   test against max(abs_tol, rel_tol |value|), or its floor, was met.
     * - **Infinite ends, opt-in:** a method that declares `static constexpr bool infinite_ends = true;` (TanhSinh does)
     *   may also receive a = -inf, b = +inf or both, always a < b after integrate()'s ordering; b - a and a/2 + b/2 are
     *   then not finite: it must not form them, nor call fn at an infinite end. Without the member, or with it false,
     *   integrate() returns non_finite for an infinite end and the method never sees one. A combinator around a method
     *   must forward the member to pass the inner method's opt-in on (observer() does). The member must be a constant
     *   expression: a static member that is not constexpr fails to compile rather than reading as false.
     *
     * @tparam M The method.
     * @tparam T The real type.
     * @tparam F fn's type, reference removed; by default a pointer to a plain function of T, as integrate() passes one.
     */
    template<typename M, typename T, typename F = T (*)(T const &)>
    concept estimating_method = requires(M const &m, F &fn, T const &a, T const &b, Criteria<T> const &criteria) {
        { m.estimate(fn, a, b, criteria) } -> std::same_as<std::expected<Integral<T>, Failure<T>>>;
    };

    /**
     * @brief The state type of an iterative method M, for a combinator that keeps an inner method's state in its own:
     * what `m.start(fn, a, b)` returns inside its std::expected.
     * @tparam M The method.
     * @tparam T The real type.
     * @tparam F fn's type, reference removed; by default a pointer to a plain function of T, as integrate() passes one.
     */
    template<typename M, typename T, typename F = T (*)(T const &)>
    using state_t = typename detail::start_result<M, T, F>::value_type;

    /**
     * @brief What a type must provide to work with integrate() as a method that improves its estimate step by step:
     * the method computes and reports, integrate() runs the loop and decides when to stop, and observer() can watch
     * the steps.
     *
     * GaussKronrod<15>, TanhSinh and Romberg are iterative methods; a user's type can be another. Splitting the work
     * this way means a method of one's own needs no loop, no stop test and no knowledge of the criteria, and is stopped
     * by the same test as the library's.
     *
     * Three const members, each returning `std::expected<..., Failure<T>>`:
     * - `m.start(fn, a, b)` returns a state S of the method's own type, movable, whose member `step` is the Step<T>
     *   after the start;
     * - `m.advance(fn, s)` takes one more step (a level, a bisection) and updates `s.step`;
     * - `m.settle(s)` refines `s.step` before integrate() decides on it (Gauss-Kronrod forms its sums afresh); a method
     *   with nothing to refine returns `{}`.
     *
     * integrate() runs the loop: it checks each step (a current result as integrate() checks results, an attainable
     * error finite and not negative; else non_finite), tests `ready` and
     * `abs_error <= max(abs_tol, rel_tol |value|, attainable)`, settles before any decision, and advances while the
     * step is not exhausted. The method never sees the criteria.
     *
     * The contract: as for fixed_rule, and
     * - **Given:** advance is never called on an exhausted step; nothing is called after a Failure; settle may be
     *   called twice between steps.
     * - **Required:** the state keeps no reference, pointer or copy of fn (fn is passed to every call and need not be
     *   the same object each time); start and every advance leave `step` describing the state; the method reports
     *   `exhausted` within finitely many steps (one that never does loops forever, and in constant evaluation reaches
     *   the compiler's limit); settle refines `step` without taking a step, leaves `ready` and `exhausted` as they
     *   are, and calling it twice gives the step calling it once gives; `attainable` is an error the method can
     *   certify.
     * - **Infinite ends, opt-in:** as for estimating_method.
     *
     * A type that is both an iterative method and a method of another kind is outside the contract. Check a method
     * where it is defined: `static_assert(nxx::quadrature::iterative_method<MyMethod, double>)`; a failure there names
     * the cause (a member that is not const, settle missing, a return type, the state's `step` member).
     *
     * @tparam M The method.
     * @tparam T The real type.
     * @tparam F fn's type, reference removed; by default a pointer to a plain function of T, as integrate() passes one.
     */
    template<typename M, typename T, typename F = T (*)(T const &)>
    concept iterative_method = requires(M const &m, F &fn, T const &a, T const &b) {
        { m.start(fn, a, b) } -> detail::expected_state<T>;
    } && requires(M const &m, F &fn, state_t<M, T, F> &s) {
        { m.advance(fn, s) } -> std::same_as<std::expected<void, Failure<T>>>;
        { m.settle(s) } -> std::same_as<std::expected<void, Failure<T>>>;
        { std::as_const(s).step } -> std::same_as<Step<T> const &>;
    };

    /**
     * @brief Marks a callable as the observer of an integration, so that integrate() tells it from fn by its type and
     * never by trying to call an argument: made by observer(), passed to integrate() before fn.
     * @tparam O The callable's type as given to observer(): an lvalue reference, kept by reference, or a value.
     */
    template<typename O>
    struct Observer
    {
        /** @brief The callable: a std::reference_wrapper for an lvalue, the callable itself for an rvalue. */
        detail::kept<O> callable;
    };

    /**
     * @brief Makes @p o the observer of an integration by an iterative method, to watch it converge, record its steps
     * or find out why it failed: `o(step)` is called with the Step<T> after the method's start and after each step.
     *
     * The step is for the ends in order (a < b), so for b < a it has the opposite sign to the result; nothing is
     * reported after settle, whose values the result shows. The observed run computes exactly what the unobserved run
     * computes. o must be callable as `o(Step<T> const &)` and return void (a static_assert says so where it is not);
     * it is called as a non-const callable, so a mutable lambda can count; an lvalue is kept by reference, so a
     * recorder keeps what it records; its exceptions pass through, which is also the way out of a method that never
     * ends. Count the evaluations in fn, and the steps in o.
     *
     * @code
     * namespace qd = nxx::quadrature;
     * auto const print = [](qd::Step<double> const &s) {
     *     if (s.current) std::println("{:.15f} +- {:.3g}", s.current->value, s.current->abs_error);
     * };
     * auto const r = qd::integrate<qd::Romberg>(qd::observer(print), fn, 0.0, 1.0, {.rel_tol = 1e-12});
     * @endcode
     *
     * @tparam O The callable's type.
     * @param o A callable taking Step<T> const &, returning void.
     * @return The observer, for integrate().
     */
    template<typename O>
    constexpr Observer<O> observer(O &&o)
    { return Observer<O>{ detail::kept<O>(std::forward<O>(o)) }; }

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief The method integrate() runs when it is given an observer: M wrapped so that each step is reported,
         * itself an iterative method with M's state, so that the observed run is the same computation and integrate()
         * needs no second loop (design.md, "Observing").
         * @tparam M The method.
         * @tparam O The observer's callable, as given to observer().
         */
        template<typename M, typename O>
        struct Observed
        {
            /**
             * @brief The method, held for the call to integrate() that observes it, so that a method that cannot be
             * copied can be observed.
             */
            M const &method;

            /** @brief The observer, called as a non-const callable: a mutable lambda can count its calls. */
            mutable kept<O> observer;

            /** @brief The method's opt-in to infinite ends, forwarded, so that an observed TanhSinh takes them. */
            static constexpr bool infinite_ends = takes_infinite_ends<M>;

            // Self = M makes each declaration dependent, so that an observer on a method without these members fails
            // the concept softly instead of where the class is instantiated.

            /**
             * @brief The method's start, then a report of its first step.
             * @tparam T The real type.
             * @tparam F fn's type.
             * @tparam Self The method, as a dependent type.
             * @param fn The function.
             * @param a The lower end.
             * @param b The upper end.
             * @return The method's state, or its Failure.
             */
            template<real T, typename F, typename Self = M>
            requires iterative_method<Self, T, F>
            constexpr std::expected<state_t<Self, T, F>, Failure<T>> start(F &fn, T const &a, T const &b) const
            {
                auto started = method.start(fn, a, b);
                if (started) report(started->step);
                return started;
            }

            /**
             * @brief The method's step, then a report of it.
             * @tparam F fn's type.
             * @tparam S The method's state.
             * @tparam Self The method, as a dependent type.
             * @param fn The function.
             * @param s The state.
             * @return Nothing, or the method's Failure.
             */
            template<typename F, typename S, typename Self = M>
            constexpr auto advance(F &fn, S &s) const -> decltype(std::declval<Self const &>().advance(fn, s))
            {
                auto const r = method.advance(fn, s);
                if (r) report(s.step);
                return r;
            }

            /**
             * @brief The method's settle, not reported: it refines the step integrate() decides on, and the result
             * shows its values.
             * @tparam S The method's state.
             * @tparam Self The method, as a dependent type.
             * @param s The state.
             * @return Nothing, or the method's Failure.
             */
            template<typename S, typename Self = M>
            constexpr auto settle(S &s) const -> decltype(std::declval<Self const &>().settle(s))
            { return method.settle(s); }

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
                constexpr bool callable = std::invocable<kept<O> &, Step<T> const &>;
                static_assert(callable, "nxx::quadrature::observer(o): o must be callable as o(Step<T> const &)");
                if constexpr (callable) {
                    constexpr bool returns_void = std::is_void_v<std::invoke_result_t<kept<O> &, Step<T> const &>>;
                    static_assert(returns_void,
                        "nxx::quadrature::observer(o): o must return void (a result is reserved for stopping the run)");
                    if constexpr (returns_void) std::invoke(observer, step);
                }
            }
        };

        /**
         * @brief The iterative method m observed by o, for the integrate() forms that take an observer.
         * @param m The method, held by reference: it must outlive the result (a temporary does, for the full
         * expression that integrates with it).
         * @param o The observer.
         * @return The observed method.
         */
        template<typename M, typename O>
        constexpr Observed<M, O> observed(M const &m, Observer<O> &&o)
        { return Observed<M, O>{ m, std::move(o.callable) }; }

        /**
         * @brief Whether an Integral is one a caller may be given (a finite value, a finite, non-negative estimate), so
         * that integrate() turns a broken result of any method, a user's included, into non_finite (design.md,
         * "Checked results").
         * @param r The result.
         * @return true if it is.
         */
        template<real T>
        constexpr bool sound(Integral<T> const &r)
        {
            return nxx::math::isfinite(r.value) && nxx::math::isfinite(r.abs_error) && !(r.abs_error < T(0));
        }

        /**
         * @brief Whether a step can be tested against the request: its current result, if any, sound, and its
         * attainable error finite and not negative.
         *
         * Checked before any arithmetic on the step, which in constant evaluation a NaN would make ill-formed; a step
         * that fails it ends the integration with non_finite.
         *
         * @param p The step.
         * @return true if it can.
         */
        template<real T>
        constexpr bool testable(Step<T> const &p)
        {
            return (!p.current || detail::sound(*p.current)) && nxx::math::isfinite(p.attainable)
                   && !(p.attainable < T(0));
        }

        /**
         * @brief The stop test, written once for every iterative method (BRIEF.md D5): whether a step is ready and has
         * a current result whose estimate is at most max(abs_tol, rel_tol |value|, attainable).
         * @param p The step, testable.
         * @param criteria The validated criteria.
         * @return true if it does.
         */
        template<real T>
        constexpr bool meets(Step<T> const &p, Criteria<T> const &criteria)
        {
            return p.ready && p.current
                   && p.current->abs_error
                          <= detail::wanted(criteria.abs_tol, criteria.rel_tol, p.current->value, p.attainable);
        }

        /**
         * @brief The loop of every iterative method, on a < b: it owns the decision to stop, so that each method only
         * computes and reports its steps (PHASE3A_API.md B, C).
         *
         * Each step is checked, then tested against the request. On a pass, or where the method is exhausted, the step
         * is settled and the decision taken on the settled step; otherwise the method takes one more step. A step that
         * passes before settle and fails after it goes on from the settled step (Gauss-Kronrod's running sums can drift
         * from the sums formed afresh).
         *
         * @param m The method.
         * @param fn The function.
         * @param a The lower end.
         * @param b The upper end.
         * @param criteria The validated criteria.
         * @return The settled result; not_converged with the step's result as the best (or none); or a Failure.
         */
        template<typename M, real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            drive(M const &m, F &fn, T const &a, T const &b, Criteria<T> const &criteria)
        {
            Failure<T> const invalid{ error::non_finite, std::nullopt };
            auto started = m.start(fn, a, b);
            if (!started) return std::unexpected(started.error());
            auto &s = *started;
            while (true) {
                if (!detail::testable(s.step)) return std::unexpected(invalid);
                bool const passed = detail::meets(s.step, criteria);
                if (passed || s.step.exhausted) {
                    // Settle before any decision, and only then: Gauss-Kronrod checks its running sums against the
                    // sums formed afresh before every decision to stop (PHASE3A_API.md, "settle").
                    std::expected<void, Failure<T>> const settled = m.settle(s);
                    if (!settled) return std::unexpected(settled.error());
                    if (!detail::testable(s.step)) return std::unexpected(invalid);
                    if (detail::meets(s.step, criteria)) {
                        // A method that knows its estimate is unreliable may not succeed, whatever the estimate says.
                        if (!s.step.trusted) return std::unexpected(Failure<T>{ error::not_converged, s.step.current });
                        return *s.step.current;
                    }
                    // Exhausted and short of the request: the method can do no more.
                    if (!passed) return std::unexpected(Failure<T>{ error::not_converged, s.step.current });
                    // The pass was the running sums' drift: the settled step is not exhausted, so go on from it.
                    continue;
                }
                std::expected<void, Failure<T>> const advanced = m.advance(fn, s);
                if (!advanced) return std::unexpected(advanced.error());
            }
        }

        /**
         * @brief A method with an estimate on a < b, so that integrate()'s validation and ordering are written once for
         * both kinds: an iterative method's loop, or an estimating method's estimate.
         * @param m The method.
         * @param fn The function, as passed.
         * @param a The lower end.
         * @param b The upper end.
         * @param criteria The validated criteria.
         * @return The method's result.
         */
        template<typename M, real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            estimate_on(M const &m, F &fn, T const &a, T const &b, Criteria<T> const &criteria)
        {
            if constexpr (iterative_method<M, T, F>) return detail::drive(m, fn, a, b, criteria);
            else return m.estimate(fn, a, b, criteria);
        }

        /**
         * @brief The integral from b to a, given the one from a to b: integrate() runs every method on ordered ends and
         * negates for b < a, so that integrate(fn, b, a) is exactly -integrate(fn, a, b).
         * @param r The integral from a to b.
         * @return -r, with the same estimate.
         */
        template<real T>
        constexpr Integral<T> negated(Integral<T> r)
        {
            r.value = T(-r.value);
            return r;
        }

    } // namespace detail
    /** @endcond */

    /**
     * @brief The N-point Gauss-Legendre rule, for an integral that must be a smooth function of the parameters inside
     * fn (an integral inside a root solve, an optimiser or another integral), or for an fn too costly for an adaptive
     * method: exactly N evaluations, no error estimate and no tolerance.
     *
     * An adaptive method chooses its points from fn, so its result jumps as a parameter of fn changes which panels it
     * bisects; a root finder or an optimiser on such a function can stall (Lyness 1983). A fixed rule samples the same
     * points whatever fn, so its result is as smooth as fn. Without an estimate, check the accuracy by comparing two
     * values of N, or one result with a method that estimates its error.
     *
     * The rule evaluates fn at the zeros of the Legendre polynomial P_N mapped to (a, b), and is exact for polynomials
     * of degree 2N - 1, up to the rounding of its nodes and weights and of the sum. Its error on [a, b] is
     * (b - a)^(2N+1) (N!)^4 / ((2N + 1) ((2N)!)^3) f^(2N)(xi) for some xi in (a, b).
     *
     * The nodes and weights come from the three-term recurrence, + - * / in T: for float, double and long double at
     * compile time, once per program; for every other type, such as a Boost.Multiprecision number, at its first use,
     * and then kept (a function-local static, initialised once by the language; a type that changed its precision at
     * run time would see the nodes of its first).
     *
     * An infinite end is error::non_finite; use TanhSinh.
     *
     * @tparam N The number of points, 1 to 64.
     */
    template<int N>
    struct GaussLegendre
    {
        /**
         * @brief Applies the rule to fn on [a, b], for integrate(), which validates and orders the ends first
         * (fixed_rule).
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function, evaluated at the N interior nodes.
         * @param a The lower end, finite.
         * @param b The upper end, finite, a < b.
         * @return The integral, or error::non_finite if a value of fn is not finite or the sum overflowed.
         */
        template<real T, typename F>
        requires(N >= 1 && N <= 64)
        constexpr std::expected<T, error> apply(F &fn, T const &a, T const &b) const
        {
            // Qualified, so that argument-dependent lookup cannot pick a function from the namespace of T or F.
            detail::Interval<T> const interval = detail::interval_of(a, b);
            return detail::with_gauss_legendre<N, T>(
                [&](detail::Rule<T, N> const &rule) { return detail::apply(rule, fn, interval); });
        }
    };

    /**
     * @brief Romberg's method, to see and teach how extrapolating trapezoidal sums converges, or to compare with the
     * other methods; for results, prefer the default, GaussKronrod<15>.
     *
     * It needs more evaluations than a Gauss rule for the same accuracy, and its estimate overstates the error on
     * smooth fn. An iterative method, one step per level: level k is the trapezoidal rule at 2^k + 1 equally spaced
     * points, ends included, extrapolated by Richardson (Romberg 1955). Its estimate at level k is the difference of
     * the last two extrapolated values, |R(k, k) - R(k-1, k-1)|, which estimates the error of the previous value. It
     * may stop from level 4 (17 evaluations) and certifies down to its rounding floor. The cap is 16 levels, 65,537
     * evaluations, then error::not_converged with the best estimate; fn with a kink or a singular derivative reach it.
     *
     * Since it samples the ends, a singularity at an end is error::non_finite; use TanhSinh. An infinite end is
     * error::non_finite; use TanhSinh.
     */
    struct Romberg
    {
        /**
         * @brief Begins the integration, for integrate(): T_0 from the ends, then level 1 (iterative_method).
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param a The lower end, finite.
         * @param b The upper end, finite, a < b.
         * @return The state after level 1, or error::non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<detail::RombergState<T>, Failure<T>> start(F &fn, T const &a, T const &b) const
        { return detail::romberg_start(fn, a, b); }

        /**
         * @brief Takes the next level: the step halved, the new midpoints added, the row extrapolated
         * (iterative_method).
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param s The state, not exhausted.
         * @return Nothing, or error::non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> advance(F &fn, detail::RombergState<T> &s) const
        { return detail::romberg_advance(fn, s); }

        /**
         * @brief Does nothing: Romberg's step is final as reported (iterative_method).
         * @tparam T The real type.
         * @param s The state.
         * @return Nothing.
         */
        template<real T>
        constexpr std::expected<void, Failure<T>> settle([[maybe_unused]] detail::RombergState<T> &s) const
        { return {}; }
    };

    /**
     * @brief The default method, adaptive Gauss-Kronrod (G7-K15), for most integrals on a finite interval: it
     * subdivides where its estimate is largest, so it copes with kinks, steps, peaks and mild endpoint singularities
     * without being told where they are.
     *
     * Choose another method for an endpoint singularity stronger than about (x - a)^-0.65, or many digits in a
     * multiprecision type (TanhSinh); for an infinite end (TanhSinh, alone or on the infinite part of a split
     * interval); or for an fn that must cost fewer than its 150 evaluations (GaussLegendre<N>, which has no estimate).
     *
     * An iterative method, one step per bisection. The interval starts as 10 equal panels, 150 evaluations of fn; then
     * the panel with the largest estimate is bisected, 30 evaluations a step, until the sum of the estimates meets the
     * criteria. A panel's estimate is |K15 - G7| on its 15 points, or its floor if larger: it estimates the 7-point
     * rule's error, not the returned 15-point result's, so it usually overstates the error by orders of magnitude, but
     * it is not a bound.
     *
     * What it misses: a feature narrower than the spacing of the points (a narrow peak between them), and a kink or a
     * step near a panel's ends, outside the outermost points, which neither rule samples; split the interval at a known
     * feature. An endpoint singularity stronger than about (x - a)^-0.65 is underestimated by a ratio that bisection
     * does not change; substitute x = a + u^k with k (1 + alpha) >= 1 for (x - a)^alpha (and u = b - x for the other
     * end), or use TanhSinh. Full precision in a 50-digit type is costly: the estimate converges at the 7-point rule's
     * order. An infinite end is error::non_finite.
     *
     * @tparam N The number of Kronrod points: only 15 exists.
     */
    template<int N>
    requires(N == 15)
    struct GaussKronrod
    {
        /**
         * @brief Begins the integration, for integrate(): the initial mesh of 10 equal panels (iterative_method).
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param a The lower end, finite.
         * @param b The upper end, finite, a < b.
         * @return The state after the mesh, or error::non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<detail::KronrodState<T>, Failure<T>> start(F &fn, T const &a, T const &b) const
        {
            return detail::with_gauss_kronrod<T>(
                [&](detail::Kronrod<T> const &rule) { return detail::kronrod_start(rule, fn, a, b); });
        }

        /**
         * @brief Bisects the panel with the largest estimate, so that the work goes where fn is hardest
         * (iterative_method).
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param s The state, not exhausted.
         * @return Nothing, or error::non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> advance(F &fn, detail::KronrodState<T> &s) const
        { return detail::kronrod_advance(fn, s); }

        /**
         * @brief Forms the sums afresh from the panels before integrate() decides, so that a decision never rests on
         * running sums that rounding has drifted (iterative_method).
         * @tparam T The real type.
         * @param s The state.
         * @return Nothing, or error::non_finite.
         */
        template<real T>
        constexpr std::expected<void, Failure<T>> settle(detail::KronrodState<T> &s) const
        { return detail::kronrod_settle(s); }
    };

    /**
     * @brief The tanh-sinh rule and its forms for infinite intervals, for an integrand singular at an end (log x,
     * 1 / sqrt x), for many digits in a multiprecision type, and for an infinite interval, which no other method takes.
     *
     * Not for kinks, steps or narrow peaks inside the interval: there it usually fails honestly, with not_converged,
     * where the default, GaussKronrod<15>, subdivides; split the interval at a known feature, or fall back to the
     * default (see integrate()). On smooth fn its cost is comparable to the default's; with an endpoint singularity it
     * is several times cheaper.
     *
     * The trapezoidal rule after the substitution x = tanh(sinh t) (Takahasi and Mori 1974), the step halved level by
     * level; an iterative method, one step per level. Level 0 samples fn at the centre and at t = 1, 2, ..., level k at
     * the odd multiples of 2^-k, so each level reuses every earlier evaluation. The points crowd towards the ends
     * double exponentially but never reach them, so log x or 1 / sqrt x at an end is integrated as it is. Each side
     * stops at its first point that rounds onto the end, at a negligible term, or at a value that is not finite within
     * eps times the half-width of the end and beyond its kept points (a singular end); a value that is not finite
     * elsewhere is error::non_finite. The estimate is the difference of the last two levels plus a bound on what lies
     * beyond each side's last point; it may stop from level 4 and certifies down to twice its rounding floor. The cap
     * is 10 levels, then error::not_converged, without a best where the bound beyond an end is unbounded (the terms
     * grow towards it, as for 1/x).
     *
     * Its steps, observed: report i is level i; level 0 has no current; a later report without one means the estimate
     * is unbounded; `ready` from report 4, `exhausted` only at report 9 (or at report 0, with no point strictly
     * inside). The evaluations roughly double per level. A failure at level 0, such as non_finite for sqrt(x - 0.5),
     * gives no report: record fn's arguments to find it. A run to the cap can be stopped by an observer that throws.
     *
     * The limits on a finite interval, each with its remedy:
     * - A singularity at an end other than 0 is limited by the spacing of T near that end: substitute u = b - x (or
     *   u = x - a) and write fn in u, since fn(b - u) rounds b - u and loses what the substitution gains.
     * - A success's abs_error can exceed the request on an interval a few million ulps wide or less, where the slivers
     *   between the last points and the ends, each up to about an ulp of its end, dominate the floor.
     * - Compile time: its node tables are instantiated wherever integrate<TanhSinh> is for float, double or long
     *   double, both maps' tables whatever the ends, more for long double; instantiate it in one translation unit.
     * - Constant evaluation: fn rational, or built from one's own constexpr functions (std::sqrt and the like are not
     *   constexpr on Clang and MSVC). A class type has no compile-time table, so each level reached is generated in
     *   the same evaluation: only a run that stops within the first few levels fits the default limit.
     * - High precision: the first call in a type without a table builds the levels it reaches, which takes a noticeable
     *   fraction of a second for the whole table, then keeps them; a type that changes its precision at run time sees
     *   its first.
     *
     * Infinite ends (PHASE4_API.md): exp-sinh x = a + E on [a, inf) (b - E on (-inf, b]), sinh-sinh x = +-sinh u on the
     * line, with u = c sinh t, E = e^u and c = 3/2; the side towards a finite end takes x = a + e^-u. Level 0 samples
     * the centre (a + 1, b - 1, or 0) and t = 1, 2, .... It calls fn at very large |x|, up to about 1e305 in double on
     * a side that is not pruned. A value that is not finite towards infinity, beyond the kept points, ends the side as
     * a tail where |x| exceeds 1/eps or the last kept term is already negligible, and is error::non_finite elsewhere; a
     * value below min() there is not taken as negligible, since its weight reaches max / 2. The limits, each with its
     * remedy:
     * - The maps centre at 0 (the line) or at the finite end, at scale 1. A peak far from 0 on the line is missed,
     *   first as an honest not_converged, farther out as a false success (ok, 0): split at the peak, or substitute
     *   x = x0 + u. A half line far from 0 is cut by rounding near a, its abs_error growing with ulp(a), until from
     *   |a| about 1/eps no point is inside: substitute x = a + u with fn written in u. A narrow feature on the line is
     *   not_converged: substitute to widen it. A wide one costs evaluations, not accuracy.
     * - A value not finite nearer than the tail rule allows is error::non_finite (x^k e^-x for large k, where x^k
     *   overflows before e^-x underflows): write fn so that it stays finite where its factors do not, as
     *   x^k e^-x = exp(k log x - x).
     * - Values below min() towards infinity are not used, so an fn that small is not_converged: scale fn.
     * - Slow decay (x^-1.01 on [1, inf)) is not_converged in double: use a type with more digits.
     * - Oscillation and no integral (sin x / x, sin x, 1/x, 1 on [0, inf)): not_converged; their bests, where any,
     *   are not reliable.
     * - A singular end other than 0 is limited by the spacing near it, as on a finite interval: substitute.
     * - A large number written for infinity integrates that finite interval, which may not converge where the infinite
     *   one would: write the infinity.
     *
     * The measured figures behind these limits are in docs/quadrature/design.md, "TanhSinh's measured limits".
     */
    struct TanhSinh
    {
        /** @brief Opts in to infinite ends: exp-sinh on a half line, sinh-sinh on the line (estimating_method). */
        static constexpr bool infinite_ends = true;

        /**
         * @brief Begins the integration, for integrate(): level 0 of the map the ends choose (iterative_method).
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param a The lower end, finite or -inf.
         * @param b The upper end, finite or inf, a < b.
         * @return The state after level 0, or error::non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<detail::TanhSinhState<T>, Failure<T>> start(F &fn, T const &a, T const &b) const
        {
            // Both ends finite: tanh-sinh; otherwise exp-sinh or sinh-sinh, which exp_sinh_start chooses by the ends.
            if (nxx::math::isfinite(a) && nxx::math::isfinite(b)) return detail::tanh_sinh_start(fn, a, b);
            return detail::exp_sinh_start(fn, a, b);
        }

        /**
         * @brief Takes the next level, the step halved, on the map start chose (iterative_method).
         * @tparam T The real type.
         * @tparam F fn's type.
         * @param fn The function.
         * @param s The state, not exhausted.
         * @return Nothing, or error::non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> advance(F &fn, detail::TanhSinhState<T> &s) const
        {
            if (s.map == detail::mapping::interval) return detail::tanh_sinh_advance(fn, s);
            return detail::exp_sinh_advance(fn, s);
        }

        /**
         * @brief Does nothing: TanhSinh's step is final as reported (iterative_method).
         * @tparam T The real type.
         * @param s The state.
         * @return Nothing.
         */
        template<real T>
        constexpr std::expected<void, Failure<T>> settle([[maybe_unused]] detail::TanhSinhState<T> &s) const
        { return {}; }
    };

    /**
     * @brief Integrates @p fn from @p a to @p b with a fixed rule given as a value, for a rule that carries parameters
     * (a user's rule with a number of panels): a plain number, no error estimate and no criteria.
     *
     * The ends are validated (NaN or infinite: error::non_finite, since no fixed rule takes an infinite end) and put
     * in order; a == b is 0 without calling fn; b < a negates. The rule's value is checked: one that is not finite is
     * error::non_finite. fn's exceptions pass through unchanged.
     *
     * @code
     * namespace qd = nxx::quadrature;
     * auto const g = qd::integrate(qd::GaussLegendre<20>{}, [](double x) { return 1.0 / (1.0 + x * x); }, 0.0, 1.0);
     * auto const m = qd::integrate(Midpoint{.panels = 1000}, fn, 0.0, 1.0); // a user's rule with a parameter
     * @endcode
     *
     * @tparam M The rule: a fixed_rule, such as `GaussLegendre<N>` with 1 <= N <= 64, or a user's.
     * @tparam T The real type, deduced from both ends, which must have the same type.
     * @tparam F A callable taking T const & and returning T (for a class type T, T itself or a built-in arithmetic
     *           type, not an expression template; for a built-in T anything convertible to T, so a double result for
     *           a float T narrows silently). The same for every method.
     * @param m The rule.
     * @param fn The function.
     * @param a The lower end of the interval.
     * @param b The upper end.
     * @return The integral; error::non_finite if an end or a value of fn is not finite, or the sum overflowed.
     */
    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T> && fixed_rule<M, T, detail::passed<F>>
    [[nodiscard]] constexpr std::expected<T, error> integrate(M const &m, F &&fn, T const &a, T const &b)
    {
        if (!nxx::math::isfinite(a) || !nxx::math::isfinite(b)) return std::unexpected(error::non_finite);
        if (a == b) return T(0);
        auto &&f = detail::as_object(fn);
        std::expected<T, error> const r =
            b < a ? m.apply(f, b, a).transform([](T const &v) { return T(-v); }) : m.apply(f, a, b);
        // A user's rule may return a value that is not finite; the contract promises the caller one that is.
        if (r && !nxx::math::isfinite(*r)) return std::unexpected(error::non_finite);
        return r;
    }

    /**
     * @brief Integrates @p fn from @p a to @p b with a method given as a value, for a method that carries parameters (a
     * user's method, a combinator around another method): the integral and an estimate of its error, or why there is
     * none.
     *
     * Every other form with an estimate calls this one. The method is an iterative_method, whose loop integrate() runs,
     * or an estimating_method. The criteria are validated (negative or not finite: error::non_finite) and rel_tol
     * clamped to 1; the ends are validated (NaN: error::non_finite; infinite: error::non_finite unless the method takes
     * infinite ends, which TanhSinh alone of the library's methods does) and put in order, (inf, a) as (a, inf) and
     * (inf, -inf) as (-inf, inf); a == b, the same infinity included, is {0, 0} without calling fn; b < a negates the
     * value and a best's value. An iterative method stops when its step is ready and its estimate meets
     * max(criteria.abs_tol, criteria.rel_tol |value|) or the error the method can certify (for the library's methods,
     * the rounding floor; see Criteria). The result is checked: a value or an estimate that is not finite, or a
     * negative estimate, is error::non_finite. The estimate is not a bound (see Integral). fn's exceptions pass through
     * unchanged.
     *
     * @code
     * namespace qd = nxx::quadrature;
     * auto const r = qd::integrate(qd::Romberg{}, [](double x) { return std::exp(x); }, 0.0, 1.0, {.rel_tol = 1e-12});
     * auto const s = qd::integrate(SquareAtLowerEnd<qd::GaussKronrod<15>>{}, fn, 0.0, 1.0); // a user's combinator
     * @endcode
     *
     * @tparam M The method: an iterative_method (`GaussKronrod<15>`, `TanhSinh`, `Romberg`, or a user's) or an
     *           estimating_method.
     * @tparam T The real type, deduced from both ends, which must have the same type.
     * @tparam F A callable taking T const & and returning T, as for the fixed rule.
     * @param m The method.
     * @param fn The function.
     * @param a The lower end of the interval.
     * @param b The upper end.
     * @param criteria When to stop; it converts to T, so `{.rel_tol = 1e-12}` works for any T. Name the member: a
     *        braced `{1e-12}` sets abs_tol.
     * @return The integral and an estimate of its error, returned when the request is met or the estimate is at the
     *         error the method can certify, so a caller who needs the request compares abs_error with it; or a
     *         Failure: error::non_finite if an end, a tolerance or a value of fn is not finite, a tolerance is negative
     *         or a sum overflowed; error::not_converged, with the best estimate where there is one, if the tolerance
     *         was not met (see error).
     */
    // criteria is a non-deduced context (type_identity_t), so T comes from the ends alone and a braced list of doubles
    // converts to a multiprecision T.
    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T>
             && (iterative_method<M, T, detail::passed<F>> || estimating_method<M, T, detail::passed<F>>)
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(M const &m, F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        Failure<T> const invalid{ error::non_finite, std::nullopt };
        // An assert would vanish under NDEBUG, and a NaN tolerance would then never stop: validated here.
        if (!nxx::math::isfinite(criteria.abs_tol) || !nxx::math::isfinite(criteria.rel_tol)
            || criteria.abs_tol < T(0) || criteria.rel_tol < T(0))
            return std::unexpected(invalid);
        // rel_tol |value| stays finite (constant evaluation rejects an overflow); abs_tol enters only a max of finite
        // values, so it needs no clamp, and the floor is the effective minimum of both.
        Criteria<T> const validated{ criteria.abs_tol, criteria.rel_tol > T(1) ? T(1) : criteria.rel_tol };
        if (detail::is_nan(a) || detail::is_nan(b)) return std::unexpected(invalid);
        if constexpr (!detail::takes_infinite_ends<M>) {
            if (!nxx::math::isfinite(a) || !nxx::math::isfinite(b)) return std::unexpected(invalid);
        }
        // Equal infinities too: the same infinity is an empty interval, and fn is never called at an end.
        if (a == b) return Integral<T>{ T(0), T(0) };
        auto &&f = detail::as_object(fn);
        std::expected<Integral<T>, Failure<T>> const r = b < a ? detail::estimate_on(m, f, b, a, validated)
                                                                     .transform(detail::negated<T>)
                                                                     .transform_error([](Failure<T> failure) {
                                                                         failure.best =
                                                                             failure.best.transform(detail::negated<T>);
                                                                         return failure;
                                                                     })
                                                               : detail::estimate_on(m, f, a, b, validated);
        // What every method must return: a sound result, or a Failure whose best, if any, is sound and not_converged's.
        bool const broken = r ? !detail::sound(*r)
                              : r.error().best
                                    && (r.error().reason != error::not_converged || !detail::sound(*r.error().best));
        if (broken) return std::unexpected(invalid);
        return r;
    }

    /**
     * @brief Integrates @p fn from @p a to @p b with the fixed rule M, named, for exactly N evaluations and a plain
     * number: integrate<GaussLegendre<N>>(fn, a, b), the same as integrate(M{}, fn, a, b).
     *
     * fn is evaluated at the N interior nodes, never at an end, except that on an interval a few ulps wide a node can
     * round onto an end. fn's exceptions pass through unchanged.
     *
     * @code
     * namespace qd = nxx::quadrature;
     * auto const g = qd::integrate<qd::GaussLegendre<20>>([](double x) { return 1.0 / (1.0 + x * x); }, 0.0, 1.0);
     * constexpr double nine = qd::integrate<qd::GaussLegendre<7>>([](double x) { return x * x; }, 0.0, 3.0).value();
     * @endcode
     *
     * @tparam M The rule: a fixed_rule that value-initialises, such as `GaussLegendre<N>` with 1 <= N <= 64.
     * @tparam T The real type, deduced from both ends, which must have the same type.
     * @tparam F A callable taking T const & and returning T, as for the rule given as a value.
     * @param fn The function.
     * @param a The lower end of the interval.
     * @param b The upper end.
     * @return The integral; error::non_finite if an end or a value of fn is not finite, or the sum overflowed.
     */
    template<typename M, real T, typename F>
    requires std::default_initializable<M> && nxx::detail::callable<F, T>
             && fixed_rule<M, T, detail::passed<F>>
    [[nodiscard]] constexpr std::expected<T, error> integrate(F &&fn, T const &a, T const &b)
    { return nxx::quadrature::integrate(M{}, fn, a, b); }

    /**
     * @brief Integrates @p fn from @p a to @p b with the method M, named, or the default GaussKronrod<15>: the form for
     * most calls, returning the integral and an estimate of its error, or why there is none.
     *
     * The same as integrate(M{}, fn, a, b, criteria).
     *
     * @code
     * namespace qd = nxx::quadrature;
     * auto const r = qd::integrate([](double x) { return std::exp(x); }, 0.0, 1.0);
     * auto const s = qd::integrate<qd::Romberg>([](double x) { return std::exp(x); }, 0.0, 1.0, {.rel_tol = 1e-12});
     * auto const t = qd::integrate<qd::TanhSinh>([](double x) { return std::log(x); }, 0.0, 1.0); // -1
     * auto const e = qd::integrate<qd::TanhSinh>([](double x) { return x * x * std::exp(-x); }, 0.0, inf); // 2
     * if (r) use(r->value, r->abs_error);
     * else if (r.error().best) report(r.error().best->value, r.error().best->abs_error);
     * @endcode
     *
     * Choosing a method (each method's comment says more):
     * - An endpoint singularity (log x, x^-0.9), or high precision in a multiprecision type: TanhSinh, several times
     *   cheaper than the default with an endpoint singularity, comparable on smooth fn.
     * - Kinks, steps, narrow peaks inside the interval, or an fn whose shape is unknown: the default, which
     *   subdivides. TanhSinh usually fails honestly there, with not_converged, which is why it is not the default.
     * - An fn too expensive for about 100 evaluations, or an integral inside another algorithm: GaussLegendre<N>,
     *   which has no estimate.
     * - An infinite interval: TanhSinh, the only method that takes an infinite end (see TanhSinh for its limits).
     *
     * A fallback from TanhSinh to the default for interior features. It accepts the default's success only where it
     * agrees with TanhSinh's best within the two estimates, since at an endpoint singularity the default can succeed
     * falsely where TanhSinh failed honestly. It cannot catch an error both methods share, as at a singular end other
     * than 0, where both lose what lies within an ulp of the end: substitute (see TanhSinh) instead. It is for finite
     * intervals: the default rejects an infinite end, so there it returns TanhSinh's failure unchanged. On a half line
     * with interior features, split at a finite c: the default on [a, c], TanhSinh on [c, inf).
     *
     * @code
     * auto const t = qd::integrate<qd::TanhSinh>(g, a, b);
     * auto const r = t.or_else([&](qd::Failure<double> const &f) -> decltype(t) {
     *     if (f.reason != qd::error::not_converged || !f.best) return t; // non_finite, or no estimate at all
     *     auto const k = qd::integrate(g, a, b);                          // interior features: the default
     *     if (k && std::abs(k->value - f.best->value) <= f.best->abs_error + k->abs_error) return k;
     *     return t;
     * });
     * @endcode
     *
     * @tparam M The method: an iterative_method or an estimating_method that value-initialises; `GaussKronrod<15>` by
     *           default.
     * @tparam T The real type, deduced from both ends, which must have the same type.
     * @tparam F A callable taking T const & and returning T, as for the fixed rule.
     * @param fn The function.
     * @param a The lower end of the interval.
     * @param b The upper end.
     * @param criteria When to stop; it converts to T, so `{.rel_tol = 1e-12}` works for any T. Name the member: a
     *        braced `{1e-12}` sets abs_tol.
     * @return As for the method given as a value.
     */
    template<typename M = GaussKronrod<15>, real T, typename F>
    requires std::default_initializable<M> && nxx::detail::callable<F, T>
             && (iterative_method<M, T, detail::passed<F>> || estimating_method<M, T, detail::passed<F>>)
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return nxx::quadrature::integrate(M{}, fn, a, b, criteria); }

    /**
     * @brief Integrates @p fn from @p a to @p b with the iterative method M, named, showing each step to the observer
     * @p o, to watch the method converge or find why it failed: the same result as integrate<M>(fn, a, b, criteria).
     *
     * @code
     * namespace qd = nxx::quadrature;
     * auto const r = qd::integrate<qd::Romberg>(qd::observer(print), fn, 0.0, 1.0, {.rel_tol = 1e-12});
     * auto const s = qd::integrate(qd::observer(print), fn, 0.0, 1.0); // GaussKronrod<15>
     * @endcode
     *
     * @tparam M The method: an iterative_method that value-initialises; `GaussKronrod<15>` by default.
     * @tparam O The observer's callable (see observer()).
     * @tparam T The real type, deduced from both ends, which must have the same type.
     * @tparam F A callable taking T const & and returning T, as for the fixed rule.
     * @param o The observer, made by observer().
     * @param fn The function.
     * @param a The lower end of the interval.
     * @param b The upper end.
     * @param criteria When to stop.
     * @return As for the method given as a value.
     */
    template<typename M = GaussKronrod<15>, typename O, real T, typename F>
    requires std::default_initializable<M> && nxx::detail::callable<F, T>
             && iterative_method<M, T, detail::passed<F>>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>> integrate(Observer<O> o, F &&fn, T const &a,
        T const &b, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return nxx::quadrature::integrate(detail::observed(M{}, std::move(o)), fn, a, b, criteria); }

    /**
     * @brief Integrates @p fn from @p a to @p b with the iterative method @p m, given as a value, showing each step to
     * the observer @p o: the same result as integrate(m, fn, a, b, criteria).
     *
     * @code
     * namespace qd = nxx::quadrature;
     * auto const r = qd::integrate(Trapezoid{}, qd::observer(print), fn, 0.0, 1.0); // a user's iterative method
     * @endcode
     *
     * @tparam M The method: an iterative_method.
     * @tparam O The observer's callable (see observer()).
     * @tparam T The real type, deduced from both ends, which must have the same type.
     * @tparam F A callable taking T const & and returning T, as for the fixed rule.
     * @param m The method.
     * @param o The observer, made by observer().
     * @param fn The function.
     * @param a The lower end of the interval.
     * @param b The upper end.
     * @param criteria When to stop.
     * @return As for the method given as a value.
     */
    template<typename M, typename O, real T, typename F>
    requires nxx::detail::callable<F, T> && iterative_method<M, T, detail::passed<F>>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>> integrate(M const &m, Observer<O> o, F &&fn,
        T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return nxx::quadrature::integrate(detail::observed(m, std::move(o)), fn, a, b, criteria); }

} // namespace nxx::quadrature

/**
 * @brief The standard library, for the std::formatter specialisations that let std::format and std::print write
 * nxx::quadrature::error and nxx::quadrature::Failure.
 */
namespace std {

    /**
     * @brief Lets std::format and std::print write an nxx::quadrature::error, as the name of its enumerator, such as
     * `not_converged`, with the format specification of a string.
     *
     * @code
     * if (!g) std::println("no integral: {}", g.error());   // no integral: non_finite
     * @endcode
     */
    template<>
    struct formatter<nxx::quadrature::error> : formatter<string_view>
    {
        /**
         * @brief Writes the name of @p e, for std::format.
         *
         * @tparam FormatContext The format context.
         * @param e The error.
         * @param ctx The format context.
         * @return The end of the output.
         */
        template<typename FormatContext>
        auto format(nxx::quadrature::error e, FormatContext &ctx) const
        { return formatter<string_view>::format(nxx::quadrature::detail::name(e), ctx); }
    };

    /**
     * @brief Lets std::format and std::print write an nxx::quadrature::Failure, as the name of its reason, such as
     * `not_converged`, with the format specification of a string.
     *
     * The best estimate is not written: that would need a formatter for T, which Boost.Multiprecision's types do not
     * have; format `f.best->value` where T has one.
     *
     * @code
     * if (!r) std::println("no integral: {}", r.error());   // no integral: not_converged
     * @endcode
     *
     * @tparam T The real type.
     */
    template<nxx::real T>
    struct formatter<nxx::quadrature::Failure<T>> : formatter<string_view>
    {
        /**
         * @brief Writes the name of @p f's reason, for std::format.
         *
         * @tparam FormatContext The format context.
         * @param f The failure.
         * @param ctx The format context.
         * @return The end of the output.
         */
        template<typename FormatContext>
        auto format(nxx::quadrature::Failure<T> const &f, FormatContext &ctx) const
        { return formatter<string_view>::format(nxx::quadrature::detail::name(f.reason), ctx); }
    };

} // namespace std
