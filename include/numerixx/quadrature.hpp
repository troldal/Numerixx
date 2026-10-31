/**
 * @file
 * @brief Integrals of a real function of one variable over a finite interval: adaptive Gauss-Kronrod, the default;
 * Romberg's method; and the fixed Gauss-Legendre rule.
 *
 * integrate(fn, a, b) integrates fn from a to b by adaptive Gauss-Kronrod, the default method, and integrate<M>(fn, a,
 * b) by a method M: `GaussKronrod<15>`, the default named; `Romberg`; or `GaussLegendre<N>`, a fixed rule that returns
 * a plain number. The methods with an estimate return an Integral, the value and an estimate of its error, or a
 * Failure. Every method returns a std::expected; Numerixx itself never throws, but the result's own value() throws
 * std::bad_expected_access on a failure, so test the result and read it with `->` or `*`.
 *
 * @code
 * namespace qd = nxx::quadrature;
 * auto const r = qd::integrate([](double x) { return std::exp(x); }, 0.0, 1.0); // e - 1, 150 evaluations
 * if (r) use(r->value, r->abs_error);
 * @endcode
 *
 * The ends are oriented: b < a gives the negated integral, and a == b gives 0 without calling fn. An interval narrow
 * against its position (|a + b| far above |b - a|) rounds each point at which fn is evaluated by about eps |a + b| / 2,
 * so the result errs by about that times the integral of |f'|: fn's conditioning, which no method sees.
 */

#pragma once

#include <numerixx/math.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <expected>
#include <format>
#include <functional>
#include <limits>
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

/** @brief Integrals of real functions of one variable. */
namespace nxx::quadrature {

    /**
     * @brief Why integrate() returned no result.
     *
     * More values may be added, and a listed cause may move to a new value: a switch on it should have a default. It
     * is formattable: std::format("{}", e) writes the enumerator's name, such as `not_converged`.
     */
    enum class error {
        /**
         * An end is NaN or infinite (an infinite interval is for a later method, TanhSinh), a tolerance is negative or
         * not finite, fn returned a value that is not finite (for GaussKronrod: on a panel of its initial mesh that no
         * short bisection avoids, such as a region where fn is NaN), or a sum overflowed.
         */
        non_finite,

        /**
         * The tolerance was not met: within the method's cap; or, for GaussKronrod, by stagnation, when its
         * bisections stop lowering the estimate's excess over the floor (fn noisier than the request, or a divergent
         * integral); or with every panel final and their errors above the tolerance (a strong endpoint singularity,
         * or a piece where fn is not finite that a short bisection did not avoid). Failure::best then holds the
         * estimate reached and its error estimate.
         */
        not_converged
    };

    /**
     * @brief An integral, and an estimate of its error.
     *
     * abs_error estimates |value - exact|; it is not a bound. No method that samples fn at finitely many points can
     * bound the error: a narrow peak between the points is missed, and an fn that vanishes at every point sampled
     * gives 0 with the floor as its estimate, whatever its integral. abs_error is never below the rounding floor (see
     * Criteria), so an integral that is exactly 0 reports the floor, not 0. Read the members by name: members may be
     * added in a later version, which would break a structured binding.
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Integral
    {
        /** @brief The integral. */
        T value;

        /** @brief An estimate of |value - exact|, not a bound; never below the rounding floor. */
        T abs_error;

        /** @brief Compares both members. */
        friend constexpr bool operator==(Integral const &, Integral const &) = default;
    };

    /**
     * @brief Why a method with an error estimate returned no Integral, and, for error::not_converged, the best
     * estimate it reached.
     *
     * Formattable: std::format("{}", f) writes the reason's name, such as `not_converged`.
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Failure
    {
        /** @brief Why there is no result. */
        error reason;

        /** @brief For error::not_converged the last estimate and its error estimate; empty for error::non_finite. */
        std::optional<Integral<T>> best;

        /** @brief Compares both members. */
        friend constexpr bool operator==(Failure const &, Failure const &) = default;
    };

    /**
     * @brief When to stop: when the error estimate is at most max(abs_tol, rel_tol |value|), or at the rounding floor.
     *
     * The rounding floor is 50 max(eps L1, denorm_min), L1 the integral of |fn| by the same method on the same
     * evaluations (for GaussKronrod the sum of every panel's): the scale of the sum's own rounding. It is the effective
     * minimum of every request (about 1e-14 relative in double): GaussKronrod accepts at twice the floor, Romberg at
     * the floor, and either ends an integral that cancels to 0 there, so abs_tol = 0 works. The defaults ask for
     * sqrt(eps) relative, not full precision: no estimate from samples certifies the last digits. A caller who chose a
     * 50-digit type for its digits sets rel_tol. nxx::roots::Criteria asks for full precision by default; under
     * using-directives for both namespaces this name and `error` collide with nxx::roots', so qualify them.
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
     * @brief The N-point Gauss-Legendre rule, 1 <= N <= 64: exact for polynomials of degree 2N - 1, up to the rounding
     * of its nodes and weights and of the sum.
     *
     * A fixed rule: N evaluations of fn at the zeros of the Legendre polynomial P_N mapped to (a, b), no error estimate
     * and no criteria, so that the result is a smooth function of parameters inside fn (Lyness 1983): the rule for a
     * root finder, an optimiser or an inner integral. Its error on [a, b] is
     * (b - a)^(2N+1) (N!)^4 / ((2N + 1) ((2N)!)^3) f^(2N)(xi) for some xi in (a, b).
     *
     * The nodes and weights come from the three-term recurrence, + - * / in T: for float, double and long double at
     * compile time, once per program; for every other type, such as a Boost.Multiprecision number, at its first use,
     * and then kept (a function-local static, initialised once by the language; a type that changed its precision at
     * run time would see the nodes of its first).
     *
     * @tparam N The number of points, 1 to 64.
     */
    template<int N>
    struct GaussLegendre
    {
    };

    /**
     * @brief Romberg's method: the trapezoidal rule on halved steps, extrapolated by Richardson (Romberg 1955).
     *
     * Level k samples fn at 2^k + 1 equally spaced points, ends included, so a singularity at an end is
     * error::non_finite. It stops at a level k >= 4 (17 evaluations) where the difference of the last two extrapolated
     * values, |R(k, k) - R(k-1, k-1)|, meets the criteria, and returns R(k, k); that difference estimates the error of
     * the previous value, so it overstates the error on smooth fn. The cap is 16 levels, 65,537 evaluations, then
     * error::not_converged with the best estimate; fn with a kink or a singular derivative reach it. It needs more
     * evaluations than a Gauss rule for the same accuracy; its uses are teaching and equally spaced data.
     */
    struct Romberg
    {
    };

    /**
     * @brief The G7-K15 Gauss-Kronrod pair, globally adaptive: the default method. Only N = 15 exists.
     *
     * The interval starts as 10 equal panels, 150 evaluations of fn (for an expensive fn that must cost less, use
     * GaussLegendre<N>, which has no estimate, or Romberg); then the panel with the largest estimate is bisected
     * until the sum of the estimates meets the criteria. A panel's estimate is |K15 - G7| on its 15 points, or its
     * floor if larger: it estimates the 7-point rule's error, not the returned 15-point result's, so it usually
     * overstates the error by orders of magnitude, but it is not a bound. What it misses: a feature narrower than the
     * spacing of the points (a narrow peak between them), and a kink or a step near a panel's ends, outside the
     * outermost points, which neither rule samples. An endpoint singularity stronger than about (x - a)^-0.65 is
     * underestimated by a ratio that bisection does not change; substitute x = a + u^k with k (1 + alpha) >= 1 for
     * (x - a)^alpha (and u = b - x for the other end) until TanhSinh. Full precision in a 50-digit type is costly:
     * the estimate converges at the 7-point rule's order.
     *
     * @tparam N The number of Kronrod points: 15.
     */
    template<int N>
    requires(N == 15)
    struct GaussKronrod
    {
    };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief The number of points of a GaussLegendre rule, or 0 for any other type.
         * @tparam M The method.
         */
        template<typename M>
        inline constexpr int gauss_points = 0;

        /**
         * @brief The number of points of GaussLegendre<N>: N.
         * @tparam N The number of points.
         */
        template<int N>
        inline constexpr int gauss_points<GaussLegendre<N>> = N;

        /**
         * @brief The fixed rules integrate() accepts: GaussLegendre<N> with 1 <= N <= 64.
         * @tparam M The method.
         */
        template<typename M>
        concept fixed_rule = gauss_points<M> >= 1 && gauss_points<M> <= 64;

        /**
         * @brief The methods with an estimate that integrate() accepts: GaussKronrod<15> and Romberg.
         * @tparam M The method.
         */
        template<typename M>
        concept estimating = std::same_as<M, GaussKronrod<15>> || std::same_as<M, Romberg>;

        /**
         * @brief The types whose Gauss-Legendre rules are computed at compile time: float, double and long double,
         * named, since std::is_floating_point_v also admits the unmeasured std::float16_t and std::bfloat16_t.
         * @tparam T The real type.
         */
        template<typename T>
        concept tabled = std::same_as<T, float> || std::same_as<T, double> || std::same_as<T, long double>;

        /**
         * @brief The name of an error's enumerator.
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
         * @brief The larger of @p a and @p b.
         * @param a A value.
         * @param b A value.
         * @return b if a < b, else a.
         */
        template<real T>
        constexpr T larger(T const &a, T const &b)
        { return a < b ? b : a; }

        /**
         * @brief Whether |@p x| > max / 2: the constant-evaluation guard, below which a sum of two such values, or a
         * difference, cannot overflow.
         * @param x A finite value.
         * @return |x| > max / 2.
         */
        template<real T>
        constexpr bool large(T const &x)
        { return nxx::math::abs(x) > T(T((std::numeric_limits<T>::max)()) / T(2)); }

        /**
         * @brief Whether |@p a b| <= max / 2, decided without an overflow: for |a| <= 1 the product, at most |b|, is
         * formed and compared; else |b| <= (max / 2) / |a|, which bounds the product by (max / 2)(1 + eps) < max.
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
         * @brief Whether @p a + @p b may be formed: at run time always (an overflow is found by a later check of the
         * sum); in constant evaluation, where an overflow is not a constant expression, only if neither exceeds
         * max / 2.
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
         * @brief An array of M elements, each T(0): nxx::real does not require T to be default constructible.
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
         * @brief fn at @p x, converted to T.
         * @param fn The function.
         * @param x The point, a T const lvalue as the constraint promises.
         * @return fn(x) as a T.
         */
        template<real T, typename F>
        constexpr T evaluate(F &fn, T const &x)
        { return static_cast<T>(std::invoke(fn, x)); }

        /**
         * @brief A sum compensated by Neumaier's variant of Kahan's method: the rounding error of each addition is
         * accumulated separately, so that the result errs by about eps |sum| rather than n eps sum |terms|.
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
         * @brief The interval [a, b], a != b, both finite, as the methods use it: its centre and half-width, which
         * place the points, and the scale of each term.
         * @tparam T The real type.
         */
        template<real T>
        struct Interval
        {
            /** @brief a / 2 + b / 2. */
            T c;

            /** @brief b / 2 - a / 2, signed; neither it nor c can overflow, on [-max, max] included. */
            T hw;

            /** @brief b - a, where finite. */
            T width;

            /** @brief Whether b - a is finite. */
            bool finite;

            /** @brief Whether |b - a| >= 1. */
            bool wide;

            /** @brief min(a, b). */
            T lo;

            /** @brief max(a, b). */
            T hi;

            /**
             * @brief The point @p x, clamped into [lo, hi].
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
             * @brief The term hw w y, formed before it is summed (a sum of values formed first could overflow).
             *
             * The order keeps every intermediate normal unless the term itself is subnormal or overflows. Where b - a
             * is finite, the term is (b - a) (w / 2) y: with |b - a| >= 1 the scale (b - a) (w / 2) is formed first,
             * which lies between w / 2 and max (w / 2 <= 1), so that a tiny y on a wide interval is not made subnormal
             * first; with |b - a| < 1 the value is multiplied first, so that a subnormal scale (b - a) w / 2 is never
             * formed. Where b - a overflows, |hw| > max / 2: (hw y) w for |y| <= 1, where hw y cannot overflow and is
             * normal unless y is 0, else hw (w y), where |w y| >= w is normal (hw w would overflow for w > 1). w / 2
             * is exact in radix 2.
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
             * @brief The term, or nothing where constant evaluation could not form it within max / 2.
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
         * @brief a / 2 + b / 2, which cannot overflow: the centre of an interval, and where a panel is bisected.
         * @param a An end.
         * @param b The other end.
         * @return The midpoint.
         */
        template<real T>
        constexpr T midpoint(T const &a, T const &b)
        { return T(T(a / T(2)) + T(b / T(2))); }

        /**
         * @brief The rounding floor of an estimate: 50 max(eps L1, denorm_min), QUADPACK's constant.
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
         * @brief The estimate a method accepts: max(abs_tol, rel_tol |value|, floor), with the floor its own
         * (Romberg's, or twice Gauss-Kronrod's).
         * @param abs_tol The absolute tolerance.
         * @param rel_tol The relative tolerance, at most 1.
         * @param value The integral.
         * @param floor The floor the method accepts at.
         * @return The largest estimate accepted.
         */
        template<real T>
        constexpr T wanted(T const &abs_tol, T const &rel_tol, T const &value, T const &floor)
        {
            T const relative = T(rel_tol * nxx::math::abs(value));
            return detail::larger(detail::larger(abs_tol, relative), floor);
        }

        /**
         * @brief The interval [a, b] for the methods.
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
         * @brief The number of bits of T's significand: digits times the bits of a digit of its radix, rounded up.
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
         * @brief A value of P_N and of its derivative.
         * @tparam T The real type.
         */
        template<real T>
        struct Legendre
        {
            /** @brief P_N(x). */
            T p;

            /** @brief P_N'(x). */
            T dp;
        };

        /**
         * @brief P_(k+1)(x) from P_k(x) and P_(k-1)(x): (k + 1) P_(k+1) = (2k + 1) x P_k - k P_(k-1).
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
         * @brief P'_(k+1)(x), the recurrence differentiated: (k + 1) P'_(k+1) = (2k + 1) (P_k + x P'_k) - k P'_(k-1).
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
         * @brief P_n(x) and P_n'(x) by the three-term recurrence and its derivative.
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
         * @brief The number of zeros of P_n strictly above @p x: the sign changes of P_0(x), ..., P_n(x).
         *
         * The sequence is a Sturm sequence (the recurrence of orthogonal polynomials): its sign changes fall by one
         * as x passes a zero of P_n upward, from n below -1 to 0 above 1. A zero term is skipped and the last nonzero
         * sign kept: where an intermediate P_k vanishes its neighbours have opposite signs, and where P_n vanishes at
         * x, x itself is not counted.
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
         * @brief The zero of a polynomial in (@p lower, @p upper) by Newton's method, safeguarded by bisection: P_n for
         * the Gauss nodes, the Stieltjes polynomial for the Kronrod nodes.
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
         * @brief The non-negative half of the N-point Gauss-Legendre rule: nodes ascending, and their weights.
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
         * @brief The N-point Gauss-Legendre rule, from + - * / in T alone (no pi, no cos).
         *
         * The positive zeros of P_N are isolated one by one, from the largest, by bisection on zeros_above(), then
         * found by bracketed_zero(); the negative nodes are their mirror images, and 0 is a node for odd N. The weights
         * are 2 / ((1 - x)(1 + x) P_N'(x)^2).
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
         * once per program.
         * @tparam N The number of points.
         * @tparam T float, double or long double.
         */
        template<int N, typename T>
        inline constexpr Rule<T, N> gauss_legendre_table = detail::gauss_legendre_rule<N, T>();

        /**
         * @brief The rule @p generate makes, built at its first use and kept: one function-local static per generator,
         * so per rule and T, whatever fn is (initialised once, thread-safely, by the language; never changed). Not
         * for constant evaluation.
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
         * @brief @p use applied to the rule @p generate makes for a type without a compile-time table: the generator
         * itself in constant evaluation, the cached rule at run time.
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
         * @brief @p use applied to the N-point rule in T: the compile-time table for float, double and long double,
         * with_generated() for every other T.
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
         * @brief The coefficients of the Stieltjes polynomial E_8 of the G7-K15 pair in the Legendre basis:
         * E_8 = sum_k a_k P_k, a_8 = 1, orthogonal to P_7 P_j for every j < 8.
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
         * @brief sum_k a_k P_k(x) and its derivative, by the three-term recurrence.
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
         * @brief The non-negative half of the G7-K15 pair: the 8 nodes ascending, the 7-point Gauss nodes at the even
         * indices (x[0] = 0), the Kronrod nodes at the odd ones, with both rules' weights (the Gauss weight 0 where a
         * node is Kronrod's only).
         * @tparam T The real type.
         */
        template<real T>
        struct Kronrod
        {
            /** @brief The number of non-negative nodes. */
            static constexpr std::size_t size = 8;

            /** @brief The non-negative nodes, ascending. */
            std::array<T, size> x;

            /** @brief The 15-point Kronrod weights. */
            std::array<T, size> wk;

            /** @brief The 7-point Gauss weights at the even indices, 0 at the odd ones. */
            std::array<T, size> wg;
        };

        /**
         * @brief The G7-K15 pair, from + - * / in T alone.
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
         * @brief The G7-K15 pair for float, double and long double, computed at compile time once per program.
         * @tparam T float, double or long double.
         */
        template<typename T>
        inline constexpr Kronrod<T> gauss_kronrod_table = detail::gauss_kronrod_rule<T>();

        /**
         * @brief @p use applied to the G7-K15 pair in T, as with_gauss_legendre(): the table, or with_generated().
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
         * @brief A rule applied to fn on an interval: sum over the nodes of hw w_i fn(c + hw x_i), compensated.
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
            T const result = sum.value();
            if (!nxx::math::isfinite(result)) return std::unexpected(error::non_finite);
            return result;
        }

        /**
         * @brief Romberg's method on an interval, with validated tolerances.
         * @param fn The function.
         * @param a The lower end.
         * @param b The upper end.
         * @param abs_tol The absolute tolerance, finite and non-negative.
         * @param rel_tol The relative tolerance, in [0, 1].
         * @return The integral, or the Failure.
         */
        template<real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            romberg(F &fn, T const &a, T const &b, T const &abs_tol, T const &rel_tol)
        {
            constexpr int cap = 16;      // 2^16 + 1 evaluations; fixed, not from T (docs/quadrature/design.md)
            constexpr int min_level = 4; // 17 evaluations: fewer samples of sin^2 on [0, 2 pi] all vanish
            constexpr std::size_t row_size = cap + 1;
            Failure<T> const overflow{ error::non_finite, std::nullopt };
            Interval<T> const interval = detail::interval_of(a, b);

            // h_k = hw w, w relative to the half-width: T_0's h_0 / 2 and level 1's h_1 are hw itself (w = 1), level
            // k's h_k is hw 2^(1-k). A non-finite value of fn makes L1 non-finite, so L1 is the check per level.
            // In constant evaluation every term, sum, L1 and extrapolated value is kept within max / 2 before the
            // operation that could overflow (Interval::checked_term); at run time L1 and the estimate are checked.
            T weight = T(1);
            std::optional<T> const ta = interval.checked_term(weight, detail::evaluate(fn, a));
            std::optional<T> const tb = interval.checked_term(weight, detail::evaluate(fn, b));
            if (!ta || !tb) return std::unexpected(overflow);
            T l1 = T(nxx::math::abs(*ta) + nxx::math::abs(*tb));
            if (!nxx::math::isfinite(l1)) return std::unexpected(overflow);
            std::array<T, row_size> previous = detail::zeros<T, row_size>();
            std::array<T, row_size> row = detail::zeros<T, row_size>();
            previous[0] = T(*ta + *tb);

            Integral<T> best{ T(0), T(0) };
            for (int k = 1; k <= cap; ++k) {
                if (k > 1) weight = T(weight / T(2));
                int const m = 1 << (k - 1);
                Sum<T> sum;
                T abs_sum = T(0); // a plain sum: its terms are positive, so it errs by about m eps relative
                for (int j = 1; j <= m; ++j) {
                    // The new midpoints u_j = (2j - 1 - m) / m, exact (at most 16 bits).
                    T const u = T(T(2 * j - 1 - m) / T(m));
                    T const x = interval.point(T(interval.c + T(interval.hw * u)));
                    std::optional<T> const t = interval.checked_term(weight, detail::evaluate(fn, x));
                    if (!t) return std::unexpected(overflow);
                    if consteval {
                        if (detail::large(sum.sum) || detail::large(abs_sum)) return std::unexpected(overflow);
                    }
                    sum.add(*t);
                    abs_sum = T(abs_sum + nxx::math::abs(*t));
                }
                if consteval {
                    if (detail::large(sum.value()) || detail::large(abs_sum)) return std::unexpected(overflow);
                }
                row[0] = T(T(previous[0] / T(2)) + sum.value());
                l1 = T(T(l1 / T(2)) + abs_sum);
                if (!nxx::math::isfinite(l1)) return std::unexpected(overflow);

                // 4^m - 1 is exact in radix 2 up to the cap for double, and for float up to m = 12, beyond which the
                // correction is at rounding level anyway.
                T power = T(1);
                for (int i = 1; i <= k; ++i) {
                    std::size_t const s = static_cast<std::size_t>(i);
                    if consteval {
                        if (detail::large(row[s - 1]) || detail::large(previous[s - 1]))
                            return std::unexpected(overflow);
                    }
                    power = T(power * T(4));
                    row[s] = T(row[s - 1] + T(T(row[s - 1] - previous[s - 1]) / T(power - T(1))));
                }
                std::size_t const diagonal = static_cast<std::size_t>(k);
                // The difference of the last two diagonal entries: it estimates the error of R(k-1, k-1), so it is
                // pessimistic for R(k, k) on smooth fn, and stays honest where Richardson's assumptions fail (sqrt x),
                // where the subdiagonal |R(k, k) - R(k, k-1)| would claim 4^k - 1 times too little. An overflow of
                // T_k or of the extrapolation reaches it (R(k-1, k-1) is finite), so it is checked once.
                if consteval {
                    if (detail::large(row[diagonal])) return std::unexpected(overflow);
                }
                T const estimate = nxx::math::abs(T(row[diagonal] - previous[diagonal - 1]));
                if (!nxx::math::isfinite(estimate)) return std::unexpected(overflow);
                T const floor = detail::rounding_floor(l1);
                best = Integral<T>{ row[diagonal], detail::larger(estimate, floor) };
                if (k >= min_level && estimate <= detail::wanted(abs_tol, rel_tol, row[diagonal], floor)) return best;
                previous = row;
            }
            return std::unexpected(Failure<T>{ error::not_converged, best });
        }

        /**
         * @brief A panel of the adaptive Gauss-Kronrod loop, with its 15-point value and its estimate.
         * @tparam T The real type.
         */
        template<real T>
        struct Panel
        {
            /** @brief The lower end. */
            T lo;

            /** @brief The upper end. */
            T hi;

            /** @brief The 15-point Kronrod value. */
            T value;

            /** @brief |K15 - G7|, which orders the heap. */
            T difference;

            /** @brief The rounding floor, 50 max(eps L1, denorm_min), L1 the panel's integral of |fn|. */
            T floor;

            /** @brief The estimate's excess over the floor, max(difference, floor) - floor. */
            T excess;

            /** @brief The order of creation: panels of equal difference leave the heap in it. */
            std::size_t order;
        };

        /** @brief Why a panel has no estimate. */
        enum class fault {
            not_finite, ///< L1 not finite: fn not finite at a point of the panel (or an overflow of its terms)
            overflow    ///< a finite L1 but a value or difference that overflowed, or (in constant evaluation) a value
                        ///< or a sum the checks refuse: the run is non_finite
        };

        /**
         * @brief The G7-K15 pair on [lo, hi], lo < hi: the two rules on the same 15 evaluations of fn, each a
         * compensated sum, and L1, the plain sum of |term| (positive terms).
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
            if (!nxx::math::isfinite(l1)) return std::unexpected(fault::not_finite);
            T const k = kronrod.value();
            T const g = gauss.value();
            if (!detail::fits(k, g)) return std::unexpected(fault::overflow);
            T const difference = nxx::math::abs(T(k - g));
            // Beyond the spec's check of L1: a finite L1 whose value or difference overflows (the 7-point weights
            // are larger than the 15-point ones at the Gauss nodes) ends the run non_finite.
            if (!nxx::math::isfinite(k) || !nxx::math::isfinite(difference)) return std::unexpected(fault::overflow);
            T const floor = detail::rounding_floor(l1);
            return Panel<T>{ lo, hi, k, difference, floor, T(detail::larger(difference, floor) - floor), order };
        }

        /**
         * @brief Whether [lo, hi] can be bisected: its midpoint lies strictly inside, and the outer nodes of both
         * children lie strictly inside them (else the nodes round onto the ends).
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
         * @brief A panel waiting in the descent, with its depth below the panel that started it.
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
         * @brief Adaptive Gauss-Kronrod on [a, b], a != b, both finite, with validated tolerances
         * (docs/quadrature/PHASE2_API.md C to E).
         * @param rule The G7-K15 pair.
         * @param fn The function.
         * @param a The lower end.
         * @param b The upper end.
         * @param abs_tol The absolute tolerance, finite and non-negative.
         * @param rel_tol The relative tolerance, in [0, 1].
         * @return The integral, or the Failure.
         */
        template<real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            gauss_kronrod(Kronrod<T> const &rule, F &fn, T const &a, T const &b, T const &abs_tol, T const &rel_tol)
        {
            constexpr int mesh = 10;       // Shampine's quadgk: one panel cannot accept K = G on a peak it missed
            constexpr int depth = 2;       // bisections below a panel whose L1 is not finite
            constexpr int stagnation = 32; // bisections without a new minimum of the excess X
            // Bisections and descent splits together, per bit of T: x^-0.9 on [0, 1] at cpp_bin_float_50's default
            // needs 3,550 of its 10,752.
            int const cap = 64 * detail::significand_bits<T>();
            Failure<T> const invalid{ error::non_finite, std::nullopt };
            // Integrated over [min, max] and negated, so that the panels' ends are in order.
            bool const reversed = b < a;
            T const lower = reversed ? b : a;
            T const upper = reversed ? a : b;

            // The running sums, compensated: the value K, the floor F and the excess X = sum (e_p - floor_p); the
            // estimate is E = F + X. A panel's contribution is subtracted with the same bits that were added, so X,
            // which panels at their floor leave at exactly 0, does not drift.
            Sum<T> total;
            Sum<T> floor;
            Sum<T> excess;
            bool refused = false; // a sum constant evaluation could not form within max / 2
            auto const contribute = [&](Panel<T> const &p, bool add) {
                T const v = add ? p.value : T(-p.value);
                T const f = add ? p.floor : T(-p.floor);
                T const x = add ? p.excess : T(-p.excess);
                if (!detail::fits(total.sum, v) || !detail::fits(floor.sum, f) || !detail::fits(excess.sum, x)) {
                    refused = true;
                    return;
                }
                total.add(v);
                floor.add(f);
                excess.add(x);
            };

            // The active panels, a heap by difference (the same order as the estimate, since only panels above their
            // floor enter it), growing on demand; final panels stay only in the sums.
            std::vector<Panel<T>> heap;
            auto const before = [](Panel<T> const &p, Panel<T> const &q) {
                return p.difference < q.difference || (p.difference == q.difference && q.order < p.order);
            };
            std::size_t order = 0;
            int work = 0; // bisections and descent splits
            std::vector<Panel<T>> pieces;
            std::vector<Pending<T>> stack;

            // [p, q] sampled into pieces: the panel itself, or, where its L1 is not finite, the pieces of a descent of
            // at most `depth` bisections below it, each counted against the cap, each child first its lower half.
            // Returns the fault that remains, if any: overflow ends the run; not_finite, a piece the descent could not
            // avoid.
            auto const resolve = [&](T const &p, T const &q) -> std::optional<fault> {
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
                        stack.push_back(Pending<T>{ mid, top.hi, top.level + 1 });
                        stack.push_back(Pending<T>{ top.lo, mid, top.level + 1 });
                    }
                    else missing = true;
                }
                if (missing) return fault::not_finite;
                return std::nullopt;
            };
            // A piece enters the sums, and the heap unless it is final: at its floor (splitting cannot lower it below
            // its floor, since L1 adds up), or not splittable.
            auto const commit = [&]() {
                for (Panel<T> const &piece : pieces) {
                    contribute(piece, true);
                    if (!(piece.difference <= piece.floor) && detail::splittable(rule, piece.lo, piece.hi)) {
                        heap.push_back(piece);
                        std::push_heap(heap.begin(), heap.end(), before);
                    }
                }
            };

            // The initial mesh of 10 equal panels; equal adjacent points (an interval too narrow for 10) give empty
            // panels, which are skipped. A mesh panel with no finite estimate is non_finite.
            Interval<T> const whole = detail::interval_of(lower, upper);
            T left = lower;
            for (int i = 1; i <= mesh; ++i) {
                T const u = T(T(-1) + T(T(2 * i) / T(mesh)));
                T const right = i == mesh ? upper : whole.point(T(whole.c + T(whole.hw * u)));
                if (left < right) {
                    pieces.clear();
                    if (resolve(left, right)) return std::unexpected(invalid);
                    commit();
                    if (refused) return std::unexpected(invalid);
                }
                left = right;
            }

            bool tainted = false; // a piece no descent avoided was undone: the run cannot end ok
            std::optional<T> least; // the least excess so far
            int since = 0;          // bisections since the excess last reached a new minimum
            while (true) {
                T const value = total.value();
                T const f = floor.value();
                T const x = excess.value();
                if (!detail::fits(f, x)) return std::unexpected(invalid);
                T const estimate = T(f + x);
                if (!nxx::math::isfinite(value) || !nxx::math::isfinite(estimate) || !nxx::math::isfinite(f))
                    return std::unexpected(invalid);
                // The result {K, E}, its estimate at least F: a Neumaier residue can leave E a hair below it.
                auto const result = [&] {
                    return Integral<T>{ reversed ? T(-value) : value, detail::larger(estimate, f) };
                };
                // Twice the floor: a kink or a jump never brings every panel to its floor (X = 0), and on a smooth fn
                // the test cannot accept an error above max(request, 2 F).
                if (estimate <= detail::wanted(abs_tol, rel_tol, value, T(T(2) * f))) {
                    if (tainted) return std::unexpected(Failure<T>{ error::not_converged, result() });
                    return result();
                }
                if (heap.empty() || work >= cap || since >= stagnation)
                    return std::unexpected(Failure<T>{ error::not_converged, result() });

                std::pop_heap(heap.begin(), heap.end(), before);
                Panel<T> const parent = heap.back();
                heap.pop_back();
                T const mid = detail::midpoint(parent.lo, parent.hi);
                ++work;
                pieces.clear();
                std::optional<fault> const lower_half = resolve(parent.lo, mid);
                if (lower_half == fault::overflow) return std::unexpected(invalid);
                std::optional<fault> const upper_half = resolve(mid, parent.hi);
                if (upper_half == fault::overflow) return std::unexpected(invalid);
                if (lower_half || upper_half) {
                    // Undone: the parent stays, final, with its own estimate; the descent's pieces are dropped.
                    tainted = true;
                }
                else {
                    contribute(parent, false);
                    commit();
                    if (refused) return std::unexpected(invalid);
                }
                T const now = excess.value();
                if (!least || now < *least) {
                    least = now;
                    since = 0;
                }
                else ++since;
            }
        }

    } // namespace detail
    /** @endcond */

    /**
     * @brief The integral of @p fn from @p a to @p b by the fixed rule @p M, GaussLegendre<N>: N evaluations of fn, a
     * plain number, no error estimate and no criteria.
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
     * @tparam M The rule, `GaussLegendre<N>` with 1 <= N <= 64.
     * @tparam T The real type, deduced from both ends, which must have the same type.
     * @tparam F A callable taking T const & and returning T (for a class type T, T itself or a built-in arithmetic
     *           type, not an expression template; for a built-in T anything convertible to T, so a double result for
     *           a float T narrows silently). The same for every method.
     * @param fn The function.
     * @param a The lower end of the interval.
     * @param b The upper end.
     * @return The integral; error::non_finite if an end or a value of fn is not finite, or the sum overflowed.
     */
    template<typename M, real T, typename F>
    requires detail::fixed_rule<M> && nxx::detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<T, error> integrate(F &&fn, T const &a, T const &b)
    {
        if (!nxx::math::isfinite(a) || !nxx::math::isfinite(b)) return std::unexpected(error::non_finite);
        if (a == b) return T(0);
        constexpr int N = detail::gauss_points<M>;
        // Qualified, so that argument-dependent lookup cannot pick a function from the namespace of T or F.
        detail::Interval<T> const interval = detail::interval_of(a, b);
        return detail::with_gauss_legendre<N, T>(
            [&](detail::Rule<T, N> const &rule) { return detail::apply(rule, fn, interval); });
    }

    /**
     * @brief The integral of @p fn from @p a to @p b by a method with an error estimate, @p M: GaussKronrod<15>, the
     * default, or Romberg.
     *
     * The method stops when its error estimate meets max(criteria.abs_tol, criteria.rel_tol |value|), or the rounding
     * floor (see Criteria). The estimate is not a bound (see Integral). fn's exceptions pass through unchanged.
     *
     * @code
     * namespace qd = nxx::quadrature;
     * auto const r = qd::integrate([](double x) { return std::exp(x); }, 0.0, 1.0);
     * auto const s = qd::integrate<qd::Romberg>([](double x) { return std::exp(x); }, 0.0, 1.0, {.rel_tol = 1e-12});
     * if (r) use(r->value, r->abs_error);
     * else if (r.error().best) report(r.error().best->value, r.error().best->abs_error);
     * @endcode
     *
     * @tparam M The method: `GaussKronrod<15>` (the default) or `Romberg`.
     * @tparam T The real type, deduced from both ends, which must have the same type.
     * @tparam F A callable taking T const & and returning T, as for the fixed rule.
     * @param fn The function.
     * @param a The lower end of the interval.
     * @param b The upper end.
     * @param criteria When to stop; it converts to T, so `{.rel_tol = 1e-12}` works for any T. Name the member: a
     *        braced `{1e-12}` sets abs_tol.
     * @return The integral and an estimate of its error, returned when the request is met or the estimate is at the
     *         floor, so a caller who needs the request compares abs_error with it; or a Failure: error::non_finite if
     *         an end, a tolerance or a value of fn is not finite, a tolerance is negative or a sum overflowed;
     *         error::not_converged, with the best estimate, if the tolerance was not met (see error).
     */
    // criteria is a non-deduced context (type_identity_t), so T comes from the ends alone and a braced list of doubles
    // converts to a multiprecision T.
    template<typename M = GaussKronrod<15>, real T, typename F>
    requires detail::estimating<M> && nxx::detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        Failure<T> const invalid{ error::non_finite, std::nullopt };
        // An assert would vanish under NDEBUG, and a NaN tolerance would then never stop: validated here.
        if (!nxx::math::isfinite(criteria.abs_tol) || !nxx::math::isfinite(criteria.rel_tol)
            || criteria.abs_tol < T(0) || criteria.rel_tol < T(0))
            return std::unexpected(invalid);
        // rel_tol |value| stays finite (constant evaluation rejects an overflow); abs_tol enters only a max of finite
        // values, so it needs no clamp, and the floor is the effective minimum of both.
        T const rel_tol = criteria.rel_tol > T(1) ? T(1) : criteria.rel_tol;
        if (!nxx::math::isfinite(a) || !nxx::math::isfinite(b)) return std::unexpected(invalid);
        if (a == b) return Integral<T>{ T(0), T(0) };
        if constexpr (std::same_as<M, Romberg>) return detail::romberg(fn, a, b, criteria.abs_tol, rel_tol);
        else
            return detail::with_gauss_kronrod<T>([&](detail::Kronrod<T> const &rule) {
                return detail::gauss_kronrod(rule, fn, a, b, criteria.abs_tol, rel_tol);
            });
    }

} // namespace nxx::quadrature

/** @brief The standard library's formatter, specialised for nxx::quadrature::error and nxx::quadrature::Failure. */
namespace std {

    /**
     * @brief Formats an nxx::quadrature::error as the name of its enumerator, such as `not_converged`, with the format
     * specification of a string.
     *
     * @code
     * if (!g) std::println("no integral: {}", g.error());   // no integral: non_finite
     * @endcode
     */
    template<>
    struct formatter<nxx::quadrature::error> : formatter<string_view>
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
        auto format(nxx::quadrature::error e, FormatContext &ctx) const
        { return formatter<string_view>::format(nxx::quadrature::detail::name(e), ctx); }
    };

    /**
     * @brief Formats an nxx::quadrature::Failure as the name of its reason, such as `not_converged`, with the format
     * specification of a string. The best estimate is not written: that would need a formatter for T, which
     * Boost.Multiprecision's types do not have; format `f.best->value` where T has one.
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
         * @brief Writes the name of @p f's reason.
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
