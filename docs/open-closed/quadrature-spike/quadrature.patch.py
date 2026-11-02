# Draft 2 of the quadrature method protocol, written into a copy of quadrature.hpp.
# Usage: patch2.py <in> <out> [normalize]   (normalize: integrate() orders the ends and negates; Q4's alternative)
import sys
src = open(sys.argv[1]).read()
normalize = len(sys.argv) > 3 and sys.argv[3] == "normalize"

def cut(start_marker, end_marker):
    global src
    i = src.index(start_marker); j = src.index(end_marker, i)
    block = src[i:j]; src = src[:i] + src[j:]
    return block

def sub(old, new):
    global src
    assert src.count(old) == 1, (old[:80], src.count(old))
    src = src.replace(old, new)

# 1. The three methods move from the top of the header to just before integrate(), with their members inline.
methods = cut("    /**\n     * @brief The N-point Gauss-Legendre rule", "    /** @cond INTERNAL */\n    namespace detail {")
methods = methods.replace("""    template<int N>
    struct GaussLegendre
    {
    };""", """    template<int N>
    struct GaussLegendre
    {
        /**
         * @brief The rule applied to fn from a to b (fixed_rule); integrate() calls it.
         * @param fn The function.
         * @param a One end, finite.
         * @param b The other end, finite, a != b.
         * @return The integral from a to b, or error::non_finite.
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
    };""")
methods = methods.replace("""    struct Romberg
    {
    };""", """    struct Romberg
    {
        /**
         * @brief Romberg's method from a to b (estimating_method); integrate() calls it.
         * @param fn The function.
         * @param a One end, finite.
         * @param b The other end, finite, a != b.
         * @param criteria The validated criteria.
         * @return The integral from a to b and its estimate, or the Failure.
         */
        template<real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            estimate(F &fn, T const &a, T const &b, Criteria<T> const &criteria) const
        { return detail::romberg(fn, a, b, criteria.abs_tol, criteria.rel_tol); }
    };""")
methods = methods.replace("""    template<int N>
    requires(N == 15)
    struct GaussKronrod
    {
    };""", """    template<int N>
    requires(N == 15)
    struct GaussKronrod
    {
        /**
         * @brief Adaptive Gauss-Kronrod from a to b (estimating_method); integrate() calls it.
         * @param fn The function.
         * @param a One end, finite.
         * @param b The other end, finite, a != b.
         * @param criteria The validated criteria.
         * @return The integral from a to b and its estimate, or the Failure.
         */
        template<real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            estimate(F &fn, T const &a, T const &b, Criteria<T> const &criteria) const
        {
            return detail::with_gauss_kronrod<T>([&](detail::Kronrod<T> const &rule) {
                return detail::gauss_kronrod(rule, fn, a, b, criteria.abs_tol, criteria.rel_tol);
            });
        }
    };""")

# 2. The closed lists go.
cut("        /**\n         * @brief The number of points of a GaussLegendre rule, or 0 for any other type.",
    "        /**\n         * @brief The types whose Gauss-Legendre rules are computed at compile time")

concepts = """    /**
     * @brief A fixed rule: M has a const member apply(fn, a, b) returning std::expected<T, error>.
     * @tparam M The method. @tparam T The real type. @tparam F fn's type, reference removed.
     */
    template<typename M, typename T, typename F = T(T const &)>
    concept fixed_rule = requires(M const &m, F &fn, T const &a, T const &b) {
        { m.apply(fn, a, b) } -> std::same_as<std::expected<T, error>>;
    };

    /**
     * @brief A method with an error estimate: M has a const member estimate(fn, a, b, criteria) returning
     * std::expected<Integral<T>, Failure<T>>.
     * @tparam M The method. @tparam T The real type. @tparam F fn's type, reference removed.
     */
    template<typename M, typename T, typename F = T(T const &)>
    concept estimating_method = requires(M const &m, F &fn, T const &a, T const &b, Criteria<T> const &criteria) {
        { m.estimate(fn, a, b, criteria) } -> std::same_as<std::expected<Integral<T>, Failure<T>>>;
    };

    /** @cond INTERNAL */
    namespace detail {

        /** @brief Whether an Integral is finite with a non-negative estimate (the protocol's checked part). */
        template<real T>
        constexpr bool sound(Integral<T> const &r)
        {
            return nxx::math::isfinite(r.value) && nxx::math::isfinite(r.abs_error) && !(r.abs_error < T(0));
        }

    } // namespace detail
    /** @endcond */

"""

if normalize:
    call_fixed = """        std::expected<T, error> const r =
            b < a ? m.apply(fn, b, a).transform([](T const &v) { return T(-v); }) : m.apply(fn, a, b);"""
    call_est = """        Criteria<T> const validated{ criteria.abs_tol, rel_tol };
        auto const negated = [](Integral<T> v) {
            v.value = T(-v.value);
            return v;
        };
        std::expected<Integral<T>, Failure<T>> const r = b < a
            ? m.estimate(fn, b, a, validated).transform(negated).transform_error([&](Failure<T> f) {
                  if (f.best) f.best = negated(*f.best);
                  return f;
              })
            : m.estimate(fn, a, b, validated);"""
else:
    call_fixed = """        std::expected<T, error> const r = m.apply(fn, a, b);"""
    call_est = """        std::expected<Integral<T>, Failure<T>> const r = m.estimate(fn, a, b, Criteria<T>{ criteria.abs_tol, rel_tol });"""

entries = concepts + methods + """    /** @brief A fixed rule given as a value: integrate(m, fn, a, b). */
    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T> && fixed_rule<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<T, error> integrate(M const &m, F &&fn, T const &a, T const &b)
    {
        if (!nxx::math::isfinite(a) || !nxx::math::isfinite(b)) return std::unexpected(error::non_finite);
        if (a == b) return T(0);
""" + call_fixed + """
        // The protocol's checked part: a value is finite.
        if (r && !nxx::math::isfinite(*r)) return std::unexpected(error::non_finite);
        return r;
    }

    /** @brief A method with an estimate given as a value: integrate(m, fn, a, b[, criteria]). */
    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T> && estimating_method<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(M const &m, F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        Failure<T> const invalid{ error::non_finite, std::nullopt };
        if (!nxx::math::isfinite(criteria.abs_tol) || !nxx::math::isfinite(criteria.rel_tol)
            || criteria.abs_tol < T(0) || criteria.rel_tol < T(0))
            return std::unexpected(invalid);
        T const rel_tol = criteria.rel_tol > T(1) ? T(1) : criteria.rel_tol;
        if (!nxx::math::isfinite(a) || !nxx::math::isfinite(b)) return std::unexpected(invalid);
        if (a == b) return Integral<T>{ T(0), T(0) };
""" + call_est + """
        // The protocol's checked part: finite values and estimates, a best only with not_converged.
        if (r ? !detail::sound(*r)
              : r.error().best && (r.error().reason != error::not_converged || !detail::sound(*r.error().best)))
            return std::unexpected(invalid);
        return r;
    }

    /** @brief The fixed rule M, named: integrate<M>(fn, a, b). */
    template<typename M, real T, typename F>
    requires std::default_initializable<M> && nxx::detail::callable<F, T>
             && fixed_rule<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<T, error> integrate(F &&fn, T const &a, T const &b)
    { return nxx::quadrature::integrate(M{}, fn, a, b); }

    /** @brief The method M with an estimate, named: integrate<M>(fn, a, b[, criteria]); GaussKronrod<15> by default. */
    template<typename M = GaussKronrod<15>, real T, typename F>
    requires std::default_initializable<M> && nxx::detail::callable<F, T>
             && estimating_method<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return nxx::quadrature::integrate(M{}, fn, a, b, criteria); }

} // namespace nxx::quadrature
"""
i = src.index("    /**\n     * @brief The integral of @p fn from @p a to @p b by the fixed rule @p M")
j = src.index("} // namespace nxx::quadrature\n", i) + len("} // namespace nxx::quadrature\n")
src = src[:i] + entries + src[j:]
open(sys.argv[2], "w").write(src)
