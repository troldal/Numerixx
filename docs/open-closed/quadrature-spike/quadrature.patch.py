import sys, re
src = open(sys.argv[1]).read()

def sub(old, new, count=1):
    global src
    assert src.count(old) == count, (old[:80], src.count(old))
    src = src.replace(old, new)

# 1. The methods declare their member.
sub("""    template<int N>
    struct GaussLegendre
    {
    };""", """    template<int N>
    struct GaussLegendre
    {
        /**
         * @brief The rule applied to fn on [a, b] (the fixed_rule protocol); integrate() calls it.
         * @param fn The function.
         * @param a One end, finite.
         * @param b The other end, finite, a != b.
         * @return The integral from a to b, or error::non_finite.
         */
        template<real T, typename F>
        requires(N >= 1 && N <= 64)
        constexpr std::expected<T, error> apply(F &fn, T const &a, T const &b) const;
    };""")
sub("""    struct Romberg
    {
    };""", """    struct Romberg
    {
        /**
         * @brief Romberg's method on [a, b] (the estimating_method protocol); integrate() calls it.
         * @param fn The function.
         * @param a One end, finite.
         * @param b The other end, finite, a != b.
         * @param criteria The validated criteria: both finite and non-negative, rel_tol at most 1.
         * @return The integral from a to b and its estimate, or the Failure.
         */
        template<real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            estimate(F &fn, T const &a, T const &b, Criteria<T> const &criteria) const;
    };""")
sub("""    template<int N>
    requires(N == 15)
    struct GaussKronrod
    {
    };""", """    template<int N>
    requires(N == 15)
    struct GaussKronrod
    {
        /**
         * @brief Adaptive Gauss-Kronrod on [a, b] (the estimating_method protocol); integrate() calls it.
         * @param fn The function.
         * @param a One end, finite.
         * @param b The other end, finite, a != b.
         * @param criteria The validated criteria: both finite and non-negative, rel_tol at most 1.
         * @return The integral from a to b and its estimate, or the Failure.
         */
        template<real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            estimate(F &fn, T const &a, T const &b, Criteria<T> const &criteria) const;
    };

    /**
     * @brief A fixed rule: a method M whose const member apply(fn, a, b) returns std::expected<T, error>, the
     * integral from a to b. integrate() validates the ends first (finite, a != b) and passes them in the caller's
     * order; the rule returns a finite value or error::non_finite.
     */
    template<typename M, typename F, typename T>
    concept fixed_rule = std::copy_constructible<M> && requires(M const &m, F &fn, T const &a, T const &b) {
        { m.apply(fn, a, b) } -> std::same_as<std::expected<T, error>>;
    };

    /**
     * @brief A method with an error estimate: a method M whose const member estimate(fn, a, b, criteria) returns
     * std::expected<Integral<T>, Failure<T>>. integrate() validates the ends (finite, a != b, in the caller's order)
     * and the criteria (finite, non-negative, rel_tol clamped to 1) first.
     */
    template<typename M, typename F, typename T>
    concept estimating_method = std::copy_constructible<M>
                                && requires(M const &m, F &fn, T const &a, T const &b, Criteria<T> const &c) {
                                       { m.estimate(fn, a, b, c) } -> std::same_as<std::expected<Integral<T>, Failure<T>>>;
                                   };""")

# 2. The closed lists go.
start = src.index("        /**\n         * @brief The number of points of a GaussLegendre rule, or 0 for any other type.")
end = src.index("        /**\n         * @brief The types whose Gauss-Legendre rules are computed at compile time")
src = src[:start] + src[end:]

# 3. The entry points dispatch through the protocol; the members are defined after detail.
old_fixed = src[src.index("    template<typename M, real T, typename F>\n    requires detail::fixed_rule<M>"):]
old_fixed = old_fixed[:old_fixed.index("    /**")]
sub(old_fixed, """    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T> && fixed_rule<M, std::remove_reference_t<F>, T>
    [[nodiscard]] constexpr std::expected<T, error> integrate(F &&fn, T const &a, T const &b, M const &m = M{})
    {
        if (!nxx::math::isfinite(a) || !nxx::math::isfinite(b)) return std::unexpected(error::non_finite);
        if (a == b) return T(0);
        std::expected<T, error> const r = m.apply(fn, a, b);
        // The protocol's postcondition: a finite value or error::non_finite.
        assert(r ? nxx::math::isfinite(*r) : r.error() == error::non_finite);
        return r;
    }

""")
old_est = src[src.index("    template<typename M = GaussKronrod<15>, real T, typename F>\n    requires detail::estimating<M>"):]
old_est = old_est[:old_est.index("} // namespace nxx::quadrature")]
sub(old_est, """    template<typename M = GaussKronrod<15>, real T, typename F>
    requires nxx::detail::callable<F, T> && estimating_method<M, std::remove_reference_t<F>, T>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>> integrate(F &&fn,
        T const &a,
        T const &b,
        Criteria<std::type_identity_t<T>> const &criteria = {},
        M const &m = M{})
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
        std::expected<Integral<T>, Failure<T>> const r = m.estimate(fn, a, b, Criteria<T>{ criteria.abs_tol, rel_tol });
        // The protocol's postcondition: a finite value and a finite, non-negative estimate, in the result and in a
        // best; no best with error::non_finite.
        assert(r ? detail::sound(*r)
                 : (r.error().best ? r.error().reason == error::not_converged && detail::sound(*r.error().best)
                                   : true));
        return r;
    }

""")

# 4. The members' definitions and the postcondition helper, after detail and before the entry points.
anchor = "    /**\n     * @brief The integral of @p fn from @p a to @p b by the fixed rule @p M"
sub(anchor, """    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief Whether an Integral meets the protocol: a finite value, a finite and non-negative abs_error.
         * @param r The result.
         * @return true if it does.
         */
        template<real T>
        constexpr bool sound(Integral<T> const &r)
        {
            return nxx::math::isfinite(r.value) && nxx::math::isfinite(r.abs_error) && !(r.abs_error < T(0));
        }

    } // namespace detail
    /** @endcond */

    template<int N>
    template<real T, typename F>
    requires(N >= 1 && N <= 64)
    constexpr std::expected<T, error> GaussLegendre<N>::apply(F &fn, T const &a, T const &b) const
    {
        // Qualified, so that argument-dependent lookup cannot pick a function from the namespace of T or F.
        detail::Interval<T> const interval = detail::interval_of(a, b);
        return detail::with_gauss_legendre<N, T>(
            [&](detail::Rule<T, N> const &rule) { return detail::apply(rule, fn, interval); });
    }

    template<real T, typename F>
    constexpr std::expected<Integral<T>, Failure<T>>
        Romberg::estimate(F &fn, T const &a, T const &b, Criteria<T> const &criteria) const
    { return detail::romberg(fn, a, b, criteria.abs_tol, criteria.rel_tol); }

    template<int N>
    requires(N == 15)
    template<real T, typename F>
    constexpr std::expected<Integral<T>, Failure<T>>
        GaussKronrod<N>::estimate(F &fn, T const &a, T const &b, Criteria<T> const &criteria) const
    {
        return detail::with_gauss_kronrod<T>([&](detail::Kronrod<T> const &rule) {
            return detail::gauss_kronrod(rule, fn, a, b, criteria.abs_tol, criteria.rel_tol);
        });
    }

""" + anchor)
open(sys.argv[2], "w").write(src)
