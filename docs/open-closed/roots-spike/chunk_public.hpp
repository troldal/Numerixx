
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
         * @return absolute + relative |x|, or max if that would overflow.
         */
        constexpr T at(T const &x) const
        {
            // scaled is at most |x|, since relative <= 1. Above max / 2 it saturates: max - absolute would round,
            // and could admit a sum that overflows; below it the sum is at most 3/4 max, since absolute <= max / 4.
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
     * bracket around 0, at 0 (zero first). A method that needs fn's derivatives (values_of<M> > 1) receives them too:
     * at_lower[0] is fn's value, at_lower[1] f', at_lower[2] f''. The derivatives are not checked.
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
     * For a bracket, lower < upper with fn's values of opposite signs, inside the caller's bracket, and strictly inside
     * the previous bracket. solve() returns point once the bracket is at most 2 tol(point) wide or point is within
     * tol(point) of both ends, with abs_error the larger distance from point to an end, and point with abs_error
     * tol(point) where point is an end at which fn is exactly 0. From a guess, before a sign change, lower and upper are
     * the last two iterates (equal after the start), and point the latest. A method that cannot go on returns
     * error::not_converged from advance.
     *
     * @tparam T The real type.
     */
    template<real T>
    struct Step
    {
        /** @brief The lower point. */
        T lower;

        /** @brief fn's value at lower. */
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

    /**
     * @brief How many values a method needs from fn at a point: M::values where M declares it (Newton 2, Halley 3),
     * else 1, fn's value alone.
     * @tparam M The method.
     */
    // Only a static integral constant counts: a method with a data member that happens to be called values (a table of
    // coefficients) needs one value (numerics review of draft 1, finding 7).
    template<typename M>
    inline constexpr std::size_t values_of = [] {
        if constexpr (requires {
                          requires std::integral<std::remove_cvref_t<decltype(M::values)>>;
                          typename std::integral_constant<std::size_t, M::values>;
                      })
            return std::size_t{ M::values };
        else
            return std::size_t{ 1 };
    }();
