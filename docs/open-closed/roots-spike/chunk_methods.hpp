
    /**
     * @brief Bisection: halves the bracket at every step, one bit per evaluation of fn, guaranteed.
     *
     * About log2(width / (2 tol)) evaluations, 52 for [1, 2] in double, whatever fn is.
     */
    struct Bisection
    {
        /** @brief The state on a validated bracket. @param s The bracket. @return The state. */
        template<real T, typename F>
        constexpr std::expected<detail::BisectionState<T>, error> start(F &, Start<T> const &s, Tolerance<T> const &) const
        {
            T const fa = s.at_lower[0];
            T const fb = s.at_upper[0];
            return detail::BisectionState<T>{ s.lower, fa, s.upper, fb, detail::bisection_step(s.lower, fa, s.upper, fb) };
        }

        /** @brief One halving. @param fn The function. @param s The state. @return Nothing, or error::non_finite. */
        template<real T, typename F>
        constexpr std::expected<void, error> advance(F &fn, detail::BisectionState<T> &s, Tolerance<T> const &) const
        {
            T const m = s.step.point;
            // Not reached where the stop test failed (PHASE7_API.md B).
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

    /** @brief Brent's method (Brent 1973). */
    struct Brent
    {
        /** @brief The state on a validated bracket. @param s The bracket. @return The state. */
        template<real T, typename F>
        constexpr std::expected<detail::SafeguardedState<T, 1>, error>
            start(F &, Start<T> const &s, Tolerance<T> const &) const
        { return detail::safeguarded_start(s); }

        /** @brief One iteration. @param fn The function. @param s The state. @param t The tolerances. @return Nothing, or an error. */
        template<real T, typename F>
        constexpr std::expected<void, error>
            advance(F &fn, detail::SafeguardedState<T, 1> &s, Tolerance<T> const &t) const
        { return detail::safeguarded_advance<detail::Safeguard::interpolation>(fn, s, t); }
    };

    /** @brief Alefeld, Potra and Shi's Algorithm 748 (1995, Algorithm 4.2). A bracket only. */
    struct Toms748
    {
        /** @brief The state on a validated bracket. @param s The bracket. @return The state. */
        template<real T, typename F>
        constexpr std::expected<detail::Toms748State<T>, error> start(F &, Start<T> const &s, Tolerance<T> const &) const
        {
            T const fa = s.at_lower[0];
            T const fb = s.at_upper[0];
            // d and e are read only once steps have set them (the first iteration uses Newton-Quadratic, without e).
            return detail::Toms748State<T>{ detail::Enclosure<T>{ s.lower, fa, s.upper, fb, s.upper, fb, s.upper, fb },
                s.lower, s.upper, 0, true, detail::bracket_step(s.lower, fa, s.upper, fb) };
        }

        /** @brief One evaluation. @param fn The function. @param s The state. @param t The tolerances. @return Nothing, or an error. */
        template<real T, typename F>
        constexpr std::expected<void, error> advance(F &fn, detail::Toms748State<T> &s, Tolerance<T> const &t) const
        { return detail::toms748_advance(fn, s, t); }
    };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief False position with a scaling of the retained end's value: the members of AndersonBjorck, Illinois
         * and RegulaFalsi.
         * @tparam S The scaling.
         */
        template<Scaling S>
        struct FalsePositionMethod
        {
            /** @brief The state on a validated bracket. @param s The bracket. @return The state. */
            template<real T, typename F>
            constexpr std::expected<FalsePositionState<T>, error> start(F &, Start<T> const &s, Tolerance<T> const &) const
            {
                T const fa = s.at_lower[0];
                T const fb = s.at_upper[0];
                return FalsePositionState<T>{ s.lower, fa, fa, s.upper, fb, detail::distance_bound(s.lower, s.upper), 0,
                    detail::bracket_step(s.lower, fa, s.upper, fb) };
            }

            /** @brief One evaluation. @param fn The function. @param s The state. @param t The tolerances. @return Nothing, or an error. */
            template<real T, typename F>
            constexpr std::expected<void, error> advance(F &fn, FalsePositionState<T> &s, Tolerance<T> const &t) const
            { return detail::false_position_advance<S>(fn, s, t); }
        };

        /**
         * @brief A method with a derivative step (Newton's or Halley's): on a bracket inside Brent's safeguard, and
         * from a guess until a sign change, then on that bracket.
         * @tparam G The step.
         * @tparam K The values it needs.
         */
        template<Safeguard G, std::size_t K>
        struct DerivativeMethod
        {
            /** @brief fn's values the method needs: the value and K - 1 derivatives. */
            static constexpr std::size_t values = K;

            /** @brief The state on a validated bracket. @param s The bracket. @return The state. */
            template<real T, typename F>
            constexpr std::expected<SafeguardedState<T, K>, error> start(F &, Start<T, K> const &s, Tolerance<T> const &) const
            { return detail::safeguarded_start(s); }

            /** @brief One iteration on a bracket. @param fn The function. @param s The state. @param t The tolerances. @return Nothing, or an error. */
            template<real T, typename F>
            constexpr std::expected<void, error> advance(F &fn, SafeguardedState<T, K> &s, Tolerance<T> const &t) const
            { return detail::safeguarded_advance<G>(fn, s, t); }

            /** @brief The state from a guess. @param fn The function. @param guess The guess, finite. @return The state, or error::non_finite. */
            template<real T, typename F>
            constexpr std::expected<OpenDerivativeState<T, K>, error> start(F &fn, T const &guess, Tolerance<T> const &) const
            {
                std::array<T, K> const v = detail::evaluate<K>(fn, guess);
                if (!math::isfinite(v[0])) return std::unexpected(error::non_finite);
                return OpenDerivativeState<T, K>{ guess, v, std::nullopt, Step<T>{ guess, v[0], guess, v[0], guess } };
            }

            /** @brief One iteration from a guess. @param fn The function. @param s The state. @param t The tolerances. @return Nothing, or an error. */
            template<real T, typename F>
            constexpr std::expected<void, error> advance(F &fn, OpenDerivativeState<T, K> &s, Tolerance<T> const &t) const
            { return detail::open_derivative_advance<G>(fn, s, t); }
        };

    } // namespace detail
    /** @endcond */

    /** @brief Regula falsi with Anderson and Bjorck's rescaling (1973). A bracket only. */
    struct AndersonBjorck : detail::FalsePositionMethod<detail::Scaling::anderson_bjorck>
    {
    };

    /** @brief False position, guarded. A bracket. */
    struct RegulaFalsi : detail::FalsePositionMethod<detail::Scaling::none>
    {
    };

    /** @brief The Illinois method (Dowell and Jarratt 1971). A bracket. */
    struct Illinois : detail::FalsePositionMethod<detail::Scaling::illinois>
    {
    };

    /** @brief Ridders' method (1979). A bracket. */
    struct Ridders
    {
        /** @brief The state on a validated bracket. @param s The bracket. @return The state. */
        template<real T, typename F>
        constexpr std::expected<detail::RiddersState<T>, error> start(F &, Start<T> const &s, Tolerance<T> const &) const
        {
            T const fa = s.at_lower[0];
            T const fb = s.at_upper[0];
            return detail::RiddersState<T>{ detail::Enclosure<T>{ s.lower, fa, s.upper, fb, s.upper, fb, s.upper, fb },
                s.lower, fa, 0, detail::bracket_step(s.lower, fa, s.upper, fb) };
        }

        /** @brief One evaluation. @param fn The function. @param s The state. @param t The tolerances. @return Nothing, or an error. */
        template<real T, typename F>
        constexpr std::expected<void, error> advance(F &fn, detail::RiddersState<T> &s, Tolerance<T> const &t) const
        { return detail::ridders_advance(fn, s, t); }
    };

    /** @brief Newton's method: the step -f / f'. Needs fn and f'. */
    struct Newton : detail::DerivativeMethod<detail::Safeguard::newton, 2>
    {
    };

    /** @brief Halley's method. Needs fn, f' and f''. */
    struct Halley : detail::DerivativeMethod<detail::Safeguard::halley, 3>
    {
    };

    /** @brief The secant method. Needs fn only; from a guess. */
    struct Secant
    {
        /** @brief The state from a guess: fn at the guess and at the second point. @return The state, or error::non_finite. */
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

        /** @brief One iteration. @param fn The function. @param s The state. @param t The tolerances. @return Nothing, or an error. */
        template<real T, typename F>
        constexpr std::expected<void, error> advance(F &fn, detail::SecantState<T> &s, Tolerance<T> const &t) const
        { return detail::secant_advance(fn, s, t); }
    };

    /** @brief Steffensen's method. Needs fn only; from a guess. */
    struct Steffensen
    {
        /** @brief The state from a guess. @return The state, or error::non_finite. */
        template<real T, typename F>
        constexpr std::expected<detail::SteffensenState<T>, error>
            start(F &fn, T const &guess, Tolerance<T> const &) const
        {
            T const fx = detail::evaluate<1>(fn, guess)[0];
            if (!math::isfinite(fx)) return std::unexpected(error::non_finite);
            return detail::SteffensenState<T>{ guess, fx, std::nullopt, Step<T>{ guess, fx, guess, fx, guess } };
        }

        /** @brief One iteration. @param fn The function. @param s The state. @param t The tolerances. @return Nothing, or an error. */
        template<real T, typename F>
        constexpr std::expected<void, error> advance(F &fn, detail::SteffensenState<T> &s, Tolerance<T> const &t) const
        { return detail::steffensen_advance(fn, s, t); }
    };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief The type of fn the concepts assume by default: a plain function of T returning the values a method
         * needs (T, or a std::array of K values of T).
         */
        template<typename T, std::size_t K>
        using archetype = std::conditional_t<K == 1, T (*)(T const &), std::array<T, K> (*)(T const &)>;

        /** @brief What a method's start returns from a start X (a Start or a guess). */
        template<typename M, typename T, typename F, typename X>
        using start_result = decltype(std::declval<M const &>().start(
            std::declval<F &>(), std::declval<X const &>(), std::declval<Tolerance<T> const &>()));

        /** @brief Whether R is std::expected<S, error> for some state S. */
        template<typename R>
        inline constexpr bool is_state_result = false;

        /** @brief It is. */
        template<typename S>
        inline constexpr bool is_state_result<std::expected<S, error>> = !std::is_reference_v<S>;

        /** @brief A start's result. */
        template<typename R>
        concept state_result = is_state_result<R>;

    } // namespace detail
    /** @endcond */

    /** @brief A bracketing method's state for fn of type F: what its start returns, inside its std::expected. */
    template<typename M, typename T, typename F = detail::archetype<T, values_of<M>>>
    using bracketing_state_t = typename detail::start_result<M, T, F, Start<T, values_of<M>>>::value_type;

    /** @brief An open method's state for fn of type F: what its start from a guess returns, inside its std::expected. */
    template<typename M, typename T, typename F = detail::archetype<T, values_of<M>>>
    using open_state_t = typename detail::start_result<M, T, F, T>::value_type;

    /**
     * @brief What a type must provide to work with solve() on a bracket, one step at a time: `start` and `advance`,
     * each reporting a Step; solve() checks each step and decides when to stop.
     */
    template<typename M, typename T, typename F = detail::archetype<T, values_of<M>>>
    concept bracketing_method =
        requires(M const &m, F &fn, Start<T, values_of<M>> const &s, Tolerance<T> const &t) {
            { m.start(fn, s, t) } -> detail::state_result;
        } && requires(M const &m, F &fn, bracketing_state_t<M, T, F> &st, Tolerance<T> const &t) {
            { m.advance(fn, st, t) } -> std::same_as<std::expected<void, error>>;
            { std::as_const(st).step } -> std::same_as<Step<T> const &>;
            // advance updates the state: one taken by value or by const reference would loop for ever on a copy
            // (usage review of draft 1, finding 2).
            requires !requires { m.advance(fn, std::move(st), t); };
        };

    /** @brief What a type must provide to work with solve() from a guess, one step at a time. */
    template<typename M, typename T, typename F = detail::archetype<T, values_of<M>>>
    concept open_method = requires(M const &m, F &fn, T const &guess, Tolerance<T> const &t) {
        { m.start(fn, guess, t) } -> detail::state_result;
    } && requires(M const &m, F &fn, open_state_t<M, T, F> &st, Tolerance<T> const &t) {
        { m.advance(fn, st, t) } -> std::same_as<std::expected<void, error>>;
        { std::as_const(st).step } -> std::same_as<Step<T> const &>;
        // advance updates the state: one taken by value or by const reference would loop for ever on a copy
        // (usage review of draft 1, finding 2).
        requires !requires { m.advance(fn, std::move(st), t); };
    };

    /** @brief A method that solves on a bracket in one call: `m.solve(fn, start, tolerance)`. */
    template<typename M, typename T, typename F = detail::archetype<T, values_of<M>>>
    concept bracketing_solver = requires(M const &m, F &fn, Start<T, values_of<M>> const &s, Tolerance<T> const &t) {
        { m.solve(fn, s, t) } -> std::same_as<std::expected<Root<T>, error>>;
    };

    /** @brief A method that solves from a guess in one call: `m.solve(fn, guess, tolerance)`. */
    template<typename M, typename T, typename F = detail::archetype<T, values_of<M>>>
    concept open_solver = requires(M const &m, F &fn, T const &guess, Tolerance<T> const &t) {
        { m.solve(fn, guess, t) } -> std::same_as<std::expected<Root<T>, error>>;
    };

    /** @cond INTERNAL */
    namespace detail {

        /** @brief M works on a bracket, stepped or in one call. */
        template<typename M, typename T, typename F>
        concept from_bracket = bracketing_method<M, T, F> || bracketing_solver<M, T, F>;

        /** @brief M works from a guess, stepped or in one call. */
        template<typename M, typename T, typename F>
        concept from_guess = open_method<M, T, F> || open_solver<M, T, F>;

        /** @brief Whether O is an Observer. */
        template<typename O>
        inline constexpr bool is_observer = false;

        /** @brief It is. */
        template<typename O>
        inline constexpr bool is_observer<Observer<O>> = true;

        /**
         * @brief The method solve() runs when it is given an observer: M wrapped so that each step is reported.
         * @tparam M The method.
         * @tparam O The observer's callable.
         */
        template<typename M, typename O>
        struct Observed
        {
            /** @brief The method. */
            M const &method;

            /** @brief The observer. */
            mutable nxx::detail::kept<O> observer;

            /** @brief The method's values. */
            static constexpr std::size_t values = values_of<M>;

            /** @brief The method's start on a bracket, then a report. */
            template<real T, std::size_t K, typename F, typename Self = M>
            requires bracketing_method<Self, T, F>
            constexpr auto start(F &fn, Start<T, K> const &s, Tolerance<T> const &t) const
                -> decltype(std::declval<Self const &>().start(fn, s, t))
            {
                auto started = method.start(fn, s, t);
                if (started) report(started->step);
                return started;
            }

            /** @brief The method's start from a guess, then a report. */
            template<real T, typename F, typename Self = M>
            requires open_method<Self, T, F>
            constexpr auto start(F &fn, T const &guess, Tolerance<T> const &t) const
                -> decltype(std::declval<Self const &>().start(fn, guess, t))
            {
                auto started = method.start(fn, guess, t);
                if (started) report(started->step);
                return started;
            }

            /** @brief The method's step, then a report. */
            template<typename F, typename S, real T, typename Self = M>
            constexpr auto advance(F &fn, S &s, Tolerance<T> const &t) const
                -> decltype(std::declval<Self const &>().advance(fn, s, t))
            {
                auto const r = method.advance(fn, s, t);
                if (r) report(s.step);
                return r;
            }

            /** @brief Calls the observer with a step. */
            template<real T>
            constexpr void report(Step<T> const &step) const
            {
                constexpr bool callable = std::invocable<nxx::detail::kept<O> &, Step<T> const &>;
                static_assert(callable, "nxx::roots: observer(o): o must be callable as o(Step<T> const &)");
                if constexpr (callable) {
                    constexpr bool returns_void =
                        std::is_void_v<std::invoke_result_t<nxx::detail::kept<O> &, Step<T> const &>>;
                    static_assert(returns_void,
                        "nxx::roots: observer(o): o must return void (a result is reserved for stopping the run)");
                    if constexpr (returns_void) std::invoke(observer, step);
                }
            }
        };

        /** @brief m observed by o. */
        template<typename M, typename O>
        constexpr Observed<M, O> observed(M const &m, Observer<O> &&o)
        { return Observed<M, O>{ m, std::move(o.callable) }; }

        /** @brief The caller's bracket after validation, inside which every step and result must lie. */
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
         * advance. Before an open method's first sign change, the steps count against the iteration cap.
         * @param m The method.
         * @param fn The function.
         * @param state The method's state, started.
         * @param tolerance The tolerances.
         * @param within The caller's bracket, for a bracketing method; nothing from a guess.
         * @return The root, or an error.
         */
        template<typename M, real T, typename F, typename S>
        constexpr std::expected<Root<T>, error> run(M const &m,
            F &fn,
            S &state,
            Tolerance<T> const &tolerance,
            std::optional<Within<T>> const &within)
        {
            std::optional<Within<T>> previous; // the last step's bracket, where it was one
            int approaches = 0;
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
                    // Strictly nested in the last bracket: finitely many such steps, so every method ends
                    // (numerics review of draft 1, finding 3).
                    if (previous
                        && (s.lower < previous->lower || previous->upper < s.upper
                            || (s.lower == previous->lower && s.upper == previous->upper)))
                        return std::unexpected(error::invalid_step);
                    previous = Within<T>{ s.lower, s.upper };
                    T const tol = tolerance.at(s.point);
                    T const below = detail::distance_bound(s.point, s.lower);
                    T const above = detail::distance_bound(s.upper, s.point);
                    T const larger = below < above ? above : below;
                    // The bracket at most 2 tol wide (Criteria's rule), or the point within tol of both ends, which
                    // the width can miss by the unit distance_bound adds (numerics review of draft 1, finding 2).
                    if (detail::distance_bound(s.lower, s.upper) <= detail::twice(tol) || larger <= tol)
                        return Root<T>{ s.point, larger };
                } else {
                    // A bracketing method must keep its sign change; from a guess, a step without one counts against
                    // the cap, so that a method that finds and loses sign changes still ends.
                    if (within) return std::unexpected(error::invalid_step);
                    previous.reset();
                    if (approaches == detail::iteration_cap<T>()) return std::unexpected(error::not_converged);
                    ++approaches;
                }
                std::expected<void, error> const r = m.advance(fn, state, tolerance);
                if (!r) return std::unexpected(r.error());
            }
        }

        /**
         * @brief A one-call method's result, checked: a finite value inside the caller's bracket and a finite,
         * non-negative abs_error.
         * @param r The result.
         * @param within The caller's bracket, or nothing.
         * @param tolerance The tolerances.
         * @return r, or error::non_finite or error::invalid_step.
         */
        template<real T>
        constexpr std::expected<Root<T>, error> checked(std::expected<Root<T>, error> const &r,
            Tolerance<T> const &tolerance,
            std::optional<Within<T>> const &within)
        {
            if (!r) return r;
            if (!math::isfinite(r->value) || !math::isfinite(r->abs_error)) return std::unexpected(error::non_finite);
            if (r->abs_error < T(0) || detail::twice(tolerance.at(r->value)) < r->abs_error)
                return std::unexpected(error::invalid_step); // Root's bound (numerics review of draft 1, finding 4)
            if (within && (r->value < within->lower || within->upper < r->value))
                return std::unexpected(error::invalid_step);
            return r;
        }

        /** @brief solve() on a bracket: validation, zero first, then the method. */
        template<typename M, real T, typename F>
        constexpr std::expected<Root<T>, error>
            on_bracket(M const &m, F &fn, Bracket<T> const &bracket, Criteria<T> const &criteria)
        {
            constexpr std::size_t K = values_of<M>;
            std::expected<Tolerance<T>, error> const tolerance = detail::checked_tolerance(criteria);
            if (!tolerance) return std::unexpected(tolerance.error());
            std::expected<Prepared<T, K>, error> const prepared = detail::prepare<K>(fn, bracket, *tolerance);
            if (!prepared) return std::unexpected(prepared.error());
            if (prepared->root) return *prepared->root;
            Start<T, K> const &start = prepared->start;
            std::optional<Within<T>> const within = Within<T>{ start.lower, start.upper };
            if constexpr (bracketing_method<M, T, F>) {
                auto state = m.start(fn, start, *tolerance);
                if (!state) return std::unexpected(state.error());
                return detail::run(m, fn, *state, *tolerance, within);
            } else {
                return detail::checked(m.solve(fn, start, *tolerance), *tolerance, within);
            }
        }

        /** @brief solve() from a guess: validation, then the method. */
        template<typename M, real T, typename F>
        constexpr std::expected<Root<T>, error>
            on_guess(M const &m, F &fn, T const &guess, Criteria<T> const &criteria)
        {
            std::expected<Tolerance<T>, error> const tolerance = detail::checked_tolerance(criteria);
            if (!tolerance) return std::unexpected(tolerance.error());
            if (!math::isfinite(guess)) return std::unexpected(error::non_finite);
            if constexpr (open_method<M, T, F>) {
                auto state = m.start(fn, guess, *tolerance);
                if (!state) return std::unexpected(state.error());
                return detail::run(m, fn, *state, *tolerance, std::optional<Within<T>>());
            } else {
                return detail::checked(m.solve(fn, guess, *tolerance), *tolerance, std::optional<Within<T>>());
            }
        }

    } // namespace detail
    /** @endcond */

    /** @brief A root of fn in bracket by the method M, named. */
    template<typename M, real T, typename F>
    requires std::default_initializable<M> && detail::from_bracket<M, T, std::remove_reference_t<F>>
             && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_bracket(M{}, fn, bracket, criteria); }

    /** @brief A root of fn from guess by the method M, named. */
    template<typename M, real T, typename F>
    requires std::default_initializable<M> && detail::from_guess<M, T, std::remove_reference_t<F>>
             && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_guess(M{}, fn, guess, criteria); }

    /** @brief A root of fn in bracket by the method m, given as a value. */
    template<typename M, real T, typename F>
    requires(!detail::is_observer<M>) && detail::from_bracket<M, T, std::remove_reference_t<F>>
            && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(M const &m, F &&fn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_bracket(m, fn, bracket, criteria); }

    /** @brief A root of fn from guess by the method m, given as a value. */
    template<typename M, real T, typename F>
    requires(!detail::is_observer<M>) && detail::from_guess<M, T, std::remove_reference_t<F>>
            && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(M const &m, F &&fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_guess(m, fn, guess, criteria); }

    /** @brief As solve<M>(fn, bracket), each step reported to the observer. */
    template<typename M, typename O, real T, typename F>
    requires std::default_initializable<M> && bracketing_method<M, T, std::remove_reference_t<F>>
             && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(Observer<O> o, F &&fn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_bracket(detail::observed(M{}, std::move(o)), fn, bracket, criteria); }

    /** @brief As solve<M>(fn, guess), each step reported to the observer. */
    template<typename M, typename O, real T, typename F>
    requires std::default_initializable<M> && open_method<M, T, std::remove_reference_t<F>>
             && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(Observer<O> o, F &&fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_guess(detail::observed(M{}, std::move(o)), fn, guess, criteria); }

    /** @brief As solve(m, fn, bracket), each step reported to the observer. */
    template<typename M, typename O, real T, typename F>
    requires bracketing_method<M, T, std::remove_reference_t<F>> && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error> solve(M const &m,
        Observer<O> o,
        F &&fn,
        Bracket<T> const &bracket,
        Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_bracket(detail::observed(m, std::move(o)), fn, bracket, criteria); }

    /** @brief As solve(m, fn, guess), each step reported to the observer. */
    template<typename M, typename O, real T, typename F>
    requires open_method<M, T, std::remove_reference_t<F>> && detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error> solve(M const &m,
        Observer<O> o,
        F &&fn,
        T const &guess,
        Criteria<std::type_identity_t<T>> const &criteria = {})
    { return detail::on_guess(detail::observed(m, std::move(o)), fn, guess, criteria); }

    /** @brief A method with values_of<M> == 2 on bracket, the derivative as a separate callable. */
    template<typename M, real T, typename F, typename D>
    requires(values_of<M> == 2) && detail::value_valued<F, T> && detail::derivative_valued<D, T>
            && std::default_initializable<M>
            && detail::from_bracket<M, T, detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, D &&dfn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const joined{ std::addressof(fn),
            { std::addressof(dfn) } };
        return nxx::roots::solve<M>(joined, bracket, criteria);
    }

    /** @brief A method with values_of<M> == 2 from guess, the derivative as a separate callable. */
    template<typename M, real T, typename F, typename D>
    requires(values_of<M> == 2) && detail::value_valued<F, T> && detail::derivative_valued<D, T>
            && std::default_initializable<M>
            && detail::from_guess<M, T, detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, D &&dfn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const joined{ std::addressof(fn),
            { std::addressof(dfn) } };
        return nxx::roots::solve<M>(joined, guess, criteria);
    }

    /** @brief A method with values_of<M> == 3 on bracket, the derivatives as separate callables. */
    template<typename M, real T, typename F, typename D, typename D2>
    requires(values_of<M> == 3) && detail::value_valued<F, T> && detail::derivative_valued<D, T>
            && detail::derivative_valued<D2, T> && std::default_initializable<M>
            && detail::from_bracket<M, T,
                detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>, std::remove_reference_t<D2>> const>
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

    /** @brief A method with values_of<M> == 3 from guess, the derivatives as separate callables. */
    template<typename M, real T, typename F, typename D, typename D2>
    requires(values_of<M> == 3) && detail::value_valued<F, T> && detail::derivative_valued<D, T>
            && detail::derivative_valued<D2, T> && std::default_initializable<M>
            && detail::from_guess<M, T,
                detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>, std::remove_reference_t<D2>> const>
    [[nodiscard]] constexpr std::expected<Root<T>, error>
        solve(F &&fn, D &&dfn, D2 &&d2fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
    {
        detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>, std::remove_reference_t<D2>> const
            joined{ std::addressof(fn), { std::addressof(dfn), std::addressof(d2fn) } };
        return nxx::roots::solve<M>(joined, guess, criteria);
    }

