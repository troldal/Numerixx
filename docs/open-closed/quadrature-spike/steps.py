# Draft 3 of the quadrature method protocol: draft 2 (quadrature.patch.py ... normalize), then the stepping protocol.
# Usage: steps.py <draft-2 header> <out>
import sys
src = open(sys.argv[1]).read()

def sub(old, new):
    global src
    assert src.count(old) == 1, (old[:100], src.count(old))
    src = src.replace(old, new)

def replace_between(start, end, new, keep_end=True):
    global src
    i = src.index(start); j = src.index(end, i)
    src = src[:i] + new + (src[j:] if keep_end else src[j + len(end):])

# 1. Step<T>, public, after Criteria (before the first detail block).
sub("    /** @cond INTERNAL */\n    namespace detail {\n\n        /**\n         * @brief The types whose Gauss-Legendre rules",
"""    /**
     * @brief Where an iterative method stands after a step: its current result, the least error it can certify, and
     * whether it may stop, can go on, and can succeed (docs/quadrature/PHASE3A_API.md, draft 3).
     * @tparam T The real type.
     */
    template<real T>
    struct Step
    {
        /** @brief The value and its estimate now, as integrate() would return them; empty where the method has no
         * usable estimate (TanhSinh: an unbounded end term, or no point strictly inside). */
        std::optional<Integral<T>> current;

        /** @brief The least error the method can certify: Romberg its rounding floor F, Gauss-Kronrod 2 F. */
        T attainable;

        /** @brief Whether the method has done the minimum work to stop (Romberg: level 4). */
        bool ready;

        /** @brief Whether the method cannot go on: its cap, stagnation, or every panel final. */
        bool exhausted;

        /** @brief Whether a result that meets the request may be returned as a success. */
        bool trusted;
    };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief The types whose Gauss-Legendre rules""")

# 2. Romberg as a state and a step per level.
replace_between("        /**\n         * @brief Romberg's method on an interval, with validated tolerances.",
"""            return std::unexpected(Failure<T>{ error::not_converged, best });
        }
""", keep_end=False, new="""        /**
         * @brief Romberg's method, stepped: the state after level k.
         * @tparam T The real type.
         */
        template<real T>
        struct RombergState
        {
            static constexpr int cap = 16;      // 2^16 + 1 evaluations; fixed, not from T (docs/quadrature/design.md)
            static constexpr int min_level = 4; // 17 evaluations: fewer samples of sin^2 on [0, 2 pi] all vanish
            static constexpr std::size_t row_size = cap + 1;

            Interval<T> interval;
            T weight;
            T l1;
            std::array<T, row_size> previous;
            std::array<T, row_size> row;
            int k;
            Step<T> step;
        };

        /**
         * @brief Level s.k of Romberg's method: T_k, the extrapolated row, the estimate and the floor; the body of
         * Phase 1's loop, unchanged, with the stop test moved to integrate().
         * @param fn The function.
         * @param s The state, with level k - 1 done (previous) and k set.
         * @return Nothing, or the Failure (non_finite).
         */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> romberg_level(F &fn, RombergState<T> &s)
        {
            Failure<T> const overflow{ error::non_finite, std::nullopt };
            int const k = s.k;
            Interval<T> const &interval = s.interval;
            if (k > 1) s.weight = T(s.weight / T(2));
            int const m = 1 << (k - 1);
            Sum<T> sum;
            T abs_sum = T(0); // a plain sum: its terms are positive, so it errs by about m eps relative
            for (int j = 1; j <= m; ++j) {
                // The new midpoints u_j = (2j - 1 - m) / m, exact (at most 16 bits).
                T const u = T(T(2 * j - 1 - m) / T(m));
                T const x = interval.point(T(interval.c + T(interval.hw * u)));
                std::optional<T> const t = interval.checked_term(s.weight, detail::evaluate(fn, x));
                if (!t) return std::unexpected(overflow);
                if consteval {
                    if (detail::large(sum.sum) || detail::large(abs_sum)) return std::unexpected(overflow);
                }
                sum.add(*t);
                abs_sum = T(abs_sum + nxx::math::abs(*t));
            }
            std::optional<T> const level = detail::compensated(sum);
            if (!level) return std::unexpected(overflow);
            if consteval {
                if (detail::large(*level) || detail::large(abs_sum)) return std::unexpected(overflow);
            }
            s.row[0] = T(T(s.previous[0] / T(2)) + *level);
            s.l1 = T(T(s.l1 / T(2)) + abs_sum);
            if (!nxx::math::isfinite(s.l1)) return std::unexpected(overflow);
            T power = T(1);
            for (int i = 1; i <= k; ++i) {
                std::size_t const index = static_cast<std::size_t>(i);
                if consteval {
                    if (detail::large(s.row[index - 1]) || detail::large(s.previous[index - 1]))
                        return std::unexpected(overflow);
                }
                power = T(power * T(4));
                s.row[index] = T(s.row[index - 1] + T(T(s.row[index - 1] - s.previous[index - 1]) / T(power - T(1))));
            }
            std::size_t const diagonal = static_cast<std::size_t>(k);
            if consteval {
                if (detail::large(s.row[diagonal])) return std::unexpected(overflow);
            }
            T const estimate = nxx::math::abs(T(s.row[diagonal] - s.previous[diagonal - 1]));
            if (!nxx::math::isfinite(estimate)) return std::unexpected(overflow);
            T const floor = detail::rounding_floor(s.l1);
            s.step = Step<T>{ Integral<T>{ s.row[diagonal], detail::larger(estimate, floor) }, floor,
                k >= RombergState<T>::min_level, k == RombergState<T>::cap, true };
            return {};
        }

        /**
         * @brief Romberg's start: T_0 from the ends, then level 1.
         * @return The state, or the Failure.
         */
        template<real T, typename F>
        constexpr std::expected<RombergState<T>, Failure<T>> romberg_start(F &fn, T const &a, T const &b)
        {
            Failure<T> const overflow{ error::non_finite, std::nullopt };
            RombergState<T> s{ detail::interval_of(a, b), T(1), T(0), detail::zeros<T, RombergState<T>::row_size>(),
                detail::zeros<T, RombergState<T>::row_size>(), 1,
                Step<T>{ std::nullopt, T(0), false, false, true } };
            std::optional<T> const ta = s.interval.checked_term(s.weight, detail::evaluate(fn, a));
            std::optional<T> const tb = s.interval.checked_term(s.weight, detail::evaluate(fn, b));
            if (!ta || !tb) return std::unexpected(overflow);
            s.l1 = T(nxx::math::abs(*ta) + nxx::math::abs(*tb));
            if (!nxx::math::isfinite(s.l1)) return std::unexpected(overflow);
            s.previous[0] = T(*ta + *tb);
            std::expected<void, Failure<T>> const r = detail::romberg_level(fn, s);
            if (!r) return std::unexpected(r.error());
            return s;
        }

        /**
         * @brief Romberg's next level.
         * @return Nothing, or the Failure.
         */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> romberg_advance(F &fn, RombergState<T> &s)
        {
            s.previous = s.row;
            ++s.k;
            return detail::romberg_level(fn, s);
        }
""")

# 3. Gauss-Kronrod as a state machine: the start (the mesh), a step (one bisection), and settle (the fresh sums).
replace_between("        /**\n         * @brief Adaptive Gauss-Kronrod on [a, b], a != b", "\n    } // namespace detail\n    /** @endcond */", """        /** @brief The running sums of the adaptive loop: the value K, the floor F and the excess X. */
        template<real T>
        struct KronrodTotals
        {
            Sum<T> value;
            Sum<T> floor;
            Sum<T> excess;
        };

        /** @brief The values of a set of sums; the estimate is E = F + X. */
        template<real T>
        struct KronrodValues
        {
            T value;
            T floor;
            T excess;
            T estimate;
        };

        /**
         * @brief Adaptive Gauss-Kronrod, stepped: Phase 2's loop (docs/quadrature/PHASE2_API.md C to E) with its locals
         * as members, a step per bisection, and the stop test moved to integrate(). The rule is copied: for a type
         * without a table, in constant evaluation, with_gauss_kronrod passes a temporary.
         * @tparam T The real type.
         */
        template<real T>
        struct KronrodState
        {
            static constexpr int mesh = 10;       // Shampine's quadgk: one panel cannot accept K = G on a peak it missed
            static constexpr int depth = 2;       // bisections below a panel whose L1 is not finite
            static constexpr int stagnation = 32; // bisections without a new minimum of the excess X

            Kronrod<T> rule;
            // Bisections and descent splits together, per bit of T: x^-0.9 on [0, 1] at cpp_bin_float_50's default
            // needs 3,550 of its 10,752.
            int cap;
            bool reversed; // integrated over [min, max] and negated, so that the panels' ends are in order
            // Every current panel, active and final, in the order it was committed; a bisected parent is marked
            // removed. Every decision to stop is checked against sums formed afresh from them (settle).
            std::vector<Panel<T>> panels;
            std::vector<bool> removed;
            KronrodTotals<T> running;
            bool refused; // a sum constant evaluation could not form within max / 2
            std::vector<std::size_t> heap;
            std::size_t order;
            int work; // bisections and descent splits
            std::vector<Panel<T>> pieces;
            std::vector<Pending<T>> stack;
            bool tainted; // a piece no descent avoided was undone: the run cannot end ok
            T least;      // the least excess so far
            int since;    // bisections since the excess last reached a new minimum
            Step<T> step;

            constexpr KronrodState(Kronrod<T> const &r, bool rev)
                : rule(r), cap(64 * detail::significand_bits<T>()), reversed(rev), panels(), removed(), running(),
                  refused(false), heap(), order(0), work(0), pieces(), stack(), tainted(false), least(T(0)), since(0),
                  step{ std::nullopt, T(0), false, false, true }
            {}

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

            // The three sums afresh, from the current panels.
            constexpr KronrodTotals<T> fresh()
            {
                KronrodTotals<T> t;
                for (std::size_t i = 0; i < panels.size(); ++i)
                    if (!removed[i]) contribute(t, panels[i], true);
                return t;
            }

            // The heap's order: by difference, the same order as the estimate.
            constexpr bool precedes(std::size_t p, std::size_t q) const
            {
                Panel<T> const &u = panels[p];
                Panel<T> const &v = panels[q];
                return u.difference < v.difference || (u.difference == v.difference && v.order < u.order);
            }

            // [p, q] sampled into pieces, or a descent of at most depth bisections below it.
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
                        stack.push_back(Pending<T>{ mid, top.hi, top.level + 1 });
                        stack.push_back(Pending<T>{ top.lo, mid, top.level + 1 });
                    }
                    else missing = true;
                }
                if (missing) return fault::not_finite;
                return std::nullopt;
            }

            // The pieces enter the panels and the running sums, and the heap unless final.
            constexpr void commit()
            {
                auto const before = [this](std::size_t p, std::size_t q) { return precedes(p, q); };
                for (Panel<T> const &piece : pieces) {
                    panels.push_back(piece);
                    removed.push_back(false);
                    contribute(running, piece, true);
                    if (!(piece.difference <= piece.floor) && detail::splittable(rule, piece.lo, piece.hi)) {
                        heap.push_back(panels.size() - 1);
                        std::push_heap(heap.begin(), heap.end(), before);
                    }
                }
            }

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

            // The step from a set of values: the result {K, max(E, F)}, and twice the floor as the attainable error
            // (a kink or a jump never brings every panel to its floor).
            constexpr void report(KronrodValues<T> const &v)
            {
                step = Step<T>{ Integral<T>{ reversed ? T(-v.value) : v.value, detail::larger(v.estimate, v.floor) },
                    T(T(2) * v.floor), true, heap.empty() || work >= cap || since >= stagnation, !tainted };
            }
        };

        /**
         * @brief Gauss-Kronrod's start: the initial mesh of 10 equal panels (empty ones skipped); a mesh panel with no
         * finite estimate is non_finite.
         */
        template<real T, typename F>
        constexpr std::expected<KronrodState<T>, Failure<T>>
            kronrod_start(Kronrod<T> const &rule, F &fn, T const &a, T const &b)
        {
            Failure<T> const invalid{ error::non_finite, std::nullopt };
            bool const reversed = b < a;
            T const lower = reversed ? b : a;
            T const upper = reversed ? a : b;
            KronrodState<T> s(rule, reversed);
            Interval<T> const whole = detail::interval_of(lower, upper);
            T left = lower;
            for (int i = 1; i <= KronrodState<T>::mesh; ++i) {
                T const u = T(T(-1) + T(T(2 * i) / T(KronrodState<T>::mesh)));
                T const right = i == KronrodState<T>::mesh ? upper : whole.point(T(whole.c + T(whole.hw * u)));
                if (left < right) {
                    s.pieces.clear();
                    if (s.resolve(fn, left, right)) return std::unexpected(invalid);
                    s.commit();
                    if (s.refused) return std::unexpected(invalid);
                }
                left = right;
            }
            std::optional<KronrodValues<T>> const after_mesh = KronrodState<T>::values_of(s.running);
            if (!after_mesh) return std::unexpected(invalid);
            s.least = after_mesh->excess;
            s.since = 0;
            s.report(*after_mesh);
            return s;
        }

        /**
         * @brief Gauss-Kronrod's step: the panel with the largest estimate bisected (precondition: not exhausted).
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
                // Undone: the parent stays, final, with its own estimate; the descent's pieces are dropped.
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
            if (current->excess < s.least) {
                s.least = current->excess;
                s.since = 0;
            }
            else ++s.since;
            s.report(*current);
            return {};
        }

        /**
         * @brief Gauss-Kronrod's settle: the step from sums formed afresh, which every decision to stop uses.
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
""")

# 4. The methods' members: start, advance (and settle) instead of estimate.
sub("""        template<real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            estimate(F &fn, T const &a, T const &b, Criteria<T> const &criteria) const
        { return detail::romberg(fn, a, b, criteria.abs_tol, criteria.rel_tol); }""",
"""        template<real T, typename F>
        constexpr std::expected<detail::RombergState<T>, Failure<T>> start(F &fn, T const &a, T const &b) const
        { return detail::romberg_start(fn, a, b); }

        /** @brief The next level (iterative_method). */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> advance(F &fn, detail::RombergState<T> &s) const
        { return detail::romberg_advance(fn, s); }""")
sub("""        template<real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            estimate(F &fn, T const &a, T const &b, Criteria<T> const &criteria) const
        {
            return detail::with_gauss_kronrod<T>([&](detail::Kronrod<T> const &rule) {
                return detail::gauss_kronrod(rule, fn, a, b, criteria.abs_tol, criteria.rel_tol);
            });
        }""",
"""        template<real T, typename F>
        constexpr std::expected<detail::KronrodState<T>, Failure<T>> start(F &fn, T const &a, T const &b) const
        {
            return detail::with_gauss_kronrod<T>(
                [&](detail::Kronrod<T> const &rule) { return detail::kronrod_start(rule, fn, a, b); });
        }

        /** @brief One bisection (iterative_method). */
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> advance(F &fn, detail::KronrodState<T> &s) const
        { return detail::kronrod_advance(fn, s); }

        /** @brief The current step from sums formed afresh, before any decision to stop. */
        template<real T>
        constexpr std::expected<void, Failure<T>> settle(detail::KronrodState<T> &s) const
        { return detail::kronrod_settle(s); }""")

# 5. The iterative protocol, the driver, Progress and observed(), after estimating_method.
sub("""    /** @cond INTERNAL */
    namespace detail {

        /** @brief Whether an Integral is finite with a non-negative estimate (the protocol's checked part). */""",
"""    /** @cond INTERNAL */
    namespace detail {

        /** @brief R is std::expected<S, Failure<T>> for some state S. */
        template<typename R, typename T>
        concept expected_state = requires { typename R::value_type; }
                                 && std::same_as<R, std::expected<typename R::value_type, Failure<T>>>;

        /** @brief The state type of an iterative method. */
        template<typename M, typename T, typename F>
        using state_of = typename decltype(std::declval<M const &>().start(
            std::declval<F &>(), std::declval<T const &>(), std::declval<T const &>()))::value_type;

    } // namespace detail
    /** @endcond */

    /**
     * @brief An iterative method: start(fn, a, b) returns std::expected<S, Failure<T>>, a state S whose member step
     * is the current Step<T>; advance(fn, s) takes one more step; an optional settle(s) refines the step before any
     * decision to stop. integrate() owns the loop and the test against the criteria.
     * @tparam M The method. @tparam T The real type. @tparam F fn's type, reference removed.
     */
    template<typename M, typename T, typename F = T(T const &)>
    concept iterative_method = requires(M const &m, F &fn, T const &a, T const &b) {
        { m.start(fn, a, b) } -> detail::expected_state<T>;
    } && requires(M const &m, F &fn, detail::state_of<M, T, F> &s) {
        { m.advance(fn, s) } -> std::same_as<std::expected<void, Failure<T>>>;
        { std::as_const(s).step } -> std::same_as<Step<T> const &>;
    };

    /**
     * @brief What an observer of an iterative method sees after its start (iteration 0) and after each step.
     * @tparam T The real type.
     */
    template<real T>
    struct Progress
    {
        /** @brief 0 after the start, then one per step. */
        int iteration;

        /** @brief The evaluations of fn so far. */
        std::size_t evaluations;

        /** @brief The method's step. */
        Step<T> step;
    };

    /** @cond INTERNAL */
    namespace detail {

        /** @brief fn, counting its evaluations. */
        template<typename F>
        struct Counted
        {
            F *fn;
            std::size_t *count;

            template<typename X>
            constexpr decltype(auto) operator()(X const &x) const
            {
                ++*count;
                return std::invoke(*fn, x);
            }
        };

        /** @brief The state of an observed method: the inner state, the count, the iteration, the step. */
        template<typename S, real T>
        struct ObservedState
        {
            S inner;
            std::size_t evaluations;
            int iteration;
            Step<T> step;
        };

    } // namespace detail
    /** @endcond */

    /**
     * @brief An iterative method M whose progress is reported: observer(Progress<T>) after the start and after each
     * step. Itself an iterative method; made by observed().
     */
    template<typename M, typename Observer>
    struct Observed
    {
        M method;
        Observer observer;

        /** @brief M's start, counted and reported. */
        template<real T, typename F>
        requires iterative_method<M, T, detail::Counted<F>>
        constexpr std::expected<detail::ObservedState<detail::state_of<M, T, detail::Counted<F>>, T>, Failure<T>>
            start(F &fn, T const &a, T const &b) const
        {
            std::size_t count = 0;
            detail::Counted<F> counted{ std::addressof(fn), &count };
            auto started = method.start(counted, a, b);
            if (!started) return std::unexpected(started.error());
            Step<T> const first = started->step;
            detail::ObservedState<detail::state_of<M, T, detail::Counted<F>>, T> s{ std::move(*started), count, 0,
                first };
            std::invoke(observer, Progress<T>{ 0, s.evaluations, s.step });
            return s;
        }

        /** @brief M's step, counted and reported. */
        template<real T, typename F, typename S>
        constexpr std::expected<void, Failure<T>> advance(F &fn, detail::ObservedState<S, T> &s) const
        {
            detail::Counted<F> counted{ std::addressof(fn), &s.evaluations };
            std::expected<void, Failure<T>> const r = method.advance(counted, s.inner);
            if (!r) return r;
            s.step = s.inner.step;
            ++s.iteration;
            std::invoke(observer, Progress<T>{ s.iteration, s.evaluations, s.step });
            return {};
        }

        /** @brief M's settle, where it has one. */
        template<real T, typename S>
        constexpr std::expected<void, Failure<T>> settle(detail::ObservedState<S, T> &s) const
        {
            if constexpr (requires { method.settle(s.inner); }) {
                std::expected<void, Failure<T>> const r = method.settle(s.inner);
                if (!r) return r;
                s.step = s.inner.step;
            }
            return {};
        }
    };

    /**
     * @brief An iterative method that reports its progress to observer: integrate(observed(m, observer), fn, a, b).
     * @param m The method.
     * @param observer A callable taking Progress<T> const &.
     * @return The observed method.
     */
    template<typename M, typename Observer>
    constexpr Observed<M, std::decay_t<Observer>> observed(M const &m, Observer &&observer)
    { return Observed<M, std::decay_t<Observer>>{ m, std::forward<Observer>(observer) }; }

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief The loop of an iterative method on a < b: test the step against the criteria; on a pass, or where the
         * method is exhausted, settle and decide on the settled step; else take a step.
         */
        template<typename M, real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            drive(M const &m, F &fn, T const &a, T const &b, Criteria<T> const &criteria)
        {
            auto started = m.start(fn, a, b);
            if (!started) return std::unexpected(started.error());
            auto &s = *started;
            // The request: abs_error <= max(abs_tol, rel_tol |value|, attainable), once the method is ready.
            auto const met = [&criteria](Step<T> const &p) {
                return p.ready && p.current
                       && p.current->abs_error
                              <= detail::wanted(criteria.abs_tol, criteria.rel_tol, p.current->value, p.attainable);
            };
            while (true) {
                bool const passed = met(s.step);
                if (passed || s.step.exhausted) {
                    if constexpr (requires { m.settle(s); }) {
                        std::expected<void, Failure<T>> const r = m.settle(s);
                        if (!r) return std::unexpected(r.error());
                    }
                    if (met(s.step)) {
                        if (!s.step.trusted) return std::unexpected(Failure<T>{ error::not_converged, s.step.current });
                        return *s.step.current;
                    }
                    if (!passed) return std::unexpected(Failure<T>{ error::not_converged, s.step.current });
                    continue; // the settled step no longer passes: go on from it
                }
                std::expected<void, Failure<T>> const r = m.advance(fn, s);
                if (!r) return std::unexpected(r.error());
            }
        }

        /** @brief Whether an Integral is finite with a non-negative estimate (the protocol's checked part). */""")

# 6. integrate takes either kind of estimating method; an iterative one is driven.
sub("""    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T> && estimating_method<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(M const &m, F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {})""",
"""    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T>
             && (iterative_method<M, T, std::remove_reference_t<F>> || estimating_method<M, T, std::remove_reference_t<F>>)
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(M const &m, F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {})""")
sub("""        std::expected<Integral<T>, Failure<T>> const r = b < a
            ? m.estimate(fn, b, a, validated).transform(negated).transform_error([&](Failure<T> f) {
                  if (f.best) f.best = negated(*f.best);
                  return f;
              })
            : m.estimate(fn, a, b, validated);""",
"""        auto const run = [&](T const &lo, T const &hi) -> std::expected<Integral<T>, Failure<T>> {
            if constexpr (iterative_method<M, T, std::remove_reference_t<F>>) return detail::drive(m, fn, lo, hi, validated);
            else return m.estimate(fn, lo, hi, validated);
        };
        std::expected<Integral<T>, Failure<T>> const r = b < a
            ? run(b, a).transform(negated).transform_error([&](Failure<T> f) {
                  if (f.best) f.best = negated(*f.best);
                  return f;
              })
            : run(a, b);""")
sub("""    template<typename M = GaussKronrod<15>, real T, typename F>
    requires std::default_initializable<M> && nxx::detail::callable<F, T>
             && estimating_method<M, T, std::remove_reference_t<F>>""",
"""    template<typename M = GaussKronrod<15>, real T, typename F>
    requires std::default_initializable<M> && nxx::detail::callable<F, T>
             && (iterative_method<M, T, std::remove_reference_t<F>> || estimating_method<M, T, std::remove_reference_t<F>>)""")
sub("#include <type_traits>\n", "#include <type_traits>\n#include <utility>\n") if src.count("#include <utility>\n") == 0 else None
open(sys.argv[2], "w").write(src)
