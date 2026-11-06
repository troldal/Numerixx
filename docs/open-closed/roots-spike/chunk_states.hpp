
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
         * of the interpolation overflows. A derivative step is accepted by Brent's test written for a step s: towards c,
         * |s| < 3 |m| / 2 - tol / 2, and |s| < |e| / 2; one that cannot be formed is a bisection step.
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
                if (T(T(2) * p) < T(T(T(3) * m * q) - math::abs(T(tol * q))) && p < math::abs(T(before_last * q / T(2)))) {
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
            // arithmetic where the stop test failed (PHASE7_API.md B); not_converged rather than a root outside the bound.
            if (!((s.b < next && next < s.c) || (s.c < next && next < s.b))) return std::unexpected(error::not_converged);
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
            // Not reached: a midpoint not strictly inside means adjacent ends, which met the stop test before (midpoint()
            // falls back on (a + b) / 2 where a + (b - a) / 2 flushes, numerics review of draft 1, finding 1).
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
         * @return Nothing, or error::non_finite.
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
         * @return Nothing, or error::non_finite.
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
         * @return Nothing, or error::non_finite.
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
         * @return Nothing, or error::non_finite.
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
         * @return Nothing, or error::non_finite.
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
         * @return Nothing, or error::non_finite.
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
