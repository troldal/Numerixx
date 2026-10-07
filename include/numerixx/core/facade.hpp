// Solver configuration and the family facades (DESIGN §6.6).
//
// Every solver holds one options aggregate: stop criterion, iteration budget, derivative source, projection and
// observer. The builders with_stop, with_budget, with_derivative, with_projection and with_observer are generic: they
// rebind that aggregate and return a new solver value, so configuration is order-independent and the old value is
// unchanged. User callables are held in copyable boxes, so a solver that holds a lambda stays copy-assignable.
//
// The family facades (bracketing, open, search) are deducing-this bases, no CRTP. Each constrains operator() and .on()
// on the inputs its solvers accept and deletes everything else with a reason, so std::is_invocable_v is false (not a
// hard error) for a solver given the wrong kind of input. operator() is also constrained on the solver protocol
// (detail::runnable_v), so a solver that lacks a protocol member gives false too, and a reason. A member of the wrong
// type is still a hard error inside detail::run or nxx::iterate (DESIGN §6.6).
#pragma once

#include <numerixx/config.hpp>
#include <numerixx/core/callable.hpp>
#include <numerixx/core/criteria.hpp>
#include <numerixx/core/iterate.hpp>
#include <numerixx/core/refined.hpp>
#include <numerixx/core/scalar.hpp>

#include <cstddef>
#include <expected>
#include <functional>
#include <type_traits>
#include <utility>

NXX_BEGIN_HEADER

namespace nxx
{
    // No derivative source configured (Newton then needs a callable with .derivative()).
    struct no_derivative
    {
        friend constexpr bool operator==(no_derivative, no_derivative) = default;
    };

    // The identity projection.
    struct no_projection
    {
        template<class X>
        constexpr X operator()(const X& x) const noexcept
        { return x; }
    };

    // An observer that ignores every iterate.
    struct no_observer
    {
        template<class V>
        constexpr void operator()(const V&) const noexcept
        {}
    };

    namespace detail
    {
        // Key for the constructor that takes a finished options aggregate (used by the builders).
        struct from_options_t
        {
            explicit constexpr from_options_t() = default;
        };
        inline constexpr from_options_t from_options {};

        // Whether stop criterion C would try to tighten solver S's own tolerance: S has one (internal_tolerance) and C
        // contains a width criterion. Such a solver's intrinsic test runs before the stop criterion and reports success
        // once its own tolerance holds, so with_stop can only add an early exit or a failure (DESIGN §6.8).
        template<class S, class C>
        inline constexpr bool tightens_tolerance_v = [] {
            if constexpr (S::internal_tolerance)
                return contains_width_v<C>;
            else
                return false;
        }();

        // Whether solver S's with_stop and public rebuild accept stop criterion C: C can stop S alone (a bare
        // min_iterations guard cannot) and does not try to tighten S's own tolerance. with_stop and every solver's
        // rebuild use this one predicate, so a solver that sets internal_tolerance gets both halves of the rule; such a
        // solver also constrains its from_options constructor with it (DESIGN §6.8). False, not a hard error, for a type
        // without views, and on GCC, Clang and clang-cl also for a views that is not a view_kind constant (an int, a data
        // member, a function): forming S::views as the template argument is part of the nested requirement's
        // substitution. cl 19.51 treats that invalid argument as a hard error; it rejects such a malformed solver at
        // with_stop's constraint anyway. tightens_tolerance_v is asked only after it, because an && in a variable
        // template's initializer does not stop the instantiation of its later operands.
        template<class S, class C>
        inline constexpr bool stop_allowed_v = [] {
            if constexpr (!requires { requires stop_criterion_for_v<C, S::views>; })
                return false;
            else
                return !tightens_tolerance_v<S, C>;
        }();
    }    // namespace detail

    template<class Stop, class Deriv = no_derivative, class Proj = no_projection, class Obs = no_observer>
    struct options
    {
        using stop_type       = Stop;
        using derivative_type = Deriv;
        using projection_type = Proj;
        using observer_type   = Obs;

        Stop                  stop;
        max_iterations        budget;
        NXX_NO_UNIQUE_ADDRESS detail::copyable_box<Deriv> derivative {};
        NXX_NO_UNIQUE_ADDRESS detail::copyable_box<Proj> projection {};
        NXX_NO_UNIQUE_ADDRESS detail::copyable_box<Obs> observer {};

        template<class V>
        constexpr void observe(const V& view) const
        { std::invoke(*observer, view); }

        template<class X>
        constexpr X project(const X& x) const
        { return static_cast<X>(std::invoke(*projection, x)); }
    };

    // The builders shared by every solver. A solver provides options(), rebuild(options) (constrained on
    // detail::stop_allowed_v<Self, typename O2::stop_type>, as with_stop is) and the kind of its views,
    // static constexpr view_kind views (named so because view(s) is the protocol function).
    struct solver_facade
    {
        // Whether the solver is ready to be called with a function of type F (Newton: it has a derivative source).
        template<class F>
        static constexpr bool ready_v = true;

        // Whether F can be called with the scalar type of an accepted input In; solvers specialise it.
        template<class F, class In>
        static constexpr bool callable_v = true;

        // Whether the solver uses a derivative source, and whether it projects its iterates. A builder that would be
        // ignored is deleted with a reason instead (DESIGN §3.3: solver values never hold a configuration that means
        // nothing).
        static constexpr bool uses_derivative = false;
        static constexpr bool projects        = false;

        // Whether the solver has its own tolerance, which its intrinsic test applies (brent; later golden, brent_min,
        // toms748, itp). Its constructor takes the width criterion; with_stop rejects one, at any depth of || and &&,
        // because the intrinsic test reports success before the stop criterion is asked (DESIGN §6.8).
        static constexpr bool internal_tolerance = false;

        template<class Self>
        constexpr auto with_budget(this const Self& self, max_iterations budget)
        {
            auto o   = self.options();
            o.budget = budget;
            return self.rebuild(std::move(o));
        }

        template<class Self, class C>
            requires detail::stop_allowed_v<Self, C>
        constexpr auto with_stop(this const Self& self, C stop)
        {
            const auto& o = self.options();
            using O       = std::remove_cvref_t<decltype(o)>;
            using O2      = options<C, typename O::derivative_type, typename O::projection_type, typename O::observer_type>;
            return self.rebuild(O2 { std::move(stop), o.budget, o.derivative, o.projection, o.observer });
        }

        // Not a criterion at all: a number, a validated tolerance (tolerance<T>) or a part (abs_tolerance<T>,
        // rel_tolerance<T>). The text names the tests in x first and says what each one bounds (DESIGN §6.6, §9.3). It
        // and the criterion catch-all below are disjoint, so cl reports no ambiguity (C2668), and this header needs no
        // tolerance trait.
        template<class Self, class C>
            requires(!is_criterion_v<C>)
        void with_stop(this const Self&, C) NXX_DELETE("with_stop takes a criterion, not a number or a validated tolerance: wrap "
                                                       "it in the test you mean: width_tol{*tol} (bisection; brent takes its "
                                                       "width in its constructor) bounds the error in x; x_tol{*tol} (open "
                                                       "methods) bounds only the last step in x; f_tol{*tol} bounds only "
                                                       "|f(x)|; a part (abs_tolerance, rel_tolerance) is not a criterion either");

        template<class Self, class C>
            requires(is_criterion_v<C> && !criterion_for_v<C, Self::views>)
        void with_stop(this const Self&, C) NXX_DELETE("this criterion does not apply to this solver: bracketing methods "
                                                       "converge on the enclosure (width_tol, floored_width), open methods on "
                                                       "successive iterates (x_tol, step_tol)");

        template<class Self, class C>
            requires(criterion_for_v<C, Self::views> && detail::guard_only_v<C>)
        void with_stop(this const Self&, C) NXX_DELETE("min_iterations only guards another criterion: combine it with a "
                                                       "convergence test using && (your_test && min_iterations{n})");

        template<class Self, class C>
            requires(stop_criterion_for_v<C, Self::views> && detail::tightens_tolerance_v<Self, C>)
        void with_stop(this const Self&, C) NXX_DELETE("this solver has its own tolerance: pass the width criterion to its "
                                                       "constructor (brent{nxx::width_tol{1e-10}}); with_stop adds an early "
                                                       "exit or a failure (f_tol, max_evaluations) and cannot tighten that "
                                                       "tolerance");

        template<class Self, class D>
            requires Self::uses_derivative
        constexpr auto with_derivative(this const Self& self, D derivative)
        {
            const auto& o = self.options();
            using O       = std::remove_cvref_t<decltype(o)>;
            using O2      = options<typename O::stop_type, D, typename O::projection_type, typename O::observer_type>;
            return self.rebuild(O2 { o.stop, o.budget, detail::copyable_box<D> { std::move(derivative) }, o.projection, o.observer });
        }

        template<class Self, class D>
            requires(!Self::uses_derivative)
        void with_derivative(this const Self&, D) NXX_DELETE("this solver does not use a derivative (newton does)");

        // Applied to every proposed iterate before the function is evaluated (D20).
        template<class Self, class P>
            requires Self::projects
        constexpr auto with_projection(this const Self& self, P projection)
        {
            const auto& o = self.options();
            using O       = std::remove_cvref_t<decltype(o)>;
            using O2      = options<typename O::stop_type, typename O::derivative_type, P, typename O::observer_type>;
            return self.rebuild(O2 { o.stop, o.budget, o.derivative, detail::copyable_box<P> { std::move(projection) }, o.observer });
        }

        template<class Self, class P>
            requires(!Self::projects)
        void with_projection(this const Self&, P) NXX_DELETE("projection applies to open methods (secant, newton): a "
                                                             "bracketing method keeps every iterate inside its "
                                                             "bracket");

        // Called with the view of every iterate: logging without polluting the stop criteria.
        template<class Self, class Ob>
        constexpr auto with_observer(this const Self& self, Ob observer)
        {
            const auto& o = self.options();
            using O       = std::remove_cvref_t<decltype(o)>;
            using O2      = options<typename O::stop_type, typename O::derivative_type, typename O::projection_type, Ob>;
            return self.rebuild(O2 { o.stop, o.budget, o.derivative, o.projection, detail::copyable_box<Ob> { std::move(observer) } });
        }
    };

    // solver.on(input): a curried solver, a value f -> result. Copy-assignable whenever the solver and input are.
    template<class S, class In>
    class bound
    {
        S  solver_;
        In in_;

    public:
        constexpr bound(S solver, In in) : solver_(std::move(solver)), in_(std::move(in)) {}

        // Constrained, so std::is_invocable_v on a curried solver is false (not a hard error) for a function it cannot
        // take. The deleted sibling keeps it false and gives a reason; the solver's own call states the cause (DESIGN §6.6).
        template<class F>
            requires std::is_invocable_v<const S&, const F&, const In&>
        constexpr auto operator()(const F& fn) const
        { return solver_(fn, in_); }

        template<class F>
            requires(!std::is_invocable_v<const S&, const F&, const In&>)
        void operator()(const F&) const NXX_DELETE("this solver cannot take this function with its bound input: call "
                                                   "solver(f, input) for the reason");

        constexpr const S&  solver() const noexcept { return solver_; }
        constexpr const In& input() const noexcept { return in_; }
    };

    namespace detail
    {
        // prepare, then iterate: the one-shot path of every facade. The function is held by reference for the call.
        template<class S, class F, class In>
        constexpr auto run(const S& solver, const F& fn, const In& in)
        {
            auto p  = solver.prepare(std::cref(fn), in);
            using R = decltype(nxx::iterate(solver, *p));
            if (!p) return R { std::unexpect, std::move(p).error() };
            return nxx::iterate(solver, *p);
        }

        // Whether S states the inputs it takes: a static bool variable template accepts_v<In> (DESIGN §6.6). An accepts_v
        // that is not a template is ruled out first, because GCC 16 makes `S::template accepts_v<In>` a hard error for
        // it, even inside a requires-expression.
        template<class S, class In>
        inline constexpr bool states_inputs_v = [] {
            if constexpr (requires { S::accepts_v; })
                return false;
            else
                return requires { std::bool_constant<S::template accepts_v<std::remove_cvref_t<In>>> {}; };
        }();
        // False, not a hard error, for a solver without accepts_v. Each test is asked only once the one before it holds:
        // an && in a variable template's initializer does not stop the instantiation of its later operands.
        template<class S, class In>
        inline constexpr bool accepts_v = [] {
            if constexpr (states_inputs_v<S, In>)
                return bool(S::template accepts_v<std::remove_cvref_t<In>>);
            else
                return false;
        }();
        template<class S, class F>
        inline constexpr bool ready_v = S::template ready_v<F>;
        // Whether F can be called with the scalar type of the input (only asked once accepts_v holds). F decays: for a
        // plain function F is a function type, and the solver's `const F&` would then qualify a function type, which MSVC
        // warns about (C4180); a function and a pointer to it are invocable alike.
        template<class S, class F, class In>
        inline constexpr bool input_callable_v = S::template callable_v<std::decay_t<F>, std::remove_cvref_t<In>>;
        // An input that is not accepted, and not a C array (MSVC cannot order the array overload against this one). A
        // solver that does not state its inputs at all is incomplete instead.
        template<class S, class In>
        inline constexpr bool rejected_v = states_inputs_v<S, In> && !accepts_v<S, In> && !std::is_array_v<std::remove_cvref_t<In>>;

        // What prepare returns on run's path: the function by std::cref (F& rather than const F&, which would qualify a
        // function type), then the input.
        template<class S, class F, class In>
        using prepared_t = decltype(std::declval<const S&>().prepare(std::cref(std::declval<F&>()), std::declval<const In&>()));

        // Whether run can call S (DESIGN §6.6): prepare(std::cref(f), in) returns a std::expected problem, and S is an
        // iterative_solver_for that problem (id, options, init, step, view, estimate, best, intrinsic, nfev, and
        // better_than for the estimate_type of init's failure). Without it, a solver that lacks a member made
        // std::is_invocable_v a hard error inside run. Members are checked for presence, not for every type the driver
        // needs: a prepare error that init's failure cannot hold, an options() that is not an options aggregate or an
        // init error that has an estimate_type but is not a failure is still a hard error in run or iterate (DESIGN
        // §6.6). Asked only once accepts_v, ready_v and input_callable_v hold.
        template<class S, class F, class In>
        inline constexpr bool runnable_v = [] {
            if constexpr (requires {
                              typename prepared_t<S, F, In>::value_type;
                              typename prepared_t<S, F, In>::error_type;
                          })
                return iterative_solver_for<S, typename prepared_t<S, F, In>::value_type>;
            else
                return false;
        }();

        // A solver that cannot take part in the call: it does not state its inputs, or it accepts this one and cannot
        // run it. Arrays are left to the array overloads, which ask runnable_v themselves.
        template<class S, class F, class In>
        inline constexpr bool incomplete_v = [] {
            if constexpr (std::is_array_v<std::remove_cvref_t<In>>)
                return false;
            else if constexpr (!states_inputs_v<S, In>)
                return true;
            else if constexpr (!accepts_v<S, In>)
                return false;
            else if constexpr (!ready_v<S, F>)
                return false;
            else if constexpr (!input_callable_v<S, F, In>)
                return false;
            else
                return !runnable_v<S, F, In>;
        }();
    }    // namespace detail

    // Bracketing methods: a bracket<T>, a braced {lo, hi}, a std::pair, the result of bracket<T>::make, a sign_bracket
    // or a search result.
    struct bracketing_facade : solver_facade
    {
        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && detail::input_callable_v<Self, F, In> &&
                     detail::runnable_v<Self, F, In>)
        constexpr auto operator()(this const Self& self, const F& fn, const In& in)
        { return detail::run(self, fn, in); }

        // A braced list or a C array: N is deduced, so {x} is not taken as {x, 0} and {a, b, c} is not cut short.
        template<class Self, class F, real T, std::size_t N>
            requires(N == 2 && detail::ready_v<Self, F> && detail::input_callable_v<Self, F, std::pair<T, T>> &&
                     detail::runnable_v<Self, F, std::pair<T, T>>)
        constexpr auto operator()(this const Self& self, const F& fn, const T (&lo_hi)[N])
        { return detail::run(self, fn, std::pair<T, T> { lo_hi[0], lo_hi[1] }); }

        // A solver that does not implement the protocol (DESIGN §6.6) for {lo, hi}, and below for any other input.
        template<class Self, class F, real T, std::size_t N>
            requires(N == 2 && detail::ready_v<Self, F> && detail::input_callable_v<Self, F, std::pair<T, T>> &&
                     !detail::runnable_v<Self, F, std::pair<T, T>>)
        void operator()(this const Self&, const F&, const T (&)[N]) NXX_DELETE("this solver does not implement the solver protocol "
                                                                               "(DESIGN 6.6): it needs accepts_v, prepare(f, in), "
                                                                               "and id, options(), init, step, view, estimate, best "
                                                                               "and intrinsic for the problem prepare returns, and "
                                                                               "better_than(const Est&, const Est&) for its estimate type");

        template<class Self, class F, real T, std::size_t N>
            requires(N == 2 && detail::ready_v<Self, F> && !detail::input_callable_v<Self, F, std::pair<T, T>>)
        void operator()(this const Self&, const F&, const T (&)[N]) NXX_DELETE("the function cannot be called with the "
                                                                               "scalar type of the bracket");

        template<class Self, class F, class T, std::size_t N>
            requires(N != 2 || !real<T>)
        void operator()(this const Self&, const F&, const T (&)[N]) NXX_DELETE("a bracket has two ends of a real type: "
                                                                               "write {lo, hi}, for example {1.0, 2.0}");

        template<class Self, class F, class In>
            requires detail::rejected_v<Self, In>
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("bracketing solvers need a bracket: pass "
                                                                          "{lo, hi}, nxx::bracket<T>::make(a, b), or a "
                                                                          "search result");

        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && !detail::input_callable_v<Self, F, In>)
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("the function cannot be called with the "
                                                                          "scalar type of the bracket");

        template<class Self, class F, class In>
            requires detail::incomplete_v<Self, F, In>
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("this solver does not implement the solver protocol "
                                                                          "(DESIGN 6.6): it needs accepts_v, prepare(f, in), "
                                                                          "and id, options(), init, step, view, estimate, best "
                                                                          "and intrinsic for the problem prepare returns, and "
                                                                          "better_than(const Est&, const Est&) for its estimate type");

        template<class Self, class In>
            requires detail::accepts_v<Self, In>
        constexpr auto on(this const Self& self, In in)
        { return bound<Self, In> { self, std::move(in) }; }

        template<class Self, real T, std::size_t N>
            requires(N == 2)
        constexpr auto on(this const Self& self, const T (&lo_hi)[N])
        { return bound<Self, std::pair<T, T>> { self, std::pair<T, T> { lo_hi[0], lo_hi[1] } }; }

        template<class Self, class T, std::size_t N>
            requires(N != 2 || !real<T>)
        void on(this const Self&, const T (&)[N]) NXX_DELETE("a bracket has two ends of a real type: write {lo, hi}, "
                                                             "for example {1.0, 2.0}");

        // A forwarding reference: by value, a C array would decay to a pointer, which cl cannot order against the array
        // overloads; rejected_v leaves arrays to them, so a pointer still gets the reason.
        template<class Self, class In>
            requires detail::rejected_v<Self, std::remove_cvref_t<In>>
        void on(this const Self&, In&&) NXX_DELETE("bracketing solvers need a bracket: pass {lo, hi}, "
                                                   "nxx::bracket<T>::make(a, b), or a search result");
    };

    // Open methods: a guess of a real type, or a root estimate (seeded with its x and f(x): no re-evaluation).
    struct open_facade : solver_facade
    {
        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && detail::input_callable_v<Self, F, In> &&
                     detail::runnable_v<Self, F, In>)
        constexpr auto operator()(this const Self& self, const F& fn, const In& in)
        { return detail::run(self, fn, in); }

        template<class Self, class F, class In>
            requires(detail::rejected_v<Self, In> && !std::is_integral_v<std::remove_cvref_t<In>>)
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("open methods take a guess of a real type or "
                                                                          "a root estimate; bracketing solvers take "
                                                                          "{lo, hi}");

        // A solver that does not state its inputs gets the protocol's reason below instead (not an ambiguous call).
        template<class Self, class F, class I>
            requires(std::is_integral_v<std::remove_cvref_t<I>> && detail::states_inputs_v<Self, I>)
        void operator()(this const Self&, const F&, const I&) NXX_DELETE("open methods take a guess of a real type: write 1.0, not 1");

        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && !detail::input_callable_v<Self, F, In>)
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("the function cannot be called with the type of the guess");

        // Newton's reason, kept here so that no solver declares an operator() of its own. A solver that did would have
        // to repeat `using open_facade::operator();`. Its operator()'s template-head would differ from each facade
        // overload's in the requires-clause, so the two would not correspond ([basic.scope.scope]/4, [temp.over.link]/6)
        // and the using-declaration would hide nothing ([namespace.udecl]/11). CLion's ReSharper C++ engine (2026.2)
        // treats the solver's overload as hiding the facade's, however, and marks every valid call as an error.
        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && !detail::ready_v<Self, F>)
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("newton needs a derivative: "
                                                                          ".with_derivative(df), "
                                                                          ".with_derivative(deriv::numeric{}), a "
                                                                          "callable with .derivative(), or use secant");

        template<class Self, class F, class In>
            requires detail::incomplete_v<Self, F, In>
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("this solver does not implement the solver protocol "
                                                                          "(DESIGN 6.6): it needs accepts_v, prepare(f, in), "
                                                                          "and id, options(), init, step, view, estimate, best "
                                                                          "and intrinsic for the problem prepare returns, and "
                                                                          "better_than(const Est&, const Est&) for its estimate type");

        // A braced list or a C array: open methods start from one point, so {lo, hi} (or {x}) gets a reason rather than a
        // bare "no matching function". Arrays are left out of the catch-alls above, so only this overload takes them.
        template<class Self, class F, class T, std::size_t N>
        void operator()(this const Self&, const F&, const T (&)[N]) NXX_DELETE("open methods take one guess of a real type "
                                                                               "or a root estimate, not a braced list: write "
                                                                               "1.0, or pass {lo, hi} to a bracketing solver");

        template<class Self, class In>
            requires detail::accepts_v<Self, In>
        constexpr auto on(this const Self& self, In in)
        { return bound<Self, In> { self, std::move(in) }; }

        template<class Self, class T, std::size_t N>
        void on(this const Self&, const T (&)[N]) NXX_DELETE("open methods take one guess of a real type or a root "
                                                             "estimate, not a braced list: write 1.0, or pass {lo, hi} to a "
                                                             "bracketing solver");

        // A forwarding reference, for the reason given in bracketing_facade: by value, a C array would decay to a
        // pointer, which cl cannot order against the array overload above (C2668).
        template<class Self, class In>
            requires(detail::rejected_v<Self, std::remove_cvref_t<In>> && !std::is_integral_v<std::remove_cvref_t<In>>)
        void on(this const Self&, In&&) NXX_DELETE("open methods take a guess of a real type or a root estimate; "
                                                   "bracketing solvers take {lo, hi}");

        template<class Self, class I>
            requires std::is_integral_v<std::remove_cvref_t<I>>
        void on(this const Self&, I) NXX_DELETE("open methods take a guess of a real type: write 1.0, not 1");
    };

    // Searchers: a start window (bracket<T>, braced {lo, hi}, std::pair, the result of bracket<T>::make). They stop at
    // the first sign change or when the budget runs out, so they have no configurable stop criterion.
    struct search_facade : solver_facade
    {
        template<class Self, class C>
        void with_stop(this const Self&, C) NXX_DELETE("searchers have no configurable stop criterion: they stop at "
                                                       "the first sign change or when the budget runs out");

        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && detail::input_callable_v<Self, F, In> &&
                     detail::runnable_v<Self, F, In>)
        constexpr auto operator()(this const Self& self, const F& fn, const In& in)
        { return detail::run(self, fn, in); }

        // A braced list or a C array: N is deduced, so {x} is not taken as {x, 0} and {a, b, c} is not cut short.
        template<class Self, class F, real T, std::size_t N>
            requires(N == 2 && detail::ready_v<Self, F> && detail::input_callable_v<Self, F, std::pair<T, T>> &&
                     detail::runnable_v<Self, F, std::pair<T, T>>)
        constexpr auto operator()(this const Self& self, const F& fn, const T (&lo_hi)[N])
        { return detail::run(self, fn, std::pair<T, T> { lo_hi[0], lo_hi[1] }); }

        // A solver that does not implement the protocol (DESIGN §6.6) for {lo, hi}, and below for any other input.
        template<class Self, class F, real T, std::size_t N>
            requires(N == 2 && detail::ready_v<Self, F> && detail::input_callable_v<Self, F, std::pair<T, T>> &&
                     !detail::runnable_v<Self, F, std::pair<T, T>>)
        void operator()(this const Self&, const F&, const T (&)[N]) NXX_DELETE("this solver does not implement the solver protocol "
                                                                               "(DESIGN 6.6): it needs accepts_v, prepare(f, in), "
                                                                               "and id, options(), init, step, view, estimate, best "
                                                                               "and intrinsic for the problem prepare returns, and "
                                                                               "better_than(const Est&, const Est&) for its estimate type");

        template<class Self, class F, real T, std::size_t N>
            requires(N == 2 && detail::ready_v<Self, F> && !detail::input_callable_v<Self, F, std::pair<T, T>>)
        void operator()(this const Self&, const F&, const T (&)[N]) NXX_DELETE("the function cannot be called with the "
                                                                               "scalar type of the window");

        template<class Self, class F, class T, std::size_t N>
            requires(N != 2 || !real<T>)
        void operator()(this const Self&, const F&, const T (&)[N]) NXX_DELETE("a window has two ends of a real type: "
                                                                               "write {lo, hi}, for example {1.0, 2.0}");

        template<class Self, class F, class In>
            requires detail::rejected_v<Self, In>
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("searchers take a start window or a guess");

        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && !detail::input_callable_v<Self, F, In>)
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("the function cannot be called with the "
                                                                          "scalar type of the window");

        template<class Self, class F, class In>
            requires detail::incomplete_v<Self, F, In>
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("this solver does not implement the solver protocol "
                                                                          "(DESIGN 6.6): it needs accepts_v, prepare(f, in), "
                                                                          "and id, options(), init, step, view, estimate, best "
                                                                          "and intrinsic for the problem prepare returns, and "
                                                                          "better_than(const Est&, const Est&) for its estimate type");

        template<class Self, class In>
            requires detail::accepts_v<Self, In>
        constexpr auto on(this const Self& self, In in)
        { return bound<Self, In> { self, std::move(in) }; }

        template<class Self, real T, std::size_t N>
            requires(N == 2)
        constexpr auto on(this const Self& self, const T (&lo_hi)[N])
        { return bound<Self, std::pair<T, T>> { self, std::pair<T, T> { lo_hi[0], lo_hi[1] } }; }

        template<class Self, class T, std::size_t N>
            requires(N != 2 || !real<T>)
        void on(this const Self&, const T (&)[N]) NXX_DELETE("a window has two ends of a real type: write {lo, hi}, "
                                                             "for example {1.0, 2.0}");

        // A forwarding reference, for the reason given in bracketing_facade.
        template<class Self, class In>
            requires detail::rejected_v<Self, std::remove_cvref_t<In>>
        void on(this const Self&, In&&) NXX_DELETE("searchers take a start window or a guess");
    };
}    // namespace nxx

NXX_END_HEADER
