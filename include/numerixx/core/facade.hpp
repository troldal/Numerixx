// Solver configuration and the family facades (DESIGN §6.6).
//
// Every solver holds one options aggregate: stop criterion, iteration budget, derivative source, projection and
// observer. The builders with_stop, with_budget, with_derivative, with_projection and with_observer are generic: they
// rebind that aggregate and return a new solver value, so configuration is order-independent and the old value is
// unchanged. User callables are held in copyable boxes, so a solver that holds a lambda stays copy-assignable.
//
// The family facades (bracketing, open, search) are deducing-this bases, no CRTP. Each constrains operator() and .on()
// on the inputs its solvers accept and deletes everything else with a reason, so std::is_invocable_v is false (not a
// hard error) for a solver given the wrong kind of input.
#pragma once

#include <numerixx/config.hpp>
#include <numerixx/core/callable.hpp>
#include <numerixx/core/criteria.hpp>
#include <numerixx/core/iterate.hpp>
#include <numerixx/core/refined.hpp>
#include <numerixx/core/scalar.hpp>

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

    // The builders shared by every solver. A solver provides options(), rebuild(options) and the kind of its views,
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

        template<class Self>
        constexpr auto with_budget(this const Self& self, max_iterations budget)
        {
            auto o   = self.options();
            o.budget = budget;
            return self.rebuild(std::move(o));
        }

        template<class Self, class C>
            requires criterion_for_v<C, Self::views>
        constexpr auto with_stop(this const Self& self, C stop)
        {
            const auto& o = self.options();
            using O       = std::remove_cvref_t<decltype(o)>;
            using O2      = options<C, typename O::derivative_type, typename O::projection_type, typename O::observer_type>;
            return self.rebuild(O2 { std::move(stop), o.budget, o.derivative, o.projection, o.observer });
        }

        template<class Self, class C>
            requires(!criterion_for_v<C, Self::views>)
        void with_stop(this const Self&, C)
            NXX_DELETE("this criterion does not apply to this solver: bracketing methods converge on the enclosure "
                       "(width_tol, floored_width), open methods on successive iterates (x_tol, step_tol)");

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
        void with_projection(this const Self&, P)
            NXX_DELETE("projection applies to open methods (secant, newton): a bracketing method keeps every iterate inside its bracket");

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
        // take.
        template<class F>
            requires std::is_invocable_v<const S&, const F&, const In&>
        constexpr auto operator()(const F& fn) const
        { return solver_(fn, in_); }

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

        template<class S, class In>
        inline constexpr bool accepts_v = S::template accepts_v<std::remove_cvref_t<In>>;
        template<class S, class F>
        inline constexpr bool ready_v = S::template ready_v<F>;
        // Whether F can be called with the scalar type of the input (only asked once accepts_v holds). F decays: for a
        // plain function F is a function type, and the solver's `const F&` would then qualify a function type, which MSVC
        // warns about (C4180); a function and a pointer to it are invocable alike.
        template<class S, class F, class In>
        inline constexpr bool input_callable_v = S::template callable_v<std::decay_t<F>, std::remove_cvref_t<In>>;
        // An input that is not accepted, and not a C array (MSVC cannot order the array overload against this one).
        template<class S, class In>
        inline constexpr bool rejected_v = !accepts_v<S, In> && !std::is_array_v<std::remove_cvref_t<In>>;
    }    // namespace detail

    // Bracketing methods: a bracket<T>, a braced {lo, hi}, a std::pair, the result of bracket<T>::make, a sign_bracket
    // or a search result.
    struct bracketing_facade : solver_facade
    {
        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && detail::input_callable_v<Self, F, In>)
        constexpr auto operator()(this const Self& self, const F& fn, const In& in)
        { return detail::run(self, fn, in); }

        template<class Self, class F, real T>
            requires(detail::ready_v<Self, F> && detail::input_callable_v<Self, F, std::pair<T, T>>)
        constexpr auto operator()(this const Self& self, const F& fn, const T (&lo_hi)[2])
        { return detail::run(self, fn, std::pair<T, T> { lo_hi[0], lo_hi[1] }); }

        template<class Self, class F, class In>
            requires detail::rejected_v<Self, In>
        void operator()(this const Self&, const F&, const In&)
            NXX_DELETE("bracketing solvers need a bracket: pass {lo, hi}, nxx::bracket<T>::make(a, b), or a search result");

        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && !detail::input_callable_v<Self, F, In>)
        void operator()(this const Self&, const F&, const In&)
            NXX_DELETE("the function cannot be called with the scalar type of the bracket");

        template<class Self, class In>
            requires detail::accepts_v<Self, In>
        constexpr auto on(this const Self& self, In in)
        { return bound<Self, In> { self, std::move(in) }; }

        template<class Self, real T>
        constexpr auto on(this const Self& self, const T (&lo_hi)[2])
        { return bound<Self, std::pair<T, T>> { self, std::pair<T, T> { lo_hi[0], lo_hi[1] } }; }

        template<class Self, class In>
            requires detail::rejected_v<Self, In>
        void on(this const Self&, In)
            NXX_DELETE("bracketing solvers need a bracket: pass {lo, hi}, nxx::bracket<T>::make(a, b), or a search result");
    };

    // Open methods: a guess of a real type, or a root estimate (seeded with its x and f(x): no re-evaluation).
    struct open_facade : solver_facade
    {
        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && detail::input_callable_v<Self, F, In>)
        constexpr auto operator()(this const Self& self, const F& fn, const In& in)
        { return detail::run(self, fn, in); }

        template<class Self, class F, class In>
            requires(detail::rejected_v<Self, In> && !std::is_integral_v<std::remove_cvref_t<In>>)
        void operator()(this const Self&, const F&, const In&)
            NXX_DELETE("open methods take a guess of a real type or a root estimate; bracketing solvers take {lo, hi}");

        template<class Self, class F, class I>
            requires std::is_integral_v<std::remove_cvref_t<I>>
        void operator()(this const Self&, const F&, const I&) NXX_DELETE("open methods take a guess of a real type: write 1.0, not 1");

        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && !detail::input_callable_v<Self, F, In>)
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("the function cannot be called with the type of the guess");

        template<class Self, class In>
            requires detail::accepts_v<Self, In>
        constexpr auto on(this const Self& self, In in)
        { return bound<Self, In> { self, std::move(in) }; }

        template<class Self, class In>
            requires(detail::rejected_v<Self, In> && !std::is_integral_v<std::remove_cvref_t<In>>)
        void on(this const Self&, In)
            NXX_DELETE("open methods take a guess of a real type or a root estimate; bracketing solvers take {lo, hi}");

        template<class Self, class I>
            requires std::is_integral_v<std::remove_cvref_t<I>>
        void on(this const Self&, I) NXX_DELETE("open methods take a guess of a real type: write 1.0, not 1");
    };

    // Searchers: a start window (bracket<T>, braced {lo, hi}, std::pair, the result of bracket<T>::make). They stop at
    // the first sign change or when the budget runs out, so they have no configurable stop criterion.
    struct search_facade : solver_facade
    {
        template<class Self, class C>
        void with_stop(this const Self&, C)
            NXX_DELETE("searchers have no configurable stop criterion: they stop at the first sign change or when the budget runs out");

        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && detail::input_callable_v<Self, F, In>)
        constexpr auto operator()(this const Self& self, const F& fn, const In& in)
        { return detail::run(self, fn, in); }

        template<class Self, class F, real T>
            requires(detail::ready_v<Self, F> && detail::input_callable_v<Self, F, std::pair<T, T>>)
        constexpr auto operator()(this const Self& self, const F& fn, const T (&lo_hi)[2])
        { return detail::run(self, fn, std::pair<T, T> { lo_hi[0], lo_hi[1] }); }

        template<class Self, class F, class In>
            requires detail::rejected_v<Self, In>
        void operator()(this const Self&, const F&, const In&) NXX_DELETE("searchers take a start window or a guess");

        template<class Self, class F, class In>
            requires(detail::accepts_v<Self, In> && detail::ready_v<Self, F> && !detail::input_callable_v<Self, F, In>)
        void operator()(this const Self&, const F&, const In&)
            NXX_DELETE("the function cannot be called with the scalar type of the window");

        template<class Self, class In>
            requires detail::accepts_v<Self, In>
        constexpr auto on(this const Self& self, In in)
        { return bound<Self, In> { self, std::move(in) }; }

        template<class Self, real T>
        constexpr auto on(this const Self& self, const T (&lo_hi)[2])
        { return bound<Self, std::pair<T, T>> { self, std::pair<T, T> { lo_hi[0], lo_hi[1] } }; }

        template<class Self, class In>
            requires detail::rejected_v<Self, In>
        void on(this const Self&, In) NXX_DELETE("searchers take a start window or a guess");
    };
}    // namespace nxx

NXX_END_HEADER
