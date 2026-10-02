// Derivatives as functions (DESIGN §6.12, D21):
//   derivative_of(f[, stencil[, step]])  a derivative_fn: x -> expected<T, fault<UE>>, so it plugs into a solver as a
//                                         callback without nesting errors, keeping f's cause (Newton's UE is f's UE).
//                                         cost_of(it) = the stencil's non-zero points times cost_of(f).
//   numeric{stencil, step}                a derivative policy: bind(f) -> derivative_fn, for Newton in a curried chain
//                                         (the function is not known when the chain is built).
// derivative_fn holds f in a copyable box: copy-assignable even when f is a lambda with captures.
#pragma once

#include <numerixx/core/callable.hpp>
#include <numerixx/core/scalar.hpp>
#include <numerixx/deriv/diff.hpp>
#include <numerixx/deriv/stencil.hpp>
#include <numerixx/deriv/step.hpp>

#include <concepts>
#include <cstdint>
#include <type_traits>
#include <utility>

NXX_BEGIN_HEADER

namespace nxx::deriv
{
    template<class F, class S = std::remove_cvref_t<decltype(central_1_2)>, class H = optimal>
    class derivative_fn
    {
        nxx::detail::copyable_box<F> fn_;
        S                            stencil_;
        H                            step_;

    public:
        constexpr derivative_fn(F fn, S s, H h) : fn_(std::move(fn)), stencil_(s), step_(h) {}

        template<real T>
        constexpr auto operator()(const T& x) const
        { return nxx::deriv::diff(*fn_, x, stencil_, step_); }

        constexpr std::uint32_t evaluation_cost() const noexcept { return stencil_.nonzero_points() * nxx::cost_of(*fn_); }
        constexpr const F&      function() const noexcept { return *fn_; }
    };

    template<class F, class S = std::remove_cvref_t<decltype(central_1_2)>, class H = optimal>
    constexpr auto derivative_of(F fn, S s = central_1_2, H h = {})
    { return derivative_fn<F, S, H> { std::move(fn), s, h }; }

    // Constructors rather than default member initialisers: numeric<S, H> is default-constructible only when S is
    // central_1_2's type and H is default-constructible, and asking (std::default_initializable, a copyable box) is not a
    // hard error for other S and H.
    template<class S = std::remove_cvref_t<decltype(central_1_2)>, class H = optimal>
    struct numeric
    {
        S s;
        H h;

        constexpr numeric()
            requires(std::same_as<S, std::remove_cvref_t<decltype(central_1_2)>> && std::default_initializable<H>)
            : s(central_1_2),
              h()
        {}
        constexpr explicit numeric(S st)
            requires std::default_initializable<H>
            : s(st),
              h()
        {}
        constexpr numeric(S st, H hh) : s(st), h(hh) {}

        template<class F>
        constexpr auto bind(const F& fn) const
        { return nxx::deriv::derivative_of(fn, s, h); }
    };
}    // namespace nxx::deriv

NXX_END_HEADER
