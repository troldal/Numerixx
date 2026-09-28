// One-shot differentiation (DESIGN §7.1): diff(f, x[, stencil[, step]]) -> expected<T, fault<UE>>, where UE is the
// callback's own error type. The default is central_1_2 with the optimal relative step.
//
// Deferred to phase 2: forward/backward/second conveniences, diff_with_error, ridders and the mixed partials.
#pragma once

#include <numerixx/core/callable.hpp>
#include <numerixx/core/error.hpp>
#include <numerixx/core/math.hpp>
#include <numerixx/core/scalar.hpp>
#include <numerixx/deriv/stencil.hpp>
#include <numerixx/deriv/step.hpp>

#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>

NXX_BEGIN_HEADER

namespace nxx::deriv
{
    // Default template arguments let diff(f, x) deduce (a defaulted stencil parameter alone cannot, DESIGN §7.1).
    template<class F, real T, int O = 1, int A = 2, std::size_t N = 2, class H = optimal>
    constexpr auto diff(const F& fn, T x, const stencil<O, A, N>& s = central_1_2, H h = {})
        -> std::expected<T, fault<callback_error_t<F, T>>>
    {
        using UE   = callback_error_t<F, T>;
        const T hh = detail::resolve<O, A>(h, x);
        if (!(hh > T(0)) || !math::isfinite(hh)) return std::unexpected(fault<UE> { errc::invalid_input, 0, {} });
        T             acc(0);
        std::uint32_t evals = 0;
        for (std::size_t i = 0; i < N; ++i) {
            if (s.weight[i] == 0) continue;
            auto y = nxx::evaluate(fn, x + T(s.offset[i]) * hh);
            if (!y) {
                fault<UE> e = y.error();
                e.evals += evals;
                return std::unexpected(e);
            }
            evals += cost_of(fn);
            acc += T(s.weight[i]) * *y;
        }
        // acc / (denominator h^O) in one division where h^O is representable; otherwise one division per power, which
        // keeps every intermediate representable whenever the derivative is (h^2 overflows for |x| > 1e158 in double,
        // and would turn a finite f'' into 0; it underflows below 1e-155).
        T den = T(s.denominator);
        for (int k = 0; k < O; ++k) den *= hh;
        T d;
        if (math::isfinite(den) && den >= (std::numeric_limits<T>::min)())
            d = acc / den;
        else {
            d = acc / T(s.denominator);
            for (int k = 0; k < O; ++k) d /= hh;
        }
        if (!math::isfinite(d)) return std::unexpected(fault<UE> { errc::non_finite_value, evals, {} });
        return d;
    }

    // f'(x) with central_1_2.
    template<class F, real T, class H = optimal>
    constexpr auto central(const F& fn, T x, H h = {})
    { return nxx::deriv::diff(fn, x, central_1_2, h); }
}    // namespace nxx::deriv

NXX_END_HEADER
