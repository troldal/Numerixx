// The program of docs/sphinx/how-to/roots-derivatives.md: giving f' and f'' to Newton and Halley. The page shows the
// parts between the markers; main returns non-zero unless what the page says holds. Counts are the independent
// prototypes' (docs/roots/design.md, Evidence, Phases 2 and 3: .scratch/roots2/proto2.py, .scratch/roots3/proto3.py),
// for fn written as here and builds without FMA contraction.

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>

#include <print>
#include <tuple>
#include <utility>
// </includes>
// clang-format on

#include <cmath>
#include <expected>
#include <limits>
#include <type_traits>

namespace {

    namespace rt = nxx::roots;

    // fn wrapped to count its calls; it returns what fn returns.
    template<typename F>
    struct Counted
    {
        F fn;
        int *calls;
        auto operator()(double x) const
        {
            ++*calls;
            return fn(x);
        }
    };

    // Within abs_error + c of the Wallis root 2.0945514815423265915 (mpmath): c = 2e-16 for fn's rounding over
    // |f'| = 11.2, plus half an ulp of the root for the reference in a double, 2.3e-16; and abs_error <= 2 tol(value).
    bool accurate(std::expected<rt::Root<double>, rt::error> const &r)
    {
        double const eps = std::numeric_limits<double>::epsilon();
        return r && std::abs(r->value - 2.0945514815423265915) <= r->abs_error + 2e-16 + 2.3e-16
               && r->abs_error <= 2 * (std::numeric_limits<double>::denorm_min() + 2 * eps * r->value);
    }

} // namespace

int main()
{
    // clang-format off
    // <pair>
    // x^3 - 2x - 5 and its derivative, in one call
    auto const fdf = [](double x) { return std::pair{(x * x - 2.0) * x - 5.0, 3.0 * x * x - 2.0}; };
    auto const newton = rt::solve<rt::Newton>(fdf, rt::Bracket{2.0, 3.0});
    auto const newton_guess = rt::solve<rt::Newton>(fdf, 2.0);
    // </pair>

    // <tuple>
    // with f'' for Halley
    auto const fdf2 = [](double x) {
        return std::tuple{(x * x - 2.0) * x - 5.0, 3.0 * x * x - 2.0, 6.0 * x};
    };
    auto const halley = rt::solve<rt::Halley>(fdf2, 2.0);
    // </tuple>

    // <separate>
    auto const f = [](double x) { return (x * x - 2.0) * x - 5.0; };
    auto const df = [](double x) { return 3.0 * x * x - 2.0; };
    auto const separate = rt::solve<rt::Newton>(f, df, rt::Bracket{2.0, 3.0});
    // </separate>

    // <numerical>
    auto const numerical = rt::solve<rt::Newton>(f, nxx::deriv::derivative_of(f), rt::Bracket{2.0, 3.0});
    // </numerical>
    // clang-format on

    if (newton) std::println("Newton {:.17g}, Halley {:.17g}", newton->value, halley ? halley->value : 0.0);

    if (!accurate(newton) || !accurate(newton_guess) || !accurate(halley) || !accurate(separate)) return 1;
    if (!accurate(numerical)) return 2;
    // The separate callables are the pair, joined: the same result.
    if (!(separate == newton)) return 3;

    // Calls of fn: one per point for the pair and the tuple; with derivative_of (Central<2>) f at the point and at
    // its two nodes, 3 per point.
    int calls = 0;
    auto const counted = [&](auto const &fn) {
        calls = 0;
        return Counted<std::decay_t<decltype(fn)>>{fn, &calls};
    };
    (void)rt::solve<rt::Newton>(counted(fdf), rt::Bracket{2.0, 3.0});
    if (!(calls == 7)) return 4;
    (void)rt::solve<rt::Newton>(counted(fdf), 2.0);
    if (!(calls == 6)) return 5;
    (void)rt::solve<rt::Halley>(counted(fdf2), 2.0);
    if (!(calls == 5)) return 6;
    auto const counted_f = counted(f);
    (void)rt::solve<rt::Newton>(counted_f, nxx::deriv::derivative_of(counted_f), rt::Bracket{2.0, 3.0});
    if (!(calls == 21)) return 7;
    calls = 0;
    (void)rt::solve<rt::Newton>(counted_f, nxx::deriv::derivative_of(counted_f), 2.0);
    if (!(calls == 18)) return 8;
    // Halley with f' and f'' both numerical, Central<2>: 6 calls per point (3 for f', 3 for f''), 30 from 2.
    calls = 0;
    (void)rt::solve<rt::Halley>(counted_f, nxx::deriv::derivative_of(counted_f),
                                nxx::deriv::derivative_of<nxx::deriv::Central<2>, 2>(counted_f), 2.0);
    if (!(calls == 30)) return 10;
    // Secant, value-only, from the same guess: 8 (Phase 2's prototype).
    calls = 0;
    (void)rt::solve<rt::Secant>(counted_f, 2.0);
    if (!(calls == 8)) return 9;
    return 0;
}
