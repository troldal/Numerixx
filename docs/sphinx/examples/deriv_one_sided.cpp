// The program of docs/sphinx/deriv/one-sided.md: a function defined for x >= 0 only, 1 / (1 + x) with NaN below 0,
// whose derivative is -1 / (1 + x)^2: at the edge 0, and at 1e-2, 1e-6 and 1e-10 from it (the near-edge table). The
// page shows the parts between the markers and the output in the raw strings; main returns non-zero unless the
// program prints that output and what the page says holds. fn uses only + - * /, so the output is the same on every
// target that does not contract a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>

#include <cmath>
#include <limits>
#include <print>

namespace dv = nxx::deriv;
// </includes>
// clang-format on

#include "printed.hpp"

#include <cstddef>
#include <initializer_list>
#include <string>
#include <vector>

namespace {

    // The relative error of a result, or -1 for a failure.
    template<typename R>
    double relative_error(R const &d, double exact)
    { return d ? std::abs(d->value - exact) / std::abs(exact) : -1.0; }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <edge>
    // 1 / (1 + x), defined for x >= 0 only: NaN below 0.
    auto const fn = [](double x) {
        if (x < 0.0) return std::numeric_limits<double>::quiet_NaN();
        return 1.0 / (1.0 + x);
    };
    auto const central = dv::derivative(fn, 0.0);
    auto const forward = dv::derivative<dv::Forward<2>>(fn, 0.0);
    auto const forward_ridders = dv::derivative<dv::ForwardRidders>(fn, 0.0);
    if (!central) std::println("Central<2>:     {}", central.error());
    if (forward) std::println("Forward<2>:     {:.15f}", forward->value);
    if (forward_ridders) std::println("ForwardRidders: {:.15f}, abs_error {:.1e}", forward_ridders->value,
                                      *forward_ridders->abs_error);
    // </edge>

    // <kink>
    // A formula that goes on below 0 with another slope: no NaN tells the methods where the edge is.
    auto const kink = [](double x) { return x >= 0.0 ? 1.0 / (1.0 + x) : 1.0 - 3.0 * x; };
    auto const both_sides = dv::derivative<dv::Ridders>(kink, 0.0);
    auto const one_side = dv::derivative<dv::ForwardRidders>(kink, 0.0);
    if (both_sides)
        std::println("Ridders:        {:.9f}, abs_error {:.1e}", both_sides->value, *both_sides->abs_error);
    if (one_side)
        std::println("ForwardRidders: {:.9f}, abs_error {:.1e}", one_side->value, *one_side->abs_error);
    // </kink>
    // clang-format on

    // The near-edge table: Ridders' method and ForwardRidders at 1e-2, 1e-6 and 1e-10 from the edge, each with its
    // relative error, its estimate relative to |f'|, and the evaluations of fn it made.
    struct Row {
        double distance;
        char const *method;
        double value;
        double error;
        double estimate;
        int evaluations;
        int nans; // the evaluations at which fn returned NaN: points past the edge
    };
    std::vector<Row> rows;
    for (double const d : {1e-2, 1e-6, 1e-10}) {
        // -1 / (1 + d)^2 in three roundings, within 1.5 ulps: a relative 3.3e-16, a hundred times below the smallest
        // error in the table, 4e-14, so it cannot change the printed digit by more than a unit in rare ties, and the
        // checks below allow for it.
        double const exact = -1.0 / ((1.0 + d) * (1.0 + d));
        auto const add = [&]<typename M>(char const *name) {
            int calls = 0;
            int nans = 0;
            auto const counted = [&](double t) {
                ++calls;
                double const y = fn(t);
                if (std::isnan(y)) ++nans;
                return y;
            };
            auto const r = dv::derivative<M>(counted, d);
            if (r && r->abs_error)
                rows.push_back({.distance = d, .method = name, .value = r->value, .error = relative_error(r, exact),
                                .estimate = *r->abs_error / std::abs(exact), .evaluations = calls, .nans = nans});
        };
        add.operator()<dv::Ridders>("Ridders");
        add.operator()<dv::ForwardRidders>("ForwardRidders");
    }
    std::println("distance  method           error  estimate  evaluations");
    for (Row const &r : rows)
        std::println("{:<8.0e}  {:<15} {:>6.0e}  {:>8.0e}  {:>11}", r.distance, r.method, r.error, r.estimate,
                     r.evaluations);

    constexpr std::string_view edge_output = R"(
Central<2>:     non_finite
Forward<2>:     -0.999999999886047
ForwardRidders: -0.999999999999835, abs_error 4.8e-12
)";
    constexpr std::string_view kink_output = R"(
Ridders:        -1.999999991, abs_error 1.2e-08
ForwardRidders: -1.000000000, abs_error 4.8e-12
)";
    constexpr std::string_view distances_output = R"(
distance  method           error  estimate  evaluations
1e-02     Ridders          4e-14     2e-13           16
1e-02     ForwardRidders   5e-13     5e-12           10
1e-06     Ridders          2e-11     9e-10           36
1e-06     ForwardRidders   3e-13     5e-12           10
1e-10     Ridders          2e-10     8e-06           64
1e-10     ForwardRidders   1e-13     5e-12           10
)";
    std::string const expected = std::string(docs::quoted(edge_output)) + std::string(docs::quoted(kink_output))
                                 + std::string(docs::quoted(distances_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN fails it.
    if (!forward || !forward_ridders || !forward_ridders->abs_error) return 2;
    // At the edge, the central methods fail: their points reach below 0, where fn is NaN.
    if (!(central == std::unexpected(dv::error::non_finite))) return 3;
    if (!(dv::derivative<dv::Ridders>(fn, 0.0) == std::unexpected(dv::error::non_finite))) return 4;
    // Forward<2>: f'(0) = -1 to its model's accuracy, h^2 / 3 |f'''| = h^2 / 3 * 6 with h = eps^(1/3) = 6.06e-6, about
    // 7e-11, plus rounding about S eps / h = 4 * 2.2e-16 / 6e-6 = 1.5e-10; the page says "about 10 digits".
    if (!(relative_error(forward, -1.0) < 3e-10)) return 5;
    // ForwardRidders: its estimate covers the actual error, below 1e-12 (the page: "about 13 digits").
    if (!(std::abs(forward_ridders->value + 1.0) <= *forward_ridders->abs_error && *forward_ridders->abs_error < 1e-11
          && std::abs(forward_ridders->value + 1.0) < 1e-12))
        return 6;
    // At 1e-6, the central difference's default step, 6.06e-6, still reaches below 0.
    if (!(dv::derivative(fn, 1e-6) == std::unexpected(dv::error::non_finite))) return 8;
    // "On a function that changes fast near its edge, such as sqrt(1 - x) near 1, Ridders' method is the more accurate
    // one", and a one-sided stencil with its default step fails there: g(x) = sqrt(1 - x), defined up to 1,
    // g'(x) = -1 / (2 sqrt(1 - x)). At 0.999 Ridders beats BackwardRidders; at 1 - 1e-7, Backward<2>'s step, 6e-6,
    // spans sixty times the distance to the edge, and its result is more than 10% off, with no estimate. 1 - x is
    // exact (Sterbenz), so each reference is within about 2 ulps, far below the errors compared.
    auto const g = [](double x) {
        if (x > 1.0) return std::numeric_limits<double>::quiet_NaN();
        return std::sqrt(1.0 - x);
    };
    double const at = 0.999;
    double const g_exact = -0.5 / std::sqrt(1.0 - at);
    auto const g_central = dv::derivative<dv::Ridders>(g, at);
    auto const g_backward = dv::derivative<dv::BackwardRidders>(g, at);
    if (!g_central || !g_backward
        || !(std::abs(g_central->value - g_exact) < std::abs(g_backward->value - g_exact)))
        return 13;
    double const close = 1.0 - 1e-7;
    auto const g_stencil = dv::derivative<dv::Backward<2>>(g, close);
    double const g_close = -0.5 / std::sqrt(1.0 - close);
    if (!g_stencil || !(std::abs(g_stencil->value - g_close) > 0.1 * std::abs(g_close))) return 14;
    // "On a function whose values happen to vanish at the points a method tries, [the estimate] can even be 0 for a
    // wrong result" (Codex's case): u (1 - (64 u)^2) (1 - (128 u)^2), u = t - 1/64, NaN below 0, at 1/64, where the
    // derivative is 1. Ridders' first step 1/32 reaches below 0; its next steps, 1/64 and 1/128, land on zeros of the
    // polynomial, so it returns 0 with an estimate of 0. ForwardRidders returns 1.
    auto const vanishing = [](double t) {
        if (t < 0.0) return std::numeric_limits<double>::quiet_NaN();
        double const u = t - 1.0 / 64.0;
        return u * (1.0 - (64.0 * u) * (64.0 * u)) * (1.0 - (128.0 * u) * (128.0 * u));
    };
    auto const fooled = dv::derivative<dv::Ridders>(vanishing, 1.0 / 64.0);
    auto const right = dv::derivative<dv::ForwardRidders>(vanishing, 1.0 / 64.0);
    if (!fooled || !fooled->abs_error || !(fooled->value == 0.0 && *fooled->abs_error == 0.0)) return 15;
    if (!right || !right->abs_error || !(std::abs(right->value - 1.0) <= *right->abs_error)) return 16;
    // The kink: Ridders returns about the mean of the slopes -3 and -1, wrong, with a small estimate; ForwardRidders
    // sees only the right-hand formula and returns -1 within its estimate.
    if (!both_sides || !both_sides->abs_error || !one_side || !one_side->abs_error) return 9;
    if (!(std::abs(both_sides->value + 2.0) < 1e-6 && *both_sides->abs_error < 1e-6)) return 10;
    if (!(std::abs(one_side->value + 1.0) <= *one_side->abs_error)) return 11;

    // The near-edge table: six rows, Ridders' method then ForwardRidders at each distance.
    if (!(rows.size() == 6)) return 17;
    for (std::size_t i = 0; i < rows.size(); i += 2) {
        Row const &central_row = rows[i];
        Row const &one_sided_row = rows[i + 1];
        // Both estimates cover their errors, allowing 3.3e-16 for the reference's rounding.
        if (!(central_row.error <= central_row.estimate + 3.3e-16
              && one_sided_row.error <= one_sided_row.estimate + 3.3e-16))
            return 18;
        // ForwardRidders stays below 1e-12 in error and 1e-11 in estimate, at 10 evaluations or fewer.
        if (!(one_sided_row.error < 1e-12 && one_sided_row.estimate < 1e-11 && one_sided_row.evaluations <= 10))
            return 19;
    }
    // Far from the edge, at 1e-2, both work, and Ridders' method is the more accurate.
    if (!(rows[0].error < 1e-12 && rows[0].error < rows[1].error)) return 20;
    // Near it, at 1e-10, Ridders' method restarts below each step that met NaN until it has made 64 evaluations, the
    // most it may make for f' (deriv.hpp, Ridders), with an estimate of about 1e-5; ForwardRidders is the more
    // accurate.
    if (!(rows[4].evaluations == 64 && rows[4].estimate > 3e-6 && rows[4].estimate < 3e-5)) return 21;
    if (!(rows[5].error < rows[4].error)) return 22;
    // The page's numbers at 1e-10: Ridders' error between 1e-10 and 1e-9; the one-sided method's below 1e-12, at 10
    // evaluations.
    if (!(rows[4].error > 1e-10 && rows[4].error < 1e-9 && rows[5].error < 1e-12 && rows[5].evaluations == 10))
        return 23;
    // "It restarts below every step that reaches past the edge, and the closer x is, the more often": at every distance
    // Ridders' method met NaN, and more often the closer it was.
    if (!(rows[0].nans > 0 && rows[2].nans > rows[0].nans && rows[4].nans > rows[2].nans)) return 24;
    // The rule: "keep the one with the smaller estimate if the two differ by less than the sum of their estimates". At
    // 1e-10 they do (estimates relative to |f'| = 1 within 2e-10, so compared as absolute values), the smaller
    // estimate is ForwardRidders', and its error is the smaller.
    if (!(std::abs(rows[4].value - rows[5].value) < rows[4].estimate + rows[5].estimate
          && rows[5].estimate < rows[4].estimate && rows[5].error < rows[4].error))
        return 25;
    // A formula that cancels (the user's decision U3): sqrt(1 - x^2) at 0.999, NaN above 1, loses up to about 125 eps
    // relative in 1 - x^2 (about 60 eps after the square root), and BackwardRidders' estimate falls below its error
    // (Ridders' happens to cover it, hence the page's "makes an estimate too small", not both); written as
    // sqrt((1 - x)(1 + x)), both Ridders' and BackwardRidders' estimates cover their errors. The reference,
    // -x / sqrt((1 - x)(1 + x)), has 1 - x exact and a few roundings, about 1e-15 relative, far below the errors
    // compared.
    auto const cancelling = [](double x) {
        if (x > 1.0) return std::numeric_limits<double>::quiet_NaN();
        return std::sqrt(1.0 - x * x);
    };
    auto const factored = [](double x) {
        if (x > 1.0) return std::numeric_limits<double>::quiet_NaN();
        return std::sqrt((1.0 - x) * (1.0 + x));
    };
    double const c_exact = -at / std::sqrt((1.0 - at) * (1.0 + at));
    auto const missed = dv::derivative<dv::BackwardRidders>(cancelling, at);
    if (!missed || !missed->abs_error || !(std::abs(missed->value - c_exact) > *missed->abs_error)) return 26;
    auto const cured_central = dv::derivative<dv::Ridders>(factored, at);
    auto const cured_one_sided = dv::derivative<dv::BackwardRidders>(factored, at);
    if (!cured_central || !cured_central->abs_error || !cured_one_sided || !cured_one_sided->abs_error
        || !(std::abs(cured_central->value - c_exact) <= *cured_central->abs_error)
        || !(std::abs(cured_one_sided->value - c_exact) <= *cured_one_sided->abs_error))
        return 27;
    return 0;
}
