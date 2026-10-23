// E6 (Chapter 4): differentiating at the edge of a domain. The Shomate fit is valid on [500, 2000] K; a wrapper that
// returns NaN outside that range makes "evaluable" and "valid" the same thing, so the library can tell. At the lower
// end, Central<2> samples 500 K - h and fails with non_finite; Forward<2> and Forward<4> never sample below x and
// succeed, with relative errors at most 1e-8 ((1/3 + 8) eps^(2/3) is about 3e-10 for Forward<2>). Quoted with
// \snippet{04-stencils-edge}{validity-edge}.

#include "shomate.hpp"

#include <numerixx/deriv.hpp>

#include <cmath>
#include <limits>
#include <print>

namespace {

    double rel_error(std::expected<nxx::deriv::Estimate<double>, nxx::deriv::error> const &d, double exact)
    { return d ? std::abs(d->value - exact) / std::abs(exact) : 1.0; }

} // namespace

int main()
{
    using namespace nxx::deriv;
    // clang-format off
// <validity-edge>
// The fit is valid on [T_min, T_max] only. Outside it the
// wrapper returns NaN, so a rule that samples there fails
// instead of returning a plausible number.
auto const H = [](double T) {
    using book::nitrogen;
    if (T < nitrogen.T_min || T > nitrogen.T_max)
        return std::numeric_limits<double>::quiet_NaN();
    return book::shomate_H(nitrogen, T);
};
double const T_edge = book::nitrogen.T_min; // 500 K
auto const central = derivative<Central<2>>(H, T_edge);
auto const fwd2 = derivative<Forward<2>>(H, T_edge);
auto const fwd4 = derivative<Forward<4>>(H, T_edge);
// </validity-edge>
    // clang-format on
    double const cp = book::shomate_cp(book::nitrogen, T_edge);
    bool const central_failed = !central && central.error() == error::non_finite;
    double const e2 = rel_error(fwd2, cp);
    double const e4 = rel_error(fwd4, cp);
    std::println("cp exact at {} K  {:.12f} J/(mol K)", T_edge, cp);
    std::println("Central<2>   {}", central ? "value" : "non_finite (samples below the range)");
    if (fwd2) std::println("Forward<2>   {:.12f}  rel. error {:.1e}", fwd2->value, e2);
    if (fwd4) std::println("Forward<4>   {:.12f}  rel. error {:.1e}", fwd4->value, e4);
    bool const ok = central_failed && fwd2 && fwd4 && e2 <= 1e-8 && e4 <= 1e-8;
    return ok ? 0 : 1;
}
