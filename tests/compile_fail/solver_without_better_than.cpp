// Compile-fail (DESIGN §6.6): a solver whose failure estimate type has no order. nxx::better_than is the customisation
// point for that order, found by ADL next to the estimate type, and iterative_solver_for requires it, so the facade's
// protocol reason names the missing better_than. Before, the driver fell back to merit_of(e) < merit_of(best), and
// without a merit_of either the call was a hard error inside nxx::iterate.
#include <numerixx/roots.hpp>

#include <cstdint>
#include <expected>
#include <optional>
#include <type_traits>

namespace mine
{
    struct halving_estimate
    {
        double x;
        double fx;
    };

#ifdef NUMERIXX_CF_CONTROL
    constexpr bool better_than(const halving_estimate& a, const halving_estimate& b) noexcept
    { return nxx::math::abs(a.fx) < nxx::math::abs(b.fx); }
#endif

    struct halving_state
    {
        double        x;
        double        fx;
        std::uint32_t nfev;
    };

    // A complete solver protocol over an open facade: x halves at every step.
    struct halving : nxx::open_facade
    {
        static constexpr nxx::algo      id    = nxx::algo::user_first;
        static constexpr nxx::view_kind views = nxx::view_kind::point;
        template<class In>
        static constexpr bool accepts_v = std::is_same_v<In, double>;

        nxx::options<nxx::never> opt { {}, nxx::max_iterations { 8 } };

        constexpr const nxx::options<nxx::never>& options() const noexcept { return opt; }

        template<class F>
        constexpr auto prepare(const F& fn, double x0) const -> std::expected<nxx::problem<F, double>, nxx::failure<halving_estimate>>
        { return nxx::problem<F, double> { fn, x0, 0 }; }

        template<class P>
        constexpr auto init(const P& p) const -> std::expected<halving_state, nxx::failure<halving_estimate>>
        { return halving_state { p.in, p.f(p.in), 1 }; }

        template<class P>
        constexpr auto step(const P& p, const halving_state& s) const -> std::expected<halving_state, nxx::fault<>>
        { return halving_state { s.x / 2.0, p.f(s.x / 2.0), s.nfev + 1 }; }

        constexpr nxx::roots::point_view<double>  view(const halving_state& s) const noexcept { return { s.x, s.fx }; }
        constexpr halving_estimate                estimate(const halving_state& s) const noexcept { return { s.x, s.fx }; }
        constexpr halving_estimate                best(const halving_state& s) const noexcept { return { s.x, s.fx }; }
        constexpr std::optional<nxx::stop_reason> intrinsic(const halving_state&) const noexcept { return std::nullopt; }
    };
}    // namespace mine

int main()
{
    auto f = [](double x) { return x - 1.0; };
    (void)mine::halving {}(f, 1.0);
}
