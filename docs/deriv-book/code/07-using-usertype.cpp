// E23 (Chapter 7, Sidebar 7.1; Chapter 8, section 8.7): a real type of one's own, built to the requirements the book
// lists and nothing more, so that the list is tested. Real wraps a double; std::numeric_limits<Real> is specialised
// with exactly the members that deriv.hpp and math.hpp use (is_specialized, is_integer, is_bounded, is_exact, radix,
// digits, epsilon(), max(), denorm_min()) and does not inherit numeric_limits<double>, whose members would hide a
// missing one. Real has no default constructor: the library does not need one. The calls differentiate e^x at 1.

#include <numerixx/deriv.hpp>

#include <cmath>
#include <compare>
#include <limits>
#include <numbers>
#include <print>

using nxx::deriv::derivative;
using nxx::deriv::Ridders;

// clang-format off
// <usertype-type>
// A real type of your own, here a wrapper of double: the
// arithmetic, a total order and construction from int.
struct Real
{
    double v;
    constexpr Real(int i) : v(i) {}
    constexpr explicit Real(double d) : v(d) {}
    friend constexpr Real operator+(Real a, Real b)
    { return Real(a.v + b.v); }
    friend constexpr Real operator-(Real a, Real b)
    { return Real(a.v - b.v); }
    friend constexpr Real operator*(Real a, Real b)
    { return Real(a.v * b.v); }
    friend constexpr Real operator/(Real a, Real b)
    { return Real(a.v / b.v); }
    friend constexpr Real operator-(Real a)
    { return Real(-a.v); }
    friend constexpr auto operator<=>(Real, Real) = default;
};
// </usertype-type>

// <usertype-limits>
// numeric_limits with exactly the members derivative()
// uses, none inherited from numeric_limits<double>.
template<>
struct std::numeric_limits<Real>
{
    using D = std::numeric_limits<double>;
    static constexpr bool is_specialized = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_bounded = true;
    static constexpr bool is_exact = false;
    static constexpr int radix = D::radix;   // 2
    static constexpr int digits = D::digits; // 53
    static constexpr Real epsilon()
    { return Real(D::epsilon()); }
    static constexpr Real max() { return Real(D::max()); }
    static constexpr Real denorm_min()
    { return Real(D::denorm_min()); }
};
// </usertype-limits>
// clang-format on

int main()
{
    // clang-format off
// <usertype-call>
static_assert(nxx::real<Real>);
auto const f = [](Real x) { return Real(std::exp(x.v)); };
auto const d = derivative(f, Real(1));
auto const r = derivative<Ridders>(f, Real(1));
if (d && r) {
    std::println("{} {}", d->value.v, r->value.v);
}
// </usertype-call>
    // clang-format on
    // Central<2> at the default step: about eps^(2/3) relative (DS7's median 5.8e-12 for e^x); Ridders about 1e-14.
    // The text claims both within 1e-9 of e.
    double const e = std::numbers::e;
    bool const ok = d && r && std::abs(d->value.v - e) <= 1e-9 * e && std::abs(r->value.v - e) <= 1e-9 * e;
    if (!ok) std::println("FAILED: a claim of the text does not hold");
    return ok ? 0 : 1;
}
