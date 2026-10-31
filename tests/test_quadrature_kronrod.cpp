// GaussKronrod<15> of nxx::quadrature (Phase 2, docs/quadrature/PHASE2_API.md), the default method; and the edge cases
// of both methods with an estimate, Romberg and GaussKronrod<15>, as one template.
//
// The rule is checked against Boost.Math's 115-digit G7-K15 table (quadrature_reference.hpp) within an a-posteriori
// bound, as Phase 1's Gauss-Legendre rules. Smooth integrals are checked within a derived bound: the truncation of K15
// on panels at most a tenth of the interval wide plus the rounding (below, kronrod_bound). Kahaner's 21 (Kahaner 1971,
// as Gander and Gautschi 2000 list them) at the default against references from closed forms or mpmath at 60 digits,
// with the honesty check |value - ref| <= abs_error; the non-smooth ones also within rel_tol |ref|. Evaluation counts
// are asserted exactly in double for fn of + - * /, sqrt and floor (exact), against an independent Python
// implementation of the specification that mirrors the header's operation order, with the C++ generator's double table
// loaded into it (docs/quadrature/design.md records its counts; the prototype itself is not kept). No count is
// asserted where a libm function's last bits could move a decision.

#include "quadrature_check.hpp"
#include "quadrature_reference.hpp"

#include <doctest/doctest.h>

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace quadrature_check;

namespace {

    using gk15 = qd::GaussKronrod<15>;

    // The pair as integrate uses it in T: the compile-time table for float, double and long double, the generator
    // for every other T.
    template<typename T>
    qd::detail::Kronrod<T> kronrod_of()
    {
        if constexpr (qd::detail::tabled<T>) return qd::detail::gauss_kronrod_table<T>;
        else return qd::detail::gauss_kronrod_rule<T>();
    }

    // The exact Stieltjes coefficients a_k (quadrature_reference.hpp) in bin100.
    bin100 exact_coefficient(std::size_t k)
    {
        return bin100(quadrature_reference::stieltjes_numerator[k])
               / bin100(quadrature_reference::stieltjes_denominator[k]);
    }

    // A first-order bound on the rounding errors of the computed E_8(x) and E_8'(x) (detail::legendre_series): the
    // coefficients' own errors, |a_k - fl(a_k)| against the exact rationals, times |P_k| and |P_k'|; the errors of the
    // computed P_k and P_k' (running_bound, Phase 1's running error analysis) times |a_k|; and the series' own
    // operations, a product and an addition per degree, u (|a_k P_k| + |partial sum|).
    struct SeriesBound
    {
        double e;
        double de;
        double d2e;
        double error_e;
        double error_de;
    };

    template<typename T>
    SeriesBound series_bound(std::array<T, 9> const &computed, double x, double u)
    {
        SeriesBound r{ 0, 0, 0, 0, 0 };
        double partial = 0.0;
        double partial_d = 0.0;
        for (std::size_t k = 0; k <= 8; ++k) {
            int const n = static_cast<int>(k);
            RunningBound const b = n == 0 ? RunningBound{ 1.0, 0.0, 0.0, 0.0 } : running_bound(n, x, u);
            double const a = static_cast<double>(exact_coefficient(k));
            double const da = static_cast<double>(bin100(abs(bin100(wide(computed[k]) - exact_coefficient(k)))));
            // P_k'' from the Legendre equation, (1 - x^2) P'' = 2x P' - k (k + 1) P.
            double const d2 = (2.0 * x * b.dp - n * (n + 1.0) * b.p) / ((1.0 - x) * (1.0 + x));
            partial += a * b.p;
            partial_d += a * b.dp;
            r.e += a * b.p;
            r.de += a * b.dp;
            r.d2e += a * d2;
            r.error_e += da * std::abs(b.p) + std::abs(a) * b.error_p + u * (std::abs(a * b.p) + std::abs(partial));
            r.error_de +=
                da * std::abs(b.dp) + std::abs(a) * b.error_dp + u * (std::abs(a * b.dp) + std::abs(partial_d));
        }
        return r;
    }

    // The bounds on a computed node and its Kronrod weight (docs/quadrature/design.md, "The tests' error model"):
    // - a Kronrod node x: |dx| <= (|E_8(x) computed| + error_e) / |E_8'(x)|, the residual the generator's own;
    //   its weight c / (P_7 E_8'), c = 1/4: the node's error times |w'| (w' / w = -(P_7' / P_7 + E_8'' / E_8')), plus
    //   the relative errors of the computed P_7 and E_8' and 3u for the product and the quotient;
    // - a Gauss node: Phase 1's model (model_of) for the node and the Gauss weight; the Kronrod weight w_G + h,
    //   h = c / (P_7' E_8), adds |h| times the relative errors of P_7' and E_8 and 3u, the node's error times |h'|
    //   (h' / h = -(P_7'' / P_7' + E_8' / E_8)), and u |w_K| for the addition.
    template<typename T>
    Model kronrod_model(qd::detail::Kronrod<T> const &rule, std::size_t i)
    {
        double const u = eps_of<T>() / 2.0;
        double const x = static_cast<double>(wide(rule.x[i]));
        std::array<T, 9> const coefficients = qd::detail::stieltjes<T>();
        SeriesBound const e = series_bound(coefficients, x, u);
        RunningBound const p = running_bound(7, x, u);
        double const c = 0.25;
        if (i % 2 == 1) {
            double const residual =
                static_cast<double>(wide(T(qd::detail::legendre_series(coefficients, rule.x[i]).p)));
            double const dx = (std::abs(residual) + e.error_e) / std::abs(e.de);
            double const w = c / (p.p * e.de);
            double const slope = std::abs(w) * std::abs(p.dp / p.p + e.d2e / e.de);
            double const relative = p.error_p / std::abs(p.p) + e.error_de / std::abs(e.de) + 3 * u;
            double const weight = slope * dx + std::abs(w) * relative;
            return { dx, weight };
        }
        Model const g = model_of(7, rule.x[i]);
        double const d2p = (2.0 * x * p.dp - 56.0 * p.p) / ((1.0 - x) * (1.0 + x));
        double const h = c / (p.dp * e.e);
        double const slope = std::abs(h) * std::abs(d2p / p.dp + e.de / e.e);
        double const weight = g.w + slope * g.x
                              + std::abs(h) * (p.error_dp / std::abs(p.dp) + e.error_e / std::abs(e.e) + 3 * u)
                              + u * std::abs(static_cast<double>(wide(rule.wk[i])));
        return { g.x, weight };
    }

    // The computed pair's actual errors against Boost's table, in units of T's eps: the sum of |dw_K| over the 15
    // points and the largest |dx|, for the bounds below.
    struct RuleErrors
    {
        double sum_dw = 0.0;
        double max_dx = 0.0;
    };

    template<typename T>
    RuleErrors measure_kronrod()
    {
        qd::detail::Kronrod<T> const rule = kronrod_of<T>();
        double const eps = eps_of<T>();
        RuleErrors r;
        for (std::size_t i = 0; i < qd::detail::Kronrod<T>::size; ++i) {
            r.sum_dw += (i == 0 ? 1.0 : 2.0) * error_of(rule.wk[i], quadrature_reference::boost_k15_w[i]) / eps;
            r.max_dx = std::max(r.max_dx, error_of(rule.x[i], quadrature_reference::boost_k15_x[i]) / eps);
        }
        return r;
    }

    // Every node and weight of the pair within its bound (the largest ratios are in docs/quadrature/design.md), the
    // Gauss part GaussLegendre<7>'s with Boost's G7 weights, the layout integrate relies on.
    template<typename T>
    void check_kronrod()
    {
        qd::detail::Kronrod<T> const rule = kronrod_of<T>();
        for (std::size_t i = 0; i < qd::detail::Kronrod<T>::size; ++i) {
            INFO("node " << i);
            Model const m = kronrod_model(rule, i);
            CHECK(error_of(rule.x[i], quadrature_reference::boost_k15_x[i]) <= m.x);
            CHECK(error_of(rule.wk[i], quadrature_reference::boost_k15_w[i]) <= m.w);
            if (i % 2 == 0) {
                CHECK(rule.x[i] == qd::detail::gauss_legendre_rule<7, T>().x[i / 2]);
                CHECK(error_of(rule.wg[i], quadrature_reference::boost7_w[i / 2]) <= model_of(7, rule.x[i]).w);
            }
            else CHECK(rule.wg[i] == T(0));
        }
        for (std::size_t i = 1; i < qd::detail::Kronrod<T>::size; ++i) CHECK(rule.x[i - 1] < rule.x[i]);
    }

    // A bound, required finite: a bound that overflowed would pass every check.
    double finite(double bound)
    {
        REQUIRE(std::isfinite(bound));
        return bound;
    }

    // The bound on the adaptive K15 result on [a, b] (docs/quadrature/design.md):
    // - truncation: every panel lies inside a mesh panel, at most (b - a) / 10 wide. On a panel of width L, K15, exact
    //   to degree 23 with positive weights summing to L, errs by at most 2 L E_23(f) (I and K15 both bounded by L times
    //   the sup of f - p), and E_23(f) <= (L / 2)^24 max |f^(24)| / (2^23 24!) (Chebyshev interpolation); summed over
    //   the panels, 2 (b - a) ((b - a) / 20)^24 max |f^(24)| / (2^23 24!);
    // - the weights' errors: sum over panels of h_p sum_i |dw_i| |f_i| <= sum |dw| max |f| (b - a) / 2;
    // - the nodes' errors and the rounding of fn's argument: w_i |f'| (h_p |dx_i| + eps (|c_p| + 1.5 h_p |x_i|)),
    //   summed: max |f'| (b - a) ((b - a) / 20 max |dx| + eps (max(|a|, |b|) + 1.5 (b - a) / 20));
    // - the terms, three operations each, fn's own rounding rho, the panels' compensated sums and the running sum:
    //   (4 eps + rho) L1 (the running sum's O(n eps^2) terms over at most the cap's 64 b updates stay below eps L1).
    struct Smooth
    {
        double log_max_24; // log max |f^(24)| on [a, b], or -inf for a polynomial of degree at most 23
        double max_f;
        double max_df;
        double rho;
        double l1;
    };

    // 1 / (1 + x^2) on [0, 1] = pi / 4: |f^(24)| <= 24! (f = Im 1 / (x - i), |x - i| >= 1), |f'| <= 0.65, three
    // roundings whose relative errors add (rho = 2 eps), L1 = pi / 4.
    Smooth runge_smooth(double eps)
    { return Smooth{ std::lgamma(25.0), 1.0, 0.65, 2.0 * eps, 0.7854 }; }

    double kronrod_bound(Smooth const &s, double a, double b, double eps, RuleErrors const &rule)
    {
        double const width = std::abs(b - a);
        double const truncation = std::exp(std::log(2.0 * width) + 24.0 * std::log(width / 20.0) + s.log_max_24
                                           - 23.0 * std::log(2.0) - std::lgamma(25.0));
        double const weights = rule.sum_dw * eps * s.max_f * width / 2.0;
        double const nodes = s.max_df * width
                             * (width / 20.0 * rule.max_dx * eps
                                + eps * (std::max(std::abs(a), std::abs(b)) + 1.5 * width / 20.0));
        return finite(truncation + weights + nodes + (4.0 * eps + s.rho) * s.l1);
    }

    // A run with its evaluations counted.
    template<typename M = gk15, typename T, typename F>
    std::pair<std::expected<qd::Integral<T>, qd::Failure<T>>, int> counted(F fn, T a, T b, qd::Criteria<T> c = {})
    {
        int n = 0;
        auto const r = qd::integrate<M>(Counted<F>{ fn, &n }, a, b, c);
        return { r, n };
    }

    // A hash of x's bits in [-1, 1]: noise that never resolves, however fine the panels.
    double hash(double x)
    {
        std::uint64_t h = std::bit_cast<std::uint64_t>(x) * 0x9E3779B97F4A7C15ull;
        h ^= h >> 29;
        h *= 0xBF58476D1CE4E5B9ull;
        h ^= h >> 32;
        return static_cast<double>(h >> 11) / 4503599627370496.0 - 1.0;
    }

    double root(double x)
    { return std::sqrt(x); }

    constexpr auto runge = [](double x) { return 1.0 / (1.0 + x * x); };

} // namespace

// The type is constrained on N, and the default method is GaussKronrod<15>.
template<int N>
concept kronrod_exists = requires { typename qd::GaussKronrod<N>; };
static_assert(kronrod_exists<15> && !kronrod_exists<21> && !kronrod_exists<7>);

// Constant evaluation: polynomial and rational cases within their bounds (the run-time twins in "quadrature: kronrod:
// constant evaluation" check the tolerances against kronrod_bound, 137.7 and 6.85 eps), and a NaN fn non_finite there.
// x^2 on [0, 3]: both rules exact, every panel at its floor at once (150 evaluations); 1 / (1 + x^2) on [0, 1] (150),
// against the double nearest pi / 4 (its half ulp inside the margin).
static_assert(nxx::math::abs(qd::integrate(square, 0.0, 3.0).value().value - 9.0) <= 14.0 * 9.0 * deps);
static_assert(nxx::math::abs(qd::integrate(runge, 0.0, 1.0).value().value - 0.78539816339744831) <= 6.0 * deps);
static_assert(qd::integrate([](double) { return std::numeric_limits<double>::quiet_NaN(); }, 0.0, 1.0)
              == std::unexpected(qd::Failure<double>{ qd::error::non_finite, std::nullopt }));
// float and long double, the other compile-time tables, within the same bound in their eps.
namespace {
    constexpr auto square_f = [](float x) { return x * x; };
    constexpr auto square_l = [](long double x) { return x * x; };
} // namespace
static_assert(nxx::math::abs(qd::integrate(square_f, 0.0f, 3.0f).value().value - 9.0f)
              <= 14.0f * 9.0f * std::numeric_limits<float>::epsilon());
static_assert(nxx::math::abs(qd::integrate(square_l, 0.0L, 3.0L).value().value - 9.0L)
              <= 14.0L * 9.0L * std::numeric_limits<long double>::epsilon());
// The bisection loop in constant evaluation: a kink at 0.37 (7 bisections, 360 evaluations), its value within its own
// estimate of the reference (0.37^2 + 0.63^2) / 2 = 0.2669 (0.37 as fn's double differs by 2e-17, far below the
// estimate, 3e-9); "quadrature: kronrod: constant evaluation" checks that run time gives the same result.
namespace {
    constexpr auto kink_37 = [](double x) { return x < 0.37 ? 0.37 - x : x - 0.37; };
    constexpr auto kink_result = qd::integrate(kink_37, 0.0, 1.0);
} // namespace
static_assert(kink_result && nxx::math::abs(kink_result->value - 0.2669) <= kink_result->abs_error);

TEST_CASE("quadrature: kronrod: the rule against Boost.Math")
{
    check_kronrod<float>();
    check_kronrod<double>();
    check_kronrod<long double>();
    check_kronrod<bin50>();
    check_kronrod<bin50_et>();
    check_kronrod<dec50>();
    // The weights are positive (the truncation bound relies on it) and sum to 2.
    qd::detail::Kronrod<double> const rule = qd::detail::gauss_kronrod_table<double>;
    bin100 sum = 0;
    for (std::size_t i = 0; i < rule.wk.size(); ++i) {
        CHECK(rule.wk[i] > 0.0);
        sum += bin100(i == 0 ? 1 : 2) * bin100(rule.wk[i]);
    }
    CHECK(error_of(sum, bin100(2)) <= measure_kronrod<double>().sum_dw * deps);
}

TEST_CASE("quadrature: kronrod: accuracy on smooth cases")
{
    // Within kronrod_bound, and, separately, the honesty check.
    auto const check = []<typename T>(T, auto fn, double a, double b, Smooth const &s, bin100 const &ref,
                                      qd::Criteria<T> const &c) {
        auto const r = qd::integrate(fn, T(a), T(b), c);
        REQUIRE(r.has_value());
        double const error = error_of(r->value, ref);
        CHECK(error <= kronrod_bound(s, a, b, eps_of<T>(), measure_kronrod<T>()));
        CHECK(error <= static_cast<double>(wide(r->abs_error))); // the honesty check
    };
    // e^x on [0, 1]: |f^(24)| = |f| = |f'| <= e; std::exp within 2 ulps (rho = 2 eps).
    Smooth const exp_case{ 1.0, 2.7183, 2.7183, 2.0 * deps, 1.7183 };
    auto const exp_fn = [](double x) { return std::exp(x); };
    for (double const rel : { 1.5e-8, 1e-14, 0.0 })
        check(0.0, exp_fn, 0.0, 1.0, exp_case, bin100(e_minus_1), { .rel_tol = rel });
    for (double const rel : { 1.5e-8, 1e-14, 0.0 })
        check(0.0, runge, 0.0, 1.0, runge_smooth(deps), bin100(pi_over_4), { .rel_tol = rel });
    // The 50-digit types, where the generator runs once and is cached: 1 / (1 + x^2) at the default, where the
    // truncation term, 1.4e-38, makes the bound loose (the error measured in cpp_bin_float_50 is 5.5e-45, 2.6e6 times
    // below it: the panels the run ends with are narrower than the mesh's tenth), and x^20 on [0, 1] = 1/21, exact for
    // K15 (degree <= 23), so the rounding alone: |f| <= 1, |f'| <= 20, 19 products (rho = 10 eps), L1 = 1/21.
    auto const mp_case = []<typename T>(T) {
        double const eps = eps_of<T>();
        auto const fn = [](T const &x) { return T(T(1) / T(T(1) + T(x * x))); };
        auto const r = qd::integrate(fn, T(0), T(1));
        REQUIRE(r.has_value());
        CHECK(error_of(r->value, pi_over_4) <= kronrod_bound(runge_smooth(eps), 0.0, 1.0, eps, measure_kronrod<T>()));
        CHECK(error_of(r->value, pi_over_4) <= static_cast<double>(wide(r->abs_error))); // the honesty check
        auto const power = [](T const &x) {
            T y = T(1);
            for (int k = 0; k < 20; ++k) y = T(y * x);
            return y;
        };
        auto const p = qd::integrate(power, T(0), T(1));
        REQUIRE(p.has_value());
        Smooth const polynomial{ -std::numeric_limits<double>::infinity(), 1.0, 20.0, 10.0 * eps, 1.0 / 21 };
        CHECK(error_of(p->value, bin100(1) / 21) <= kronrod_bound(polynomial, 0.0, 1.0, eps, measure_kronrod<T>()));
        CHECK(error_of(p->value, bin100(1) / 21) <= static_cast<double>(wide(p->abs_error))); // the honesty check
    };
    mp_case(bin50{});
    mp_case(bin50_et{});
    mp_case(dec50{});
}

TEST_CASE("quadrature: kronrod: Kahaner's 21")
{
    // At the default criteria, every one ends ok with its estimate covering its error (the honesty check). The
    // non-smooth ones (a step, sqrt x, x^1.5, 1/sqrt x, log x) also within rel_tol |ref|. References: closed forms, or
    // mpmath at 60 digits split at the features (8, 12, 13, 17, 18); the constants as the doubles fn uses (0.3, pi).
    double const pi = 3.141592653589793;
    auto const sech = [](double u) { return std::abs(u) > 710.0 ? 0.0 : 1.0 / std::cosh(u); };
    struct Case
    {
        char const *name;
        std::function<double(double)> fn;
        double a;
        double b;
        char const *ref;
        bool smooth;
    };
    Case const cases[] = {
        { "K01 e^x", [](double x) { return std::exp(x); }, 0.0, 1.0, e_minus_1, true },
        { "K02 step at 0.3", [](double x) { return x >= 0.3 ? 1.0 : 0.0; }, 0.0, 1.0,
          "0.700000000000000011102230246251565404236316681", false },
        { "K03 sqrt x", root, 0.0, 1.0, "0.666666666666666666666666666666666666666666667", false },
        { "K04 23/25 cosh x - cos x", [](double x) { return 23.0 / 25 * std::cosh(x) - std::cos(x); }, -1.0, 1.0,
          "0.479428226688801667358577961835307500641394964", true },
        { "K05 1/(x^4 + x^2 + 0.9)", [](double x) { return 1.0 / (x * x * x * x + x * x + 0.9); }, -1.0, 1.0,
          "1.58223296372967290254422882158697662278341583", true },
        { "K06 x^1.5", [](double x) { return x * std::sqrt(x); }, 0.0, 1.0, "0.4", false },
        { "K07 1/sqrt x", [](double x) { return 1.0 / std::sqrt(x); }, 0.0, 1.0, "2", false },
        { "K08 1/(1 + x^4)", [](double x) { return 1.0 / (1.0 + x * x * x * x); }, 0.0, 1.0,
          "0.866972987339911037573995163882870713652175367", true },
        { "K09 2/(2 + sin 10 pi x)", [pi](double x) { return 2.0 / (2.0 + std::sin(10.0 * pi * x)); }, 0.0, 1.0,
          "1.15470053837925153504879037285971226033032707", true },
        { "K10 1/(1 + x)", [](double x) { return 1.0 / (1.0 + x); }, 0.0, 1.0,
          "0.693147180559945309417232121458176568075500134", true },
        { "K11 1/(1 + e^x)", [](double x) { return 1.0 / (1.0 + std::exp(x)); }, 0.0, 1.0,
          "0.379885493041722475368236626490320926160220049", true },
        { "K12 x/(e^x - 1)", [](double x) { return x != 0.0 ? x / std::expm1(x) : 1.0; }, 0.0, 1.0,
          "0.777504634112248276417586545425710507192477296", true },
        { "K13 sin(100 pi x)/(pi x)", [pi](double x) { return std::sin(100.0 * pi * x) / (pi * x); }, 0.1, 1.0,
          "0.00909863753916684327023835633605198863744240186", true },
        { "K14 sqrt(50) e^(-50 pi x^2)", [pi](double x) { return std::sqrt(50.0) * std::exp(-50.0 * pi * x * x); },
          0.0, 10.0, "0.500000000000000009745429581298438901331102966", true },
        { "K15 25 e^(-25 x)", [](double x) { return 25.0 * std::exp(-25.0 * x); }, 0.0, 10.0, "1", true },
        { "K16 50/(pi (2500 x^2 + 1))", [pi](double x) { return 50.0 / (pi * (2500.0 * x * x + 1.0)); }, 0.0, 10.0,
          "0.499363381076456764102291181350590803222697435", true },
        { "K17 50 sinc^2(50 pi x)",
          [pi](double x) {
              double const s = std::sin(50.0 * pi * x) / (50.0 * pi * x);
              return 50.0 * s * s;
          },
          0.01, 1.0, "0.11213930374163741832311858996970147464243012", true },
        { "K18 cos(cos x + 3 sin x + ...)",
          [](double x) {
              return std::cos(std::cos(x) + 3.0 * std::sin(x) + 2.0 * std::cos(2.0 * x) + 3.0 * std::sin(2.0 * x)
                              + 3.0 * std::cos(3.0 * x));
          },
          0.0, pi, "0.83867634269442966550584383486134618746223133", true },
        { "K19 log x", [](double x) { return std::log(x); }, 0.0, 1.0, "-1", false },
        { "K20 1/(1.005 + x^2)", [](double x) { return 1.0 / (1.005 + x * x); }, -1.0, 1.0,
          "1.56439644406904990893786151733458442945248431", true },
        { "K21 sum sech(20^i (x - i/5))",
          [sech](double x) {
              double s = 0.0;
              double scale = 1.0;
              for (int i = 1; i <= 3; ++i) {
                  scale *= 20.0;
                  s += sech(scale * (x - 2.0 * i / 10));
              }
              return s;
          },
          0.0, 1.0, "0.163494943018637226588196396112688150878648606", true },
    };
    double const rel_tol = qd::Criteria<double>{}.rel_tol;
    for (Case const &c : cases) {
        INFO(c.name);
        auto const r = qd::integrate(c.fn, c.a, c.b);
        REQUIRE(r.has_value());
        double const error = error_of(r->value, c.ref);
        CHECK(error <= r->abs_error); // the honesty check
        if (!c.smooth) CHECK(error <= rel_tol * std::abs(static_cast<double>(bin100(c.ref))));
    }
}

TEST_CASE("quadrature: kronrod: evaluation counts and exits in double")
{
    // The prototype's counts with the C++ table (docs/quadrature/design.md, counts3.py), for fn of + - * /, sqrt and
    // floor, with the status and whether a failure carries a best estimate. The prototype names the exits the C++
    // does not report: the mesh, a bisection, the 2F test, the cap, stagnation, the heap empty (final panels), the
    // taint, the bounded descent.
    double const nan = std::numeric_limits<double>::quiet_NaN();
    enum class status { ok, not_converged, non_finite };
    struct Case
    {
        char const *name;
        std::pair<std::expected<qd::Integral<double>, qd::Failure<double>>, int> run;
        status expected;
        int evaluations;
    };
    auto const kink = [](double x) { return x < 0.37 ? 0.37 - x : x - 0.37; };
    auto const power_75 = [](double x) { return 1.0 / (root(x) * root(root(x))); };
    auto const power_875 = [](double x) { return 1.0 / (root(x) * root(root(x)) * root(root(root(x)))); };
    // x^-(127/128): 1 / (x^(1/2) x^(1/4) ... x^(1/128)), 7 nested square roots; beyond the x^-0.65 limit, its end
    // panels become final with errors above the tolerance.
    auto const power_127 = [](double x) {
        double r = x;
        double product = 1.0;
        for (int k = 0; k < 7; ++k) {
            r = root(r);
            product *= r;
        }
        return 1.0 / product;
    };
    // A NaN micro-region that a descent of 2 levels does not avoid, found by the prototype's search: at 1e-10 the stop
    // test passes, and only the taint makes the run not_converged (without it, ok).
    auto const holed = [nan](double x) {
        return std::abs(x - 0.48753003597833966) < 7.857561861803385e-05 ? nan : 1.0 / (1e-4 + (x - 0.42) * (x - 0.42));
    };
    auto const comb = [](double x) { return std::abs(x * 127.0 - std::floor(x * 127.0) - 0.5); };
    Case const cases[] = {
        { "1/(1+x^2) [0,1] (the mesh alone)", counted(runge, 0.0, 1.0), status::ok, 150 },
        { "1/(1+x^2) [0,1] 1e-14", counted(runge, 0.0, 1.0, { .rel_tol = 1e-14 }), status::ok, 150 },
        { "|x - 0.37| [0,1]", counted(kink, 0.0, 1.0), status::ok, 360 },
        { "|x - 0.37| [0,1] rel 0 (the 2F test)", counted(kink, 0.0, 1.0, { .rel_tol = 0.0 }), status::ok, 660 },
        // D14: the centre of the mesh panel [-0.2, 0.2] of [-1, 3] is exactly 0: one descent split.
        { "1/sqrt|x| [-1,3]", counted([](double x) { return 1.0 / root(std::abs(x)); }, -1.0, 3.0), status::ok, 2550 },
        { "sqrt(x - 0.5) [0,1] (the bounded descent)", counted([](double x) { return root(x - 0.5); }, 0.0, 1.0),
          status::non_finite, 105 },
        { "sqrt x [0,1]", counted(root, 0.0, 1.0), status::ok, 360 },
        { "x^-0.75 [0,1]", counted(power_75, 0.0, 1.0), status::ok, 2790 },
        { "x^-0.875 [0,1]", counted(power_875, 0.0, 1.0), status::ok, 5520 },
        { "x^-0.875 [0,1] 1e-12", counted(power_875, 0.0, 1.0, { .rel_tol = 1e-12 }), status::ok, 10710 },
        { "x^-(127/128) [0,1] (the heap empty)", counted(power_127, 0.0, 1.0), status::not_converged, 92130 },
        { "1/x^2 [0,1] (stagnation)", counted([](double x) { return 1.0 / (x * x); }, 0.0, 1.0), status::not_converged,
          1140 },
        { "1/(x - 0.5)^2 [0,1] (stagnation)",
          counted([](double x) { return 1.0 / ((x - 0.5) * (x - 0.5)); }, 0.0, 1.0), status::not_converged, 1140 },
        { "an undone NaN region, 1e-10 (the taint)", counted(holed, 0.0, 1.0, { .rel_tol = 1e-10 }),
          status::not_converged, 630 },
        // The cap: 3,392 bisections in double, 150 + 30 3,392 evaluations.
        { "a comb of 254 kinks, rel 0 (the cap)", counted(comb, 0.0, 1.0, { .rel_tol = 0.0 }), status::not_converged,
          101910 },
        // The check per panel beyond L1's: max on |x - 2.5| < 0.25, 0 elsewhere, on [0, 50]: the mesh panel [0, 5]
        // samples it at its centre only, where the 15-point term is finite (L1 = 0.52 max) and the 7-point one
        // overflows.
        { "max on |x - 2.5| < 0.25, [0,50] (G overflows)",
          counted([](double x) { return std::abs(x - 2.5) < 0.25 ? dmax : 0.0; }, 0.0, 50.0), status::non_finite, 15 },
    };
    for (Case const &c : cases) {
        INFO(c.name);
        auto const &r = c.run.first;
        status const got = r ? status::ok : (r.error().reason == qd::error::not_converged ? status::not_converged
                                                                                          : status::non_finite);
        CHECK(got == c.expected);
        CHECK(c.run.second == c.evaluations);
        // A failure carries its best estimate exactly when it is not_converged.
        if (!r) CHECK(r.error().best.has_value() == (got == status::not_converged));
    }
    // The taint's case passed the stop test: its estimate is within the request.
    auto const &taint = cases[13].run.first;
    REQUIRE(!taint.has_value());
    CHECK(taint.error().best->abs_error <= 1e-10 * std::abs(taint.error().best->value));
    // The cap's best estimate is honest (the comb's integral is exactly 1/4).
    auto const &cap = cases[14].run.first;
    REQUIRE(!cap.has_value());
    CHECK(error_of(cap.error().best->value, "0.25") <= cap.error().best->abs_error);
}

TEST_CASE("quadrature: kronrod: termination")
{
    // Stagnation on noise that never resolves: e^x (1 + 1e-10 h(x)), h a hash of x's bits. The best estimate is within
    // its error estimate plus the noise's amplitude, 1e-10 (e - 1), of e - 1.
    auto const noisy = qd::integrate([](double x) { return std::exp(x) * (1.0 + 1e-10 * hash(x)); }, 0.0, 1.0,
                                     { .rel_tol = 1e-12 });
    REQUIRE(!noisy.has_value());
    CHECK(noisy.error().reason == qd::error::not_converged);
    REQUIRE(noisy.error().best.has_value());
    CHECK(error_of(noisy.error().best->value, e_minus_1) <= noisy.error().best->abs_error + 1e-10 * 1.7183);

    // D14 on [-1, 3]: log |x| = 3 ln 3 - 4, the singularity at the centre node of a mesh panel; ok, honest and within
    // rel_tol.
    auto const logarithm = qd::integrate([](double x) { return std::log(std::abs(x)); }, -1.0, 3.0);
    REQUIRE(logarithm.has_value());
    char const *log_ref = "-0.704163133995670925814264289232422886057528327";
    CHECK(error_of(logarithm->value, log_ref) <= logarithm->abs_error);
    CHECK(error_of(logarithm->value, log_ref) <= 1.5e-8 * 0.7042);

    // The bounded descent in cpp_bin_float_50, where an unbounded one would walk about 2^31 levels: non_finite within
    // the 105 evaluations of double.
    auto const [deep, evaluations] =
        counted([](bin50 const &x) { return bin50(sqrt(bin50(x - bin50(0.5)))); }, bin50(0), bin50(1));
    CHECK(deep == std::unexpected(qd::Failure<bin50>{ qd::error::non_finite, std::nullopt }));
    CHECK(evaluations <= 105);

    // The cap's justification in cpp_bin_float_50 (64 b = 10,752): x^-0.875 on [0, 1] at its default (1e-25) ends ok,
    // after 2,820 bisections (84,750 evaluations); x^-0.9 by pow takes 3,550 (PHASE2_API.md E: 3,472 in the
    // prototype).
    auto const [strong, strong_evaluations] = counted(
        [](bin50 const &x) { return bin50(1 / bin50(sqrt(x) * sqrt(sqrt(x)) * sqrt(sqrt(sqrt(x))))); }, bin50(0),
        bin50(1));
    REQUIRE(strong.has_value());
    CHECK(error_of(strong->value, "8") <= strong->abs_error);
    CHECK(strong_evaluations < 150 + 30 * 64 * 168);
}

TEST_CASE("quadrature: kronrod: the substitution for a strong endpoint singularity")
{
    // x^-0.9 cos x on [0, 1] = sum_n (-1)^n / ((2n)! (2n + 0.1)) (mpmath). The singularity is stronger than about
    // x^-0.65, where |K - G| understates the end panel's error; with x = u^k, k (1 + alpha) >= 1, here k = 10,
    // dx = 10 u^9 du and x^-0.9 = u^-9, so the integrand becomes 10 cos(u^10): smooth, ok, honest and within rel_tol.
    char const *ref = "9.77184271163928998089477500829800531756618798";
    auto const substituted = qd::integrate(
        [](double u) {
            double u10 = u;
            for (int k = 1; k < 10; ++k) u10 *= u;
            return 10.0 * std::cos(u10);
        },
        0.0, 1.0);
    REQUIRE(substituted.has_value());
    CHECK(error_of(substituted->value, ref) <= substituted->abs_error);
    CHECK(error_of(substituted->value, ref) <= 1.5e-8 * 9.772);
}

TEST_CASE("quadrature: kronrod: constant evaluation")
{
    // The run-time twins of the static_asserts, and their tolerances against kronrod_bound.
    RuleErrors const rule = measure_kronrod<double>();
    Smooth const square_case{ -std::numeric_limits<double>::infinity(), 9.0, 6.0, deps / 2, 9.0 };
    CHECK(14.0 * 9.0 * deps <= kronrod_bound(square_case, 0.0, 3.0, deps, rule));
    CHECK(6.0 * deps <= kronrod_bound(runge_smooth(deps), 0.0, 1.0, deps, rule));
    // The bisection loop: the same result at run time as in constant evaluation, bit for bit.
    CHECK(qd::integrate(kink_37, 0.0, 1.0) == kink_result);
    // float and long double: their tolerances against the bound, and the run-time values.
    CHECK(14.0 * 9.0 * eps_of<long double>()
          <= kronrod_bound(Smooth{ -std::numeric_limits<double>::infinity(), 9.0, 6.0, eps_of<long double>() / 2, 9.0 },
                           0.0, 3.0, eps_of<long double>(), measure_kronrod<long double>()));
    CHECK(std::abs(qd::integrate(square_f, 0.0f, 3.0f)->value - 9.0f)
          <= 14.0f * 9.0f * std::numeric_limits<float>::epsilon());
    RuleErrors const rule_f = measure_kronrod<float>();
    Smooth const square_float{ -std::numeric_limits<double>::infinity(), 9.0, 6.0, eps_of<float>() / 2, 9.0 };
    CHECK(14.0 * 9.0 * eps_of<float>() <= kronrod_bound(square_float, 0.0, 3.0, eps_of<float>(), rule_f));
}

// The edge cases of both methods with an estimate, Romberg and GaussKronrod<15> (Phase 1's, moved here).
namespace {

    // The number of panels a method's floor sums over on a mesh where every panel is final at its floor.
    template<typename M>
    constexpr double floors = std::same_as<M, qd::Romberg> ? 1.0 : 10.0;

    // The rounding allowance of a result, by method: Romberg's 14 eps L1 and fn's 1.97 rho L1 (Phase 1);
    // Gauss-Kronrod's kronrod_bound without truncation (the cases are constants and e^x, whose truncation is
    // negligible: below).
    template<typename M>
    double allowance(double rho, double l1, double max_f, double max_df, double a, double b)
    {
        if constexpr (std::same_as<M, qd::Romberg>) return finite((14.0 * deps + 1.97 * rho) * l1);
        else {
            // On [-max, max] the width overflows double: its terms are the half-widths', doubled.
            double const half = std::abs(b / 2 - a / 2);
            RuleErrors const rule = measure_kronrod<double>();
            double const weights = rule.sum_dw * deps * max_f * half;
            // A constant has no node term (and 0 times the overflowing width would be NaN).
            double const argument = half / 10.0 * rule.max_dx * deps
                                    + deps * (std::max(std::abs(a), std::abs(b)) + 1.5 * half / 10.0);
            double const nodes = max_df == 0.0 ? 0.0 : max_df * 2.0 * half * argument;
            return finite(weights + nodes + (4.0 * deps + rho) * l1);
        }
    }

    constexpr auto zero_fn = [](double) { return 0.0; };
    constexpr auto nan_fn = [](double) { return dnan; };
    constexpr auto inf_fn = [](double x) { return x > 0.45 && x < 0.55 ? dinf : x; };
    constexpr auto max_fn = [](double) { return dmax; };
    constexpr auto large_fn = [](double) { return 0.8 * dmax; };
    constexpr auto wrapped_one = [](Wrapped const &) { return Wrapped(1); };

    // The constant-evaluation edges of the independent review of Phase 1, for one method: a zero fn on [-max, max]
    // (b - a not formed) is 0 with its floor; NaN and infinite values and a sum that overflows are non_finite; a sum
    // above max / 2 is non_finite there, conservatively, where run time gives 0.8 max; a type without a default
    // constructor integrates.
    template<typename M>
    constexpr bool constant_edges()
    {
        qd::Failure<double> const invalid{ qd::error::non_finite, std::nullopt };
        auto const zero = qd::integrate<M>(zero_fn, -dmax, dmax);
        auto const one = qd::integrate<M>(wrapped_one, Wrapped(0), Wrapped(1));
        return zero && zero->value == 0.0 && zero->abs_error == floors<M> * 50.0 * dmin
               && qd::integrate<M>(nan_fn, 0.0, 1.0) == std::unexpected(invalid)
               && qd::integrate<M>(inf_fn, 0.0, 1.0) == std::unexpected(invalid)
               && qd::integrate<M>(max_fn, 0.0, 3.0) == std::unexpected(invalid)
               && qd::integrate<M>(large_fn, 0.0, 1.0) == std::unexpected(invalid)
               && qd::integrate<M>(square, 2.0, 2.0) == qd::Integral<double>{ 0.0, 0.0 } && one
               && nxx::math::abs(one->value.v - 1.0) <= 4.0 * deps;
    }

} // namespace

static_assert(constant_edges<qd::Romberg>());
static_assert(constant_edges<gk15>());

TEST_CASE_TEMPLATE("quadrature: estimating methods: edges", M, qd::Romberg, gk15)
{
    qd::Failure<double> const invalid{ qd::error::non_finite, std::nullopt };

    // a == b: {0, 0}, fn not called; NaN and infinite ends and invalid tolerances: non_finite, fn not called.
    int n = 0;
    auto const c = Counted{ square, &n };
    auto const same = qd::integrate<M>(c, 2.0, 2.0);
    REQUIRE(same.has_value());
    CHECK(same->value == 0.0); // defined as 0 (A3), not computed
    CHECK(same->abs_error == 0.0);
    for (auto const &[a, b] : { std::pair{ dnan, 1.0 }, std::pair{ 0.0, dnan }, std::pair{ -dinf, 1.0 },
                                std::pair{ 0.0, dinf } })
        CHECK(qd::integrate<M>(c, a, b) == std::unexpected(invalid));
    for (qd::Criteria<double> const criteria :
         { qd::Criteria<double>{ .abs_tol = -1e-12 }, qd::Criteria<double>{ .rel_tol = -1e-12 },
           qd::Criteria<double>{ .abs_tol = dnan }, qd::Criteria<double>{ .rel_tol = dnan },
           qd::Criteria<double>{ .abs_tol = dinf }, qd::Criteria<double>{ .rel_tol = dinf } })
        CHECK(qd::integrate<M>(c, 0.0, 1.0, criteria) == std::unexpected(invalid));
    CHECK(n == 0);

    // The run-time twins of constant_edges, where run time gives the conservative case's value: 0.8 max, a constant.
    CHECK(qd::integrate<M>(zero_fn, -dmax, dmax) == qd::Integral<double>{ 0.0, floors<M> * 50.0 * dmin });
    CHECK(qd::integrate<M>(nan_fn, 0.0, 1.0) == std::unexpected(invalid));
    CHECK(qd::integrate<M>(inf_fn, 0.0, 1.0) == std::unexpected(invalid));
    // A sum that overflows: max on [0, 3].
    CHECK(qd::integrate<M>(max_fn, 0.0, 3.0) == std::unexpected(invalid));
    auto const large = qd::integrate<M>(large_fn, 0.0, 1.0);
    REQUIRE(large.has_value());
    CHECK(error_of(large->value, bin100(0.8 * dmax)) <= allowance<M>(0.0, 0.8 * dmax, 0.8 * dmax, 0.0, 0.0, 1.0));
    auto const one = qd::integrate<M>(wrapped_one, Wrapped(0), Wrapped(1));
    REQUIRE(one.has_value());
    CHECK(std::abs(one->value.v - 1.0) <= allowance<M>(0.0, 1.0, 1.0, 0.0, 0.0, 1.0));

    // Orientation: reversed ends give the negated integral, within the forward case's allowance. e^x on [0, 1] and
    // [1, 0]: |f| = |f'| <= e, std::exp within 2 ulps; the truncation is below 1e-60 for both (Phase 1's bound for
    // Romberg at its level; kronrod_bound's for the panels).
    auto const exp_fn = [](double x) { return std::exp(x); };
    auto const forward = qd::integrate<M>(exp_fn, 0.0, 1.0);
    auto const backward = qd::integrate<M>(exp_fn, 1.0, 0.0);
    REQUIRE(forward.has_value());
    REQUIRE(backward.has_value());
    CHECK(error_of(-backward->value, e_minus_1) <= backward->abs_error); // the honesty check
    CHECK(error_of(-backward->value, bin100(forward->value))
          <= 2.0 * allowance<M>(2.0 * deps, 1.7183, 2.7183, 2.7183, 0.0, 1.0));
    // Gauss-Kronrod integrates over [min, max] and negates: exactly.
    if constexpr (std::same_as<M, gk15>) {
        CHECK(backward->value == -forward->value);
        CHECK(backward->abs_error == forward->abs_error);
    }

    // The scales of fn: 1e-300 e^x and 1e300 e^x on [0, 1], honest (the relative errors as at scale 1).
    for (double const scale : { 1e-300, 1e300 }) {
        INFO("scale " << scale);
        auto const r = qd::integrate<M>([scale](double x) { return scale * std::exp(x); }, 0.0, 1.0);
        REQUIRE(r.has_value());
        CHECK(error_of(r->value, bin100(scale) * bin100(e_minus_1)) <= r->abs_error + 3.0 * deps * scale * 1.7183);
    }

    // An odd fn: x^3 on [-1, 1] = 0, its abs_error the floor (abs_tol 0): 50 eps L1, L1 between the exact 1/2
    // (Gauss-Kronrod: |x|^3 is a polynomial on each panel, the mesh having a point at 0) and Romberg's trapezoid sum
    // 0.5078125, each formed with a few roundings.
    auto const odd = qd::integrate<M>(cube, -1.0, 1.0, { .abs_tol = 0.0, .rel_tol = 1e-14 });
    REQUIRE(odd.has_value());
    CHECK(error_of(odd->value, "0") <= allowance<M>(deps, 0.5, 1.0, 3.0, -1.0, 1.0));
    CHECK(odd->abs_error >= 50.0 * deps * 0.5 * (1.0 - 1e-12));
    CHECK(odd->abs_error <= 50.0 * deps * 0.5078125 * (1.0 + 1e-12));

    // The order of each term's products (detail::Interval::term), on constant fn: the rounding alone (rho = 0).
    // fn = 1e308 on [0, 1e-10]: 1e298 (1e308 fl(1e-10)), no overflow, each term scaled before the sum.
    auto const huge = qd::integrate<M>([](double) { return 1e308; }, 0.0, 1e-10);
    REQUIRE(huge.has_value());
    CHECK(error_of(huge->value, "1.00000000000000004741126094494e298")
          <= allowance<M>(0.0, 1e298, 1e308, 0.0, 0.0, 1e-10));
    // fn = 1e-310 on [0, 1e10]: fl(1e-310) 1e10, the scale formed first, so nothing is subnormal.
    auto const tiny = qd::integrate<M>([](double) { return 1e-310; }, 0.0, 1e10);
    REQUIRE(tiny.has_value());
    CHECK(error_of(tiny->value, "9.99999999999996944932750289769e-301")
          <= allowance<M>(0.0, 1e-300, 1e-310, 0.0, 0.0, 1e10));
    // [-max, max], b - a overflows: fn = 1e-300 gives 2 max fl(1e-300); the constant denorm_min gives
    // 2 max denorm_min = 1.78e-15, with (hw y) w formed for |y| <= 1.
    auto const wide_tiny = qd::integrate<M>([](double) { return 1e-300; }, -dmax, dmax);
    REQUIRE(wide_tiny.has_value());
    CHECK(error_of(wide_tiny->value, "359538626.972463141629054847463")
          <= allowance<M>(0.0, 3.6e8, 1e-300, 0.0, -dmax, dmax));
    auto const smallest = qd::integrate<M>([](double) { return dmin; }, -dmax, dmax);
    REQUIRE(smallest.has_value());
    CHECK(error_of(smallest->value, bin100(2) * bin100(dmax) * bin100(dmin))
          <= allowance<M>(0.0, 1.78e-15, dmin, 0.0, -dmax, dmax));
    // [0, 3 denorm_min] with fn = 1e300: 1.48e-23, fn's value multiplied first.
    auto const subnormal = qd::integrate<M>([](double) { return 1e300; }, 0.0, 3.0 * dmin);
    REQUIRE(subnormal.has_value());
    CHECK(error_of(subnormal->value, "1.48219693752373971035210123429e-23")
          <= allowance<M>(0.0, 1.5e-23, 1e300, 0.0, 0.0, 3.0 * dmin));

    // The domain: every sample inside [min(a, b), max(a, b)] on subnormal widths, where c + hw u would leave the
    // interval without the clamp (on [D, 3D], D = denorm_min, c = hw = 2D), and on [-max, max]; both orientations.
    double const D = dmin;
    for (auto const &[a, b] : { std::pair{ D, 3.0 * D }, std::pair{ 3.0 * D, D }, std::pair{ 0.0, 3.0 * D },
                                std::pair{ D, 2.0 * D }, std::pair{ -3.0 * D, D }, std::pair{ -dmax, dmax },
                                std::pair{ dmax, -dmax }, std::pair{ 0.0, dmax } }) {
        INFO("[" << a << ", " << b << "]");
        int outside = 0;
        auto const fn = Confined{ [](double x) { return x / 2; }, std::min(a, b), std::max(a, b), &outside };
        (void)qd::integrate<M>(fn, a, b);
        (void)qd::integrate<M>(fn, a, b, { .rel_tol = 0.0 });
        CHECK(outside == 0);
    }
    // Phase 1's review's counterexample: 0 on [D, 3D], 1e300 outside it; the integral is 0, with a floor.
    auto const r = qd::integrate<M>([D](double x) { return x < D || x > 3.0 * D ? 1e300 : 0.0; }, D, 3.0 * D);
    REQUIRE(r.has_value());
    CHECK(error_of(r->value, "0") <= 0.0);
    CHECK(r->abs_error >= 50.0 * dmin);
    CHECK(r->abs_error <= 500.0 * dmin);
}

// Emscripten's default build compiles throw but aborts on it rather than unwinding to a catch, so the test runs
// everywhere else.
#ifndef __EMSCRIPTEN__
TEST_CASE_TEMPLATE("quadrature: estimating methods: fn's exceptions pass through", M, qd::Romberg, gk15)
{
    auto const throwing = [](double x) -> double {
        if (x > 0.6) throw std::runtime_error("fn failed");
        return x;
    };
    CHECK_THROWS_AS((void)qd::integrate<M>(throwing, 0.0, 1.0), std::runtime_error);
}
#endif
