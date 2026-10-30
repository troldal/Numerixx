// GaussLegendre<N> of nxx::quadrature (Phase 1, docs/quadrature/PHASE1_API.md D): the generator, the fixed rule's
// accuracy, its edges and its constant evaluation.
//
// The generator is checked against independent nodes and weights, never against Numerixx (quadrature_reference.hpp:
// Boost.Math's and mpmath's). Each node and weight must lie within an a-posteriori error bound: a first-order running
// error analysis of the recurrences at the computed node (below, and docs/quadrature/design.md). At every N from 1 to
// 64 the weights must sum to 2 and x^(2N-2) and x^(2N-1) on [0, 1] be integrated within the per-node rounding bound of
// D. Integrals are checked against closed forms within the truncation bound of D plus that rounding bound.

#include "quadrature_check.hpp"
#include "quadrature_reference.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace quadrature_check;
using quadrature_reference::Reference;
using quadrature_reference::reference_for;
using quadrature_reference::references;

namespace {

    // The computed rule's non-negative half in T, as integrate uses it: the compile-time table for float, double and
    // long double, the generator for every other T.
    template<typename T, int N>
    qd::detail::Rule<T, N> computed()
    {
        if constexpr (qd::detail::tabled<T>) return qd::detail::gauss_legendre_table<N, T>;
        else return qd::detail::gauss_legendre_rule<N, T>();
    }

    // A first-order bound on the rounding errors of the computed P_n(x) and P_n'(x) (detail::legendre), by a running
    // error analysis in the manner of Higham (Accuracy and Stability of Numerical Algorithms, 2002, ch. 3): the local
    // error of each step, from its own operations under the standard model fl(a op b) = (a op b)(1 + d), |d| <= u,
    // times its exact propagation to degree n through the linear recurrences. For a step k -> k + 1:
    //   p:  t1 = (2k + 1) x, t2 = t1 p_k, t3 = k p_(k-1), t4 = t2 - t3, p_(k+1) = t4 / (k + 1),
    //       local error <= u ((2 |t2| + |t3| + |t4|) / (k + 1) + |p_(k+1)|);
    //   p': s1 = x p'_k, s2 = p_k + s1, s3 = (2k + 1) s2, s4 = k p'_(k-1), s5 = s3 - s4, p'_(k+1) = s5 / (k + 1),
    //       local error <= u (((2k + 1) (|s1| + |s2|) + |s3| + |s4| + |s5|) / (k + 1) + |p'_(k+1)|).
    // An error d injected into p_j propagates by the homogeneous recurrence d_(k+1) = a_k d_k - b_k d_(k-1),
    // a_k = (2k + 1) x / (k + 1), b_k = k / (k + 1), and into p' through the term c_k d_k, c_k = (2k + 1) / (k + 1); an
    // error in p'_j by the same homogeneous recurrence. The responses at n are computed exactly (in double, enough for
    // a bound's first digits) by running the linearised recurrences from each j: O(n^2). This is a bound to first order
    // in u, not an estimate of typical errors: the local errors enter with their absolute values, the propagation with
    // its signs.
    struct RunningBound
    {
        double p;
        double dp;
        double error_p;
        double error_dp;
    };

    RunningBound running_bound(int n, double x, double u)
    {
        std::vector<double> p(static_cast<std::size_t>(n) + 1), d(p.size()), local_p(p.size()), local_d(p.size());
        p[0] = 1.0;
        p[1] = x;
        d[0] = 0.0;
        d[1] = 1.0;
        for (int k = 1; k < n; ++k) {
            std::size_t const i = static_cast<std::size_t>(k);
            double const kk = k;
            double const t2 = (2.0 * kk + 1.0) * x * p[i];
            double const t3 = kk * p[i - 1];
            double const t4 = t2 - t3;
            p[i + 1] = t4 / (kk + 1.0);
            local_p[i + 1] = u * ((2.0 * std::abs(t2) + std::abs(t3) + std::abs(t4)) / (kk + 1.0) + std::abs(p[i + 1]));
            double const s1 = x * d[i];
            double const s2 = p[i] + s1;
            double const s3 = (2.0 * kk + 1.0) * s2;
            double const s4 = kk * d[i - 1];
            double const s5 = s3 - s4;
            d[i + 1] = s5 / (kk + 1.0);
            local_d[i + 1] = u * (((2.0 * kk + 1.0) * (std::abs(s1) + std::abs(s2)) + std::abs(s3) + std::abs(s4)
                                   + std::abs(s5)) / (kk + 1.0)
                                  + std::abs(d[i + 1]));
        }
        double error_p = 0.0;
        double error_dp = 0.0;
        for (int j = 2; j <= n; ++j) {
            // Unit errors at p_j and at p'_j, carried to degree n.
            double dp0 = 0.0, dp1 = 1.0, dd0 = 0.0, dd1 = 0.0; // from p_j
            double ed0 = 0.0, ed1 = 1.0;                       // from p'_j
            for (int k = j; k < n; ++k) {
                double const kk = k;
                double const a = (2.0 * kk + 1.0) * x / (kk + 1.0);
                double const b = kk / (kk + 1.0);
                double const c = (2.0 * kk + 1.0) / (kk + 1.0);
                double const dp2 = a * dp1 - b * dp0;
                double const dd2 = a * dd1 - b * dd0 + c * dp1;
                double const ed2 = a * ed1 - b * ed0;
                dp0 = dp1;
                dp1 = dp2;
                dd0 = dd1;
                dd1 = dd2;
                ed0 = ed1;
                ed1 = ed2;
            }
            std::size_t const i = static_cast<std::size_t>(j);
            error_p += std::abs(dp1) * local_p[i];
            error_dp += std::abs(dd1) * local_p[i] + std::abs(ed1) * local_d[i];
        }
        std::size_t const last = static_cast<std::size_t>(n);
        return { p[last], d[last], error_p, error_dp };
    }

    // The error model of a computed node x and weight w of the N-point rule in T (docs/quadrature/design.md):
    // - the node: |x - x*| = |P_N(x)| / |P_N'| to first order, and the exact P_N(x) differs from the computed one,
    //   residual, by at most error_p: |dx| <= (|residual| + error_p) / |P_N'(x)|. The residual is the generator's own
    //   value at its result, computed in T.
    // - the weight 2 / ((1 - x)(1 + x) P_N'(x)^2), relative: the node's error times 2 |x| / (1 - x^2) (|d ln w / dx|
    //   at a zero, from the Legendre equation), 6u for the formula's six roundings, and twice the relative error of the
    //   computed P_N', 2 error_dp / |P_N'(x)|.
    struct Model
    {
        double x;
        double w;
    };

    template<typename T>
    Model model_of(int n, T const &node)
    {
        double const u = eps_of<T>() / 2.0;
        double const x = static_cast<double>(wide(node));
        RunningBound const r = running_bound(n, x, u);
        double const residual = static_cast<double>(wide(T(qd::detail::legendre(n, node).p)));
        double const dp = std::abs(r.dp);
        double const dx = (std::abs(residual) + r.error_p) / dp;
        double const w = 2.0 / ((1.0 - x) * (1.0 + x) * dp * dp);
        double const relative = 2.0 * std::abs(x) / ((1.0 - x) * (1.0 + x)) * dx + 6.0 * u + 2.0 * r.error_dp / dp;
        return { dx, w * relative };
    }

    // A node of the full rule on [-1, 1], with the errors of its node and weight.
    struct Node
    {
        double x;
        double w;
        double dx;
        double dw;
    };

    // The full N-point rule: the reference values and the computed rule's actual errors where a reference exists,
    // else the computed values and the model above.
    template<typename T, int N>
    std::vector<Node> nodes_of()
    {
        auto const rule = computed<T, N>();
        Reference const *ref = reference_for(N);
        std::vector<Node> full;
        for (std::size_t i = 0; i < rule.x.size(); ++i) {
            Node node{};
            if (ref != nullptr) {
                bin100 const x(ref->x[i]);
                bin100 const w(ref->w[i]);
                node = { static_cast<double>(x), static_cast<double>(w), error_of(rule.x[i], x),
                         error_of(rule.w[i], w) };
            }
            else {
                Model const m = model_of(N, rule.x[i]);
                node = { static_cast<double>(wide(rule.x[i])), static_cast<double>(wide(rule.w[i])), m.x, m.w };
            }
            full.push_back(node);
            if (node.x != 0.0) full.push_back({ -node.x, node.w, node.dx, node.dw });
        }
        return full;
    }

    // The per-node rounding bound of D for a rule on [a, b], in double:
    //   |hw| sum_i (|dw_i| |f_i| + w_i |f'_i| (|hw| |dx_i| + eps (|c| + 1.5 |hw x_i|))),
    // the errors of the weights and of the nodes, and the rounding of fn's argument (of c, hw, hw x_i and c + hw x_i,
    // u each); plus the rounding of each term t_i, three operations (w_i / 2 and two products), and of fn's value, rho
    // relative; plus the compensated sum's, u |I| + 2 N u^2 sum |t_i| (Neumaier): (2 eps + rho + 2 N eps^2) sum |t_i|.
    template<typename F, typename D>
    double rounding_bound(std::vector<Node> const &nodes, F f, D df, double a, double b, double eps, double rho)
    {
        double const hw = std::abs(b / 2 - a / 2);
        double const c = a / 2 + b / 2;
        double weights = 0.0;
        double terms = 0.0;
        for (Node const &n : nodes) {
            double const x = c + (b / 2 - a / 2) * n.x;
            double const fx = std::abs(f(x));
            // 1.5 eps first: 1.5 hw |x_i| can overflow on [-max, max].
            double const argument = hw * n.dx + eps * std::abs(c) + 1.5 * eps * (hw * std::abs(n.x));
            weights += n.dw * fx + n.w * std::abs(df(x)) * argument;
            terms += hw * n.w * fx;
        }
        double const count = static_cast<double>(nodes.size());
        return hw * weights + (2.0 * eps + rho + 2.0 * count * eps * eps) * terms;
    }

    // The truncation bound of D, (b - a)^(2N+1) (N!)^4 / ((2N + 1) ((2N)!)^3) max |f^(2N)|, through logarithms.
    double truncation(int n, double width, double log_max_derivative)
    {
        double const m = static_cast<double>(n);
        return std::exp((2.0 * m + 1.0) * std::log(width) + 4.0 * std::lgamma(m + 1.0)
                        - 3.0 * std::lgamma(2.0 * m + 1.0) - std::log(2.0 * m + 1.0) + log_max_derivative);
    }

    // x^p by p - 1 products, so fn rounds by at most about (p / 2) eps relative.
    template<typename T>
    T power(T const &x, int p)
    {
        T y = T(1);
        for (int k = 0; k < p; ++k) y = T(y * x);
        return y;
    }

    // Regression baselines, per referenced rule (CLAUDE.md's exception for a quality bound): the model above bounds
    // every node and weight, but for the weights it is loose (the running bound on P' sums its local errors in
    // absolute value: about 50 times the measured sum of |dw|), so these pin the generator's measured quality. The
    // largest |dx_i| / eps and the sum over the N nodes of |dw_i| / eps, each the maximum over float, double, long
    // double (80-bit and binary128; MSVC's is double) and cpp_bin_float_50, recorded on 2026-10-30 at the commit that
    // added this file ("Add nxx::quadrature Phase 1"), with a margin of 1.5. cpp_dec_float_50 errs below the
    // references' resolution.
    struct Baseline
    {
        int n;
        double max_dx;
        double sum_dw;
    };

    constexpr Baseline baselines[] = {
        { 1, 0.0, 0.0 },      { 2, 0.244, 0.0 },    { 3, 0.275, 2.444 },  { 5, 0.375, 3.678 },  { 7, 0.436, 3.933 },
        { 10, 0.316, 5.5 },   { 15, 0.404, 5.09 },  { 20, 0.319, 5.426 }, { 25, 0.331, 6.756 }, { 30, 0.28, 8.38 },
        { 40, 0.35, 8.632 },  { 50, 0.428, 7.785 }, { 64, 0.333, 10.372 },
    };

    constexpr Baseline baseline_for(int n)
    {
        for (Baseline const &b : baselines)
            if (b.n == n) return b;
        return { n, 0.0, 0.0 };
    }

    // The largest ratio of an actual error to its model, over the nodes and over the weights checked so far.
    struct Ratios
    {
        double x = 0.0;
        double w = 0.0;
    };

    template<typename T, int N>
    void check_generator(Ratios &ratios)
    {
        Reference const *ref = reference_for(N);
        REQUIRE(ref != nullptr);
        auto const rule = computed<T, N>();
        REQUIRE(ref->x.size() == rule.x.size());
        double const eps = eps_of<T>();
        double max_dx = 0.0;
        double sum_dw = 0.0;
        for (std::size_t i = 0; i < rule.x.size(); ++i) {
            INFO("N = " << N << ", node " << i);
            double const dx = error_of(rule.x[i], ref->x[i]);
            double const dw = error_of(rule.w[i], ref->w[i]);
            Model const m = model_of(N, rule.x[i]);
            CHECK(dx <= m.x);
            CHECK(dw <= m.w);
            ratios.x = std::max(ratios.x, dx / m.x);
            ratios.w = std::max(ratios.w, dw / m.w);
            max_dx = std::max(max_dx, dx / eps);
            sum_dw += ((N % 2 == 1 && i == 0) ? 1.0 : 2.0) * dw / eps;
        }
        INFO("N = " << N << ", the regression baseline");
        CHECK(max_dx <= 1.5 * baseline_for(N).max_dx);
        CHECK(sum_dw <= 1.5 * baseline_for(N).sum_dw);
        // The layout integrate relies on: ascending, 0 first for odd N.
        for (std::size_t i = 1; i < rule.x.size(); ++i) CHECK(rule.x[i - 1] < rule.x[i]);
        if constexpr (N % 2 == 1) CHECK(rule.x[0] == T(0));
        else CHECK(rule.x[0] > T(0));
    }

    // Every reference, for one type.
    template<typename T, std::size_t... I>
    Ratios check_generators(std::index_sequence<I...>)
    {
        Ratios ratios;
        (check_generator<T, references[I].n>(ratios), ...);
        return ratios;
    }

    // At one N: the weights sum to 2, and x^(2N-2) and x^(2N-1) on [0, 1], which the rule integrates exactly, are
    // within the rounding bound. The weight sum is formed in bin100 (exact for these terms), so |sum - 2| <= sum |dw_i|
    // with the model's dw_i: checked where no reference exists. At a referenced N, whose dw_i are the actual errors, it
    // would hold by the triangle inequality alone; check_generator covers that N.
    template<typename T, int N>
    void check_every_rule()
    {
        INFO("N = " << N);
        auto const rule = computed<T, N>();
        std::vector<Node> const nodes = nodes_of<T, N>();
        double const eps = eps_of<T>();
        bin100 sum = 0;
        for (std::size_t i = 0; i < rule.w.size(); ++i)
            sum += bin100((N % 2 == 1 && i == 0) ? 1 : 2) * wide(rule.w[i]);
        if (reference_for(N) == nullptr) {
            double dw = 0.0;
            for (Node const &n : nodes) dw += n.dw;
            CHECK(error_of(sum, bin100(2)) <= dw);
        }
        for (int const p : { 2 * N - 2, 2 * N - 1 }) {
            INFO("x^" << p);
            auto const r = qd::integrate<qd::GaussLegendre<N>>([p](T const &x) { return power(x, p); }, T(0), T(1));
            REQUIRE(r.has_value());
            double const bound = rounding_bound(
                nodes, [p](double x) { return std::pow(x, p); },
                [p](double x) { return p == 0 ? 0.0 : p * std::pow(x, p - 1); }, 0.0, 1.0, eps, p * eps / 2);
            CHECK(error_of(*r, bin100(1) / (p + 1)) <= bound);
        }
    }

    template<typename T, int... N>
    void check_every_rules(std::integer_sequence<int, N...>)
    { (check_every_rule<T, N + 1>(), ...); }

    // 1 / (1 + x^2) on [0, 1] = pi / 4: |f^(2N)| <= (2N)! on the real line (f = Im 1 / (x - i), and |x - i| >= 1), so
    // the truncation bound is (N!)^4 / ((2N + 1) ((2N)!)^2). |f'| = 2x / (1 + x^2)^2. fn rounds three times on values
    // whose relative errors add (x^2, 1 + x^2, the quotient): rho = 2 eps.
    template<typename T, int N>
    void check_runge()
    {
        INFO("N = " << N);
        auto const fn = [](T const &x) { return T(T(1) / T(T(1) + T(x * x))); };
        auto const r = qd::integrate<qd::GaussLegendre<N>>(fn, T(0), T(1));
        REQUIRE(r.has_value());
        double const eps = eps_of<T>();
        double const bound = truncation(N, 1.0, std::lgamma(2.0 * N + 1.0))
                             + rounding_bound(
                                 nodes_of<T, N>(), [](double x) { return 1.0 / (1.0 + x * x); },
                                 [](double x) { return 2.0 * x / ((1.0 + x * x) * (1.0 + x * x)); }, 0.0, 1.0, eps,
                                 2.0 * eps);
        CHECK(error_of(*r, pi_over_4) <= bound);
    }

    // The per-node rounding bound of x^p on [a, b] with the N-point rule in T, for the constant-evaluation cases.
    template<int N, typename T = double>
    double polynomial_bound(int p, double a, double b)
    {
        double const eps = eps_of<T>();
        return rounding_bound(
            nodes_of<T, N>(), [p](double x) { return std::pow(x, p); },
            [p](double x) { return p == 0 ? 0.0 : p * std::pow(x, p - 1); }, a, b, eps, p * eps / 2);
    }

    constexpr auto identity = [](double x) { return x; };
    constexpr auto power127 = [](double x) {
        double y = 1.0;
        for (int k = 0; k < 127; ++k) y *= x;
        return y;
    };

    // The constant-evaluation tolerances, each at most the rounding bound of its case, which
    // "quadrature: gauss: constant evaluation" checks at run time (the bounds are 15.7, 1.88, 58.7 and 2.24 eps), so a
    // static_assert that passes is within the bound. GL7's is the summary's (PHASE1_API.md).
    constexpr double tol_gl1 = 15.0 * deps;
    constexpr double tol_gl2 = 1.5 * deps;
    constexpr double tol_gl7 = 4.0 * 9.0 * deps;
    constexpr double tol_gl64 = 1.9 * deps;

    // Neumaier's sum of 1, 1e100, 1, -1e100 is 2, where a plain sum gives 0.
    constexpr double neumaier_example()
    {
        qd::detail::Sum<double> s;
        for (double const t : { 1.0, 1e100, 1.0, -1e100 }) s.add(t);
        return s.value();
    }

} // namespace

// Constant evaluation of the smallest rule, the first even one, the summary's case and the largest (the budget's
// limit), within their rounding bounds; .value(), whose throw on an empty result is a compile error (Q6).
static_assert(nxx::math::abs(qd::integrate<qd::GaussLegendre<1>>(identity, 0.0, 3.0).value() - 4.5) <= tol_gl1);
static_assert(nxx::math::abs(qd::integrate<qd::GaussLegendre<2>>(cube, 0.0, 1.0).value() - 0.25) <= tol_gl2);
static_assert(nxx::math::abs(qd::integrate<qd::GaussLegendre<7>>(square, 0.0, 3.0).value() - 9.0) <= tol_gl7);
static_assert(nxx::math::abs(qd::integrate<qd::GaussLegendre<64>>(power127, 0.0, 1.0).value() - 1.0 / 128) <= tol_gl64);
// float and long double, the other compile-time tables: the summary's case within 36 eps of the type, below its bound
// in each (checked at run time too).
static_assert(nxx::math::abs(qd::integrate<qd::GaussLegendre<7>>([](float x) { return x * x; }, 0.0f, 3.0f).value()
                             - 9.0f)
              <= 36.0f * std::numeric_limits<float>::epsilon());
static_assert(
    nxx::math::abs(qd::integrate<qd::GaussLegendre<7>>([](long double x) { return x * x; }, 0.0L, 3.0L).value() - 9.0L)
    <= 36.0L * std::numeric_limits<long double>::epsilon());
// The edges are values in constant expressions too.
static_assert(qd::integrate<qd::GaussLegendre<10>>(square, 1.0, 1.0) == 0.0); // a == b: 0, fn not called
static_assert(qd::integrate<qd::GaussLegendre<10>>(square, dnan, 1.0) == std::unexpected(qd::error::non_finite));
static_assert(qd::integrate<qd::GaussLegendre<10>>(square, 0.0, dinf) == std::unexpected(qd::error::non_finite));
static_assert(nxx::math::abs(qd::integrate<qd::GaussLegendre<7>>(square, 3.0, 0.0).value() + 9.0) <= tol_gl7);

// The internals the rule rests on, in constant expressions (run-time twins in "quadrature: gauss: internals").
static_assert(qd::detail::zeros_above(7, 0.0) == 3);  // three positive zeros; 0, itself a zero, is not counted
static_assert(qd::detail::zeros_above(7, -1.0) == 7); // every zero lies in (-1, 1)
static_assert(qd::detail::zeros_above(7, 1.0) == 0);
static_assert(qd::detail::zeros_above(8, 0.5) == 3); // P_8's zeros above 0.5: 0.525532, 0.796666, 0.960290
static_assert(neumaier_example() == 2.0);

// The overloads: the fixed rule takes no criteria, both ends must have one real type, N is 1 to 64, and fn returns T.
template<typename M, typename A, typename B>
concept integrable = requires(A a, B b) { qd::integrate<M>(square, a, b); };
template<typename M>
concept with_criteria = requires { qd::integrate<M>(square, 0.0, 1.0, { .rel_tol = 1e-12 }); };
static_assert(integrable<qd::GaussLegendre<1>, double, double> && integrable<qd::GaussLegendre<64>, double, double>);
static_assert(!integrable<qd::GaussLegendre<0>, double, double> && !integrable<qd::GaussLegendre<65>, double, double>);
static_assert(!integrable<qd::GaussLegendre<-3>, double, double>);
static_assert(!integrable<qd::GaussLegendre<10>, double, int>); // an integer end: both ends deduce T
static_assert(!integrable<int, double, double>);                // not a method
static_assert(!with_criteria<qd::GaussLegendre<10>> && with_criteria<qd::Romberg>);
template<typename F>
concept integrable_et = requires(F fn) { qd::integrate<qd::GaussLegendre<5>>(fn, bin50_et(0), bin50_et(1)); };
static_assert(!integrable_et<decltype([](bin50_et const &x) { return x * x; })>); // an expression template
static_assert(integrable_et<decltype([](bin50_et const &x) { return bin50_et(x * x); })>);

TEST_CASE("quadrature: gauss: nodes and weights against Boost.Math and mpmath")
{
    auto const all = std::make_index_sequence<std::size(references)>{};
    // The largest ratios of actual error to model over the 13 references, recorded in docs/quadrature/design.md: nodes
    // 0.27 (float), 0.26 (double), 0.24 (80-bit long double), 0.44 (cpp_bin_float_50), weights 0.15, 0.11, 0.13, 0.09.
    for (Ratios const r : { check_generators<float>(all), check_generators<double>(all),
                            check_generators<long double>(all), check_generators<bin50>(all),
                            check_generators<bin50_et>(all), check_generators<dec50>(all) }) {
        CHECK(r.x <= 1.0);
        CHECK(r.w <= 1.0);
    }
}

TEST_CASE("quadrature: gauss: every N from 1 to 64")
{
    check_every_rules<float>(std::make_integer_sequence<int, 64>{});
    check_every_rules<double>(std::make_integer_sequence<int, 64>{});
    check_every_rules<long double>(std::make_integer_sequence<int, 64>{});
    // The 50-digit types, whose rules the generator computes at each call, on a sample of N.
    check_every_rule<bin50, 7>();
    check_every_rule<bin50, 20>();
    check_every_rule<bin50_et, 20>();
    check_every_rule<dec50, 5>();
    check_every_rule<dec50, 20>();
}

TEST_CASE("quadrature: gauss: accuracy")
{
    // e^x on [0, 1] = e - 1, with std::exp (a budget of 2 ulps, rho = 2 eps); |f^(2N)| <= e. N = 1 is the midpoint
    // rule (error 0.070 against a bound of e / 24 = 0.113); from N = 7 on the truncation is below 1e-16.
    auto const exp_fn = [](double x) { return std::exp(x); };
    auto const exp_case = [&]<int N>(qd::GaussLegendre<N>) {
        INFO("N = " << N);
        auto const r = qd::integrate<qd::GaussLegendre<N>>(exp_fn, 0.0, 1.0);
        REQUIRE(r.has_value());
        double const bound = truncation(N, 1.0, 1.0)
                             + rounding_bound(nodes_of<double, N>(), exp_fn, exp_fn, 0.0, 1.0, deps, 2.0 * deps);
        CHECK(error_of(*r, e_minus_1) <= bound);
    };
    exp_case(qd::GaussLegendre<1>{});
    exp_case(qd::GaussLegendre<3>{});
    exp_case(qd::GaussLegendre<7>{});
    exp_case(qd::GaussLegendre<10>{});
    exp_case(qd::GaussLegendre<20>{});

    check_runge<double, 20>();
    check_runge<double, 64>();
    check_runge<float, 10>();
    check_runge<float, 64>();
    check_runge<long double, 30>();
    check_runge<long double, 64>();
    // The 50-digit types, whose rules the generator computes at each call: truncation 1e-78 at N = 64.
    check_runge<bin50, 64>();
    check_runge<bin50_et, 64>();
    check_runge<dec50, 64>();
    check_runge<bin50, 40>(); // truncation 4e-49: a case the truncation dominates
}

TEST_CASE("quadrature: gauss: constant evaluation")
{
    // The run-time twins of the static_asserts above, and their tolerances against the rounding bound.
    CHECK(tol_gl1 <= polynomial_bound<1>(1, 0.0, 3.0));
    CHECK(tol_gl2 <= polynomial_bound<2>(3, 0.0, 1.0));
    CHECK(tol_gl7 <= polynomial_bound<7>(2, 0.0, 3.0));
    CHECK(tol_gl64 <= polynomial_bound<64>(127, 0.0, 1.0));
    CHECK(36.0 * eps_of<float>() <= polynomial_bound<7, float>(2, 0.0, 3.0));
    CHECK(36.0 * eps_of<long double>() <= polynomial_bound<7, long double>(2, 0.0, 3.0));
    CHECK(std::abs(*qd::integrate<qd::GaussLegendre<1>>(identity, 0.0, 3.0) - 4.5) <= tol_gl1);
    CHECK(std::abs(*qd::integrate<qd::GaussLegendre<2>>(cube, 0.0, 1.0) - 0.25) <= tol_gl2);
    CHECK(std::abs(*qd::integrate<qd::GaussLegendre<7>>(square, 0.0, 3.0) - 9.0) <= tol_gl7);
    CHECK(std::abs(*qd::integrate<qd::GaussLegendre<64>>(power127, 0.0, 1.0) - 1.0 / 128) <= tol_gl64);
    CHECK(std::abs(*qd::integrate<qd::GaussLegendre<7>>(square, 3.0, 0.0) + 9.0) <= tol_gl7);
}

TEST_CASE("quadrature: gauss: internals")
{
    CHECK(qd::detail::zeros_above(7, 0.0) == 3);
    CHECK(qd::detail::zeros_above(7, -1.0) == 7);
    CHECK(qd::detail::zeros_above(7, 1.0) == 0);
    CHECK(qd::detail::zeros_above(8, 0.5) == 3);
    CHECK(neumaier_example() == 2.0);
    // A non-finite term, or an overflow of the running sum, makes the sum non-finite: the one check per sum.
    qd::detail::Sum<double> inf;
    for (double const t : { 1.0, dinf, 1.0 }) inf.add(t);
    CHECK(!nxx::math::isfinite(inf.value()));
    qd::detail::Sum<double> overflow;
    for (double const t : { dmax, dmax, -dmax }) overflow.add(t);
    CHECK(!nxx::math::isfinite(overflow.value()));
}

TEST_CASE("quadrature: gauss: edges")
{
    int n = 0;
    auto const counted_exp = Counted{ [](double x) { return std::exp(x); }, &n };

    // N evaluations of fn, and none for a == b, a NaN or an infinite end.
    REQUIRE(qd::integrate<qd::GaussLegendre<10>>(counted_exp, 0.0, 1.0).has_value());
    CHECK(n == 10);
    n = 0;
    REQUIRE(qd::integrate<qd::GaussLegendre<7>>(counted_exp, 0.0, 1.0).has_value());
    CHECK(n == 7);
    n = 0;
    auto const same = qd::integrate<qd::GaussLegendre<10>>(counted_exp, 2.0, 2.0);
    REQUIRE(same.has_value());
    CHECK(*same == 0.0); // defined as 0 (A3), not computed
    for (auto const &[a, b] : { std::pair{ dnan, 1.0 }, std::pair{ 0.0, dnan }, std::pair{ -dinf, 1.0 },
                                std::pair{ 0.0, dinf } })
        CHECK(qd::integrate<qd::GaussLegendre<10>>(counted_exp, a, b) == std::unexpected(qd::error::non_finite));
    CHECK(n == 0);

    // Reversed ends negate every term exactly (hw and b - a change sign, the nodes are the same), so the negated
    // result is within the forward case's bound.
    auto const exp_fn = [](double x) { return std::exp(x); };
    auto const reversed = qd::integrate<qd::GaussLegendre<10>>(exp_fn, 1.0, 0.0);
    REQUIRE(reversed.has_value());
    CHECK(error_of(-*reversed, e_minus_1)
          <= rounding_bound(nodes_of<double, 10>(), exp_fn, exp_fn, 0.0, 1.0, deps, 2.0 * deps));

    // fn not finite at a node: non_finite, from the one check of the sum. The centre of an odd rule; an outer node.
    CHECK(qd::integrate<qd::GaussLegendre<7>>([](double x) { return 1.0 / (x - 0.5); }, 0.0, 1.0)
          == std::unexpected(qd::error::non_finite));
    CHECK(qd::integrate<qd::GaussLegendre<10>>([](double x) { return x > 0.9 ? dnan : x; }, 0.0, 1.0)
          == std::unexpected(qd::error::non_finite));
    CHECK(qd::integrate<qd::GaussLegendre<10>>([](double x) { return x > 0.9 ? -dinf : x; }, 0.0, 1.0)
          == std::unexpected(qd::error::non_finite));

    // The order of each term's products (detail::Interval::term). fn is constant in every case: f' = 0, rho = 0.
    auto const zero = [](double) { return 0.0; };
    // fn = 1e308 on [0, 1e-10]: 1e298 (1e308 fl(1e-10)); each term is scaled before the sum, which does not overflow.
    auto const huge = [](double) { return 1e308; };
    auto const h = qd::integrate<qd::GaussLegendre<10>>(huge, 0.0, 1e-10);
    REQUIRE(h.has_value());
    CHECK(error_of(*h, "1.00000000000000004741126094494e298")
          <= rounding_bound(nodes_of<double, 10>(), huge, zero, 0.0, 1e-10, deps, 0.0));

    // fn = 1e-310 (subnormal) on [0, 1e10]: fl(1e-310) 1e10 = 1e-300, with the scale (b - a)(w_i / 2) formed first, so
    // that no intermediate is subnormal; fn's value first would lose about 1000 eps.
    auto const subnormal_fn = [](double) { return 1e-310; };
    auto const s = qd::integrate<qd::GaussLegendre<10>>(subnormal_fn, 0.0, 1e10);
    REQUIRE(s.has_value());
    CHECK(error_of(*s, "9.99999999999996944932750289769e-301")
          <= rounding_bound(nodes_of<double, 10>(), subnormal_fn, zero, 0.0, 1e10, deps, 0.0));

    // [-max, max], where b - a overflows and the terms are hw (w_i f_i): fn = 0 gives 0 (the bound is 0), and
    // fn = 1e-300 gives 2 max fl(1e-300) = 3.5953862697246314e8, for N = 1 too, whose weight 2 would overflow hw w_1.
    auto const tiny = [](double) { return 1e-300; };
    auto const z = qd::integrate<qd::GaussLegendre<10>>(zero, -dmax, dmax);
    REQUIRE(z.has_value());
    CHECK(error_of(*z, "0") <= rounding_bound(nodes_of<double, 10>(), zero, zero, -dmax, dmax, deps, 0.0));
    char const *wide_ref = "359538626.972463141629054847463";
    auto const t10 = qd::integrate<qd::GaussLegendre<10>>(tiny, -dmax, dmax);
    auto const t1 = qd::integrate<qd::GaussLegendre<1>>(tiny, -dmax, dmax);
    REQUIRE(t10.has_value());
    REQUIRE(t1.has_value());
    CHECK(error_of(*t10, wide_ref) <= rounding_bound(nodes_of<double, 10>(), tiny, zero, -dmax, dmax, deps, 0.0));
    CHECK(error_of(*t1, wide_ref) <= rounding_bound(nodes_of<double, 1>(), tiny, zero, -dmax, dmax, deps, 0.0));

    // [0, 3 denorm_min] with fn = 1e300: 3 denorm_min 1e300 = 1.48e-23, with fn's value multiplied first, where
    // (b - a) w_i / 2 would round in the subnormals (by up to 33 %). The bound is linear in the scale of a constant fn:
    // that of [0, 3] with fn = 1, times denorm_min 1e300 (formed in double, where the subnormal width would round).
    auto const big = [](double) { return 1e300; };
    auto const one = [](double) { return 1.0; };
    for (auto const &[g, bound] :
         { std::pair{ qd::integrate<qd::GaussLegendre<10>>(big, 0.0, 3.0 * dmin),
                      rounding_bound(nodes_of<double, 10>(), one, zero, 0.0, 3.0, deps, 0.0) },
           std::pair{ qd::integrate<qd::GaussLegendre<2>>(big, 0.0, 3.0 * dmin),
                      rounding_bound(nodes_of<double, 2>(), one, zero, 0.0, 3.0, deps, 0.0) } }) {
        REQUIRE(g.has_value());
        CHECK(error_of(*g, "1.48219693752373971035210123429e-23") <= bound * (1e300 * dmin));
    }

    // A double result for a float T narrows silently, as documented: 2x is exact in double and in float at the nodes.
    auto const narrowing = qd::integrate<qd::GaussLegendre<7>>([](float x) { return 2.0 * x; }, 0.0f, 1.0f);
    REQUIRE(narrowing.has_value());
    CHECK(error_of(*narrowing, "1") <= rounding_bound(nodes_of<float, 7>(), [](double x) { return 2.0 * x; },
                                                      [](double) { return 2.0; }, 0.0, 1.0, eps_of<float>(), 0.0));
}

// Emscripten's default build compiles throw but aborts on it rather than unwinding to a catch, so the test runs
// everywhere else.
#ifndef __EMSCRIPTEN__
TEST_CASE("quadrature: gauss: fn's exceptions pass through")
{
    auto const throwing = [](double x) -> double {
        if (x > 0.9) throw std::runtime_error("fn failed");
        return x;
    };
    CHECK_THROWS_AS((void)qd::integrate<qd::GaussLegendre<10>>(throwing, 0.0, 1.0), std::runtime_error);
}
#endif
