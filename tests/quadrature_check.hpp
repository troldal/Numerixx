// What the nxx::quadrature tests share: the types, the references, the error measurement, an evaluation counter and a
// domain check, a real type without a default constructor, and the running error analysis of the Legendre recurrence
// with the error model of a computed Gauss node and weight (Phase 1's, which the Kronrod rule's builds on).
//
// References are closed forms at 70 digits (mpmath at 80), never Numerixx. Every tolerance comes from an error
// estimate stated where it is used: for GaussLegendre the truncation bound and the per-node rounding bound of
// docs/quadrature/PHASE1_API.md D, for Romberg the bound of Bauer, Rutishauser and Stiefel plus its rounding (C), for
// GaussKronrod the bound of docs/quadrature/design.md ("The tests' error model of the K15 rule").

#pragma once

#include <numerixx/quadrature.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ios>
#include <limits>
#include <type_traits>
#include <vector>

namespace quadrature_check {

    namespace mp = boost::multiprecision;
    namespace qd = nxx::quadrature;

    using bin50 = mp::cpp_bin_float_50;
    using bin50_et = mp::number<mp::cpp_bin_float<50>, mp::et_on>;
    using dec50 = mp::cpp_dec_float_50;
    // The type in which errors are measured: it holds every value of the types tested, and their differences, exactly
    // enough (100 digits against cpp_dec_float_50's at most about 60).
    using bin100 = mp::cpp_bin_float_100;

    inline constexpr double dmax = std::numeric_limits<double>::max();
    inline constexpr double dmin = std::numeric_limits<double>::denorm_min();
    inline constexpr double deps = std::numeric_limits<double>::epsilon();
    inline constexpr double dnan = std::numeric_limits<double>::quiet_NaN();
    inline constexpr double dinf = std::numeric_limits<double>::infinity();

    // 70-digit references (mpmath at 80 digits).
    inline constexpr char const *e_minus_1 = "1.718281828459045235360287471352662497757247093699959574966967627724077";
    inline constexpr char const *pi_over_4 = "0.7853981633974483096156608458198757210492923498437764552437361480769541";

    // A value of a tested type in bin100. cpp_dec_float_50 goes through its decimal string, exact at 70 significant
    // digits (more than it carries).
    template<typename T>
    bin100 wide(T const &v)
    {
        if constexpr (std::is_same_v<T, dec50>) return bin100(v.str(70, std::ios_base::scientific));
        else return bin100(v);
    }

    // |value - ref| in double, measured in bin100 (a difference formed in T would round).
    template<typename T>
    double error_of(T const &value, bin100 const &ref)
    { return static_cast<double>(bin100(abs(bin100(wide(value) - ref)))); }

    template<typename T>
    double error_of(T const &value, char const *ref)
    { return error_of(value, bin100(ref)); }

    // T's epsilon in double (cpp_bin_float_50's, 2^-167, is a double too).
    template<typename T>
    double eps_of()
    { return static_cast<double>(wide(T(std::numeric_limits<T>::epsilon()))); }

    // fn wrapped to count its evaluations.
    template<typename F>
    struct Counted
    {
        F fn;
        int *count;
        template<typename T>
        auto operator()(T const &x) const
        {
            ++*count;
            return fn(x);
        }
    };

    // fn wrapped to count its samples outside [lo, hi], the interval's ends in order, or outside (lo, hi) where open
    // (TanhSinh never samples an end).
    template<typename F>
    struct Confined
    {
        F fn;
        double lo;
        double hi;
        int *outside;
        bool open = false;
        double operator()(double x) const
        {
            if (x < lo || hi < x || (open && (x == lo || x == hi))) ++*outside;
            return fn(x);
        }
    };

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

    inline RunningBound running_bound(int n, double x, double u)
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

    // The a-priori bound of docs/quadrature/PHASE3_API.md B on the relative error of a tanh-sinh node's complement and
    // weight, in eps (exp's error taken as 1.5 ulp): for t >= 1/2, 2.5 + 4 sinh t + 3.5 e^-t and 6.2 + 4 sinh t +
    // 3.5 e^-t; for t < 1/2, where sinh and cosh come from their series of N additions, 2.5 + (N + 0.5) sinh t and
    // 6.2 + (N + 0.5) sinh t. The error of s = sinh t moves q = e^(-2s) by 2s times it, which dominates the outer
    // nodes.
    struct NodeBound
    {
        double comp;
        double w;
    };

    inline NodeBound node_bound(double t, int series_additions)
    {
        if (t >= 0.5) {
            double const common = 4.0 * std::sinh(t) + 3.5 * std::exp(-t);
            return { 2.5 + common, 6.2 + common };
        }
        double const series = (series_additions + 0.5) * std::sinh(t);
        return { 2.5 + series, 6.2 + series };
    }

    // Sums over the tanh-sinh double table, both sides, at the spacing h = 2^-k of a level k: sum h w (about 2, the
    // integral of dx/dt over the line), sum h w B_w and sum h w comp B_comp (B: node_bound, 8 series additions in
    // double); the largest of each over the levels 4 to 9, where a run can stop.
    struct TableSums
    {
        double weight;
        double w_bound;
        double comp_bound;
    };

    inline TableSums tanh_sinh_sums()
    {
        TableSums largest{ 0.0, 0.0, 0.0 };
        for (int k = 4; k <= 9; ++k) {
            double const h = 1.0 / static_cast<double>(1 << k);
            TableSums sums{ 0.0, 0.0, 0.0 };
            for (int level = 0; level <= k; ++level) {
                qd::detail::with_tanh_sinh_level<double>(level, [&sums, h, level](auto const &nodes) {
                    for (int i = 0; i < nodes.size; ++i) {
                        double const t = qd::detail::position(level, i) / 512.0;
                        NodeBound const b = node_bound(t, 8);
                        double const sides = t == 0.0 ? 1.0 : 2.0;
                        double const w = nodes.weight(i);
                        sums.weight += sides * h * w;
                        sums.w_bound += sides * h * w * b.w;
                        sums.comp_bound += sides * h * w * nodes.complement(i) * b.comp;
                    }
                });
            }
            largest = { std::max(largest.weight, sums.weight), std::max(largest.w_bound, sums.w_bound),
                        std::max(largest.comp_bound, sums.comp_bound) };
        }
        return largest;
    }

    constexpr auto square = [](double x) { return x * x; };
    constexpr auto cube = [](double x) { return x * x * x; };

    // A literal real type without a default constructor (the independent review's case): nxx::real does not require
    // one, so neither may the methods.
    struct Wrapped
    {
        double v;
        Wrapped() = delete;
        constexpr Wrapped(int i) : v(i) {}
        constexpr Wrapped(double d) : v(d) {}
        friend constexpr Wrapped operator+(Wrapped x, Wrapped y) { return Wrapped(x.v + y.v); }
        friend constexpr Wrapped operator-(Wrapped x, Wrapped y) { return Wrapped(x.v - y.v); }
        friend constexpr Wrapped operator*(Wrapped x, Wrapped y) { return Wrapped(x.v * y.v); }
        friend constexpr Wrapped operator/(Wrapped x, Wrapped y) { return Wrapped(x.v / y.v); }
        friend constexpr Wrapped operator-(Wrapped x) { return Wrapped(-x.v); }
        friend constexpr bool operator==(Wrapped, Wrapped) = default;
        friend constexpr auto operator<=>(Wrapped x, Wrapped y) { return x.v <=> y.v; }
    };

} // namespace quadrature_check

template<>
struct std::numeric_limits<quadrature_check::Wrapped> : std::numeric_limits<double>
{
    static constexpr quadrature_check::Wrapped epsilon() { return std::numeric_limits<double>::epsilon(); }
    static constexpr quadrature_check::Wrapped max() { return (std::numeric_limits<double>::max)(); }
    static constexpr quadrature_check::Wrapped denorm_min() { return std::numeric_limits<double>::denorm_min(); }
};

static_assert(nxx::real<quadrature_check::Wrapped>);
