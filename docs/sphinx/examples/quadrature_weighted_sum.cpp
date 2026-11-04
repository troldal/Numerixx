// The data of the figure in docs/sphinx/quadrature/index.md ("A weighted sum, and an error estimate"): the Gauss rule
// of five points, GaussLegendre<5>, on 4 / (1 + x^2) over [0, 1], drawn as five bars side by side, the i-th of width
// w_i and height f(x_i), above the curve. The program prints the bars and the curve as CSV and compares them with the
// table in the raw string, which docs/sphinx/plots/quadrature_weighted_sum.py plots; main returns non-zero unless the
// two agree and what the page says of them holds (PHASE8_PLAN.md). The integrand and the rule use only + - * /, which
// IEEE 754 rounds correctly, so the table is the same on every target that does not contract a * b + c into a fused
// multiply-add (x86-64 without -mfma, WebAssembly).
//
// The rule does not report its nodes and weights, so the program asks it: the nodes are the points fn receives, and
// the weight of a node is the rule applied to the function that is 1 there and 0 at the other nodes.

#include <numerixx/quadrature.hpp>

#include "printed.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <print>
#include <string_view>
#include <vector>

namespace {

    namespace qd = nxx::quadrature;
    using Rule = qd::GaussLegendre<5>;

    double integrand(double x) { return 4.0 / (1.0 + x * x); }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    std::vector<double> nodes;
    auto const sum = qd::integrate<Rule>([&nodes](double x) {
        nodes.push_back(x);
        return integrand(x);
    }, 0.0, 1.0);
    if (!(sum && nodes.size() == 5)) return 2;
    std::ranges::sort(nodes);

    std::vector<double> weights;
    for (double const node : nodes) {
        auto const w = qd::integrate<Rule>([node](double x) { return x == node ? 1.0 : 0.0; }, 0.0, 1.0);
        if (!w) return 3;
        weights.push_back(*w);
    }

    // The bars side by side from 0, the i-th from the sum of the weights before it.
    std::println("node,weight,left,right,f");
    double left = 0.0;
    std::vector<double> lefts;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        lefts.push_back(left);
        std::println("{:.6g},{:.6g},{:.6g},{:.6g},{:.6g}", nodes[i], weights[i], left, left + weights[i],
                     integrand(nodes[i]));
        left += weights[i];
    }
    std::println("");
    std::println("x,f");
    for (int k = 0; k <= 50; ++k) std::println("{:.4g},{:.6g}", k / 50.0, integrand(k / 50.0));

    constexpr std::string_view weighted_sum_table = R"(
node,weight,left,right,f
0.0469101,0.118463,0,0.118463,3.99122
0.230765,0.239314,0.118463,0.357778,3.79776
0.5,0.284444,0.357778,0.642222,3.2
0.769235,0.239314,0.642222,0.881537,2.513
0.95309,0.118463,0.881537,1,2.09602

x,f
0,4
0.02,3.9984
0.04,3.99361
0.06,3.98565
0.08,3.97456
0.1,3.9604
0.12,3.94322
0.14,3.92311
0.16,3.90016
0.18,3.87447
0.2,3.84615
0.22,3.81534
0.24,3.78215
0.26,3.74672
0.28,3.7092
0.3,3.66972
0.32,3.62845
0.34,3.58551
0.36,3.54108
0.38,3.49528
0.4,3.44828
0.42,3.4002
0.44,3.35121
0.46,3.30142
0.48,3.25098
0.5,3.2
0.52,3.14861
0.54,3.09693
0.56,3.04507
0.58,2.99312
0.6,2.94118
0.62,2.88934
0.64,2.83768
0.66,2.78629
0.68,2.73523
0.7,2.68456
0.72,2.63435
0.74,2.58465
0.76,2.5355
0.78,2.48694
0.8,2.43902
0.82,2.39177
0.84,2.34522
0.86,2.29938
0.88,2.25428
0.9,2.20994
0.92,2.16638
0.94,2.12359
0.96,2.0816
0.98,2.0404
1,2
)";
    if (!(docs::printed() == docs::quoted(weighted_sum_table))) return 1;

    double const eps = std::numeric_limits<double>::epsilon();
    // What the page says. The bars' total area, the sum of w_i f(x_i), is the rule's result: the two sums differ only
    // in the order and grouping of their terms, a few rounding errors of a sum near 3.
    double area = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) area += weights[i] * integrand(nodes[i]);
    if (!(std::abs(area - *sum) <= 8.0 * eps * 4.0)) return 4;
    // The widths add up to the interval's length, 1: the rule is exact for a constant (to rounding).
    if (!(std::abs(left - 1.0) <= 8.0 * eps)) return 5;
    // Each bar contains its node, strictly inside: the Chebyshev-Markov-Stieltjes separation theorem (Szego,
    // Orthogonal Polynomials, Theorem 3.41.1), for the positive weights and the measure dx, whose distribution x is
    // strictly increasing; the margins are 0.047 or more.
    for (std::size_t i = 0; i < nodes.size(); ++i)
        if (!(lefts[i] < nodes[i] && nodes[i] < lefts[i] + weights[i])) return 6;
    // Five evaluations give pi to within 1.4e-8, to two digits: the rule's error is -1.37e-8 (numpy's leggauss(5),
    // and the closed-form nodes and weights in 50 digits, both independent of Numerixx), far above the rounding of the
    // sum, 1e-15.
    double const error = *sum - std::numbers::pi;
    if (!(error >= -1.4e-8 && error <= -1.35e-8)) return 7;
    return 0;
}
