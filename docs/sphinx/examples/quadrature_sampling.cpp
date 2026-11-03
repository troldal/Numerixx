// The data of the figure in docs/sphinx/quadrature/choose-a-method.md ("Measure the cost on your own function"): where
// the default and TanhSinh evaluate 1 / sqrt x on [0, 1], counted per decade of x. The program prints the counts as CSV
// and compares them with the table in the raw string, which docs/sphinx/plots/quadrature_sampling.py plots; main
// returns non-zero unless the two agree and what the page says of them holds. 1 / sqrt x and the methods use only
// + - * / and std::sqrt, which IEEE 754 rounds correctly, so the points are the same on every target that does not
// contract a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly); the decade of x is found by
// multiplying by 10, not by a libm's log10, for the same reason.

#include <numerixx/quadrature.hpp>

#include "printed.hpp"

#include <cmath>
#include <cstddef>
#include <map>
#include <print>
#include <string_view>
#include <utility>
#include <vector>

namespace {

    namespace qd = nxx::quadrature;

    // The decade d of x in (0, 1), 10^d <= x < 10^(d + 1) up to the rounding of the products: minus the number of times
    // x is multiplied by 10 to reach 1. Only for x > 0, which main checks first: 0 would never reach 1.
    int decade_of(double x)
    {
        int d = -1;
        for (double y = x * 10.0; y < 1.0; y *= 10.0) --d;
        return d;
    }

    // Every x the method M evaluates fn at, in order; empty if the integral fails.
    template<typename M>
    std::vector<double> points_of()
    {
        std::vector<double> xs;
        auto const fn = [&xs](double x) {
            xs.push_back(x);
            return 1.0 / std::sqrt(x);
        };
        if (!qd::integrate<M>(fn, 0.0, 1.0)) xs.clear();
        return xs;
    }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    auto const kronrod = points_of<qd::GaussKronrod<15>>();
    auto const tanh_sinh = points_of<qd::TanhSinh>();
    // Every point inside (0, 1), before decade_of sees it.
    for (auto const *points : {&kronrod, &tanh_sinh})
        for (double const x : *points)
            if (!(x > 0.0 && x < 1.0)) return 3;
    // decade -> evaluations of the default and of TanhSinh
    std::map<int, std::pair<int, int>> counts;
    for (double const x : kronrod) ++counts[decade_of(x)].first;
    for (double const x : tanh_sinh) ++counts[decade_of(x)].second;
    std::println("decade,default,tanh_sinh");
    for (auto const &[d, n] : counts) std::println("{},{},{}", d, n.first, n.second);

    constexpr std::string_view sampling_table = R"(
decade,default,tanh_sinh
-176,0,1
-65,0,1
-40,0,1
-35,0,1
-33,0,1
-31,0,1
-29,0,1
-27,0,1
-26,0,1
-24,0,1
-23,0,1
-21,0,1
-20,0,1
-19,0,1
-18,0,1
-17,0,1
-16,2,1
-15,7,1
-14,34,1
-13,95,1
-12,99,2
-11,104,1
-10,97,2
-9,98,2
-8,104,2
-7,98,2
-6,97,3
-5,103,4
-4,98,4
-3,92,6
-2,87,10
-1,135,73
)";
    if (!(docs::printed() == docs::quoted(sampling_table))) return 1;

    // What the page says: 1350 and 131 evaluations.
    if (!(kronrod.size() == 1350 && tanh_sinh.size() == 131)) return 2;
    // The default: its initial mesh of ten panels (150 evaluations), then 40 cuts of the panel at 0, 30 evaluations
    // each, the k-th inside [0, 0.1 / 2^k], so about 3.3 cuts in a decade (2^3.3 is about 10); 135 evaluations in
    // [0.1, 1), the nine panels it never cuts; about 100 (87 to 104) in every decade from 1e-13 to 0.1; none below
    // 1e-16.
    for (std::size_t k = 0; k < 40; ++k)
        for (std::size_t i = 150 + 30 * k; i < 150 + 30 * (k + 1); ++i)
            if (!(kronrod[i] <= 0.1 / static_cast<double>(1ULL << k))) return 11;
    if (!(counts[-1].first == 135)) return 12;
    for (int d = -13; d <= -2; ++d)
        if (!(counts[d].first >= 87 && counts[d].first <= 104)) return 4;
    if (!(counts[-16].first > 0)) return 5;
    for (auto const &[d, n] : counts)
        if (d < -16 && n.first != 0) return 6;
    // TanhSinh: 73 of its 131 points in [0.1, 1), one or two in every decade from 1e-20 to 1e-6, then 12 points
    // spread over the decades below, the smallest below 1e-175.
    if (!(counts[-1].second == 73)) return 7;
    for (int d = -20; d <= -7; ++d)
        if (!(counts[d].second >= 1 && counts[d].second <= 2)) return 8;
    int below_20 = 0;
    for (auto const &[d, n] : counts)
        if (d < -20) below_20 += n.second;
    if (!(below_20 == 12)) return 9;
    double smallest = 1.0;
    for (double const x : tanh_sinh) smallest = x < smallest ? x : smallest;
    if (!(smallest < 1e-175)) return 10;
    return 0;
}
