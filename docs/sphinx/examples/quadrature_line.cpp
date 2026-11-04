// The data of the figure in docs/sphinx/quadrature/infinite-intervals.md ("Where the maps look"): the points TanhSinh
// places on the whole line in its first five levels, which every run takes, here its run on 1 / (1 + x^2), and the peak
// 1 / (1 + (x - x0)^2) of the page's example at x0 = 0, 10 and 40. The program prints both as CSV and compares them
// with the table in the raw string, which docs/sphinx/plots/quadrature_line.py plots; main returns non-zero unless the
// two agree and what the page says of them holds (PHASE8_PLAN.md). The integrand, the peaks and the method use only
// + - * /, which IEEE 754 rounds correctly, so the table is the same on every target that does not contract a * b + c
// into a fused multiply-add (x86-64 without -mfma, WebAssembly).

#include <numerixx/quadrature.hpp>

#include "printed.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <iterator>
#include <limits>
#include <print>
#include <string_view>
#include <vector>

namespace {

    namespace qd = nxx::quadrature;

    double peak(double x, double x0) { return 1.0 / (1.0 + (x - x0) * (x - x0)); }

    // The points within the peak's half-width, |x - x0| <= 1, where it is at least half its height.
    std::ptrdiff_t within(std::vector<double> const &points, double x0)
    {
        return std::ranges::count_if(points, [x0](double x) { return std::abs(x - x0) <= 1.0; });
    }

    // The distance between the two points on either side of x0.
    double gap_at(std::vector<double> const &points, double x0)
    {
        auto const above = std::ranges::upper_bound(points, x0);
        return *above - *(above - 1);
    }

    // The points within the figure's range, |x| <= 60.
    std::vector<double> in_range(std::vector<double> const &points)
    {
        std::vector<double> kept;
        std::ranges::copy_if(points, std::back_inserter(kept), [](double x) { return std::abs(x) <= 60.0; });
        return kept;
    }

    // A run of TanhSinh on the peak at x0: the points fn has received after each level, sorted, and the result.
    struct Levels {
        std::vector<std::vector<double>> after;
        std::expected<qd::Integral<double>, qd::Failure<double>> result;
    };

    Levels levels_of(double x0)
    {
        std::vector<double> xs;
        std::vector<std::vector<double>> after;
        auto const watch = [&xs, &after](qd::Step<double> const &) {
            std::vector<double> sorted = xs;
            std::ranges::sort(sorted);
            after.push_back(std::move(sorted));
        };
        double const inf = std::numeric_limits<double>::infinity();
        auto result = qd::integrate<qd::TanhSinh>(qd::observer(watch), [&xs, x0](double x) {
            xs.push_back(x);
            return peak(x, x0);
        }, -inf, inf);
        return {.after = std::move(after), .result = std::move(result)};
    }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // The run on the peak at 0, 1 / (1 + x^2), ends after the five levels every run takes.
    Levels const centred = levels_of(0.0);
    if (!(centred.result && centred.after.size() == 5)) return 2;
    std::vector<double> const &points = centred.after.back();

    // The points within the figure's range, |x| <= 60; then each peak within 5 of its centre, every 0.25.
    std::println("x");
    for (double const x : points)
        if (std::abs(x) <= 60.0) std::println("{:.6g}", x);
    std::println("");
    std::println("x0,x,f");
    for (double const x0 : {0.0, 10.0, 40.0})
        for (int k = -20; k <= 20; ++k) std::println("{},{},{:.4g}", x0, x0 + k / 4.0, peak(x0 + k / 4.0, x0));

    constexpr std::string_view line_table = R"(
x
-59.2864
-43.7545
-32.8596
-25.0828
-19.4398
-15.2812
-12.1714
-9.81335
-8.00166
-6.59213
-5.48225
-4.59814
-3.88595
-3.30597
-2.82859
-2.4315
-2.09774
-1.81425
-1.57093
-1.35986
-1.1748
-1.01075
-0.863703
-0.730391
-0.608122
-0.494645
-0.388051
-0.28669
-0.189098
-0.0939487
0
0.0939487
0.189098
0.28669
0.388051
0.494645
0.608122
0.730391
0.863703
1.01075
1.1748
1.35986
1.57093
1.81425
2.09774
2.4315
2.82859
3.30597
3.88595
4.59814
5.48225
6.59213
8.00166
9.81335
12.1714
15.2812
19.4398
25.0828
32.8596
43.7545
59.2864

x0,x,f
0,-5,0.03846
0,-4.75,0.04244
0,-4.5,0.04706
0,-4.25,0.05246
0,-4,0.05882
0,-3.75,0.06639
0,-3.5,0.07547
0,-3.25,0.08649
0,-3,0.1
0,-2.75,0.1168
0,-2.5,0.1379
0,-2.25,0.1649
0,-2,0.2
0,-1.75,0.2462
0,-1.5,0.3077
0,-1.25,0.3902
0,-1,0.5
0,-0.75,0.64
0,-0.5,0.8
0,-0.25,0.9412
0,0,1
0,0.25,0.9412
0,0.5,0.8
0,0.75,0.64
0,1,0.5
0,1.25,0.3902
0,1.5,0.3077
0,1.75,0.2462
0,2,0.2
0,2.25,0.1649
0,2.5,0.1379
0,2.75,0.1168
0,3,0.1
0,3.25,0.08649
0,3.5,0.07547
0,3.75,0.06639
0,4,0.05882
0,4.25,0.05246
0,4.5,0.04706
0,4.75,0.04244
0,5,0.03846
10,5,0.03846
10,5.25,0.04244
10,5.5,0.04706
10,5.75,0.05246
10,6,0.05882
10,6.25,0.06639
10,6.5,0.07547
10,6.75,0.08649
10,7,0.1
10,7.25,0.1168
10,7.5,0.1379
10,7.75,0.1649
10,8,0.2
10,8.25,0.2462
10,8.5,0.3077
10,8.75,0.3902
10,9,0.5
10,9.25,0.64
10,9.5,0.8
10,9.75,0.9412
10,10,1
10,10.25,0.9412
10,10.5,0.8
10,10.75,0.64
10,11,0.5
10,11.25,0.3902
10,11.5,0.3077
10,11.75,0.2462
10,12,0.2
10,12.25,0.1649
10,12.5,0.1379
10,12.75,0.1168
10,13,0.1
10,13.25,0.08649
10,13.5,0.07547
10,13.75,0.06639
10,14,0.05882
10,14.25,0.05246
10,14.5,0.04706
10,14.75,0.04244
10,15,0.03846
40,35,0.03846
40,35.25,0.04244
40,35.5,0.04706
40,35.75,0.05246
40,36,0.05882
40,36.25,0.06639
40,36.5,0.07547
40,36.75,0.08649
40,37,0.1
40,37.25,0.1168
40,37.5,0.1379
40,37.75,0.1649
40,38,0.2
40,38.25,0.2462
40,38.5,0.3077
40,38.75,0.3902
40,39,0.5
40,39.25,0.64
40,39.5,0.8
40,39.75,0.9412
40,40,1
40,40.25,0.9412
40,40.5,0.8
40,40.75,0.64
40,41,0.5
40,41.25,0.3902
40,41.5,0.3077
40,41.75,0.2462
40,42,0.2
40,42.25,0.1649
40,42.5,0.1379
40,42.75,0.1168
40,43,0.1
40,43.25,0.08649
40,43.5,0.07547
40,43.75,0.06639
40,44,0.05882
40,44.25,0.05246
40,44.5,0.04706
40,44.75,0.04244
40,45,0.03846
)";
    if (!(docs::printed() == docs::quoted(line_table))) return 1;

    // What the page says. The first five levels: 133 points in all, 61 of them within the figure's range; within the
    // half-width of the peak at 0, where it is at least half its height, 17; at 10 one, the points about 2 apart; at
    // 40 none, the two nearest about 11 apart, about five times the peak's width at half its height, 2.
    if (!(points.size() == 133 && in_range(points).size() == 61)) return 3;
    if (!(within(points, 0.0) == 17 && within(points, 10.0) == 1 && within(points, 40.0) == 0)) return 4;
    if (!(gap_at(points, 10.0) > 1.5 && gap_at(points, 10.0) < 2.5)) return 5;
    if (!(gap_at(points, 40.0) > 10.5 && gap_at(points, 40.0) < 11.5 && gap_at(points, 40.0) / 2.0 > 5.0 &&
          gap_at(points, 40.0) / 2.0 < 6.0))
        return 6;
    // Within that range the map puts its points in the same places whatever the integrand: the first five levels of
    // the runs on the peaks at 10 and 40 hold the same points there (how far out a run goes depends on fn).
    Levels const at_10 = levels_of(10.0);
    Levels const at_40 = levels_of(40.0);
    if (!(at_10.after.size() > 5 && at_40.after.size() > 5)) return 7;
    if (!(in_range(at_10.after[4]) == in_range(points) && in_range(at_40.after[4]) == in_range(points))) return 8;
    // Each further level halves the spacing. At 10, four more levels bring the points 0.13 apart, 15 within the
    // half-width, and the run succeeds after 2,069 evaluations.
    if (!(at_10.result && at_10.after.size() == 9 && at_10.after.back().size() == 2069)) return 9;
    if (!(gap_at(at_10.after.back(), 10.0) > 0.125 && gap_at(at_10.after.back(), 10.0) < 0.135 &&
          within(at_10.after.back(), 10.0) == 15))
        return 10;
    // At 40 the tenth and last level leaves them 0.36 apart, 5 within the half-width, and the run ends not_converged
    // after 4,119 evaluations.
    if (!(!at_40.result && at_40.result.error().reason == qd::error::not_converged && at_40.after.size() == 10 &&
          at_40.after.back().size() == 4119))
        return 11;
    if (!(gap_at(at_40.after.back(), 40.0) > 0.355 && gap_at(at_40.after.back(), 40.0) < 0.365 &&
          within(at_40.after.back(), 40.0) == 5))
        return 12;
    return 0;
}
