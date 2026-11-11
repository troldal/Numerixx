// The data of the first figure in docs/sphinx/roots/index.md ("Two ways to start"): the width of the bracket after
// each evaluation of fn, for Bisection and Brent's method on x^2 - 2 over [1, 2], from the steps observer(o) reports.
// The program prints them as CSV and compares them with the table in the raw string, which
// docs/sphinx/plots/roots_closing_in.py plots; main returns non-zero unless the two agree and what the page says of
// them holds (PHASE10_PLAN.md), against sqrt(2). fn and the methods use only + - * /, which IEEE 754 rounds correctly,
// so the table is the same on every target that does not contract a * b + c into a fused multiply-add (x86-64 without
// -mfma, WebAssembly).
//
// Each point is one evaluation, so there is no sampling to choose. The flat step of Brent's width at the 5th and 6th
// evaluations is the method's: two interpolated points in a row land below the root, the second next to the first, so
// the far end stays and the bracket narrows little; the next point lands just above the root. The same step shows on
// [1.1, 2], [1, 1.9], [0.7, 2.3] and [1, 3] (a probe), and on Wallis's cubic over [2, 3], which main checks.

#include <numerixx/roots.hpp>

#include "printed.hpp"

#include <cmath>
#include <cstddef>
#include <expected>
#include <print>
#include <string_view>
#include <utility>
#include <vector>

namespace rt = nxx::roots;

namespace {

    struct Point {
        int evaluations; // of fn, when the step was reported
        double width;    // upper - lower
    };

    struct Run {
        std::vector<rt::Step<double>> steps;
        std::vector<Point> points;
        std::expected<rt::Root<double>, rt::error> root;
        int evaluations;
    };

    // A solve by M, each step recorded with the evaluations of fn so far. upper - lower is exact for the brackets here,
    // whose ends lie within a factor 2 of each other (Sterbenz).
    template<typename M, typename F>
    Run run(F const &f, rt::Bracket<double> const &bracket)
    {
        int n = 0;
        std::vector<rt::Step<double>> steps;
        std::vector<Point> points;
        auto const fn = [&n, &f](double x) {
            ++n;
            return f(x);
        };
        auto const record = [&](rt::Step<double> const &s) {
            steps.push_back(s);
            points.push_back({.evaluations = n, .width = s.upper - s.lower});
        };
        auto const root = rt::solve<M>(rt::observer(record), fn, bracket);
        return {.steps = steps, .points = points, .root = root, .evaluations = n};
    }

    // A step that narrows the bracket by less than 5%, followed by one that narrows it more than 10,000-fold.
    bool flat_then_cliff(std::vector<Point> const &points)
    {
        for (std::size_t i = 2; i < points.size(); ++i)
            if (points[i - 1].width > 0.95 * points[i - 2].width && points[i].width < 1e-4 * points[i - 1].width)
                return true;
        return false;
    }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    auto const square = [](double x) { return x * x - 2.0; };
    Run const bisection = run<rt::Bisection>(square, rt::Bracket{1.0, 2.0});
    Run const brent = run<rt::Brent>(square, rt::Bracket{1.0, 2.0});
    if (!bisection.root || !brent.root) return 2;

    std::println("method,evaluations,width");
    for (auto const &[name, r] : {std::pair{"bisection", &bisection}, std::pair{"brent", &brent}})
        for (Point const &p : r->points) std::println("{},{},{:.3g}", name, p.evaluations, p.width);

    constexpr std::string_view closing_in_table = R"(
method,evaluations,width
bisection,2,1
bisection,3,0.5
bisection,4,0.25
bisection,5,0.125
bisection,6,0.0625
bisection,7,0.0312
bisection,8,0.0156
bisection,9,0.00781
bisection,10,0.00391
bisection,11,0.00195
bisection,12,0.000977
bisection,13,0.000488
bisection,14,0.000244
bisection,15,0.000122
bisection,16,6.1e-05
bisection,17,3.05e-05
bisection,18,1.53e-05
bisection,19,7.63e-06
bisection,20,3.81e-06
bisection,21,1.91e-06
bisection,22,9.54e-07
bisection,23,4.77e-07
bisection,24,2.38e-07
bisection,25,1.19e-07
bisection,26,5.96e-08
bisection,27,2.98e-08
bisection,28,1.49e-08
bisection,29,7.45e-09
bisection,30,3.73e-09
bisection,31,1.86e-09
bisection,32,9.31e-10
bisection,33,4.66e-10
bisection,34,2.33e-10
bisection,35,1.16e-10
bisection,36,5.82e-11
bisection,37,2.91e-11
bisection,38,1.46e-11
bisection,39,7.28e-12
bisection,40,3.64e-12
bisection,41,1.82e-12
bisection,42,9.09e-13
bisection,43,4.55e-13
bisection,44,2.27e-13
bisection,45,1.14e-13
bisection,46,5.68e-14
bisection,47,2.84e-14
bisection,48,1.42e-14
bisection,49,7.11e-15
bisection,50,3.55e-15
bisection,51,1.78e-15
bisection,52,8.88e-16
brent,2,1
brent,3,0.667
brent,4,0.0857
brent,5,0.00498
brent,6,0.00483
brent,7,2.42e-07
brent,8,4.15e-14
brent,9,6.66e-16
)";
    if (!(docs::printed() == docs::quoted(closing_in_table))) return 1;

    // What the page says. Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    // "52 and 9" evaluations, the last step reported after the last evaluation, the first after the two ends.
    if (!(bisection.evaluations == 52 && brent.evaluations == 9)) return 3;
    for (Run const *r : {&bisection, &brent})
        if (!(r->points.front().evaluations == 2 && r->points.front().width == 1.0
              && r->points.back().evaluations == r->evaluations))
            return 4;
    // "Bisection's width halves with every evaluation": exactly, since every end is a dyadic number of [1, 2], from 1
    // at the second evaluation to 2^-50 = 8.9e-16 at the 52nd.
    for (std::size_t i = 1; i < bisection.points.size(); ++i)
        if (!(bisection.points[i].width == bisection.points[i - 1].width / 2.0
              && bisection.points[i].evaluations == bisection.points[i - 1].evaluations + 1))
            return 5;
    if (!(bisection.points.back().width == 0x1p-50)) return 6;
    // "Brent's last width is the abs_error it returns", below 1e-15. Its point is an end of the last bracket, so the
    // larger distance from the point to an end is the width, and abs_error, the smaller of the two, is the width.
    if (!(brent.root->abs_error == brent.points.back().width && brent.root->abs_error < 1e-15)) return 7;
    // Both roots within their bounds of sqrt(2) = 1.41421356237309504880 (20 digits), plus c: x * x rounds by at most
    // half an ulp of 2 and the subtraction is exact, so the computed sign change lies within 2.2e-16 / |f'| =
    // 2.2e-16 / 2.83 < 8e-17 of sqrt(2); and half an ulp of sqrt(2), 1.2e-16, for the reference held in a double.
    double const sqrt2 = 1.41421356237309504880;
    for (Run const *r : {&bisection, &brent})
        if (!(std::abs(r->root->value - sqrt2) <= r->root->abs_error + 8e-17 + 1.2e-16)) return 8;
    // "A step can land next to the end it replaces and narrow the bracket little; the next makes up for it": at the
    // 6th evaluation the upper end stays and the new lower end lies between the old one and the root, so the width
    // falls by less than 5%; at the 7th it falls more than 10,000-fold.
    rt::Step<double> const &fifth = brent.steps[3];
    rt::Step<double> const &sixth = brent.steps[4];
    if (!(brent.points[3].evaluations == 5 && sixth.upper == fifth.upper && fifth.lower < sixth.lower
          && sixth.lower < sqrt2 && brent.points[4].width > 0.95 * brent.points[3].width
          && brent.points[5].width < 1e-4 * brent.points[4].width))
        return 9;
    // The same on Wallis's cubic over [2, 3]: the step is the method's, not this function's.
    Run const wallis = run<rt::Brent>([](double x) { return (x * x - 2.0) * x - 5.0; }, rt::Bracket{2.0, 3.0});
    if (!(wallis.root && flat_then_cliff(brent.points) && flat_then_cliff(wallis.points))) return 10;
    return 0;
}
