// The data of the second figure in docs/sphinx/roots/index.md ("What you get back"): (x - 1)^3 computed two ways near
// its root, written out as x^3 - 3x^2 + 3x - 1 by Horner's rule and as a product, at the 201 points x = 1 + k 1e-7,
// k = -100 ... 100, across [1 - 1e-5, 1 + 1e-5]. The program prints them as CSV and compares them with the table in
// the raw string, which docs/sphinx/plots/roots_cancellation.py plots; main returns non-zero unless the two agree and
// what the page says of them holds (PHASE10_PLAN.md), against the exact root 1 and the rounding analysis of
// roots_overview.cpp. Both forms use only + - * /, which IEEE 754 rounds correctly, so the table is the same on every
// target that does not contract a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly).
//
// A decimal step, not a dyadic one: with x = 1 + k 2^-24 every point would be exact and Horner's rule would round on a
// lattice, the same few errors repeating; the rounding a reader meets is better sampled by points that are not. Each x
// is the double nearest (10^7 + k) / 10^7, a quotient of exact integers, correctly rounded, so no fused multiply-add
// can change the grid.

#include "printed.hpp"

#include <algorithm>
#include <cmath>
#include <print>
#include <string_view>
#include <vector>

namespace {

    double written_out(double x) { return ((x - 3.0) * x + 3.0) * x - 1.0; } // x^3 - 3x^2 + 3x - 1
    double product(double x) { return (x - 1.0) * (x - 1.0) * (x - 1.0); }

    struct Point {
        double offset;      // x - 1, exact (Sterbenz)
        double written_out; // fn written out
        double product;     // fn as a product
    };

    // The sign changes of a series between successive nonzero values; a zero between them does not count as one.
    int sign_changes(std::vector<Point> const &points, double Point::*value)
    {
        int changes = 0;
        double last = 0.0;
        for (Point const &p : points) {
            double const y = p.*value;
            if (y == 0.0) continue;
            if (last != 0.0 && (y > 0.0) != (last > 0.0)) ++changes;
            last = y;
        }
        return changes;
    }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    std::vector<Point> points;
    for (int k = -100; k <= 100; ++k) {
        double const x = (1e7 + k) / 1e7;
        points.push_back({.offset = x - 1.0, .written_out = written_out(x), .product = product(x)});
    }

    std::println("offset,written_out,product");
    for (Point const &p : points) std::println("{:.3g},{:.3g},{:.3g}", p.offset, p.written_out, p.product);

    constexpr std::string_view cancellation_table = R"(
offset,written_out,product
-1e-05,-1.11e-15,-1e-15
-9.9e-06,-1.11e-15,-9.7e-16
-9.8e-06,-8.88e-16,-9.41e-16
-9.7e-06,-9.99e-16,-9.13e-16
-9.6e-06,-9.99e-16,-8.85e-16
-9.5e-06,-7.77e-16,-8.57e-16
-9.4e-06,-9.99e-16,-8.31e-16
-9.3e-06,-7.77e-16,-8.04e-16
-9.2e-06,-6.66e-16,-7.79e-16
-9.1e-06,-5.55e-16,-7.54e-16
-9e-06,-5.55e-16,-7.29e-16
-8.9e-06,-8.88e-16,-7.05e-16
-8.8e-06,-4.44e-16,-6.81e-16
-8.7e-06,-7.77e-16,-6.59e-16
-8.6e-06,-7.77e-16,-6.36e-16
-8.5e-06,-5.55e-16,-6.14e-16
-8.4e-06,-5.55e-16,-5.93e-16
-8.3e-06,-7.77e-16,-5.72e-16
-8.2e-06,-4.44e-16,-5.51e-16
-8.1e-06,-6.66e-16,-5.31e-16
-8e-06,-7.77e-16,-5.12e-16
-7.9e-06,-4.44e-16,-4.93e-16
-7.8e-06,-5.55e-16,-4.75e-16
-7.7e-06,-7.77e-16,-4.57e-16
-7.6e-06,-4.44e-16,-4.39e-16
-7.5e-06,-6.66e-16,-4.22e-16
-7.4e-06,-3.33e-16,-4.05e-16
-7.3e-06,-3.33e-16,-3.89e-16
-7.2e-06,-1.11e-16,-3.73e-16
-7.1e-06,-3.33e-16,-3.58e-16
-7e-06,-5.55e-16,-3.43e-16
-6.9e-06,-2.22e-16,-3.29e-16
-6.8e-06,-4.44e-16,-3.14e-16
-6.7e-06,-3.33e-16,-3.01e-16
-6.6e-06,-1.11e-16,-2.87e-16
-6.5e-06,-3.33e-16,-2.75e-16
-6.4e-06,-5.55e-16,-2.62e-16
-6.3e-06,-1.11e-16,-2.5e-16
-6.2e-06,-4.44e-16,-2.38e-16
-6.1e-06,-4.44e-16,-2.27e-16
-6e-06,-1.11e-16,-2.16e-16
-5.9e-06,-3.33e-16,-2.05e-16
-5.8e-06,-1.11e-16,-1.95e-16
-5.7e-06,-2.22e-16,-1.85e-16
-5.6e-06,-2.22e-16,-1.76e-16
-5.5e-06,0,-1.66e-16
-5.4e-06,-1.11e-16,-1.57e-16
-5.3e-06,0,-1.49e-16
-5.2e-06,-1.11e-16,-1.41e-16
-5.1e-06,-1.11e-16,-1.33e-16
-5e-06,0,-1.25e-16
-4.9e-06,-2.22e-16,-1.18e-16
-4.8e-06,-2.22e-16,-1.11e-16
-4.7e-06,2.22e-16,-1.04e-16
-4.6e-06,-1.11e-16,-9.73e-17
-4.5e-06,-3.33e-16,-9.11e-17
-4.4e-06,0,-8.52e-17
-4.3e-06,-1.11e-16,-7.95e-17
-4.2e-06,-4.44e-16,-7.41e-17
-4.1e-06,0,-6.89e-17
-4e-06,-1.11e-16,-6.4e-17
-3.9e-06,0,-5.93e-17
-3.8e-06,0,-5.49e-17
-3.7e-06,-2.22e-16,-5.07e-17
-3.6e-06,2.22e-16,-4.67e-17
-3.5e-06,-1.11e-16,-4.29e-17
-3.4e-06,2.22e-16,-3.93e-17
-3.3e-06,0,-3.59e-17
-3.2e-06,-1.11e-16,-3.28e-17
-3.1e-06,2.22e-16,-2.98e-17
-3e-06,-1.11e-16,-2.7e-17
-2.9e-06,-2.22e-16,-2.44e-17
-2.8e-06,0,-2.2e-17
-2.7e-06,-1.11e-16,-1.97e-17
-2.6e-06,-3.33e-16,-1.76e-17
-2.5e-06,0,-1.56e-17
-2.4e-06,-1.11e-16,-1.38e-17
-2.3e-06,-2.22e-16,-1.22e-17
-2.2e-06,-1.11e-16,-1.06e-17
-2.1e-06,-1.11e-16,-9.26e-18
-2e-06,0,-8e-18
-1.9e-06,0,-6.86e-18
-1.8e-06,0,-5.83e-18
-1.7e-06,0,-4.91e-18
-1.6e-06,0,-4.1e-18
-1.5e-06,2.22e-16,-3.37e-18
-1.4e-06,0,-2.74e-18
-1.3e-06,-1.11e-16,-2.2e-18
-1.2e-06,2.22e-16,-1.73e-18
-1.1e-06,-1.11e-16,-1.33e-18
-1e-06,0,-1e-18
-9e-07,2.22e-16,-7.29e-19
-8e-07,-1.11e-16,-5.12e-19
-7e-07,-2.22e-16,-3.43e-19
-6e-07,0,-2.16e-19
-5e-07,-1.11e-16,-1.25e-19
-4e-07,-1.11e-16,-6.4e-20
-3e-07,-1.11e-16,-2.7e-20
-2e-07,-1.11e-16,-8e-21
-1e-07,0,-1e-21
0,0,0
1e-07,0,1e-21
2e-07,2.22e-16,8e-21
3e-07,2.22e-16,2.7e-20
4e-07,0,6.4e-20
5e-07,0,1.25e-19
6e-07,2.22e-16,2.16e-19
7e-07,0,3.43e-19
8e-07,-1.11e-16,5.12e-19
9e-07,-2.22e-16,7.29e-19
1e-06,-1.11e-16,1e-18
1.1e-06,2.22e-16,1.33e-18
1.2e-06,0,1.73e-18
1.3e-06,2.22e-16,2.2e-18
1.4e-06,2.22e-16,2.74e-18
1.5e-06,0,3.37e-18
1.6e-06,2.22e-16,4.1e-18
1.7e-06,-1.11e-16,4.91e-18
1.8e-06,-1.11e-16,5.83e-18
1.9e-06,0,6.86e-18
2e-06,2.22e-16,8e-18
2.1e-06,0,9.26e-18
2.2e-06,2.22e-16,1.06e-17
2.3e-06,2.22e-16,1.22e-17
2.4e-06,-1.11e-16,1.38e-17
2.5e-06,2.22e-16,1.56e-17
2.6e-06,2.22e-16,1.76e-17
2.7e-06,2.22e-16,1.97e-17
2.8e-06,2.22e-16,2.2e-17
2.9e-06,0,2.44e-17
3e-06,-1.11e-16,2.7e-17
3.1e-06,-1.11e-16,2.98e-17
3.2e-06,-1.11e-16,3.28e-17
3.3e-06,0,3.59e-17
3.4e-06,-1.11e-16,3.93e-17
3.5e-06,2.22e-16,4.29e-17
3.6e-06,2.22e-16,4.67e-17
3.7e-06,2.22e-16,5.07e-17
3.8e-06,0,5.49e-17
3.9e-06,-1.11e-16,5.93e-17
4e-06,-1.11e-16,6.4e-17
4.1e-06,2.22e-16,6.89e-17
4.2e-06,0,7.41e-17
4.3e-06,2.22e-16,7.95e-17
4.4e-06,2.22e-16,8.52e-17
4.5e-06,-1.11e-16,9.11e-17
4.6e-06,0,9.73e-17
4.7e-06,2.22e-16,1.04e-16
4.8e-06,2.22e-16,1.11e-16
4.9e-06,2.22e-16,1.18e-16
5e-06,-1.11e-16,1.25e-16
5.1e-06,0,1.33e-16
5.2e-06,2.22e-16,1.41e-16
5.3e-06,4.44e-16,1.49e-16
5.4e-06,-1.11e-16,1.57e-16
5.5e-06,0,1.66e-16
5.6e-06,2.22e-16,1.76e-16
5.7e-06,2.22e-16,1.85e-16
5.8e-06,2.22e-16,1.95e-16
5.9e-06,4.44e-16,2.05e-16
6e-06,2.22e-16,2.16e-16
6.1e-06,2.22e-16,2.27e-16
6.2e-06,2.22e-16,2.38e-16
6.3e-06,2.22e-16,2.5e-16
6.4e-06,2.22e-16,2.62e-16
6.5e-06,4.44e-16,2.75e-16
6.6e-06,2.22e-16,2.87e-16
6.7e-06,4.44e-16,3.01e-16
6.8e-06,2.22e-16,3.14e-16
6.9e-06,4.44e-16,3.29e-16
7e-06,2.22e-16,3.43e-16
7.1e-06,2.22e-16,3.58e-16
7.2e-06,4.44e-16,3.73e-16
7.3e-06,2.22e-16,3.89e-16
7.4e-06,4.44e-16,4.05e-16
7.5e-06,2.22e-16,4.22e-16
7.6e-06,4.44e-16,4.39e-16
7.7e-06,6.66e-16,4.57e-16
7.8e-06,4.44e-16,4.75e-16
7.9e-06,6.66e-16,4.93e-16
8e-06,6.66e-16,5.12e-16
8.1e-06,6.66e-16,5.31e-16
8.2e-06,4.44e-16,5.51e-16
8.3e-06,6.66e-16,5.72e-16
8.4e-06,6.66e-16,5.93e-16
8.5e-06,6.66e-16,6.14e-16
8.6e-06,8.88e-16,6.36e-16
8.7e-06,8.88e-16,6.59e-16
8.8e-06,6.66e-16,6.81e-16
8.9e-06,6.66e-16,7.05e-16
9e-06,8.88e-16,7.29e-16
9.1e-06,6.66e-16,7.54e-16
9.2e-06,6.66e-16,7.79e-16
9.3e-06,8.88e-16,8.04e-16
9.4e-06,1.11e-15,8.31e-16
9.5e-06,8.88e-16,8.57e-16
9.6e-06,6.66e-16,8.85e-16
9.7e-06,6.66e-16,9.13e-16
9.8e-06,1.11e-15,9.41e-16
9.9e-06,8.88e-16,9.7e-16
1e-05,1.11e-15,1e-15
)";
    if (!(docs::printed() == docs::quoted(cancellation_table))) return 1;

    // What the page says. Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    // The product computes its sign exactly (roots_overview.cpp): one sign change, at 1, where it is 0.
    if (!(sign_changes(points, &Point::product) == 1 && points[100].offset == 0.0 && points[100].product == 0.0))
        return 2;
    // The written-out form changes sign 31 times between nonzero values, more than 30, and is exactly 0 at 35 points.
    int const zeros = static_cast<int>(std::ranges::count(points, 0.0, &Point::written_out));
    if (!(sign_changes(points, &Point::written_out) == 31 && zeros == 35)) return 3;
    // Its farthest zero or wrong sign lies beyond 4e-6 from 1, and, as the rounding analysis requires, within 8.2e-6.
    double farthest = 0.0;
    for (Point const &p : points)
        if (p.written_out == 0.0 || (p.written_out > 0.0) != (p.offset > 0.0))
            farthest = std::max(farthest, std::abs(p.offset));
    if (!(farthest > 4e-6 && farthest <= 8.2e-6)) return 4;
    // "0 or a few multiples of 1.1e-16 of either sign": near 1, t * x - 1 is exact (Sterbenz) and t * x is a double
    // within a factor 2 of 1, so every written-out value is a whole multiple of 2^-53 = 1.1e-16; and by the running
    // error analysis of roots_overview.cpp it lies within 5.5e-16 of (x - 1)^3, which the product computes to a
    // relative 2.3e-16 (two roundings of an exact x - 1), under 1e-30 here.
    for (Point const &p : points) {
        double const units = p.written_out * 0x1p53;
        if (!(units == std::round(units) && std::abs(p.written_out - p.product) <= 5.5e-16 + 1e-30)) return 5;
    }
    return 0;
}
