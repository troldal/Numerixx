// The data of the figure in docs/sphinx/quadrature/rough-integrands.md ("A singularity at an end"): the substitution
// TanhSinh makes on [0, 1], x = (1 + tanh(sinh t)) / 2, and the integrand 1 / sqrt x after it, f(x(t)) dx/dt, as
// functions of t. The program prints both as CSV and compares them with the table in the raw string, which
// docs/sphinx/plots/quadrature_tanh_sinh_map.py plots; main returns non-zero unless the two agree, the map is the one
// the library uses, and what the page says holds (PHASE8_PLAN.md).
//
// The map needs e^y, and a libm's exp may differ in its last bits between platforms. So the program computes e^y
// itself, from + - * / alone (halving y, a Taylor series, squaring back), which IEEE 754 rounds the same everywhere;
// the program has no a * b + c that a compiler could contract into a fused multiply-add. The map's relative error stays
// below 1e-12 (largest, 8.8e-13 near t = -3.8, mostly the error of sinh t multiplied by |2 sinh t|); the table prints
// six digits.

#include <numerixx/quadrature.hpp>

#include "printed.hpp"

#include <algorithm>
#include <cmath>
#include <print>
#include <string_view>
#include <vector>

namespace {

    namespace qd = nxx::quadrature;

    // e^y: y halved k times to |r| <= 1/16 (exact), e^r from its Taylor series to r^12 / 12!, far below rounding for
    // |r| <= 1/16, then squared k times.
    double exp_of(double y)
    {
        double r = y;
        int k = 0;
        for (; r > 0.0625 || r < -0.0625; ++k) r /= 2.0;
        double term = 1.0;
        double sum = 1.0;
        for (int n = 1; n <= 12; ++n) {
            term = term * r / n;
            sum += term;
        }
        for (; k > 0; --k) sum *= sum;
        return sum;
    }

    double sinh_of(double t) { return (exp_of(t) - exp_of(-t)) / 2.0; }
    double cosh_of(double t) { return (exp_of(t) + exp_of(-t)) / 2.0; }

    // The map onto [0, 1] and its complement, each formed directly so that neither is a difference of nearly equal
    // numbers: x = 1 / (1 + e^(-2 sinh t)) and 1 - x = 1 / (1 + e^(2 sinh t)).
    double x_of(double t) { return 1.0 / (1.0 + exp_of(-2.0 * sinh_of(t))); }
    double one_minus_x_of(double t) { return 1.0 / (1.0 + exp_of(2.0 * sinh_of(t))); }

    // The integrand in t for f = 1 / sqrt x: f(x) dx/dt, with dx/dt = 2 x (1 - x) cosh t, so 2 sqrt(x) (1 - x) cosh t.
    double g_of(double t) { return 2.0 * std::sqrt(x_of(t)) * one_minus_x_of(t) * cosh_of(t); }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // t from -4 to 4 in steps of 1/16; the points of the step 1/2, which the figure marks, flagged.
    std::println("t,x,g,marked");
    for (int k = -64; k <= 64; ++k) {
        double const t = k / 16.0;
        std::println("{},{:.6g},{:.6g},{}", t, x_of(t), g_of(t), k % 8 == 0 ? 1 : 0);
    }

    constexpr std::string_view tanh_sinh_map_table = R"(
t,x,g,marked
-4,1.97824e-24,7.6818e-11,1
-3.9375,5.41265e-23,3.77489e-10,0
-3.875,1.21209e-21,1.6782e-09,0
-3.8125,2.24867e-20,6.79082e-09,0
-3.75,3.49577e-19,2.51545e-08,0
-3.6875,4.60299e-18,8.57537e-08,0
-3.625,5.18553e-17,2.7041e-07,0
-3.5625,5.04559e-16,7.92464e-07,0
-3.5,4.27815e-15,2.16798e-06,1
-3.4375,3.18754e-14,5.55985e-06,0
-3.375,2.10337e-13,1.34187e-05,0
-3.3125,1.23834e-12,3.05912e-05,0
-3.25,6.54994e-12,6.6104e-05,0
-3.1875,3.1328e-11,0.000135837,0
-3.125,1.36328e-10,0.000266256,0
-3.0625,5.42859e-10,0.000499251,0
-3,1.98877e-09,0.000897947,1
-2.9375,6.73723e-09,0.0015531,0
-2.875,2.12054e-08,0.00258941,0
-2.8125,6.22915e-08,0.00417091,0
-2.75,1.71496e-07,0.00650442,0
-2.6875,4.44266e-07,0.00983998,0
-2.625,1.08695e-06,0.0144677,0
-2.5625,2.52041e-06,0.0207104,0
-2.5,5.55721e-06,0.028912,1
-2.4375,1.16871e-05,0.0394225,0
-2.375,2.35115e-05,0.05258,0
-2.3125,4.53693e-05,0.0686917,0
-2.25,8.41911e-05,0.0880151,0
-2.1875,0.000150604,0.11074,0
-2.125,0.000260287,0.136975,0
-2.0625,0.000435546,0.166734,0
-2,0.000707037,0.199933,1
-1.9375,0.00111553,0.236386,0
-1.875,0.00171357,0.275806,0
-1.8125,0.00256687,0.317809,0
-1.75,0.00375529,0.361929,0
-1.6875,0.00537307,0.40762,0
-1.625,0.00752846,0.454277,0
-1.5625,0.0103422,0.501247,0
-1.5,0.0139454,0.547848,1
-1.4375,0.018476,0.593387,0
-1.375,0.0240744,0.637179,0
-1.3125,0.0308793,0.678572,0
-1.25,0.0390215,0.716961,0
-1.1875,0.0486194,0.751813,0
-1.125,0.0597732,0.782684,0
-1.0625,0.0725606,0.80923,0
-1,0.0870338,0.831222,1
-0.9375,0.103217,0.848548,0
-0.875,0.121104,0.861211,0
-0.8125,0.140664,0.869322,0
-0.75,0.161836,0.873091,0
-0.6875,0.184537,0.872808,0
-0.625,0.208666,0.868823,0
-0.5625,0.234106,0.861524,0
-0.5,0.260728,0.85132,1
-0.4375,0.288397,0.838619,0
-0.375,0.316976,0.823807,0
-0.3125,0.346327,0.807241,0
-0.25,0.376314,0.789231,0
-0.1875,0.406802,0.770036,0
-0.125,0.437663,0.749861,0
-0.0625,0.46877,0.728853,0
0,0.5,0.707107,1
0.0625,0.53123,0.684666,0
0.125,0.562337,0.661534,0
0.1875,0.593198,0.63768,0
0.25,0.623686,0.61305,0
0.3125,0.653673,0.587579,0
0.375,0.683024,0.561205,0
0.4375,0.711603,0.533877,0
0.5,0.739272,0.505573,1
0.5625,0.765894,0.476309,0
0.625,0.791334,0.446146,0
0.6875,0.815463,0.415201,0
0.75,0.838164,0.383647,0
0.8125,0.859336,0.351715,0
0.875,0.878896,0.319684,0
0.9375,0.896783,0.287877,0
1,0.912966,0.256646,1
1.0625,0.927439,0.22635,0
1.125,0.940227,0.197344,0
1.1875,0.951381,0.169956,0
1.25,0.960978,0.144474,0
1.3125,0.969121,0.121127,0
1.375,0.975926,0.100076,0
1.4375,0.981524,0.0814124,0
1.5,0.986055,0.0651517,1
1.5625,0.989658,0.0512409,0
1.625,0.992472,0.0395654,0
1.6875,0.994627,0.0299597,0
1.75,0.996245,0.0222209,0
1.8125,0.997433,0.0161223,0
1.875,0.998286,0.0114268,0
1.9375,0.998884,0.00789959,0
2,0.999293,0.00531814,1
2.0625,0.999564,0.00348045,0
2.125,0.99974,0.00221016,0
2.1875,0.999849,0.00135911,0
2.25,0.999916,0.000807623,0
2.3125,0.999955,0.000462696,0
2.375,0.999976,0.000254956,0
2.4375,0.999988,0.000134772,0
2.5,0.999994,6.81566e-05,1
2.5625,0.999997,3.28795e-05,0
2.625,0.999999,1.50836e-05,0
2.6875,1,6.55867e-06,0
2.75,1,2.69362e-06,0
2.8125,1,1.04099e-06,0
2.875,1,3.77072e-07,0
2.9375,1,1.27479e-07,0
3,1,4.00445e-08,1
3.0625,1,1.16322e-08,0
3.125,1,3.1088e-09,0
3.1875,1,7.60302e-10,0
3.25,1,1.69179e-10,0
3.3125,1,3.40421e-11,0
3.375,1,6.15414e-12,0
3.4375,1,9.92638e-13,0
3.5,1,1.41802e-13,1
3.5625,1,1.78006e-14,0
3.625,1,1.94724e-15,0
3.6875,1,1.83981e-16,0
3.75,1,1.48726e-17,0
3.8125,1,1.01832e-18,0
3.875,1,5.84267e-20,0
3.9375,1,2.77721e-21,0
4,1,1.08045e-22,1
)";
    if (!(docs::printed() == docs::quoted(tanh_sinh_map_table))) return 1;

    // The map is the library's: among the points TanhSinh evaluates, those of its first level, the step 1, are x(t)
    // for t = -4 to 3 (and beyond -4; t = 4 rounds onto 1). Near 0 a point may differ from x(t) by the map's error,
    // below 1e-12 of x, which keeps the library's full relative accuracy there; near 1 also by two ulps of the numbers
    // just below 1, 2.2e-16, since both are rounded there.
    std::vector<double> points;
    if (!qd::integrate<qd::TanhSinh>([&points](double x) {
            points.push_back(x);
            return 1.0 / std::sqrt(x);
        }, 0.0, 1.0))
        return 2;
    auto const evaluated = [&points](double x, double tolerance) {
        return std::ranges::any_of(points, [=](double p) { return std::abs(p - x) <= tolerance; });
    };
    for (int t = -4; t <= 3; ++t)
        if (!evaluated(x_of(t), t <= 0 ? 1e-12 * x_of(t) : 1e-12 * one_minus_x_of(t) + 2.3e-16)) return 3;

    // What the page says: at t = -3 the point is about 2e-9; the integrand in t has fallen below 1e-10 at t = -4 and
    // t = 4, the ends of the figure.
    if (!(x_of(-3.0) > 1.5e-9 && x_of(-3.0) < 2.5e-9)) return 4;
    // Each further step of 1/2 towards -4 multiplies the distance to 0 by a factor that itself shrinks.
    for (double t = -1.0; t > -3.75; t -= 0.5)
        if (!(x_of(t - 1.0) / x_of(t - 0.5) < x_of(t - 0.5) / x_of(t))) return 8;
    if (!(g_of(-4.0) < 1e-10 && g_of(4.0) < 1e-10)) return 5;
    // On the singular side it falls more slowly: g(-t) > g(t) for t from 1 to 4, where 1 / sqrt x ~ e^(-sinh t) takes
    // back half of the exponent of dx/dt ~ 2 cosh t e^(-2 sinh t).
    for (int k = 16; k <= 64; ++k)
        if (!(g_of(-k / 16.0) > g_of(k / 16.0))) return 9;
    // And what the alt text says: on the figure's grid the integrand is largest, about 0.87, at t = -0.75, and at t = 4
    // it is about 1e-22.
    int top = -64;
    for (int k = -64; k <= 64; ++k)
        if (g_of(k / 16.0) > g_of(top / 16.0)) top = k;
    if (!(top == -12 && g_of(-0.75) > 0.86 && g_of(-0.75) < 0.88)) return 6;
    if (!(g_of(4.0) > 0.9e-22 && g_of(4.0) < 1.2e-22)) return 7;
    return 0;
}
