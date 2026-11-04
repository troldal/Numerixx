// The data of the figure in docs/sphinx/quadrature/high-precision.md ("Ask for the digits you want"): how many
// evaluations the default method and TanhSinh need for 1 / (1 + x^2) over [0, 1], which is pi / 4, in cpp_bin_float_50
// at requests of 4 to 48 digits. The program prints the counts as CSV and compares them with the table in the raw
// string, which docs/sphinx/plots/quadrature_digits.py plots; main returns non-zero unless the two agree, each result
// meets its request, and what the page says holds (PHASE8_PLAN.md). The integrand and the methods use only + - * /,
// which Boost rounds correctly, so the counts are the same everywhere.

#include <numerixx/quadrature.hpp>

#include "printed.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>

#include <cmath>
#include <cstddef>
#include <print>
#include <string>
#include <string_view>
#include <vector>

namespace {

    namespace qd = nxx::quadrature;
    using bin50 = boost::multiprecision::cpp_bin_float_50;

    struct Count {
        int digits;
        std::size_t kronrod;
        std::size_t tanh_sinh;
    };

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // pi / 4 to 60 digits, from mpmath: an independent value to check each result against.
    bin50 const exact("0.785398163397448309615660845819875721049292349843776455243736");

    std::vector<Count> counts;
    for (int d = 4; d <= 48; ++d) {
        bin50 const rel_tol("1e-" + std::to_string(d));
        std::size_t n = 0;
        auto const fn = [&n](bin50 const &x) {
            ++n;
            return bin50(bin50(1) / (bin50(1) + x * x));
        };
        auto const kronrod = qd::integrate(fn, bin50(0), bin50(1), {.rel_tol = rel_tol});
        std::size_t const kronrod_n = n;
        n = 0;
        auto const tanh_sinh = qd::integrate<qd::TanhSinh>(fn, bin50(0), bin50(1), {.rel_tol = rel_tol});
        if (!(kronrod && tanh_sinh)) return 2;
        // Each result meets its request: its actual error is below rel_tol times the integral.
        using boost::multiprecision::abs;
        if (!(abs(kronrod->value - exact) <= rel_tol * exact && abs(tanh_sinh->value - exact) <= rel_tol * exact))
            return 3;
        counts.push_back({.digits = d, .kronrod = kronrod_n, .tanh_sinh = n});
    }

    std::println("digits,default,tanh_sinh");
    for (Count const &c : counts) std::println("{},{},{}", c.digits, c.kronrod, c.tanh_sinh);

    constexpr std::string_view digits_table = R"(
digits,default,tanh_sinh
4,150,172
5,150,172
6,150,172
7,150,172
8,150,172
9,150,172
10,150,172
11,150,172
12,150,172
13,150,172
14,150,172
15,150,172
16,150,172
17,150,172
18,150,172
19,150,172
20,150,172
21,150,172
22,150,172
23,210,172
24,300,172
25,390,172
26,450,172
27,450,172
28,720,172
29,870,325
30,990,325
31,1020,325
32,1470,325
33,1890,325
34,2100,325
35,2250,325
36,2850,325
37,3750,325
38,4230,325
39,4560,325
40,5430,325
41,7230,325
42,8460,325
43,9150,325
44,10050,325
45,13830,325
46,16680,325
47,18240,325
48,18990,325
)";
    if (!(docs::printed() == docs::quoted(digits_table))) return 1;

    // What the page says. Up to 22 digits the default needs only its first 150 evaluations, fewer than TanhSinh's
    // 172; at 23 it begins to cut; its count never falls as the request tightens, and reaches 18,990 at 48 digits.
    for (Count const &c : counts)
        if (!(c.digits > 22 || (c.kronrod == 150 && c.tanh_sinh == 172))) return 4;
    if (!(counts[23 - 4].kronrod > 150)) return 5;
    for (std::size_t i = 1; i < counts.size(); ++i)
        if (!(counts[i].kronrod >= counts[i - 1].kronrod)) return 6;
    if (!(counts.back().digits == 48 && counts.back().kronrod == 18990)) return 7;
    // From about 30 digits its cost doubles with about every four more digits: from 30 to 48, 18 digits over
    // log2(18990 / 990), 4.2 (halving every panel doubles the evaluations, once there are many, and divides the 7-point
    // rule's error by 2^14, 4.2 digits; from 23 to 30 a bisection's 30 new evaluations for a panel of 15 make it
    // about three). log2 from a libm is fine here: only the range is asserted.
    double const per_doubling = 18.0 / std::log2(static_cast<double>(counts.back().kronrod) /
                                                 static_cast<double>(counts[30 - 4].kronrod));
    if (!(per_doubling > 3.8 && per_doubling < 4.6)) return 8;
    // TanhSinh's cost rises in one step: 172 evaluations up to 28 digits, 325 from 29 to 48.
    for (Count const &c : counts)
        if (!(c.tanh_sinh == (c.digits <= 28 ? 172u : 325u))) return 9;
    // At the default request, sqrt(eps), about 25 digits: 420 evaluations for the default and 172 for TanhSinh.
    {
        std::size_t n = 0;
        auto const fn = [&n](bin50 const &x) {
            ++n;
            return bin50(bin50(1) / (bin50(1) + x * x));
        };
        if (!(qd::integrate(fn, bin50(0), bin50(1)) && n == 420)) return 10;
        n = 0;
        if (!(qd::integrate<qd::TanhSinh>(fn, bin50(0), bin50(1)) && n == 172)) return 11;
    }
    return 0;
}
