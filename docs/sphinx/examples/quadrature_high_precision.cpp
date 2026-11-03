// The program of docs/sphinx/quadrature/high-precision.md: integrals in cpp_bin_float_50, a Boost.Multiprecision type
// with 50 decimal digits: the default request against one for 40 digits, Gauss-Kronrod against TanhSinh at that
// request, and the pendulum's K of the tutorial to 40 digits, checked against the arithmetic-geometric mean in the
// same type. The page shows the parts between the markers and the output in the raw strings; main returns non-zero
// unless the program prints that output and what the page says holds. The first integrand uses only + - * /, which
// Boost rounds correctly, so its output is the same everywhere; the pendulum's goes through Boost's sqrt, so only its
// first 40 digits and its count are quoted, and its error is checked by a bound.

// clang-format off
// <includes>
#include <numerixx/quadrature.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>

#include <cmath>
#include <limits>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <string>

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <request>
    namespace qd = nxx::quadrature;
    using bin50 = boost::multiprecision::cpp_bin_float_50;
    int evaluations = 0;
    auto const fn = [&evaluations](bin50 const &x) {
        ++evaluations;
        return bin50(bin50(1) / (bin50(1) + x * x));
    };
    auto const by_default = qd::integrate(fn, bin50(0), bin50(1));
    if (by_default)
        std::println("default:      abs_error {}, {} evaluations", by_default->abs_error.str(2), evaluations);
    evaluations = 0;
    auto const forty = qd::integrate(fn, bin50(0), bin50(1), {.rel_tol = bin50("1e-40")});
    if (forty) std::println("rel_tol 1e-40: abs_error {}, {} evaluations", forty->abs_error.str(2), evaluations);
    evaluations = 0;
    auto const tanh_sinh = qd::integrate<qd::TanhSinh>(fn, bin50(0), bin50(1), {.rel_tol = bin50("1e-40")});
    if (tanh_sinh) std::println("TanhSinh:     abs_error {}, {} evaluations", tanh_sinh->abs_error.str(2), evaluations);
    // </request>

    // <pendulum>
    using std::sqrt;
    bin50 const kp = sqrt(bin50(1) / 2); // k' = cos 45 degrees, correctly rounded in the type
    bin50 const k2 = bin50(1) - kp * kp;
    int count = 0;
    auto const integrand = [kp, k2, &count](bin50 const &u) {
        ++count;
        bin50 const w = u * (bin50(2) - u);
        return bin50(bin50(1) / (sqrt(w) * sqrt(kp * kp + k2 * w)));
    };
    auto const K = qd::integrate<qd::TanhSinh>(integrand, bin50(0), bin50(1), {.rel_tol = bin50("1e-40")});
    if (K) std::println("K = {}, {} evaluations", K->value.str(40), count);
    // </pendulum>

    // <reference>
    using std::atan;
    bin50 const pi = bin50(4) * atan(bin50(1));
    bin50 const eps = std::numeric_limits<bin50>::epsilon();
    bin50 a = 1;
    bin50 b = kp;
    while (a - b > bin50(4) * eps * a) {
        bin50 const mean = (a + b) / 2;
        b = sqrt(a * b);
        a = mean;
    }
    bin50 const exact = pi / (a + b); // K = pi / (2 AGM(1, k'))
    // </reference>
    // clang-format on

    constexpr std::string_view request_output = R"(
default:      abs_error 1.8e-26, 420 evaluations
rel_tol 1e-40: abs_error 7.6e-41, 5430 evaluations
TanhSinh:     abs_error 2.1e-49, 325 evaluations
)";
    constexpr std::string_view pendulum_output = R"(
K = 1.854074677301371918433850347195260046218, 345 evaluations
)";
    std::string const expected = std::string(docs::quoted(request_output)) + std::string(docs::quoted(pendulum_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    if (!by_default || !forty || !tanh_sinh || !K) return 2;
    using std::abs;
    bin50 const quarter_pi = pi / 4; // the first integral
    // "The default asks for about 25 digits": sqrt(eps), here 2^-83.5 = 7.3e-26, since eps is 2^-167 for
    // cpp_bin_float_50's 168 bits (the quoted abs_error 2.1e-49 is 50 2^-167 pi / 4, TanhSinh's floor).
    bin50 const half = qd::Criteria<bin50>{}.rel_tol;
    if (!(half > bin50("1e-26") && half < bin50("1e-25"))) return 3;
    // "Its value is far more accurate than its abs_error says" (5.5e-45 against 1.8e-26), and every estimate covers
    // its error. The reference pi / 4 errs by a few ulps of the type (atan and the product), 1e-49, which the margins
    // below leave room for.
    bin50 const margin = bin50("1e-48");
    if (!(abs(by_default->value - quarter_pi) < bin50("1e-40") && by_default->abs_error > bin50("1e-27"))) return 4;
    if (!(abs(forty->value - quarter_pi) <= forty->abs_error + margin
          && abs(tanh_sinh->value - quarter_pi) <= tanh_sinh->abs_error + margin))
        return 5;
    // "Both meet the request", and "TanhSinh needs about a sixteenth of the evaluations" (325 against 5,430, above).
    if (!(forty->abs_error <= bin50("1e-40") * quarter_pi && tanh_sinh->abs_error <= bin50("1e-40") * quarter_pi))
        return 6;
    // The pendulum: K within 1e-48 of the AGM (its own error a few ulps, the AGM's at most about 14 ulps of the type,
    // as in the tutorial's analysis), far inside the 40 digits printed.
    if (!(abs(K->value - exact) <= bin50("1e-48"))) return 7;
    return 0;
}
