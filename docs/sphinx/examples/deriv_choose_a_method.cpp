// The program of docs/sphinx/deriv/choose-a-method.md: evaluations and accuracy of the methods on 1/(1 + x^2) at 0.5,
// whose derivative is -2x / (1 + x^2)^2 = -0.64. The page shows the parts between the markers and the output in the
// raw string; main returns non-zero unless the program prints that output and what the page says holds. fn uses only
// + - * /, so the output is the same on every target that does not contract a * b + c into a fused multiply-add
// (x86-64 without -mfma, WebAssembly); the cpp_bin_float_50 part is Boost's own arithmetic.

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>

#include <cmath>
#include <print>
// </includes>

// <count>
// The relative error of derivative<M> on fn at x, and the number of calls of fn it made.
template<typename M, typename F>
void compare(char const *name, F const &fn, double x, double exact)
{
    int calls = 0;
    auto const counted = [&](double t) {
        ++calls;
        return fn(t);
    };
    auto const d = nxx::deriv::derivative<M>(counted, x);
    if (d) std::println("{:<12} {:>2} calls   relative error {:.1e}", name, calls,
                        std::abs(d->value - exact) / std::abs(exact));
}
// </count>
// clang-format on

#include "printed.hpp"

#include <limits>
#include <string>

namespace {

    // The relative error of derivative<M>.
    template<typename M, typename F>
    double relative_error(F const &fn, double x, double exact)
    {
        auto const d = nxx::deriv::derivative<M>(fn, x);
        return d ? std::abs(d->value - exact) / std::abs(exact) : std::numeric_limits<double>::quiet_NaN();
    }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <compare>
    namespace dv = nxx::deriv;
    auto const fn = [](double x) { return 1.0 / (1.0 + x * x); };
    double const exact = 0.64; // f'(0.5) = -0.64
    compare<dv::Central<2>>("Central<2>", fn, 0.5, -exact);
    compare<dv::Central<4>>("Central<4>", fn, 0.5, -exact);
    compare<dv::Central<6>>("Central<6>", fn, 0.5, -exact);
    compare<dv::Ridders>("Ridders", fn, 0.5, -exact);
    // </compare>

    // <more-digits>
    using bin50 = boost::multiprecision::cpp_bin_float_50; // 50 decimal digits
    auto const root = [](bin50 const &x) {
        using std::sqrt;
        return bin50(sqrt(x));
    };
    auto const d50 = dv::derivative<dv::Ridders>(root, bin50(2));
    if (d50) std::println("{}, abs_error {}", d50->value.str(40), d50->abs_error->str(2));
    // </more-digits>
    // clang-format on

    constexpr std::string_view more_digits_output = R"(
0.3535533905932737622004221810524245196424, abs_error 6.5e-47
)";
    constexpr std::string_view output = R"(
Central<2>    2 calls   relative error 3.4e-11
Central<4>    4 calls   relative error 2.0e-13
Central<6>    6 calls   relative error 4.9e-13
Ridders      14 calls   relative error 5.3e-14
)";
    if (!(docs::printed() == std::string(docs::quoted(output)) + std::string(docs::quoted(more_digits_output))))
        return 1;

    // Each check states what must hold, so that a NaN fails it. The bounds, truncation plus rounding, for 1/(1 + x^2)
    // at 0.5, with f^(n) = (-1)^n n! Im[(x - i)^-(n+1)]: Central<2>, f''' = 3.69, 3.5e-11 + 4.6e-11; Central<4>,
    // f^(5) = -21.6, 3.4e-13 + 5.6e-13; Central<6>, f^(7) = -1110, 4.7e-13 + 9e-14. Ridders reaches about 1e-13 or
    // better.
    double const e2 = relative_error<dv::Central<2>>(fn, 0.5, -exact);
    double const e4 = relative_error<dv::Central<4>>(fn, 0.5, -exact);
    double const e6 = relative_error<dv::Central<6>>(fn, 0.5, -exact);
    double const er = relative_error<dv::Ridders>(fn, 0.5, -exact);
    if (!(e2 < 2e-10 && e4 < 3e-12 && e6 < 2e-12 && er < 1e-13)) return 2;
    // "Each higher order is more accurate than Central<2>", and Ridders more than any stencil.
    if (!(e4 < e2 && e6 < e2 && er < e4 && er < e6)) return 3;
    // "About two digits over Central<2>": at least a factor 30.
    if (!(e4 < e2 / 30 && e6 < e2 / 30)) return 5;
    auto const r = dv::derivative<dv::Ridders>(fn, 0.5);
    if (!r || !r->abs_error || !(std::abs(r->value + exact) <= *r->abs_error)) return 4;
    // In cpp_bin_float_50: within its estimate of 1 / (2 sqrt 2) (Boost's sqrt, correctly rounded in the type), and
    // the estimate near 1e-46, against about 1e-14 in double.
    bin50 const exact50 = bin50(1) / (2 * sqrt(bin50(2)));
    if (!d50 || !d50->abs_error || !(abs(bin50(d50->value - exact50)) <= *d50->abs_error) || !(*d50->abs_error < 1e-45))
        return 6;
    // "A stencil gains less ... about two thirds of its digits": the central difference in cpp_bin_float_50, with its
    // default step eps^(1/3), within 1e-30 relative but not within 1e-40 (model: eps^(2/3), about 1e-33).
    auto const c50 = dv::derivative(root, bin50(2));
    bin50 const c50_error = c50 ? bin50(abs(bin50(c50->value - exact50))) : bin50(1);
    if (!c50 || !(c50_error < 1e-30 * exact50 && c50_error > 1e-40 * exact50)) return 7;
    return 0;
}
