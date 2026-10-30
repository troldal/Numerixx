// The program of docs/sphinx/deriv/accuracy.md: evaluations and relative errors of every method for f', with the
// default step, on three functions with closed-form derivatives: 1/x at 2 (-1/4), sqrt x at 2 (1 / (2 sqrt 2)) and the
// peak x / (x^2 + a^2), a = 1e-4, at 2e-4 ((a^2 - x^2) / (x^2 + a^2)^2). The page shows the output in the raw string;
// main returns non-zero unless the program prints it and the Ridders methods' estimates cover their errors. The
// functions use only + - * / and std::sqrt (correctly rounded), so the output is the same on every target that does not
// contract a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly).

#include <numerixx/deriv.hpp>

#include "printed.hpp"

#include <cmath>
#include <format>
#include <print>
#include <string>

namespace {

    namespace dv = nxx::deriv;

    constexpr double width = 1e-4; // the peak's, a in the page
    constexpr double peak_x = 2e-4;

    double inverse(double x) { return 1.0 / x; }
    double root(double x) { return std::sqrt(x); }
    double peak(double x) { return x / (x * x + width * width); }

    bool covered = true; // whether every estimate covered its error

    // One cell: the calls of fn and the relative error of derivative<M>(fn, x).
    template<typename M>
    std::string cell(double (*fn)(double), double x, double exact)
    {
        int calls = 0;
        auto const counted = [&](double t) {
            ++calls;
            return fn(t);
        };
        auto const d = dv::derivative<M>(counted, x);
        if (!d) return "failed";
        if (d->abs_error && !(std::abs(d->value - exact) <= *d->abs_error)) covered = false;
        return std::format("{:>3} {:.0e}", calls, std::abs(d->value - exact) / std::abs(exact));
    }

    template<typename M>
    void row(char const *name)
    {
        double const peak_exact = (width * width - peak_x * peak_x)
                                  / ((peak_x * peak_x + width * width) * (peak_x * peak_x + width * width));
        std::println("{:<16} {:>10} {:>10} {:>10}", name, cell<M>(inverse, 2.0, -0.25),
                     cell<M>(root, 2.0, std::sqrt(2.0) / 4.0), cell<M>(peak, peak_x, peak_exact));
    }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    std::println("{:<16} {:>10} {:>10} {:>10}", "method", "1/x", "sqrt x", "peak");
    row<dv::Central<2>>("Central<2>");
    row<dv::Central<4>>("Central<4>");
    row<dv::Central<6>>("Central<6>");
    row<dv::Central<10>>("Central<10>");
    row<dv::Forward<2>>("Forward<2>");
    row<dv::Backward<2>>("Backward<2>");
    row<dv::Ridders>("Ridders");
    row<dv::ForwardRidders>("ForwardRidders");
    row<dv::BackwardRidders>("BackwardRidders");

    constexpr std::string_view output = R"(
method                  1/x     sqrt x       peak
Central<2>          2 4e-11    2 6e-12    2 3e-04
Central<4>          4 1e-12    4 2e-14    4 1e+00
Central<6>          6 1e-12    6 7e-14    6 1e+00
Central<10>        10 9e-11   10 1e-12   10 1e+00
Forward<2>          3 1e-10    3 2e-10    3 6e-04
Backward<2>         3 5e-11    3 1e-10    3 8e-04
Ridders            14 2e-15   14 6e-14   36 3e-14
ForwardRidders     10 2e-13    9 2e-13   21 2e-13
BackwardRidders    10 2e-13    9 1e-13   21 8e-13
)";
    if (!(docs::printed() == docs::quoted(output))) return 1;
    // The Ridders methods' estimates cover their actual errors on these functions.
    if (!covered) return 2;
    return 0;
}
