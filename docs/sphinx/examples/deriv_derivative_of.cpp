// The program of docs/sphinx/deriv/derivative-of.md: a derivative as a callable, called at several points, and its copy
// of fn. The page shows the parts between the markers and the output in the raw strings; main returns non-zero unless
// the program prints that output and what the page says holds, against the closed form (sqrt x)' = 1 / (2 sqrt x).
// std::sqrt is correctly rounded (IEEE 754), so the output is the same on every target that does not contract
// a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>

#include <cmath>
#include <functional>
#include <initializer_list>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <format>
#include <string>

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;
    namespace dv = nxx::deriv;

    // clang-format off
    // <callable>
    auto const fn = [](double x) { return std::sqrt(x); };
    auto const slope = dv::derivative_of<dv::Ridders>(fn); // a function of x
    for (double const x : {1.0, 2.0, 4.0}) {
        auto const d = slope(x); // a std::expected, as from derivative<Ridders>(fn, x)
        if (d) std::println("f'({}) = {:.15f}, abs_error {:.1e}", x, d->value, *d->abs_error);
    }
    // </callable>

    // <state>
    // A function object that counts its own calls. derivative_of keeps a copy, whose count is its own;
    // std::ref(counter) passes the original instead.
    struct Counter
    {
        int calls = 0;
        double operator()(double x)
        {
            ++calls;
            return std::sqrt(x);
        }
    };
    Counter counter;
    auto copied = dv::derivative_of(counter);
    auto shared = dv::derivative_of(std::ref(counter));
    (void)copied(2.0);
    (void)shared(2.0);
    std::println("the original counted {} calls", counter.calls);
    // </state>
    // clang-format on

    constexpr std::string_view callable_output = R"(
f'(1) = 0.499999999999952, abs_error 1.1e-13
f'(2) = 0.353553390593294, abs_error 8.0e-14
f'(4) = 0.249999999999976, abs_error 5.7e-14
)";
    constexpr std::string_view state_output = R"(
the original counted 2 calls
)";
    if (!(docs::printed() == std::string(docs::quoted(callable_output)) + std::string(docs::quoted(state_output))))
        return 1;

    // "The exact values are 0.5, 0.353553390593274 and 0.25": 1 / (2 sqrt 2) to 15 decimals.
    if (!(std::format("{:.15f}", std::sqrt(2.0) / 4.0) == "0.353553390593274")) return 5;
    // The values: within their estimates of 1 / (2 sqrt x), and bit for bit what derivative<Ridders> gives.
    for (double const x : {1.0, 2.0, 4.0}) {
        auto const d = slope(x);
        if (!d || !d->abs_error || !(std::abs(d->value - 0.5 / std::sqrt(x)) <= *d->abs_error + 1.2e-16)) return 2;
        if (!(d == dv::derivative<dv::Ridders>(fn, x))) return 3;
    }
    // The copy: Central<2> calls fn twice for f', so the copy counted 2 and the original only the 2 calls made through
    // std::ref.
    if (!(counter.calls == 2)) return 4;
    return 0;
}
