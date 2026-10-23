// TEMPLATE example, used by the layout test. Every listing in the book is a file like this one: compiled and run
// as a test (CMakeLists.txt in the book folder), and quoted between snippet markers with \snippet{<file>}{<name>}.
// A marker is a comment line "// <name>" before the snippet and "// </name>" after it; the markers do not appear in
// the book. Snippet code and markers start in column 0 (listings cannot strip indentation from an input file):
// either a function or declaration at namespace scope, or statements inside a function between
// "// clang-format off" and "// clang-format on". The checks in main make the example fail its test if the library
// stops doing what the text says.

#include <numerixx/deriv.hpp>

#include <cmath>
#include <cstdio>

int main()
{
    // clang-format off
// <quickstart>
auto const f = [](double x) { return std::exp(x); };
auto const d = nxx::deriv::derivative<nxx::deriv::Ridders>(f, 1.0);
if (d) std::printf("f'(1) = %.15f +- %.1e\n", d->value, *d->abs_error);
// </quickstart>
    // clang-format on
    return d && std::abs(d->value - std::exp(1.0)) <= *d->abs_error ? 0 : 1;
}
