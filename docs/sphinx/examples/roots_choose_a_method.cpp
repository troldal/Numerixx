// The program of docs/sphinx/roots/choose-a-method.md: measuring the evaluations of fn per method. The page
// shows the parts between the markers and the output in the raw string; main returns non-zero unless the program
// prints that output, whose counts are the independent prototypes' (docs/roots/design.md, Evidence, Phases 1 and 4),
// valid for fn written as here and builds without FMA contraction.

// clang-format off
// <includes>
#include <numerixx/roots.hpp>

#include <print>
// </includes>

// <count>
// The number of calls of fn that solve<M> makes on bracket, or -1 where it finds no root.
template<typename M, typename F>
int evaluations(F const &fn, nxx::roots::Bracket<double> const &bracket)
{
    int calls = 0;
    auto const counted = [&](double x) {
        ++calls;
        return fn(x);
    };
    return nxx::roots::solve<M>(counted, bracket) ? calls : -1;
}
// </count>
// clang-format on

#include "printed.hpp"

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <compare>
    namespace rt = nxx::roots;
    auto const wallis = [](double x) { return (x * x - 2.0) * x - 5.0; };
    rt::Bracket const bracket{2.0, 3.0};
    int const bisection = evaluations<rt::Bisection>(wallis, bracket);
    int const brent = evaluations<rt::Brent>(wallis, bracket);
    int const toms748 = evaluations<rt::Toms748>(wallis, bracket);
    int const anderson_bjorck = evaluations<rt::AndersonBjorck>(wallis, bracket);
    std::println("Bisection {}, Brent {}, Toms748 {}, AndersonBjorck {}",
                 bisection, brent, toms748, anderson_bjorck);
    // </compare>
    // clang-format on

    constexpr std::string_view output = R"(
Bisection 51, Brent 8, Toms748 9, AndersonBjorck 9
)";
    if (!(docs::printed() == docs::quoted(output))) return 1;
    return 0;
}
