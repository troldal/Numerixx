// The program of docs/sphinx/how-to/roots-find-a-bracket.md: search for a bracket from a guess, then solve in it. The
// page shows the parts between the markers; main returns non-zero unless what the page says holds. The count and the
// bracket are the independent prototype's (.scratch/roots5/proto5d.py, as in tests/test_roots_search.cpp), for fn
// written as here and builds without FMA contraction.

// clang-format off
// <includes>
#include <numerixx/roots.hpp>

#include <print>
// </includes>
// clang-format on

#include <cmath>
#include <limits>

namespace {

    constexpr auto wallis_fn = [](double x) { return (x * x - 2.0) * x - 5.0; };

    // search takes T from its guess, so search(fn, 0) asks for an int and fails the constraint real<T>.
    template<typename G>
    constexpr bool searchable = requires(G const &guess) { nxx::roots::search(wallis_fn, guess); };
    static_assert(searchable<double> && !searchable<int>);

} // namespace

int main()
{
    // clang-format off
    // <search>
    namespace rt = nxx::roots;
    auto const wallis = [](double x) { return (x * x - 2.0) * x - 5.0; };
    auto const bracket = rt::search(wallis, 0.0);
    if (bracket) std::println("[{}, {}]", bracket->lower, bracket->upper);
    // </search>

    // <chain>
    auto const root = rt::search(wallis, 0.0).and_then([&](rt::Bracket<double> const &b) {
        return rt::solve<rt::Brent>(wallis, b);
    });
    if (root) std::println("{:.17g} +- {:.1e}", root->value, root->abs_error);
    else std::println("no root: {}", root.error());
    // </chain>
    // clang-format on

    if (!bracket || !root) return 1;

    // The prototype's bracket and its 27 evaluations.
    int calls = 0;
    auto const counted = rt::search(
        [&](double x) {
            ++calls;
            return wallis(x);
        },
        0.0);
    if (!(counted == rt::Bracket{0x1.bade147ae147cp+0, 0x1.4c268f5c28f5dp+1} && calls == 27)) return 2;
    if (!(*bracket == *counted)) return 3;

    // The root 2.0945514815423265915 (mpmath): within abs_error + c, c = 2e-16 for fn's rounding (x * x, - 2, * x and
    // - 5 near the root, over |f'| = 11.2), plus half an ulp of the root for the reference in a double, 2.3e-16.
    double const reference = 2.0945514815423265915;
    double const eps = std::numeric_limits<double>::epsilon();
    if (!(std::abs(root->value - reference) <= root->abs_error + 2e-16 + 2.3e-16)) return 4;
    if (!(root->abs_error <= 2 * (std::numeric_limits<double>::denorm_min() + 2 * eps * root->value))) return 5;
    return 0;
}
