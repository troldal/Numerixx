// The program of docs/sphinx/basics.md, its sections on real types and compile time: a solve in a Boost.Multiprecision
// type with expression templates, long double's precision by platform, the concept nxx::real, and a constant
// expression. The page shows the parts between the markers and the output in the raw string; main returns non-zero
// unless the program prints that output and what the page says holds, against the root of cos x = x to 80 digits
// (mpmath).

// clang-format off
// <includes>
#include <numerixx/roots.hpp>

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <cmath>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <numerixx/deriv.hpp>

#include <limits>

namespace {

    namespace rt = nxx::roots;

    // An int is not a real type, and a bracket of ints does not compile: the constraint real<T> fails.
    template<typename T>
    constexpr bool bracketable = requires(T const a) { rt::Bracket{a, a}; };
    static_assert(nxx::real<float> && nxx::real<long double> && !nxx::real<int>);
    static_assert(bracketable<double> && !bracketable<int>);
    // Nor does a derivative at an int point: derivative(fn, 2) deduces T = int.
    template<typename T>
    constexpr bool differentiable = requires(T const x) { nxx::deriv::derivative([](double t) { return t; }, x); };
    static_assert(differentiable<double> && !differentiable<int>);

    // long double's precision as the page states it: double's with MSVC and clang-cl, x87 extended (a 64-bit
    // significand) with GCC and Clang on x86-64, quadruple (113 bits) with Emscripten.
#if defined(_MSC_VER)
    static_assert(std::numeric_limits<long double>::digits == 53);
#elif defined(__EMSCRIPTEN__)
    static_assert(std::numeric_limits<long double>::digits == 113);
#elif defined(__x86_64__)
    static_assert(std::numeric_limits<long double>::digits == 64);
#endif

    // clang-format off
    // <constexpr>
    constexpr auto square = [](double x) { return x * x - 2.0; };
    constexpr auto root = rt::solve<rt::Brent>(square, rt::Bracket{1.0, 2.0});
    static_assert(root && root->value > 1.414 && root->value < 1.415);
    // </constexpr>
    // clang-format on

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <multiprecision>
    using dec50 = boost::multiprecision::cpp_dec_float_50; // 50 decimal digits, expression templates on
    auto const fn = [](dec50 const &x) {
        using std::cos;           // so that cos(x) finds Boost's cos for dec50
        return dec50(cos(x) - x); // return the type, not the expression
    };
    auto const r = rt::solve<rt::Brent>(fn, rt::Bracket{dec50(0), dec50(1)});
    if (r) std::println("{}", r->value.str(50));
    // </multiprecision>
    // clang-format on

    constexpr std::string_view output = R"(
0.73908513321516064165531208767387340401341175890076
)";
    if (!(docs::printed() == docs::quoted(output))) return 1;

    // The root of cos x = x within abs_error + c, with abs_error <= 2 tol(value): c = eps for Boost's cos near the
    // root (accurate to about an ulp) over |f'| = 1 + sin x = 1.67, plus half an ulp of the reference held in dec50.
    dec50 const reference("0.73908513321516064165531208767387340401341175890075746496568063577328465488354759");
    dec50 const eps = std::numeric_limits<dec50>::epsilon();
    if (!r || !(abs(dec50(r->value - reference)) <= dec50(r->abs_error + eps))) return 2;
    if (!(r->abs_error <= dec50(2 * dec50(dec50(std::numeric_limits<dec50>::denorm_min()) + 2 * eps * r->value))))
        return 3;

    // The constant expression gives what a run gives (on builds without FMA contraction), and sqrt(2): x * x rounds by
    // half an ulp of 2 over |f'| = 2.83, 8e-17, plus half an ulp of sqrt(2) for the reference in a double, 1.2e-16.
    if (!(root == rt::solve<rt::Brent>(square, rt::Bracket{1.0, 2.0}))) return 4;
    if (!(std::abs(root->value - 1.4142135623730950488) <= root->abs_error + 8e-17 + 1.2e-16)) return 5;
    return 0;
}
