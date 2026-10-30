// The program of docs/sphinx/roots/other-types.md: float, long double and Boost.Multiprecision. The page shows the
// parts between the markers and the output in the raw string; main returns non-zero unless the program prints that
// output and what the page says holds, against sqrt(2) to 70 digits (mpmath).

// clang-format off
// <includes>
#include <numerixx/roots.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <limits>
#include <print>
#include <utility>
// </includes>
// clang-format on

#include "printed.hpp"

#include <cmath>
#include <expected>

namespace {

    namespace rt = nxx::roots;

    constexpr char const *sqrt2_70 = "1.414213562373095048801688724209698078569671875376948073176679737990732";

    // Within abs_error + c of sqrt(2), and abs_error <= 2 tol(value), in T. c = eps sqrt(2): x * x rounds by half an
    // ulp of 2 (eps) over |f'| = 2.83, plus half an ulp of sqrt(2) for a reference held in T (0.7 eps).
    template<typename T>
    bool accurate(std::expected<rt::Root<T>, rt::error> const &r, T const &reference)
    {
        using std::abs;
        if (!r) return false;
        T const eps = std::numeric_limits<T>::epsilon();
        T const bound = T(2) * T(T(std::numeric_limits<T>::denorm_min()) + T(T(2) * eps * abs(r->value)));
        return abs(T(r->value - reference)) <= T(r->abs_error + T(eps * reference)) && r->abs_error <= bound;
    }


    // The walk to 0 in cpp_bin_float_50: from a bracket of width about 1 down to denorm_min = 2^min_exponent (about),
    // one halving per bit, about 2e9 halvings.
    static_assert(std::numeric_limits<boost::multiprecision::cpp_bin_float_50>::min_exponent < -2000000000
                  && std::numeric_limits<boost::multiprecision::cpp_bin_float_50>::min_exponent > -2200000000);

    // A pair deduced from cpp_dec_float_50 arithmetic holds expressions, not numbers, and solve rejects it.
    template<typename T>
    constexpr bool solvable_with_deduced_pair = requires(T const &x0) {
        rt::solve<rt::Newton>([](T const &x) { return std::pair{x * x - 2, 2 * x}; }, x0);
    };
    static_assert(!solvable_with_deduced_pair<boost::multiprecision::cpp_dec_float_50>);
    static_assert(solvable_with_deduced_pair<double>);

    // A fn returning x * x - 2 compiles in cpp_bin_float_50, where that is a number, not in cpp_dec_float_50, where
    // it is an expression; nor does a separate derivative returning 2 * x there.
    template<typename T>
    constexpr bool solvable_returning_arithmetic = requires(T const &x0) {
        rt::solve<rt::Brent>([](T const &x) { return x * x - 2; }, rt::Bracket<T>{x0, x0});
    };
    template<typename T>
    constexpr bool solvable_with_derivative_returning_arithmetic = requires(T const &x0) {
        rt::solve<rt::Newton>([](T const &x) { return T(x * x - 2); }, [](T const &x) { return 2 * x; }, x0);
    };
    static_assert(solvable_returning_arithmetic<boost::multiprecision::cpp_bin_float_50>
                  && !solvable_returning_arithmetic<boost::multiprecision::cpp_dec_float_50>);
    static_assert(solvable_with_derivative_returning_arithmetic<boost::multiprecision::cpp_bin_float_50>
                  && !solvable_with_derivative_returning_arithmetic<boost::multiprecision::cpp_dec_float_50>);

    // cpp_dec_float_50 uses expression templates, cpp_bin_float_50 does not.
    static_assert(boost::multiprecision::cpp_dec_float_50::et == boost::multiprecision::et_on
                  && boost::multiprecision::cpp_bin_float_50::et == boost::multiprecision::et_off);

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <float>
    auto const in_float = rt::solve<rt::Brent>([](float x) { return x * x - 2.0f; }, rt::Bracket{1.0f, 2.0f});
    auto const in_long_double =
        rt::solve<rt::Brent>([](long double x) { return x * x - 2.0L; }, rt::Bracket{1.0L, 2.0L});
    // </float>

    // <multiprecision>
    using bin50 = boost::multiprecision::cpp_bin_float_50;
    auto const in_bin50 = rt::solve<rt::Brent>([](bin50 const &x) { return bin50(x * x - 2); },
                                               rt::Bracket{bin50(1), bin50(2)});
    if (in_bin50) std::println("{}", in_bin50->value.str(50));
    // </multiprecision>

    // <expression-templates>
    using dec50 = boost::multiprecision::cpp_dec_float_50; // expression templates on
    auto const fdf = [](dec50 const &x) { return std::pair<dec50, dec50>{x * x - 2, 2 * x}; };
    auto const in_dec50 = rt::solve<rt::Newton>(fdf, dec50(1));
    // </expression-templates>
    // clang-format on

    constexpr std::string_view output = R"(
1.4142135623730950488016887242096980785696718753769
)";
    if (!(docs::printed() == docs::quoted(output))) return 7;

    if (!accurate(in_float, 1.41421356f)) return 1;
    if (!accurate(in_long_double, 1.414213562373095048801688724209698078570L)) return 2;
    if (!accurate(in_bin50, bin50(sqrt2_70))) return 3;
    if (!accurate(in_dec50, dec50(sqrt2_70))) return 4;

    // A float fn with a double bracket compiles, but its root is only float-accurate while abs_error is a double's
    // bound: x * x - 2.0f rounds x to float, so the computed sign change lies at the float nearest sqrt(2) (or a
    // neighbour), 3.5e-8 from sqrt(2), far beyond abs_error <= 2 tol, about 9e-16.
    // The lambda rounds x to float as a float parameter would, without the conversion warning MSVC gives for one.
    auto const mixed = rt::solve<rt::Brent>(
        [](double x) {
            float const y = static_cast<float>(x);
            return y * y - 2.0f;
        },
        rt::Bracket{1.0, 2.0});
    if (!mixed || !(std::abs(mixed->value - 1.4142135623730950488) > 2e-8
                    && std::abs(mixed->value - 1.4142135623730950488) < 5e-8 && mixed->abs_error < 1e-15))
        return 5;

    // The open methods' iteration cap, 4 times the bits of T: 212 steps for a double. Newton from 0 on x^3 - 2x + 2
    // cycles between 0 and 1 and never brackets the root, so it evaluates fn at the guess and at each of 212 steps.
    int calls = 0;
    auto const cycle = [&](double x) {
        ++calls;
        return std::pair{(x * x - 2.0) * x + 2.0, 3.0 * x * x - 2.0};
    };
    if (!(rt::solve<rt::Newton>(cycle, 0.0) == std::unexpected(rt::error::not_converged))) return 6;
    if (!(calls == 1 + 4 * std::numeric_limits<double>::digits && calls == 213)) return 8;
    return 0;
}
