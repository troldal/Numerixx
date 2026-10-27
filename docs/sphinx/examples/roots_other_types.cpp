// The program of docs/sphinx/how-to/roots-other-types.md: float, long double, Boost.Multiprecision and constant
// expressions. The page shows the parts between the markers; main returns non-zero unless what the page says holds,
// against sqrt(2) to 70 digits (mpmath).

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

#include <cmath>
#include <expected>

namespace {

    namespace rt = nxx::roots;

    constexpr char const *sqrt2_70 = "1.414213562373095048801688724209698078569671875376948073176679737990732";

    // clang-format off
    // <constexpr>
    constexpr auto square = [](double x) { return x * x - 2.0; };
    constexpr auto root = rt::solve<rt::Brent>(square, rt::Bracket{1.0, 2.0});
    static_assert(root.has_value());
    // </constexpr>
    // clang-format on

    // long double's precision as the page states it: double's with MSVC and clang-cl, x87 extended with GCC and Clang
    // on x86-64, quadruple with Emscripten.
#if defined(_MSC_VER)
    static_assert(std::numeric_limits<long double>::digits == 53);
#elif defined(__EMSCRIPTEN__)
    static_assert(std::numeric_limits<long double>::digits == 113);
#elif defined(__x86_64__)
    static_assert(std::numeric_limits<long double>::digits == 64);
#endif

    // Within abs_error + c of sqrt(2), and abs_error <= 2 tol(value), in T. c = eps sqrt(2): x * x rounds by half an
    // ulp of 2 (eps) over |f'| = 2.83, plus half an ulp of sqrt(2) for a reference held in T (0.7 eps).
    template<typename T>
    bool accurate(std::expected<rt::Root<T>, rt::error> const &r, T const &reference)
    {
        using std::abs;
        T const eps = std::numeric_limits<T>::epsilon();
        T const bound = T(2) * T(T(std::numeric_limits<T>::denorm_min()) + T(T(2) * eps * abs(r->value)));
        return r.has_value() && abs(T(r->value - reference)) <= T(r->abs_error + T(eps * reference))
               && r->abs_error <= bound;
    }

} // namespace

int main()
{
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

    if (!accurate(in_float, 1.41421356f)) return 1;
    if (!accurate(in_long_double, 1.414213562373095048801688724209698078570L)) return 2;
    if (!accurate(in_bin50, bin50(sqrt2_70))) return 3;
    if (!accurate(in_dec50, dec50(sqrt2_70))) return 4;
    // The constant expression gives what a run gives (on builds without FMA contraction), and sqrt(2).
    if (!(root == rt::solve<rt::Brent>(square, rt::Bracket{1.0, 2.0}))) return 5;
    if (!accurate(root, 1.4142135623730950488)) return 6;
    return 0;
}
