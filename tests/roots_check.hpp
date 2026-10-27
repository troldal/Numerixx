// What the nxx::roots tests share: the types, the references, the accuracy checks and an evaluation counter.
//
// A result must lie within abs_error of the reference plus a conditioning term c: the shift of fn's sign change from
// the exact root, fn's error budget near the root divided by |f'(root)|, derived per case (an fn of + - * /: its own
// rounding; a libm fn: a budget of 2 ulps). References are closed forms and published constants, cross-checked with
// mpmath, never Numerixx.

#pragma once

#include <numerixx/roots.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <cstdlib>
#include <expected>
#include <limits>
#include <type_traits>
#include <utility>

namespace roots_check {

    namespace mp = boost::multiprecision;
    namespace rt = nxx::roots;

    using bin50 = mp::cpp_bin_float_50;
    using bin50_et = mp::number<mp::cpp_bin_float<50>, mp::et_on>;
    using dec50 = mp::cpp_dec_float_50;

    inline constexpr double dmax = std::numeric_limits<double>::max();
    inline constexpr double dmin = std::numeric_limits<double>::denorm_min();
    inline constexpr double deps = std::numeric_limits<double>::epsilon();
    inline constexpr double dnan = std::numeric_limits<double>::quiet_NaN();
    inline constexpr double dinf = std::numeric_limits<double>::infinity();

    // 70-digit references (mpmath at 80 digits): more digits than any type tested (cpp_bin_float_50 carries about
    // 50.6 decimal digits, so a 50-digit reference would err by up to 10 of its epsilons).
    inline constexpr char const *sqrt2_70 = "1.414213562373095048801688724209698078569671875376948073176679737990732";
    inline constexpr char const *wallis_70 =
        "2.094551481542326591482386540579302963857306105628239180304128529045312"; // x^3 - 2x - 5
    inline constexpr double sqrt2 = 1.4142135623730950488;
    inline constexpr double wallis = 2.0945514815423265915;
    inline constexpr double dottie = 0.73908513321516064166; // cos x = x, OEIS A003957
    inline constexpr double omega = 0.56714329040978387300; // x e^x = 1, W(1), OEIS A030178
    inline constexpr double ln2 = 0.69314718055994530942;

    // A reference in T: a built-in type rounds the decimal string, a multiprecision type parses it.
    template<typename T>
    T from_string(char const *digits)
    {
        if constexpr (std::is_floating_point_v<T>)
            return static_cast<T>(std::strtold(digits, nullptr));
        else
            return T(digits);
    }

    // tol(x) with the default criteria, and the bound every result keeps: abs_error <= 2 tol(value).
    constexpr double default_tol(double x) { return dmin + 2.0 * deps * nxx::math::abs(x); }

    // |value - ref| <= abs_error + c, and abs_error <= 2 tol(value): the two checks of every accuracy case.
    template<typename T>
    bool brackets(std::expected<rt::Root<T>, rt::error> const &r, T const &ref, T const &c, T const &tol2)
    {
        using std::abs;
        return r.has_value() && abs(T(r->value - ref)) <= T(r->abs_error + c) && r->abs_error <= tol2;
    }

    inline bool brackets(std::expected<rt::Root<double>, rt::error> const &r, double ref, double c)
    { return r.has_value() && brackets(r, ref, c, 2.0 * default_tol(r->value)); }

    // The same check in cpp_bin_float_50, which holds every double and the difference of any two exactly: for
    // references and distances that double cannot represent (2.5 denorm_min; 2^53 + 0.5).
    inline bool brackets_exactly(std::expected<rt::Root<double>, rt::error> const &r, bin50 const &ref)
    { return r.has_value() && bin50(abs(bin50(bin50(r->value) - ref))) <= bin50(r->abs_error); }

    // fn wrapped to count its evaluations; it returns what fn returns, a value or a tuple of values.
    template<typename F>
    struct Counted
    {
        F fn;
        int *count;
        auto operator()(double x) const
        {
            ++*count;
            return fn(x);
        }
    };

    // solve<M> on fn from a start (a Bracket or a guess), with the number of evaluations of fn.
    template<typename M, typename F, typename S>
    std::pair<std::expected<rt::Root<double>, rt::error>, int> counted(F fn, S start, rt::Criteria<double> c = {})
    {
        int n = 0;
        auto const r = rt::solve<M>(Counted<F>{ fn, &n }, start, c);
        return { r, n };
    }

} // namespace roots_check
