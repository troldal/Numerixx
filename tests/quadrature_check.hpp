// What the nxx::quadrature tests share: the types, the references, the error measurement and an evaluation counter.
//
// References are closed forms at 70 digits (mpmath at 80), never Numerixx. Every tolerance comes from an error
// estimate stated where it is used: for GaussLegendre the truncation bound and the per-node rounding bound of
// docs/quadrature/PHASE1_API.md D, for Romberg the bound of Bauer, Rutishauser and Stiefel plus its rounding (C).

#pragma once

#include <numerixx/quadrature.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <ios>
#include <limits>
#include <type_traits>

namespace quadrature_check {

    namespace mp = boost::multiprecision;
    namespace qd = nxx::quadrature;

    using bin50 = mp::cpp_bin_float_50;
    using bin50_et = mp::number<mp::cpp_bin_float<50>, mp::et_on>;
    using dec50 = mp::cpp_dec_float_50;
    // The type in which errors are measured: it holds every value of the types tested, and their differences, exactly
    // enough (100 digits against cpp_dec_float_50's at most about 60).
    using bin100 = mp::cpp_bin_float_100;

    inline constexpr double dmax = std::numeric_limits<double>::max();
    inline constexpr double dmin = std::numeric_limits<double>::denorm_min();
    inline constexpr double deps = std::numeric_limits<double>::epsilon();
    inline constexpr double dnan = std::numeric_limits<double>::quiet_NaN();
    inline constexpr double dinf = std::numeric_limits<double>::infinity();

    // 70-digit references (mpmath at 80 digits).
    inline constexpr char const *e_minus_1 = "1.718281828459045235360287471352662497757247093699959574966967627724077";
    inline constexpr char const *pi_over_4 = "0.7853981633974483096156608458198757210492923498437764552437361480769541";

    // A value of a tested type in bin100. cpp_dec_float_50 goes through its decimal string, exact at 70 significant
    // digits (more than it carries).
    template<typename T>
    bin100 wide(T const &v)
    {
        if constexpr (std::is_same_v<T, dec50>) return bin100(v.str(70, std::ios_base::scientific));
        else return bin100(v);
    }

    // |value - ref| in double, measured in bin100 (a difference formed in T would round).
    template<typename T>
    double error_of(T const &value, bin100 const &ref)
    { return static_cast<double>(bin100(abs(bin100(wide(value) - ref)))); }

    template<typename T>
    double error_of(T const &value, char const *ref)
    { return error_of(value, bin100(ref)); }

    // T's epsilon in double (cpp_bin_float_50's, 2^-167, is a double too).
    template<typename T>
    double eps_of()
    { return static_cast<double>(wide(T(std::numeric_limits<T>::epsilon()))); }

    // fn wrapped to count its evaluations.
    template<typename F>
    struct Counted
    {
        F fn;
        int *count;
        template<typename T>
        auto operator()(T const &x) const
        {
            ++*count;
            return fn(x);
        }
    };

    constexpr auto square = [](double x) { return x * x; };
    constexpr auto cube = [](double x) { return x * x * x; };

} // namespace quadrature_check
