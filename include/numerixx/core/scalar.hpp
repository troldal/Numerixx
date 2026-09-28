// Scalars (DESIGN §6.1): the open trait nxx::scalar_traits and the concept nxx::real.
#pragma once

#include <concepts>
#include <limits>
#include <type_traits>

namespace nxx
{
    // Primary template: "not a scalar". Users specialise it (is_real = true, epsilon()) for a real type that has no
    // std::numeric_limits specialisation.
    template<class T>
    struct scalar_traits
    {
    };

    // Every type with a specialised numeric_limits that is neither integral nor exact: float, double, long double,
    // and multiprecision floating-point types such as cpp_bin_float_50. Rejects int, bool, cpp_int and cpp_rational.
    template<class T>
        requires(std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_integer && !std::numeric_limits<T>::is_exact)
    struct scalar_traits<T>
    {
        static constexpr bool is_real = true;
        static constexpr T    epsilon() noexcept { return std::numeric_limits<T>::epsilon(); }
    };

    // A bool variable template rather than a concept used in folds: clang-cl safe (DESIGN §5.3).
    template<class T>
    inline constexpr bool is_real_v = requires { requires scalar_traits<std::remove_cvref_t<T>>::is_real; };

    template<class T>
    concept real = is_real_v<T> && std::regular<T> && std::totally_ordered<T> && requires(const T a, const T b) {
        { a + b } -> std::convertible_to<T>;
        { a - b } -> std::convertible_to<T>;
        { a * b } -> std::convertible_to<T>;
        { a / b } -> std::convertible_to<T>;
        { -a } -> std::convertible_to<T>;
    };
}    // namespace nxx
