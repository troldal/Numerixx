//
// Created by Kenneth Balslev on 20/03/2023.
//

#ifndef NUMERIXX_CONCEPTS_HPP
#define NUMERIXX_CONCEPTS_HPP

#include <complex>

namespace nxx {

//    template <typename T>
//    struct is_complex : std::false_type {};
//
//    template <std::floating_point T>
//    struct is_complex<std::complex<T>> : std::true_type {};

//    template <typename T>
//    concept IsComplex = is_complex<T>::value;

    template <typename T>
    concept IsComplex = std::same_as<T, std::complex<typename T::value_type>>;

}

#endif    // NUMERIXX_CONCEPTS_HPP
