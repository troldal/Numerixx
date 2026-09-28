// Boost.Multiprecision scalars in Eigen, for numerixx::linalg and numerixx::multiroots (DESIGN §3.5, §7.5).
// Pulls in Boost's Eigen glue (NumTraits for the multiprecision types). Do not stream Eigen matrices of
// cpp_bin_float values: Eigen's operator<< does not compile with them on Boost 1.92.
//
// Status: build skeleton (roadmap phase 0). Implemented in phase 9.
#pragma once

#include <numerixx/adapters/multiprecision.hpp>
#include <numerixx/linalg.hpp>

#include <Eigen/Core>
#include <boost/multiprecision/eigen.hpp>

namespace nxx::adapters
{
}    // namespace nxx::adapters
