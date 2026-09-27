// Prototype of the Eigen-backed linalg facade: lu_solve returning std::expected, plus vector_traits for
// Eigen column vectors so N-D solvers run on Eigen storage. Concrete return types only (no `auto` expressions).
#pragma once
#include "linalg.hpp"
#include <Eigen/Core>
#include <Eigen/LU>

namespace nxx::linalg::eigen {

template<class T, int N>
auto lu_solve(const Eigen::Matrix<T, N, N>& A, const Eigen::Matrix<T, N, 1>& b) -> std::expected<Eigen::Matrix<T, N, 1>, errc> {
    if (A.rows() != A.cols() || A.rows() != b.size()) return std::unexpected(errc::dimension_mismatch);
    if (!A.allFinite() || !b.allFinite()) return std::unexpected(errc::non_finite_input);
    const Eigen::PartialPivLU<Eigen::Matrix<T, N, N>> lu(A);
    const T rc = lu.rcond();   // PartialPivLU never reports singularity itself: check rcond and the result
    if (!(rc > T(A.rows()) * std::numeric_limits<T>::epsilon())) return std::unexpected(errc::singular);
    Eigen::Matrix<T, N, 1> x = lu.solve(b);
    if (!x.allFinite()) return std::unexpected(errc::singular);
    return x;
}

}    // namespace nxx::linalg::eigen

template<class T, int N> struct nxx::linalg::vector_traits<Eigen::Matrix<T, N, 1>> {
    using scalar = T;
    using vector = Eigen::Matrix<T, N, 1>;
    using matrix = Eigen::Matrix<T, N, N>;
    static std::size_t size(const vector& v) noexcept { return static_cast<std::size_t>(v.size()); }
    static matrix zero_matrix(const vector& v) { return matrix::Zero(v.size(), v.size()); }
    static T norm_inf(const vector& v) noexcept { return v.size() ? v.cwiseAbs().maxCoeff() : T(0); }
    static T sq_norm(const vector& v) noexcept { return v.squaredNorm(); }
    static bool all_finite(const vector& v) noexcept { return v.allFinite(); }
    static auto solve(const matrix& J, const vector& b) { return nxx::linalg::eigen::lu_solve(J, b); }
};
