// Prototype of PLAN_v1 7.5 linalg: a tiny in-house fixed-size vec/mat + pivoted LU (constexpr, zero heap),
// plus the vector_traits customisation point every storage (in-house, Eigen) must provide.
#pragma once
#include "core.hpp"

namespace nxx::linalg {

// ---------------------------------------------------------------- the adapter surface N-D solvers use
template<class V> struct vector_traits;   // scalar, matrix, size, norm_inf, sq_norm, all_finite, solve, zero_matrix
template<class V> concept dense_vector = requires { typename vector_traits<V>::scalar; typename vector_traits<V>::matrix; };

// ---------------------------------------------------------------- in-house storage
template<real T, std::size_t N> struct vec {
    std::array<T, N> v{};
    constexpr T& operator[](std::size_t i) noexcept { return v[i]; }
    constexpr const T& operator[](std::size_t i) const noexcept { return v[i]; }
    static constexpr std::size_t size() noexcept { return N; }
    friend constexpr vec operator+(vec a, const vec& b) noexcept { for (std::size_t i = 0; i < N; ++i) a[i] += b[i]; return a; }
    friend constexpr vec operator-(vec a, const vec& b) noexcept { for (std::size_t i = 0; i < N; ++i) a[i] -= b[i]; return a; }
    friend constexpr vec operator-(vec a) noexcept { for (auto& x : a.v) x = -x; return a; }
    friend constexpr vec operator*(const T& s, vec a) noexcept { for (auto& x : a.v) x *= s; return a; }
    friend constexpr bool operator==(const vec&, const vec&) = default;
};
template<real T, std::size_t N> struct mat {
    std::array<T, N * N> a{};   // row-major
    constexpr T& operator()(std::size_t i, std::size_t j) noexcept { return a[i * N + j]; }
    constexpr const T& operator()(std::size_t i, std::size_t j) const noexcept { return a[i * N + j]; }
};

// pivoted LU solve with a RELATIVE pivot threshold n*eps*max|a|; singular -> errc, never NaN
template<real T, std::size_t N>
constexpr auto lu_solve(mat<T, N> A, vec<T, N> b) noexcept -> std::expected<vec<T, N>, errc> {
    T amax(0);
    for (const T& x : A.a) { if (!math::isfinite(x)) return std::unexpected(errc::non_finite_input); amax = math::max(amax, math::abs(x)); }
    const T thresh = T(N) * std::numeric_limits<T>::epsilon() * amax;
    for (std::size_t k = 0; k < N; ++k) {
        std::size_t piv = k;
        for (std::size_t i = k + 1; i < N; ++i) if (math::abs(A(i, k)) > math::abs(A(piv, k))) piv = i;
        if (!(math::abs(A(piv, k)) > thresh)) return std::unexpected(errc::singular);
        if (piv != k) { for (std::size_t j = 0; j < N; ++j) std::swap(A(k, j), A(piv, j)); std::swap(b[k], b[piv]); }
        for (std::size_t i = k + 1; i < N; ++i) {
            const T l = A(i, k) / A(k, k);
            for (std::size_t j = k; j < N; ++j) A(i, j) -= l * A(k, j);
            b[i] -= l * b[k];
        }
    }
    vec<T, N> x{};
    for (std::size_t i = N; i-- > 0;) {
        T s = b[i];
        for (std::size_t j = i + 1; j < N; ++j) s -= A(i, j) * x[j];
        x[i] = s / A(i, i);
    }
    return x;
}

template<real T, std::size_t N> struct vector_traits<vec<T, N>> {
    using scalar = T;
    using vector = vec<T, N>;
    using matrix = mat<T, N>;
    static constexpr std::size_t size(const vector&) noexcept { return N; }
    static constexpr matrix zero_matrix(const vector&) noexcept { return {}; }
    static constexpr T norm_inf(const vector& v) noexcept { T m(0); for (const T& x : v.v) m = math::max(m, math::abs(x)); return m; }
    static constexpr T sq_norm(const vector& v) noexcept { T s(0); for (const T& x : v.v) s += x * x; return s; }
    static constexpr bool all_finite(const vector& v) noexcept { for (const T& x : v.v) if (!math::isfinite(x)) return false; return true; }
    static constexpr auto solve(const matrix& J, const vector& b) noexcept { return lu_solve(J, b); }
};

}    // namespace nxx::linalg
