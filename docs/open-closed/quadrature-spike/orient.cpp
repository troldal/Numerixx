// Orientation with today's header: is integrate(fn, b, a) exactly -integrate(fn, a, b), the whole result (status,
// value, abs_error, best), for each method? Intervals of both signs, across 0, and wide; integrands that cancel.
#include <numerixx/quadrature.hpp>
#include <cmath>
#include <cstdio>
namespace qd = nxx::quadrature;

template<typename T>
bool negated(std::expected<qd::Integral<T>, qd::Failure<T>> const &r, std::expected<qd::Integral<T>, qd::Failure<T>> const &s)
{
    auto const same = [](qd::Integral<T> const &x, qd::Integral<T> const &y) { return x.value == -y.value && x.abs_error == y.abs_error; };
    if (r.has_value() != s.has_value()) return false;
    if (r) return same(*r, *s);
    if (r.error().reason != s.error().reason || r.error().best.has_value() != s.error().best.has_value()) return false;
    return !r.error().best || same(*r.error().best, *s.error().best);
}

int main()
{
    int n = 0, gl = 0, ro = 0, ro_tight = 0, gk = 0, gk0 = 0;
    for (int i = 1; i <= 3000; ++i) {
        double const a = (i % 4 == 0) ? -0.013 * i : (i % 4 == 1) ? 0.013 * i : (i % 4 == 2) ? -0.5 - 0.001 * i : -1e3 * i;
        double const b = (i % 4 == 2) ? 0.3 + 0.002 * i : (i % 4 == 3) ? 1e3 * i * 0.7 : a + 0.7 + 0.01 * i;
        auto f = [i](double x) {
            return i % 3 == 0 ? std::sin(37.0 * i * x) * 1e3 - x
                 : i % 3 == 1 ? std::exp(std::sin(i * x)) / (1 + x * x)
                              : std::cos(i * x) - 0.999;
        };
        ++n;
        auto const g1 = qd::integrate<qd::GaussLegendre<64>>(f, a, b), g2 = qd::integrate<qd::GaussLegendre<64>>(f, b, a);
        if (g1.has_value() != g2.has_value() || (g1 && *g1 != -*g2)) ++gl;
        if (!negated(qd::integrate<qd::Romberg>(f, a, b), qd::integrate<qd::Romberg>(f, b, a))) ++ro;
        if (!negated(qd::integrate<qd::Romberg>(f, a, b, { .rel_tol = 0 }), qd::integrate<qd::Romberg>(f, b, a, { .rel_tol = 0 }))) ++ro_tight;
        if (!negated(qd::integrate(f, a, b), qd::integrate(f, b, a))) ++gk;
        if (!negated(qd::integrate(f, a, b, { .rel_tol = 0 }), qd::integrate(f, b, a, { .rel_tol = 0 }))) ++gk0;
    }
    std::printf("of %d intervals, results that are not the exact negation when reversed:\n"
                "GaussLegendre<64> %d, Romberg %d (rel_tol 0: %d), GaussKronrod<15> %d (rel_tol 0: %d)\n",
        n, gl, ro, ro_tight, gk, gk0);
}
