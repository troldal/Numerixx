#include <numerixx/quadrature.hpp>
#include <cmath>
#include <cstdio>
namespace qd = nxx::quadrature;
int main()
{
    int differ_gl = 0, differ_r = 0, differ_gk = 0, n = 0;
    for (int i = 1; i <= 3000; ++i) {
        double const a = 0.013 * i, b = a + 0.7 + 0.01 * i;
        auto f = [i](double x) { return (i % 3 == 0 ? std::sin(37.0 * i * x) * 1e3 - x : i % 3 == 1 ? std::exp(std::sin(i * x)) / (1 + x * x) : std::cos(i * x) - 0.999); };
        ++n;
        auto const gl = qd::integrate<qd::GaussLegendre<64>>(f, a, b), glr = qd::integrate<qd::GaussLegendre<64>>(f, b, a);
        if (*gl != -*glr) ++differ_gl;
        auto const r = qd::integrate<qd::Romberg>(f, a, b), rr = qd::integrate<qd::Romberg>(f, b, a);
        if (r && rr && r->value != -rr->value) ++differ_r;
        auto const k = qd::integrate(f, a, b), kr = qd::integrate(f, b, a);
        if (k && kr && k->value != -kr->value) ++differ_gk;
    }
    std::printf("of %d: integrate(b, a) != -integrate(a, b) for GaussLegendre<64> %d, Romberg %d, GaussKronrod %d\n", n, differ_gl, differ_r, differ_gk);
}
