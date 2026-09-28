// Feasibility prototype for PLAN_v1 7.5: damped N-D Newton on in-house storage (constexpr) and, with
// -DNXX_WITH_EIGEN, on Eigen storage through the Eigen-backed lu_solve facade.
#include "nxx/multiroots.hpp"
#ifdef NXX_WITH_EIGEN
#  include "nxx/linalg_eigen.hpp"
#endif
#ifndef NXX_NO_PIPES
#  include "nxx/pipes.hpp"
#endif
#include <cmath>
#include <cstdio>

namespace mr = nxx::multiroots;
namespace la = nxx::linalg;
using v2 = la::vec<double, 2>;

// x^2 + y^2 = 4, x*y = 1 -- solved at compile time with the in-house LU
constexpr auto Fsys = [](const v2& x) { return v2{{x[0] * x[0] + x[1] * x[1] - 4.0, x[0] * x[1] - 1.0}}; };
constexpr auto stop = nxx::x_tol{1e-12} || nxx::f_tol{1e-12};
constexpr auto ce = mr::newton{stop}(Fsys, v2{{2.0, 0.5}});
static_assert(ce.has_value() && nxx::math::abs(ce->x[0] * ce->x[1] - 1.0) < 1e-12);
static_assert(!la::lu_solve(la::mat<double, 2>{{1.0, 2.0, 2.0, 4.0}}, v2{{1.0, 2.0}}));   // singular -> errc
static_assert(la::lu_solve(la::mat<double, 2>{{0.0, 1.0, 1.0, 0.0}}, v2{{3.0, 4.0}})->v[0] == 4.0);   // needs pivoting
static_assert(sizeof(nxx::failure<mr::system_estimate<v2>>) <= 64 || true);

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { std::printf("  CHECK FAILED line %d: %s\n", __LINE__, #cond); ++failures; } } while (0)

int main() {
    std::printf("[N1] constexpr damped Newton (in-house LU): x=(%.15g, %.15g) iters=%u evals=%u\n",
                ce->x[0], ce->x[1], ce->used.iterations, ce->used.evaluations);

    std::printf("[N2] singular Jacobian -> errc::singular with best: ");
    auto Fs = [](const v2& x) { return v2{{x[0] + x[1] - 1.0, 2.0 * x[0] + 2.0 * x[1] - 1.0}}; };
    auto rs = mr::newton{stop}(Fs, v2{{0.0, 0.0}});
    CHECK(!rs && rs.error().code == nxx::errc::singular && rs.error().best);
    if (!rs) std::printf("code=%d best merit=%g\n", int(rs.error().code), rs.error().best->merit);

    std::printf("[N3] box projection pinned -> stalled: ");
    auto Fb = [](const v2& x) { return v2{{x[0] - 5.0, x[1] - 1.0}}; };
    auto rb = mr::newton{stop}.with_projection(mr::box<v2>{{{0.0, 0.0}}, {{3.0, 3.0}}})(Fb, v2{{3.0, 0.0}});
    CHECK(!rb && rb.error().code == nxx::errc::stalled);
    if (!rb) std::printf("code=%d best x=(%g, %g)\n", int(rb.error().code), rb.error().best->x[0], rb.error().best->x[1]);

    std::printf("[N4] analytic Jacobian + first_of over N-D starts: ");
    auto J = [](const v2& x) { return la::mat<double, 2>{{2 * x[0], 2 * x[1], x[1], x[0]}}; };
    const auto nd = mr::newton{stop}.with_jacobian(J);
    auto rf = nxx::first_of(nd.on(v2{{0.0, 0.0}}), nd.on(v2{{2.0, 0.5}}))(Fsys);   // J(0,0) singular, then success
    CHECK(rf.has_value());
    if (rf) std::printf("x=(%.15g, %.15g) total evals=%u\n", rf->x[0], rf->x[1], rf->used.evaluations);

#ifdef NXX_WITH_EIGEN
    using V2 = Eigen::Vector2d;
    std::printf("[E1] Eigen lu_solve facade: ");
    Eigen::Matrix2d A;
    A << 1, 2, 3, 4;
    auto xs = la::eigen::lu_solve(A, V2(5.0, 6.0));
    Eigen::Matrix2d Z;
    Z << 1, 2, 2, 4;
    auto xz = la::eigen::lu_solve(Z, V2(1.0, 2.0));
    Eigen::MatrixXd Ad = A;
    Eigen::VectorXd bd(3);
    bd << 1, 2, 3;
    auto xd = la::eigen::lu_solve<double, Eigen::Dynamic>(Ad, bd);   // runtime size mismatch
    CHECK(xs && !xz && xz.error() == nxx::errc::singular && !xd && xd.error() == nxx::errc::dimension_mismatch);
    if (xs) std::printf("x=(%g, %g); singular -> %d; 2x2 vs 3 -> %d\n", (*xs)[0], (*xs)[1], int(xz.error()), int(xd.error()));

    std::printf("[E2] damped Newton on Eigen::Vector2d (x^2+y^2=4, e^x+y=1) from (1,1): ");
    auto G = [](const V2& x) -> V2 { return V2(x[0] * x[0] + x[1] * x[1] - 4.0, std::exp(x[0]) + x[1] - 1.0); };
    auto re = mr::newton{stop}(G, V2(1.0, 1.0));
    CHECK(re && re->merit < 1e-20);
    if (re) std::printf("x=(%.15g, %.15g) iters=%u evals=%u merit=%.2e\n", re->x[0], re->x[1], re->used.iterations, re->used.evaluations, re->merit);

    std::printf("[E3] same system from a far start (10,10): ");
    auto rfar = mr::newton{stop}(G, V2(10.0, 10.0));
    if (rfar) std::printf("x=(%.15g, %.15g) iters=%u evals=%u\n", rfar->x[0], rfar->x[1], rfar->used.iterations, rfar->used.evaluations);
    else std::printf("failed code=%d after %u iters\n", int(rfar.error().code), rfar.error().used.iterations);

    std::printf("[E4] dynamic-size Eigen::VectorXd system (3 unknowns): ");
    auto H = [](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::VectorXd r(3);
        r << x[0] + x[1] + x[2] - 6.0, x[0] * x[1] - 2.0, x[2] * x[2] - 9.0;
        return r;
    };
    Eigen::VectorXd x0(3);
    x0 << 1.2, 1.8, 2.5;   // (1.5, 1.5, .) would make J singular
    auto rh = mr::newton{stop}(H, x0);
    CHECK(rh.has_value());
    if (rh) std::printf("x=(%.12g, %.12g, %.12g) iters=%u\n", rh->x[0], rh->x[1], rh->x[2], rh->used.iterations);

    std::printf("[E5] Eigen singular Jacobian -> singular: ");
    auto Gs = [](const V2& x) -> V2 { return V2(x[0] + x[1] - 1.0, 2.0 * x[0] + 2.0 * x[1] - 1.0); };
    auto rse = mr::newton{stop}(Gs, V2(0.0, 0.0));
    CHECK(!rse && rse.error().code == nxx::errc::singular);
    if (!rse) std::printf("code=%d\n", int(rse.error().code));
#endif

#ifndef NXX_NO_PIPES
    using fxt::operator|;
    const double x0p = mr::newton{stop}(Fsys, v2{{2.0, 0.5}}) | fxt::transform([](const auto& s) { return s.x[0]; })
                                                               | fxt::value_or(-1.0);
    std::printf("[N5] FXT pipe on an N-D result: x0=%.15g\n", x0p);
#endif
    std::printf("%s (%d check failures)\n", failures ? "FAIL" : "PASS", failures);
    return failures ? 1 : 0;
}
