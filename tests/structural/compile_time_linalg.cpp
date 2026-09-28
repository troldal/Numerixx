// Compile-time record (DESIGN §3.6): a linalg TU that instantiates Eigen solves, fixed-size and dynamic, as the
// multiroots solvers will. structural.compile_time.linalg rebuilds this TU and records its time; it is not gated.
#include <numerixx/linalg.hpp>

#include <Eigen/LU>

namespace compile_time
{
    Eigen::Vector3d solve_fixed(const Eigen::Matrix3d& a, const Eigen::Vector3d& b);
    Eigen::VectorXd solve_dynamic(const Eigen::MatrixXd& a, const Eigen::VectorXd& b);

    Eigen::Vector3d solve_fixed(const Eigen::Matrix3d& a, const Eigen::Vector3d& b) { return a.partialPivLu().solve(b); }

    Eigen::VectorXd solve_dynamic(const Eigen::MatrixXd& a, const Eigen::VectorXd& b) { return a.partialPivLu().solve(b); }
}    // namespace compile_time
