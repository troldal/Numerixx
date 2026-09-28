// A consumer of every non-adapter module: the umbrella header, the FXT pipes and the Eigen-backed linalg.
#include <numerixx/linalg.hpp>
#include <numerixx/multiroots.hpp>
#include <numerixx/numerixx.hpp>
#include <numerixx/pipes.hpp>

#include <Eigen/Core>

#include <cstdio>
#include <expected>

int main()
{
    using nxx::operator|;
    const auto            doubled = std::expected<int, int> { nxx::version.major } | fxt::transform([](int v) { return 2 * v; });
    const Eigen::Vector2d v { 3.0, 4.0 };
    if (!doubled || *doubled != 2 * NUMERIXX_VERSION_MAJOR || v.norm() != 5.0) return 1;
    std::printf("Numerixx %s: pipes and linalg OK\n", NUMERIXX_VERSION_STRING);
    return 0;
}
