// A consumer of the scalar modules only. It must build without FXT and Eigen on the include path.
#include <numerixx/numerixx.hpp>

#include <cstdio>

#if __has_include(<fxt/monads/Expected.hpp>) || __has_include(<Eigen/Core>)
#    error "a scalar-only consumer must not see FXT or Eigen"
#endif

int main()
{
    std::printf("Numerixx %s: scalar modules OK\n", NUMERIXX_VERSION_STRING);
    return nxx::version.major == NUMERIXX_VERSION_MAJOR ? 0 : 1;
}
