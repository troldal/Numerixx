// Compile-fail (DESIGN §6.2, §9.1): the mixed literal constructor of width_tol with run-time parts. It is consteval, so
// run-time (already validated) parts are not a constant expression, and the diagnostic does not mention make(): the docs
// carry the run-time path, width_tol<T>::make(abs_tolerance, rel_tolerance) or make(a, *rel).
#include <numerixx/core.hpp>

int main(int argc, char**)
{
    const double a = static_cast<double>(argc) * 1e-10;
    const double r = static_cast<double>(argc) * 1e-8;
    const auto   A = nxx::abs_tolerance<double>::make(a);
    const auto   R = nxx::rel_tolerance<double>::make(r);
    if (A && R) {
#ifdef NUMERIXX_CF_CONTROL
        const auto w = nxx::width_tol<double>::make(*A, *R);
        (void)w;
#else
        const auto w = nxx::width_tol { *A, *R };
        (void)w;
#endif
    }
}
