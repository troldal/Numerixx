// A consumer translation unit compiled with strict warnings as errors (DESIGN §9.1, §5.3). The global `f` makes
// MSVC's C4459 ("declaration hides global declaration") fire in any Numerixx header that names a parameter `f`,
// so the library names callback parameters `fn` or `func`.
#include <numerixx/numerixx.hpp>

double f(double x);
double f(double x) { return x * x - 2.0; }

namespace consumer
{
    int uses_numerixx() { return nxx::version.major + static_cast<int>(f(2.0)); }
}    // namespace consumer
