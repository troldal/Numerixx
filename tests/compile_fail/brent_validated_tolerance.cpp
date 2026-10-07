// Compile-fail (DESIGN §6.6, §7.2): a validated tolerance is not a criterion. brent{*tol} was a long CTAD error with no
// reason; brent's bare-number deletion now takes validated tolerances and their parts too, and its explicit deduction
// guide sends them to brent<>, so the deletion reports, not CTAD.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
    if (auto tol = nxx::tolerance<double>::make(1e-10)) {
#ifdef NUMERIXX_CF_CONTROL
        auto s = r::brent { nxx::width_tol { *tol } };
#else
        auto s = r::brent { *tol };
#endif
        (void)s;
    }
}
