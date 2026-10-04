// Compile-fail (DESIGN §3.3, §6.8): a bare number as Brent's tolerance. brent's width_tolerance_v asked
// W::applies_to of every W (an && in a variable template's initializer does not stop the instantiation of its later
// operands), so brent{1e-10} was a hard error inside brent.hpp ("'applies_to' is not a member of double"), and
// std::is_constructible_v and CTAD probes were hard errors too. A deleted constructor now says what to write.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    auto s = r::brent { nxx::width_tol { 1e-10 } };
#else
    auto s = r::brent { 1e-10 };
#endif
    (void)s;
}
