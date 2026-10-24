// Compiled by the warning tests in CMakeLists.txt, never linked. A derivative above the first in float selects a
// deprecated overload of derivative(): with warnings as errors the "warns_step" and "warns_default" cases (one per
// overload, with a step and with the default step) must fail to build, while the "silenced" case (a pragma at the
// calls) and the "macro" case (NUMERIXX_NO_PRECISION_WARNING) must build both. A first derivative in float and every
// order in double never warn.

#include <numerixx/deriv.hpp>

namespace {
    constexpr auto cube = [](auto x) { return x * x * x; };
} // namespace

#if defined(NUMERIXX_WARNING_SILENCED)
#if defined(_MSC_VER) && !defined(__clang__)
#pragma warning(push)
#pragma warning(disable : 4996)
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#endif
#if !defined(NUMERIXX_WARNING_DEFAULT_ONLY)
bool second_float_step() { return nxx::deriv::derivative<nxx::deriv::Central<2>, 2>(cube, 1.0f, 0.25f).has_value(); }
#endif
#if !defined(NUMERIXX_WARNING_STEP_ONLY)
bool second_float_default() { return nxx::deriv::derivative<nxx::deriv::ForwardRidders, 2>(cube, 1.0f).has_value(); }
#endif
#if defined(NUMERIXX_WARNING_SILENCED)
#if defined(_MSC_VER) && !defined(__clang__)
#pragma warning(pop)
#else
#pragma GCC diagnostic pop
#endif
#endif

bool first_float()
{
    using namespace nxx::deriv;
    return derivative<Central<2>>(cube, 1.0f).has_value() && derivative<Ridders>(cube, 1.0f).has_value()
           && derivative<BackwardRidders, order::first>(cube, 1.0f, 0.25f).has_value();
}
bool fourth_double() { return nxx::deriv::derivative<nxx::deriv::Ridders, 4>(cube, 1.0).has_value(); }
