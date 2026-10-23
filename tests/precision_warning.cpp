// Compiled by the warning tests in CMakeLists.txt, never linked. A 4th derivative in float selects the deprecated
// overload of derivative(): with warnings as errors the "warns" case must fail to build, while the "silenced" case (a
// pragma at the call) and the "macro" case (NUMERIXX_NO_PRECISION_WARNING) must build. double never warns.

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
bool fourth_float() { return nxx::deriv::derivative<nxx::deriv::Central<2>, 4>(cube, 1.0f).has_value(); }
#if defined(NUMERIXX_WARNING_SILENCED)
#if defined(_MSC_VER) && !defined(__clang__)
#pragma warning(pop)
#else
#pragma GCC diagnostic pop
#endif
#endif

bool fourth_double() { return nxx::deriv::derivative<nxx::deriv::Ridders, 4>(cube, 1.0).has_value(); }
bool third_float() { return nxx::deriv::derivative<nxx::deriv::Ridders, 3>(cube, 1.0f).has_value(); }
