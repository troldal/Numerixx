#include "nxx/roots.hpp"
#include "nxx/deriv.hpp"
namespace r = nxx::roots;
namespace d = nxx::deriv;
inline constexpr auto f  = [](double x) { return x * x - 2.0; };
inline constexpr auto df = [](double x) { return 2.0 * x; };
