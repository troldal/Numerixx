// Plan 6.2 declares `using max_evaluations = detail::refined<tag::max_evaluations, std::uint32_t>` while 6.8
// uses `max_evaluations{n}` as a stop criterion. Both cannot be the same name.
#include "prelude.hpp"
namespace plan { struct max_evals_tag { static constexpr const char* message = "x"; static constexpr bool check(std::uint32_t v) { return v > 0; } };
                 using max_evaluations = nxx::detail::refined<max_evals_tag, std::uint32_t>; }
int main() { auto c = nxx::never{} || plan::max_evaluations{20}; (void)c; }
