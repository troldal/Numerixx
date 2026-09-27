// Are solver values and combinator results regular values (copy AND assign)?
#include "prelude.hpp"
#include <cstdio>
#include <optional>
int main() {
    auto chain = nxx::first_of(r::newton{}.with_derivative(df).on(0.0), r::bisection{}.on(nxx::bracket{0.0, 2.0}));
    auto piped = nxx::then(r::bisection{nxx::x_tol{1e-4}}.on(nxx::bracket{0.0, 2.0}), r::newton{}.with_derivative(df));
    auto solver = r::newton{}.with_derivative(df);
    auto curried = r::bisection{}.on(nxx::bracket{0.0, 2.0});
    auto dfn = d::derivative_of([k = 2.0](double x) { return k * x * x; });   // capturing lambda inside
    std::printf("copy-constructible: chain=%d then=%d solver=%d curried=%d derivative_of=%d\n",
        std::is_copy_constructible_v<decltype(chain)>, std::is_copy_constructible_v<decltype(piped)>,
        std::is_copy_constructible_v<decltype(solver)>, std::is_copy_constructible_v<decltype(curried)>,
        std::is_copy_constructible_v<decltype(dfn)>);
    std::printf("copy-assignable:    chain=%d then=%d solver=%d curried=%d derivative_of=%d\n",
        std::is_copy_assignable_v<decltype(chain)>, std::is_copy_assignable_v<decltype(piped)>,
        std::is_copy_assignable_v<decltype(solver)>, std::is_copy_assignable_v<decltype(curried)>,
        std::is_copy_assignable_v<decltype(dfn)>);
    std::printf("std::regular / semiregular: solver=%d/%d chain=%d/%d\n", std::regular<decltype(solver)>, std::semiregular<decltype(solver)>,
        std::regular<decltype(chain)>, std::semiregular<decltype(chain)>);
    std::optional<decltype(chain)> slot;
    slot.emplace(chain);            // fine
    // slot = chain;                // would not compile: closure copy-assignment is deleted
    std::printf("optional<chain>.emplace works: %d\n", slot.has_value());
}
