// Plan 6.9 engine-held loop, verbatim shape: `s = stepper.step(*problem, *s)` where s came from init().
// init() returns expected<S, failure<Est,UE>>, step() returns expected<S, fault<UE>>: not assignable.
#include "prelude.hpp"
int main() {
    const r::brent<> stepper{};
    auto problem = stepper.prepare(std::cref(f), nxx::bracket{1.0, 2.0});
    auto s = stepper.init(*problem);
    int budget = 10;
    while (s && !stepper.intrinsic(*s) && budget-- > 0) s = stepper.step(*problem, *s);
    return s ? 0 : 1;
}
