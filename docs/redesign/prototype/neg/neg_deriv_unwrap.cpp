// Plan 6.12: derivative_of returns x -> expected<T, failure<T, UE>>; 6.4 only unwraps expected<T, fault<UE>>.
// Without the extra unwrap rule, Newton's UE becomes failure<double, none> and the chain no longer type-checks.
#include "prelude.hpp"
int main() {
    auto res = nxx::first_of(r::newton{}.with_derivative(d::derivative_of(f)).on(1.0),
                             r::brent{}.on(nxx::bracket{0.0, 2.0}))(f);
    (void)res;
}
