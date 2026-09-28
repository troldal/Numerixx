// f and f' with different user error types: the plan does not say what Newton's UE is.
#include "prelude.hpp"
enum class e1 { a }; enum class e2 { b };
int main() {
    auto g  = [](double x) -> std::expected<double, e1> { return x * x - 2.0; };
    auto dg = [](double x) -> std::expected<double, e2> { return 2.0 * x; };
    auto res = r::newton{}.with_derivative(dg)(g, 1.0);
    (void)res;
}
