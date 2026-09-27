#include "prelude.hpp"
int main() { auto res = nxx::then(r::newton{}.with_derivative(df).on(1.0), r::bisection{})(f); (void)res; }
