#include "prelude.hpp"
int main() { auto res = r::newton{nxx::width_tol{1e-6}}.with_derivative(df)(f, 1.0); (void)res; }
