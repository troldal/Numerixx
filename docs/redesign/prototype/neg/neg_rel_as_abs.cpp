#include "prelude.hpp"
int main() { nxx::x_tol<double> c{nxx::rel_tolerance<double>{0.1}}; (void)c; }   // roles are not interchangeable
