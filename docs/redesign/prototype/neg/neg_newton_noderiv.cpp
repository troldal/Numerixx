#include "prelude.hpp"
int main() { auto res = r::newton{}(f, 1.0); (void)res; }        // D13: deleted with reason
