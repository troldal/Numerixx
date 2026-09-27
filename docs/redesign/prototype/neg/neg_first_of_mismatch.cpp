#include "prelude.hpp"
int main() {   // search result type != root result type
    auto res = nxx::first_of(r::bisection{}.on(nxx::bracket{0.0, 2.0}), r::expand_out{}.on(nxx::bracket{0.0, 2.0}))(f);
    (void)res;
}
