// Plan 7.1 signature verbatim: the defaulted stencil parameter cannot deduce O, A, N.
#include "prelude.hpp"
template<class F, nxx::real T, int O, int A, std::size_t N>
constexpr auto diff_plan(const F& fn, T x, const d::stencil<O, A, N>& s = d::central_1_2) { return d::diff(fn, x, s); }
int main() { auto v = diff_plan(f, 1.0); (void)v; }
