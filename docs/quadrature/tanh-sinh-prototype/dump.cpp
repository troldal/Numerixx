// The double tanh-sinh table of quadrature.hpp as lines "level j comp w", t = j / 512, for the prototype (ts.py):
// g++ -std=c++23 -I include docs/quadrature/tanh-sinh-prototype/dump.cpp -o dump && ./dump > tab_double.txt
#include <numerixx/quadrature.hpp>

#include <cstdio>

namespace detail = nxx::quadrature::detail;

int main()
{
    for (int level = 0; level < detail::tanh_sinh_levels; ++level)
        detail::with_tanh_sinh_level<double>(level, [level](detail::TanhSinhLevel<double> const &nodes) {
            for (int i = 0; i < nodes.size; ++i) {
                int const j = level == 0 ? i << 9 : (2 * i + 1) << (9 - level);
                std::printf("%d %d %.17g %.17g\n", level, j, nodes.complement(i), nodes.weight(i));
            }
        });
}
