// The double tables of quadrature.hpp for the prototypes, t = j / 512, each value to 17 digits:
// - the tanh-sinh table, as lines "level j comp w", for ts.py:
//   g++ -std=c++23 -I include docs/quadrature/tanh-sinh-prototype/dump.cpp -o dump && ./dump > tab_double.txt
// - with the argument es, the exp-sinh table (c = 3/2, u = c sinh t), as lines "level j q f" for a node toward the
//   finite end only and "level j q f E x" for one toward infinity too (q = exp(-u), f = c cosh t, E = 1/q, x the
//   line's node sinh u; the weights are formed at use, f q, f E and (f E + f q) / 2):
//   ./dump es > es_double.txt
#include <numerixx/quadrature.hpp>

#include <cstdio>
#include <string_view>

namespace detail = nxx::quadrature::detail;

int main(int argc, char **argv)
{
    bool const exp_sinh = argc > 1 && std::string_view(argv[1]) == "es";
    for (int level = 0; level < detail::tanh_sinh_levels; ++level) {
        if (exp_sinh) {
            detail::with_exp_sinh_level<double>(level, [level](detail::ExpSinhLevel<double> const &nodes) {
                for (int i = 0; i < nodes.inward; ++i) {
                    std::printf("%d %d %.17g %.17g", level, detail::position(level, i), nodes.q(i), nodes.f(i));
                    if (i < nodes.outward) std::printf(" %.17g %.17g", nodes.e(i), nodes.x(i));
                    std::printf("\n");
                }
            });
            continue;
        }
        detail::with_tanh_sinh_level<double>(level, [level](detail::TanhSinhLevel<double> const &nodes) {
            for (int i = 0; i < nodes.size; ++i) {
                int const j = level == 0 ? i << 9 : (2 * i + 1) << (9 - level);
                std::printf("%d %d %.17g %.17g\n", level, j, nodes.complement(i), nodes.weight(i));
            }
        });
    }
}
