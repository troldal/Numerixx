// The data of the figure in docs/sphinx/quadrature/index.md ("A weighted sum, and an error estimate"): the panels the
// default method ends with on the peak 1 / ((x - 0.3)^2 + 1e-4) over [0, 1], at the default request and at a relative
// tolerance of 1e-12, and the curve above them. The program prints both as CSV and compares them with the tables in the
// raw strings, which docs/sphinx/plots/quadrature_panels.py plots; main returns non-zero unless the two agree and what
// the page says of them holds (PHASE8_PLAN.md). The integrand and the method use only + - * /, which IEEE 754 rounds
// correctly, so the tables are the same on every target that does not contract a * b + c into a fused multiply-add
// (x86-64 without -mfma, WebAssembly).
//
// The method does not report its panels, so the program rebuilds them from the points fn receives: the first 150 are
// the ten initial panels, 15 points each, and every further 30 are the two halves of the panel it bisected, which is
// the one panel that contains all 30. Kronrod's nodes lie strictly inside a panel, so a panel is told by which points
// it contains, without the nodes' values.

#include <numerixx/quadrature.hpp>

#include "printed.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <print>
#include <span>
#include <string_view>
#include <vector>

namespace {

    namespace qd = nxx::quadrature;

    double peak(double x) { return 1.0 / ((x - 0.3) * (x - 0.3) + 1e-4); }

    struct Panel {
        double left;
        double right;
        [[nodiscard]] bool holds(std::span<double const> xs) const
        {
            return std::ranges::all_of(xs, [this](double x) { return left < x && x < right; });
        }
    };

    struct Run {
        std::vector<Panel> panels; // in order of x
        std::vector<Panel> cuts;   // the panels bisected, in order
        std::size_t evaluations;
        std::size_t bisections;
    };

    // The final panels of the default method at the request rel_tol, or nothing if the integral fails or the points do
    // not come in the blocks described above.
    std::optional<Run> run(double rel_tol)
    {
        std::vector<double> xs;
        auto const fn = [&xs](double x) {
            xs.push_back(x);
            return peak(x);
        };
        if (!qd::integrate(fn, 0.0, 1.0, {.rel_tol = rel_tol})) return std::nullopt;
        if (xs.size() < 150 || (xs.size() - 150) % 30 != 0) return std::nullopt;

        std::vector<Panel> panels;
        std::vector<Panel> cuts;
        for (std::size_t k = 0; k < 10; ++k) {
            Panel const p{.left = static_cast<double>(k) / 10.0, .right = static_cast<double>(k + 1) / 10.0};
            if (!p.holds(std::span(xs).subspan(15 * k, 15))) return std::nullopt;
            panels.push_back(p);
        }
        for (std::size_t i = 150; i < xs.size(); i += 30) {
            auto const both = std::span(xs).subspan(i, 30);
            auto const cut = std::ranges::find_if(panels, [both](Panel const &p) { return p.holds(both); });
            if (cut == panels.end()) return std::nullopt;
            double const middle = (cut->left + cut->right) / 2.0;
            Panel const low{.left = cut->left, .right = middle};
            Panel const high{.left = middle, .right = cut->right};
            // The two halves of 15 points each, in either order.
            if (!(low.holds(both.first(15)) && high.holds(both.last(15))) &&
                !(high.holds(both.first(15)) && low.holds(both.last(15))))
                return std::nullopt;
            cuts.push_back(*cut);
            *cut = high;
            panels.insert(cut, low);
        }
        return Run{.panels = panels, .cuts = cuts, .evaluations = xs.size(), .bisections = (xs.size() - 150) / 30};
    }

    double width(Panel const &p) { return p.right - p.left; }

    bool near(double x, double y) { return std::abs(x - y) <= 1e-15; } // panel ends a few halvings from k / 10

    // The narrowest panels of a run: n of them, of width narrowest, which together cover [0.275, 0.325], the 0.05
    // around the peak.
    bool narrowest_cover(Run const &r, double narrowest, long n)
    {
        auto const narrow = [narrowest](Panel const &p) { return near(width(p), narrowest); };
        auto const around = [](Panel const &p) { return p.left > 0.275 - 1e-15 && p.right < 0.325 + 1e-15; };
        return std::ranges::all_of(r.panels, [narrowest](Panel const &p) { return width(p) > narrowest - 1e-15; }) &&
               std::ranges::count_if(r.panels, narrow) == n && near(static_cast<double>(n) * narrowest, 0.05) &&
               std::ranges::all_of(r.panels, [&](Panel const &p) { return !narrow(p) || around(p); });
    }

    // [0, 0.1] and [0.5, 1] still hold their six first panels, of width 0.1.
    bool far_whole(Run const &r)
    {
        auto const far = [](Panel const &p) { return p.right <= 0.1 + 1e-15 || p.left >= 0.5 - 1e-15; };
        return std::ranges::count_if(r.panels, far) == 6 &&
               std::ranges::all_of(r.panels, [far](Panel const &p) { return !far(p) || near(width(p), 0.1); });
    }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    auto const loose = run(qd::Criteria<double>{}.rel_tol); // the default request
    auto const tight = run(1e-12);
    // Every run is a prefix of one sequence of bisections, since which panel is cut next depends on fn alone and the
    // request decides only when the run stops; a request of 0 stops at the rounding floor, the longest run of all.
    auto const longest = run(0.0);
    if (!(loose && tight && longest)) return 2;

    std::println("request,left,right");
    for (auto const &[name, r] : {std::pair{"default", &*loose}, std::pair{"1e-12", &*tight}})
        for (Panel const &p : r->panels) std::println("{},{:.6g},{:.6g}", name, p.left, p.right);
    std::println("");
    // The curve: every 0.02, and every 0.001 within 0.05 of the peak, where it changes fast.
    std::println("x,f");
    for (int k = 0; k <= 1000; ++k)
        if (k % 20 == 0 || (k >= 250 && k <= 350)) std::println("{:.4g},{:.4g}", k / 1000.0, peak(k / 1000.0));

    constexpr std::string_view panels_table = R"(
request,left,right
default,0,0.1
default,0.1,0.2
default,0.2,0.25
default,0.25,0.275
default,0.275,0.2875
default,0.2875,0.3
default,0.3,0.3125
default,0.3125,0.325
default,0.325,0.35
default,0.35,0.4
default,0.4,0.5
default,0.5,0.6
default,0.6,0.7
default,0.7,0.8
default,0.8,0.9
default,0.9,1
1e-12,0,0.1
1e-12,0.1,0.15
1e-12,0.15,0.2
1e-12,0.2,0.225
1e-12,0.225,0.25
1e-12,0.25,0.2625
1e-12,0.2625,0.275
1e-12,0.275,0.28125
1e-12,0.28125,0.2875
1e-12,0.2875,0.29375
1e-12,0.29375,0.3
1e-12,0.3,0.30625
1e-12,0.30625,0.3125
1e-12,0.3125,0.31875
1e-12,0.31875,0.325
1e-12,0.325,0.3375
1e-12,0.3375,0.35
1e-12,0.35,0.375
1e-12,0.375,0.4
1e-12,0.4,0.45
1e-12,0.45,0.5
1e-12,0.5,0.6
1e-12,0.6,0.7
1e-12,0.7,0.8
1e-12,0.8,0.9
1e-12,0.9,1

x,f
0,11.1
0.02,12.74
0.04,14.77
0.06,17.33
0.08,20.62
0.1,24.94
0.12,30.77
0.14,38.91
0.16,50.76
0.18,68.97
0.2,99.01
0.22,153.8
0.24,270.3
0.25,384.6
0.251,399.8
0.252,416
0.253,433.1
0.254,451.3
0.255,470.6
0.256,491.2
0.257,513.1
0.258,536.5
0.259,561.5
0.26,588.2
0.261,616.9
0.262,647.7
0.263,680.7
0.264,716.3
0.265,754.7
0.266,796.2
0.267,841
0.268,889.7
0.269,942.5
0.27,1000
0.271,1063
0.272,1131
0.273,1206
0.274,1289
0.275,1379
0.276,1479
0.277,1590
0.278,1712
0.279,1848
0.28,2000
0.281,2169
0.282,2358
0.283,2571
0.284,2809
0.285,3077
0.286,3378
0.287,3717
0.288,4098
0.289,4525
0.29,5000
0.291,5525
0.292,6098
0.293,6711
0.294,7353
0.295,8000
0.296,8621
0.297,9174
0.298,9615
0.299,9901
0.3,1e+04
0.301,9901
0.302,9615
0.303,9174
0.304,8621
0.305,8000
0.306,7353
0.307,6711
0.308,6098
0.309,5525
0.31,5000
0.311,4525
0.312,4098
0.313,3717
0.314,3378
0.315,3077
0.316,2809
0.317,2571
0.318,2358
0.319,2169
0.32,2000
0.321,1848
0.322,1712
0.323,1590
0.324,1479
0.325,1379
0.326,1289
0.327,1206
0.328,1131
0.329,1063
0.33,1000
0.331,942.5
0.332,889.7
0.333,841
0.334,796.2
0.335,754.7
0.336,716.3
0.337,680.7
0.338,647.7
0.339,616.9
0.34,588.2
0.341,561.5
0.342,536.5
0.343,513.1
0.344,491.2
0.345,470.6
0.346,451.3
0.347,433.1
0.348,416
0.349,399.8
0.35,384.6
0.36,270.3
0.38,153.8
0.4,99.01
0.42,68.97
0.44,50.76
0.46,38.91
0.48,30.77
0.5,24.94
0.52,20.62
0.54,17.33
0.56,14.77
0.58,12.74
0.6,11.1
0.62,9.756
0.64,8.643
0.66,7.71
0.68,6.92
0.7,6.246
0.72,5.666
0.74,5.163
0.76,4.724
0.78,4.338
0.8,3.998
0.82,3.697
0.84,3.428
0.86,3.188
0.88,2.972
0.9,2.777
0.92,2.601
0.94,2.441
0.96,2.295
0.98,2.162
1,2.04
)";
    if (!(docs::printed() == docs::quoted(panels_table))) return 1;

    // What the page says. At the default request: 6 cuts, 16 panels (330 evaluations), the narrowest four of 1/80
    // covering [0.275, 0.325], each cut a panel that ends at 0.3, the peak.
    if (!(loose->evaluations == 330 && loose->bisections == 6 && loose->panels.size() == 16)) return 3;
    if (!narrowest_cover(*loose, 1.0 / 80.0, 4)) return 4;
    if (!std::ranges::all_of(loose->cuts, [](Panel const &p) { return near(p.left, 0.3) || near(p.right, 0.3); }))
        return 5;
    // Asked for 1e-12: 16 cuts, 26 panels (630 evaluations), the narrowest eight of 1/160 over the same stretch.
    if (!(tight->evaluations == 630 && tight->bisections == 16 && tight->panels.size() == 26)) return 6;
    if (!narrowest_cover(*tight, 1.0 / 160.0, 8)) return 7;
    // [0, 0.1] and [0.5, 1] keep their first panels whatever the request: in both runs, and in the longest (690
    // evaluations, 18 cuts), of which every other run is a prefix.
    if (!(longest->evaluations == 690 && longest->bisections == 18)) return 8;
    if (!(far_whole(*loose) && far_whole(*tight) && far_whole(*longest))) return 9;
    // The peak reaches 10,000 at 0.3.
    if (!(std::abs(peak(0.3) - 1e4) < 1e-8)) return 10;
    return 0;
}
