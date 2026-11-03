// The program of docs/sphinx/quadrature/costs.md: the evaluations, error, estimate and status per method in the page's
// tables, at the default request, in double. It asserts each entry as the page types it. The Python prototypes in
// docs/quadrature/costs-prototype/ (python -I costs.py prints them, run 2026-11-03) reproduce: Gauss-Kronrod's counts,
// from a K15 rule of their own (values within two ulps, estimates to the digits printed); TanhSinh's counts, values and
// estimates bit for bit, through docs/quadrature/tanh-sinh-prototype/ts.py and es.py, which read the C++ node tables
// (so they check the loop and the estimate, not the nodes); Romberg's counts, values and estimates bit for bit on the
// first two integrands. Romberg on 1 / sqrt x: the prototype raises at fn(0), while the C++ evaluates fn(a) then fn(b)
// before checking, 2 evaluations by its definition; GaussLegendre<N> evaluates fn N times by definition, and a
// method without infinite ends refuses one before calling fn. The errors are against the closed forms pi, 2/3 and 2,
// formed in double (pi to half an ulp, 2/3 to half an ulp of 0.67: both far below the two digits printed). The
// integrands use only + - * / and
// std::sqrt, so the entries hold on every target that does not contract a * b + c into a fused multiply-add. The page
// shows no code from this program.

#include <numerixx/quadrature.hpp>

#include <cmath>
#include <format>
#include <limits>
#include <numbers>
#include <print>
#include <string>
#include <string_view>

namespace {

    namespace qd = nxx::quadrature;

    int failures = 0;

    // One row of a table: the method's evaluations, then the actual error, the estimate and the status as the page
    // writes them ("-" for no estimate; a status of "number" for a fixed rule's plain result).
    struct Row
    {
        std::string_view method;
        int evaluations;
        std::string_view error;
        std::string_view estimate;
        std::string_view status;
    };

    std::string two_digits(double x) { return x == 0.0 ? std::string("0") : std::format("{:.1e}", x); }

    void check(Row const &expected, int evaluations, std::string const &error, std::string const &estimate,
               std::string const &status)
    {
        if (evaluations == expected.evaluations && error == expected.error && estimate == expected.estimate
            && status == expected.status)
            return;
        ++failures;
        std::println(stderr, "{}: got {} evaluations, error {}, estimate {}, {}; the page says {}, {}, {}, {}",
                     expected.method, evaluations, error, estimate, status, expected.evaluations, expected.error,
                     expected.estimate, expected.status);
    }

    // A method with an estimate, on fn from a to b, against the exact integral.
    template<typename M, typename F>
    void with_estimate(Row const &expected, F const &f, double a, double b, double exact)
    {
        int evaluations = 0;
        auto const r = qd::integrate<M>([&](double x) { ++evaluations; return f(x); }, a, b);
        if (r) {
            check(expected, evaluations, two_digits(std::abs(r->value - exact)), two_digits(r->abs_error), "ok");
            // The page: "every estimate covers its error".
            if (!(std::abs(r->value - exact) <= r->abs_error)) ++failures;
        }
        else {
            check(expected, evaluations, "-", "-", std::format("{}", r.error()));
        }
    }

    // A fixed rule.
    template<int N, typename F>
    void fixed(Row const &expected, F const &f, double a, double b, double exact)
    {
        int evaluations = 0;
        auto const r = qd::integrate<qd::GaussLegendre<N>>([&](double x) { ++evaluations; return f(x); }, a, b);
        if (r) check(expected, evaluations, two_digits(std::abs(*r - exact)), "-", "number");
        else check(expected, evaluations, "-", "-", std::format("{}", r.error()));
    }

} // namespace

int main()
{
    double const pi = std::numbers::pi;
    double const inf = std::numeric_limits<double>::infinity();
    auto const arctan = [](double x) { return 4.0 / (1.0 + x * x); };
    auto const root = [](double x) { return std::sqrt(x); };
    auto const inverse_root = [](double x) { return 1.0 / std::sqrt(x); };
    auto const lorentz = [](double x) { return 1.0 / (1.0 + x * x); };

    // 4 / (1 + x^2) on [0, 1] = pi
    with_estimate<qd::GaussKronrod<15>>({"GaussKronrod<15>", 150, "0", "3.5e-14", "ok"}, arctan, 0.0, 1.0, pi);
    with_estimate<qd::TanhSinh>({"TanhSinh", 121, "4.4e-16", "3.5e-14", "ok"}, arctan, 0.0, 1.0, pi);
    with_estimate<qd::Romberg>({"Romberg", 33, "4.8e-11", "1.2e-08", "ok"}, arctan, 0.0, 1.0, pi);
    fixed<10>({"GaussLegendre<10>", 10, "2.5e-13", "-", "number"}, arctan, 0.0, 1.0, pi);
    fixed<20>({"GaussLegendre<20>", 20, "8.9e-16", "-", "number"}, arctan, 0.0, 1.0, pi);

    // sqrt x on [0, 1] = 2/3: a singular derivative at 0
    with_estimate<qd::GaussKronrod<15>>({"GaussKronrod<15>", 360, "2.9e-10", "5.1e-09", "ok"}, root, 0.0, 1.0,
                                        2.0 / 3.0);
    with_estimate<qd::TanhSinh>({"TanhSinh", 116, "0", "7.5e-15", "ok"}, root, 0.0, 1.0, 2.0 / 3.0);
    with_estimate<qd::Romberg>({"Romberg", 65537, "4.1e-09", "7.5e-09", "ok"}, root, 0.0, 1.0, 2.0 / 3.0);
    fixed<10>({"GaussLegendre<10>", 10, "8.9e-05", "-", "number"}, root, 0.0, 1.0, 2.0 / 3.0);
    fixed<20>({"GaussLegendre<20>", 20, "1.2e-05", "-", "number"}, root, 0.0, 1.0, 2.0 / 3.0);

    // 1 / sqrt x on [0, 1] = 2: a singularity at 0
    with_estimate<qd::GaussKronrod<15>>({"GaussKronrod<15>", 1350, "1.4e-08", "2.1e-08", "ok"}, inverse_root, 0.0,
                                        1.0, 2.0);
    with_estimate<qd::TanhSinh>({"TanhSinh", 131, "0", "2.2e-14", "ok"}, inverse_root, 0.0, 1.0, 2.0);
    with_estimate<qd::Romberg>({"Romberg", 2, "-", "-", "non_finite"}, inverse_root, 0.0, 1.0, 2.0);
    fixed<10>({"GaussLegendre<10>", 10, "8.3e-02", "-", "number"}, inverse_root, 0.0, 1.0, 2.0);
    fixed<20>({"GaussLegendre<20>", 20, "4.2e-02", "-", "number"}, inverse_root, 0.0, 1.0, 2.0);

    // 1 / (1 + x^2) on the line = pi: an infinite interval, TanhSinh's alone
    with_estimate<qd::GaussKronrod<15>>({"GaussKronrod<15>", 0, "-", "-", "non_finite"}, lorentz, -inf, inf, pi);
    with_estimate<qd::TanhSinh>({"TanhSinh", 133, "0", "3.5e-14", "ok"}, lorentz, -inf, inf, pi);
    with_estimate<qd::Romberg>({"Romberg", 0, "-", "-", "non_finite"}, lorentz, -inf, inf, pi);
    fixed<10>({"GaussLegendre<10>", 0, "-", "-", "non_finite"}, lorentz, -inf, inf, pi);

    return failures == 0 ? 0 : 1;
}
