// Spike (draft 4): print the Step of every iteration, the example in PHASE3A_API.md, "Example: printing each step".
#include <numerixx/quadrature.hpp>

#include <cmath>
#include <cstdio>

namespace qd = nxx::quadrature;

// Prints one Step: the iteration, the evaluations of fn so far, and every field of the step.
struct Printer
{
    std::size_t const *evaluations; // counted by fn
    int iteration = 0;

    void operator()(qd::Step<double> const &s)
    {
        std::printf("%4d %6zu  ", iteration++, *evaluations);
        if (s.current)
            std::printf("%22.17f %10.3g", s.current->value, s.current->abs_error);
        else
            std::printf("%22s %10s", "-", "-");
        std::printf(" %10.3g  %-5s %-9s %s\n", s.attainable, s.ready ? "yes" : "no", s.exhausted ? "yes" : "no",
            s.trusted ? "yes" : "no");
    }
};

void header(char const *title)
{
    std::printf("%s\n%4s %6s  %22s %10s %10s  %-5s %-9s %s\n", title, "step", "evals", "value", "abs_error",
        "attainable", "ready", "exhausted", "trusted");
}

template<typename R>
void footer(R const &r)
{
    if (r) std::printf("result: %.17f +- %.3g\n\n", r->value, r->abs_error);
    else std::printf("no result: reason %d\n\n", static_cast<int>(r.error().reason));
}

int main()
{
    std::size_t evaluations = 0;

    // Romberg on e^x over [0, 1] at rel_tol 1e-12: one step per level.
    auto const e = [&evaluations](double x) {
        ++evaluations;
        return std::exp(x);
    };
    header("Romberg, e^x on [0, 1], rel_tol 1e-12");
    Printer romberg{ &evaluations };
    footer(qd::integrate(qd::observed(qd::Romberg{}, romberg), e, 0.0, 1.0, { .rel_tol = 1e-12 }));

    // Gauss-Kronrod on a peak of width 1e-3 at 0.3 over [0, 1] at rel_tol 1e-10: the mesh, then one step per bisection.
    evaluations = 0;
    auto const peak = [&evaluations](double x) {
        ++evaluations;
        double const d = x - 0.3;
        return 1.0 / (d * d + 1e-6);
    };
    header("GaussKronrod<15>, 1 / ((x - 0.3)^2 + 1e-6) on [0, 1], rel_tol 1e-10");
    Printer kronrod{ &evaluations };
    footer(qd::integrate(qd::observed(qd::GaussKronrod<15>{}, kronrod), peak, 0.0, 1.0, { .rel_tol = 1e-10 }));
}
