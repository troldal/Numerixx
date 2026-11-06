// Sketch 1: today's calls, unchanged, and the method given as a value.
#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>
#include <print>
#include <utility>

namespace rt = nxx::roots;

int main()
{
    auto const f = [](double x) { return x * x - 2.0; };
    auto const fdf = [](double x) { return std::pair{ x * x - 2.0, 2.0 * x }; };
    auto const df = [](double x) { return 2.0 * x; };

    auto const a = rt::solve<rt::Brent>(f, rt::Bracket{ 1.0, 2.0 });
    auto const b = rt::solve<rt::Secant>(f, 1.0);
    auto const c = rt::solve<rt::Newton>(fdf, 1.0);
    auto const d = rt::solve<rt::Newton>(f, df, rt::Bracket{ 1.0, 2.0 });
    auto const e = rt::solve<rt::Newton>(f, nxx::deriv::derivative_of(f), 1.0);
    auto const g = rt::solve<rt::Bisection>(f, rt::Bracket{ 2.0, 1.0 }, { .abs_tol = 1e-12 });
    auto const h = rt::solve<rt::Toms748>([](auto x) { return x * x - 2.0; }, rt::Bracket{ 1.0, 2.0 }); // generic lambda

    // New: the method as a value.
    auto const v1 = rt::solve(rt::Toms748{}, f, rt::Bracket{ 1.0, 2.0 }, { .rel_tol = 1e-12 });
    auto const v2 = rt::solve(rt::Secant{}, f, 1.0);
    auto const v3 = rt::solve(rt::Newton{}, fdf, rt::Bracket{ 1.0, 2.0 });
    rt::Brent const brent;
    auto const v4 = rt::solve(brent, f, rt::Bracket{ 1.0, 2.0 });

    for (auto const &[name, r] : { std::pair{ "brent", a }, { "secant", b }, { "newton", c }, { "newton sep", d },
                                   { "newton num", e }, { "bisection", g }, { "toms748 generic", h },
                                   { "value toms748", v1 }, { "value secant", v2 }, { "value newton", v3 },
                                   { "value brent", v4 } }) {
        if (r) std::println("{:16} {:.17g} +- {:.3g}", name, r->value, r->abs_error);
        else std::println("{:16} error {}", name, r.error());
    }
    // Same result named and as a value:
    std::println("named == value (Toms748): {}", rt::solve<rt::Toms748>(f, rt::Bracket{ 1.0, 2.0 }, { .rel_tol = 1e-12 }) == v1);
    std::println("named == value (Brent): {}", a == v4);
}
