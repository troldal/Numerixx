// Sketch 2: each step printed, on a bracket and from a guess; a recorder kept by reference.
#include <numerixx/roots.hpp>
#include <print>
#include <vector>

namespace rt = nxx::roots;

int main()
{
    auto const f = [](double x) { return x * x * x - 2.0 * x - 5.0; }; // root 2.0945514815423265
    auto const fdf = [](double x) { return std::pair{ x * x * x - 2.0 * x - 5.0, 3.0 * x * x - 2.0 }; };

    auto const print = [](rt::Step<double> const &s) {
        std::println("  [{:.17g}, {:.17g}] f=({:.3g}, {:.3g}) point {:.17g}{}", s.lower, s.upper, s.f_lower, s.f_upper,
                     s.point, s.exhausted ? " exhausted" : "");
    };

    std::println("Brent on [2, 3]:");
    auto const r1 = rt::solve<rt::Brent>(rt::observer(print), f, rt::Bracket{ 2.0, 3.0 });
    std::println(" -> {:.17g} +- {:.3g}", r1->value, r1->abs_error);

    std::println("Newton from 2:");
    auto const r2 = rt::solve<rt::Newton>(rt::observer(print), fdf, 2.0);
    std::println(" -> {:.17g} +- {:.3g}", r2->value, r2->abs_error);

    std::println("Secant from 2 (method as a value):");
    auto const r3 = rt::solve(rt::Secant{}, rt::observer(print), f, 2.0);
    std::println(" -> {:.17g} +- {:.3g}", r3->value, r3->abs_error);

    // A recorder the caller owns, kept by reference.
    std::vector<rt::Step<double>> steps;
    auto record = [&steps](rt::Step<double> const &s) { steps.push_back(s); };
    auto const r4 = rt::solve<rt::Bisection>(rt::observer(record), f, rt::Bracket{ 2.0, 3.0 });
    std::println("Bisection: {} steps recorded, last point {:.17g}, result {:.17g}", steps.size(), steps.back().point,
                 r4->value);

    // A stateful functor by lvalue: is it kept by reference? (count survives the call)
    struct Counter
    {
        int n = 0;
        void operator()(rt::Step<double> const &) { ++n; }
    } counter;
    (void)rt::solve<rt::Toms748>(rt::observer(counter), f, rt::Bracket{ 2.0, 3.0 });
    std::println("Toms748 counter (lvalue functor): {}", counter.n);

    // The same functor as an rvalue: a copy, so the caller's count does not change.
    Counter c2;
    (void)rt::solve<rt::Toms748>(rt::observer(Counter{}), f, rt::Bracket{ 2.0, 3.0 });
    std::println("rvalue functor: caller's copy stays {}", c2.n);

    // Observed == unobserved
    std::println("observed == unobserved: {} {} {}", r1 == rt::solve<rt::Brent>(f, rt::Bracket{ 2.0, 3.0 }),
                 r2 == rt::solve<rt::Newton>(fdf, 2.0), r3 == rt::solve<rt::Secant>(f, 2.0));

    // A generic observer: auto parameter.
    int generic = 0;
    (void)rt::solve<rt::Ridders>(rt::observer([&](auto const &) { ++generic; }), f, rt::Bracket{ 2.0, 3.0 });
    std::println("Ridders, generic observer: {} reports", generic);

    // Newton on a bracket, observed, through the pair-returning fn.
    int nb = 0;
    (void)rt::solve<rt::Newton>(rt::observer([&](rt::Step<double> const &) { ++nb; }), fdf, rt::Bracket{ 2.0, 3.0 });
    std::println("Newton on a bracket: {} reports", nb);
}
