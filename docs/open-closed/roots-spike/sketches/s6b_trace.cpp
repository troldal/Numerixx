// Sketch 6b: trace of NewtonBisect on the cubic.
#define main s6_main
#include "s6_values2.cpp"
#undef main
int main()
{
    auto const fdf = [](double x) { return std::pair{ x * x * x - 2.0 * x - 5.0, 3.0 * x * x - 2.0 }; };
    (void)rt::solve<user::NewtonBisect>(rt::observer([](rt::Step<double> const &s) {
        std::println("  [{:.17g}, {:.17g}] w {:.3g} f=({:.3g}, {:.3g})", s.lower, s.upper, s.upper - s.lower, s.f_lower, s.f_upper);
    }), fdf, rt::Bracket{ 2.0, 3.0 });
    return 0;
}
