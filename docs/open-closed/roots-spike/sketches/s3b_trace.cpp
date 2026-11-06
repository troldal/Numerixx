// Sketch 3b: why the naive Pegasus takes 33 calls on the cubic; the naive sign test; overflow in the user's arithmetic.
#define main s3_main
#include "s3_bracketing.cpp"
#undef main
int main()
{
    auto const cubic = [](double x) { return x * x * x - 2.0 * x - 5.0; };
    auto const print = [](rt::Step<double> const &s) {
        std::println("  [{:.17g}, {:.17g}] width {:.3g} f=({:.3g}, {:.3g}) point {:.17g}", s.lower, s.upper,
                     s.upper - s.lower, s.f_lower, s.f_upper, s.point);
    };
    std::println("Pegasus on the cubic:");
    (void)rt::solve<user::Pegasus<>>(rt::observer(print), cubic, rt::Bracket{ 2.0, 3.0 });

    // fa * fb underflows: values about 1e-200.
    auto const tiny = [](double x) { return 1e-200 * (x * x * x - 2.0 * x - 5.0); };
    std::println("-- tiny scale, nonlinear");
    run<user::Pegasus<>>("Pegasus", tiny, { 2.0, 3.0 }, 2.0945514815423265);
    run<user::Pegasus<true>>("Pegasus (naive sign)", tiny, { 2.0, 3.0 }, 2.0945514815423265);
    run<rt::AndersonBjorck>("AndersonBjorck (lib)", tiny, { 2.0, 3.0 }, 2.0945514815423265);

    // fn finite at the ends of a bracket whose width overflows.
    auto const lin = [](double x) { return x - 1.5; };
    auto const atan_ = [](double x) { return std::atan(x) - 0.5; }; // root tan(0.5)
    std::println("-- huge bracket [-1e308, 1e308]");
    run<user::Pegasus<>>("Pegasus", lin, { -1e308, 1e308 }, 1.5);
    run<user::Trisection>("Trisection", lin, { -1e308, 1e308 }, 1.5);
    run<user::Pegasus<>>("Pegasus atan", atan_, { -1e308, 1e308 }, std::tan(0.5));
    run<user::Trisection>("Trisection atan", atan_, { -1e308, 1e308 }, std::tan(0.5));
    run<rt::AndersonBjorck>("AndersonBjorck atan", atan_, { -1e308, 1e308 }, std::tan(0.5));
}
