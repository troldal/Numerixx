// Sketch 9: user methods in constant evaluation, observed too; and in cpp_bin_float_50 at run time.
#define main s3c_main
#include "s3c_pegasus_fixed.cpp"
#undef main
#define main s4_main
#define run s4_run
#include "s4_open.cpp"
#undef run
#undef main
#define main s7_main
#include "s7_parameter.cpp"
#undef main
#include <boost/multiprecision/cpp_bin_float.hpp>

constexpr double cubic(double x) { return x * x * x - 2.0 * x - 5.0; }
constexpr double exact = 2.0945514815423265;

constexpr auto pegasus = rt::solve<user::Pegasus>(cubic, rt::Bracket{ 2.0, 3.0 });
static_assert(pegasus && pegasus->value == exact && pegasus->abs_error <= 1e-15);

constexpr auto section = rt::solve(user::Section<>{ .fraction = 0.25 }, cubic, rt::Bracket{ 2.0, 3.0 });
static_assert(section && section->value - exact <= section->abs_error && exact - section->value <= section->abs_error);

constexpr auto chebyshev = rt::solve<user::Chebyshev<>>(
    [](double x) { return std::array{ x * x * x - 2.0 * x - 5.0, 3.0 * x * x - 2.0, 6.0 * x }; }, 2.0);
static_assert(chebyshev && chebyshev->value == exact);

// Observed in constant evaluation: count the steps.
constexpr int pegasus_steps = [] {
    int n = 0;
    (void)rt::solve<user::Pegasus>(rt::observer([&n](rt::Step<double> const &) { ++n; }), cubic, rt::Bracket{ 2.0, 3.0 });
    return n;
}();
static_assert(pegasus_steps == 9); // the start and 8 advances: s3c counted 10 calls of fn, 2 of them at the ends

int main()
{
    std::println("constexpr: pegasus {:.17g}, section {:.17g}, chebyshev {:.17g}, pegasus steps {}", pegasus->value,
                 section->value, chebyshev->value, pegasus_steps);
    using big = boost::multiprecision::cpp_bin_float_50;
    big const exact50("2.0945514815423265914823865405793029638573061056282391803041285");
    auto const f50 = [](big const &x) -> big { return x * x * x - 2 * x - 5; };
    for (auto const &[name, r] : { std::pair{ "Pegasus 50", rt::solve<user::Pegasus>(f50, rt::Bracket<big>{ 2, 3 }) },
                                   { "Section 50", rt::solve(user::Section<>{}, f50, rt::Bracket<big>{ 2, 3 }) },
                                   { "AndersonBjorck 50", rt::solve<rt::AndersonBjorck>(f50, rt::Bracket<big>{ 2, 3 }) } })
        std::println("{:18} |x - exact| = {}  abs_error = {}", name, static_cast<double>(abs(r->value - exact50)),
                     static_cast<double>(r->abs_error));
}
