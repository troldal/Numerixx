// Battery for docs/roots/PHASE7_API.md: every method on many functions, starts and tolerances, in five types. Prints
// every result with its count of evaluations, so that the output of today's roots.hpp and of the spike's can be
// compared line by line. Uses only the forms of solve() that exist today.
#include <numerixx/roots.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <array>
#include <cmath>
#include <format>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace rt = nxx::roots;
using bin50 = boost::multiprecision::cpp_bin_float_50;
using dec50 = boost::multiprecision::cpp_dec_float_50;

template<typename T>
std::string show(T const &x)
{
    std::ostringstream o;
    o << std::setprecision(std::numeric_limits<T>::max_digits10) << x;
    return o.str();
}

template<typename T>
std::string show(std::expected<rt::Root<T>, rt::error> const &r)
{
    if (!r) return std::format("error {}", r.error());
    return show(r->value) + " +- " + show(r->abs_error);
}

// A function with its first two derivatives, the brackets and the guesses it is solved from.
template<typename T>
struct Case
{
    char const *name;
    std::array<T, 3> (*f)(T const &);
    std::vector<std::array<double, 2>> brackets;
    std::vector<double> guesses;
};

template<typename T>
std::vector<Case<T>> cases()
{
    using std::atan;
    using std::cos;
    using std::exp;
    using std::log;
    using std::sin;
    using std::tan;
    return {
        { "x^2-2", [](T const &x) { return std::array<T, 3>{ T(x * x - T(2)), T(T(2) * x), T(2) }; },
            { { 1, 2 }, { 2, 1 }, { 0, 10 }, { -1, 1e6 } }, { 1, 10, 1e-3, 0, -3 } },
        { "(x-1)^3",
            [](T const &x) {
                T const d = T(x - T(1));
                return std::array<T, 3>{ T(d * d * d), T(T(3) * d * d), T(T(6) * d) };
            },
            { { 0, 3 }, { -10, 1.5 } }, { 2, 0 } },
        { "exp(x)-2",
            [](T const &x) {
                T const e = T(exp(x));
                return std::array<T, 3>{ T(e - T(2)), e, e };
            },
            { { 0, 1 }, { -5, 5 }, { -700, 700 } }, { 0, 5, -20 } },
        { "x^20-1",
            [](T const &x) {
                T p = T(1);
                for (int i = 0; i < 19; ++i) p = T(p * x);
                return std::array<T, 3>{ T(p * x - T(1)), T(T(20) * p), T(T(380) * T(p / x)) };
            },
            { { 0, 5 }, { 0.5, 1.5 } }, { 0.5, 2 } },
        { "log(x/1e-12)",
            [](T const &x) { return std::array<T, 3>{ T(log(T(x / T(1e-12)))), T(T(1) / x), T(T(-1) / T(x * x)) }; },
            { { 1e-20, 1 } }, { 1e-6 } },
        { "tan",
            [](T const &x) {
                T const t = T(tan(x));
                return std::array<T, 3>{ t, T(T(1) + t * t), T(T(2) * t * T(T(1) + t * t)) };
            },
            { { 1, 2 }, { -1, 1 } }, { 1.4, 0.5 } },
        { "x", [](T const &x) { return std::array<T, 3>{ x, T(1), T(0) }; }, { { -1, 2 }, { -1e-300, 1 } }, { 1, -2 } },
        { "x-1e-30", [](T const &x) { return std::array<T, 3>{ T(x - T(1e-30)), T(1), T(0) }; }, { { -1, 1 }, { 0, 1 } },
            { 1 } },
        { "1e300(x-1)", [](T const &x) { return std::array<T, 3>{ T(T(1e300) * T(x - T(1))), T(1e300), T(0) }; },
            { { 0, 1e10 }, { -1e308, 1e308 } }, { 1e7, 0 } },
        { "step", [](T const &x) { return std::array<T, 3>{ x < T(0.3) ? T(-1) : T(1), T(0), T(0) }; }, { { 0, 1 } },
            { 0 } },
        { "cos(x)-x", [](T const &x) { return std::array<T, 3>{ T(cos(x) - x), T(-sin(x) - T(1)), T(-cos(x)) }; },
            { { 0, 1 }, { -2, 3 } }, { 0, 1, 100 } },
        { "x^3-2x-5",
            [](T const &x) {
                return std::array<T, 3>{ T(x * x * x - T(2) * x - T(5)), T(T(3) * x * x - T(2)), T(T(6) * x) };
            },
            { { 2, 3 }, { -10, 10 } }, { 2, -2, 0 } },
        { "atan",
            [](T const &x) {
                T const q = T(T(1) + x * x);
                return std::array<T, 3>{ T(atan(x)), T(T(1) / q), T(T(-2) * x / T(q * q)) };
            },
            { { -1, 3 } }, { 3, 1, 0.5 } },
        { "1/x-1",
            [](T const &x) { return std::array<T, 3>{ T(T(1) / x - T(1)), T(T(-1) / T(x * x)), T(T(2) / T(x * x * x)) }; },
            { { 0.5, 3 }, { -1, 2 } }, { 0.5, 3 } },
        { "x^2", [](T const &x) { return std::array<T, 3>{ T(x * x), T(T(2) * x), T(2) }; }, { { -1, 2 } },
            { 1, -1e-200 } },
        { "sin(x)/x-0.5",
            [](T const &x) {
                T const s = x == T(0) ? T(1) : T(sin(x) / x);
                return std::array<T, 3>{ T(s - T(0.5)), T(0), T(0) };
            },
            { { -1, 3 }, { 0, 3 } }, { 1 } },
        { "1e-300x-1e-310", [](T const &x) { return std::array<T, 3>{ T(T(1e-300) * x - T(1e-310)), T(1e-300), T(0) }; },
            { { -1, 1 } }, { 1 } },
    };
}

template<typename T>
std::vector<rt::Criteria<T>> criteria()
{
    return { {}, { .rel_tol = T(1e-8) }, { .abs_tol = T(1e-6), .rel_tol = T(0) },
        { .abs_tol = T(1e300), .rel_tol = T(1) }, { .abs_tol = T(0), .rel_tol = T(0.25) } };
}

template<typename M, typename T>
void on_bracket(char const *method, Case<T> const &f, std::array<double, 2> const &b, std::size_t c,
    rt::Criteria<T> const &cr)
{
    int n = 0;
    auto const counted = [&](T const &x) {
        ++n;
        return f.f(x);
    };
    auto const r = rt::solve<M>(counted, rt::Bracket<T>{ T(b[0]), T(b[1]) }, cr);
    std::cout << f.name << " [" << b[0] << ", " << b[1] << "] " << method << " c" << c << ": " << show(r)
              << " n=" << n << '\n';
}

template<typename M, typename T>
void from_guess(char const *method, Case<T> const &f, double g, std::size_t c, rt::Criteria<T> const &cr)
{
    int n = 0;
    auto const counted = [&](T const &x) {
        ++n;
        return f.f(x);
    };
    auto const r = rt::solve<M>(counted, T(g), cr);
    std::cout << f.name << " from " << g << " " << method << " c" << c << ": " << show(r) << " n=" << n << '\n';
}

// The separate-derivative forms.
template<typename T>
void separate(Case<T> const &f, std::array<double, 2> const &b, double g)
{
    int n = 0;
    auto const v = [&](T const &x) {
        ++n;
        return f.f(x)[0];
    };
    auto const d = [&](T const &x) { return f.f(x)[1]; };
    auto const d2 = [&](T const &x) { return f.f(x)[2]; };
    auto const r = rt::solve<rt::Newton>(v, d, rt::Bracket<T>{ T(b[0]), T(b[1]) });
    std::cout << f.name << " separate Newton bracket: " << show(r) << " n=" << n << '\n';
    n = 0;
    auto const s = rt::solve<rt::Halley>(v, d, d2, T(g));
    std::cout << f.name << " separate Halley guess: " << show(s) << " n=" << n << '\n';
}

template<typename T>
void battery(char const *type)
{
    std::cout << "=== " << type << '\n';
    double const max = double((std::numeric_limits<T>::max)());
    auto const cs = criteria<T>();
    for (Case<T> const &f : cases<T>()) {
        for (std::size_t c = 0; c < cs.size(); ++c) {
            for (auto const &b : f.brackets) {
                if (!(std::abs(b[0]) <= max) || !(std::abs(b[1]) <= max)) continue; // float: no 1e308
                on_bracket<rt::Bisection>("Bisection", f, b, c, cs[c]);
                on_bracket<rt::Brent>("Brent", f, b, c, cs[c]);
                on_bracket<rt::Toms748>("Toms748", f, b, c, cs[c]);
                on_bracket<rt::AndersonBjorck>("AndersonBjorck", f, b, c, cs[c]);
                on_bracket<rt::Ridders>("Ridders", f, b, c, cs[c]);
                on_bracket<rt::RegulaFalsi>("RegulaFalsi", f, b, c, cs[c]);
                on_bracket<rt::Illinois>("Illinois", f, b, c, cs[c]);
                on_bracket<rt::Newton>("Newton", f, b, c, cs[c]);
                on_bracket<rt::Halley>("Halley", f, b, c, cs[c]);
            }
            for (double g : f.guesses) {
                from_guess<rt::Newton>("Newton", f, g, c, cs[c]);
                from_guess<rt::Halley>("Halley", f, g, c, cs[c]);
                from_guess<rt::Secant>("Secant", f, g, c, cs[c]);
                from_guess<rt::Steffensen>("Steffensen", f, g, c, cs[c]);
            }
        }
        separate(f, f.brackets[0], f.guesses[0]);
    }
}

int main()
{
    battery<float>("float");
    battery<double>("double");
    battery<long double>("long double");
    battery<bin50>("cpp_bin_float_50");
    battery<dec50>("cpp_dec_float_50");
}
