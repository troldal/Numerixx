// Prints every result of a battery bit for bit, so that two builds (old and new header) can be diffed.
#include <numerixx/quadrature.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <iomanip>

namespace qd = nxx::quadrature;
using bin50 = boost::multiprecision::cpp_bin_float_50;
using dec50 = boost::multiprecision::cpp_dec_float_50;

template<typename T> std::string show(T const &x)
{
    std::ostringstream o;
    if constexpr (std::is_floating_point_v<T>) o << std::hexfloat << x;
    else o << std::setprecision(std::numeric_limits<T>::max_digits10) << x;
    return o.str();
}
template<typename T> void out(char const *what, std::expected<T, qd::error> const &r)
{ std::printf("%s %s\n", what, r ? show(*r).c_str() : "error"); }
template<typename T> void out(char const *what, std::expected<qd::Integral<T>, qd::Failure<T>> const &r)
{
    if (r) std::printf("%s ok %s %s\n", what, show(r->value).c_str(), show(r->abs_error).c_str());
    else if (r.error().best) std::printf("%s fail %d best %s %s\n", what, int(r.error().reason), show(r.error().best->value).c_str(), show(r.error().best->abs_error).c_str());
    else std::printf("%s fail %d\n", what, int(r.error().reason));
}

template<typename T> void battery(char const *type)
{
    using std::exp; using std::sqrt; using std::log; using std::sin; using std::abs;
    using boost::multiprecision::exp; using boost::multiprecision::sqrt; using boost::multiprecision::log;
    using boost::multiprecision::sin; using boost::multiprecision::abs;
    auto f_exp = [](T const &x) -> T { return T(exp(x)); };
    auto f_sqrt = [](T const &x) -> T { return T(sqrt(x)); };
    auto f_inv_sqrt = [](T const &x) -> T { return T(T(1) / T(sqrt(x))); };
    auto f_log = [](T const &x) -> T { return T(log(x)); };
    auto f_kink = [](T const &x) -> T { return T(abs(T(x - T(1) / T(3)))); };
    auto f_peak = [](T const &x) -> T { T d = T(x - T(3) / T(10)); return T(T(1) / T(d * d + T(1) / T(1000000))); };
    auto f_sin = [](T const &x) -> T { return T(sin(T(T(50) * x))); };
    auto f_nan = [](T const &x) -> T { return x < T(1) / T(2) ? T(std::numeric_limits<double>::quiet_NaN()) : x; };
    std::array<std::pair<T, T>, 6> const ends{ { { T(0), T(1) }, { T(1), T(0) }, { T(-2), T(3) }, { T(3), T(-2) },
        { T(1) / T(1000), T(7) }, { T(7), T(1) / T(1000) } } };
    char buf[128];
    int i = 0;
    auto run = [&](char const *name, auto const &f) {
        for (auto const &[a, b] : ends) {
            ++i;
            std::snprintf(buf, sizeof buf, "%s %s #%d GL1", type, name, i); out(buf, qd::integrate<qd::GaussLegendre<1>>(f, a, b));
            std::snprintf(buf, sizeof buf, "%s %s #%d GL7", type, name, i); out(buf, qd::integrate<qd::GaussLegendre<7>>(f, a, b));
            std::snprintf(buf, sizeof buf, "%s %s #%d GL20", type, name, i); out(buf, qd::integrate<qd::GaussLegendre<20>>(f, a, b));
            std::snprintf(buf, sizeof buf, "%s %s #%d GL64", type, name, i); out(buf, qd::integrate<qd::GaussLegendre<64>>(f, a, b));
            std::snprintf(buf, sizeof buf, "%s %s #%d GK", type, name, i); out(buf, qd::integrate(f, a, b));
            std::snprintf(buf, sizeof buf, "%s %s #%d GK-12", type, name, i); out(buf, qd::integrate<qd::GaussKronrod<15>>(f, a, b, { .rel_tol = 1e-12 }));
            std::snprintf(buf, sizeof buf, "%s %s #%d GK-0", type, name, i); out(buf, qd::integrate(f, a, b, { .abs_tol = 0, .rel_tol = 0 }));
            std::snprintf(buf, sizeof buf, "%s %s #%d R", type, name, i); out(buf, qd::integrate<qd::Romberg>(f, a, b));
            std::snprintf(buf, sizeof buf, "%s %s #%d R-10", type, name, i); out(buf, qd::integrate<qd::Romberg>(f, a, b, { .rel_tol = 1e-10 }));
            std::snprintf(buf, sizeof buf, "%s %s #%d R-big", type, name, i); out(buf, qd::integrate<qd::Romberg>(f, a, b, { .abs_tol = 1e-3, .rel_tol = 5 }));
        }
    };
    run("exp", f_exp); run("sqrt", f_sqrt); run("inv_sqrt", f_inv_sqrt); run("log", f_log); run("kink", f_kink);
    run("peak", f_peak); run("sin", f_sin); run("nan", f_nan);
    // Edges.
    T const max = (std::numeric_limits<T>::max)();
    auto one = [](T const &) -> T { return T(1); };
    auto quarter = [](T const &) -> T { return T(1) / T(4); };
    out("edge 1 on [-max,max] GK", qd::integrate(one, T(-max), max));
    out("edge 1/4 on [-max,max] GK", qd::integrate(quarter, T(-max), max));
    out("edge 1/4 on [-max,max] GL7", qd::integrate<qd::GaussLegendre<7>>(quarter, T(-max), max));
    out("edge 1/4 on [-max,max] R", qd::integrate<qd::Romberg>(quarter, T(-max), max));
    out("edge a==b", qd::integrate(one, T(2), T(2)));
    out("edge nan end", qd::integrate(one, T(0), T(std::numeric_limits<double>::quiet_NaN())));
    out("edge neg tol", qd::integrate(one, T(0), T(1), { .abs_tol = -1 }));
}

int main()
{
    battery<float>("float");
    battery<double>("double");
    battery<long double>("ldouble");
    battery<bin50>("bin50");
    battery<dec50>("dec50");
    constexpr auto c1 = qd::integrate<qd::GaussLegendre<7>>([](double x) { return x * x; }, 0.0, 3.0);
    constexpr auto c2 = qd::integrate([](double x) { return x * x * x - x; }, -1.0, 2.0);
    constexpr auto c3 = qd::integrate<qd::Romberg>([](double x) { return 1.0 / (1.0 + x * x); }, 0.0, 1.0);
    out("constexpr GL7", c1); out("constexpr GK", c2); out("constexpr R", c3);
}
