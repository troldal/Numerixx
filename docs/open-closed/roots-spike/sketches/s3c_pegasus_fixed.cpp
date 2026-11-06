// Sketch 3c: Pegasus as a user must write it to converge like the library's methods: the secant's quotient formed
// first (no overflow of f * dx), and each new point kept tol inside the ends, so that the bracket closes from both sides.
#include <numerixx/roots.hpp>
#include <cmath>
#include <print>

namespace rt = nxx::roots;

namespace user {

    struct Pegasus
    {
        template<nxx::real T>
        struct State
        {
            T x0, f0, g0; // retained point, fn there, scaled value
            T x1, f1;     // latest point
            rt::Step<T> step;
        };

        template<nxx::real T>
        static constexpr bool opposite(T const &a, T const &b) { return (a < T(0)) != (b < T(0)); }

        template<nxx::real T>
        static constexpr rt::Step<T> report(State<T> const &s)
        {
            bool const x0_low = s.x0 < s.x1;
            T const point = nxx::math::abs(s.f1) <= nxx::math::abs(s.f0) ? s.x1 : s.x0;
            return x0_low ? rt::Step<T>{ s.x0, s.f0, s.x1, s.f1, point } : rt::Step<T>{ s.x1, s.f1, s.x0, s.f0, point };
        }

        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
        {
            State<T> st{ s.lower, s.at_lower[0], s.at_lower[0], s.upper, s.at_upper[0], {} };
            st.step = report(st);
            return st;
        }

        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &t) const
        {
            T const lo = s.step.lower, hi = s.step.upper;
            T const mid = T(T(lo / T(2)) + T(hi / T(2)));
            // The secant through (x1, f1) and (x0, g0): q in (0, 1) since f1 and g0 have opposite signs.
            T const q = T(s.f1 / T(s.f1 - s.g0));
            T x = T(s.x1 - T(q * T(s.x1 - s.x0)));
            if (!nxx::math::isfinite(x)) x = mid;
            // Kept tol(point) inside the ends, so that a point converging from one side lands on the other.
            T const margin = t.at(s.step.point);
            if (T(T(hi - lo) / T(4)) <= margin) x = mid;
            else if (x < T(lo + margin)) x = T(lo + margin);
            else if (x > T(hi - margin)) x = T(hi - margin);
            if (!(lo < x && x < hi)) {
                s.step.exhausted = true;
                return {};
            }
            T const fx = T(std::invoke(fn, x));
            if (!nxx::math::isfinite(fx)) return std::unexpected(rt::error::non_finite);
            if (fx == T(0)) {
                s.step = rt::Step<T>{ x, fx, x, fx, x };
                return {};
            }
            if (opposite(fx, s.f1)) {
                s.x0 = s.x1;
                s.f0 = s.f1;
                s.g0 = s.f1;
            } else {
                s.g0 = T(s.g0 * T(s.f1 / T(s.f1 + fx)));
            }
            s.x1 = x;
            s.f1 = fx;
            s.step = report(s);
            return {};
        }
    };

} // namespace user

template<typename M, typename F>
void run(char const *name, F f, rt::Bracket<double> b, double exact, rt::Criteria<double> c = {})
{
    int calls = 0;
    auto counted = [&](double x) { ++calls; return f(x); };
    auto const r = rt::solve<M>(counted, b, c);
    if (r) std::println("{:22} {:.17g} +- {:.3g}  |x - exact| = {:.3g}  calls {:4}  {}", name, r->value, r->abs_error,
                        std::abs(r->value - exact), calls, std::abs(r->value - exact) <= r->abs_error ? "ok" : "BOUND FAILS");
    else std::println("{:22} error {} (calls {})", name, r.error(), calls);
}

int main()
{
    auto const sq = [](double x) { return x * x - 2.0; };
    auto const cubic = [](double x) { return x * x * x - 2.0 * x - 5.0; };
    auto const tiny = [](double x) { return 1e-200 * (x * x * x - 2.0 * x - 5.0); };
    auto const x20 = [](double x) { return std::pow(x, 20) - 1.0; };
    auto const flat = [](double x) { return std::exp(-1.0 / (x * x)) - 0.5; };
    auto const lin = [](double x) { return x - 1.5; };
    auto const atan_ = [](double x) { return std::atan(x) - 0.5; };
    auto const cosx = [](double x) { return std::cos(x) - x; };
    struct Case { char const *name; double (*f)(double); rt::Bracket<double> b; double exact; };
    Case const cases[] = {
        { "sqrt2 [1,2]", +sq, { 1.0, 2.0 }, std::sqrt(2.0) },
        { "cubic [2,3]", +cubic, { 2.0, 3.0 }, 2.0945514815423265 },
        { "tiny cubic [2,3]", +tiny, { 2.0, 3.0 }, 2.0945514815423265 },
        { "x^20-1 [0,1.5]", +x20, { 0.0, 1.5 }, 1.0 },
        { "flat [0.1,10]", +flat, { 0.1, 10.0 }, 1.0 / std::sqrt(std::log(2.0)) },
        { "cos x - x [0,1]", +cosx, { 0.0, 1.0 }, 0.73908513321516064 },
        { "x-1.5 [+-1e308]", +lin, { -1e308, 1e308 }, 1.5 },
        { "atan-0.5 [+-1e308]", +atan_, { -1e308, 1e308 }, std::tan(0.5) },
    };
    for (auto const &k : cases) {
        std::println("-- {}", k.name);
        run<user::Pegasus>("Pegasus (fixed)", k.f, k.b, k.exact);
        run<rt::AndersonBjorck>("AndersonBjorck (lib)", k.f, k.b, k.exact);
        run<rt::Illinois>("Illinois (lib)", k.f, k.b, k.exact);
        run<rt::Brent>("Brent (lib)", k.f, k.b, k.exact);
    }
    std::println("-- cubic, rel_tol 1e-6");
    run<user::Pegasus>("Pegasus (fixed)", cubic, { 2.0, 3.0 }, 2.0945514815423265, { .rel_tol = 1e-6 });
    run<rt::AndersonBjorck>("AndersonBjorck (lib)", cubic, { 2.0, 3.0 }, 2.0945514815423265, { .rel_tol = 1e-6 });
}
