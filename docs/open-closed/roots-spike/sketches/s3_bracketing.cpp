// Sketch 3: bracketing methods of one's own, against the public protocol only (no nxx::detail).
#include <numerixx/roots.hpp>
#include <cmath>
#include <print>

namespace rt = nxx::roots;

namespace user {

    // ---- Pegasus (Dowell and Jarratt 1972): false position whose retained value is scaled by f1 / (f1 + f2). ----
    // Naive = true: the sign test a user writes first, fa * fb < 0, which underflows.
    template<bool Naive = false>
    struct Pegasus
    {
        template<nxx::real T>
        struct State
        {
            T x0, f0, g0; // the retained point, fn there, and the scaled value the secant uses
            T x1, f1;     // the latest point
            rt::Step<T> step;
        };

        template<nxx::real T>
        static constexpr bool opposite(T const &a, T const &b)
        {
            if constexpr (Naive) return a * b < T(0);
            else return (a < T(0)) != (b < T(0));
        }

        // The report: the two points in order, fn's TRUE values there (not g0), the end with the smaller |f|.
        template<nxx::real T>
        static constexpr rt::Step<T> report(State<T> const &s)
        {
            bool const x0_low = s.x0 < s.x1;
            T const lo = x0_low ? s.x0 : s.x1, flo = x0_low ? s.f0 : s.f1;
            T const hi = x0_low ? s.x1 : s.x0, fhi = x0_low ? s.f1 : s.f0;
            T const point = nxx::math::abs(s.f1) <= nxx::math::abs(s.f0) ? s.x1 : s.x0;
            return rt::Step<T>{ lo, flo, hi, fhi, point };
        }

        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
        {
            State<T> st{ s.lower, s.at_lower[0], s.at_lower[0], s.upper, s.at_upper[0], {} };
            st.step = report(st);
            return st;
        }

        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &) const
        {
            T const lo = s.step.lower, hi = s.step.upper;
            T x = T(s.x1 - T(s.f1 * T(T(s.x1 - s.x0) / T(s.f1 - s.g0))));
            if (!(lo < x && x < hi)) {
                x = T(T(lo / T(2)) + T(hi / T(2))); // overflow-free midpoint
                if (!(lo < x && x < hi)) {
                    s.step.exhausted = true;
                    return {};
                }
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

    // ---- Plain trisection: the bracket cut in three, the third with the sign change kept. ----
    struct Trisection
    {
        template<nxx::real T>
        struct State
        {
            T a, fa, b, fb;
            rt::Step<T> step;
        };

        template<nxx::real T>
        static constexpr rt::Step<T> report(T const &a, T const &fa, T const &b, T const &fb)
        { return rt::Step<T>{ a, fa, b, fb, nxx::math::abs(fa) <= nxx::math::abs(fb) ? a : b }; }

        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
        {
            return State<T>{ s.lower, s.at_lower[0], s.upper, s.at_upper[0],
                             report(s.lower, s.at_lower[0], s.upper, s.at_upper[0]) };
        }

        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &) const
        {
            T const third = T(T(s.b / T(3)) - T(s.a / T(3)));
            T const x1 = T(s.a + third), x2 = T(s.b - third);
            if (!(s.a < x1 && x1 < x2 && x2 < s.b)) {
                s.step.exhausted = true;
                return {};
            }
            auto const negative = [](T const &v) { return v < T(0); };
            T const f1 = T(std::invoke(fn, x1));
            if (!nxx::math::isfinite(f1)) return std::unexpected(rt::error::non_finite);
            if (f1 == T(0)) { s.step = rt::Step<T>{ x1, f1, x1, f1, x1 }; return {}; }
            if (negative(f1) != negative(s.fa)) { s.b = x1; s.fb = f1; }
            else {
                T const f2 = T(std::invoke(fn, x2));
                if (!nxx::math::isfinite(f2)) return std::unexpected(rt::error::non_finite);
                if (f2 == T(0)) { s.step = rt::Step<T>{ x2, f2, x2, f2, x2 }; return {}; }
                if (negative(f2) != negative(f1)) { s.a = x1; s.fa = f1; s.b = x2; s.fb = f2; }
                else { s.a = x2; s.fa = f2; }
            }
            s.step = report(s.a, s.fa, s.b, s.fb);
            return {};
        }
    };

} // namespace user

static_assert(rt::bracketing_method<user::Pegasus<>, double>);
static_assert(rt::bracketing_method<user::Trisection, double>);
static_assert(!rt::open_method<user::Trisection, double>);

template<typename M, typename F>
void run(char const *name, F f, rt::Bracket<double> b, double exact, rt::Criteria<double> c = {})
{
    int calls = 0, steps = 0;
    auto counted = [&](double x) { ++calls; return f(x); };
    auto const r = rt::solve<M>(rt::observer([&](rt::Step<double> const &) {
        if (++steps > 10000) { std::println("  runaway: stopping"); std::exit(9); }
    }), counted, b, c);
    if (r) std::println("{:22} {:.17g} +- {:.3g}  |x - exact| = {:.3g}  calls {:3}  {}", name, r->value, r->abs_error,
                        std::abs(r->value - exact), calls, std::abs(r->value - exact) <= r->abs_error ? "ok" : "BOUND FAILS");
    else std::println("{:22} error {} (calls {})", name, r.error(), calls);
}

int main()
{
    auto const sq = [](double x) { return x * x - 2.0; };
    auto const cubic = [](double x) { return x * x * x - 2.0 * x - 5.0; };
    auto const tiny = [](double x) { return 1e-200 * (x - 1.4142135623730951); };
    auto const x20 = [](double x) { return std::pow(x, 20) - 1.0; };
    auto const flat = [](double x) { return std::exp(-1.0 / (x * x)) - 0.5; }; // root 1/sqrt(ln 2)
    double const flat_root = 1.0 / std::sqrt(std::log(2.0));
    double const s2 = std::sqrt(2.0), c = 2.0945514815423265;

    for (auto const &[name, fn_root] : { std::pair{ "sqrt2", 0 }, { "cubic", 1 }, { "tiny scale", 2 }, { "x^20 - 1", 3 },
                                        { "flat", 4 }, { "huge bracket", 5 } }) {
        std::println("-- {}", name);
        auto go = [&]<typename M>(char const *mname) {
            switch (fn_root) {
            case 0: run<M>(mname, sq, { 1.0, 2.0 }, s2); break;
            case 1: run<M>(mname, cubic, { 2.0, 3.0 }, c); break;
            case 2: run<M>(mname, tiny, { 1.0, 2.0 }, s2); break;
            case 3: run<M>(mname, x20, { 0.0, 1.5 }, 1.0); break;
            case 4: run<M>(mname, flat, { 0.1, 10.0 }, flat_root); break;
            case 5: run<M>(mname, sq, { 0.0, 1e300 }, s2); break;
            }
        };
        go.template operator()<user::Pegasus<>>("Pegasus");
        go.template operator()<user::Pegasus<true>>("Pegasus (naive sign)");
        go.template operator()<user::Trisection>("Trisection");
        go.template operator()<rt::AndersonBjorck>("AndersonBjorck (lib)");
        go.template operator()<rt::Illinois>("Illinois (lib)");
        go.template operator()<rt::Bisection>("Bisection (lib)");
    }
    std::println("-- abs_tol 1e-6 on the cubic");
    run<user::Pegasus<>>("Pegasus", cubic, { 2.0, 3.0 }, c, { .abs_tol = 1e-6, .rel_tol = 0 });
    run<rt::AndersonBjorck>("AndersonBjorck (lib)", cubic, { 2.0, 3.0 }, c, { .abs_tol = 1e-6, .rel_tol = 0 });
    std::println("-- as a value: {}", rt::solve(user::Trisection{}, sq, rt::Bracket{ 1.0, 2.0 })->value);
}
