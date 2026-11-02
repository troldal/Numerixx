// Spike (draft 3): the iterative protocol and observed(), against the header steps.py writes.
#include <numerixx/quadrature.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

namespace qd = nxx::quadrature;

namespace user {

    // A method of one's own, iterative: the trapezoidal rule, halving the step. It computes and reports; integrate()
    // owns the loop, the test against the criteria and (through observed()) the progress. Its own limits: ready from
    // level 3, exhausted at level 20; its attainable error a floor from its own sum of |terms| (Criteria's recipe).
    struct Trapezoid
    {
        template<nxx::real T>
        struct State
        {
            T a, b, sum, abs_sum, previous;
            int level;
            qd::Step<T> step;
        };

        template<nxx::real T>
        static constexpr void report(State<T> &s)
        {
            T const h = T(T(s.b - s.a) / T(1 << s.level));
            T const value = T(h * s.sum);
            T const floor = T(T(50) * T(std::numeric_limits<T>::epsilon()) * T(h * s.abs_sum));
            T const change = nxx::math::abs(T(value - s.previous));
            s.step = qd::Step<T>{ qd::Integral<T>{ value, change > floor ? change : floor }, floor, s.level >= 3,
                s.level >= 20, true };
            s.previous = value;
        }

        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, qd::Failure<T>> start(F &fn, T const &a, T const &b) const
        {
            T const fa = T(std::invoke(fn, a));
            T const fb = T(std::invoke(fn, b));
            T const mid = T(std::invoke(fn, T(T(a + b) / T(2))));
            State<T> s{ a, b, T(T(T(fa + fb) / T(2)) + mid), T(T(T(nxx::math::abs(fa) + nxx::math::abs(fb)) / T(2)) + nxx::math::abs(mid)),
                T(T(T(b - a) / T(2)) * T(fa + fb)), 1,
                qd::Step<T>{ std::nullopt, T(0), false, false, true } };
            report(s);
            if (!nxx::math::isfinite(s.step.current->value)) return std::unexpected(qd::Failure<T>{ qd::error::non_finite, std::nullopt });
            return s;
        }

        template<nxx::real T, typename F>
        constexpr std::expected<void, qd::Failure<T>> advance(F &fn, State<T> &s) const
        {
            ++s.level;
            int const n = 1 << s.level;
            T const h = T(T(s.b - s.a) / T(n));
            for (int i = 1; i < n; i += 2) {
                T const y = T(std::invoke(fn, T(s.a + T(T(i) * h))));
                s.sum = T(s.sum + y);
                s.abs_sum = T(s.abs_sum + nxx::math::abs(y));
            }
            report(s);
            if (!nxx::math::isfinite(s.step.current->value)) return std::unexpected(qd::Failure<T>{ qd::error::non_finite, std::nullopt });
            return {};
        }
    };
    static_assert(qd::iterative_method<Trapezoid, double>);

} // namespace user

constexpr double cubic(double x) { return x * x * x - x; }
constexpr auto quiet = [](auto const &) {};

// The library's iterative methods, observed or not, give the same result, in constant expressions too.
static_assert(qd::iterative_method<qd::Romberg, double> && qd::iterative_method<qd::GaussKronrod<15>, double>);
static_assert(!qd::estimating_method<qd::Romberg, double>); // one form each: they step
static_assert(qd::integrate(qd::observed(qd::GaussKronrod<15>{}, quiet), cubic, -1.0, 2.0)
              == qd::integrate(cubic, -1.0, 2.0));
static_assert(qd::integrate(qd::observed(qd::Romberg{}, quiet), cubic, 2.0, -1.0, { .rel_tol = 1e-10 })
              == qd::integrate<qd::Romberg>(cubic, 2.0, -1.0, { .rel_tol = 1e-10 }));
// A user's iterative method: x^3 - x on [-1, 2] = 2.25; the trapezoid's error (b - a) h^2 f'' / 12 shrinks by 4 per level.
static_assert(nxx::math::abs(qd::integrate<user::Trapezoid>(cubic, -1.0, 2.0, { .rel_tol = 1e-6 })->value - 2.25) < 1e-5);

int main()
{
    // Progress of Gauss-Kronrod on a peak of width 1e-3 at 0.3 over [0, 1].
    std::size_t independent = 0;
    auto const peak = [&independent](double x) {
        ++independent;
        double const d = x - 0.3;
        return 1.0 / (d * d + 1e-6);
    };
    std::vector<qd::Progress<double>> seen;
    auto const record = [&seen](qd::Progress<double> const &p) { seen.push_back(p); };
    auto const observed = qd::integrate(qd::observed(qd::GaussKronrod<15>{}, record), peak, 0.0, 1.0, { .rel_tol = 1e-10 });
    std::size_t const observed_calls = independent;
    independent = 0;
    auto const plain = qd::integrate(peak, 0.0, 1.0, { .rel_tol = 1e-10 });
    std::printf("GaussKronrod<15> on a peak: %zu progress reports; same result as unobserved: %s; evaluations reported %zu, "
                "counted by fn %zu (unobserved run: %zu)\n",
        seen.size(), observed == plain ? "yes" : "NO", seen.back().evaluations, observed_calls, independent);
    for (std::size_t i = 0; i < seen.size(); i += seen.size() / 6) {
        qd::Progress<double> const &p = seen[i];
        std::printf("  iteration %3d  evaluations %5zu  value %.15f  abs_error %.3g  attainable %.3g  exhausted %d\n",
            p.iteration, p.evaluations, p.step.current->value, p.step.current->abs_error, p.step.attainable,
            p.step.exhausted);
    }
    std::printf("  result: %.15f +- %.3g\n", observed->value, observed->abs_error);

    // Romberg on e^x over [0, 1], every level.
    auto const exp_fn = [](double x) { return std::exp(x); };
    std::printf("Romberg on e^x:\n");
    auto const r = qd::integrate(qd::observed(qd::Romberg{}, [](qd::Progress<double> const &p) {
        std::printf("  level %2d  evaluations %5zu  value %.17g  abs_error %.3g  ready %d\n", p.iteration + 1,
            p.evaluations, p.step.current->value, p.step.current->abs_error, p.step.ready);
    }), exp_fn, 0.0, 1.0, { .rel_tol = 1e-12 });
    std::printf("  result: %.17g +- %.3g (e - 1 = 1.7182818284590452)\n", r->value, r->abs_error);

    // A user's iterative method, observed.
    std::printf("Trapezoid (a user's method) on e^x:\n");
    auto const t = qd::integrate(qd::observed(user::Trapezoid{}, [](qd::Progress<double> const &p) {
        if (p.iteration % 4 == 0)
            std::printf("  level %2d  evaluations %6zu  value %.17g  abs_error %.3g\n", p.iteration + 1, p.evaluations,
                p.step.current->value, p.step.current->abs_error);
    }), exp_fn, 0.0, 1.0, { .rel_tol = 1e-10 });
    if (t) std::printf("  result: %.17g +- %.3g\n", t->value, t->abs_error);
    else std::printf("  failed: reason %d\n", static_cast<int>(t.error().reason));
}
