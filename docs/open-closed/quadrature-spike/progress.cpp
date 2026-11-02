// Spike (draft 5): the iterative protocol and observer(), against the header steps.py writes.
#include <numerixx/quadrature.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace qd = nxx::quadrature;

namespace user {

    // A method of one's own, iterative: the trapezoidal rule, halving the step. It computes and reports; integrate()
    // owns the loop and the test against the criteria; observer() reports its steps. Its own limits: ready from
    // level 3, exhausted at level 20, the attainable error a floor from its own sum of |terms| (Criteria's recipe).
    // fn only in [a, b] (a < b, from integrate()): the midpoint as a/2 + b/2, a width that overflows refused, a
    // point a + i h clamped to b. It refines only while its next step (b - a) / 2^k is exact (in the subnormal range a
    // step rounds, and the sums would be of the wrong width).
    struct Trapezoid
    {
        template<nxx::real T>
        struct State
        {
            T a, b, h, sum, abs_sum, previous;
            int level;
            qd::Step<T> step;
        };

        // Whether level k's step (b - a) / 2^k is exact: in radix 2 it is until it is subnormal.
        template<nxx::real T>
        static constexpr bool exact_step(T const &width, int k)
        {
            T const n = T(1 << k);
            T const h = T(width / n);
            return h > T(0) && T(h * n) == width;
        }

        template<nxx::real T>
        static constexpr void report(State<T> &s)
        {
            T const value = T(s.h * s.sum);
            T const floor = T(T(50) * T(std::numeric_limits<T>::epsilon()) * T(s.h * s.abs_sum));
            T const change = nxx::math::abs(T(value - s.previous));
            s.step = qd::Step<T>{ .current = qd::Integral<T>{ value, change > floor ? change : floor },
                .attainable = floor, .ready = s.level >= 3, .exhausted = s.level >= 20 || !exact_step(T(s.b - s.a), s.level + 1) };
            s.previous = value;
        }

        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, qd::Failure<T>> start(F &fn, T const &a, T const &b) const
        {
            qd::Failure<T> const invalid{ qd::error::non_finite, std::nullopt };
            T const width = T(b - a);
            if (!nxx::math::isfinite(width)) return std::unexpected(invalid);
            if (!exact_step(width, 1)) return std::unexpected(qd::Failure<T>{ qd::error::not_converged, std::nullopt });
            T const fa = T(std::invoke(fn, a));
            T const fb = T(std::invoke(fn, b));
            T const mid = T(std::invoke(fn, T(T(a / T(2)) + T(b / T(2)))));
            State<T> s{ a, b, T(width / T(2)), T(T(T(fa + fb) / T(2)) + mid),
                T(T(T(nxx::math::abs(fa) + nxx::math::abs(fb)) / T(2)) + nxx::math::abs(mid)),
                T(T(width / T(2)) * T(fa + fb)), 1, {} };
            report(s);
            if (!nxx::math::isfinite(s.step.current->value)) return std::unexpected(invalid);
            return s;
        }

        template<nxx::real T, typename F>
        constexpr std::expected<void, qd::Failure<T>> advance(F &fn, State<T> &s) const
        {
            ++s.level;
            int const n = 1 << s.level;
            s.h = T(T(s.b - s.a) / T(n));
            for (int i = 1; i < n; i += 2) {
                T x = T(s.a + T(T(i) * s.h));
                if (s.b < x) x = s.b;
                T const y = T(std::invoke(fn, x));
                s.sum = T(s.sum + y);
                s.abs_sum = T(s.abs_sum + nxx::math::abs(y));
            }
            report(s);
            if (!nxx::math::isfinite(s.step.current->value))
                return std::unexpected(qd::Failure<T>{ qd::error::non_finite, std::nullopt });
            return {};
        }

        template<nxx::real T>
        constexpr std::expected<void, qd::Failure<T>> settle(State<T> &) const { return {}; }
    };
    static_assert(qd::iterative_method<Trapezoid, double>);

    // A scripted method: its steps come from a list, so that every branch of integrate()'s loop can be driven. settle
    // replaces the step at index settle_at by `settled`.
    struct Scripted
    {
        std::array<qd::Step<double>, 4> steps;
        int settle_at = -1;
        qd::Step<double> settled{};

        struct State
        {
            int i;
            qd::Step<double> step;
        };

        template<typename F>
        constexpr std::expected<State, qd::Failure<double>> start(F &, double const &, double const &) const
        { return State{ 0, steps[0] }; }

        template<typename F>
        constexpr std::expected<void, qd::Failure<double>> advance(F &, State &s) const
        {
            // Never called on an exhausted step (the loop's promise): a call there would be a non_finite here.
            if (s.step.exhausted) return std::unexpected(qd::Failure<double>{ qd::error::non_finite, std::nullopt });
            s.step = steps[static_cast<std::size_t>(++s.i)];
            return {};
        }

        constexpr std::expected<void, qd::Failure<double>> settle(State &s) const
        {
            if (s.i == settle_at) s.step = settled;
            return {};
        }
    };
    static_assert(qd::iterative_method<Scripted, double>);

} // namespace user

constexpr double cubic(double x) { return x * x * x - x; }
constexpr double one(double) { return 1.0; }
constexpr auto quiet = [](auto const &) {};
using Integral = qd::Integral<double>;
using Step = qd::Step<double>;

// The library's iterative methods, observed or not, give the same result, in constant expressions too: both observed
// forms, the method named (Gauss-Kronrod by default) or given as a value, reversed ends included.
static_assert(qd::iterative_method<qd::Romberg, double> && qd::iterative_method<qd::GaussKronrod<15>, double>);
static_assert(qd::integrate(qd::observer(quiet), cubic, -1.0, 2.0) == qd::integrate(cubic, -1.0, 2.0));
static_assert(qd::integrate<qd::Romberg>(qd::observer(quiet), cubic, 2.0, -1.0, { .rel_tol = 1e-10 })
              == qd::integrate<qd::Romberg>(cubic, 2.0, -1.0, { .rel_tol = 1e-10 }));
static_assert(qd::integrate(qd::Romberg{}, qd::observer(quiet), cubic, 2.0, -1.0, { .rel_tol = 1e-10 })
              == qd::integrate(qd::Romberg{}, cubic, 2.0, -1.0, { .rel_tol = 1e-10 }));
static_assert(qd::integrate(user::Trapezoid{}, qd::observer(quiet), cubic, -1.0, 2.0, { .rel_tol = 1e-6 })
              == qd::integrate<user::Trapezoid>(cubic, -1.0, 2.0, { .rel_tol = 1e-6 }));
// An observer on what does not step is rejected softly, in both forms.
template<typename M>
concept observable = requires { qd::integrate<M>(qd::observer(quiet), cubic, 0.0, 1.0); }
                     && requires { qd::integrate(M{}, qd::observer(quiet), cubic, 0.0, 1.0); };
static_assert(observable<qd::Romberg> && observable<qd::GaussKronrod<15>> && !observable<qd::GaussLegendre<7>>);
// The calls of today, with a generic lambda as fn, still resolve: no overload calls an argument on trial.
constexpr auto generic_square = [](auto const &x) { return x * x; };
static_assert(qd::integrate<qd::Romberg>(generic_square, 0.0, 1.0, { .rel_tol = 1e-12 }).has_value());
static_assert(qd::integrate(qd::Romberg{}, generic_square, 0.0, 1.0, { .rel_tol = 1e-12 }).has_value());
static_assert(qd::integrate(generic_square, 0.0, 1.0).has_value());

// A user's iterative method against the closed form: the trapezoidal rule on x^3 - x over [-1, 2] errs by exactly
// (h^2 / 12) (f'(2) - f'(-1)) = 3 h^2 / 4 (Euler-Maclaurin ends after one term: f''' is constant); its change between
// levels is 9 h^2 / 4, at most 1e-6 * 9/4 first at level 12, h = 3 / 4096, where every sum is exact in double.
static_assert(qd::integrate<user::Trapezoid>(cubic, -1.0, 2.0, { .rel_tol = 1e-6 })->value == 2.25 + 27.0 / 67108864.0);

// The loop's branches, by the scripted method (each in a constant expression).
constexpr Integral two{ 2.0, 1e-3 };
constexpr Step pass{ .current = two };                     // meets rel_tol 1e-2 on 2
constexpr Step miss{ .current = Integral{ 2.0, 1.0 } };    // does not
constexpr qd::Criteria<double> loose{ .rel_tol = 1e-2 };
constexpr qd::Failure<double> non_finite{ qd::error::non_finite, std::nullopt };
// Not ready: it goes on; ready and passing: success.
static_assert(qd::integrate(user::Scripted{ { Step{ .current = two, .ready = false }, pass, miss, miss } }, one, 0.0,
                  1.0, loose)
              == two);
// Passing but untrusted: not_converged with it as the best.
static_assert(qd::integrate(user::Scripted{ { Step{ .current = two, .trusted = false }, miss, miss, miss } }, one, 0.0,
                  1.0, loose)
              == std::unexpected(qd::Failure<double>{ qd::error::not_converged, two }));
// Exhausted and failing: not_converged with the best; without one where current is empty; advance never called there.
static_assert(qd::integrate(user::Scripted{ { miss, Step{ .current = Integral{ 2.0, 1.0 }, .exhausted = true }, miss,
                                miss } },
                  one, 0.0, 1.0, loose)
              == std::unexpected(qd::Failure<double>{ qd::error::not_converged, Integral{ 2.0, 1.0 } }));
static_assert(qd::integrate(user::Scripted{ { miss, Step{ .exhausted = true }, miss, miss } }, one, 0.0, 1.0, loose)
              == std::unexpected(qd::Failure<double>{ qd::error::not_converged, std::nullopt }));
// Passing before settle, failing after: it goes on from the settled step (here to a success at the next step).
static_assert(qd::integrate(user::Scripted{ { pass, pass, miss, miss }, 0, miss }, one, 0.0, 1.0, loose) == two);
// A current that is not finite, or an attainable error that is not finite or negative: non_finite, before any
// arithmetic on it (a NaN would not be a constant expression).
static_assert(qd::integrate(user::Scripted{ { Step{ .current = Integral{ std::numeric_limits<double>::quiet_NaN(), 0.0 } },
                                miss, miss, miss } },
                  one, 0.0, 1.0, { .rel_tol = 0.0 })
              == std::unexpected(non_finite));
static_assert(qd::integrate(user::Scripted{ { Step{ .current = two, .attainable = std::numeric_limits<double>::infinity() },
                                miss, miss, miss } },
                  one, 0.0, 1.0, loose)
              == std::unexpected(non_finite));
static_assert(qd::integrate(user::Scripted{ { Step{ .current = two, .attainable = -1.0 }, miss, miss, miss } }, one,
                  0.0, 1.0, loose)
              == std::unexpected(non_finite));

int main()
{
    // Gauss-Kronrod on a peak of width 1e-3 at 0.3 over [0, 1]: the user counts evaluations in fn, steps in the observer.
    std::size_t evaluations = 0;
    auto const peak = [&evaluations](double x) {
        ++evaluations;
        double const d = x - 0.3;
        return 1.0 / (d * d + 1e-6);
    };
    struct Row
    {
        std::size_t evaluations;
        Step step;
    };
    std::vector<Row> rows;
    auto const record = [&](Step const &s) { rows.push_back(Row{ evaluations, s }); };
    auto const observed = qd::integrate<qd::GaussKronrod<15>>(qd::observer(record), peak, 0.0, 1.0, { .rel_tol = 1e-10 });
    std::size_t const observed_calls = evaluations;
    evaluations = 0;
    auto const plain = qd::integrate(peak, 0.0, 1.0, { .rel_tol = 1e-10 });
    bool exact_counts = true;
    for (std::size_t i = 0; i < rows.size(); ++i) exact_counts = exact_counts && rows[i].evaluations == 150 + 30 * i;
    std::printf("GaussKronrod<15> on a peak: %zu steps; same result as unobserved: %s; evaluations 150 + 30 i at step "
                "i: %s; %zu in all (unobserved: %zu)\n",
        rows.size(), observed == plain ? "yes" : "NO", exact_counts ? "yes" : "NO", observed_calls, evaluations);
    for (std::size_t i = 0; i < rows.size(); i += rows.size() / 5) {
        Step const &s = rows[i].step;
        if (s.current)
            std::printf("  step %2zu  evaluations %4zu  value %.15f  abs_error %.3g  attainable %.3g\n", i,
                rows[i].evaluations, s.current->value, s.current->abs_error, s.attainable);
    }

    // Romberg on e^x over [0, 1], every level: a mutable lambda counts its calls.
    auto const exp_fn = [](double x) { return std::exp(x); };
    std::printf("Romberg on e^x:\n");
    auto const r = qd::integrate<qd::Romberg>(qd::observer([level = 0](Step const &s) mutable {
        ++level;
        if (s.current)
            std::printf("  level %2d  evaluations %5d  value %.17g  abs_error %.3g  ready %d\n", level, (1 << level) + 1,
                s.current->value, s.current->abs_error, s.ready);
    }), exp_fn, 0.0, 1.0, { .rel_tol = 1e-12 });
    std::printf("  result: %.17g +- %.3g (e - 1 = 1.7182818284590452)\n", r->value, r->abs_error);

    // An lvalue observer is kept by reference: what it records is there afterwards.
    struct Recorder
    {
        std::vector<Step> steps;
        void operator()(Step const &s) { steps.push_back(s); }
    } recorder;
    auto const t = qd::integrate(user::Trapezoid{}, qd::observer(recorder), exp_fn, 0.0, 1.0, { .rel_tol = 1e-10 });
    std::printf("Trapezoid (a user's method) on e^x: %zu steps recorded; result %.17g +- %.3g\n",
        recorder.steps.size(), t->value, t->abs_error);
}
