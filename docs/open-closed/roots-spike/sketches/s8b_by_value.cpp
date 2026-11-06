// Sketch 8b: a forgotten & on the state in advance: compiles, and loops forever on a copy. Could the concept see it?
#include <numerixx/roots.hpp>
#include <print>
namespace rt = nxx::roots;

namespace user {
    struct Halves
    {
        template<nxx::real T> struct State { T a, fa, b, fb; rt::Step<T> step; };
        template<nxx::real T, typename F>
        constexpr std::expected<State<T>, rt::error> start(F &, rt::Start<T> const &s, rt::Tolerance<T> const &) const
        { return State<T>{ s.lower, s.at_lower[0], s.upper, s.at_upper[0], { s.lower, s.at_lower[0], s.upper, s.at_upper[0], s.lower } }; }
    };
    struct ByRef : Halves
    {
        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &) const
        {
            T const m = T(T(s.a / T(2)) + T(s.b / T(2)));
            T const fm = T(std::invoke(fn, m));
            if ((fm < T(0)) == (s.fa < T(0))) { s.a = m; s.fa = fm; } else { s.b = m; s.fb = fm; }
            s.step = rt::Step<T>{ s.a, s.fa, s.b, s.fb, s.a };
            return {};
        }
    };
    struct ByValue : Halves // the same body, & forgotten
    {
        template<nxx::real T, typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State<T> s, rt::Tolerance<T> const &t) const
        { return ByRef{}.advance(fn, s, t); }
    };
}

// A clause the concept could add: advance must not accept an rvalue state (so it takes it by non-const reference).
template<typename M, typename T, typename F = T (*)(T const &)>
concept takes_state_by_reference = !requires(M const &m, F &fn, rt::bracketing_state_t<M, T, F> &&st, rt::Tolerance<T> const &t) {
    m.advance(fn, std::move(st), t);
};

static_assert(rt::bracketing_method<user::ByValue, double>); // accepted today
static_assert(takes_state_by_reference<user::ByRef, double>);
static_assert(!takes_state_by_reference<user::ByValue, double>);
static_assert(takes_state_by_reference<rt::Brent, double> && takes_state_by_reference<rt::Toms748, double>
              && takes_state_by_reference<rt::Bisection, double> && takes_state_by_reference<rt::Ridders, double>
              && takes_state_by_reference<rt::Illinois, double>);

int main()
{
    auto const f = [](double x) { return x * x - 2.0; };
    std::println("ByRef: {:.17g}", rt::solve<user::ByRef>(f, rt::Bracket{ 1.0, 2.0 })->value);
    int n = 0;
    (void)rt::solve<user::ByValue>(rt::observer([&](rt::Step<double> const &) {
        if (++n == 1000000) { std::println("ByValue: 1e6 steps, the bracket never changes"); std::exit(0); }
    }), f, rt::Bracket{ 1.0, 2.0 });
}
