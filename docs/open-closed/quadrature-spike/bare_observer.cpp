// Spike (draft 5): why observer(o) wraps the callable. Built as is, against the header steps.py writes, it compiles and
// prints the integral. With -DBARE it adds what the wrapper avoids, an overload integrate<M>(o, fn, a, b[, criteria])
// for a bare callable o, told from the other forms by whether o can be called with a Step<T>; and then a call that
// compiles today fails: to answer that question for a generic lambda the compiler instantiates its body, and the
// integrand's body does not accept a Step (GCC 13: "no match for 'operator*' (operand types are 'const
// nxx::quadrature::Step<double>' and 'const nxx::quadrature::Step<double>')"). Ordering the overload's constraints so
// that fn and the method are checked first avoids it for today's overloads, but leaves correctness to that order.
#include <numerixx/quadrature.hpp>

#include <cstdio>

#ifdef BARE
namespace nxx::quadrature {
    template<typename M = GaussKronrod<15>, typename O, real T, typename F>
    requires std::default_initializable<M> && std::invocable<O &, Step<T> const &> && nxx::detail::callable<F, T>
             && iterative_method<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(O &&o, F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {})
    { return nxx::quadrature::integrate(M{}, observer(std::forward<O>(o)), fn, a, b, criteria); }
} // namespace nxx::quadrature
#endif

namespace qd = nxx::quadrature;

int main()
{
    auto fn = [](auto const &x) { return x * x; }; // generic: a natural way to write an integrand
    auto const r = qd::integrate<qd::Romberg>(fn, 0.0, 1.0, { .rel_tol = 1e-12 });
    std::printf("%.17g\n", r->value);
}
