// Sidebar 7.1 (Chapter 7, section 7.5): for a class-type T, fn must return T, not an expression. With Boost's
// expression templates (et_on, the default of cpp_dec_float_50) x * x - 2 is an expression object that refers to its
// operands; a fn returning it does not compile (nxx::detail::value_of in math.hpp). The same body with -> T, or with
// T(...) around it, is accepted. Built only when the target boost_multiprecision exists (PLAN.md, decision A2).

#include <numerixx/deriv.hpp>

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <functional>
#include <print>
#include <type_traits>

using nxx::deriv::derivative;

// clang-format off
// <et-rule>
using dec50 = boost::multiprecision::cpp_dec_float_50;

// x * x - 2 is an expression object, not a dec50: it
// refers to x, which dies when bare returns. Only the
// checks below use bare
[[maybe_unused]] auto const bare =
    [](dec50 x) { return x * x - 2; };
// the cures: -> dec50, or dec50(...)
auto const typed = [](dec50 const &x) -> dec50 {
    return x * x - 2;
};
auto const wrapped = [](dec50 const &x) {
    return dec50(x * x - 2);
};

// whether derivative(fn, x) compiles
template<typename F>
concept differentiable =
    requires(F fn, dec50 x) { derivative(fn, x); };

static_assert(!differentiable<decltype(bare)>);
static_assert(differentiable<decltype(typed)>);
static_assert(differentiable<decltype(wrapped)>);
// </et-rule>
// clang-format on

// The claim of the text: bare's result is one of Boost's expression types, not dec50.
static_assert(boost::multiprecision::is_number_expression<
              std::invoke_result_t<decltype(bare) const &, dec50 const &>>::value);
static_assert(std::same_as<std::invoke_result_t<decltype(typed) const &, dec50 const &>, dec50>);

int main()
{
    // d/dx (x^2 - 2) = 2x = 3 at x = 1.5. Central<2> is exact on a quadratic apart from rounding, about
    // eps / h with eps about 1e-49 and the default step h about eps^(1/3) |x|, so about 1e-33 relative; 1e-30 leaves
    // a margin of three digits.
    dec50 const x = dec50(3) / 2;
    auto const a = derivative(typed, x);
    auto const b = derivative(wrapped, x);
    dec50 const tol = dec50(3) * dec50("1e-30");
    bool ok = a.has_value() && b.has_value();
    ok = ok && abs(dec50(a->value - 3)) <= tol && abs(dec50(b->value - 3)) <= tol;
    if (ok) std::println("-> dec50: {}", a->value.str(40));
    if (!ok) std::println("FAILED: a claim of the text does not hold");
    return ok ? 0 : 1;
}
