// Compile-fail (DESIGN §6.6): a solver that does not implement the protocol, here one without estimate(). The facade
// called it anyway, so the call (and std::is_invocable_v) was a hard error inside detail::run ("no matching function
// for call to iterate"); the facade now checks the protocol and a deleted overload names it.
#include <numerixx/roots.hpp>

namespace r = nxx::roots;

namespace
{
#ifdef NUMERIXX_CF_CONTROL
    // A complete solver: bisection's protocol under an id of its own.
    struct own_bisection : r::bisection<>
    {
        static constexpr nxx::algo id = nxx::algo::user_first;
    };
#else
    // bisection's protocol with estimate() hidden: what a solver that lacks it looks like to the facade.
    struct own_bisection : r::bisection<>
    {
        void estimate() const = delete;
    };
#endif
}    // namespace

int main()
{
    auto f = [](double x) { return x * x - 2.0; };
    (void)own_bisection {}(f, { 1.0, 2.0 });
}
