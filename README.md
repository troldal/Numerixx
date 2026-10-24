# Numerixx

A header-only C++23 library of numerical methods, in the spirit of GSL: small algorithms that combine, work with any
real type (including Boost.Multiprecision) and run in constant expressions. Numerixx itself never throws or catches:
failures are returned as `std::expected`, and an exception thrown by your callback or your number type passes through
unchanged.

Numerixx 2 is in development on the `dev` branch and currently provides numerical derivatives: central and one-sided
finite differences of orders up to 10, and Ridders' extrapolation with an error estimate, central or one-sided, for
the first to the fourth derivative. Numerixx 1.x is at the tag `v1.1.0-legacy`.

```cpp
#include <numerixx/deriv.hpp>
#include <cmath>

auto d = nxx::deriv::derivative([](double x) { return std::sin(x); }, 1.0); // central difference
if (d) { double const slope = *d; /* cos(1) to about 1e-10 */ }

using nxx::deriv::Central, nxx::deriv::Forward, nxx::deriv::Ridders;
auto d4 = nxx::deriv::derivative<Central<4>>([](double x) { return std::sin(x); }, 1.0); // about 1e-13
auto df = nxx::deriv::derivative<Forward<2>>([](double x) { return std::sqrt(x); }, 0.5); // never evaluates below x
auto dr = nxx::deriv::derivative<Ridders>([](double x) { return std::sin(x); }, 1.0);     // adapts its step
if (dr) { /* dr->value is cos(1) to about 1e-14, *dr->abs_error estimates its error */ }
auto d2 = nxx::deriv::derivative<Ridders, 2>([](double x) { return std::sin(x); }, 1.0);  // f'', about 1e-12
auto de = nxx::deriv::derivative<nxx::deriv::ForwardRidders>([](double x) { return std::sqrt(x); }, 0.0625);
// x and above only, with an error estimate: for a domain edge, or the end of a correlation's range
auto d3 = nxx::deriv::derivative<Central<2>, nxx::deriv::order::third>([](double x) { return std::sin(x); }, 1.0);

static_assert(*nxx::deriv::derivative([](double x) { return x * x; }, 3.0, 0.5) == 6.0);
```

## Using it

With [CPM](https://github.com/cpm-cmake/CPM.cmake):

```cmake
CPMAddPackage("gh:troldal/Numerixx#dev")
target_link_libraries(app PRIVATE numerixx::numerixx)
```

The API documentation is generated with Doxygen: build the target `numerixx_docs`, or run `doxygen docs/Doxyfile` in
the project root.

## License

MIT, except third-party data: the Shomate coefficients of nitrogen in `docs/deriv-book/code/shomate.hpp`
are NIST Standard Reference Data (NIST Chemistry WebBook, SRD 69), (c) U.S. Secretary of Commerce, reproduced with
attribution and not covered by the MIT licence.
