# Numerixx

A header-only C++23 library of numerical methods, in the spirit of GSL: small algorithms that combine, work with any
real type (including Boost.Multiprecision) and run in constant expressions. Numerixx itself never throws or catches:
failures are returned as `std::expected`, and an exception thrown by your callback or your number type passes through
unchanged.

Numerixx 2 is in development on the `dev` branch and currently provides numerical derivatives. Numerixx 1.x is at the
tag `v1.1.0-legacy`.

```cpp
#include <numerixx/deriv.hpp>
#include <cmath>

auto d = nxx::deriv::derivative([](double x) { return std::sin(x); }, 1.0); // central difference
if (d) { /* *d is cos(1) to about 1e-10 */ }

static_assert(*nxx::deriv::derivative([](double x) { return x * x; }, 3.0, 0.5) == 6.0);
```

## Using it

With [CPM](https://github.com/cpm-cmake/CPM.cmake):

```cmake
CPMAddPackage("gh:troldal/Numerixx#dev")
target_link_libraries(app PRIVATE numerixx::numerixx)
```

## License

MIT
