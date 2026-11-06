# Spike: the roots method protocol (`docs/roots/PHASE7_API.md`)

Evidence for the proposal, not part of the build. `assemble.py` writes the proposal's headers from today's
`include/numerixx/` and three chunks; the programs here are built against them. From the project root:

```sh
# The headers: roots.hpp with the protocol, math.hpp with nxx::observer, quadrature.hpp using it.
python3 docs/open-closed/roots-spike/assemble.py build/roots-spike/include

# The existing tests, unchanged, against the spike's headers (the roots tests, and the quadrature tests that observe).
printf '#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN\n#include <doctest/doctest.h>\n' > build/roots-spike/main.cpp
docs/open-closed/roots-spike/build.sh build/roots-spike/include build/roots-spike/obj \
    tests/test_roots_*.cpp tests/test_quadrature_protocol.cpp tests/test_quadrature_examples.cpp \
    tests/test_quadrature_romberg.cpp tests/test_deriv_of.cpp build/roots-spike/main.cpp
g++ build/roots-spike/obj/*.o -o build/roots-spike/tests && build/roots-spike/tests

# Bit for bit: the battery against today's header and the spike's, every result and evaluation count compared.
SYS="-isystem $(echo ~/.cache/cpm/boost_multiprecision/*/include) -isystem $(echo ~/.cache/cpm/boost_config/*/include)"
g++ -std=c++23 -O1 -DBOOST_MP_STANDALONE -Iinclude $SYS docs/open-closed/roots-spike/dump.cpp -o build/roots-spike/dump_old
g++ -std=c++23 -O1 -DBOOST_MP_STANDALONE -Ibuild/roots-spike/include -Iinclude $SYS \
    docs/open-closed/roots-spike/dump.cpp -o build/roots-spike/dump_new
build/roots-spike/dump_old > build/roots-spike/old.txt
build/roots-spike/dump_new > build/roots-spike/new.txt
cmp build/roots-spike/old.txt build/roots-spike/new.txt
```

- `chunk_public.hpp`: `Tolerance` (promoted from `detail`), `Start`, `Step`, `values_of`.
- `chunk_states.hpp`: each library method's state and one step, transcribed from today's `detail` loops, the stop
  test removed (it is `solve`'s).
- `chunk_methods.hpp`: the eleven method types, the four concepts, the observer wrapper, the loop (`detail::run`),
  validation (`on_bracket`, `on_guess`) and the forms of `solve`.
- `build.sh`: compiles sources against a chosen include directory with the tests' warnings, two at a time.
- `dump.cpp`: the battery (17 functions, their brackets and guesses, five tolerances, every method, the
  separate-derivative forms; float, double, long double, `cpp_bin_float_50`, `cpp_dec_float_50`).
- `sketches/`: the usage review's programs.

Built with GCC 16.2 in the Codespace, the only toolchain there; the phase itself builds and tests on all five.
