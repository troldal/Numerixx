// Compile-time guard (DESIGN §3.6): the umbrella header alone, instantiating nothing. structural.compile_time.umbrella
// rebuilds this TU and fails on GCC when it takes longer than 2 s.
#include <numerixx/numerixx.hpp>
