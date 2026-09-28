// A consumer of the multiprecision adapter.
#include <numerixx/adapters/multiprecision.hpp>
#include <numerixx/numerixx.hpp>

#include <cstdio>

int main()
{
    using mp50         = boost::multiprecision::cpp_bin_float_50;
    const mp50 residue = mp50 { 1 } / 3 * 3 - 1;
    if (abs(residue) > mp50 { 1e-45 }) return 1;
    std::printf("Numerixx %s: multiprecision adapter OK\n", NUMERIXX_VERSION_STRING);
    return 0;
}
