// Prints the Numerixx version. quick_tour.cpp shows the spike's modules; the others arrive with their phases (DESIGN §10.3).
//
// With MinGW's libstdc++, std::print needs libstdc++exp at link time, which examples/CMakeLists.txt adds where the
// toolchain needs it.
#include <numerixx/numerixx.hpp>

#include <print>

int main()
{
    std::println("Numerixx {}.{}.{}", nxx::version.major, nxx::version.minor, nxx::version.patch);
    return 0;
}
