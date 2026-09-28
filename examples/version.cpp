// Prints the Numerixx version. The module examples arrive with their phases (DESIGN §10.3).
//
// std::printf rather than std::println: with MinGW's libstdc++, std::print needs -lstdc++exp at link time.
#include <numerixx/numerixx.hpp>

#include <cstdio>

int main()
{
    std::printf("Numerixx %d.%d.%d\n", nxx::version.major, nxx::version.minor, nxx::version.patch);
    return 0;
}
