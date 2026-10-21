include(${PROJECT_SOURCE_DIR}/cmake/CPM.cmake)

CPMAddPackage(
        NAME doctest
        VERSION 2.5.3
        GITHUB_REPOSITORY doctest/doctest
        SYSTEM YES
        OPTIONS "DOCTEST_WITH_TESTS OFF" "DOCTEST_NO_INSTALL ON")

# Boost.Multiprecision supplies the user-defined real types the tests run on. Its standalone mode is header-only and
# needs only Boost.Config.
CPMAddPackage(
        NAME boost_config
        VERSION 1.92.0
        URL https://github.com/boostorg/config/archive/refs/tags/boost-1.92.0.tar.gz
        URL_HASH SHA256=b4171037f13373203ba79cbc141d612982052283e696a315185ab5bea46102a0
        DOWNLOAD_ONLY YES)
CPMAddPackage(
        NAME boost_multiprecision
        VERSION 1.92.0
        URL https://github.com/boostorg/multiprecision/archive/refs/tags/boost-1.92.0.tar.gz
        URL_HASH SHA256=9da997843edd802f9d5e25d53872a7676330f31a52a5e47a418d1c94ba3b01d0
        DOWNLOAD_ONLY YES)
add_library(boost_multiprecision INTERFACE)
target_include_directories(boost_multiprecision SYSTEM INTERFACE
        ${boost_multiprecision_SOURCE_DIR}/include ${boost_config_SOURCE_DIR}/include)
target_compile_definitions(boost_multiprecision INTERFACE BOOST_MP_STANDALONE)

add_executable(numerixx_tests test_math.cpp test_deriv.cpp test_stencils.cpp test_callables.cpp test_ridders.cpp)
target_link_libraries(numerixx_tests PRIVATE numerixx::numerixx doctest::doctest_with_main boost_multiprecision)
if (MSVC)
    target_compile_options(numerixx_tests PRIVATE /W4 /WX /permissive-)
else ()
    target_compile_options(numerixx_tests PRIVATE
            -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Werror)
endif ()

add_test(NAME numerixx_tests COMMAND numerixx_tests)
