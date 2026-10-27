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

add_executable(numerixx_tests test_math.cpp test_deriv.cpp test_stencils.cpp test_callables.cpp test_ridders.cpp
        test_orders.cpp test_onesided.cpp test_deriv_of.cpp test_roots_bracketing.cpp test_roots_open.cpp
        test_roots_numerical.cpp test_roots_toms_anderson.cpp test_roots_search.cpp test_roots_classic.cpp)
target_link_libraries(numerixx_tests PRIVATE numerixx::numerixx doctest::doctest_with_main boost_multiprecision)
if (MSVC)
    target_compile_options(numerixx_tests PRIVATE /W4 /WX /permissive-)
else ()
    target_compile_options(numerixx_tests PRIVATE
            -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Werror)
endif ()

add_test(NAME numerixx_tests COMMAND numerixx_tests)

# The low-precision warning of a derivative above the first in float. Each case is an object library outside the normal
# build, built by its test, since two of them must fail: with warnings as errors, "warns_step" and "warns_default" (one
# per overload of derivative()) fail with the warning's text in the output, while "silenced" (a pragma at the calls)
# and "macro" (NUMERIXX_NO_PRECISION_WARNING) build.
foreach (case IN ITEMS warns_step warns_default silenced macro)
    set(target numerixx_precision_warning_${case})
    add_library(${target} OBJECT EXCLUDE_FROM_ALL precision_warning.cpp)
    target_link_libraries(${target} PRIVATE numerixx::numerixx)
    if (MSVC)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive-)
    else ()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Werror)
    endif ()
    add_test(NAME ${target} COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${target} --config $<CONFIG>)
endforeach ()
target_compile_definitions(numerixx_precision_warning_warns_step PRIVATE NUMERIXX_WARNING_STEP_ONLY)
target_compile_definitions(numerixx_precision_warning_warns_default PRIVATE NUMERIXX_WARNING_DEFAULT_ONLY)
target_compile_definitions(numerixx_precision_warning_silenced PRIVATE NUMERIXX_WARNING_SILENCED)
target_compile_definitions(numerixx_precision_warning_macro PRIVATE NUMERIXX_NO_PRECISION_WARNING)
set_tests_properties(numerixx_precision_warning_warns_step numerixx_precision_warning_warns_default
        PROPERTIES PASS_REGULAR_EXPRESSION "keep few correct digits")
# The four build the same tree, so they must not run concurrently under ctest -j.
set_tests_properties(numerixx_precision_warning_warns_step numerixx_precision_warning_warns_default
        numerixx_precision_warning_silenced numerixx_precision_warning_macro
        PROPERTIES RESOURCE_LOCK numerixx_build_tree)
