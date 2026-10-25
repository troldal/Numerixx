# One demo program per module, built with the library's warnings as errors.
function(numerixx_add_demo module)
    set(target numerixx_${module}_demo)
    add_executable(${target} ${module}_demo.cpp)
    target_link_libraries(${target} PRIVATE numerixx::numerixx)
    # libstdc++ on MinGW keeps the Windows console part of std::print in libstdc++exp.
    if (MINGW AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_link_libraries(${target} PRIVATE stdc++exp)
    endif ()
    if (MSVC)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive-)
    else ()
        target_compile_options(${target} PRIVATE
                -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Werror)
    endif ()

    # Builds and runs the demo: under Emscripten through node (CMAKE_CROSSCOMPILING_EMULATOR); natively it is empty.
    add_custom_target(run_${module}_demo
            COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:${target}>
            DEPENDS ${target}
            USES_TERMINAL)
endfunction()

numerixx_add_demo(deriv)
numerixx_add_demo(roots)
