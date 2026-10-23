# The book "Numerical Derivatives with Numerixx" (docs/deriv-book/CLAUDE.md):
# - every code example in code/ is an executable and a test, part of the normal build, so a listing cannot rot;
# - numerixx_book_data runs the generator, which writes data/ (values.tex and the figures' CSV files) in the source
#   tree, from the library itself;
# - numerixx_book regenerates data/, builds the PDF with latexmk into this build directory and runs the checks.
# The data and PDF targets are on demand only and need a native toolchain (and latexmk); not in CI.

# The tests' warning flags; MSVC is also true for clang-cl, whose -Wall would mean -Weverything.
if (MSVC)
    set(book_warnings /W4 /WX /permissive-)
else ()
    set(book_warnings -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Werror)
endif ()

file(GLOB book_examples CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/code/*.cpp)
foreach (source IN LISTS book_examples)
    get_filename_component(name ${source} NAME_WE)
    set(target numerixx_book_${name})
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE numerixx::numerixx)
    target_compile_options(${target} PRIVATE ${book_warnings})
    # libstdc++ on MinGW keeps the Windows console part of std::print in libstdc++exp.
    if (MINGW AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_link_libraries(${target} PRIVATE stdc++exp)
    endif ()
    add_test(NAME ${target} COMMAND ${target})
endforeach ()

if (NOT EMSCRIPTEN)
    add_executable(numerixx_book_generator EXCLUDE_FROM_ALL tools/generate_data.cpp)
    target_link_libraries(numerixx_book_generator PRIVATE numerixx::numerixx)
    target_compile_options(numerixx_book_generator PRIVATE ${book_warnings})
    add_custom_target(numerixx_book_data
            COMMAND numerixx_book_generator ${CMAKE_CURRENT_SOURCE_DIR}/data
            COMMENT "Generating docs/deriv-book/data from the library"
            VERBATIM)

    find_program(NUMERIXX_LATEXMK latexmk)
    find_package(Python3 COMPONENTS Interpreter)
    if (NUMERIXX_LATEXMK AND Python3_FOUND)
        set(book_pdf ${CMAKE_CURRENT_BINARY_DIR}/pdf)
        add_custom_target(numerixx_book
                COMMAND ${CMAKE_COMMAND} -E make_directory ${book_pdf}
                COMMAND ${NUMERIXX_LATEXMK} -interaction=nonstopmode -halt-on-error -outdir=${book_pdf} main.tex
                COMMAND ${Python3_EXECUTABLE} tools/check.py ${book_pdf}/main.log
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                DEPENDS numerixx_book_data
                COMMENT "Building ${book_pdf}/main.pdf"
                USES_TERMINAL
                VERBATIM)
    else ()
        message(STATUS "latexmk or Python not found: no target numerixx_book")
    endif ()
endif ()
