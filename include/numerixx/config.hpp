// Configuration macros shared by every Numerixx header (DESIGN §5.3).
#pragma once

#include <cassert>

// [[no_unique_address]]: MSVC (and clang-cl) accept the attribute only in their own spelling.
#if defined(_MSC_VER)
#    define NXX_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#else
#    define NXX_NO_UNIQUE_ADDRESS [[no_unique_address]]
#endif

// NXX_DELETE(reason): a deleted overload that states why, where the compiler supports `= delete("reason")` (C++26,
// P2573). Clang 19+ (the floor) accepts it as an extension, silenced inside the macro. GCC 15+ accepts it too, but rejects
// _Pragma inside a declaration, so headers that use NXX_DELETE bracket their body with NXX_BEGIN_HEADER and
// NXX_END_HEADER. Elsewhere (MSVC, GCC 14) the overload is deleted without a reason; the diagnostic names
// the file and line of the deleted declaration, so NXX_DELETE always starts on the declarator's own line.
#if defined(__clang__) && __clang_major__ >= 19
#    define NXX_DELETE(reason)                                                                        \
        _Pragma("clang diagnostic push") _Pragma("clang diagnostic ignored \"-Wc++26-extensions\"") = \
            delete (reason)_Pragma("clang diagnostic pop")
#elif defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 15
#    define NXX_DELETE(reason) = delete (reason)
#else
#    define NXX_DELETE(reason) = delete
#endif

// NXX_BEGIN_HEADER / NXX_END_HEADER bracket the body of every header with arithmetic:
//   - Clang (also clang-cl and em++): no floating-point contraction in library code (on the targets listed below). Clang contracts a * b +
//   c into a
//     fused multiply-add by default; on FMA hardware (ARM64, x86 with -march=haswell) and in constant folding that
//     rounds once instead of twice, so a solver would take a different path on different platforms. With contraction
//     off, the library's arithmetic is the same everywhere (DESIGN §6.1): given the same values of f, a solver takes
//     the same path. The setting is saved and restored, so the user's code keeps its own.
//   - GCC 15+: the -Wc++26-extensions diagnostic of NXX_DELETE is silenced (see above). Contraction is NOT turned off:
//     GCC contracts a * b + c by default in C++, ISO mode included (only ISO C defaults to -ffp-contract=off), and it
//     has no pragma for a region (#pragma GCC optimize would stop inlining). On FMA targets (AArch64, x86 with -mfma
//     or -march=haswell) a GCC build may therefore differ in the last ulp; build with -ffp-contract=off for
//     bit-identical results there (DESIGN §5.3). The presets build baseline x86-64, which has no FMA.
//   - Other Clang targets: #pragma float_control is supported where the target has strict floating point (x86,
//     x86-64, AArch64, RISC-V, PowerPC, SystemZ; checked with Clang 22); elsewhere (WebAssembly, 32-bit ARM, MIPS,
//     LoongArch, ...) it is ignored with a warning, so the setting cannot be saved and restored there. Contraction is
//     turned off and stays off after the header; WebAssembly has no scalar FMA, so there this only makes constant
//     folding agree with run time.
//   - MSVC: C4459 (a declaration hides a global) is silenced in the library's templates, so a consumer's global named
//     like one of the library's locals (e, n, s, ...) does not break a /W4 /WX build. Consumers that include Numerixx
//     through its CMake targets see its headers as external anyway (DESIGN D24).
#if defined(__clang__) && (defined(__x86_64__) || defined(__i386__) || defined(__aarch64__) || defined(_M_X64) || defined(_M_IX86) || \
                           defined(_M_ARM64) || defined(__riscv) || defined(__powerpc__) || defined(__s390__))
#    define NXX_BEGIN_HEADER _Pragma("float_control(push)") _Pragma("clang fp contract(off)")
#    define NXX_END_HEADER _Pragma("float_control(pop)")
#elif defined(__clang__)
#    define NXX_BEGIN_HEADER _Pragma("clang fp contract(off)")
#    define NXX_END_HEADER
#elif defined(__GNUC__) && __GNUC__ >= 15
#    define NXX_BEGIN_HEADER _Pragma("GCC diagnostic push") _Pragma("GCC diagnostic ignored \"-Wc++26-extensions\"")
#    define NXX_END_HEADER _Pragma("GCC diagnostic pop")
#elif defined(_MSC_VER)
#    define NXX_BEGIN_HEADER __pragma(warning(push)) __pragma(warning(disable : 4459))
#    define NXX_END_HEADER __pragma(warning(pop))
#else
#    define NXX_BEGIN_HEADER
#    define NXX_END_HEADER
#endif

// Contracts of the low-level protocol (DESIGN §3.4): calling init/step with a problem that prepare() did not produce,
// narrowing a bracket at a point outside it. Checked in debug builds; C++26 contracts later.
#define NXX_EXPECTS(cond) assert(cond)
#define NXX_ASSERT(cond) assert(cond)
