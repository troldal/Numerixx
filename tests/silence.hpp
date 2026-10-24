// Silences the low-precision warning of derivative() around a deliberate call, for the tests: float derivatives
// above the first select a deprecated overload. The warning is reported at the call, so the call itself must be
// between the two macros, not inside a template that a silenced line instantiates (GCC and MSVC report it at the line
// in the template).

#pragma once

#if defined(_MSC_VER) && !defined(__clang__)
#define NUMERIXX_SILENCE_BEGIN __pragma(warning(push)) __pragma(warning(disable : 4996))
#define NUMERIXX_SILENCE_END __pragma(warning(pop))
#else
#define NUMERIXX_SILENCE_BEGIN \
    _Pragma("GCC diagnostic push") _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#define NUMERIXX_SILENCE_END _Pragma("GCC diagnostic pop")
#endif
