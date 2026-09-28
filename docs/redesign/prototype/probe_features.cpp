#include <version>
#include <cstdio>
#define SHOW(m) std::printf("%-34s %s\n", #m, show_(#m, STR(m)))
#define STR(x) STR2(x)
#define STR2(x) #x
static const char* show_(const char* name, const char* val) { return (val[0]=='_' ) ? "(undefined)" : val; }
int main() {
  SHOW(__cplusplus);
  SHOW(__cpp_deleted_function);
  SHOW(__cpp_explicit_this_parameter);
  SHOW(__cpp_consteval);
  SHOW(__cpp_if_consteval);
  SHOW(__cpp_exceptions);
  SHOW(__cpp_lib_expected);
  SHOW(__cpp_lib_ranges);
  SHOW(__cpp_lib_move_only_function);
  SHOW(__cpp_lib_generator);
  SHOW(__cpp_lib_constexpr_cmath);
  SHOW(__cpp_lib_forward_like);
  SHOW(__cpp_lib_unreachable);
  SHOW(__cpp_lib_print);
  SHOW(__cpp_lib_stacktrace);
  SHOW(__cpp_lib_function_ref);
  SHOW(__cpp_lib_optional);
  SHOW(__cpp_lib_ranges_iota);
  SHOW(__cpp_lib_ranges_zip);
  SHOW(__cpp_static_call_operator);
}
