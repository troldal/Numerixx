// Does the compiler accept a deletion reason, `= delete("...")` (C++26)? neg_msvc.ps1 compiles this with cl, which
// rejects the syntax (MSVC 19.51).
struct S { void f(int) = delete("reason text here"); };
int main() { S{}.f(1); }
