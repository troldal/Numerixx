// Self-test of the compile-fail harness (cmake/RunCompileFail.cmake, DESIGN §9.1). The EXPECT text of this case sits in
// a comment on the deleted declaration: the compiler quotes that line in its "declared here" note, but never says it in
// its own message. The harness must reject the case, or a missing deletion reason could pass on the quoted line.
void refuse(int) = delete;    // this reason is visible only in quoted source

int main()
{
#ifdef NUMERIXX_CF_CONTROL
    return 0;
#else
    refuse(1);
#endif
}
