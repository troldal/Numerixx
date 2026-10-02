// Compile-fail (DESIGN §6.4, §9.1): Newton with a function and a derivative that report different error types. The
// common-cause rule rejects it and says how to make them agree.
#include <numerixx/roots.hpp>

#include <expected>

namespace r = nxx::roots;

namespace
{
    enum class f_error { domain };
    enum class df_error { domain };
}    // namespace

int main()
{
    auto f = [](double x) -> std::expected<double, f_error> { return x * x - 2.0; };
#ifdef NUMERIXX_CF_CONTROL
    auto df = [](double x) -> std::expected<double, f_error> { return 2.0 * x; };
#else
    auto df = [](double x) -> std::expected<double, df_error> { return 2.0 * x; };
#endif
    (void)r::newton {}.with_derivative(df)(f, 1.0);
}
