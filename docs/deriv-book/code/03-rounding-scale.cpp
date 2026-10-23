// E4 (Chapter 3, section 3.5): when the default step is wrong. The default step eps^(1/3) max(|x|, 1) assumes that
// f varies on the scale max(|x|, 1). Two cases where it does not, with the same fix, a step eps^(1/3) w built from
// the scale w on which f varies:
// - narrow-feature: tanh((x - x0) / w) with x0 = 1000 and w = 0.5, at x0 + 0.2 (DS4, F4);
// - small-x: y ln y at y = 1e-6, where the default step exceeds y and samples below zero (DS4b, F4b).
// The parameters are those of the generator (tools/generate_data.cpp, ds4_step_scale and ds4b_step_scale_small).

#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

int main()
{
    using nxx::deriv::derivative;
    using nxx::deriv::Ridders;

    // clang-format off
// <narrow-feature>
// a step from 0 to 1 over a width of about w, at x0 = 1000
double const x0 = 1000.0, w = 0.5;
auto const f = [=](double x) {
    return std::tanh((x - x0) / w);
};
double const x = x0 + 0.2;

auto const d = derivative(f, x); // step eps^(1/3) * |x|
using nxx::math::root_epsilon;
double const cbrt_eps = root_epsilon<3, double>();
auto const d_w = derivative(f, x, cbrt_eps * w);
auto const r = derivative<Ridders>(f, x, w / 4);
// </narrow-feature>

// <small-x>
// y ln y is defined for y > 0 only; the default step is
// eps^(1/3) * 1, far larger than y itself
double const y = 1e-6;
auto const g = [](double t) { return t * std::log(t); };
auto const e = derivative(g, y);           // samples y < 0
auto const e_y = derivative(g, y, cbrt_eps * y);
// </small-x>
    // clang-format on

    if (!d || !d_w || !r || !e_y) return 1;
    double const th = std::tanh((x - x0) / w);
    double const f_exact = (1 - th * th) / w;
    double const g_exact = std::log(y) + 1;
    double const rel_d = std::abs(d->value - f_exact) / f_exact;
    double const rel_w = std::abs(d_w->value - f_exact) / f_exact;
    double const r_actual = std::abs(r->value - f_exact);
    double const rel_r = r_actual / f_exact;
    double const rel_e_y = std::abs(e_y->value - g_exact) / std::abs(g_exact);
    std::println("narrow feature, default step {:.1e}: relative error {:.1e}", cbrt_eps * x, rel_d);
    std::println("narrow feature, step {:.1e}:         relative error {:.1e}", cbrt_eps * w, rel_w);
    std::println("narrow feature, Ridders from w/4: relative error {:.1e}, estimate {:.1e} (relative {:.1e})", rel_r,
        *r->abs_error, *r->abs_error / f_exact);
    std::println("small y, default step {:.1e}: {}", cbrt_eps, e ? "a value" : "non_finite");
    std::println("small y, step {:.1e}: relative error {:.1e}", cbrt_eps * y, rel_e_y);

    // The model predicts the fitted step about 1e4 times better than the default; 100 is a safe margin. Ridders:
    // covered within the tests' factor 6 and accurate to 1e-9 (measured 4.8e-15). The small-y default step fails;
    // the fitted step is accurate to about eps^(2/3) (measured 7.9e-12, asserted 1e-9).
    bool const ok = rel_w * 100 <= rel_d && r_actual <= 6.0 * *r->abs_error && rel_r <= 1e-9 && !e
                    && e.error() == nxx::deriv::error::non_finite && rel_e_y <= 1e-9;
    return ok ? 0 : 1;
}
