// TanhSinh on infinite intervals (Phase 4 step 2, docs/quadrature/PHASE4_API.md): the opt-in of integrate(), exp-sinh
// on [a, inf) and (-inf, b], sinh-sinh on (-inf, inf). The edge cases every method shares are in
// test_quadrature_kronrod.cpp's template (an infinite end non_finite for the methods that do not take one), a method of
// one's own that takes infinite ends in test_quadrature_protocol.cpp, the node table in test_quadrature_exp_sinh.cpp.
//
// Integrals are checked against closed forms or mpmath at 60 digits (an exponent or a scale written as a double literal
// taken as that double, exactly) with the honesty check |value - ref| <= abs_error, and within the request where it
// succeeds. Evaluation counts are asserted exactly in double, in total and per level through the observer, for fn of
// + - * / and sqrt (correctly rounded, so bit for bit in any language), against the kept prototype
// (docs/quadrature/tanh-sinh-prototype/es.py, run with cpp=True, which loads the C++ double table es_double.txt); no
// count is asserted where a libm function's last bits could move a decision. How each side ended is read from the
// method's state after a level (detail::TanhSinhState), the prototype's record of the same run.

#include "quadrature_check.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <cstdlib>
#include <expected>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

using namespace quadrature_check;

namespace {

    using tanh_sinh = qd::TanhSinh;
    using gk15 = qd::GaussKronrod<15>;
    using cut_by = qd::detail::cut_by;

    // References (mpmath at 80 digits, printed to 60).
    constexpr char const *sqrt_pi = "1.77245385090551602729816748334114518279754945612238712821381";
    constexpr char const *half_sqrt_pi = "0.886226925452758013649083741670572591398774728061193564106904";
    constexpr char const *pi = "3.14159265358979323846264338327950288419716939937510582097494";
    constexpr char const *half_pi = "1.57079632679489661923132169163975144209858469968755291048747";

    enum class status { ok, not_converged, non_finite };

    template<typename T>
    status status_of(std::expected<qd::Integral<T>, qd::Failure<T>> const &r)
    {
        if (r) return status::ok;
        return r.error().reason == qd::error::not_converged ? status::not_converged : status::non_finite;
    }

    // An observed run in double: its result, and the evaluations counted after each report.
    struct Run
    {
        std::expected<qd::Integral<double>, qd::Failure<double>> result;
        std::vector<int> counts;
        std::vector<qd::Step<double>> steps;
    };

    template<typename F>
    Run run_observed(F fn, double a, double b, qd::Criteria<double> c = {})
    {
        int n = 0;
        std::vector<int> counts;
        std::vector<qd::Step<double>> steps;
        auto const record = [&](qd::Step<double> const &s) {
            counts.push_back(n);
            steps.push_back(s);
        };
        auto const r = qd::integrate<tanh_sinh>(qd::observer(record), Counted<F>{ fn, &n }, a, b, c);
        return { r, counts, steps };
    }

    // The method's state after level `level`, driven as integrate() drives it (no criteria), for the sides' causes.
    template<typename F>
    qd::detail::TanhSinhState<double> state_at(F fn, double a, double b, int level)
    {
        auto s = tanh_sinh{}.start(fn, a, b);
        REQUIRE(s.has_value());
        while (s->k < level) {
            REQUIRE(!s->step.exhausted);
            REQUIRE(tanh_sinh{}.advance(fn, *s).has_value());
        }
        return *s;
    }

    double root_of(double x)
    { return std::sqrt(x); }

    // x^(1/64) by six square roots: powers by sqrt alone, for exact counts.
    double root64(double x)
    {
        for (int k = 0; k < 6; ++k) x = root_of(x);
        return x;
    }

    // sqrt where defined, else NaN (as the prototype's).
    double root_or_nan(double x)
    { return x >= 0.0 ? root_of(x) : dnan; }

    // The level-k tail cases of the prototype (r4_widek3.py, r4_widek.py): NaN on (X1, X2), between level-0 nodes.
    constexpr double tail_x1 = 598066702.4976885; // 2 1.05^400
    constexpr double tail_x2 = 1196133404.995377; // 2 X1
    double level_k_tail(double x)
    { return root_or_nan((x - tail_x1) * (x - tail_x2)) / (x * x) / x; }

    constexpr auto runge = [](double x) { return 1.0 / (1.0 + x * x); };
    constexpr auto inverse_square = [](double x) { return 1.0 / ((1.0 + x) * (1.0 + x)); };
    constexpr auto step_at_2 = [](double x) { return (x < 2.0 ? 1.0 : 0.0) / (1.0 + x * x); };

} // namespace

static_assert(qd::detail::takes_infinite_ends<tanh_sinh>);
static_assert(!qd::detail::takes_infinite_ends<qd::Romberg> && !qd::detail::takes_infinite_ends<gk15>
              && !qd::detail::takes_infinite_ends<qd::GaussLegendre<7>>);
static_assert(qd::detail::takes_infinite_ends<qd::detail::Observed<tanh_sinh, decltype(runge)>>
              && !qd::detail::takes_infinite_ends<qd::detail::Observed<qd::Romberg, decltype(runge)>>);

TEST_CASE("quadrature: infinite: the opt-in")
{
    qd::Failure<double> const invalid{ qd::error::non_finite, std::nullopt };
    int n = 0;
    auto const counted = Counted{ +runge, &n };
    std::pair<double, double> const infinite[] = { { 0.0, dinf }, { -dinf, 0.0 }, { -dinf, dinf }, { dinf, 0.0 },
                                                   { dinf, -dinf }, { dinf, dinf } };
    // Every method but TanhSinh: non_finite, fn never called; the fixed rule takes no infinite end in any form.
    for (auto const &[a, b] : infinite) {
        INFO(a << ", " << b);
        CHECK(qd::integrate(counted, a, b) == std::unexpected(invalid));
        CHECK(qd::integrate<qd::Romberg>(counted, a, b) == std::unexpected(invalid));
        CHECK(qd::integrate<qd::Romberg>(qd::observer([](qd::Step<double> const &) {}), counted, a, b)
              == std::unexpected(invalid));
        CHECK(qd::integrate(gk15{}, counted, a, b) == std::unexpected(invalid));
        CHECK(qd::integrate<qd::GaussLegendre<7>>(counted, a, b) == std::unexpected(qd::error::non_finite));
    }
    CHECK(n == 0);
    // TanhSinh: a NaN end is non_finite, an equal infinity {0, 0}, both without calling fn.
    CHECK(qd::integrate<tanh_sinh>(counted, dnan, dinf) == std::unexpected(invalid));
    CHECK(qd::integrate<tanh_sinh>(counted, -dinf, dnan) == std::unexpected(invalid));
    CHECK(qd::integrate<tanh_sinh>(counted, dinf, dinf) == qd::Integral<double>{ 0.0, 0.0 });
    CHECK(qd::integrate<tanh_sinh>(counted, -dinf, -dinf) == qd::Integral<double>{ 0.0, 0.0 });
    CHECK(n == 0);

    // Observed: the same run, its steps reported; the last report is the result (ends in order).
    std::vector<qd::Step<double>> steps;
    auto const r = qd::integrate<tanh_sinh>(qd::observer([&steps](qd::Step<double> const &s) { steps.push_back(s); }),
                                            runge, -dinf, dinf);
    REQUIRE(r.has_value());
    CHECK(r == qd::integrate<tanh_sinh>(runge, -dinf, dinf));
    CHECK(r == qd::integrate(tanh_sinh{}, qd::observer([](qd::Step<double> const &) {}), runge, -dinf, dinf));
    REQUIRE(steps.size() == 5);
    CHECK(steps.back().current == *r);

    // Reversed ends are the exactly negated result, value and best (the run is on the ends in order).
    auto const negated = [](std::expected<qd::Integral<double>, qd::Failure<double>> x) {
        auto const flip = [](qd::Integral<double> i) {
            i.value = -i.value;
            return i;
        };
        if (x) return decltype(x)(flip(*x));
        qd::Failure<double> f = x.error();
        if (f.best) f.best = flip(*f.best);
        return decltype(x)(std::unexpected(f));
    };
    for (double const tol : { 1.5e-8, 0.0 }) {
        qd::Criteria<double> const c{ .rel_tol = tol };
        CHECK(qd::integrate<tanh_sinh>(inverse_square, dinf, 0.0, c)
              == negated(qd::integrate<tanh_sinh>(inverse_square, 0.0, dinf, c)));
        CHECK(qd::integrate<tanh_sinh>(inverse_square, 0.0, -dinf, c)
              == negated(qd::integrate<tanh_sinh>(inverse_square, -dinf, 0.0, c)));
        CHECK(qd::integrate<tanh_sinh>(runge, dinf, -dinf, c)
              == negated(qd::integrate<tanh_sinh>(runge, -dinf, dinf, c)));
    }
    // A best too: the step at 2 runs to the cap with one.
    auto const capped = qd::integrate<tanh_sinh>(step_at_2, 0.0, dinf);
    REQUIRE(!capped.has_value());
    REQUIRE(capped.error().best.has_value());
    CHECK(qd::integrate<tanh_sinh>(step_at_2, dinf, 0.0) == negated(capped));
}

TEST_CASE("quadrature: infinite: accuracy on S in double")
{
    // PHASE4_API.md F's 26 smooth, endpoint and decay cases on half lines and the line: ok at 1e-6, the default, 1e-12
    // and 0, the honesty check, and the error within the request rel_tol |ref| (where rel_tol > 0: a success's
    // estimate is at most max(rel_tol |value|, 2 F), and 2 F, 100 eps L1, is far below 1e-12 |ref| for every case,
    // L1 / |ref| being at most about 2). References: closed forms (mpmath at 80 digits).
    struct Case
    {
        char const *name;
        double a;
        double b;
        double (*fn)(double);
        char const *ref;
    };
    Case const cases[] = {
        { "e^-x", 0.0, dinf, [](double x) { return std::exp(-x); }, "1" },
        { "x e^-x", 0.0, dinf, [](double x) { return x * std::exp(-x); }, "1" },
        { "x^2 e^-x", 0.0, dinf, [](double x) { return std::pow(x, 2) * std::exp(-x); }, "2" },
        { "x^5 e^-x", 0.0, dinf, [](double x) { return std::pow(x, 5) * std::exp(-x); }, "120" },
        { "x^10 e^-x", 0.0, dinf, [](double x) { return std::pow(x, 10) * std::exp(-x); }, "3628800" },
        { "e^-x^2 [0, inf)", 0.0, dinf, [](double x) { return std::exp(-x * x); }, half_sqrt_pi },
        { "1/(1 + x^2) [0, inf)", 0.0, dinf, [](double x) { return 1.0 / (1.0 + x * x); }, half_pi },
        { "1/(1 + x)^2", 0.0, dinf, [](double x) { return 1.0 / ((1.0 + x) * (1.0 + x)); }, "1" },
        { "x^-0.5 e^-x", 0.0, dinf, [](double x) { return std::exp(-x) / std::sqrt(x); }, sqrt_pi },
        { "log x e^-x", 0.0, dinf, [](double x) { return std::log(x) * std::exp(-x); },
          "-0.577215664901532860606512090082402431042159335939923598805767" },
        { "e^-x cos x", 0.0, dinf, [](double x) { return std::exp(-x) * std::cos(x); }, "0.5" },
        // Gamma(1 - p), p the double 0.9.
        { "x^-0.9 e^-x", 0.0, dinf, [](double x) { return std::exp(-x) * std::pow(x, -0.9); },
          "9.51350769866873403823051630302669684791590698432035078689421" },
        // (pi / 1.5) / sin(pi / 1.5).
        { "1/(1 + x^1.5)", 0.0, dinf, [](double x) { return 1.0 / (1.0 + std::pow(x, 1.5)); },
          "2.41839915231229046745877101018954097637875499745698743409318" },
        { "e^-sqrt x", 0.0, dinf, [](double x) { return std::exp(-std::sqrt(x)); }, "2" },
        { "x^-1.5 [1, inf)", 1.0, dinf, [](double x) { return std::pow(x, -1.5); }, "2" },
        { "x^-2 [1, inf)", 1.0, dinf, [](double x) { return std::pow(x, -2.0); }, "1" },
        { "log x / x^2 [1, inf)", 1.0, dinf, [](double x) { return std::log(x) / (x * x); }, "1" },
        { "e^x (-inf, 0]", -dinf, 0.0, [](double x) { return std::exp(x); }, "1" },
        { "x^-2 (-inf, -1]", -dinf, -1.0, [](double x) { return std::pow(-x, -2.0); }, "1" },
        { "e^-x^2", -dinf, dinf, [](double x) { return std::exp(-x * x); }, sqrt_pi },
        { "1/(1 + x^2)", -dinf, dinf, [](double x) { return 1.0 / (1.0 + x * x); }, pi },
        { "sech x", -dinf, dinf, [](double x) { return 1.0 / std::cosh(x); }, pi },
        // pi / sqrt 2.
        { "1/(1 + x^4)", -dinf, dinf, [](double x) { return 1.0 / (1.0 + std::pow(x, 4)); },
          "2.22144146907918312350794049503034684930731084468784511154270" },
        { "x^2 e^-x^2", -dinf, dinf, [](double x) { return x * x * std::exp(-x * x); }, half_sqrt_pi },
        // sqrt(pi) e^(-1/4).
        { "e^-x^2 cos x", -dinf, dinf, [](double x) { return std::exp(-x * x) * std::cos(x); },
          "1.38038844704314297477341524672559127427077246556221079845025" },
        // sqrt(pi) Gamma(p - 1/2) / Gamma(p), p the double 0.6.
        { "1/(1 + x^2)^0.6", -dinf, dinf, [](double x) { return std::pow(1.0 + x * x, -0.6); },
          "11.3230869752157559548809108990077018370939290765966283668498" },
    };
    for (Case const &c : cases) {
        for (double const tol : { 1e-6, 1.4901161193847656e-8, 1e-12, 0.0 }) {
            INFO(c.name << ", rel_tol " << tol);
            auto const r = qd::integrate<tanh_sinh>(c.fn, c.a, c.b, { .rel_tol = tol });
            REQUIRE(r.has_value());
            double const error = error_of(r->value, c.ref);
            CHECK(error <= r->abs_error); // the honesty check
            if (tol > 0.0) CHECK(error <= tol * std::abs(static_cast<double>(bin100(c.ref))));
        }
    }
}

TEST_CASE("quadrature: infinite: accuracy in the 50-digit types")
{
    // Two or three cases per map at 1e-25 and 1e-45, in cpp_bin_float_50 with and without expression templates and
    // cpp_dec_float_50: ok, the honesty check, and the error within rel_tol |ref| (2 F, 100 eps L1, is at most about
    // 1e-47 relative for each, cpp_dec_float_50's eps being 1e-49 and L1 about |ref|). References: closed forms.
    auto const check = []<typename T>(T, char const *rel) {
        T const inf = std::numeric_limits<T>::infinity();
        T const rel_tol(rel);
        auto const run = [&](char const *name, auto fn, T const &a, T const &b, char const *ref) {
            INFO(name << ", rel_tol " << rel);
            auto const r = qd::integrate<tanh_sinh>(fn, a, b, { .rel_tol = rel_tol });
            REQUIRE(r.has_value());
            double const error = error_of(r->value, ref);
            CHECK(error <= static_cast<double>(wide(r->abs_error))); // the honesty check
            CHECK(error <= static_cast<double>(wide(rel_tol)) * std::abs(static_cast<double>(bin100(ref))));
        };
        using std::exp, std::sqrt;
        run("e^-x", [](T const &x) { return T(exp(T(-x))); }, T(0), inf, "1");
        run("x^-0.5 e^-x", [](T const &x) { return T(exp(T(-x)) / sqrt(x)); }, T(0), inf, sqrt_pi);
        run("x^-1.5 [1, inf)", [](T const &x) { return T(T(1) / T(x * sqrt(x))); }, T(1), inf, "2");
        run("e^x (-inf, 0]", [](T const &x) { return T(exp(x)); }, T(-inf), T(0), "1");
        run("x^-2 (-inf, -1]", [](T const &x) { return T(T(1) / T(x * x)); }, T(-inf), T(-1), "1");
        run("e^-x^2", [](T const &x) { return T(exp(T(-T(x * x)))); }, T(-inf), inf, sqrt_pi);
        run("1/(1 + x^2)", [](T const &x) { return T(T(1) / T(T(1) + T(x * x))); }, T(-inf), inf, pi);
    };
    for (char const *rel : { "1e-25", "1e-45" }) {
        check(bin50{}, rel);
        check(bin50_et{}, rel);
        check(dec50{}, rel);
    }
}

TEST_CASE("quadrature: infinite: evaluation counts in double")
{
    // The kept prototype's results with the C++ table (es.py, cpp=True), per level, with the status and whether a
    // failure carries a best: absent only where the estimate is unbounded or nothing can be sampled.
    struct Case
    {
        char const *name;
        Run run;
        status expected;
        std::vector<int> counts;
        bool best;
    };
    Case const cases[] = {
        // [a, inf): pruning on both sides, the minimum level.
        { "1/(1 + x)^2 [0, inf)", run_observed(+inverse_square, 0.0, dinf), status::ok, { 13, 21, 37, 69, 133 },
          false },
        { "x^-0.5 / (1 + x) [0, inf)", run_observed([](double x) { return 1.0 / (root_of(x) * (1.0 + x)); }, 0.0, dinf),
          status::ok, { 13, 23, 43, 81, 155 }, false },
        // A rounding cut at a = 1 on the side toward it.
        { "x^-1.5 [1, inf)", run_observed([](double x) { return 1.0 / (x * root_of(x)); }, 1.0, dinf), status::ok,
          { 10, 19, 37, 72, 141 }, false },
        { "1/(1 + (x - 1e10))^2 [1e10, inf)",
          run_observed([](double x) { return 1.0 / ((1.0 + (x - 1e10)) * (1.0 + (x - 1e10))); }, 1e10, dinf),
          status::ok, { 9, 16, 30, 58, 113 }, false },
        // The table's end toward infinity and the underflow rule at level 7; at level 0 and k at 1e-250.
        { "x^-(1 + 1/64) [1, inf)", run_observed([](double x) { return 1.0 / (x * root64(x)); }, 1.0, dinf),
          status::not_converged, { 10, 21, 43, 86, 172, 342, 682, 1361, 2717, 5428 }, true },
        { "1e-250 x^-(1 + 1/64) [1, inf) 1e-3",
          run_observed([](double x) { return 1e-250 / (x * root64(x)); }, 1.0, dinf, { .rel_tol = 1e-3 }),
          status::not_converged, { 10, 20, 39, 76, 149, 293, 580, 1152, 2294, 4577 }, true },
        // The tail rule's second clause at level 4.
        { "x^-2 sqrt((x - X1)(x - X2)) / x [1, inf)", run_observed(level_k_tail, 1.0, dinf), status::ok,
          { 10, 18, 33, 63, 121 }, false },
        // A singular end other than 0: the spacing near 1.
        { "(x - 1)^-0.5 / x^2 [1, inf)", run_observed([](double x) { return 1.0 / (root_or_nan(x - 1.0) * x * x); },
                                                      1.0, dinf),
          status::not_converged, { 10, 18, 34, 65, 125, 244, 482, 956, 1903, 3796 }, true },
        // (-inf, b]: the side toward -inf first.
        { "1/(1 + x^2)^2 (-inf, 0]",
          run_observed([](double x) { return 1.0 / ((1.0 + x * x) * (1.0 + x * x)); }, -dinf, 0.0), status::ok,
          { 13, 20, 34, 62, 117 }, false },
        { "(-x)^-1.5 (-inf, -1]", run_observed([](double x) { return 1.0 / (-x * root_of(-x)); }, -dinf, -1.0),
          status::ok, { 10, 19, 37, 72, 141 }, false },
        { "1/(1 + (x - 3)^2) (-inf, -2]",
          run_observed([](double x) { return 1.0 / (1.0 + (x - 3.0) * (x - 3.0)); }, -dinf, -2.0), status::ok,
          { 10, 19, 36, 68, 132 }, false },
        // The line.
        { "1/(1 + x^2)", run_observed(+runge, -dinf, dinf), status::ok, { 13, 21, 37, 69, 133 }, false },
        { "(1 + x^2)^-1.5", run_observed([](double x) { return 1.0 / ((1.0 + x * x) * root_of(1.0 + x * x)); }, -dinf,
                                         dinf),
          status::ok, { 13, 21, 35, 63, 117 }, false },
        { "1/(1 + x^2)^3 rel 0",
          run_observed([](double x) { return 1.0 / ((1.0 + x * x) * (1.0 + x * x) * (1.0 + x * x)); }, -dinf, dinf,
                       { .rel_tol = 0.0 }),
          status::ok, { 13, 19, 29, 49, 89 }, false },
        // The cap, with no best: the estimate unbounded.
        { "1/x [0, inf)", run_observed([](double x) { return 1.0 / x; }, 0.0, dinf), status::not_converged,
          { 13, 27, 55, 109, 219, 438, 876, 1752, 3505, 7011 }, false },
        { "1/x [1, inf)", run_observed([](double x) { return 1.0 / x; }, 1.0, dinf), status::not_converged,
          { 10, 21, 43, 86, 172, 342, 682, 1361, 2718, 5431 }, false },
        { "1 [0, inf)", run_observed([](double) { return 1.0; }, 0.0, dinf), status::not_converged,
          { 13, 21, 36, 64, 120, 230, 450, 889, 1766, 3519 }, false },
        // An interior step: the cap with an honest best (the constant-evaluation case).
        { "the step at 2 [0, inf)", run_observed(+step_at_2, 0.0, dinf), status::not_converged,
          { 13, 18, 27, 45, 81, 152, 292, 570, 1125, 2233 }, true },
        // Level 0: a value not finite toward infinity, pending, then non_finite (its last term above the threshold).
        { "sqrt(10 - x) / (1 + x)^2 [0, inf)",
          run_observed([](double x) { return root_or_nan(10.0 - x) / ((1.0 + x) * (1.0 + x)); }, 0.0, dinf),
          status::non_finite, {}, false },
        // No point strictly inside: a + 1 rounds onto a.
        { "on [-1e17, inf)", run_observed(+runge, -1e17, dinf), status::not_converged, { 0 }, false },
    };
    for (Case const &c : cases) {
        INFO(c.name);
        auto const &r = c.run.result;
        CHECK(status_of(r) == c.expected);
        CHECK(c.run.counts == c.counts);
        if (!r) CHECK(r.error().best.has_value() == c.best);
    }
    // The failure at level 0 evaluated 9 points (the centre, 2 toward infinity, 6 toward 0), and reported none.
    int n = 0;
    auto const shifted = qd::integrate<tanh_sinh>(
        Counted{ [](double x) { return root_or_nan(10.0 - x) / ((1.0 + x) * (1.0 + x)); }, &n }, 0.0, dinf);
    CHECK(status_of(shifted) == status::non_finite);
    CHECK(n == 9);
    // The steps as reported: level 0 without a current, ready from level 4, exhausted only at the cap.
    Run const &cap = cases[14].run;
    REQUIRE(cap.steps.size() == 10);
    for (std::size_t i = 0; i < cap.steps.size(); ++i) {
        INFO("step " << i);
        CHECK(cap.steps[i].ready == (i >= 4));
        CHECK(cap.steps[i].exhausted == (i == 9));
    }
    CHECK(!cases[0].run.steps.front().current.has_value());
    CHECK(cases[0].run.steps.back().current == *cases[0].run.result);
    // The values and bests against closed forms or mpmath (the honesty check): 1, pi, 2, 1; 64 and 64 1e-250 (the
    // integral of x^-(1 + 1/64) on [1, inf) is 64); the level-k tail's integral over [1, X1] and [X2, inf) (mpmath);
    // B(1/2, 3/2) = pi / 2; pi / 4, 2, pi / 2 - atan 5; pi, 2, 3 pi / 8; atan 2.
    auto const honest = [&](std::size_t i, bin100 const &ref) {
        INFO(cases[i].name);
        auto const &r = cases[i].run.result;
        qd::Integral<double> const got = r ? *r : *r.error().best;
        CHECK(error_of(got.value, ref) <= got.abs_error);
    };
    honest(0, bin100(1));
    honest(1, bin100(pi));
    honest(2, bin100(2));
    honest(3, bin100(1));
    honest(4, bin100(64));
    honest(5, bin100(1e-250) * 64);
    honest(6, bin100("422897019.877332888190726155703294944220443073249716393198446")); // mpmath
    honest(7, bin100(half_pi));
    honest(8, bin100(pi) / 4);
    honest(9, bin100(2));
    honest(10, bin100("0.197395559849880758370049765194790293447585103787852101517689"));
    honest(11, bin100(pi));
    honest(12, bin100(2));
    honest(13, bin100(pi) * 3 / 8);
    honest(17, bin100("1.10714871779409050301706546017853704007004764540143264667654"));
}

TEST_CASE("quadrature: infinite: how the sides end")
{
    using state = qd::detail::TanhSinhState<double>;
    using mapping = qd::detail::mapping;
    auto const e_minus_x = [](double x) { return std::exp(-x); };

    SUBCASE("toward infinity: pruning, the table's end, the tail rule, the underflow rule")
    {
        // e^-x on [0, inf): both sides pruned by level 4, where it succeeds.
        state const pruned = state_at(e_minus_x, 0.0, dinf, 4);
        CHECK(pruned.map == mapping::upper);
        CHECK(pruned.right.cause == cut_by::pruning);
        CHECK(pruned.left.cause == cut_by::pruning);
        auto const e = qd::integrate<tanh_sinh>(e_minus_x, 0.0, dinf);
        REQUIRE(e.has_value());
        CHECK(error_of(e->value, "1") <= e->abs_error);

        // x^-1.01 on [1, inf): the side toward infinity ends as at the table's end (there, or by the underflow rule,
        // which cuts it the same way where x^-1.01 falls below min(), beyond about 4e304); an honest best,
        // 1 / (p - 1) for p the double 1.01.
        state const end = state_at([](double x) { return std::pow(x, -1.01); }, 1.0, dinf, 9);
        CHECK(end.right.cause == cut_by::table);
        auto const slow = qd::integrate<tanh_sinh>([](double x) { return std::pow(x, -1.01); }, 1.0, dinf);
        REQUIRE(!slow.has_value());
        REQUIRE(slow.error().best.has_value());
        CHECK(error_of(slow.error().best->value, "99.9999999999999111821580299875556521999886538448512496402875")
              <= slow.error().best->abs_error);

        // The tail rule at level 0: x^50 e^-x is NaN (inf 0) from t = 3, below eps L1 at t = 2: a tail (then t = 2
        // pruned), ok (50!); sqrt(10 - x) e^-x is NaN from t = 2 with a large term at t = 1: non_finite. x^51 e^-x:
        // non_finite (PHASE4_API.md F: ok up to k = 50).
        auto const power = [](int k) { return [k](double x) { return std::pow(x, k) * std::exp(-x); }; };
        CHECK(state_at(power(50), 0.0, dinf, 0).right.cut == 2 * 512);
        auto const k50 = qd::integrate<tanh_sinh>(power(50), 0.0, dinf);
        REQUIRE(k50.has_value());
        CHECK(error_of(k50->value, "3.04140932017133780436126081660647688443776415689605120000000e+64")
              <= k50->abs_error);
        qd::Failure<double> const invalid{ qd::error::non_finite, std::nullopt };
        CHECK(qd::integrate<tanh_sinh>(power(51), 0.0, dinf) == std::unexpected(invalid));
        CHECK(qd::integrate<tanh_sinh>([](double x) { return root_or_nan(10.0 - x) * std::exp(-x); }, 0.0, dinf)
              == std::unexpected(invalid));

        // The tail rule at level k: x^-2 sqrt((x - X1)(x - X2)) / x on [1, inf) is NaN between X1 and X2, first seen at
        // level 4 beyond a kept term below the threshold: a tail (the integral over [1, X1] and [X2, inf), mpmath);
        // e^-x sqrt((x - Y1)(x - Y2)) at level 1 beyond a large one: non_finite (r4_widek.txt).
        CHECK(state_at(level_k_tail, 1.0, dinf, 3).right.cause == cut_by::pruning);
        CHECK(state_at(level_k_tail, 1.0, dinf, 4).right.cause == cut_by::value);
        constexpr double y1 = 20.489948203756917; // 3 1.03^65
        constexpr double y2 = 24.5879378445083;   // 1.2 Y1
        CHECK(qd::integrate<tanh_sinh>([](double x) { return std::exp(-x) * root_or_nan((x - y1) * (x - y2)); }, 0.0,
                                       dinf)
              == std::unexpected(invalid));

        // The underflow rule: at level 0 (t = 6, where 1e-250 x^-(1 + 1/64) falls below min() with a weight of about
        // 8e133) and at level k (x^-(1 + 1/64) at level 7), each as the table's end, the node not kept.
        state const level0 = state_at([](double x) { return 1e-250 / (x * root64(x)); }, 1.0, dinf, 0);
        CHECK(level0.right.cause == cut_by::table);
        CHECK(level0.right.cut == 6 * 512);
        state const before = state_at([](double x) { return 1.0 / (x * root64(x)); }, 1.0, dinf, 6);
        CHECK(before.right.cut == -1);
        state const levelk = state_at([](double x) { return 1.0 / (x * root64(x)); }, 1.0, dinf, 7);
        CHECK(levelk.right.cause == cut_by::table);
        CHECK(levelk.right.cut > 0);
        // Its honest failures (PHASE4_API.md Q4): without the rule 1e-250 x^-1.01 at 1e-3 and 1e-15 x^-1.02 ended ok
        // with errors of 19 % and 8.9e-7. References: the scale and 1 / (p - 1) for the doubles written.
        auto const tiny = qd::integrate<tanh_sinh>([](double x) { return 1e-250 * std::pow(x, -1.01); }, 1.0, dinf,
                                                   { .rel_tol = 1e-3 });
        REQUIRE(!tiny.has_value());
        REQUIRE(tiny.error().best.has_value());
        CHECK(error_of(tiny.error().best->value, "9.99999999999999165821117553759408541979405345279455584478349e-249")
              <= tiny.error().best->abs_error);
        auto const small = qd::integrate<tanh_sinh>([](double x) { return 1e-15 * std::pow(x, -1.02); }, 1.0, dinf);
        REQUIRE(!small.has_value());
        CHECK(small.error().reason == qd::error::not_converged);
        REQUIRE(small.error().best.has_value());
        CHECK(error_of(small.error().best->value, "4.99999999999999594763490088268283372023903760280663522849167e-14")
              <= small.error().best->abs_error);
    }
    SUBCASE("toward the finite end")
    {
        // A table's end: x^-0.99 e^-x on [0, inf) at 1e-4 (fn finite below min()); an honest best, Gamma(1 - p).
        auto const strong = [](double x) { return std::pow(x, -0.99) * std::exp(-x); };
        CHECK(state_at(strong, 0.0, dinf, 9).left.cause == cut_by::table);
        auto const r = qd::integrate<tanh_sinh>(strong, 0.0, dinf, { .rel_tol = 1e-4 });
        REQUIRE(!r.has_value());
        CHECK(r.error().reason == qd::error::not_converged);
        REQUIRE(r.error().best.has_value());
        CHECK(error_of(r.error().best->value, "99.4325851191505149043169814984083123892644333785657256533160")
              <= r.error().best->abs_error);

        // A singular end other than 0: (x - 1)^-0.5 e^-x on [1, inf), sqrt(pi) / e, limited by the spacing of double
        // near 1 (a rounding cut): ok at 1e-6, an honest not_converged at the default; with u = x - 1 written in fn,
        // e^-(1 + u) / sqrt u on [0, inf), ok within the default.
        constexpr char const *ref = "0.652049332173292183059158613247067249185040697619428215807973";
        auto const shifted = [](double x) { return std::exp(-x) / root_or_nan(x - 1.0); };
        CHECK(state_at(shifted, 1.0, dinf, 4).left.cause == cut_by::rounding);
        auto const loose = qd::integrate<tanh_sinh>(shifted, 1.0, dinf, { .rel_tol = 1e-6 });
        REQUIRE(loose.has_value());
        CHECK(error_of(loose->value, ref) <= loose->abs_error);
        auto const tight = qd::integrate<tanh_sinh>(shifted, 1.0, dinf);
        REQUIRE(!tight.has_value());
        CHECK(tight.error().reason == qd::error::not_converged);
        REQUIRE(tight.error().best.has_value());
        CHECK(error_of(tight.error().best->value, ref) <= tight.error().best->abs_error);
        auto const substituted =
            qd::integrate<tanh_sinh>([](double u) { return std::exp(-(1.0 + u)) / root_of(u); }, 0.0, dinf);
        REQUIRE(substituted.has_value());
        CHECK(error_of(substituted->value, ref) <= substituted->abs_error);
        CHECK(error_of(substituted->value, ref) <= 1.5e-8 * 0.652);
    }
    SUBCASE("no integral: not_converged without a best")
    {
        qd::Failure<double> const none{ qd::error::not_converged, std::nullopt };
        CHECK(qd::integrate<tanh_sinh>([](double x) { return 1.0 / x; }, 0.0, dinf) == std::unexpected(none));
        CHECK(qd::integrate<tanh_sinh>([](double x) { return 1.0 / x; }, 1.0, dinf) == std::unexpected(none));
        CHECK(qd::integrate<tanh_sinh>([](double) { return 1.0; }, 0.0, dinf) == std::unexpected(none));
    }
}

TEST_CASE("quadrature: infinite: the limits' remedies")
{
    // A peak far from 0 on the line (e^-(x - 1000)^2, missed): split at the peak, each half ok, sqrt(pi) / 2; or
    // substitute x = 1000 + u, e^-u^2 on the line, ok.
    auto const peak = [](double x) { return std::exp(-(x - 1e3) * (x - 1e3)); };
    auto const lower = qd::integrate<tanh_sinh>(peak, -dinf, 1e3);
    auto const upper = qd::integrate<tanh_sinh>(peak, 1e3, dinf);
    REQUIRE(lower.has_value());
    REQUIRE(upper.has_value());
    for (auto const &half : { *lower, *upper }) {
        CHECK(error_of(half.value, half_sqrt_pi) <= half.abs_error);
        CHECK(half.abs_error <= 1.5e-8 * 0.887);
    }
    auto const centred = qd::integrate<tanh_sinh>([](double u) { return std::exp(-u * u); }, -dinf, dinf);
    REQUIRE(centred.has_value());
    CHECK(error_of(centred->value, sqrt_pi) <= centred->abs_error);
    CHECK(centred->abs_error <= 1.5e-8 * 1.773);

    // A half line far from 0: e^-(x - a) on [a, inf), a = 1e15, is cut by rounding near a, its estimate far above
    // the request (PHASE4_API.md F: 0.094); with x = a + u written in fn, e^-u on [0, inf), ok within it.
    auto const far = qd::integrate<tanh_sinh>([](double x) { return std::exp(-(x - 1e15)); }, 1e15, dinf);
    REQUIRE(far.has_value());
    CHECK(error_of(far->value, "1") <= far->abs_error);
    CHECK(far->abs_error > 1e-3);
    auto const near = qd::integrate<tanh_sinh>([](double u) { return std::exp(-u); }, 0.0, dinf);
    REQUIRE(near.has_value());
    CHECK(error_of(near->value, "1") <= near->abs_error);
    CHECK(near->abs_error <= 1.5e-8);
}

// Constant evaluation: 1/(1 + x^2) on the line (pi, 133 evaluations) and 1/(1 + x)^2 on [0, inf) (1, 133), with the
// honesty check against the doubles nearest pi and 1 (pi's half ulp far inside the estimate) and within the request;
// the step at 2 to the cap (2,233 evaluations; 408,612 Clang steps, 39 % of its default), not_converged with an honest
// best against the double nearest atan 2; a class type without a table, its levels 0 to 4 generated in the same
// expression (312,971 steps, 29.8 %); 1/x^2 on [1, inf) and on (-inf, -1], both 1, which put a rounding cut at the
// finite end (its sliver and power-law excess) and the reflected map with the side toward -inf first in a constant
// expression. The run-time twins check that run time gives the same results, bit for bit.
namespace {
    constexpr auto line_result = qd::integrate<tanh_sinh>(runge, -dinf, dinf);
    constexpr auto half_line_result = qd::integrate<tanh_sinh>(inverse_square, 0.0, dinf);
    constexpr auto cap_result = qd::integrate<tanh_sinh>(step_at_2, 0.0, dinf);
    constexpr auto wrapped_square = [](Wrapped const &x) { return Wrapped(1) / ((Wrapped(1) + x) * (Wrapped(1) + x)); };
    constexpr auto wrapped_result = qd::integrate<tanh_sinh>(wrapped_square, Wrapped(0), Wrapped(dinf));
    constexpr auto reciprocal_square = [](double x) { return 1.0 / (x * x); };
    constexpr auto upper_result = qd::integrate<tanh_sinh>(reciprocal_square, 1.0, dinf);
    constexpr auto lower_result = qd::integrate<tanh_sinh>(reciprocal_square, -dinf, -1.0);
} // namespace
static_assert(line_result && nxx::math::abs(line_result->value - 3.141592653589793) <= line_result->abs_error
              && line_result->abs_error <= 1.5e-8 * 3.1416);
static_assert(half_line_result && nxx::math::abs(half_line_result->value - 1.0) <= half_line_result->abs_error
              && half_line_result->abs_error <= 1.5e-8);
static_assert(!cap_result && cap_result.error().reason == qd::error::not_converged && cap_result.error().best
              && nxx::math::abs(cap_result.error().best->value - 1.1071487177940904)
                     <= cap_result.error().best->abs_error);
static_assert(wrapped_result && nxx::math::abs(wrapped_result->value.v - 1.0) <= wrapped_result->abs_error.v);
static_assert(upper_result && nxx::math::abs(upper_result->value - 1.0) <= upper_result->abs_error
              && upper_result->abs_error <= 1.5e-8);
static_assert(lower_result && nxx::math::abs(lower_result->value - 1.0) <= lower_result->abs_error
              && lower_result->abs_error <= 1.5e-8);
static_assert(qd::integrate<tanh_sinh>(runge, dinf, -dinf) == qd::integrate<tanh_sinh>(runge, -dinf, dinf).transform(
                                                                   [](qd::Integral<double> i) {
                                                                       i.value = -i.value;
                                                                       return i;
                                                                   }));
static_assert(qd::integrate<tanh_sinh>(runge, dinf, dinf) == qd::Integral<double>{ 0.0, 0.0 });
static_assert(qd::integrate<qd::Romberg>(runge, 0.0, dinf)
              == std::unexpected(qd::Failure<double>{ qd::error::non_finite, std::nullopt }));

TEST_CASE("quadrature: infinite: constant evaluation")
{
    CHECK(qd::integrate<tanh_sinh>(runge, -dinf, dinf) == line_result);
    CHECK(qd::integrate<tanh_sinh>(inverse_square, 0.0, dinf) == half_line_result);
    CHECK(qd::integrate<tanh_sinh>(step_at_2, 0.0, dinf) == cap_result);
    auto const wrapped = qd::integrate<tanh_sinh>(wrapped_square, Wrapped(0), Wrapped(dinf));
    REQUIRE(wrapped.has_value());
    CHECK(wrapped->value == wrapped_result->value);
    CHECK(wrapped->abs_error == wrapped_result->abs_error);
    CHECK(qd::integrate<tanh_sinh>(reciprocal_square, 1.0, dinf) == upper_result);
    CHECK(qd::integrate<tanh_sinh>(reciprocal_square, -dinf, -1.0) == lower_result);
    // Against the closed forms, from bin100.
    CHECK(error_of(line_result->value, pi) <= line_result->abs_error);
    CHECK(error_of(half_line_result->value, "1") <= half_line_result->abs_error);
    CHECK(error_of(cap_result.error().best->value, "1.10714871779409050301706546017853704007004764540143264667654")
          <= cap_result.error().best->abs_error);
}

// The independent review's cases (Codex, PHASE4_API.md): a type whose max() is small against 1 / eps, here 2^40 with
// 53 bits, so that a + 1 overflows on [max(), inf) and 1 / eps does not fit. As a Boost type at run time, and as a
// double wrapper in constant evaluation, whose overflow beyond 2^40 fails the constant expression.
namespace {
    using capped_boost = mp::number<mp::cpp_bin_float<53, mp::digit_base_2, void, int, -1022, 39>>;

    // Reached only on an overflow, which the library must never cause: not constexpr, so constant evaluation fails.
    [[noreturn]] void capped_overflow()
    { std::abort(); }

    struct Capped
    {
        double v;
        static constexpr double top = 1099511627776.0; // 2^40
        constexpr Capped(int i) : v(static_cast<double>(i)) {}
        constexpr Capped(double x) : v(x) {}
        static constexpr Capped checked(double x)
        {
            if (x > top || x < -top) capped_overflow();
            return Capped(x);
        }
        friend constexpr Capped operator+(Capped x, Capped y) { return checked(x.v + y.v); }
        friend constexpr Capped operator-(Capped x, Capped y) { return checked(x.v - y.v); }
        friend constexpr Capped operator*(Capped x, Capped y) { return checked(x.v * y.v); }
        friend constexpr Capped operator/(Capped x, Capped y) { return checked(x.v / y.v); }
        friend constexpr Capped operator-(Capped x) { return Capped(-x.v); }
        friend constexpr bool operator==(Capped, Capped) = default;
        friend constexpr auto operator<=>(Capped x, Capped y) { return x.v <=> y.v; }
    };

} // namespace

template<>
struct std::numeric_limits<Capped> : std::numeric_limits<double>
{
    static constexpr Capped max() noexcept { return Capped(Capped::top); }
    static constexpr Capped infinity() noexcept { return Capped(std::numeric_limits<double>::infinity()); }
};

namespace {
    constexpr Capped capped_inf = std::numeric_limits<Capped>::infinity();
    constexpr auto capped_zero = [](Capped const &) { return Capped(0); };
    // 1/(1 + x)^2 as (1/(1 + x))/(1 + x): (1 + x)^2 itself would overflow the type near its max().
    constexpr auto capped_square = [](Capped const &x) { return Capped(Capped(1) / (Capped(1) + x)) / (Capped(1) + x); };
    // 1/eps did not fit: fn = 0 on [0, inf) failed to compile; now 0 with its floor.
    constexpr auto capped_zero_result = qd::integrate<qd::TanhSinh>(capped_zero, Capped(0), capped_inf);
    constexpr auto capped_square_result = qd::integrate<qd::TanhSinh>(capped_square, Capped(0), capped_inf);
    // a + 1 overflowed on [max(), inf): now no point inside, not_converged without a best, fn not called.
    constexpr auto capped_top_result = qd::integrate<qd::TanhSinh>(capped_square, std::numeric_limits<Capped>::max(),
                                                                   capped_inf);
} // namespace
static_assert(nxx::real<Capped> && nxx::real<capped_boost>);
static_assert(capped_zero_result && capped_zero_result->value.v == 0.0);
static_assert(capped_square_result
              && nxx::math::abs(capped_square_result->value.v - 1.0) <= capped_square_result->abs_error.v);
static_assert(!capped_top_result && capped_top_result.error().reason == qd::error::not_converged
              && !capped_top_result.error().best);

TEST_CASE("quadrature: infinite: a type whose max() is small against 1 / eps")
{
    // The run-time twins of the constant expressions above.
    CHECK(qd::integrate<qd::TanhSinh>(capped_zero, Capped(0), capped_inf) == capped_zero_result);
    CHECK(qd::integrate<qd::TanhSinh>(capped_square, Capped(0), capped_inf) == capped_square_result);
    CHECK(qd::integrate<qd::TanhSinh>(capped_square, std::numeric_limits<Capped>::max(), capped_inf)
          == capped_top_result);
    // The Boost type: [max(), inf) and (-inf, -max()] sample nothing (fn was called at inf before), and an
    // ordinary integral, 1/(1 + x)^2 on [0, inf) = 1, is honest.
    capped_boost const top = (std::numeric_limits<capped_boost>::max)();
    capped_boost const inf = std::numeric_limits<capped_boost>::infinity();
    int calls = 0;
    auto const counted = [&calls](capped_boost const &x) {
        ++calls;
        return capped_boost(capped_boost(1 / (1 + x)) / (1 + x));
    };
    qd::Failure<capped_boost> const none{ qd::error::not_converged, std::nullopt };
    CHECK(qd::integrate<qd::TanhSinh>(counted, top, inf) == std::unexpected(none));
    CHECK(qd::integrate<qd::TanhSinh>(counted, -inf, capped_boost(-top)) == std::unexpected(none));
    CHECK(calls == 0);
    auto const one = qd::integrate<qd::TanhSinh>(counted, capped_boost(0), inf);
    REQUIRE(one.has_value());
    CHECK(error_of(one->value, "1") <= one->abs_error);
}
