// E22 (Chapter 7, section 7.4.2; values E22 in data/values.tex): the Jacobian of a set of balances. Two residuals in
// the unknowns x = (T, n), T in K and n in mol: an energy balance r1 = n H(T) - Q on the thread's enthalpy, and an
// ideal-gas pressure r2 = n R T / V - P, the model being one function that returns both. Parameters as in
// data/README.md (the generator uses the same): T = 1000 K, n = 2 mol, V = 1 m^3, R = 8.31446261815324 J/(mol K),
// with Q = n H(T) and P = n R T / V, so that both residuals vanish at the point. Two ways:
// - jacobian-model (the residuals) and jacobian-balances: the central difference applied by hand to the whole vector, column j from the model at
//   x_j + h_j and x_j - h_j with h_j = default_step(x_j), the step of Central<2>: two model calls per column;
// - jacobian-by-calls: one derivative() call per entry, each on a lambda that evaluates the whole model: two model
//   calls per entry with Central<2>.
// The answer key is the closed form: dr1/dT = n cp(T), dr1/dn = H(T), dr2/dT = n R / V, dr2/dn = R T / V.

#include <numerixx/deriv.hpp>

#include "shomate.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <print>

using nxx::deriv::derivative;

namespace {

    double H(double T) { return book::shomate_H(book::nitrogen, T); }

    double relative(double value, double exact) { return std::abs(value - exact) / std::abs(exact); }

} // namespace

// clang-format off
// <jacobian-model>
// The model: two residuals in the unknowns x = (T, n), an
// energy balance n H(T) - Q and a pressure n R T / V - P.
using Vec = std::array<double, 2>;
struct Balances
{
    double Q, P, V, R;
    Vec operator()(Vec const &x) const
    {
        return { x[1] * H(x[0]) - Q,
                 x[1] * R * x[0] / V - P };
    }
};
// </jacobian-model>

// J[i][j] = dr_i/dx_j by the central difference applied to
// the whole vector: two model calls per column. The book
// quotes the loop only.
template<typename Model>
std::array<Vec, 2> jacobian(Model const &r, Vec const &x)
{
    std::array<Vec, 2> J{};
// <jacobian-balances>
for (std::size_t j = 0; j < 2; ++j) {
    double const h = nxx::deriv::default_step(x[j]);
    Vec xp = x, xm = x;
    xp[j] += h;
    xm[j] -= h;
    Vec const rp = r(xp), rm = r(xm);
    for (std::size_t i = 0; i < 2; ++i) {
        J[i][j] = (rp[i] - rm[i]) / (2 * h);
    }
}
// </jacobian-balances>
    return J;
}
// clang-format on

int main()
{
    double const R = 8.31446261815324; // J/(mol K)
    double const T = 1000.0;           // K
    double const n = 2.0;              // mol
    double const V = 1.0;              // m^3
    Balances const balances{ n * H(T), n * R * T / V, V, R };
    Vec const x{ T, n };
    int calls = 0;
    auto const model = [&](Vec const &y) {
        ++calls;
        return balances(y);
    };
    double const exact[2][2] = { { n * book::shomate_cp(book::nitrogen, T), H(T) }, { n * R / V, R * T / V } };
    auto const worst = [&exact](std::array<Vec, 2> const &J) {
        double w = 0;
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j) w = std::max(w, relative(J[i][j], exact[i][j]));
        return w;
    };

    auto const J_hand = jacobian(model, x);
    int const hand_calls = calls;
    calls = 0;

    // clang-format off
// <jacobian-by-calls>
// With the library: one derivative() call per entry, each
// on the whole model with one unknown varied.
std::array<Vec, 2> J{};
for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) {
        auto const r_ij = [&](double v) {
            Vec y = x;
            y[j] = v;
            return model(y)[i];
        };
        auto const d = derivative(r_ij, x[j]);
        if (!d) {
            return 1;
        }
        J[i][j] = d->value;
    }
}
// </jacobian-by-calls>
    // clang-format on
    int const library_calls = calls;

    std::println("by hand:     [[{:.9e}, {:.9e}], [{:.9e}, {:.9e}]], {} model calls, worst relative error {:.1e}",
        J_hand[0][0], J_hand[0][1], J_hand[1][0], J_hand[1][1], hand_calls, worst(J_hand));
    std::println("by library:  [[{:.9e}, {:.9e}], [{:.9e}, {:.9e}]], {} model calls, worst relative error {:.1e}",
        J[0][0], J[0][1], J[1][0], J[1][1], library_calls, worst(J));
    std::println("closed form: [[{:.9e}, {:.9e}], [{:.9e}, {:.9e}]]", exact[0][0], exact[0][1], exact[1][0],
        exact[1][1]);
    // Every entry within 1e-8 of its closed form by both ways: Central<2> at the default step, about eps^(2/3)
    // relative times the factor |r| / (s |dr/dx|), where |r| is the size of the terms n H and n R T / V before they
    // cancel; the generator measures 2.4e-11 for both (E22). The counts: 2 per column by hand, 2 per entry by calls.
    bool const ok = worst(J_hand) <= 1e-8 && worst(J) <= 1e-8 && hand_calls == 4 && library_calls == 8;
    if (!ok) std::println("FAILED: a claim of the text does not hold");
    return ok ? 0 : 1;
}
