// E13 (Chapter 7, section 7.4; values DS18): sensitivities of an Arrhenius rate k = A exp(-Ea / (R T)) to the
// temperature and to the activation energy, each by a lambda that holds the other arguments fixed; the scaled
// sensitivity d ln k / d ln Ea that engineers report; and a loop that fills one Jacobian column per parameter and
// counts the evaluations. Parameters as in data/README.md (the generator uses the same): A = 1e10 1/s,
// Ea = 80 kJ/mol, T = 600 K, R = 8.31446261815324 J/(mol K). The answer keys are the closed forms
// dk/dT = k Ea / (R T^2), dk/dEa = -k / (R T), dk/dA = k / A and d ln k / d ln Ea = -Ea / (R T).

#include <numerixx/deriv.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <print>

using nxx::deriv::derivative;
using nxx::deriv::Ridders;

// clang-format off
// <rate>
// An Arrhenius rate constant in 1/s: A in 1/s, Ea in J/mol,
// T in K.
constexpr double R = 8.31446261815324; // J/(mol K)

double rate(double A, double Ea, double T)
{ return A * std::exp(-Ea / (R * T)); }
// </rate>
// clang-format on

namespace {

    double relative(double value, double exact) { return std::abs(value - exact) / std::abs(exact); }

    // Ridders' criterion of the tests: the actual error at most 6 times the estimate.
    template<typename E>
    bool covered(E const &r, double exact)
    { return r && r->abs_error && std::abs(r->value - exact) <= 6 * *r->abs_error; }

} // namespace

int main()
{
    bool ok = true;
    // clang-format off
// <sensitivity>
double const A = 1e10, Ea = 80e3, T = 600.0;
// k as a function of one argument, the others held fixed
auto const k_T = [&](double t) { return rate(A, Ea, t); };
auto const k_Ea = [&](double e) { return rate(A, e, T); };
// dk/dT in 1/(s K), dk/dEa in mol/(s J)
auto const dk_dT = derivative<Ridders>(k_T, T);
auto const dk_dEa = derivative<Ridders>(k_Ea, Ea);
if (dk_dT && dk_dEa) {
    std::println("dk/dT = {:.6e} +- {:.1e}", dk_dT->value,
        *dk_dT->abs_error);
}
// </sensitivity>
    // clang-format on
    double const k0 = rate(A, Ea, T);
    double const exact_T = k0 * Ea / (R * T * T);
    double const exact_Ea = -k0 / (R * T);
    if (!dk_dT || !dk_dEa) {
        std::println("FAILED: no sensitivity");
        return 1;
    }
    std::println("dk/dEa = {:.6e} +- {:.1e}", dk_dEa->value, *dk_dEa->abs_error);
    std::println("relative errors: {:.1e} (T), {:.1e} (Ea)", relative(*dk_dT, exact_T), relative(*dk_dEa, exact_Ea));
    // Ridders: the generator measures 1.1e-13 and 1.2e-13 (DS18); the text claims coverage and 1e-8.
    ok = ok && covered(dk_dT, exact_T) && covered(dk_dEa, exact_Ea);
    ok = ok && relative(*dk_dT, exact_T) <= 1e-8 && relative(*dk_dEa, exact_Ea) <= 1e-8;

    // clang-format off
// <scaled>
// The scaled sensitivity d ln k/d ln Ea = (Ea/k) dk/dEa:
// the relative change of k per relative change of Ea.
double const k = rate(A, Ea, T);
double const S = Ea / k * dk_dEa->value; // dimensionless
std::println("d ln k / d ln Ea = {:.4f}", S);
// </scaled>
    // clang-format on
    double const S_exact = -Ea / (R * T);
    std::println("exact {:.4f}, relative error {:.1e}", S_exact, relative(S, S_exact));
    ok = ok && relative(S, S_exact) <= 1e-8;

    // clang-format off
// <jacobian>
// dk/dp for the parameters p = (A, Ea): one Ridders call
// per column, the other parameter held fixed. evals counts
// the calls of the model: the cost of the Jacobian.
std::array<double, 2> const p{ A, Ea };
std::array<double, 2> J{};
std::array<int, 2> cost{};
int evals = 0;
for (std::size_t j = 0; j < p.size(); ++j) {
    auto const k_j = [&](double v) {
        auto q = p;
        q[j] = v;
        ++evals;
        return rate(q[0], q[1], T);
    };
    int const before = evals;
    auto const d = derivative<Ridders>(k_j, p[j]);
    if (!d) return 1;
    J[j] = d->value;
    cost[j] = evals - before;
}
std::println("J = [{:.6e}, {:.6e}], {} evaluations",
    J[0], J[1], evals);
// </jacobian>
    // clang-format on
    std::println("per column: {} (A), {} (Ea)", cost[0], cost[1]);
    // k is linear in A, so Ridders stops at the first exact agreement (the generator measures 4 evaluations and a
    // relative error of 0); the Ea column is the same call as above.
    ok = ok && relative(J[0], k0 / A) <= 1e-8 && relative(J[1], exact_Ea) <= 1e-8;
    // The cost is the sum over the columns, each within Ridders' cap of 64 evaluations for f'.
    ok = ok && evals == cost[0] + cost[1] && cost[0] >= 1 && cost[0] <= 64 && cost[1] >= 1 && cost[1] <= 64;

    if (!ok) std::println("FAILED: a claim of the text does not hold");
    return ok ? 0 : 1;
}
