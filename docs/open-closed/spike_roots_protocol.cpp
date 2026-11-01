// Spike: an open method protocol for nxx::roots, built on the existing header without changing it.
// Checks that (1) library methods expressed as method objects give bit-for-bit the results of nxx::roots::solve,
// (2) a policy extracted from an if-constexpr chain (false position's scaling) reproduces three library methods,
// (3) a method the library does not have (Pegasus) plugs in from outside, in constant expressions,
// (4) both call forms work: solve<M>(fn, bracket) and solve(fn, bracket, criteria, m).
#include <numerixx/roots.hpp>

#include <cstdio>

namespace proposal {

    namespace rt = nxx::roots;
    using nxx::real;

    // --- What would become public: the protocol's argument types (today rt::detail::Start and Tolerance). ---------
    template<real T>
    using Tolerance = rt::detail::Tolerance<T>;
    template<real T, std::size_t K>
    using Start = rt::detail::Start<T, K>;

    // The values a method needs from fn at a point: M::values where M declares it, else 1 (replaces values_needed<M>).
    template<typename M>
    inline constexpr std::size_t values_of = [] {
        if constexpr (requires { M::values; })
            return std::size_t{ M::values };
        else
            return std::size_t{ 1 };
    }();

    // The protocol: a value type with a const member bracketed(fn, start, tolerance) (replaces detail::bracketing<M>,
    // a list of std::same_as).
    template<typename M, typename F, typename T>
    concept bracketing_method = std::copyable<M> && requires(M const m, F &fn, Start<T, values_of<M>> const &s,
                                                       Tolerance<T> const &t) {
        { m.bracketed(fn, s, t) } -> std::same_as<std::expected<rt::Root<T>, rt::error>>;
    };

    // --- The entry point: validation once, then the method. No list of methods anywhere. ------------------------
    // M is deduced from m, or named with m defaulted: solve<M>(fn, b) and solve(fn, b, c, M{...}) are one function.
    template<typename M, real T, typename F>
    requires bracketing_method<M, std::remove_reference_t<F>, T> && rt::detail::gives<F, T, values_of<M>>
    [[nodiscard]] constexpr std::expected<rt::Root<T>, rt::error> solve(F &&fn,
        rt::Bracket<T> const &bracket,
        rt::Criteria<std::type_identity_t<T>> const &criteria = {},
        M const &m = M{})
    {
        constexpr std::size_t K = values_of<M>;
        auto const tolerance = rt::detail::checked_tolerance(criteria);
        if (!tolerance) return std::unexpected(tolerance.error());
        auto const start = rt::detail::prepare<K>(fn, bracket, *tolerance);
        if (!start) return std::unexpected(start.error());
        if (start->root) return *start->root;
        return m.bracketed(fn, *start, *tolerance);
    }

    // --- Library methods as method objects: each forwards to the detail function the if-constexpr chain chose. ---
    struct Bisection
    {
        template<std::size_t K, real T, typename F>
        constexpr std::expected<rt::Root<T>, rt::error> bracketed(F &fn, Start<T, K> const &s,
            Tolerance<T> const &t) const
        { return rt::detail::bisection(fn, s, t); }
    };

    struct Brent
    {
        template<std::size_t K, real T, typename F>
        constexpr std::expected<rt::Root<T>, rt::error> bracketed(F &fn, Start<T, K> const &s,
            Tolerance<T> const &t) const
        { return rt::detail::brent<rt::Brent>(fn, s, t); }
    };

    struct Newton
    {
        static constexpr std::size_t values = 2;

        template<std::size_t K, real T, typename F>
        constexpr std::expected<rt::Root<T>, rt::error> bracketed(F &fn, Start<T, K> const &s,
            Tolerance<T> const &t) const
        { return rt::detail::brent<rt::Newton>(fn, s, t); }
    };

    // --- A strategy point made public: false position with a scaling policy. ----------------------------------------
    // Today AndersonBjorck, Illinois and RegulaFalsi share detail::false_position<M>, which picks the factor by
    // if constexpr on M. Here the factor is a policy object: factor(f_latest, f_new) when the new point has the sign
    // of the latest, so the retained end is retained again. The loop is detail::false_position, copied, with only that
    // line changed.
    template<typename Scaling>
    struct FalsePosition
    {
        Scaling scaling{};

        template<std::size_t K, real T, typename F>
        constexpr std::expected<rt::Root<T>, rt::error> bracketed(F &fn, Start<T, K> const &start,
            Tolerance<T> const &tolerance) const
        {
            namespace d = rt::detail;
            T retained = start.lower;
            T f_retained = start.at_lower[0];
            T scaled = f_retained;
            T latest = start.upper;
            T f_latest = start.at_upper[0];
            if (std::optional<rt::Root<T>> const r = d::stop_at(retained, f_retained, latest, f_latest, tolerance))
                return *r;
            T reference = d::distance_bound(retained, latest);
            int unhalved = 0;
            while (true) {
                bool const ordered = retained < latest;
                d::Enclosure<T> const s{ ordered ? retained : latest, ordered ? f_retained : f_latest,
                    ordered ? latest : retained, ordered ? f_latest : f_retained, retained, f_retained, retained,
                    f_retained };
                bool const bisect = unhalved >= 3;
                std::optional<T> proposed;
                if (!d::wide(s)) {
                    if (bisect) {
                        proposed = d::midpoint(s.a, s.b);
                    } else if (std::optional<T> const spread = d::difference(f_latest, scaled)) {
                        if (std::optional<T> const r = d::quotient(f_latest, *spread))
                            proposed = T(latest - T(*r * T(latest - retained)));
                    }
                }
                std::expected<d::Enclosed<T>, rt::error> const step = d::enclose(fn, s, proposed, false, tolerance);
                if (!step) return std::unexpected(step.error());
                if (step->root) return *step->root;
                if (!d::same_sign(step->fc, f_latest)) {
                    retained = latest;
                    f_retained = f_latest;
                    scaled = f_latest;
                } else {
                    scaled = T(scaled * scaling.factor(f_latest, step->fc)); // the one line that differs
                }
                latest = step->c;
                f_latest = step->fc;
                T const width = d::distance_bound(retained < latest ? retained : latest,
                    retained < latest ? latest : retained);
                if (bisect || width <= T(reference / T(2))) {
                    reference = width;
                    unhalved = 0;
                } else {
                    ++unhalved;
                }
            }
        }
    };

    struct AndersonBjorckScaling
    {
        template<real T>
        constexpr T factor(T const &f_latest, T const &f_new) const
        {
            return nxx::math::abs(f_new) < nxx::math::abs(f_latest) ? T(T(1) - T(f_new / f_latest)) : T(T(1) / T(2));
        }
    };

    struct IllinoisScaling
    {
        template<real T>
        constexpr T factor(T const &, T const &) const { return T(T(1) / T(2)); }
    };

    struct NoScaling
    {
        template<real T>
        constexpr T factor(T const &, T const &) const { return T(1); }
    };

    using AndersonBjorck = FalsePosition<AndersonBjorckScaling>;
    using Illinois = FalsePosition<IllinoisScaling>;
    using RegulaFalsi = FalsePosition<NoScaling>;

} // namespace proposal

// --- A user's code, outside every nxx namespace: a method the library does not have. -----------------------------
namespace user {

    // The Pegasus method (Dowell and Jarratt 1972): the retained value scaled by f_latest / (f_latest + f_new), written
    // as 1 / (1 + f_new / f_latest) (the ratio is positive: the two values have one sign), and 1/2 where |f_new| is
    // not below |f_latest|, so that the ratio cannot overflow.
    struct PegasusScaling
    {
        template<nxx::real T>
        constexpr T factor(T const &f_latest, T const &f_new) const
        {
            if (!(nxx::math::abs(f_new) < nxx::math::abs(f_latest))) return T(T(1) / T(2));
            return T(T(1) / T(T(1) + T(f_new / f_latest)));
        }
    };
    using Pegasus = proposal::FalsePosition<PegasusScaling>;

    // A method written from scratch against the protocol alone (a deliberately trivial one): the point a third of the
    // way up the bracket, keeping the sign change. Only the protocol's types are used.
    struct Trisection
    {
        template<std::size_t K, nxx::real T, typename F>
        constexpr std::expected<nxx::roots::Root<T>, nxx::roots::error> bracketed(F &fn,
            proposal::Start<T, K> const &s, proposal::Tolerance<T> const &tol) const
        {
            T a = s.lower;
            T fa = s.at_lower[0];
            T b = s.upper;
            while (true) {
                T const c = T(a + T(T(b - a) / T(3)));
                T const half = T(T(b - a) / T(2));
                if (!(a < c && c < b) || half <= tol.at(c)) return nxx::roots::Root<T>{ c, T(b - a) };
                T const fc = T(std::invoke(fn, c));
                if (!nxx::math::isfinite(fc)) return std::unexpected(nxx::roots::error::non_finite);
                if (fc == T(0)) return nxx::roots::Root<T>{ c, tol.at(c) };
                if ((fc < T(0)) == (fa < T(0))) {
                    a = c;
                    fa = fc;
                } else {
                    b = c;
                }
            }
        }
    };

} // namespace user

namespace rt = nxx::roots;
namespace pp = proposal;

constexpr double square_minus_2(double x) { return x * x - 2.0; }
constexpr double wallis(double x) { return (x * x - 2.0) * x - 5.0; }
constexpr double power20(double x)
{
    double y = x * x;     // x^2
    y = y * y;            // x^4
    double const x16 = (y * y) * (y * y);
    return x16 * y - 1.0; // x^20 - 1
}
constexpr std::pair<double, double> fdf(double x) { return { x * x - 2.0, 2.0 * x }; }

// (1) The method objects give the library's results bit for bit (Root's operator== compares value and abs_error).
static_assert(pp::solve<pp::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 })
              == rt::solve<rt::Brent>(square_minus_2, rt::Bracket{ 1.0, 2.0 }));
static_assert(pp::solve<pp::Bisection>(wallis, rt::Bracket{ 2.0, 3.0 })
              == rt::solve<rt::Bisection>(wallis, rt::Bracket{ 2.0, 3.0 }));
static_assert(pp::solve<pp::Newton>(fdf, rt::Bracket{ 1.0, 2.0 }) == rt::solve<rt::Newton>(fdf, rt::Bracket{ 1.0, 2.0 }));

// (2) The scaling policy reproduces the three library methods that today share an if-constexpr chain.
static_assert(pp::solve<pp::AndersonBjorck>(power20, rt::Bracket{ 0.0, 5.0 })
              == rt::solve<rt::AndersonBjorck>(power20, rt::Bracket{ 0.0, 5.0 }));
static_assert(pp::solve<pp::Illinois>(power20, rt::Bracket{ 0.0, 5.0 })
              == rt::solve<rt::Illinois>(power20, rt::Bracket{ 0.0, 5.0 }));
static_assert(pp::solve<pp::RegulaFalsi>(wallis, rt::Bracket{ 2.0, 3.0 })
              == rt::solve<rt::RegulaFalsi>(wallis, rt::Bracket{ 2.0, 3.0 }));

// (3) User methods, in constant expressions, with the library's validation (zero first, no_sign_change, ...).
static_assert(pp::solve<user::Pegasus>(square_minus_2, rt::Bracket{ 1.0, 2.0 }).has_value());
static_assert(pp::solve<user::Trisection>(square_minus_2, rt::Bracket{ 1.0, 2.0 }).has_value());
static_assert(pp::solve<user::Pegasus>(square_minus_2, rt::Bracket{ 2.0, 3.0 })
              == std::unexpected(rt::error::no_sign_change));

// (4) A method object passed by value, for a method with parameters.
static_assert(pp::solve(square_minus_2, rt::Bracket{ 1.0, 2.0 }, {}, user::Pegasus{})
              == pp::solve<user::Pegasus>(square_minus_2, rt::Bracket{ 1.0, 2.0 }));

// The protocol rejects what is not a method, as detail::bracketing does today.
template<typename M>
concept solvable = requires { pp::solve<M>(square_minus_2, rt::Bracket{ 1.0, 2.0 }); };
static_assert(!solvable<int> && !solvable<rt::Brent>); // rt::Brent is a tag without the protocol
static_assert(!solvable<pp::Newton>);                  // Newton needs fn returning (f, f'): rejected by gives<>

int main()
{
    int calls = 0;
    auto const counted = [&calls](double x) {
        ++calls;
        return power20(x);
    };
    auto const report = [&](char const *name, auto const &r) {
        std::printf("%-15s calls %3d  value %.17g  abs_error %.3g\n", name, calls, r->value, r->abs_error);
        calls = 0;
    };
    report("Brent", pp::solve<pp::Brent>(counted, rt::Bracket{ 0.0, 5.0 }));
    report("AndersonBjorck", pp::solve<pp::AndersonBjorck>(counted, rt::Bracket{ 0.0, 5.0 }));
    report("Illinois", pp::solve<pp::Illinois>(counted, rt::Bracket{ 0.0, 5.0 }));
    report("Pegasus", pp::solve<user::Pegasus>(counted, rt::Bracket{ 0.0, 5.0 }));
    report("Trisection", pp::solve<user::Trisection>(counted, rt::Bracket{ 0.0, 5.0 }));
}
