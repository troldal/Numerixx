// Refined types (DESIGN §6.2): values that cannot hold an illegal state. A literal is checked at compile time by a
// consteval constructor; a run-time value goes through make(), which returns std::expected.
#pragma once

#include <numerixx/config.hpp>
#include <numerixx/core/error.hpp>
#include <numerixx/core/math.hpp>
#include <numerixx/core/scalar.hpp>

#include <compare>
#include <concepts>
#include <cstdint>
#include <expected>
#include <type_traits>

NXX_BEGIN_HEADER

namespace nxx
{
    // Validated inputs: refined values, max_iterations, bracket and the stop criteria. std::is_constructible_v is true
    // for tolerance<double> from double on every compiler, so generic code must test this instead (DESIGN §6.2).
    template<class T>
    inline constexpr bool is_refined_v = requires { typename std::remove_cvref_t<T>::nxx_refined_tag; };

    namespace detail
    {
        // Not constexpr: reaching it inside a consteval constructor is a compile error that names it.
        inline void literal_violates_invariant(const char*) noexcept {}

        // Library-internal construction key: bypasses validation where the library has already checked.
        struct trust_me
        {
            explicit constexpr trust_me() = default;
        };

        template<class Tag, class T>
        class refined
        {
            T v_;

        public:
            using value_type      = T;
            using nxx_refined_tag = void;

            constexpr refined(trust_me, T v) noexcept : v_(v) {}

            consteval refined(T v) : v_(v)
            {
                if (!Tag::check(v)) Tag::reject();    // the reason is a literal in reject(), so diagnostics quote it
            }

            template<class B>
                requires std::same_as<B, bool>
            refined(B) NXX_DELETE("a bool is not a numeric refinement");

            static constexpr auto make(T v) noexcept -> std::expected<refined, errc>
            {
                if (!Tag::check(v)) return std::unexpected(errc::invalid_input);
                return refined { trust_me {}, v };
            }

            constexpr T value() const noexcept { return v_; }

            friend constexpr auto operator<=>(const refined&, const refined&) = default;
        };
    }    // namespace detail

    namespace tag
    {
        struct positive_tolerance
        {
            static consteval void reject() { detail::literal_violates_invariant("a tolerance must be finite and > 0"); }
            template<class T>
            static constexpr bool check(const T& v)
            { return v > T(0) && math::isfinite(v); }
        };
        struct abs_tolerance
        {
            static consteval void reject() { detail::literal_violates_invariant("an absolute tolerance must be finite and >= 0"); }
            template<class T>
            static constexpr bool check(const T& v)
            { return v >= T(0) && math::isfinite(v); }
        };
        struct rel_tolerance
        {
            static consteval void reject() { detail::literal_violates_invariant("a relative tolerance must be in [0, 1)"); }
            template<class T>
            static constexpr bool check(const T& v)
            { return v >= T(0) && v < T(1); }
        };
    }    // namespace tag

    template<real T>
    using tolerance = detail::refined<tag::positive_tolerance, T>;    // finite, > 0: a criterion's only threshold
    template<real T>
    using abs_tolerance = detail::refined<tag::abs_tolerance, T>;    // finite, >= 0: the absolute part of a mixed test
    template<real T>
    using rel_tolerance = detail::refined<tag::rel_tolerance, T>;    // 0 <= r < 1; roles are not interchangeable

    // Traits for the reasoned deletions of the role-typed criteria (DESIGN §6.2) and of the solvers (§7.2).
    namespace detail
    {
        // A bare number: what a refined literal is spelled with, and what a criterion never takes as its relative part.
        template<class A>
        inline constexpr bool is_bare_number_v = std::is_arithmetic_v<std::remove_cvref_t<A>> || real<std::remove_cvref_t<A>>;
        // The criterion type that a deletion guide names for bare numbers: an integer literal names width_tol<double>.
        template<class A>
        using bare_scalar_t = std::conditional_t<real<std::remove_cvref_t<A>>, std::remove_cvref_t<A>, double>;

        // R is the relative part of a T criterion: rel_tolerance<T> exactly, so that no bare number converts into it.
        template<class R, class T>
        inline constexpr bool is_rel_v = std::same_as<std::remove_cvref_t<R>, rel_tolerance<T>>;

        // A part of a mixed tolerance, abs_tolerance<U> or rel_tolerance<U>, which alone is not a criterion.
        template<class R>
        inline constexpr bool is_tolerance_part_v = false;
        template<real U>
        inline constexpr bool is_tolerance_part_v<refined<tag::abs_tolerance, U>> = true;
        template<real U>
        inline constexpr bool is_tolerance_part_v<refined<tag::rel_tolerance, U>> = true;

        // What a solver's constructor rejects with a reason (§7.2): a bare number, a validated tolerance<U>
        // (brent{*tol}) or a part. One list, so that every solver deletion and brent's guide reject the same inputs.
        template<class R>
        inline constexpr bool not_a_criterion_v = is_bare_number_v<R> || is_tolerance_part_v<R>;
        template<real U>
        inline constexpr bool not_a_criterion_v<refined<tag::positive_tolerance, U>> = true;
    }    // namespace detail

    // 1 .. 2^32 - 1 iterations.
    class max_iterations
    {
        std::uint32_t n_;

        constexpr max_iterations(detail::trust_me, std::uint32_t n) noexcept : n_(n) {}

    public:
        using nxx_refined_tag = void;

        template<std::integral I>
            requires(!std::same_as<I, bool>)
        consteval max_iterations(I n) : n_(static_cast<std::uint32_t>(n))
        {
            if (n < I(1) || static_cast<unsigned long long>(n) > 0xFFFF'FFFFull)
                detail::literal_violates_invariant("max_iterations must be in [1, 2^32)");
        }

        template<class B>
            requires std::same_as<B, bool>
        max_iterations(B) NXX_DELETE("a bool is not an iteration count");

        // Checked in the source type, so a negative or oversized long long is rejected, never wrapped.
        static constexpr auto make(long long n) noexcept -> std::expected<max_iterations, errc>
        {
            if (n < 1 || n > 0xFFFF'FFFFll) return std::unexpected(errc::invalid_input);
            return max_iterations { detail::trust_me {}, static_cast<std::uint32_t>(n) };
        }

        constexpr std::uint32_t value() const noexcept { return n_; }

        friend constexpr bool operator==(max_iterations, max_iterations) = default;
    };
    // 1 .. 2^32 - 1 evaluations (the argument of the max_evaluations criterion). Checked in the source type, like
    // max_iterations: a negative value is rejected, never wrapped to about 4e9.
    class evaluation_budget
    {
        std::uint32_t n_;

        constexpr evaluation_budget(detail::trust_me, std::uint32_t n) noexcept : n_(n) {}

    public:
        using nxx_refined_tag = void;

        template<std::integral I>
            requires(!std::same_as<I, bool>)
        consteval evaluation_budget(I n) : n_(static_cast<std::uint32_t>(n))
        {
            if (n < I(1) || static_cast<unsigned long long>(n) > 0xFFFF'FFFFull)
                detail::literal_violates_invariant("an evaluation budget must be in [1, 2^32)");
        }

        template<class B>
            requires std::same_as<B, bool>
        evaluation_budget(B) NXX_DELETE("a bool is not an evaluation count");

        static constexpr auto make(long long n) noexcept -> std::expected<evaluation_budget, errc>
        {
            if (n < 1 || n > 0xFFFF'FFFFll) return std::unexpected(errc::invalid_input);
            return evaluation_budget { detail::trust_me {}, static_cast<std::uint32_t>(n) };
        }

        constexpr std::uint32_t value() const noexcept { return n_; }

        friend constexpr bool operator==(evaluation_budget, evaluation_budget) = default;
    };
}    // namespace nxx

NXX_END_HEADER
