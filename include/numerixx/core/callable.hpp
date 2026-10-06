// Callables and evaluation (DESIGN §6.4): one unwrapping rule for callbacks, evaluation with finiteness checks, the
// cost of an evaluation in calls of the user's function, the common cause of several callbacks, and the copyable box
// that keeps solvers and chains copy-assignable when they hold lambdas.
#pragma once

#include <numerixx/config.hpp>
#include <numerixx/core/error.hpp>
#include <numerixx/core/math.hpp>
#include <numerixx/core/scalar.hpp>

#include <concepts>
#include <cstdint>
#include <expected>
#include <functional>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>

NXX_BEGIN_HEADER

namespace nxx
{
    namespace detail
    {
        template<class T>
        inline constexpr bool is_reference_wrapper_v = false;
        template<class T>
        inline constexpr bool is_reference_wrapper_v<std::reference_wrapper<T>> = true;

        // The callable behind a std::reference_wrapper (one-shot solves hold the user's function by std::cref).
        template<class F>
        constexpr decltype(auto) unref(const F& fn) noexcept
        {
            if constexpr (is_reference_wrapper_v<F>)
                return fn.get();
            else
                return (fn);
        }

        // x -> T: cannot fail; x -> expected<T, E>: fails with E; x -> expected<T, fault<UE>>: already a Numerixx fault,
        // unwrapped so that errors never nest (derivative_of returns this form).
        template<class R>
        struct unwrap
        {
            using value               = R;
            using error               = none;
            static constexpr int kind = 0;
        };
        template<class V, class E>
        struct unwrap<std::expected<V, E>>
        {
            using value               = V;
            using error               = E;
            static constexpr int kind = 1;
        };
        template<class V, class UE>
        struct unwrap<std::expected<V, fault<UE>>>
        {
            using value               = V;
            using error               = UE;
            static constexpr int kind = 2;
        };

        template<class F, class X>
        using call_result_t = std::remove_cvref_t<std::invoke_result_t<const F&, const X&>>;
    }    // namespace detail

    // The user's error type of a callback F at X: none for x -> T, E for x -> expected<T, E>, UE for
    // x -> expected<T, fault<UE>>.
    template<class F, class X>
    using callback_error_t = typename detail::unwrap<detail::call_result_t<F, X>>::error;

    namespace detail
    {
        // Evaluation does not throw when the call does not and copying its value and its error does not.
        template<class F, class X>
        inline constexpr bool nothrow_evaluation_v =
            std::is_nothrow_invocable_v<const F&, const X&> && std::is_nothrow_copy_constructible_v<X> &&
            std::is_nothrow_copy_constructible_v<callback_error_t<F, X>>;
    }    // namespace detail

    // nxx::cost_of(fn): what one call of fn costs in calls of the user's function (D33). The default is 1; a callable
    // that evaluates the user's function several times (a finite-difference derivative) reports evaluation_cost().
    namespace detail
    {
        struct cost_of_fn
        {
            template<class F>
            constexpr std::uint32_t operator()(const F& fn) const noexcept
            {
                const auto& g = nxx::detail::unref(fn);
                if constexpr (requires {
                                  { g.evaluation_cost() } -> std::convertible_to<std::uint32_t>;
                              })
                    return static_cast<std::uint32_t>(g.evaluation_cost());
                else
                    return 1u;
            }
        };
    }    // namespace detail

    inline constexpr detail::cost_of_fn cost_of {};

    namespace detail
    {
        template<bool Sample, class X, class F>
        constexpr auto evaluate_impl(const F& fn, const X& x) noexcept(detail::nothrow_evaluation_v<F, X>)
            -> std::expected<X, fault<callback_error_t<F, X>>>
        {
            using U  = unwrap<call_result_t<F, X>>;
            using UE = typename U::error;
            using R  = std::expected<X, fault<UE>>;

            const auto rejected = [](const X& y) {
                if constexpr (!real<X>)
                    return false;    // vector values: the caller checks
                else if constexpr (Sample)
                    return math::isnan(y);    // a sample may be +-inf (a pole or a log at 0); only NaN is an error
                else
                    return !math::isfinite(y);
            };

            if constexpr (U::kind == 0) {
                const X y = std::invoke(fn, x);
                if (rejected(y)) return R { std::unexpect, fault<UE> { errc::non_finite_value, cost_of(fn), {} } };
                return R { y };
            }
            else {
                auto r = std::invoke(fn, x);
                if (!r) {
                    if constexpr (U::kind == 1)
                        return R { std::unexpect, fault<UE> { errc::callback_failed, cost_of(fn), cause_slot<UE>(r.error()) } };
                    else {
                        // Already a Numerixx fault: it reports its own cost and cause. Its input codes become
                        // non_finite_value (DESIGN §6.4, §12 item 22): an input code means "the caller's input was
                        // rejected before iterating", and the nested callable's input is not the caller's input. Other
                        // input codes (no_sign_change from a nested solve, ...) pass through unchanged.
                        auto e = r.error();
                        if (e.code == errc::invalid_input || e.code == errc::non_finite_input) e.code = errc::non_finite_value;
                        return R { std::unexpect, e };    // copies, as nothrow_evaluation_v assumes
                    }
                }
                const X y = *std::move(r);
                if (rejected(y)) return R { std::unexpect, fault<UE> { errc::non_finite_value, cost_of(fn), {} } };
                return R { y };
            }
        }
    }    // namespace detail

    // One evaluation: NaN or +-inf -> errc::non_finite_value; the callback's own error -> errc::callback_failed with the
    // error as cause; a Numerixx fault from the callback (derivative_of, a nested solve) -> that fault, with
    // invalid_input and non_finite_input turned into non_finite_value. A failed evaluation reports its cost in
    // fault::evaluations.
    template<class X, class F>
    constexpr auto evaluate(const F& fn, const X& x) noexcept(detail::nothrow_evaluation_v<F, X>)
        -> std::expected<X, fault<callback_error_t<F, X>>>
    { return detail::evaluate_impl<false>(fn, x); }

    // A sample for a bracketing method: +-inf is a signed value, only NaN is an error.
    template<class X, class F>
    constexpr auto evaluate_sample(const F& fn, const X& x) noexcept(detail::nothrow_evaluation_v<F, X>)
        -> std::expected<X, fault<callback_error_t<F, X>>>
    { return detail::evaluate_impl<true>(fn, x); }

    // The error type of a solver that calls several callbacks (f and f'): the same type, or the one that is not none.
    namespace detail
    {
        template<class A, class B>
        struct common_cause
        {
            static_assert(std::is_same_v<A, B>,
                          "nxx: the function and its derivative report different error types; "
                          "map one with .transform_error so they agree");
            using type = A;
        };
        template<class A>
        struct common_cause<A, none>
        {
            using type = A;
        };
        template<class B>
        struct common_cause<none, B>
        {
            using type = B;
        };
        template<>
        struct common_cause<none, none>
        {
            using type = none;
        };
    }    // namespace detail

    template<class A, class B>
    using common_cause_t = typename detail::common_cause<A, B>::type;

    // nxx::derivative_source(fn): the structural derivative of a callable that has one (fn.derivative(), e.g. a
    // polynomial or a spline), so Newton needs no explicit derivative for it (D13).
    namespace detail
    {
        struct derivative_source_fn
        {
            template<class F>
                requires requires(const F& fn) { nxx::detail::unref(fn).derivative(); }
            constexpr auto operator()(const F& fn) const
            { return nxx::detail::unref(fn).derivative(); }
        };
    }    // namespace detail

    inline constexpr detail::derivative_source_fn derivative_source {};

    template<class F>
    inline constexpr bool has_derivative_source_v = std::is_invocable_v<const detail::derivative_source_fn&, const F&>;

    namespace detail
    {
        // copyable_box<T> (DESIGN §3.2): holds a copy-constructible T and is copy-assignable even when T is not (a
        // lambda with captures), using the std::ranges movable-box technique: assignment is destroy + construct.
        // Kind 3 covers the common case of a capture whose copy may throw but whose move cannot (a std::vector, a
        // std::string): it copies into a temporary first, so a throwing copy leaves the box unchanged. It also keeps
        // std::optional out of that path: GCC 16.2 reports a false -Wmaybe-uninitialized for optional's reset + emplace.
        template<class T>
        inline constexpr int box_kind_v = std::is_copy_assignable_v<T>              ? 0     // T itself
                                          : std::is_empty_v<T>                      ? 1     // no state
                                          : std::is_nothrow_copy_constructible_v<T> ? 2     // in place
                                          : std::is_nothrow_move_constructible_v<T> ? 3     // copy, then move in place
                                                                                    : 4;    // optional

        template<class T, int Kind = box_kind_v<T>>
        class copyable_box;

        template<class T>
        class copyable_box<T, 0>
        {
            NXX_NO_UNIQUE_ADDRESS T v_;

        public:
            constexpr copyable_box()
                requires std::default_initializable<T>
                : v_()
            {}
            constexpr explicit copyable_box(T v) noexcept(std::is_nothrow_move_constructible_v<T>) : v_(std::move(v)) {}
            constexpr const T& operator*() const noexcept { return v_; }
            constexpr const T* operator->() const noexcept { return std::addressof(v_); }
        };

        template<class T>
        class copyable_box<T, 1>
        {
            NXX_NO_UNIQUE_ADDRESS T v_;

        public:
            constexpr copyable_box()
                requires std::default_initializable<T>
                : v_()
            {}
            constexpr explicit copyable_box(T v) noexcept(std::is_nothrow_move_constructible_v<T>) : v_(std::move(v)) {}
            constexpr copyable_box(const copyable_box&) = default;
            constexpr copyable_box& operator=(const copyable_box&) noexcept { return *this; }    // nothing to copy
            constexpr const T&      operator*() const noexcept { return v_; }
            constexpr const T*      operator->() const noexcept { return std::addressof(v_); }
        };

        template<class T>
        class copyable_box<T, 2>
        {
            T v_;

        public:
            constexpr copyable_box()
                requires std::default_initializable<T>
                : v_()
            {}
            constexpr explicit copyable_box(T v) noexcept(std::is_nothrow_move_constructible_v<T>) : v_(std::move(v)) {}
            constexpr copyable_box(const copyable_box&) = default;
            constexpr copyable_box& operator=(const copyable_box& other) noexcept
            {
                if (this != std::addressof(other)) {
                    std::destroy_at(std::addressof(v_));
                    std::construct_at(std::addressof(v_), other.v_);
                }
                return *this;
            }
            constexpr const T& operator*() const noexcept { return v_; }
            constexpr const T* operator->() const noexcept { return std::addressof(v_); }
        };

        template<class T>
        class copyable_box<T, 3>
        {
            T v_;

        public:
            constexpr copyable_box()
                requires std::default_initializable<T>
                : v_()
            {}
            constexpr explicit copyable_box(T v) noexcept : v_(std::move(v)) {}
            constexpr copyable_box(const copyable_box&) = default;
            constexpr copyable_box& operator=(const copyable_box& other)
            {
                if (this != std::addressof(other)) {
                    T copy(other.v_);    // may throw: the box is still unchanged
                    std::destroy_at(std::addressof(v_));
                    std::construct_at(std::addressof(v_), std::move(copy));    // cannot throw
                }
                return *this;
            }
            constexpr const T& operator*() const noexcept { return v_; }
            constexpr const T* operator->() const noexcept { return std::addressof(v_); }
        };

        template<class T>
        class copyable_box<T, 4>
        {
            std::optional<T> v_;    // empty only if a copy threw during assignment, as with std::ranges' movable-box

        public:
            constexpr copyable_box()
                requires std::default_initializable<T>
                : v_(std::in_place)
            {}
            constexpr explicit copyable_box(T v) : v_(std::in_place, std::move(v)) {}
            constexpr copyable_box(const copyable_box&) = default;
            constexpr copyable_box& operator=(const copyable_box& other)
            {
                if (this != std::addressof(other)) {
                    if (other.v_)
                        v_.emplace(*other.v_);    // destroys the old value first; empty if the copy throws
                    else
                        v_.reset();
                }
                return *this;
            }
            constexpr const T& operator*() const noexcept { return *v_; }
            constexpr const T* operator->() const noexcept { return std::addressof(*v_); }
        };
    }    // namespace detail
}    // namespace nxx

NXX_END_HEADER
