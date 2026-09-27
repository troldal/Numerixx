// Prototype of run-time solver chains: any_solver (a copyable, type-erased curried solver) and first_of
// over a run-time range, with the same semantics as the static first_of in core.hpp.
// Standard library only (no FXT). Heap: std::function may allocate when a solver is wrapped or copied
// (target larger than its small buffer); a call never allocates by itself.
#pragma once
#include "core.hpp"
#include <functional>
#include <ranges>
#include <vector>

namespace nxx {

// A curried solver for callables of type F: a copyable value s with s(f) -> R for const F& f.
template<class S, class F, class R>
concept curried_solver_for = std::copy_constructible<S> && std::invocable<const S&, const F&> &&
                             std::same_as<std::invoke_result_t<const S&, const F&>, R>;

// any_solver<F, Est, UE>: "const F& -> result<Est, UE>" for ONE fixed callable type F
// (e.g. std::function<double(double)>). Built on std::function (std::move_only_function is missing on
// libc++ 22 / Emscripten, and copyability is the point here).
// Invariant: never empty. There is no default constructor and there are no move operations (a move is a
// copy), because a moved-from std::function is unspecified and may be empty.
template<class F, class Est, class UE = none>
class any_solver {
public:
    using function_type = F;
    using estimate_type = Est;
    using cause_type = UE;
    using result_type = result<Est, UE>;

    template<class S>
        requires(!std::same_as<std::remove_cvref_t<S>, any_solver> &&
                 curried_solver_for<std::remove_cvref_t<S>, F, result<Est, UE>>)
    any_solver(S&& s) : impl_(std::forward<S>(s)) {}   // implicit: every matching solver value IS an any_solver

    template<class S>
        requires(!std::same_as<std::remove_cvref_t<S>, any_solver> &&
                 !curried_solver_for<std::remove_cvref_t<S>, F, result<Est, UE>>)
    any_solver(S&&) NXX_DELETE("nxx::any_solver<F, Est, UE>: needs a copyable curried solver `const F& -> "
                               "result<Est, UE>` with exactly this F, Est and UE");

    any_solver(const any_solver&) = default;
    any_solver& operator=(const any_solver&) = default;   // replaces the whole value (strong guarantee)
    ~any_solver() = default;

    [[nodiscard]] result_type operator()(const F& f) const { return impl_(f); }

private:
    std::function<result_type(const F&)> impl_;
};

template<class T> inline constexpr bool is_any_solver_v = false;
template<class F, class Est, class UE> inline constexpr bool is_any_solver_v<any_solver<F, Est, UE>> = true;

namespace detail {
    template<class F, class Est, class UE>
    class runtime_first_of {   // owns its alternatives; never changed after construction
        std::vector<any_solver<F, Est, UE>> alts_;
    public:
        explicit runtime_first_of(std::vector<any_solver<F, Est, UE>> a) : alts_(std::move(a)) {}
        [[nodiscard]] result<Est, UE> operator()(const F& f) const {
            using R = result<Est, UE>;
            using Fail = failure<Est, UE>;
            if (alts_.empty()) return R{std::unexpect, Fail{errc::invalid_input, algo::none, {}, std::nullopt, {}}};
            std::optional<Fail> failed;   // merged so far: last code/cause, BEST estimate, total cost
            for (const auto& s : alts_) {
                R r = s(f);
                if (r) { if (failed) r->used = r->used + failed->used; return r; }   // lazy; success pays for failures
                failed = failed ? nxx::detail::merge(*failed, std::move(r).error()) : std::move(r).error();
            }
            return R{std::unexpect, *std::move(failed)};
        }
    };
}    // namespace detail

// Same shape as core.hpp's single-argument first_of(S) (one by-value parameter, `class` template head) so
// that this constrained overload is more specialised and wins for ranges of any_solver.
// Accepts std::vector, std::span, std::array, ... of any_solver; the chain owns a copy (it is curried,
// so borrowing the range would dangle). An empty range fails with errc::invalid_input when called.
template<class Rng>
    requires(std::ranges::input_range<Rng> &&
             is_any_solver_v<std::remove_cvref_t<std::ranges::range_reference_t<Rng>>>)
auto first_of(Rng alts) {
    using S = std::remove_cvref_t<std::ranges::range_reference_t<Rng>>;
    using Chain = detail::runtime_first_of<typename S::function_type, typename S::estimate_type, typename S::cause_type>;
    if constexpr (std::is_same_v<Rng, std::vector<S>>) return S{Chain{std::move(alts)}};
    else return S{Chain{std::ranges::to<std::vector<S>>(alts)}};
}

}    // namespace nxx
