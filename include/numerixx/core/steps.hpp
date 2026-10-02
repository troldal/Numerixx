// Manual stepping and tracing (DESIGN §6.9): steps_view, a lazy input range over the solver protocol.
//
//     const auto solver = nxx::roots::brent{};
//     auto p = solver.prepare(std::cref(fn), nxx::bracket{1.0, 2.0});
//     if (p)
//         for (const auto& st : nxx::steps_view{solver, *p} | std::views::take(8))
//             if (st) trace(solver.estimate(*st).x);
//
// Element 0 is init(p), its failure mapped to a fault. The range ends after the first error or the first intrinsic
// stop, and is otherwise infinite: bound it with std::views::take. It applies neither the stop criterion nor the
// budget; those belong to nxx::iterate. fn must outlive the view, because the problem holds std::cref(fn).
#pragma once

#include <numerixx/core/error.hpp>
#include <numerixx/core/iterate.hpp>

#include <cstddef>
#include <expected>
#include <iterator>
#include <optional>
#include <ranges>
#include <utility>

NXX_BEGIN_HEADER

namespace nxx
{
    template<class A, class P>
    class steps_view : public std::ranges::view_interface<steps_view<A, P>>
    {
        using S  = detail::state_t<A, P>;
        using UE = typename detail::init_failure_t<A, P>::cause_type;

    public:
        using element = std::expected<S, fault<UE>>;

        constexpr steps_view(A alg, P problem) : alg_(std::move(alg)), p_(std::move(problem)) {}

        class iterator
        {
            const steps_view*      parent_ = nullptr;
            std::optional<element> cur_ {};
            bool                   done_ = false;
            friend class steps_view;

        public:
            using value_type       = element;
            using difference_type  = std::ptrdiff_t;
            using iterator_concept = std::input_iterator_tag;

            iterator() = default;

            constexpr const element& operator*() const { return *cur_; }

            constexpr iterator& operator++()
            {
                if (!*cur_ || parent_->alg_.intrinsic(**cur_))
                    done_ = true;
                else
                    cur_ = parent_->alg_.step(parent_->p_, **cur_);
                return *this;
            }

            constexpr void operator++(int) { ++*this; }

            friend constexpr bool operator==(const iterator& it, std::default_sentinel_t) { return it.done_; }
        };

        constexpr iterator begin() const
        {
            iterator it;
            it.parent_ = this;
            auto s     = alg_.init(p_);
            if (s)
                it.cur_.emplace(*std::move(s));
            else
                it.cur_.emplace(std::unexpect, fault<UE> { s.error().code, s.error().used.evaluations, s.error().cause });
            return it;
        }

        constexpr std::default_sentinel_t end() const noexcept { return {}; }

    private:
        A alg_;
        P p_;
    };
}    // namespace nxx

NXX_END_HEADER
