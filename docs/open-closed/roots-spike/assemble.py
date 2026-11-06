"""Writes the spike's headers: roots.hpp with the method protocol, math.hpp with nxx::observer, and quadrature.hpp
taking the observer from math.hpp (docs/roots/PHASE7_API.md).

    python3 docs/open-closed/roots-spike/assemble.py <out>      # from the project root; writes <out>/numerixx/*.hpp

roots.hpp keeps today's header where it is unchanged (the result types, the fn-shape concepts, validation, the
safe arithmetic, the interpolation formulas, search()) and replaces the closed lists, the dispatch and the
detail loops with chunk_public.hpp (Tolerance, Start, Step, values_of), chunk_states.hpp (each method's state and
step) and chunk_methods.hpp (the methods, the concepts, the loop, observing and the forms of solve()).
"""
import re
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
root = here.parent.parent.parent
out = Path(sys.argv[1]) / "numerixx"
out.mkdir(parents=True, exist_ok=True)

src = (root / "include/numerixx/roots.hpp").read_text(encoding="utf-8").split("\n")


def lines(a, b):  # 1-based, inclusive
    return src[a - 1:b]


def find(pattern, start=1):
    for i in range(start - 1, len(src)):
        if re.search(pattern, src[i]):
            return i + 1
    raise LookupError(pattern)


def doc_start(line):  # the line of the '/**' that opens the comment above `line`
    i = line - 1
    while not src[i - 1].strip().startswith("/**"):
        i -= 1
    return i


def body_end(line, close="        }"):
    i = line
    while src[i - 1] != close:
        i += 1
    return i


chunk = lambda name: (here / name).read_text(encoding="utf-8").split("\n")

out_lines = []
out_lines += lines(1, find(r"^    };", find(r"struct Criteria")))
out_lines += chunk("chunk_public.hpp")
out_lines += ["", "    /** @cond INTERNAL */", "    namespace detail {", ""]
# The fn-shape concepts, largest, iteration_cap.
tolerance = find(r"^        struct Tolerance")
out_lines += lines(doc_start(find(r"inline constexpr bool is_standard_tuple = false")), doc_start(tolerance) - 1)
# twice, checked_tolerance, same_sign, the safe arithmetic, evaluate, derivative_value, Joined.
out_lines += lines(body_end(tolerance, "        };") + 1, body_end(find(r"^        struct Joined"), "        };"))
out_lines += ["", "        /** @brief A validated bracket, or a root found during validation (an end, or 0). */",
              "        template<real T, std::size_t K>", "        struct Prepared", "        {",
              "            /** @brief The root, if validation found one. */", "            std::optional<Root<T>> root;",
              "            /** @brief The bracket. */", "            Start<T, K> start;", "        };", ""]
prepare = find(r"prepare\(F &fn, Bracket")
block = "\n".join(lines(doc_start(prepare - 2), body_end(prepare)))
block = block.replace("constexpr std::expected<Start<T, K>, error>", "constexpr std::expected<Prepared<T, K>, error>")
block = re.sub(r"Start<T, K>\{ (Root<T>\{[^}]*\}|std::nullopt), lower, at_lower, upper, at_upper \}",
               lambda m: "Prepared<T, K>{ " + m.group(1) + ", Start<T, K>{ lower, at_lower, upper, at_upper } }", block)
out_lines += block.split("\n")
for first, last in [
    (r"constexpr std::optional<T> newton_step", r"constexpr std::optional<T> halley_step"),
    (r"constexpr std::optional<T> difference", r"constexpr std::optional<T> sum\("),
    (r"^        struct Enclosure", r"constexpr T place\("),
    (r"constexpr std::optional<T> secant_point", r"constexpr std::optional<T> double_secant"),
    (r"constexpr T sqrt\(", r"ridders_point\(Enclosure"),
    (r"constexpr T secant_second_point", r"secant_move\(T const &x"),
]:
    out_lines += [""] + lines(doc_start(find(first)), body_end(find(last)))
out_lines += chunk("chunk_states.hpp")
out_lines += ["", "    } // namespace detail", "    /** @endcond */", "", "    using nxx::Observer;",
              "    using nxx::observer;", "",
              "    // SPIKE ONLY: one test names detail::Tolerance; the phase moves it to the public name.",
              "    namespace detail {", "        using nxx::roots::Tolerance;", "    }", ""]
out_lines += chunk("chunk_methods.hpp")
# search() and the formatter, unchanged.
i = find(r"struct SearchSide")
while src[i - 1].strip() != "/** @cond INTERNAL */":
    i -= 1
out_lines += lines(i, len(src))

text = "\n".join(out_lines)

# The numerics review of draft 1, finding 1: in a type without subnormals, (b - a) / 2 flushes to 0 where b - a is
# below 2 min, and the midpoint was an end although numbers lie between. Unchanged wherever a + (b - a) / 2 is inside.
old_midpoint = """            if ((a < T(0)) != (b < T(0))) return T(T(a + b) / T(2));
            return T(a + T(T(b - a) / T(2)));"""
assert old_midpoint in text
text = text.replace(old_midpoint, """            if ((a < T(0)) != (b < T(0))) return T(T(a + b) / T(2));
            T const m = T(a + T(T(b - a) / T(2)));
            if ((a < m && m < b) || (b < m && m < a)) return m;
            // Ends of one sign with no point formed between them: adjacent (then none exists), or, in a type without
            // subnormals, closer than 2 min, where (b - a) / 2 flushes to 0. Then both are tiny and their sum fits; the
            // guard keeps the sum from overflowing near max, which constant evaluation would reject.
            T const quarter = T(detail::largest<T>() / T(4));
            if (!(math::abs(a) <= quarter && math::abs(b) <= quarter)) return m;
            T const sum = T(T(a + b) / T(2));
            if ((a < sum && sum < b) || (b < sum && sum < a)) return sum;
            return m;""", 1)

old_not_converged = """         * An open method found no root from its guess: no sign change between its iterates within its iteration cap,
         * or no usable step"""
assert old_not_converged in text
text = text.replace(old_not_converged, """         * A method found no root. From a guess: no sign change between its iterates within its iteration cap, or no
         * usable step""", 1)
text = text.replace("""         * is as small as the derivative's own error (a multiple root): see Newton with separate derivatives.
         */
        not_converged""", """         * is as small as the derivative's own error (a multiple root): see Newton with separate derivatives. On a
         * bracket: a method that could place no point inside it, which the library's methods do not reach.
         */
        not_converged""", 1)

text = text.replace("""        not_converged
    };""", """        not_converged,

        /**
         * A method broke the protocol's contract: a step outside the caller's bracket or with its point outside its
         * points, a bracketing step without a sign change, or a result outside the caller's bracket. Never returned by
         * the library's methods.
         */
        invalid_step
    };""", 1)
text = text.replace("""            case nxx::roots::error::not_converged:
                name = "not_converged";
                break;""", """            case nxx::roots::error::not_converged:
                name = "not_converged";
                break;
            case nxx::roots::error::invalid_step:
                name = "invalid_step";
                break;""", 1)
(out / "roots.hpp").write_text(text, encoding="utf-8")

# math.hpp: nxx::Observer and nxx::observer, moved from nxx::quadrature.
m = (root / "include/numerixx/math.hpp").read_text(encoding="utf-8")
m = m.replace("#include <concepts>\n", "#include <concepts>\n#include <functional>\n", 1)
observer = '''
    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief How an Observer stores its callable: an lvalue by reference, so that a recorder the caller owns sees
         * every step; an rvalue by value, so that a temporary lambda lives as long as the Observer.
         * @tparam O The callable's type as given to observer().
         */
        template<typename O>
        using kept = std::conditional_t<std::is_lvalue_reference_v<O>,
            std::reference_wrapper<std::remove_reference_t<O>>, std::remove_cvref_t<O>>;

    } // namespace detail
    /** @endcond */

    /**
     * @brief Marks a callable as the observer of a stepped method, so that integrate() and solve() tell it from fn by
     * its type and never by trying to call an argument: made by observer(), passed before fn.
     * @tparam O The callable's type as given to observer().
     */
    template<typename O>
    struct Observer
    {
        /** @brief The callable: a std::reference_wrapper for an lvalue, the callable itself for an rvalue. */
        detail::kept<O> callable;
    };

    /**
     * @brief Makes @p o the observer of a stepped method: `o(step)` is called with the module's Step after the
     * method's start and after each step.
     * @tparam O The callable's type.
     * @param o The callable.
     * @return The observer.
     */
    template<typename O>
    constexpr Observer<O> observer(O &&o)
    { return Observer<O>{ detail::kept<O>(std::forward<O>(o)) }; }

} // namespace nxx'''
assert "    } // namespace math\n} // namespace nxx" in m
m = m.replace("    } // namespace math\n} // namespace nxx", "    } // namespace math\n" + observer, 1)
(out / "math.hpp").write_text(m, encoding="utf-8")

# quadrature.hpp: its own kept, Observer and observer() replaced by using-declarations.
q = (root / "include/numerixx/quadrature.hpp").read_text(encoding="utf-8")
i = q.index("        /**\n         * @brief How an Observer stores its callable")
end = "std::remove_cvref_t<O>>;\n"
j = q.index(end, i) + len(end)
q = q[:i] + "        using nxx::detail::kept;\n" + q[j:]
i = q.index("    /**\n     * @brief Marks a callable as the observer of an integration")
end = "    { return Observer<O>{ detail::kept<O>(std::forward<O>(o)) }; }\n"
j = q.index(end, i) + len(end)
q = q[:i] + "    using nxx::Observer;\n    using nxx::observer;\n" + q[j:]
(out / "quadrature.hpp").write_text(q, encoding="utf-8")
print(f"wrote {out}")
