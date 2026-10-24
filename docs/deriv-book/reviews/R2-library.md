# R2-library: targeted review of the revision (5432abf)

Reviewer: book-reviewer-library, round 2. Subject: commit `5432abf`, PDF `build/gcc/docs/deriv-book/pdf/main.pdf`
(88 pages). `include/` and `tests/` are unchanged since round 1, so the ground truth is the same.

Method: Chapter 8 read in full against the headers; 7.1, 7.2, 7.4, 7.6 to 7.9, 10.5, 10.6, Sidebar 3.3 and A.4 read
against the decision letter and the responses; all 25 examples compiled and run with GCC 16 and Clang 22 (llvm-mingw,
libc++), each with `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Werror`: every one builds
without a diagnostic and exits 0 on both. Scratch programs in
`C:\Users\kenne\AppData\Local\Temp\claude\C--Dev-XLThermo-Numerixx\5424a43d-ee2e-47d6-89e3-20a17882b9f6\scratchpad\r1lib`
(`ut2.cpp`, `acc.cpp`, `fstep.cpp`, `sigcheck.py`).

Summary: 21 of the 22 round-1 findings are resolved and 1 is partly resolved. New: 2 minor, 5 editorial, no major.

## Round-1 findings

| Id | Status | Note |
|---|---|---|
| 01 | resolved | Sidebar 3.3, A.4 and 8.2 (Preconditions, Failures) now say that the two-node rules fail only when both points coincide, and that one-sided movement near a power of two returns a value that is mostly rounding |
| 02 | resolved | `estimate-pitfalls` reads `*r` inside `if (r)` |
| 03 | resolved | 8.2: the assert is gone under NDEBUG, and `Forward<1>` with a negative step samples below x; 10.5 has the h > 0 hypothesis; the quick reference mentions NDEBUG |
| 04 | resolved | The E14 comment, 7.6 and Sidebar 7.1 now separate et_off from et_on |
| 05 | resolved | 7.2 and 8.2: the library throws nothing, exceptions pass through, not `noexcept` |
| 06 | resolved | 8.5 and 8.8: the weight check is a static_assert in the body, which a requires-expression does not see |
| 07 | resolved | 8.2: an untrusted result is also returned at the resolution of x |
| 08 | resolved | The budget claim is moved to 10.6, "for a selection of rules and orders, two-node rules included"; Table 10.5 says "of each family"; root_epsilon "for selected N" |
| 09 | resolved | 3.5, Sidebar 5.1 and 5.3 now share the generated factors 1.7, 5.5 and 8.0 |
| 10 | resolved | 8.3: `==` with a T compares the value only (false on an error); abs_error saturates at max() |
| 11, 12, 13, 14 | resolved | Wording as accepted; Table 4.1 drops its Ridders row and points to Table 8.2 |
| 15 | resolved (modified) | `tools/check.py` compares the signature lines with the headers; its granularity is R2-library-03 |
| 16 | partly | One two-line unbraced `if` is still in a quoted snippet: `03-rounding-step.cpp` lines 25-26, `explicit-step` (Chapter 3's first listing); the `return 1;` part was rejected by the editor |
| 17, 19, 20, 21, 22 | resolved | B.4, Sidebar 7.1, Table 9.10 with 9.5, Table 7.1, 7.2 and 8.3 as accepted |
| 18 | resolved | Table 8.1 now sits at the head of 8.1, and the lead-in ends with a full stop. Table 8.2 now floats into 8.6 instead (R2-library-06) |

Checked and correct in the revision:
- Chapter 8 statements: the precondition and NDEBUG; the resolution rule on both sides; the untrusted result at the
  resolution; static_assert against constraint; `==` with a T; saturation; `.value()`; the step conversion; member
  pointers of T only.
- The `numeric_limits` list. The members the headers read are is_specialized, is_integer, is_bounded, is_exact, radix,
  digits, epsilon, max and denorm_min (from a grep of `numeric_limits<T>::` in both headers), and that is exactly what
  `07-using-usertype.cpp` defines. I also drove that minimal type through `Central<10>` with N = 4, `Backward<6>` with
  N = 3, Ridders with N = 2 and 4, and a constexpr call (`ut2.cpp`); all compile and give e on GCC and on Clang/libc++.
- The pragma listing in 8.6 builds silently under `-Werror` on both compilers.
- 10.5's "no shared state" holds: the headers have only constexpr namespace-scope variables.
- Chapter 4's new claim checks out: Ridders exactly at the edge returns non_finite after 64 evaluations, which is 32
  failed levels of 2 each.
- `getting-started.cmake` matches the root CMakeLists: tests, demos and the book default to off when Numerixx is not
  the top-level project, there is no install rule, and the remote is right.
- The Newton listing returns non-convergence as an error.

## New findings

### R2-library-01 (minor): the Ridders acceptance test rejects every good result whose derivative is near 0

Location: `code/07-using-expected.cpp`, snippet `ridders-check` (7.2); `chapters/qr-quick-reference.tex`, "Two
checks" ("For Ridders, also require `*r->abs_error <= tol * std::abs(r->value)`").

Problem. The test is purely relative, so it fails whenever the value is 0 or tiny. That covers the book's own case
from three lines earlier, cos at 0, where the value is exactly 0. Section 7.5 tells the reader to judge such a
derivative by its absolute error; the check contradicts that.

Evidence (`acc.cpp`, the listing's `accept` with tol = 1e-8):

```
cos at 0: value 0 est 1.4e-14 accept(1e-8) false
cos at 1e-9: value -1.000e-09 exact -1.000e-09 est 2.8e-14 accept(1e-8) false
```

Suggested change. Add an absolute tolerance: `*res->abs_error <= atol + rtol * std::abs(res->value)`. Alternatively,
scale by the size of f (`tol * std::max(std::abs(res->value), fscale)`). Either way, add one sentence: near a zero
of f' use an absolute tolerance (7.5). The quick reference needs the same change.

### R2-library-02 (minor): the nested call for a mixed partial does not say how the inner failure propagates

Location: `chapters/07-using.tex` 7.9, "Mixed partials and shared evaluations".

Problem. The outer `derivative` differentiates a lambda that calls the inner one. That lambda must return a T, but
the inner call returns a `std::expected`. Every obvious way to unwrap it is wrong on a failure: `*inner` and
`inner->value` are undefined behaviour, and `inner.value()` throws, which 7.2 warns against.

Suggested change. One sentence: "return the inner value, or a NaN when the inner call fails
(`d ? d->value : std::numeric_limits<double>::quiet_NaN()`), so that the outer call reports non_finite".

### R2-library-03 (editorial): the signature check compares lines, not declarations

Location: `tools/check.py`, `signature_problems`.

Problem. The check does compare against the headers: unchanged it reports nothing, and a changed default (`auto N =
2`), a new member of `Ridders` and a new `error` value are each reported (`sigcheck.py`). But each line is matched on
its own against the whole header. A change confined to one of the two `derivative` overloads therefore passes
whenever the other overload still has the same line. Removing `[[nodiscard]] constexpr` from the overload without a
step reports nothing.

Suggested change. Join the `\code` lines of each declaration, up to the line that ends in `;` or `};`, and look up
that normalised string in the normalised header.

### R2-library-04 (editorial): "a double step is diagnosed" for a float x

Location: `chapters/07-using.tex` 7.1; `chapters/08-api.tex` 8.2, T.

Problem. Only a double that float cannot hold exactly is diagnosed, or a double variable. With `-Wconversion`, the
step `1e-3` gives 1 warning on GCC and 2 on Clang, and a double variable likewise. The step `0.25` gives none on
either compiler (`fstep.cpp`).

Suggested change. "a double step is diagnosed under -Wconversion unless float holds it exactly".

### R2-library-05 (editorial): one multi-line unbraced `if` remains in a quoted snippet

Location: `code/03-rounding-step.cpp` 25-26 (`explicit-step`, the first listing of Chapter 3).

Problem. This is the residue of R1-library-16: the editor narrowed style decision 18 to multi-line `if`s, and this one
still spans two lines.

Suggested change. Brace it.

### R2-library-06 (editorial): Table 8.2 splits a sentence of 8.6

Location: PDF p. 66 to 67 (book pp. 58-59).

Problem. The float `tab:api-methods`, referenced in 8.5, is set inside 8.6, between "bfloat16," and "not double),".

Suggested change. Use `[t]` placement, or put the table before the paragraph that ends 8.5.

### R2-library-07 (editorial): `getting-started.cmake` pins the stage-4 commit by a short hash

Location: `code/getting-started.cmake`, `GIT_TAG aeec7e6`; 7.1.

Problem. The pin is `aeec7e6`, the integration commit, not the commit that the version statement will name at
release. It is also a short hash, and CMake's documentation asks for the full hash. In addition, the library needs
CMake 3.25 (its `cmake_minimum_required`), which the listing does not say.

Suggested change. At release, write the full hash of the commit that the front matter names, and add the CMake
version to the listing's comment.

## Verdict

The revision fixed what round 1 found, and the new material is correct. Part II now states the precondition, the
resolution rule, the outcomes of Ridders' method, the constraint/static_assert split, comparisons, saturation and the
numeric_limits members as the code has them. The minimal user type compiles and runs on both compilers along every
path I exercised. Two things should be fixed before release: the Ridders acceptance test (R2-library-01), which
rejects correct results near a zero derivative, and the nested-call advice for mixed partials (R2-library-02). The
rest is editorial. With those fixes I sign off on the library content.
