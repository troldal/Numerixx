# code/

Every code listing in the book is a file here: built as `numerixx_book_<file>` and run as a CTest test with the
normal build, so a listing cannot stop compiling or drift from the library. `quickstart.cpp` shows the conventions:

- the snippet the book quotes lies between the comment lines `// <name>` and `// </name>` and is quoted with
  `\snippet{<file>}{<name>}`;
- markers and snippet start in column 0 (a function or declaration at namespace scope, or statements between
  `// clang-format off` and `// clang-format on`), because `\lstinputlisting` cannot remove indentation;
- snippet lines stay within 60 characters, the width of a listing in the 114 mm text block;
- `main` checks that the library does what the text says (exit status 0), so the test fails if it stops doing so.
