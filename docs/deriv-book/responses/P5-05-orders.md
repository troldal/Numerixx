# P5-05-orders: the Phase 5 update (PLAN 12.8)

Author: `book-author-methods`, Phase 5 update. Files: `chapters/05-orders.tex`, `figures/05-orders-order-accuracy.tex`,
`code/05-orders-float.cpp`. `tools/check_chapter.sh`: `COMPILED OK`, 5 pages standalone (unchanged).

| Item (PLAN 12.8) | What changed |
|---|---|
| 5.2 "why float cannot afford a fourth derivative at all" | "why float cannot afford a derivative above the first" |
| Sidebar 5.1, "only the fourth carries the warning" | "The second and third derivatives keep more, `\valFloatSecondCentralRelError` and `\valFloatThirdCentralRelError` with Central<2>, but far less than the `\valFloatFirstCentralRelError` of f', which is why every derivative above the first carries the warning of 5.4" |
| 5.4 The float warning ("a call with N = 4") | Rewritten: every derivative above the first, by any method, in a type with epsilon above 1e-7; the medians `\valFloat{First,Second,Third,Fourth}CentralRelError`, `\valFloatSecondRiddersRelError` to `\valFloatFourthRiddersRelError`, `\valFloatSecondForwardRiddersRelError` to `\valFloatFourthForwardRiddersRelError` (pointing to 6.7); the message quoted verbatim from `deriv.hpp` lines 1810/1831; silencing as before (8.6); Ridders' estimates `\valFloatSecondRiddersRelEstimate` and `\valFloatFourthRiddersRelEstimate` |
| Figure 5.2 and its caption | The accent box now encloses the float points for N = 2 to 4 ("float f'' to f'''': warning"); the "float, dashed" label moved to the upper left; the double labels moved right of the box. Caption: "the float points for f'' to f'''' boxed as the cases that warn" |
| E8 (coordinator, from the practice author) | `code/05-orders-float.cpp`, snippet `float-fourth`: "its warning says it keeps few correct digits" (the new message) instead of "mostly rounding noise"; the snippet keeps its 4 comment lines, each within 60 characters; the file comment updated. Rebuilt with `-Werror` and run: exit 0, output unchanged |
