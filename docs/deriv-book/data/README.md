# data/

Written only by `tools/generate_data.cpp` (CMake target `numerixx_book_data`, also run by `numerixx_book`); every file
carries the generator's header, and `tools/check.py` rejects one that does not. Never edit these files by hand: change
the generator and rerun it. The generated files are committed, so that a diff shows when a number in the book changes.
The data sets are DS1 to DS21 of `PLAN.md` section 4.3; the generator exits with status 1 if one of its self-checks
fails (the thread function against the evidence reviewer's exact values, a moment condition, Ridders' returned entry
not found in the tableau, a `\val` macro defined twice).

## Conventions

- `values.tex`: one `\val<Name>` macro per number quoted in the text, read by `macros.tex`; names are letters only.
  Errors, estimates and steps have two significant digits, counts are integers, coverage is a percentage with one
  decimal; typeset them with `\num{\valX}` or `\SI{\valX}{unit}`. Naming: `<X>RelError` is |value - exact| / |exact|;
  `<X>Actual` is |value - exact|; `<X>Estimate` is the returned `abs_error` (absolute, in the derivative's unit);
  `<X>RelEstimate` is `abs_error` / |exact|; `<X>Evals` the number of evaluations of fn; `<X>Value` the returned value
  to twelve digits; `<X>Step` a step in the argument's unit (kelvin for the thread).
- `<figure>.csv`: one file per plotted data set; `#` comment lines first (the header, then the data set's description),
  then the column names, then the rows, read with `\addplot table[col sep=comma, comment chars={\#}]`. A result that is
  not finite is written as `nan`, which pgfplots skips (set `unbounded coords=jump` to break the curve there). Steps are
  10^(k/4), four per decade. Errors are relative unless the column name says otherwise.
- `table-<name>.tex`: the body rows of a `tabular` (cells separated by `&`, every row ending in `\\`), the generator's
  header and the column list as `%` comments on the first lines. Read a body with `\tablebody{<name>}` (`macros.tex`)
  on a line of its own between `\midrule` and `\bottomrule` of a `tabular` whose header the chapter writes; not with
  `\input`, whose file hooks add tokens after the file and give "Misplaced \noalign". Numbers in scientific notation
  are wrapped in `\num{}`; a cell that could not be computed is `--`.
- Reproducibility: the generator forms powers by multiplication and the step grid from `sqrt` and exact powers of
  ten, never with libm's `pow`, whose results differed by an ulp between runs on the same machine (UCRT); the libm
  calls that remain (exp, sin, log, tanh) are shared by fn and the reference value. `data/` is generated on the gcc
  preset (GCC 16, UCRT); another libm may change the last bits and the exact-zero counts of DS9, which the committed
  files and their diff make visible.
- The thread: nitrogen's Shomate enthalpy from `code/shomate.hpp` (range 500 to 2000 K) at T0 = 1000 K; H in J/mol
  (relative to 298.15 K), T and steps in K, cp in J/(mol K). exp(x) is the function of the method figures.
- "Default-step accuracy" sets (DS6, DS7, DS15) report the median over 11 points x = 0.5, 0.55, ..., 1 (same scale,
  same step), with the smallest and largest value alongside; a single point's rounding error is luck by a factor of ten.
- The random study (DS9, DS19's coverage) uses `std::mt19937` with the seed 20261023 (plus the function's index) and
  x = 0.2 + 2.8 * u32 / 2^32 from the raw 32-bit output, so every standard library gives the same 300 points.
- Parameters an example shares with the generator (the example must use the same): E13 Arrhenius rate A = 1e10 1/s,
  Ea = 80000 J/mol, T = 600 K, R = 8.31446261815324 J/(mol K); E19 delta = 1e-3 K (inside) and 50 K (outside) from the
  500 K junction; E20 Colebrook equation 1/sqrt(f) = -2 log10(eps_r / 3.7 + 2.51 / (Re sqrt(f))) with eps_r = 1e-4,
  Re = 1e5, bisection on [1, 20] in x = 1/sqrt(f) to the relative tolerance eta = 1e-10; E17 noise eta = 1e-8 through
  `book::noisy` of `code/noise.hpp`.

## Files and macros, by figure or chapter

| Set | File(s) | Macros (`\val...`) | Used by |
|---|---|---|---|
| hand | `table-hand-cube.tex` (x^3 at 1, h = 0.1 and 0.05: f(1+h), f(1-h), Forward<1> and Central<2> values and errors to 7 digits), `table-hand-calculator.tex` (a six-digit calculator on exp at 1, h = 0.1 to 0.0001: values rounded to 6 digits, estimate, signed error) | `HandCalcExact` (e to 6 digits), `HandCalcErrorRatio`, `HandCalcLastDiff` | the pencil-and-paper tables of 2.3 and 3.2 |
| DS1 | `step-error-cp.csv` (`h, forward1, central2, central4, model_forward1, model_central2`) | | F2 `fig:fd-slopes` (ch. 2), F3 `fig:v-curve` (ch. 3) |
| DS1b | `values.tex` | `CpExact`, `CpExactLong`, `EnthalpyExact`, `CpDefaultStep`, `CpDefaultRelError`, `CpDefaultValueLong`, `CpBestStep`, `CpBestRelError`, `CpCentralFour{RelError,Step}`, `CpCentralSix{RelError,Step}`, `CpRidders{RelError,Actual,Estimate,RelEstimate,Evals,Value,InitialStep}`, `CpSecond{Exact,Step,RelError}`, `CpRiddersSecond{RelError,Actual,Estimate,RelEstimate,Evals,Value}`, `CpSecondModelRelError` (4 eps |H| / (h^2 |H''|) at the default step), `CpRelativeFactor`, `CpSecondRelativeFactor` | ch. 1 (1.4, 1.5), ch. 3 (3.2 to 3.4), ch. 4 (E5), ch. 5 (E7, 5.2), ch. 6 (E9); F3's default-step marker |
| E6 | `values.tex` | `EdgeForwardTwoRelError`, `EdgeForwardFourRelError` (Forward<2> and <4> at 500 K with the NaN wrapper of `code/04-stencils-edge.cpp`) | 4.4 (ch. 4) |
| DS2 | `step-error-exp.csv` (`h, forward1, forward2, forward4, central2 ... central10`) | | F6 `fig:order-ladder`, F7 `fig:one-sided` (ch. 4) |
| DS3 | `step-error-float.csv` (`h, error_float, error_double`) | `Float{Best,Default}{RelError,Step}`, `DoubleBest{RelError,Step}` | F5 `fig:float-double` (ch. 3, 3.7) |
| DS4 | `step-scale.csv` (`h, error`) | `Scale{Default,Fit,Best,ModelBest}Step`, `Scale{Default,Fit,Best}RelError`, `ScaleRiddersUser{RelError,Actual,Estimate,RelEstimate,Evals,Value,Step}`, `ScaleRiddersDefault{...}` | F4 `fig:step-scale`, E4 (ch. 3, 3.5) |
| DS4b | `step-scale-small.csv` (`h, error`) | `SmallDefaultStep`, `SmallDefaultStatus` (the word "finite" or "non-finite"), `SmallFitStep`, `SmallFitRelError`, `SmallRiddersUser{...}` | F4b `fig:step-scale-small`, E4 (ch. 3, 3.5) |
| DS5 | `step-error-orders.csv` (`h, n1, n2, n3, n4`) | | F10 `fig:v-curve-orders` (ch. 5) |
| DS6 | `order-accuracy.csv` (`N, c2_double, c4_double, rid_double, rid_est_double, rid_evals_double, c2_float, c4_float, rid_float, rid_est_float, rid_evals_float`; medians over 11 points) | `FloatThirdRiddersRelError`, `FloatFourth{Central,Ridders}RelError`, `FloatFourthRidders{RelEstimate,Evals}`, `DoubleFourth{Central,Ridders}RelError`, `DoubleFourthRiddersRelEstimate`; the single point x = 1: `FloatFourthRiddersAtOne{RelError,Estimate,RelEstimate}`, `FloatFourthCentralAtOneRelError` | F11 `fig:order-accuracy`, 5.4, E8 (ch. 5) |
| DS7 | `default-step-accuracy.csv` (`P, central, central_min, central_max, forward, forward_min, forward_max`; medians) | `ExpCentral{Two,Four,Six,Ten}RelError`, `ExpForward{One,Two,Six}RelError`, `ExpCentralTwoAtOneRelError` | F8 `fig:default-step-accuracy` (ch. 4, 4.2) |
| DS8 | `ridders-tableau.csv` (`i, j, h, A, error, estimate, floor`: relative error, relative estimate e_ij / cp with `nan` for j = 0, relative floor r_i / cp; 16 levels), `table-ridders-tableau.tex` (values, the 7 levels Ridders used, returned entry in `accent`), `table-ridders-tableau-errors.tex` (relative errors and the floor, same levels) | `Tableau{LevelsUsed,StopLevel,BestLevel,BestColumn,Levels,BaseBestLevel,BaseBestRelError,BaseFirstRelError,BestRelError}` | F12 `tab:ridders-tableau`, F13 `fig:ridders-convergence` (ch. 6, 6.2, 6.3) |
| DS9 | `ridders-coverage-N1.csv` .. `-N4.csv` (`fn, name, x, exact, actual, estimate, rel_actual, rel_estimate, evaluations, restarted`; absolute `actual`/`estimate`, relative `rel_*`); `ridders-coverage-plot-N1.csv` .. `-N4.csv` (`actual, estimate` only, every third case, an exact 0 clipped to eps times the exact derivative: what the scatter F14 draws, since pgfplots keeps about 480 words of TeX's memory per plotted point and 5400 points exhausted it; `\valRiddersPlotEvery`, `\valRiddersPlotCases`) | per N (`First` .. `Fourth`): `RiddersCoverage`, `RiddersUncovered`, `RiddersMedianRatio`, `RiddersWorstRatio`, `RiddersEvals{Min,Median,Max}`, `RiddersExactCases`, `RiddersRestarted`, `RiddersFailed`; `RiddersCases`, `RiddersFunctions`, `RiddersPointsPerFunction`, `RiddersExactCases` | F14 `fig:ridders-coverage`, F14b `fig:ridders-cost` (N = 1, 2; ch. 6, S11), 6.5, ch. 10 (10.2); T15 (N = 1 to 4) |
| DS10 | `ridders-edge.csv` (`k, x, ok, covered, rel_error, rel_estimate, evaluations`; last row x = 0.01 with k = log2(1/x)) | `RiddersEdge{Limit,LimitExponent,TenDigits,TenDigitsExponent,CapExponent}`, `RiddersEdgeAtOneHundredth{RelError,Actual,Estimate,RelEstimate,Evals,Value}`, `RiddersEdgeEvalsAtOneHundredth`, `EdgeCentralAtOneHundredthRelError` | F15 `fig:ridders-edge`, E10 (ch. 6, 6.4), ch. 10 (10.4) |
| DS11 | `newton.csv` (`iteration, T, residual, error`; residual in J/mol, error in K) | `Newton{Iterations,FinalError,FinalResidual,Start}` | F16 `fig:newton`, E12 (ch. 7, 7.4) |
| DS12 | `table-weights-central-N1.tex` .. `-N4.tex`, `table-weights-forward-N1.tex` .. `-N4.tex` (columns: P, nodes, numerators, d, S, C, leading error term) | `Ampl{Central,Forward}<P>`, `Errc{Central,Forward}<P>` with P in words (`Two`, `Six`) for N = 1; `AmplForwardSixFourth` (the largest S of all the tables) | T4 to T11 (ch. 9, 9.2, 9.3); S and C quoted in ch. 4 |
| DS13 | `table-default-steps.tex` (family, P, N, exponent, step for double at s = 1, model error, evaluations) | `RootEps{Three,Four,Five,Six}Double`, `RootEpsThreeFloat`, `ModelCentralTwo{First,Second,Third,Fourth}Double`, `ModelCentralTenFirstDouble`, `ModelCentralTwo{First,Fourth}Float`, `ModelForwardOne{First,Fourth}Double`, `StepConstantCentralTwo{First,Second,Fourth}`, `StepConstantCentralTenFirst`, `StepCost{CentralTwoFirst,CentralTwoSecond,CentralTwoFourth,CentralFourFourth,CentralTenFirst}` (the model's cost of dropping the constant, App. B.4) | T12 (ch. 9, 9.4); ch. 3 (3.3, 3.4, S5), ch. 5 (5.2, S8), App. A |
| DS14 | `table-ridders-constants.tex` (constant, value) | `RiddersLevelCap`, `RiddersDivisor{First,Higher}`, `RiddersMaxEvals{First,Second,Third,Fourth}` | T13 (ch. 9, 9.5), S11, ch. 10 |
| DS15 | `accuracy-matrix.csv` (`method, type, N, median, min, max`), `table-accuracy-matrix.tex` (method; N = 1..4 for float, double, cpp_bin_float_50; float N = 4 marked `*`) | `AccuracyMatrixHasMultiprecision`, `AccuracyMatrixPoints`, `Headline{First,Fourth}Digits`, `HeadlineRidders{First,Fourth}Digits` | T14 (ch. 10, 10.1); 1.5 |
| DS16 | `table-cost.tex` (method; evaluations per N; Ridders median with range, then the caps) | | T15 (ch. 10, 10.3); 1.5 |
| DS17 | `values.tex` | `Eps{Float,Double,Multiprecision}`, `Max{Float,Double}`, `DenormMin{Float,Double}`, `Min{Float,Double}`, `Digits{Float,Double,Multiprecision}`, `DigitsTen{Float,Double,Multiprecision}`, `UlpAtOneThousand`, `UlpAtOne`, `UlpBelowOne` | App. A (T16), S3, ch. 3 (3.1) |
| DS18 | `values.tex` | E13: `SensitivityRate{RelError,Actual,Estimate,RelEstimate,Evals,Value}`, `SensitivityEa{...}`, `SensitivityK`, `SensitivityScaled`, `SensitivityScaledExact`, `SensitivityScaledRelError`, `JacobianEvals`, `JacobianEvals{A,Ea}`, `Jacobian{A,Ea}RelError`; E14: `MpCentralSix{RelError,Step}`, `MpCentralTwoRelError`, `MpRidders{RelError,Estimate,RelEstimate,Evals}`; E18 (recomputed as `code/07-using-check.cpp` computes them): `Check{Steps,Exact,Order}Digits` (agreeing digits, one decimal), `CheckRiddersDiff`, `CheckStencilActual`, `CheckRiddersEstimate`, `CheckSine{Derivative,AbsError,RelError,Spread,DigitsDerivative,DigitsScale}` (sin x at x = 1.5708) | ch. 7 (7.4, 7.5, 7.6) |
| DS19 | `step-error-noise.csv` (`h, clean, noise_ten, noise_eight`) | `NoiseClean{BestStep,BestRelError}`, `Noise{Ten,Eight}{BestStep,BestRelError,FitStep,FitRelError}`, `NoiseEightDefaultRelError`, `NoiseEightModelBest{Step,RelError}`, `NoiseRidders{RelError,Actual,Estimate,RelEstimate,Evals,Value}` (x = 1), `NoiseRiddersEstimateRatio`, `Noise{Clean,Ten,Eight}Ridders{Coverage,MedianRatio,WorstRatio,MedianRelError}` over `NoiseRiddersCases` (300) points | F17 `fig:noise`, 3.6, E17 (ch. 3); 6.4 |
| DS20 | `kink.csv` (`delta, central2, central4, ridders, ridders_estimate, ridders_evals`; delta in K from the 500 K junction) | `KinkJunction`, `KinkJump{Enthalpy,Cp,CpSlope,CpSlopeRel}`, `KinkCpAtJunction`, `KinkCentralStep`, `KinkCentralFourStep`, `KinkRiddersInitialStep`, `Kink{Inside,Outside}Delta`, `Kink{Inside,Outside}{Central,CentralFour}RelError`, `Kink{Inside,Outside}Ridders{RelError,Actual,Estimate,RelEstimate,Evals,Value}`, and the plan's `KinkCentralRelError`, `KinkRiddersRelError`, `KinkRiddersEstimate` (the inside point) | F18 `fig:kink`, E19, 7.8 (ch. 7) |
| DS21 | `values.tex`, `step-error-solver.csv` (`h, bisection, fixed_point`; h in units of Re) | `Solver{Reynolds,Roughness,Tolerance,Friction,FrictionNoise,Exact,ScaleFactor,DefaultStep,FitStep}`, `Solver{Default,Fit,ModelFit}RelError`, `SolverDeriv{RelError,Actual,Estimate,RelEstimate,Evals,Value,EstimateRatio}` (Ridders on the bisection solver), `SolverFixedPoint{Noise,DefaultRelError,FitRelError}`, `SolverFixedPointRidders{...}` | E20, 7.4 (ch. 7); 3.6 |
