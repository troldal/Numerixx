// /*
//     888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8  8b        d8
//     8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P    Y8,    ,8P
//     88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8'
//     88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
//     88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b          d88b
//     88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,      ,8P  Y8,
//     88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b    d8'    `8b
//     88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8  8P        Y8
//
//     Copyright © 2022 Kenneth Troldal Balslev
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy
//     of this software and associated documentation files (the “Software”), to deal
//     in the Software without restriction, including without limitation the rights
//     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//     copies of the Software, and to permit persons to whom the Software is furnished
//     to do so, subject to the following conditions:
//
//     The above copyright notice and this permission notice shall be included in all
//     copies or substantial portions of the Software.
//
//     THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// */
//
// #ifndef NUMERIXX_MULTIROOTS_IMPL_HPP
// #define NUMERIXX_MULTIROOTS_IMPL_HPP
//
// // ===== Numerixx Includes
// #include "MultiDerivatives.hpp"
// #include <Constants.hpp>
//
// // ===== External Includes
// #include <blaze/Blaze.h>
//
// // ===== Standard Library Includes
// #include <functional>
// #include <iostream>
// #include <stdexcept>
//
// namespace nxx::multiroots
// {
//
//     namespace impl
//     {
//         template< typename... >
//         struct FunctionTraits;
//
//         template< typename RET, typename ARG >
//         struct FunctionTraits< std::function< RET(ARG) > >
//         {
//             using return_type    = RET;
//             using container_type = ARG;
//             using argument_type  = typename std::remove_cvref_t<ARG>::value_type;
//         };
//
//     }    // namespace impl
//
//     // clang-format off
//     template< typename FunctionType >
//     concept IsMultirootFunction =
//         std::same_as< typename impl::FunctionTraits< FunctionType >::return_type, typename impl::FunctionTraits< FunctionType >::argument_type > &&
//         std::same_as< FunctionType, std::function< typename impl::FunctionTraits< FunctionType >::return_type(typename impl::FunctionTraits< FunctionType >::container_type) > > &&
//         std::floating_point< typename impl::FunctionTraits< FunctionType >::return_type > &&
//         std::floating_point< typename impl::FunctionTraits< FunctionType >::argument_type >;
//     // clang-format on
//
//     template< typename T >
//     //requires IsMultirootFunction< FunctionType >
//     class DMultiNewton;
//
//     /*
//      * Forward declaration of the PolishingTraits class.
//      */
//     template< typename SOLVER >
//     struct MultirootsSolverTraits;
//
//     /*
//      * Specialization of the PolishingTraits class for Newton<FN, DFN>
//      */
//     template< typename T >
//     struct MultirootsSolverTraits< DMultiNewton< T > >
//     {
//         //using function_type  = FunctionType;
//         using return_type    = T;//typename impl::FunctionTraits< FunctionType >::return_type;
//         //using container_type = typename impl::FunctionTraits< FunctionType >::container_type;
//         //using argument_type  = typename impl::FunctionTraits< FunctionType >::argument_type;
//     };
//
//     // ========================================================================
//     // MULTIROOT-FINDING WITH DERIVATIVES
//     // ========================================================================
//
//     template< typename SUBCLASS >
//     //requires IsMultirootFunction< typename MultirootsSolverTraits< SOLVER >::function_type >
//     class MultirootBase
//     {
//         /*
//          * Friend declarations.
//          */
//         friend SUBCLASS;
//
//     public:
//         //using function_type  = typename MultirootsSolverTraits< SOLVER >::function_type;
//         using return_type    = typename MultirootsSolverTraits< SUBCLASS >::return_type;
//         //using container_type = typename MultirootsSolverTraits< SOLVER >::container_type;
//         //using argument_type  = typename MultirootsSolverTraits< SOLVER >::argument_type;
//
//     private:
//         //std::vector< function_type > m_functions {}; /**< The function object to find the root for. */
//         MultiFunctionArray<return_type> m_functions {};
//
//         using RT = blaze::DynamicVector< return_type >;
//         RT m_guess; /**< The current root estimate. */
//
//         template< typename ARR >
//         explicit MultirootBase(const MultiFunctionArray<return_type>& functions, const ARR& guess)
//             : m_functions(functions)//,
//               //m_functions(functions.begin(), functions.end()),
//               //m_guess(guess.begin(), guess.end())
//         {
//             m_guess.resize(functions.size());
//
//             size_t index = 0;
//                 for (auto& g : guess) m_guess[index++] = g;
//         }
//
//         auto size() const { return m_functions.size(); }
//
//     public:
//
//         template< typename ARR >
// //        requires std::floating_point< typename nxx::deriv::impl::VectorTraits< ARR >::value_type > ||
// //                 std::floating_point< typename ARR::value_type >
//         auto evaluate(ARR values)
//         {
//             std::vector< return_type > vals(values.begin(), values.end());
//
//             auto tmp = m_functions(vals);
//             size_t index = 0;
//                 for (auto& f : tmp) values[index++] = f;
//
//             return values;
//         }
//
//         [[nodiscard]]
//         auto result() const
//         {
//             return std::vector< return_type >(m_guess.begin(), m_guess.end());
//         }
//     };
//
//     template< typename T >
//     //requires IsMultirootFunction< FunctionType >
//     class DMultiNewton final : public MultirootBase< DMultiNewton< T > >
//     {
//     public:
//         /*
//          * Public alias declarations.
//          */
//         //using function_type  = FunctionType;
//         using return_type    = T;//typename impl::FunctionTraits< FunctionType >::return_type;
//         //using container_type = typename impl::FunctionTraits< FunctionType >::container_type;
//         //using argument_type  = typename impl::FunctionTraits< FunctionType >::argument_type;
//
//     private:
//         /*
//          * Private alias declarations.
//          */
//         using BASE = MultirootBase< DMultiNewton< T > >;
//
//     public:
//         template< typename ARR >
//         explicit DMultiNewton(const MultiFunctionArray<T>& funcArray, const ARR& guess)
//             : BASE(funcArray, guess)
//         {}
//
//         explicit DMultiNewton(const MultiFunctionArray<T>& funcArray, const std::initializer_list<T> guess)
//             : BASE(funcArray, std::vector(guess.begin(), guess.end()))
//         {}
//
//         void iterate()
//         {
//             // Solve the linear system J * dx = -f(x) (solving for dx) and update the root estimate (x_new = x_old + dx).
//             BASE::m_guess += blaze::solve(nxx::deriv::jacobian(BASE::m_functions, BASE::m_guess),
//                                           -BASE::evaluate(BASE::m_guess));
//         }
//     };
//
//     // deduction guide:
//         //template <typename T>
//         //DMultiNewton(const DynamicFunctionArray<T>&) -> DMultiNewton<std::function<T(const std::vector<T>&)>>;
//
//     template< typename SOLVER >
//     //    requires requires(SOLVER solver, typename SOLVER::function_return_type guess) {
//     //                 // clang-format off
//     //                 { solver.evaluate(0.0) } -> std::floating_point;
//     //                 { solver.init(guess) };
//     //                 { solver.iterate() };
//     //                 // clang-format on
//     //             }
//     inline auto multisolve(SOLVER                                      solver,
//                            std::vector< typename SOLVER::return_type > guess,
//                            typename SOLVER::return_type                eps     = nxx::EPS,
//                            int                                         maxiter = nxx::MAXITER)
//     {
//         // using ET = std::runtime_error;
//         // using RT = tl::expected< typename SOLVER::return_type, ET >;
//         using RT = std::vector< typename SOLVER::return_type >;
//
//         //solver.init(guess);
//         RT result = solver.result();
//
//         auto calcEPS = [&]() {
//             auto                         fvals = solver.evaluate(result);
//             typename SOLVER::return_type eps   = 0.0;
//             for (auto& fval : fvals) eps += abs(fval);
//             return eps;
//         };
//
//         // Check for NaN or Inf.
//         //        if (!std::isfinite(solver.evaluate(result.value()))) {
//         //            result = tl::make_unexpected(ET("Invalid initial guess!", RootErrorType::NumericalError, result.value()));
//         //            return result;
//         //        }
//
//         // Begin iteration loop.
//         int iter = 1;
//         while (true) {
//             result = solver.result();
//
//             // Check for NaN or Inf
//             //            if (!std::isfinite(result.value())) {
//             //                result = tl::make_unexpected(ET("Non-finite result!", RootErrorType::NumericalError, result.value(), iter));
//             //                break;
//             //            }
//
//             // Check for convergence
//             if (calcEPS() < eps) {
//                 std::cout << "Converged after " << iter << " iterations." << std::endl;
//                 break;
//             }
//
//             if (iter >= maxiter) {
//                 std::cout << "Maximum number of iterations exceeded!" << std::endl;
//                 break;
//             }
//
//             // Check for max. iterations
//             //            if (iter >= maxiter) {
//             //                result = tl::make_unexpected(
//             //                    ET("Maximum number of iterations exceeded!", RootErrorType::MaxIterationsExceeded, result.value(), iter));
//             //                break;
//             //            }
//
//             // Perform one iteration
//             ++iter;
//             solver.iterate();
//         }
//
//         return result;
//     }
//
// }    // namespace nxx::multiroots
//
// #endif    // NUMERIXX_MULTIROOTS_IMPL_HPP
