/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P    Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b          d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,      ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b    d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8  8P        Y8

    Copyright © 2023 Kenneth Troldal Balslev

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the “Software”), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is furnished
    to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef NUMERIXX_MULTIFUNCTIONARRAY_HPP
#define NUMERIXX_MULTIFUNCTIONARRAY_HPP

#include <concepts>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <Concepts.hpp>

namespace nxx::multiroots
{
    /**
     * @class DynamicFunctionArray
     * @brief A class for managing a dynamic array of functions.
     *
     * DynamicFunctionArray is designed to store a collection of std::function objects that take a std::vector<T>
     * as input and return a value of type T. It provides mechanisms to add new functions and apply all stored
     * functions to a set of inputs.
     *
     * @tparam T The type of elements processed by the functions.
     */
    template< IsFloat T >
    class DynamicFunctionArray
    {
    public:
        // Type alias for the function signature and iterator
        using FUNC_T         = std::function< T(const std::span< T >&) >;
        using const_iterator = typename std::vector< FUNC_T >::const_iterator;

        /**
         * @brief Default constructor.
         */
        DynamicFunctionArray() = default;

        /**
         * @brief Initializer list constructor.
         * @param initList An initializer list of functions.
         */
//        DynamicFunctionArray(std::initializer_list< FUNC_T > initList)
//            : functions(initList)
//        {}

        /**
         * @brief Constructor from a container of functions.
         * @param container A container of function objects.
         */
        template< typename Container >
        requires IsSpanInvocable< typename Container::value_type >
        explicit DynamicFunctionArray(const Container& container)
        {
            functions.assign(container.begin(), container.end());
        }

        /**
         * @brief Constructor for range-based initialization.
         * @param first The beginning of the range.
         * @param last The end of the range.
         */
        template< typename Iter >
        requires IsFloatInvocable< typename std::iterator_traits<Iter>::value_type >
        DynamicFunctionArray(Iter first, Iter last)
        {
            functions.assign(first, last);
        }

        /**
         * @brief Adds a new function to the array.
         * @param func The function to add.
         */
        void addFunction(const FUNC_T& func) { functions.push_back(func); }

        /**
         * @brief Function call operator for arbitrary containers.
         * @param input The input container.
         * @return A std::vector<T> containing the results of applying each function to the input.
         */
        //template< ContainerOfT< T > Container >
        template< typename Container >
        std::vector< T > operator()(const Container& input) const
        {
            return evaluate(input.begin(), input.end());
        }

        /**
         * @brief Function call operator for initializer lists.
         * @param input The initializer list input.
         * @return A std::vector<T> containing the results of applying each function to the input.
         */
        std::vector< T > operator()(std::initializer_list< T > input) const { return evaluate(input.begin(), input.end()); }

        /**
         * @brief Accesses an individual function object.
         * @param index The index of the function in the array.
         * @return A const reference to the function.
         * @throws std::out_of_range If the index is out of range.
         */
        const FUNC_T& operator[](size_t index) const
        {
            if (index >= functions.size()) {
                throw std::out_of_range("Index out of range");
            }
            return functions[index];
        }

        /**
         * @brief Returns a const iterator to the beginning of the function array.
         * @return A const iterator.
         */
        const_iterator begin() const { return functions.cbegin(); }
        const_iterator cbegin() const { return functions.cbegin(); }

        /**
         * @brief Returns a const iterator to the end of the function array.
         * @return A const iterator.
         */
        const_iterator end() const { return functions.cend(); }
        const_iterator cend() const { return functions.cend(); }

        /**
         * @brief Returns the number of functions in the array.
         * @return The number of functions.
                */
        auto size() const { return functions.size(); }

    private:
        std::vector< FUNC_T > functions;    ///< Internal storage for function objects.

        /**
         * @brief Helper function to evaluate all functions using iterators.
         * @param begin The beginning iterator of the input range.
         * @param end The end iterator of the input range.
         * @return A std::vector<T> containing the results of each function.
         */
        template< typename Iter >
        std::vector< T > evaluate(Iter begin, Iter end) const
        {
            std::vector< T > results;
            results.reserve(functions.size());
            std::vector< T > args(begin, end);
            for (const auto& func : functions) {
                results.push_back(func(args));
            }
            return results;
        }
    };








    // Concepts for checking if a type is a std::vector
    template< typename T, typename Container >
    concept IsVectorOfT = requires(Container a) {
                              {
                                  *a.begin()
                              } -> std::convertible_to< T >;
                          };

    template< typename T, typename... Functions >
    class StaticFunctionArray
    {
    public:
        // Constructor
        explicit StaticFunctionArray(Functions... funcs)
            : functions(funcs...)
        {}

        // Function call operator for arbitrary container
        template< IsVectorOfT< T > Container >
        std::vector< T > operator()(const Container& input) const
        {
            std::vector< T > results;
            results.reserve(sizeof...(Functions));
            applyFunctions< 0 >(input, results);
            return results;
        }

        // Function call operator for initializer list
        std::vector< T > operator()(std::initializer_list< T > input) const { return this->operator()(std::vector< T >(input)); }

    private:
        std::tuple< Functions... > functions;

        // Helper to recursively apply functions from the tuple
        template< std::size_t I, IsVectorOfT< T > Container >
        void applyFunctions(const Container& input, std::vector< T >& results) const
        {
            if constexpr (I < sizeof...(Functions)) {
                results.push_back(std::get< I >(functions)(input));
                applyFunctions< I + 1 >(input, results);
            }
        }
    };

}    // namespace nxx::multiroots

#endif    // NUMERIXX_MULTIFUNCTIONARRAY_HPP
