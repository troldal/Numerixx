/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P    Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b          d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,      ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b    d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8  8P        Y8

    Copyright © 2022 Kenneth Troldal Balslev

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

#ifndef NUMERIXX_ERROR_HPP
#define NUMERIXX_ERROR_HPP

// ===== External Includes
#include <boost/stacktrace.hpp>

// ===== Standard Library Includes
#include <optional>
#include <source_location>
#include <sstream>
#include <stdexcept>
#include <utility>

using JSONString = std::string;

namespace nxx
{
    enum class NumerixxErrorType { General, Poly, Polyroots, Roots, MultiRoots, Deriv, Func };

    class NumerixxError : public std::runtime_error
    {
    public:
        explicit NumerixxError(const std::string&            str,
                               NumerixxErrorType             type  = NumerixxErrorType::General,
                               //                               JSONString                    data  = {},
                               const std::source_location&   loc   = std::source_location::current(),
                               boost::stacktrace::stacktrace trace = boost::stacktrace::stacktrace())
            : std::runtime_error { str },
              m_type { type },
              //              m_data { std::move(data) },
              m_location { loc },
              m_backtrace { std::move(trace) }
        {}

        virtual ~NumerixxError() = default;

        [[nodiscard]]
        virtual NumerixxErrorType type() const noexcept
        {
            return m_type;
        }

        //        [[nodiscard]]
        //        virtual JSONString data() const noexcept
        //        {
        //            return m_data;
        //        }

        [[nodiscard]]
        virtual const std::source_location& where() const noexcept
        {
            return m_location;
        }

        [[nodiscard]]
        virtual const boost::stacktrace::stacktrace& stack() const noexcept
        {
            return m_backtrace;
        }

        [[nodiscard]]
        virtual std::string log() const
        {
            std::stringstream logStream;
            logStream << "Error: " << what() << "\n\n";
            logStream << "Occurred in:\n\t";
            logStream << "File: " << m_location.file_name() << "\n\t";
            logStream << "Function: " << m_location.function_name() << "\n\t";
            logStream << "Line: " << m_location.line() << "\n\t";
            logStream << "Column: " << m_location.column() << "\n\n";
            //            logStream << "Details:\n" << m_data << "\n\n";
            logStream << "Stacktrace:\n" << boost::stacktrace::to_string(m_backtrace) << "\n";
            return logStream.str();
        }

    private:
        NumerixxErrorType m_type { NumerixxErrorType::General };
        //        const JSONString                    m_data {};
        std::source_location          m_location;
        boost::stacktrace::stacktrace m_backtrace;
    };

    template< typename T >
    class Error : public NumerixxError
    {
    public:
        explicit Error(const std::string& str, NumerixxErrorType type = NumerixxErrorType::General, T data = {})
            : NumerixxError(str, type),
              m_data { std::move(data) }
        {}

        [[nodiscard]]
        virtual T data() const noexcept
        {
            return m_data;
        }

        [[nodiscard]]
        std::string log() const override
        {
            std::stringstream logStream;
            logStream << NumerixxError::log() << "\n";
            logStream << "Details:\n" << m_data << "\n\n";
            return logStream.str();
        }

    private:
        T m_data;
    };

}    // namespace nxx

#endif    // NUMERIXX_ERROR_HPP
