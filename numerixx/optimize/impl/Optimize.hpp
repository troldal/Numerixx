//
// Created by kenne on 13/11/2023.
//

#ifndef NUMERIXX_OPTIMIZE_HPP
#define NUMERIXX_OPTIMIZE_HPP

#include <cmath>
#include <functional>

namespace nxx::optimize
{

    class GradientDescentOptimizer
    {
    public:
        enum class Mode { Minimize, Maximize };

        GradientDescentOptimizer(double learningRate = 0.01, Mode mode = Mode::Minimize)
            : learningRate(learningRate),
              mode(mode)
        {}

        template< typename Function >
        double operator()(Function func, double initialGuess) const
        {
            double x    = initialGuess;
            double grad = 0;
            do {
                grad = derivative(func, x);
                if (mode == Mode::Maximize) {
                    grad = -grad;
                }
                x -= learningRate * grad;
            }
            while (std::abs(grad) > tolerance);

            return x;
        }

    private:
        double learningRate;
        double tolerance = 1e-6;
        Mode   mode;

        template< typename Function >
        static double derivative(Function func, double x)
        {
            const double h = 1e-6;
            return (func(x + h) - func(x)) / h;
        }
    };

    template< typename ALGO, typename Function >
    double optimize(ALGO algorithm, Function func, double initialGuess)
    {
        return algorithm(func, initialGuess);
    }

    template< typename ALGO, typename Function >
    auto optimizationOf(ALGO algorithm, Function func)
    {
        return [=](double initialGuess) { return algorithm(func, initialGuess); };
    }

}    // namespace nxx::optimize

#endif    // NUMERIXX_OPTIMIZE_HPP
