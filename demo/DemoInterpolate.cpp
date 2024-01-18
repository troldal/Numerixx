#include <Interp.hpp>

#include <cmath>

#include <iomanip>
#include <iostream>
#include <numbers>

int main()
{
    using namespace nxx::interp;

    //        std::vector<std::pair<double, double>> points = {{0, 0}, {1, 1}, {2, 4}, {3, 9}};
    //
    //        // Direct interpolation
    //        double val1 = interpolate(CubicSplineInterp{}, points, 1.5);
    //        std::cout << "Interpolated value at 1.5: " << val1 << std::endl;
    //
    //        // Using lambda function for interpolation
    //        auto interpFunc = interpolationOf(CubicSplineInterp{}, points);
    //        double val2 = interpFunc(2.5);
    //        std::cout << "Interpolated value at 2.5: " << val2 << std::endl;
    //
    //        return 0;

    std::vector< std::pair< double, double > > points = { { 0, 0 }, { 1, 2 }, { 2, 3 }, { 3, 2 }, { 4, 0 } };
    CubicSplineInterp                          splineInterp;
    SplineCoefficients                         coeffs = splineInterp(points);

    double xToEvaluate       = 3.5;
    double interpolatedValue = evaluateSpline(points, coeffs.a, coeffs.b, coeffs.c, coeffs.d, xToEvaluate);

    std::cout << "Interpolated value at " << xToEvaluate << " is " << interpolatedValue << std::endl;

    return 0;
}