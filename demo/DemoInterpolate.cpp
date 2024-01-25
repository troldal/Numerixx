#include <Interp.hpp>
#include <algorithm>

#include "_external.hpp"

int main()
{
    using namespace nxx::interp;
    using namespace sciplot;

    std::vector x_points { 0.0,         1.420735492, 2.454648713, 3.070560004, 3.621598752,
                           4.520537863, 5.860292251, 7.328493299, 8.494679123, 9.206059243 };
    std::vector y_points { 1.0f,         1.540302306f, 1.346356379f, 2.088869738f, 3.042340520f,
                           5.991202812f, 5.872036310f, 7.300592544f, 8.391857230f, 9.776685982f };

    auto linear   = interpolationOf< Linear >(x_points, y_points);
    auto lagrange = interpolationOf< Lagrange >(x_points, y_points);
    auto steffen  = interpolationOf< Steffen >(x_points, y_points);
    auto spline   = interpolationOf< Spline >(x_points, y_points);

    Vec x = linspace(0.0, 10.0, 200);

    std::vector< double > y_lin;
    std::transform(begin(x), end(x), std::back_inserter(y_lin), [&](auto val) { return linear.extrapolate(val); });

    std::vector< double > y_lag;
    std::transform(begin(x), end(x), std::back_inserter(y_lag), [&](auto val) { return lagrange.extrapolate(val); });

    std::vector< double > y_ste;
    std::transform(begin(x), end(x), std::back_inserter(y_ste), [&](auto val) { return steffen.extrapolate(val); });

    std::vector< double > y_spl;
    std::transform(begin(x), end(x), std::back_inserter(y_spl), [&](auto val) { return spline.extrapolate(val); });

    Plot2D plot1;
    plot1.palette("paired");
    plot1.drawPoints(x_points, y_points).label("Points").lineWidth(4);
    plot1.drawCurve(x, y_lin).label("Linear").lineWidth(4);
    plot1.drawCurve(x, y_lag).label("Lagrange").lineWidth(4);
    plot1.drawCurve(x, y_ste).label("Steffen").lineWidth(4);
    plot1.drawCurve(x, y_spl).label("Spline").lineWidth(4);

    Figure figure = { { plot1 } };
    Canvas canvas = { { figure } };

    canvas.defaultPalette("set1");
    canvas.size(1600, 1200);
    canvas.show();
    canvas.save("plot.svg");

    return 0;
}