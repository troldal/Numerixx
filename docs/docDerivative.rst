*********************
Numerical Derivatives
*********************

This sections describes classes and functions for finding derivatives of arbitrary one-dimensional functions. Functions for computing numerical derivatives can be found in :file:`calculus/Derivatives.hpp`.

Overview
========

Numerical derivatives are a numerical approximation of the derivative of a function. The derivative of a function represents the rate at which the function changes with respect to its independent variable. It is often useful in science and engineering to know the derivative of a function at a particular point, as it provides information about the function's behavior in the vicinity of that point.

There are several methods to compute numerical derivatives, but the most common ones are the forward difference, backward difference, and central difference methods.

The forward difference method involves computing the slope of a secant line between two points on the function, one slightly to the right of the point of interest, and one at the point of interest. The slope of this line approximates the derivative of the function at the point of interest, and can be computed as:

:math:`f'(x) ≈ \frac{f(x + h) - f(x)}{h}`

where h is a small positive number that determines the distance between the two points.

Similarly, the backward difference method involves computing the slope of a secant line between two points on the function, one slightly to the left of the point of interest, and one at the point of interest. The slope of this line approximates the derivative of the function at the point of interest, and can be computed as:

:math:`f'(x) ≈ \frac{f(x) - f(x - h)}{h}`

Finally, the central difference method involves computing the slope of a secant line between two points on the function, one slightly to the right of the point of interest, and one slightly to the left of the point of interest. The slope of this line approximates the derivative of the function at the point of interest, and can be computed as:

:math:`f'(x) ≈ \frac{f(x + h) - f(x - h)}{2h}`

where h is a small positive number that determines the distance between the two points.

These numerical methods are called "approximations" because they introduce some error compared to the exact derivative of the function. However, the error can be made arbitrarily small by decreasing the value of h, but this leads to increased numerical instability and rounding errors in the computation. Guidance for choosing the step size (h) is given by Press et al. [1]_.

Typically, central difference methods yield the best results because they evaluate the function symetrically around the point where the derivative is to be found. However, the forward and backward methods are useful in cases where the function can only be evaluated in one direction. An example of this is the function :math:`f(x) = \sqrt{x}`. Evaluation of the derivative at the point :math:`0.0` can only be done using  forward difference methods, as central and backward difference methods would require function evaluations at negative values, which would yield results in the complex domain for this particular function.

Quick Start
===========

The simplest way to compute the derivative of a function is to use one of the convenience functions in the :code:`nxx::deriv` namespace, i.e. the :code:`central`, :code:`forward`, or :code:`backward` functions. These functions will compute the 1st derivative of a function, using (unsurprisingly) central difference, forward difference, and backward difference methods. Internally, these functions use Richardson extrapolation to compute the derivatives.

All of the functions will return an :code:`tl::expected` (or :code:`std::expected` when C++23 is in widespread use) which will contain the the result or an :code:`DerivativeError` object in case of an error. For further details, please see the :ref:`error-handling` section.

The following example shows how to find the derivative of the function :math:`f(x) = ln(x) + 2x` at the value :math:`e`, using the :code:`central` function:

.. code-block:: c++

    #include <numbers>
    #include <calculus/Derivatives.hpp>

    using nxx::deriv::central;
    auto func = [](double x) { return std::log(x) + 2 * x; };
    auto result = *central(func, std::numbers::e);

In the above code example, the function is provided as a lambda. However, any callable object can be used, as long as it takes a single floating point argument, and returns a floating point value. The :code:`forward` and :code:`backward` functions works in an identical manner.

The API for using the :code:`central`, :code:`forward`, or :code:`backward` functions is given below.

.. doxygenfunction:: central
   :project: Numerixx

.. doxygenfunction:: forward
   :project: Numerixx

.. doxygenfunction:: backward
   :project: Numerixx

Advanced Use
============

Algorithms
----------

In case more fine-grained usage is required, the :code:`derivative<>` template function can be used. The :code:`derivative<>` function has the same signature as the convenience functions described above. However, in addition it also takes a single template argument with the type of the algorithm used to compute the derivative.

The following example shows how to use the function:

.. code-block:: c++

    #include <numbers>
    #include <calculus/Derivatives.hpp>

    using nxx::deriv::central;
    auto func = [](double x) { return std::log(x) + 2 * x; };
    auto result = *derivative<Order1CentralRichardson>(func, std::numbers::e);

In the example above, the :code:`derivative<>` function is passed the :code:`Order1CentralRichardson` type as the algorithm type. This is one of several algorithms that comes bundled with the Numerixx library. The following is a list of algorithms for computing the 1st derivatives:

* :code:`Order1CentralRichardson`
* :code:`Order1Central3Point`
* :code:`Order1Central5Point`
* :code:`Order1ForwardHighOrder`
* :code:`Order1Forward2Point`
* :code:`Order1Forward3Point`
* :code:`Order1BackwardHighOrder`
* :code:`Order1Backward2Point`
* :code:`Order1Backward3Point`

In addition, the following algorithms can be used for computation of the 2nd derivatives:

* :code:`Order2Central3Point`
* :code:`Order2Central5Point`
* :code:`Order2Forward3Point`
* :code:`Order2Forward4Point`
* :code:`Order2Backward3Point`
* :code:`Order2Backward4Point`

Descriptions of the algorithms listed above can be found in the references [2]_ [3]_ [4]_ [5]_.

In general, the algorithms with 3-5 function evaluations and the Richardson extrapolation algorithms provide the highest accuracy. One might think that more function evaluations would require more CPU time, but it turns out not to be significant. The following table shows benchmarks of the different algorithms, and it can be seen that the choice of algorithm does not significantly affect the computation time.

.. code-block::

    ----------------------------------------------------------------------
    Benchmark                            Time             CPU   Iterations
    ----------------------------------------------------------------------
    BM_Order1CentralRichardson       0.239 ns        0.239 ns   1000000000
    BM_Order1Central3Point           0.237 ns        0.237 ns   1000000000
    BM_Order1Central5Point           0.232 ns        0.232 ns   1000000000
    BM_Order1ForwardRichardson       0.231 ns        0.231 ns   1000000000
    BM_Order1Forward2Point           0.232 ns        0.232 ns   1000000000
    BM_Order1Forward3Point           0.229 ns        0.229 ns   1000000000
    BM_Order1BackwardRichardson      0.230 ns        0.230 ns   1000000000
    BM_Order1Backward2Point          0.231 ns        0.231 ns   1000000000
    BM_Order1Backward3Point          0.230 ns        0.230 ns   1000000000

The algorithms for computing 2nd derivatives have similar performance characteristics.

In addition to the bundled algorithms, it is also possible to provide a custom algorithm, as long as it conforms to the correct interface. Any callable object which takes as arguments a single-variable function, a floating point value and a floating point step size, and that returns a floating point result (essentially the same interface as the :code:`derivative<>` function) can be used. The following example shows how to supply a custom algorithm in the form of a lambda:

.. code-block:: c++

    #include <numbers>
    #include <calculus/Derivatives.hpp>

    using nxx::deriv::central;
    auto algo = [](auto function, double val, double stepsize) {
        return (function(val + stepsize) - function(val - stepsize)) / (2 * stepsize);
    };
    auto func = [](double x) { return std::log(x) + 2 * x; };
    auto result = *derivative<decltype(algo)>(func, std::numbers::e);

The API for using the :code:`diff<>` function is given below.

.. doxygenfunction:: diff(IsFunction auto function, ReturnType<decltype(function)> val, ReturnType<decltype(function)> stepsize)
   :project: Numerixx

Derivative Function Objects
---------------------------
The functions and algorithms described above can of course be wrapped in a function object representing the derivative
of the function in question. However, Numerixx provides a convenience function for creating function objects like this
automatically.

The :code:`derivativeOf<>` template function can be used to create a function object representing the derivative of
any input function, provided it has the correct interface. The function uses the :code:`Order1CentralRichardson`
algorithm by default, but other algorithms can be specified manually.

The following code illustrates how to use it:

.. code-block:: c++

    #include <numbers>
    #include <calculus/Derivatives.hpp>

    using namespace nxx::deriv;
    auto func = [](double x) { return std::log(x) + 2 * x; };
    auto d1func = derivativeOf(func);
    auto d2func = derivativeOf<Order2Central5Point>(func);

    double d1func_val = d1func(std::numbers::e);
    double d2func_val = d2func(std::numbers::e);

.. note::
    For objects of the :code:`Polynomial` class, a separate :code:`derivativeOf<>` function is defined in the
    :code:`nxx::poly` namespace. This function will create a function object that will compute the derivative
    analytically, rather than numerically. See the :ref:`polynomials` section for details.

The API for using the :code:`derivativeOf<>` function is given below.

.. doxygenfunction:: derivativeOf(IsFunction auto function, ReturnType< decltype(function) > stepsize)
   :project: Numerixx

.. [1] William H. Press et al. (2007). Numerical Recipes, 3rd Edition.
.. [2] Richard L. Burden et al. (2011). Numerical Analysis, 9th Edition.
.. [3] Steven C. Chapra et al. (2021). Numerical Methods for Engineers, 8th Edition.
.. [4] Sauer, T. (2012). Numerical Analysis, 2nd Edition.
.. [5] Ward Cheney et al. (2008). Numerical Mathematics and Computing, 6th Edition.

