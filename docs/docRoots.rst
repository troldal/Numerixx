****************************
One Dimensional Root-Finding
****************************

This chapter describes classes and functions for finding roots of one-dimensional functions. The algorithms can either be bracketing algorithms (not requiring function derivatives), or polishing algorithms (require computation of the function derivative).

The header files :file:`RootBracketing.hpp` and :file:`RootPolishing.hpp` contains the implementations for the root finding algorithms as well as common base classes to keep track of internal state during iteration.

Overview
========

One-dimensional root-finding methods are numerical algorithms that are used to locate the roots or zeros of a function of one variable. These methods are important in a wide range of applications in science, engineering, and mathematics, where it is often necessary to find the value of the independent variable that makes a function equal to zero.

The most commonly used one-dimensional root-finding methods include the bisection method, Newton's method, Ridders' method, and the secant method. Each method has its own strengths and weaknesses, and the choice of method depends on the specific problem being solved.

The bisection method is a robust and reliable method that is guaranteed to converge, but may be slow to converge. Newton's method is a faster method that requires the function to be differentiable, but may not converge or converge to a local minimum. Ridders' method and the secant method are faster than the bisection method and can converge faster, but may require the function to have a continuous second derivative and may not converge in some cases. As with any numerical method, it is important to understand the assumptions and limitations of each method and to carefully evaluate the accuracy and convergence properties of the solutions obtained.

Caveats
=======

It should be noted that root finding functions have the capability to search for only one root at a time. In situations where there exist multiple roots within the search range, the function will return the first root it finds, but it's difficult to determine which of the roots it will be. Typically, attempting to find a root in an area containing multiple roots will not generate any error message.

Quick Start
===========

The easiest way to use the one-dimensional root finding solvers is to create an instance of one of the solver classes (e.g. :code:`Bisection` or :code:`Newton`) and plug it in to the :code:`fsolve` or :code:`fdfsolve` function. The solver objects can be instanitated with any callable object, such a lambda or a normal C++ function.

The following example shows how to find the root of the function :math:`f(x) = x^2 - 5` (defined as a lambda) using a bracketing method:

.. code-block:: c++

    using nxx::roots::fsolve;
    using nxx::roots::Bisection;
    auto result = fsolve(Bisection([](double x){return x * x - 5;}), {0.0, 2.5});

Similarly, the same function can be solved using Newton's method (with the derivative being :math:`f'(x) = 2x`)

.. code-block:: c++

    using nxx::roots::fdfsolve;
    using nxx::roots::Newton;
    auto result = fdfsolve(Newton([](double x){return x * x - 5;}, [](double x){return 2 * x;}), 1.25);

As can be seen from the previous examples, the root-finding solvers can be invoked in a single line of code, which makes it very easy to use.

In addition to the solver object, the :code:`fsolve` function also needs the initial brackets around the root to be found. This must be provided as a :code:`std::pair` object which can be provided as an initializer list, as shown in the first example. Similarly, the :code:`fdfsolve` function must be given an initial guess which should be relatively close to the actual root. Both the :code:`fsolve` and the :code:`fdfsolve` functions also take two optional arguments for terminating the iterations, namely the relative error (epsilon) and the maximum number of iterations.

The documentation for the :code:`fsolve` and :code:`fdfsolve` functions are shown below:

.. doxygenfunction:: fsolve
   :project: Numerixx

.. doxygenfunction:: fdfsolve
   :project: Numerixx

Root-Finding Methods
====================

A number of algorithms are implemented for both the bracketing methods and polishing methods. They are all template classes, but the template parameters are automatically deduced by providing the constructor with a callable object as an argument for the function to find the root for. For the polishing methods, a callable object for the derivative is also required.

The methods are described in the following.

Bracketing methods
------------------

Bracketing methods are one-dimensional root-finding algorithms that work by first identifying an interval, or bracket, containing the root. They start with two points in the interval and iteratively narrow it down until the root is isolated within a desired tolerance. These methods are guaranteed to converge to a root as long as the function is continuous and changes sign within the interval.

Bisection
^^^^^^^^^

The bisection method is a simple algorithm for finding a root of a one-dimensional function. It works by repeatedly dividing an interval in half and determining which half contains a root, until the root is found to within a desired level of accuracy. The method is guaranteed to converge to a root as long as the function is continuous and changes sign on the interval.

This method is considered to be inefficient as it typically requires more iterations compared to more advanced methods (e.g. Ridder's method). However, due to the simple nature of the method, each iteration takes much less computational effort and the overall performance is therefore often quite good.

.. doxygenclass:: nxx::roots::Bisection
   :members:

Ridders' method
^^^^^^^^^^^^^^^

Ridders' method is a root-finding algorithm for one-dimensional functions that uses an iterative process to refine the location of the root. It works by fitting a parabola through three points and using the vertex of the parabola as the next estimate for the root. This estimate is then refined by applying a scaling factor to the distance between the estimates to reduce the error. The method is efficient and can converge faster than the bisection method, but it requires the function to be twice differentiable and have a continuous second derivative.

.. doxygenclass:: nxx::roots::Ridders
   :members:

Polishing Methods
-----------------

Polishing methods typically start with an initial guess and iteratively improve it by evaluating the function and its derivative at the guess point. By using this information, the algorithms can better estimate the location of the root and converge more quickly than non-derivative methods. These methods can be very effective, but they may be sensitive to the choice of initial guess and can converge slowly or not at all in some cases.

Newton's method
^^^^^^^^^^^^^^^

Newton's method is a popular root-finding algorithm for one-dimensional functions. It works by making a linear approximation of the function at the current estimate of the root and finding the point where this approximation crosses the x-axis. This point becomes the next estimate for the root, and the process is repeated until convergence is achieved. Newton's method is generally faster than the bisection and Ridders' methods, but it requires the function to be differentiable and the derivative to be non-zero at the estimate. Additionally, the method may fail to converge or converge to a local minimum instead of a root.

.. doxygenclass:: nxx::roots::Newton
   :members:

Discrete Newton's method
^^^^^^^^^^^^^^^^^^^^^^^^

Discrete Newton's method is a variant of Newton's method that is used for finding roots of discrete functions or numerical data. Instead of computing the derivative of the function at each estimate, the discrete derivative is computed using the available data points. This method approximates the second derivative using the difference between the first derivatives at adjacent points, and then iteratively refines the estimate of the root using a similar approach as Newton's method. Discrete Newton's method can be an effective way to find roots of numerical data, but it may be less stable than Newton's method when used on analytic functions.

.. doxygenclass:: nxx::roots::DNewton
   :members:

Design and Implementation Details
=================================

Both the bracketing and polishing algorithms are implemented using the overall architecture: a base class is defined, for keeping track of internal state between iterations, and the individual algorithms inherits from the base class. However, in order to avoid virtual functions, the architecture is implemented using static polymorphism through the Curiously Recurring Template Pattern [1]_.

The purpose of this design is that it avoids dynamic polymorphim and virtual functions, while ensuring that the solvers share a common interface. The downside of this approach is that all dependencies has to be resolved during compile-time, and it is not possible to dynamically plug in a different solver.

While the base classes (BracketingBase and PolishingBase, respectively) will not be called directly in client code, it is useful to know what the classes look like, as the individual solvers will inherit the interface of the base classes.

BracketingBase
--------------

The BracketingBase class (located in the :file:`root::impl` namespace) look as follows:

.. doxygenclass:: nxx::roots::impl::BracketingBase
    :members: m_func, m_bounds, init, evaluate, result

PolishingBase
--------------

Similarly, the PolishingBase class (located in the :file:`root::impl` namespace) look as follows:

.. doxygenclass:: nxx::roots::impl::PolishingBase
    :members: m_func, m_deriv, m_guess, init, evaluate, derivative, result

Creating Concrete Root-Finding Solvers
--------------------------------------

Blah



.. [1] Vandevoorde, D., Josuttis, N., Gregor, D. (2018). C++ Templates - The Complete Guide