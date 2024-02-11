***********************************
One-Dimensional Root-Finding
***********************************

This section delves into a range of classes and functions dedicated to the finding roots in one-dimensional functions. The featured algorithms are divided into two main types: bracketing algorithms, which operate without the need for function derivatives, and polishing algorithms, which require the computation of the function's derivative. In addition, the library offers search algorithms that can be used to find a bracket where a root exists when the initial guess is not near the actual root.

All functions are available through the :file:`Roots.hpp` header file.

Overview
========
Popular one-dimensional root-finding methods encompass the bisection method, Newton's method, Ridder's method, and the secant method. Each of these methods brings its unique advantages and potential drawbacks, and their applicability hinges on the nature of the problem at hand.

Bracketing methods such as the bisection method and Ridder's method are designed to locate a root within a specified interval. These methods are particularly useful when the function's derivative is unavailable or when the function is not differentiable. The bisection method, in particular, is known for its simplicity and reliability, as it guarantees convergence to a root as long as the function is continuous and changes sign within the interval. Ridder's method, on the other hand, can converge with fewer iterations, but each iteration is more computationally expensive.

The following bracketing methods are available in the library:

* **Bisection Method**: Known for its reliability and surefire convergence, the bisection method is a good choice when you need a method that is guaranteed to find a root. However, it can sometimes converge at a slower pace compared to other methods. (Wikipedia: `Bisection Method <https://en.wikipedia.org/wiki/Bisection_method>`_)

* **Ridder's Method**: Outpacing the bisection method, Ridder's method can converge more rapidly. However, each iteration is more computationally expensive, and the function must have a continuous second derivative. (Wikipedia: `Ridder's Method <https://en.wikipedia.org/wiki/Ridder%27s_method>`_)

* **Regula Falsi Method**: Also known as the false position method, the regula falsi method is a bracketing algorithm that uses linear interpolation to find the root of a function. It is similar to the bisection method, but it can converge faster in some cases. (Wikipedia: `Regula Falsi Method <https://en.wikipedia.org/wiki/Regula_falsi>`_)

Polishing methods, such as Newton's method and Steffensen's method, are designed to refine an initial guess to find a root. These methods are generally faster than bracketing methods, but they require the first derivative of the function to be available. Newton's method, in particular, is known for its speed and efficiency, but it may fail to converge or converge to a local minimum instead of a root.

The following polishing methods are available in the library:

* **Newton's Method**: A popular root-finding algorithm, Newton's method is known for its speed and efficiency. However, it requires the function to be differentiable and the derivative to be non-zero at the estimate. (Wikipedia: `Newton's Method <https://en.wikipedia.org/wiki/Newton%27s_method>`_)

* **Secant Method**: The secant method is a root-finding algorithm that does not require the derivative of the function. It is similar to the Newton-Raphson method, but it uses a finite difference approximation to the derivative. (Wikipedia: `Secant Method <https://en.wikipedia.org/wiki/Secant_method>`_)

* **Steffensen's Method**: Steffensen's method is an iterative root-finding algorithm that improves upon the simple fixed-point iteration by incorporating a form of Aitken's Δ² process. This method is particularly effective for functions where the derivative is difficult to compute or is not readily available. (Wikipedia: `Steffensen's Method <https://en.wikipedia.org/wiki/Steffensen%27s_method>`_)

Search methods, such as BracketSearchUp and BracketSearchDown, are designed to incrementally expand or subdivide the search bounds to find a bracket where a root exists. These methods are useful when the initial guess is not near the actual root.

The following search methods are available in the library:

* **BracketSearchUp**: A specialized search algorithm designed to incrementally expand the search bounds upwards (increasing values) to find a bracket where a root exists.

* **BracketSearchDown**: A specialized search algorithm designed to incrementally expand the search bounds downwards (decreasing values) to find a bracket where a root exists.

* **BracketExpandUp**: A specialized search algorithm designed to incrementally expand the upper bound upwards (increasing values) while keeping the lower bound fixed. This is useful for finding a bracket where a root exists when the initial guess is lower than the actual root.

* **BracketExpandDown**: A specialized search algorithm designed to incrementally expand the lower bound downwards (decreasing values) while keeping the upper bound fixed. This is useful for finding a bracket where a root exists when the initial guess is higher than the actual root.

* **BracketExpandOut**: A specialized search algorithm designed to incrementally expand both the lower and upper bounds symmetrically outwards. This is useful for finding a bracket where a root exists when the initial guess is not near the actual root.

* **BracketSubdivide**: A specialized search algorithm designed to subdivide the current search bounds into smaller segments in an attempt to find a bracket where a root exists.

Using these methods should be done through the :code:`fsolve`, :code:`fdfsolve`, and :code:`search` functions, which are designed to be easy to use and provide a consistent interface for all the solvers.

Important Considerations
========================

It's imperative to recognize that root-finding functions are designed to locate a single root at any given instance. In cases where multiple roots are present within the search range, the function will identify and return the initial root it encounters. Pinpointing which root will be found in a region with several roots is generally unpredictable. Notably, attempting to find a root in such areas usually doesn't trigger any error messages, despite the inherent challenges.


Quick Start
===========

The library provides the following functions for finding roots in one-dimensional functions:

* **fsolve**: A function for finding a root of a one-dimensional function using a bracketing method. It takes a bracketing method object and an initial bracket around the root to be found.

* **fdfsolve**: A function for finding a root of a one-dimensional function using a polishing method. It takes a polishing method object and an initial guess for the root.

* **search**: A function for finding a bracket where a root exists when the initial guess is not near the actual root. It takes a search method object and an initial guess for the root.

These functions return a proxy object that contains the result of the root-finding process. The result can be accessed using the :code:`result` member function, which returns the estimated root.

.. note::
    The proxy objects returned by the root-finding functions are not intended to be stored or passed around. They are designed to be used immediately to access the result of the root-finding process. Consequently, they cannot be copied or moved, and the member functions are only available for r-value references.

The following example shows how to find the root of the function :math:`f(x) = x^2 - 5` (defined as a lambda) using a bracketing method:

.. code-block:: c++

    using nxx::roots::fsolve;
    using nxx::roots::Bisection;
    auto result = fsolve<Bisection>([](double x){return x * x - 5;}), {0.0, 2.5}).result();

Similarly, the same function can be solved using Newton's method (with the derivative being :math:`f'(x) = 2x`)

.. code-block:: c++

    using nxx::roots::fdfsolve;
    using nxx::roots::Newton;
    auto result = fdfsolve<Newton>([](double x){return x * x - 5;},
                                   [](double x){return 2 * x;}), 1.25).result();

.. note::
    The derivative is not required to be provided, as it will be computed numerically using a finite difference approximation. However, this can be less accurate and less efficient than providing the derivative directly.

Finally, the following example shows how to find a bracket where a root exists when the initial guess is not near the actual root:

.. code-block:: c++

    using nxx::roots::search;
    using nxx::roots::BracketSearchUp;
    auto result = search<BracketSearchUp>([](double x){return x * x - 5;}, { 0.0, 1.0 }).result();

In addition to the arguments provided in the examples above, the root-finding functions also take two optional arguments for terminating the iterations, namely the relative error (epsilon) and the maximum number of iterations. Those arguments can be provided in any order; the only requirement is that epsilon must be a floating-point number and the maximum number of iterations must be an integer.

As can be seen from the previous examples, the root-finding solvers can be invoked in a single line of code, which makes them very easy to use.

Advanced Usage
==============

Although using the root-finding functions and algorithms can be as simple as shown in the previous examples, there are a number of advanced features and considerations that can be useful in more complex scenarios.

Customizing the Solvers
-----------------------

In many cases, the default settings of the solvers will be sufficient. However, in some cases it may be required to provide more fine-grained control over the solvers and the criteria for terminating the iterations. This can be done by creating a "stop token" object and passing it to the solver.

A stop token is an object that provides a way to stop the iterations of the solver based on a set of criteria. It can be any callable object, such as a lambda, that takes the solver's state as an argument and returns a boolean value indicating whether the iterations should continue or stop. The stop token can be used to check for convergence, reach a maximum number of iterations, or any other condition that should stop the iterations.

The following example shows how to create a custom stop token and use it with the :code:`fsolve` function:

.. code-block:: c++

    using nxx::roots::fsolve;
    using nxx::roots::Bisection;
    auto stop_token = [](const auto& solver) { return solver.iter() >= 100; };
    auto result = fsolve<Bisection>([](double x){return x * x - 5;}), {0.0, 2.5}, stop_token).result();

In this example, the stop token is a lambda that checks whether the number of iterations has reached 100. If the condition is met, the iterations will stop and the result will be returned.

As an alternative to providing the stop token as a function argument, it can also be provided as a second template argument, after the solver algorithm. This can be useful when the stop token is a complex type or when it needs to be reused with multiple solvers:

.. code-block:: c++

    using nxx::roots::fsolve;
    using nxx::roots::Bisection;
    auto stop_token = [](const auto& solver) { return solver.iter() >= 100; };
    using StopToken = decltype(stop_token);
    auto result = fsolve<Bisection, StopToken>([](double x){return x * x - 5;}), {0.0, 2.5}).result();

Both syntaxes are equivalent, so which one to use is a matter of taste.

However, determining when to terminate the iterations is not the only use of the stop token object; it can also be used to provide additional information about the iterations, such as the number of iterations that have been performed. This can be useful for logging or debugging purposes:

.. code-block:: c++

    using nxx::roots::fsolve;
    using nxx::roots::Bisection;
    auto stop_token = [](const auto& solver) {
        std::cout << "Iteration: " << solver.iter() << " | ";
        std::cout << "Guess:     " << solver.guess() << std::endl;
        return solver.iter() >= 100;
        return false; };
    auto result = fsolve<Bisection>([](double x){return x * x - 5;}), {0.0, 2.5}, stop_token).result();

This example shows how the stop token can be used to print the iteration number and the current guess at each iteration. This can be useful for understanding how the solver is converging and for diagnosing any issues that may arise. The output can be redirected to a file or a log for further analysis.

Finally, a stop token can also be used to deal with errors. For example, if reaching the maximum number of iterations is considered an error, the stop token can be used to throw an exception when the condition is met. Or if the current guess of the root is NaN or infinite, the stop token can be used to throw an exception as well:

.. code-block:: c++

    using nxx::roots::fsolve;
    using nxx::roots::Bisection;
    auto stop_token = [](const auto& solver) {
        if (solver.iter() >= 100) {
            throw std::runtime_error("Maximum number of iterations reached");
        }
        if (std::isnan(solver.guess()) || std::isinf(solver.guess())) {
            throw std::runtime_error("Invalid root");
        }
        return false; };
    auto result = fsolve<Bisection>([](double x){return x * x - 5;}), {0.0, 2.5}, stop_token).result();

As the previous examples show, the stop token can be used to provide fine-grained control over the solvers and the criteria for terminating the iterations. This can be useful in a wide range of scenarios, from simple cases where the default settings are sufficient to complex cases where more control is needed.

Customizing the Output
----------------------

In addition to customizing the solvers, it is also possible to customize the output of the root-finding functions. As mentioned earlier, the solver functions return a proxy object that contains the result of the root-finding process. This proxy object can be used to access the final result of the root-finding, using the :code:`result` member function.

In some cases, however, it may be useful to access additional information about the result, or to get the result in a different format. In the simplest form, a type that can be constructed from a floating point number, can be passed as a template argument to the :code:`result` function:

.. code-block:: c++

    using nxx::roots::fsolve;
    using nxx::roots::Bisection;
    auto result = fsolve<Bisection>([](double x){return x * x - 5;}), {0.0, 2.5}).result<double>();

In more complex scenarios, a callable object, such as a lambda, can be passed as an argument. This can be usefule in situations where some error condition need to be reported, but where throwing an exception is not appropriate. This could be the C++23 std::expected type, or one of the alternatives such as tl::expected from TartanLlama:

.. code-block:: c++

    using nxx::roots::fsolve;
    using nxx::roots::Bisection;
        auto outputter = [](const auto &data) -> tl::expected<decltype(data.guess), std::string> {
        auto [iter, guess, previous] = data;
        using expected = tl::expected<decltype(guess), std::string>;

        if (iter >= 5) return tl::make_unexpected(std::string("Too many iterations"));
        return expected(guess);
    };
    auto result = fsolve<Bisection>([](double x){return x * x - 5;}), {0.0, 2.5}).result(outputter);

The outputter can also be provided as a template argument. The functioanlity is the same,  so which one to use is a matter of taste:

.. code-block:: c++

    using nxx::roots::fsolve;
    using nxx::roots::Bisection;
    using Outputter = decltype(outputter);
    auto result = fsolve<Bisection>([](double x){return x * x - 5;}), {0.0, 2.5}).result<Outputter>();

The outputter can be used to provide additional information about the result, such as the number of iterations that were performed, or the reason for stopping the iterations. This can be useful for logging or debugging purposes, or for providing more context about the result of the root-finding process.

API Reference
=============

This section provides a detailed reference for the root-finding functions. The reference includes the function signature, the template parameters, the function arguments, and the return type.

For each function, the reference also includes a list of the available algorithms, along with a brief description of each algorithm and its use cases.

.. note::
    The algorithm classes are not intended to be used directly, but are provided for reference. When used in the root-finding functions, the algorithm template arguments will be deduced automatically, and should not be provided explicitly.

fsolve
------

The :code:`fsolve` functions are used to find a root of a one-dimensional function using a bracketing method. It takes a bracketing method object and an initial bracket around the root to be found. The concrete algorithm to be used is provided as a template argument, and the function to be solved is provided as a function argument. The function returns a proxy object that contains the result of the root-finding process.

Four overloads of the :code:`fsolve` function are provided, with different combinations of arguments, enabling the user to provide custom termination conditions and output formats.

.. doxygenfunction:: nxx::roots::fsolve(FN_T func, STRUCT_T bounds, ARGS... args)
    :project: Numerixx

.. doxygenfunction:: nxx::roots::fsolve(FN_T func, STRUCT_T bounds)
    :project: Numerixx

.. doxygenfunction:: nxx::roots::fsolve(FN_T func, const ARG_T (&bounds)[N], ARGS... args)
    :project: Numerixx

.. doxygenfunction:: nxx::roots::fsolve(FN_T func, const ARG_T (&bounds)[N])
    :project: Numerixx

Algorithms
^^^^^^^^^^

The following algorithms are available for use with the :code:`fsolve` functions. Note that the algorithms are not intended to be used directly, but are provided for reference. Also, when used in the :code:`fsolve` function, the algorithm template arguments will be deduced automatically, and should not be provided explicitly.

.. doxygenclass:: nxx::roots::Bisection

.. doxygenclass:: nxx::roots::Ridder

.. doxygenclass:: nxx::roots::RegulaFalsi

fdfsolve
--------

.. doxygenfunction:: nxx::roots::fdfsolve(FN_T func, DERIV_T derivative, GUESS_T guess, ARGS... args)
    :project: Numerixx

.. doxygenfunction:: nxx::roots::fdfsolve(FN_T func, GUESS_T guess, ARGS... args)
    :project: Numerixx

.. doxygenfunction:: nxx::roots::fdfsolve(FN_T func, DERIV_T derivative, GUESS_T guess)
    :project: Numerixx

.. doxygenfunction:: nxx::roots::fdfsolve(FN_T func, GUESS_T guess)
    :project: Numerixx


Algorithms
^^^^^^^^^^

The following algorithms are available for use with the :code:`fdfsolve` functions:

.. doxygenclass:: nxx::roots::Newton

.. doxygenclass:: nxx::roots::Secant

.. doxygenclass:: nxx::roots::Steffensen

search
------

.. doxygenfunction:: nxx::roots::search(FN_T function, STRUCT_T bounds, FACTOR_T ratio, ITER_T maxiter)
    :project: Numerixx

.. doxygenfunction:: nxx::roots::search(FN_T function, const ARG_T (&bounds)[N], FACTOR_T ratio, ITER_T maxiter)
    :project: Numerixx

Algorithms
^^^^^^^^^^

The following algorithms are available for use with the :code:`search` functions:

.. doxygenclass:: nxx::roots::BracketSearchUp

.. doxygenclass:: nxx::roots::BracketSearchDown

.. doxygenclass:: nxx::roots::BracketExpandUp

.. doxygenclass:: nxx::roots::BracketExpandDown

.. doxygenclass:: nxx::roots::BracketExpandOut

.. doxygenclass:: nxx::roots::BracketSubdivide



Design and Implementation Details
=================================

The solver functions documented above should be sufficient for the vast majority of use cases. However, in some cases it may be necessary to use the solvers directly, or to create custom solvers. This section provides an overview of the design and implementation details of the solvers, and how they can be used and customized.

Interface
---------

Both the bracketing and polishing algorithms are implemented using the overall architecture: a base class is defined, for keeping track of internal state between iterations, and the individual algorithms inherits from the base class. However, in order to avoid virtual functions, the architecture is implemented using static polymorphism through the Curiously Recurring Template Pattern (CRTP) [1]_.

The purpose of this design is that it avoids dynamic polymorphim and virtual functions, while ensuring that the solvers share a common interface. The downside of this approach is that all dependencies has to be resolved during compile-time, and it is not possible to dynamically plug in a different solver.

While the base classes (BracketingBase and PolishingBase, respectively) will not be called directly in client code, it is useful to know what the classes look like, as the individual solvers will inherit the interface of the base classes.

BracketingBase
^^^^^^^^^^^^^^

The BracketingBase class (located in the :code:`root::detail` namespace) look as follows:

.. doxygenclass:: nxx::roots::detail::BracketingBase
    :members: m_func, m_bounds, init, evaluate, result

PolishingBase
^^^^^^^^^^^^^

Similarly, the PolishingBase class (located in the :code:`root::detail` namespace) look as follows:

.. doxygenclass:: nxx::roots::detail::PolishingBase
    :members: m_func, m_deriv, m_guess, init, evaluate, derivative, result

Implementation
--------------

The implementation of the solvers is based on the algorithms described in the previous sections. The solvers are implemented as classes, with the individual algorithms being implemented through the function call operator (:code:`operator()`). This makes the implementation reasonably straightforward and easy to understand.


.. [1] Vandevoorde, D., Josuttis, N., Gregor, D. (2018). C++ Templates - The Complete Guide