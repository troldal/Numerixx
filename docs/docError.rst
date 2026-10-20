.. _error-handling:

***************************
A Word About Error Handling
***************************

Error handling in Numerixx is quite simple: there is none!

That may sound strange, as there may be many situations where you would like to handle errors. However, what conditions are considered errors is often a matter of opinion. For example, if you are trying to calculate the square root of a negative number, you may consider that an error. However, in some contexts, it may be perfectly acceptable to return a complex number in that case. In other cases, you may want to return a NaN or raise an exception.

For that reason, Numerixx approaches matters in a different manner. Out of the box, there is no error handling of any kind. However, the library has been designed to be customizable, so you can easily add whatever error handling you want.

Each section of the documentation will provide specific details about how a given function or class can be customized, for example to handle errors. Here follows an example from the section about :doc:`one-dimensional root finding <docRoots>`:

Finding the root of a function using the bisection method can be done like so:

.. code-block:: c++

   auto root = fsolve<Bisection>(func, { 0.0, 2.5 }).result();

This will use default settings for the accuracy and max. number of iterations. But what should happen if the max. number of iterations are exceeded? In some situations, this may be considered an exceptional conditions, and an exception should be thrown; in other cases, that may be overkill.

To handle this, you can provide a custom "stop token". This is a functor (such as a lambda), that can be provided as an additional argument to the `fsolve` function. This functor will be called at each iteration, and if it returns true, the iteration will stop. Here is an example:

.. code-block:: c++

   auto root = fsolve<Bisection>(func, { 0.0, 2.5 }, [](auto& state) {
       if (/*number of iterations exceeded*/)
           throw std::runtime_error("Too many iterations");
        If(/*required accuracy has been reached*/)
            return true;
       return false;
   }).result();

In this example, the lambda will throw an exception if the number of iterations exceeds 100. If you want to handle the error in a different way, you can do so by modifying the lambda accordingly.

As an alternative to providing the token as an argument, it can also be provided as an additional template argument to the `fsolve` function. This can be useful if you want to use the same stop token for many different calls to fsolve. Here is an example:

.. code-block:: c++

   auto root = fsolve<Bisection, MyStopToken>(func, { 0.0, 2.5 }).result();

In this case, `MyStopToken` is a class that implements the stop token interface. The interface is quite simple, and consists of a single method, which is called at each iteration. The method should return true if the iteration should stop, and false otherwise. The method is passed a reference to the state object, which contains information about the current iteration. The state object is a simple struct, and you can access its members directly.

If throwing an exception is not an option, you can also customize the return type of the fsolve function. By default, it just returns the latest guess of the root, regardless of the reason for stopping. However, you can customize the return type to include additional information, such as the reason for stopping. One option is to use `std::expected` from C++23 (or alternative implementations, such as `tl::expected` from TartanLlama), which allows you to return either a value or an error. Here is an example:

.. code-block:: c++

    auto outputter = [](const auto &data) -> tl::expected<decltype(data.guess), std::string> {
        auto [iter, guess, previous] = data;
        using expected = tl::expected<decltype(guess), std::string>;

        if (iter >= 5)
            return tl::make_unexpected(std::string("Too many iterations"));
        return expected(guess);
    };

    auto root = fsolve<Bisection>(func, { 0.0, 2.5 }).result(outputter);

The `outputter` is a functor that takes the state object as input, and returns a `tl::expected` object, that contains either the value or an error. In this case, the outputter will return an error message if the number of iterations exceeds 5. If you want to handle the error in a different way, you can modify the outputter accordingly.

Similar to the stop token, the outputter can also be provided as an additional template argument to the `result` function. This can be useful if you want to use the same outputter for many different calls to `fsolve`. Here is an example:

.. code-block:: c++

    auto root = fsolve<Bisection>(func, { 0.0, 2.5 }).result<MyOutputter>();

In this case, `MyOutputter` is a class that implements the outputter interface.

In summary, error handling in Numerixx is quite simple, but also quite flexible. You can easily customize the error handling to suit your needs, and you can do so in a way that is consistent across the entire library.