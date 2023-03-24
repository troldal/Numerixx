.. _polynomials:

***********
Polynomials
***********

This sections describes classes and functions for creating polynomial functions of arbitrary order. Declarations can be found in :file:`poly/Polynomial.hpp`.

Overview
========
A polynomial is a mathematical expression that consists of variables, coefficients, and exponents. It is a function that can be written in the form of an equation where the variables are raised to non-negative integer powers, and the coefficients are real numbers. In general, a polynomial can have one or more variables, and its degree is the highest power of the variable in the expression.

For example, the expression :math:`3x^2 - 2x + 1` is a polynomial in one variable, x, with degree 2. The coefficient of :math:`x^2` is 3, the coefficient of :math:`x` is -2, and the constant term is 1. The degree of this polynomial is 2, which means that the highest power of :math:`x` in the expression is 2.

Polynomials are a fundamental concept in mathematics and have many important applications in various fields, such as engineering, physics, computer science, and economics. They are used to model real-world problems, approximate complex functions, and solve equations.

One of the most important properties of polynomials is that they are closed under addition, subtraction, and multiplication. This means that if two polynomials are added, subtracted, or multiplied, the result is always another polynomial. For example, if we add the polynomials :math:`2x^3 + 3x^2 + 5x + 7` and :math:`x^3 - 2x^2 + 4x - 3`, we get the polynomial :math:`3x^3 + x^2 + 9x + 4`. This property makes polynomials easy to work with and manipulate algebraically.

Another important property of polynomials is that they have roots or zeros, which are the values of the variable that make the polynomial equal to zero. The number of roots of a polynomial is equal to its degree, and each root corresponds to a linear factor of the polynomial. For example, the polynomial :math:`3x^2 - 2x + 1` has two roots, which correspond to the linear factors :math:`(x - \frac{2 + \sqrt{2}}{6})` and :math:`(x - \frac{2 - \sqrt{2}}{6})`.

Polynomials also have a unique factorization property, which means that every polynomial can be factored into a product of irreducible polynomials. An irreducible polynomial is a polynomial that cannot be factored into two or more polynomials of lower degree. For example, the polynomial :math:`x^2 + 1` is irreducible over the real numbers, but it can be factored into the product of the two linear factors :math:`(x + i)` and :math:`(x - i)` over the complex numbers.

Quick Start
===========
