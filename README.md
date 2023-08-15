# Numerixx: A Header-Only C++ Library for Numerical Computations

Numerixx is a header-only C++ library that prioritizes ease of use and encapsulation, providing a collection of numerical computation algorithms for common tasks. While performance is important, the primary focus is on offering user-friendly interfaces and easy extensibility. The library covers algorithms for numerical derivatives, polynomials, one-dimensional root finding, root searching, and multi-dimensional root finding. Additional algorithms will be added in the future.

## Features

- User-friendly interfaces: Easily integrate numerical computations into your projects with intuitive function calls and classes.
- Polynomial operations: Evaluate, differentiate, integrate, and manipulate polynomials with ease.
- One-dimensional root finding: Quickly find roots of functions using both bracketing and polishing methods.
- Root searching: Identify brackets where roots may be found for further root-finding iterations.
- Multidimensional root finding: Solve systems of equations to find multidimensional roots.

## Getting Started

1. Clone the Numerixx repository:

   ```bash
   git clone https://github.com/yourusername/numerixx.git
    ```
2. Add the `include` directory to your project's include path.
3. Include the desired header files in your source code.

## Example Usage

```cpp
#include <iostream>
#include "numerixx/numerixx.hpp"

int main() {
// Calculate the derivative of a function
double result = numerixx::derivative([](double x) { return x * x; }, 2.0);

    std::cout << "Derivative at x = 2: " << result << std::endl;

    // Find a root of a function
    double root = numerixx::findRoot([](double x) { return x * x - 4.0; }, 0.0, 3.0);

    std::cout << "Root: " << root << std::endl;

    // ... More examples ...
    
    return 0;
}
```

## Documentation
For detailed usage instructions, function references, and examples, please refer to the Numerixx Documentation.

## Contributing
Contributions to Numerixx are welcome! If you find any issues or want to add new features, please submit a pull request or open an issue in the repository.

## License
This project is licensed under the MIT License.

## Acknowledgments
Numerixx is inspired by the need for efficient and reliable numerical computation tools in C++. It draws inspiration from various numerical analysis textbooks and open-source libraries.
