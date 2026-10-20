************
Introduction
************

The Numerixx library is a collection of functions and classes for numerical computing. Numerixx was developed in C++, with the purpose of providing implementations of the most common numerical algorithms. It is by no means the purpose to be as feature rich as other numerical libraries, such as the GNU Scientific Library (GSL), nor is it the purpose to be as super-optimized for performance as other numerical libraries. Rather, the intent is to provide modern C++ implementations of common numerical algorithms, and make them available as a header-only library, so that it can be integrated easily into client code.

Emphasis have been put on the architecture of the library, so that it is easy to maintain and easy to extend. In many cases, it will be possible for users to use their own custom algorithms, and plug them in to the Numerixx library. Each section will provide guidance for how to extend the library with custom algorithms.

While performance is always a concern, the primary focus of the library is to provide correct and reliable numerical algorithms, with an intuitive interface. Performance is a secondary concern, and the library is not intended to be used for high-performance computing. So if you intend to implement the next-generation weather forecast model, you should probably look elsewhere.

License and Warranty
====================

Numerixx is licensed under the MIT Software License. The MIT Software License is a permissive, open-source software license that allows users to freely use, modify, distribute, and sublicense the licensed software without any warranty. The library can also be used in commercial closed-source software; the license only requires that the license and copyright notice be included in any copies or modifications of the software.

Please note, however, that the Numerixx library has no warranty, and is provided “as is”. It is the responsibility of the users to validate the behavior and accuracy of the library.

MIT License

Copyright (c) 2024 Kenneth Troldal Balslev and KinetiQ.dev

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


