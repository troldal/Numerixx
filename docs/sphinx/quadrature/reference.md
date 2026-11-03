# nxx::quadrature

Integrals of a real function of one variable, from `numerixx/quadrature.hpp`. The [overview](index.md) introduces
the methods and when to use each; [Set tolerances](tolerances.md) explains the default request, which here asks for
half the type's digits, unlike `nxx::roots`; [Working with Numerixx](../basics.md#what-a-result-holds) says how to read
a result, and why `Criteria` and `error` must be qualified when both modules are in use. A method of your own is
checked against the concepts below ([A method of your own](method-of-your-own.md)).

```{doxygenconcept} nxx::quadrature::fixed_rule
```

```{doxygenconcept} nxx::quadrature::iterative_method
```

```{doxygenconcept} nxx::quadrature::estimating_method
```

```{doxygennamespace} nxx::quadrature
:members:
:undoc-members:
```

```{doxygenstruct} std::formatter< nxx::quadrature::error >
:members:
```

```{doxygenstruct} std::formatter< nxx::quadrature::Failure< T > >
:members:
```
