# [Matrix Conversions](@id matrix_conversions)
```@meta
CurrentModule = FiniteStateProjection
```

FiniteStateProjection.jl provides functionality for building the right-hand side of the CME as a (sparse) matrix. This provides another way to solve the CME in time via
```julia
...

A = convert(SparseMatrixCSC, sys, idxhandler, dims, p)

prob = ODEProblem((du,u,p,t) -> mul!(du, p, u), u0, tt, A)

...
```

Note that the matrix `A` has to be rebuilt for every truncation size and every set of parameters,
a restriction not shared by the `ODEFunction` returned by FiniteStateProjection.jl.

```@docs
convert(::Type{SparseMatrixCSC}, ::FSPSystem, ::AbstractIndexHandler, 
::NTuple, ps; combinatoric_ratelaw::Bool=true)
vec
```
