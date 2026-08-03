# Index Handlers

The task of an index handler is to provide a mapping between the system state and the way it is stored in memory, usually as a multidimensional array. The standard approach is to represent the states of a system with ``s`` reactions as an ``s``-dimensional array and have the index ``(i_1, \ldots, i_s)`` correspond to the state ``(n_1 = i_1, \ldots, n_s = i_s)``. This is implemented by the class [`DefaultIndexHandler`](@ref), which accepts an offset argument to deal with Julia's 1-based indexing (so the Julia index $(1,\ldots,1)$ corresponds to the state with no molecules).

## Extending Index Handlers

To define a custom index handler, subtype [`AbstractIndexHandler`](@ref) and import
the extension functions before adding methods:

```julia
import Base: LinearIndices, vec
import FiniteStateProjection: build_rhs_header, getsubstitutions, pairedindices, singleindices
```

Implement `getsubstitutions`, plus array and dimension-tuple methods for
`singleindices` and `pairedindices`. Implement `LinearIndices` and `vec` as well when
the handler supports sparse-matrix conversions.
Override `build_rhs_header` only when generated RHS functions need parameters arranged
differently from the default vector form.

```@docs
AbstractIndexHandler
DefaultIndexHandler
NaiveIndexHandler
singleindices
pairedindices
getsubstitutions
build_rhs_header
Base.LinearIndices
```
