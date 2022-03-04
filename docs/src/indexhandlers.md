# Index Handlers

The task of an index handler is to provide a mapping between the system state and the way it is stored in memory, usually as a multidimensional array. The standard approach is to represent the states of a system with ``s`` reactions as an ``s``-dimensional array and have the index ``(i_1, \ldots, i_s)`` correspond to the state ``(n_1 = i_1, \ldots, n_s = i_s)``. This is implemented by the class [`DefaultIndexHandler`](@ref), which accepts an offset argument to deal with Julia's 1-based indexing (so the Julia index $(1,\ldots,1)$ corresponds to the state with no molecules). 

See the [internal API](@ref index_handlers_internal) on how to define your own `IndexHandler` type.

```@docs
DefaultIndexHandler
```
