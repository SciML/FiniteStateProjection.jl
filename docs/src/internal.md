# Internal API
```@meta
CurrentModule = FiniteStateProjection
```

## Reaction Systems
```@docs
conservationlaws(::AbstractMatrix{Int})
```

## [Index Handlers](@id index_handlers_internal)

### Index Handler Interface

User-defined index handlers should inherit from `AbstractIndexHandler` and implement the following methods:
- [`getsubstitutions`](@ref getsubstitutions)
- [`build_rhs_header`](@ref build_rhs_header)
- [`singleindices`](@ref singleindices)
- [`pairedindices`](@ref pairedindices)

For matrix conversions they should additionally implement:
- [`LinearIndices`](@ref Base.LinearIndices)
- [`vec`](@ref Base.vec)

```@docs
singleindices
pairedindices
getsubstitutions
build_rhs_header
LinearIndices
```

### Built-in implementations
```@docs
elidedspecies(::AbstractMatrix{Int})
elisions
getsubstitutions(::NaiveIndexHandler, ::ReactionSystem)
getsubstitutions(::ReducingIndexHandler, ::ReactionSystem)
build_rhs_header(::FSPSystem{ReducingIndexHandler})
pairedindices(::ReducingIndexHandler, ::AbstractArray, ::CartesianIndex)
```

## Function Building
```@docs
build_rhs
unpackparams
build_ratefuncs
build_rhs_firstpass
build_rhs_secondpass
```
