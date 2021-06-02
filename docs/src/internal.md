# Internal API

## [Index Handlers](@id index_handlers_internal)

### Index Handler Interface

User-defined index handlers should inherit from `AbstractIndexHandler` and implement the following methods:
- [`getsubstitutions`](@ref FiniteStateProjection.getsubstitutions)
- [`build_rhs_header`](@ref FiniteStateProjection.build_rhs_header)
- [`singleindices`](@ref FiniteStateProjection.singleindices)
- [`pairedindices`](@ref FiniteStateProjection.pairedindices)

```@docs
FiniteStateProjection.singleindices(::AbstractIndexHandler, arr)
FiniteStateProjection.pairedindices
FiniteStateProjection.getsubstitutions
FiniteStateProjection.build_rhs_header(::AbstractIndexHandler, ::FSPSystem)
```

### Built-in implementations
```@docs
elidedspecies(::AbstractMatrix{Int})
FiniteStateProjection.elisions
FiniteStateProjection.getsubstitutions(::NaiveIndexHandler, ::FSPSystem)
FiniteStateProjection.getsubstitutions(::DefaultIndexHandler, ::FSPSystem)
FiniteStateProjection.build_rhs_header(::DefaultIndexHandler, ::FSPSystem)
FiniteStateProjection.pairedindices(::DefaultIndexHandler, ::AbstractArray, ::CartesianIndex)
```

## Function Building
```@docs
FiniteStateProjection.build_rhs
FiniteStateProjection.unpackparams
FiniteStateProjection.build_ratefuncs
FiniteStateProjection.build_rhs_firstpass
FiniteStateProjection.build_rhs_secondpass
```
