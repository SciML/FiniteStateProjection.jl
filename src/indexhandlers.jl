"""
    singleindices(idxhandler::AbstractIndexHandler, arr)

Returns the indices of the finite state space represented by `arr`.

# Arguments
- `idxhandler`: The index handler that defines the state-space layout.
- `arr`: Either the state array or its tuple of dimensions.

# Returns
- An iterator of Cartesian indices. The default implementation returns
  `CartesianIndices(arr)`.

# Extension Rules
Custom index handlers must import and extend this function. Implement both the
state-array and dimension-tuple forms when supporting matrix conversions.
"""
singleindices(::AbstractIndexHandler, arr::AbstractArray) = CartesianIndices(arr)
singleindices(::AbstractIndexHandler, arr::Tuple) = CartesianIndices(arr)

"""
    pairedindices(idxhandler::AbstractIndexHandler, arr, shift::CartesianIndex)

Returns pairs `(I .- shift, I)` for state transitions that remain within `arr`.

# Arguments
- `idxhandler`: The index handler that defines the state-space layout.
- `arr`: Either the state array or its tuple of dimensions.
- `shift`: The reaction stoichiometry as a Cartesian index in Catalyst species order.

# Returns
- An iterator of source and destination Cartesian-index pairs.

# Extension Rules
Custom index handlers must import and extend this function. Implement both the
state-array and dimension-tuple forms when supporting matrix conversions.
"""
function pairedindices end

"""
    getsubstitutions(idxhandler::AbstractIndexHandler, rs::ReactionSystem; state_sym::Symbol)

Returns a dictionary `S_i => f_i(state_sym)`, where each `f_i` computes a species
abundance from the symbolic state variable.

# Arguments
- `idxhandler`: The index handler that defines the state-space layout.
- `rs`: Catalyst reaction system whose species require substitutions.

# Keyword Arguments
- `state_sym`: Symbol naming the generated RHS state variable.

# Returns
- A dictionary mapping each Catalyst species to a symbolic abundance expression.

# Extension Rules
Custom index handlers must import and extend this function.
"""
function getsubstitutions end

"""
    vec(idxhandler::AbstractIndexHandler, arr)

Converts the CME state array into the vector ordering used by sparse-matrix
representations.

# Arguments
- `idxhandler`: The index handler that defines the vector ordering.
- `arr`: State array to flatten.

# Returns
- A one-dimensional state vector compatible with [`LinearIndices`](@ref Base.LinearIndices).

# Extension Rules
Custom index handlers that support matrix conversions must import and extend this
`Base` function.

See also: [`LinearIndices`](@ref Base.LinearIndices)
"""
function vec end

"""
    LinearIndices(idxhandler::AbstractIndexHandler, arr)

Returns an object `lind` converting the Cartesian indices produced by
[`singleindices`](@ref) and [`pairedindices`](@ref) to the linear indices of
[`vec`](@ref Base.vec). The required invariant is

```julia
arr[idx_cart] == vec(idxhandler, arr)[idx_lin]
```

# Arguments
- `idxhandler`: The index handler that defines the vector ordering.
- `arr`: Either the state array or its tuple of dimensions.

# Returns
- An indexable object mapping Cartesian indices to vector positions.

# Extension Rules
Custom index handlers that support matrix conversions must import and extend this
`Base` function, preserving the documented indexing invariant.

See also: [`vec`](@ref Base.vec)
"""
function LinearIndices end

##

"""
    DefaultIndexHandler{N}()
    DefaultIndexHandler{N}(offset, perm)

Default index handler for an FSP system with `N` species.

# Fields
- `offset`: Julia index corresponding to a molecular count of zero.
- `perm`: Mapping from state-array dimensions to Catalyst species order.

It represents the state as an `N`-dimensional array, maps a molecule count of zero to
`offset`, and uses `perm` to map state-array dimensions to Catalyst's species order.

The zero-argument constructor uses Julia's one-based indexing and preserves the
species order. This representation is appropriate when every state in the truncated
array is reachable; reduce conserved species before construction when possible.

# Examples
```julia
julia > DefaultIndexHandler{2}()
DefaultIndexHandler{2}(1, (1, 2))
```
"""
struct DefaultIndexHandler{N} <: AbstractIndexHandler
    offset::Int
    perm::NTuple{N, Int}
end

DefaultIndexHandler{N}() where {N} = DefaultIndexHandler{N}(1, Tuple(1:N))

"""
    NaiveIndexHandler

Deprecated constructor forwarding to [`DefaultIndexHandler`](@ref).

# Arguments
- `args...`: Positional arguments accepted by `DefaultIndexHandler`.

# Keyword Arguments
- `kwargs...`: Keyword arguments accepted by `DefaultIndexHandler`.

Use `DefaultIndexHandler` directly instead.
"""
function NaiveIndexHandler(args...; kwargs...)
    Base.depwarn(
        "`NaiveIndexHandler` is deprecated, use `DefaultIndexHandler` instead.",
        :NaiveIndexHandler
    )
    return DefaultIndexHandler(args...; kwargs...)
end

Base.vec(::DefaultIndexHandler, arr) = vec(arr)
Base.LinearIndices(::DefaultIndexHandler, arr) = LinearIndices(arr)

function pairedindices(
        ih::DefaultIndexHandler{N}, arr::AbstractArray{T, N},
        shift::CartesianIndex{N}
    ) where {T, N}
    return pairedindices(ih, axes(arr), shift)
end

# `dims` is written as `Tuple{T, Vararg{T}}` rather than `NTuple{N, T}` so that the
# element type `T` is always bound: `NTuple{0, T} === Tuple{}` leaves `T` free, which
# trips Aqua's unbound-type-parameter check. The zero-dimensional case is handled by the
# explicit method below.
function pairedindices(
        ih::DefaultIndexHandler{N}, dims::Tuple{T, Vararg{T}},
        shift::CartesianIndex{N}
    ) where {N, T <: Number}
    return pairedindices(ih, Base.OneTo.(dims), shift)
end

# Handles the degenerate zero-dimensional case (no species), where neither `T`-constrained
# method below matches because `Tuple{}` has no element type.
function pairedindices(
        ::DefaultIndexHandler{0}, ::Tuple{}, ::CartesianIndex{0}
    )
    return zip(CartesianIndices(()), CartesianIndices(()))
end

# Important: the species in `shift` are ordered according to `Catalyst.species`!
function pairedindices(
        ih::DefaultIndexHandler{N}, dims::Tuple{T, Vararg{T}},
        shift::CartesianIndex{N}
    ) where {N, T <: AbstractVector}
    ranges = tuple(
        (
            UnitRange(
                    max(first(ax), first(ax) + shift[ih.perm[i]]),
                    min(last(ax), last(ax) + shift[ih.perm[i]])
                )
                for (i, ax) in enumerate(dims)
        )...
    )

    ranges_shifted = tuple((rng .- shift[ih.perm[i]] for (i, rng) in enumerate(ranges))...)

    return zip(CartesianIndices(ranges_shifted), CartesianIndices(ranges))
end

function pairedindices(
        ::DefaultIndexHandler, dims::Tuple,
        shift::CartesianIndex
    )
    return @error "Dimension of state space ($(length(dims))) does not match number of species ($(length(shift)))"
end

"""
    getsubstitutions(sys::FSPSystem{DefaultIndexHandler}; state_sym::Symbol)::Dict

Defines the abundance of species ``S_i`` to be `state_sym[i] - offset`.
"""
function getsubstitutions(ih::DefaultIndexHandler, rs::ReactionSystem; state_sym::Symbol)
    nspecs = numspecies(rs)
    state_array = (@variables ($state_sym)[1:nspecs])[1]
    state_sym_vec = [state_array[i] for i in 1:nspecs]

    species_orig = species(rs)
    species_perm = [species_orig[ih.perm[i]] for i in 1:nspecs]

    return Dict(symbol => state_sym_vec[i] - ih.offset for (i, symbol) in enumerate(species_perm))
end

#"""
#    PermutingIndexHandler(rs::ReactionSystem, order::AbstractVector)
#
#Constructs an index handler for the reaction system in which the species appear in the order
#defined by the vector `order`.
#"""
function PermutingIndexHandler(rs::ReactionSystem, order::AbstractVector{Symbol})
    system_species = species(rs)
    resolved_order = map(order) do sym
        index = findfirst(species -> getname(species) == sym, system_species)
        isnothing(index) && error("Cannot find species $sym in reaction system")
        system_species[index]
    end
    return PermutingIndexHandler(rs, resolved_order)
end

function PermutingIndexHandler(rs::ReactionSystem, order::AbstractVector)
    spec = species(rs)
    nspec = length(spec)

    if nspec != length(order)
        @error "Length of species vector ($(length(order))) does not match number of species ($nspec)"
    end

    perm = zeros(Int, nspec)
    count = zeros(Int, nspec)

    for i in 1:nspec
        idx = findfirst(s -> isequal(s, order[i]), spec)
        if isnothing(idx)
            @error "Cannot find species $(order[i]) in reaction system"
        end

        if count[idx] > 0
            @error "Species $(order[i]) specified twice in ordering"
        end

        count[idx] += 1
        perm[i] = idx
    end

    @assert count == ones(Int, nspec)

    return DefaultIndexHandler(1, Tuple(perm))
end
