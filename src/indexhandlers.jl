"""
    singleindices(idxhandler::AbstractIndexHandler, arr)

Returns all indices `I` in `arr`. Defaults to CartesianIndices, but can
be overloaded for arbitrary index handlers.
"""
singleindices(::AbstractIndexHandler, arr::AbstractArray) = CartesianIndices(arr)
singleindices(::AbstractIndexHandler, arr::Tuple) = CartesianIndices(arr)

"""
    pairedindices(idxhandler::AbstractIndexHandler, arr, shift::CartesianIndex)

Returns all pairs of indices `(I .- shift, I)` in `arr`.
"""
function pairedindices end

"""
    getsubstitutions(idxhandler::AbstractIndexHandler, rs::ReactionSystem; state_sym::Symbol)

Returns a dict of the form `S_i => f_i(state_sym)`, where each `f_i` is an expression
for the abundance of species `S_i` in terms of the state variable `state_sym`.
"""
function getsubstitutions end

"""
    vec(idxhandler::AbstractIndexHandler, arr)

Converts the right-hand side defining the solution of the CME into a
one-dimensional vector to which a matrix can be applied.

See also: [`LinearIndices`](@ref Base.LinearIndices)
"""
function vec end

"""
    LinearIndices(idxhandler::AbstractIndexHandler, arr)

Returns an object `lind` which converts indices returned from [`singleindices`](@ref)
and [`pairedindices`](@ref) to linear indices compatible with [`vec`](@ref Base.vec)
via `lind[idx_cart] = idx_lin`. The indices are related via

```julia
arr[idx_cart] == vec(idxhandler, arr)[idx_lin]
```

See also: [`vec`](@ref Base.vec)
"""
function LinearIndices end

##

"""
    DefaultIndexHandler{N}()
    DefaultIndexHandler{N}(offset, perm)

Default index handler for an FSP system with `N` species. It represents the state
as an `N`-dimensional array, maps a molecule count of zero to `offset`, and uses
`perm` to map state-array dimensions to Catalyst's species order.

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

Deprecated alias for [`DefaultIndexHandler`](@ref). Use `DefaultIndexHandler` instead.
"""
function NaiveIndexHandler(args...; kwargs...)
    Base.depwarn(
        "`NaiveIndexHandler` is deprecated, use `DefaultIndexHandler` instead.",
        :NaiveIndexHandler
    )
    return DefaultIndexHandler(args...; kwargs...)
end
export NaiveIndexHandler

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
