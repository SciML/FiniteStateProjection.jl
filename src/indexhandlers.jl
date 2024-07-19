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
    struct DefaultIndexHandler{N} <: AbstractIndexHandler
        offset::Int
        perm::NTuple{N,Int}
    end

Basic index handler that stores the state of a system with
`s` species in an `s`-dimensional array. The `offset` parameter
denotes the offset by which the array is indexed (defaults to 1
in Julia). The order of the species is given by the tuple `perm`.

This is the simplest index handler, but it will not be optimal
if some states cannot be reached from the initial state, e.g.
due to the presence of conservation laws. In these cases one should
use try to remove redundant species where possible.

Constructors: `DefaultIndexHandler([sys::FSPSystem, offset::Int=1])`
"""
struct DefaultIndexHandler{N} <: AbstractIndexHandler
    offset::Int
    perm::NTuple{N,Int}
end

DefaultIndexHandler{N}() where {N} = DefaultIndexHandler{N}(1, Tuple(1:N))

@deprecate NaiveIndexHandler DefaultIndexHandler true

Base.vec(::DefaultIndexHandler, arr) = vec(arr)
Base.LinearIndices(::DefaultIndexHandler, arr) = LinearIndices(arr)

function pairedindices(ih::DefaultIndexHandler{N}, arr::AbstractArray{T,N},
                       shift::CartesianIndex{N}) where {T,N}
    pairedindices(ih, axes(arr), shift)
end

function pairedindices(ih::DefaultIndexHandler{N}, dims::NTuple{N,T},
                       shift::CartesianIndex{N}) where {N,T<:Number}
    pairedindices(ih, Base.OneTo.(dims), shift)
end

# Important: the species in `shift` are ordered according to `Catalyst.species`!
function pairedindices(ih::DefaultIndexHandler{N}, dims::NTuple{N,T},
                       shift::CartesianIndex{N}) where {N,T<:AbstractVector}
    ranges = tuple((UnitRange(max(first(ax), first(ax)+shift[ih.perm[i]]),
                              min(last(ax), last(ax)+shift[ih.perm[i]]))
                    for (i, ax) in enumerate(dims))...)

    ranges_shifted = tuple((rng .- shift[ih.perm[i]] for (i, rng) in enumerate(ranges))...)

    zip(CartesianIndices(ranges_shifted), CartesianIndices(ranges))
end

function pairedindices(::DefaultIndexHandler, dims::NTuple{N,T},
                       shift::CartesianIndex{M}) where {N,M,T<:AbstractVector}
    @error "Dimension of state space ($(length(dims))) does not match number of species ($(length(shift)))"
end


"""
    getsubstitutions(sys::FSPSystem{DefaultIndexHandler}; state_sym::Symbol)::Dict

Defines the abundance of species ``S_i`` to be `state_sym[i] - offset`.
"""
function getsubstitutions(ih::DefaultIndexHandler, rs::ReactionSystem; state_sym::Symbol)
    nspecs = numspecies(rs)
    state_sym_vec = ModelingToolkit.value.(ModelingToolkit.scalarize((@variables ($state_sym)[1:nspecs])[1]))

    species_orig = species(rs)
    species_perm = [ species_orig[ih.perm[i]] for i in 1:nspecs ]

    Dict(symbol => state_sym_vec[i] - ih.offset for (i, symbol) in enumerate(species_perm))
end

#"""
#    PermutingIndexHandler(rs::ReactionSystem, order::AbstractVector)
#
#Constructs an index handler for the reaction system in which the species appear in the order
#defined by the vector `order`.
#"""
function PermutingIndexHandler(rs::ReactionSystem, order::AbstractVector{Symbol})
    PermutingIndexHandler(rs, map(sym -> Catalyst._symbol_to_var(rs, sym), order))
end

function PermutingIndexHandler(rs::ReactionSystem, order::AbstractVector)
    spec = Catalyst.species(rs)
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

    DefaultIndexHandler(1, Tuple(perm))
end
