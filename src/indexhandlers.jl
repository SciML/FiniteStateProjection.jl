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
    struct NaiveIndexHandler <: AbstractIndexHandler
        offset::Int
    end

Basic index handler that stores the state of a system with
`s` species in an `s`-dimensional array. The `offset` parameter
denotes the offset by which the array is indexed (defaults to 1
in Julia).

This is the simplest index handler, but it will not be optimal
if some states cannot be reached from the initial state, e.g.
due to the presence of conservation laws. In these cases one should
use `ReducingIndexHandler`, which will automatically elide species
where possible.

Constructors: `NaiveIndexHandler([sys::FSPSystem, offset::Int=1])`

See also: [`ReducingIndexHandler`](@ref)
"""
struct NaiveIndexHandler <: AbstractIndexHandler
    offset::Int
end

NaiveIndexHandler() = NaiveIndexHandler(1)

Base.vec(::NaiveIndexHandler, arr) = vec(arr)
Base.LinearIndices(::NaiveIndexHandler, arr) = LinearIndices(arr)

function pairedindices(ih::NaiveIndexHandler, arr::AbstractArray{T,N},
                       shift::CartesianIndex{N}) where {T,N}
    pairedindices(ih, axes(arr), shift)
end

function pairedindices(ih::NaiveIndexHandler, dims::NTuple{N,T},
                       shift::CartesianIndex{N}) where {N,T<:Number}
    pairedindices(ih, Base.OneTo.(dims), shift)
end

function pairedindices(::NaiveIndexHandler, dims::NTuple{N,T},
                       shift::CartesianIndex{N}) where {N,T<:AbstractVector}
    ranges = tuple((UnitRange(max(first(ax), first(ax)+shift[i]),
                              min(last(ax), last(ax)+shift[i]))
                    for (i, ax) in enumerate(dims))...)

    ranges_shifted = tuple((rng .- shift[i] for (i, rng) in enumerate(ranges))...)

    zip(CartesianIndices(ranges_shifted), CartesianIndices(ranges))
end

"""
    getsubstitutions(sys::FSPSystem{NaiveIndexHandler}; state_sym::Symbol)::Dict

Defines the abundance of species ``S_i`` to be `state_sym[i] - offset`.
"""
function getsubstitutions(ih::NaiveIndexHandler, rs::ReactionSystem; state_sym::Symbol)
    Dict(symbol => Term(Base.getindex, (state_sym, i)) - ih.offset for (i, symbol) in enumerate(species(rs)))
end

##

"""
   struct ReducingIndexHandler <: AbstractIndexHandler

More efficient index handler that improves upon [`NaiveIndexHandler`](@ref)
by eliminating variables whose abundances can be computed from other variables
using conservation laws. Describes the system using a subset of the original
species which can be obtained via [`reducedspecies`](@ref). Reduces the
dimensionality of the FSP by the number of conservation laws in the system.

Note: This feature is currently being moved into Catalyst.jl.

Constructors: `ReducingIndexHandler(rs::ReactionSystem[, offset::Int=1])`

See also: [`reducedspecies`](@ref), [`elidedspecies`](@ref), [`NaiveIndexHandler`](@ref)
"""
struct ReducingIndexHandler <: AbstractIndexHandler
    cons_laws::Matrix{Int}
    offset::Int
    specs_red::Vector{Int}
    specs_elided::Vector{Int}
    cons_syms::Vector{Symbol}
end

function ReducingIndexHandler(rs::ReactionSystem, offset::Int=1)
    cons_laws = conservationlaws(rs)
    specs_elided = elidedspecies(cons_laws)
    specs_red = [ i for i in 1:length(species(rs)) if !(i in specs_elided) ]

    cons_syms = [ gensym("c$i") for i in 1:size(cons_laws, 1) ]

    return ReducingIndexHandler(cons_laws, offset, specs_red, specs_elided, cons_syms)
end

function FSPSystem(rs::ReactionSystem, ih::ReducingIndexHandler; kwargs...)
    rfs = create_ratefuncs(rs, ih; kwargs...)
    FSPSystem(rs, ih, ih.cons_laws, rfs)
end


"""
    reducedspecies(idxhandler::ReducingIndexHandler)

Return indices of reduced species.

See also: [`elidedspecies`](@ref)
"""
reducedspecies(idxhandler::ReducingIndexHandler) = idxhandler.specs_red

"""
    elidedspecies(idxhandler::ReducingIndexHandler)

Return indices of elided species.

See also: [`reducedspecies`](@ref)
"""
elidedspecies(idxhandler::ReducingIndexHandler) = idxhandler.specs_elided

##

"""
    elidedspecies(cons_laws::AbstractMatrix{Int})::Vector

Returns a list of species ``[ s_1, ... ]`` which can be removed from the reaction system
description using the provided matrix of conservation laws.
"""
function elidedspecies(cons_laws::AbstractMatrix{Int})::Vector{Int}
    n_cons, n_spec = size(cons_laws)
    ret = zeros(Int, n_cons)
    for i in 1:n_cons
        cons_law = cons_laws[i,:]
        cons_law_specs = filter(j -> cons_law[j] != 0, 1:n_spec)
        possible_specs = filter(j -> !(j in ret[1:i-1]), cons_law_specs)
        idx = possible_specs[end]
        ret[i] = idx
    end

    @assert length(ret) == length(unique(ret))

    return ret
end

"""
    elisions(idxhandler::ReducingIndexHandler, rs::ReactionSystem)

Replaces the symbols ``A(t)``, ``B(t)``, ... of elided species by
``N_1(t) - X(t) - Y(t)``, ``N_2(t) - U(t) - V(t)``, ..., where ``N_i(t)``
are the conserved quantities of the system.

See also: [`getsubstitutions`](@ref)
"""
function elisions(idxhandler::ReducingIndexHandler, rs::ReactionSystem)
    ret = Dict()
    spec_syms = species(rs)

    for (i, spec) in enumerate(elidedspecies(idxhandler))
        sym = spec_syms[spec]
        cons_law = idxhandler.cons_laws[i,:]

        # Does this always work? What if some of the species on the RHS
        # also end up getting elided at some point?
        rhs = (Variable(idxhandler.cons_syms[i]) - sum(cons_law[j] * spec_syms[j] for j in 1:length(spec_syms) if j != spec))
        rhs /= cons_law[i]

        ret[sym] = rhs
    end

    ret
end

##

"""
    getsubstitutions(idxhandler::ReducingIndexHandler, rs::ReactionSystem; state_sym::Symbol)::Dict

Similar to its [`NaiveIndexHandler`](@ref) variant, but computes the abundances of elided species
from the conserved quantities and the reduced species.
"""
function getsubstitutions(idxhandler::ReducingIndexHandler, rs::ReactionSystem; state_sym::Symbol)
    symbols = states(rs)

    ret = Dict{Any,Any}(symbols[spec] => Term(Base.getindex, (state_sym, i)) - idxhandler.offset
                                              for (i, spec) in enumerate(reducedspecies(idxhandler)))

    elision_dict = elisions(idxhandler, rs)
    for (spec, ex) in pairs(elision_dict)
        ret[spec] = substitute(ex, ret)
    end

    ret
end

"""
    build_rhs_header(idxhandler::ReducingIndexHandler, sys::FSPSystem)::Expr

Assumes `p` is of the form `(params, cons::AbstractVector{Int})` where `params`
are the system parameters and `cons` the conserved quantities.
"""
function build_rhs_header(sys::FSPSystem{ReducingIndexHandler})::Expr
    cons_names = Expr(:tuple, sys.ih.cons_syms...)

    quote
        (ps, $(cons_names)) = p
        $(unpackparams(sys, :ps))
    end
end

"""
    pairedindices(idxhandler::ReducingIndexHandler, arr::AbstractArray, shift::CartesianIndex)

Similar to its `NaiveIndexHandler` variant, but converts the indices into indices into
the reduced state space array.
"""
function pairedindices(idxhandler::ReducingIndexHandler, arr::AbstractArray{T,M},
                       shift::CartesianIndex{N}) where {T,M,N}
    shift_red = CartesianIndex{M}(convert(Tuple, shift)[reducedspecies(idxhandler)]...)
    pairedindices(NaiveIndexHandler(idxhandler.offset), arr, shift_red)
end

function pairedindices(idxhandler::ReducingIndexHandler, dims::NTuple{M},
                       shift::CartesianIndex{N}) where {M,N}
    shift_red = CartesianIndex{M}(convert(Tuple, shift)[reducedspecies(idxhandler)]...)
    pairedindices(NaiveIndexHandler(idxhandler.offset), dims, shift_red)
end
