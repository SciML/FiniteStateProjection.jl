""" 
    abstract type AbstractIndexHandler end

`FSP.jl` splits handling of the FSP into two parts. The first defines
how the CME is compute

See also: [`singleindices`](@ref), [`pairedindices`](@ref)
"""
abstract type AbstractIndexHandler end

""" 
    singleindices(idxhandler::AbstractIndexHandler, arr)

Returns all indices `I` in `arr`. Defaults to CartesianIndices, but can
be overloaded for arbitrary index handlers. 
"""
singleindices(::AbstractIndexHandler, arr::AbstractArray) = CartesianIndices(arr)

""" 
    pairedindices(idxhandler::AbstractIndexHandler, arr, shift::CartesianIndex)

Returns all pairs of indices `(I .- shift, I)` in `arr`.
The default implementation can be overloaded for arbitrary index handlers. 
"""
function pairedindices(::AbstractIndexHandler, arr::AbstractArray{T,N}, 
                       shift::CartesianIndex{N}) where {T,N}
    ranges = tuple((UnitRange(max(first(ax), first(ax)+shift[i]), 
                              min(last(ax), last(ax)+shift[i])) 
                    for (i, ax) in enumerate(axes(arr)))...)
    
    ranges_shifted = tuple((rng .- shift[i] for (i, rng) in enumerate(ranges))...)
   
    zip(CartesianIndices(ranges_shifted), CartesianIndices(ranges))
end

##

""" 
    struct NaiveIndexHandler <: AbstractIndexHandler
        offset::Int
    end

Basic index handler that stores the state of a system with
`s` species in an `s`-dimensional array. The `offset` parameter
denotes the offset by which the array is indexed (defaults to 1
in Julia). Use `OffsetArrays.jl` to enable 0-based indexing.

This is the simplest index handler, but it will not be optimal
if some states cannot be reached from the initial state, e.g.
due to the presence of conservation laws. It is generally better
to use `DefaultIndexHandler`, which will automatically elide species
where possible.

Constructors: `NaiveIndexHandler([sys::FSPSystem,] offset::Int)`

See also: [`DefaultIndexHandler`](@ref)
"""
struct NaiveIndexHandler <: AbstractIndexHandler
    offset::Int
end

NaiveIndexHandler(sys::FSPSystem, offset::Int) = NaiveIndexHandler(offset)

"""
    getsubstitutions(idxhandler::NaiveIndexHandler, sys::FSPSystem; state_sym::Symbol)::Dict

Defines the abundance of species ``S_i`` to be `state_sym[i] - offset`.

See also: [`NaiveIndexHandler`](@ref)
"""
function getsubstitutions(idxhandler::NaiveIndexHandler, sys::FSPSystem; state_sym::Symbol)
    Dict(symbol => Term(Base.getindex, (state_sym, i)) - idxhandler.offset for (i, symbol) in enumerate(species(sys.rs)))
end

##

""" 
    struct DefaultIndexHandler <: AbstractIndexHandler

More efficient index handler that improves upon [`NaiveIndexHandler`](@ref)
by eliminating variables whose abundances can be computed from other variables
using conservation laws. Describes the system using a subset of the original
species which can be obtained via [`reducedspecies`](@ref). Reduces the 
dimensionality of the FSP by the number of conservation laws in the system.

Constructors: `DefaultIndexHandler(sys::FSPSystem, offset::Int)`

See also: [`reducedspecies`](@ref), [`elidedspecies`](@ref), [`NaiveIndexHandler`](@ref)
"""
struct DefaultIndexHandler <: AbstractIndexHandler
    offset::Int
    specs_red::Vector{Int}
    specs_elided::Vector{Int}
    cons_syms::Vector{Symbol}
end

function DefaultIndexHandler(sys::FSPSystem, offset::Int)
    cons_laws = conservationlaws(sys)
    specs_elided = elidedspecies(cons_laws)
    specs_red = [ i for i in 1:length(species(sys.rs)) if !(i in specs_elided) ]
    
    cons_syms = [ gensym("c$i") for i in 1:size(cons_laws, 1) ]
    
    return DefaultIndexHandler(offset, specs_red, specs_elided, cons_syms)
end

"""
    reducedspecies(idxhandler::DefaultIndexHandler)

Return indices of reduced species.

See also: [`elidedspecies`](@ref)
"""
reducedspecies(idxhandler::DefaultIndexHandler) = idxhandler.specs_red

"""
    elidedspecies(idxhandler::DefaultIndexHandler)

Return indices of elided species.
"""
elidedspecies(idxhandler::DefaultIndexHandler) = idxhandler.specs_elided

##

""" 
    elidedspecies(cons_laws)::Vector

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
Replaces the symbols ``A(t)``, ``B(t)``, ... of elided species by
``N_1(t) - X(t) - Y(t)``, ``N_2(t) - U(t) - V(t)``, ..., where ``N_i(t)``
are the conserved quantities of the system.

See also: [`build_ratefuncs`](@ref)
"""
function elisions(idxhandler::DefaultIndexHandler, sys::FSPSystem)
    ret = Dict()
    spec_syms = species(sys.rs)
    
    for (i, spec) in enumerate(elidedspecies(idxhandler))
        sym = spec_syms[spec]
        cons_law = sys.cons_laws[i,:]
        
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
    getsubstitutions(idxhandler::DefaultIndexHandler, sys::FSPSystem; state_sym::Symbol)::Dict

Similar to its [`NaiveIndexHandler`](@ref) variant, but computes the abundances of elided species
from the conserved quantities and the reduced species.

See also: [`DefaultIndexHandler`](@ref)
"""
function getsubstitutions(idxhandler::DefaultIndexHandler, sys::FSPSystem; state_sym::Symbol)
    symbols = states(sys.rs)
    
    ret = Dict{Any,Any}(symbols[spec] => Term(Base.getindex, (state_sym, i)) - idxhandler.offset 
                                              for (i, spec) in enumerate(reducedspecies(idxhandler)))
    
    elision_dict = elisions(idxhandler, sys)
    for (spec, ex) in pairs(elision_dict)
        ret[spec] = substitute(ex, ret)
    end
    
    ret
end

"""
    build_rhs_header(idxhandler::DefaultIndexHandler, sys::FSPSystem)::Expr

Assumes `p` is of the form `(params, cons::AbstractVector{Int})` where `params` 
are the system parameters and `cons` the conserved quantities.

See also: [`build_rhs_header`](@ref)
"""
function build_rhs_header(idxhandler::DefaultIndexHandler, sys::FSPSystem)::Expr
    cons_names = Expr(:tuple, idxhandler.cons_syms...)
    
    quote 
        (ps::AbstractVector{Float64}, cons::AbstractVector{Int}) = p
        $(unpackparams(sys, :ps))
        $(cons_names) = cons
    end
end

"""
    pairedindices(idxhandler::DefaultIndexHandler, arr::AbstractArray, shift::CartesianIndex)

Similar to its `NaiveIndexHandler` variant, but converts the indices into indices into
the reduced state space array.

See also: ['pairedindices`](@ref)
"""
function pairedindices(idxhandler::DefaultIndexHandler, arr::AbstractArray{T,M}, 
                        shift::CartesianIndex{N}) where {T,M,N}
    shift_red = CartesianIndex{M}(convert(Tuple, shift)[reducedspecies(idxhandler)]...)
    pairedindices(NaiveIndexHandler(idxhandler.offset), arr, shift_red)
end
