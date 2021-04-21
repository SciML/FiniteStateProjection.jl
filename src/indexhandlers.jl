abstract type AbstractIndexHandler end

struct NaiveIndexHandler <: AbstractIndexHandler
    offset::Int
end

NaiveIndexHandler(sys::FSPSystem, offset::Int) = NaiveIndexHandler(offset)

function getsubstitutions(idxhandler::NaiveIndexHandler, sys::FSPSystem; state_sym::Symbol)
    Dict{Any,Any}(symbol => Term(Base.getindex, (state_sym, i)) - idxhandler.offset for (i, symbol) in enumerate(species(sys.rs)))
end

""" 
    singleindices(idxhandler::AbstractIndexHandler, arr)

Returns all indices `I` in `arr`. Defaults to CartesianIndices, but can
be overloaded for arbitrary index handlers. 
"""
singleindices(::AbstractIndexHandler, arr::AbstractArray) = CartesianIndices(arr)

""" 
    pairedindices(idxhandler::AbstractIndexHandler, arr, shift::CartesianIndex)

Returns all pairs of indices (`I .- shift`, `I`) such that in `arr`.
The default implementation can be overloaded for arbitrary index handlers. 
"""
function pairedindices(::AbstractIndexHandler, arr::AbstractArray{T,N}, 
                       shift::CartesianIndex{N}) where {T,N}#::Base.Iterators.Zip{NTuple{2, CartesianIndices{N, NTuple{N, UnitRange{Int}}}}} where {T,N}
    ranges = tuple((UnitRange(max(first(ax), first(ax)+shift[i]), 
                              min(last(ax), last(ax)+shift[i])) 
                    for (i, ax) in enumerate(axes(arr)))...)
    
    ranges_shifted = tuple((rng .- shift[i] for (i, rng) in enumerate(ranges))...)
    
    zip(CartesianIndices(ranges_shifted), CartesianIndices(ranges))
end

####

struct DefaultIndexHandler <: AbstractIndexHandler
    offset::Int
    specs_red::Vector{Int}
    specs_elided::Vector{Int}
    cons_syms::Vector{Symbol}
end

function DefaultIndexHandler(sys::FSPSystem, offset::Int)
    cons_laws = conservationlaws(sys)
    specs_elided = _elidedspecies(cons_laws)
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

""" 
    _elidedspecies(cons_laws)::Vector

Returns a list of species ``[ s_1, ... ]`` which can be removed from the reaction system
description using the provided matrix of conservation laws.
"""
function _elidedspecies(cons_laws::AbstractMatrix{Int})::Vector{Int}
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

function getsubstitutions(idxhandler::DefaultIndexHandler, sys::FSPSystem; state_sym::Symbol)
    symbols = states(sys.rs)
    
    ret = Dict{Any,Any}(symbols[spec] => Term(Base.getindex, (state_sym, i)) - idxhandler.offset for (i, spec) in enumerate(reducedspecies(idxhandler)))
    
    elision_dict = elisions(idxhandler, sys)
    for (spec, ex) in pairs(elision_dict)
        ret[spec] = substitute(ex, ret)
    end
    
    ret
end

"""
    build_rhs_header(sys::FSPSystem)::Expr

Return initialisation code for the RHS function. Unpacks 
conserved quantities and parameters.

See also: [`build_rhs`](@ref)
"""
function build_rhs_header(idxhandler::DefaultIndexHandler, sys::FSPSystem)::Expr
    cons_names = Expr(:tuple, idxhandler.cons_syms...)
    
    quote 
        (ps::AbstractVector{Float64}, cons::AbstractVector{Int64}) = p
        $(unpackparams(sys, :ps))
        $(cons_names) = cons
    end
end

function pairedindices(idxhandler::DefaultIndexHandler, arr::AbstractArray{T,M}, 
                        shift::CartesianIndex{N}) where {T,M,N}#::Base.Iterators.Zip{NTuple{2, CartesianIndices{N, NTuple{N, UnitRange{Int}}}}} where {T,N}
    shift_red = CartesianIndex{M}(convert(Tuple, shift)[reducedspecies(idxhandler)]...)
    pairedindices(NaiveIndexHandler(idxhandler.offset), arr, shift_red)
end