abstract type AbstractIndexHandler end

struct NaiveIndexHandler <: AbstractIndexHandler
    offset::Int
end
    
""" 
    pairedindices(idxhandler, arr, shift::CartesianIndex)::CartesianIndices

Returns all Cartesian indices `I` such that `I` and `I .+ shift` are in `arr`.
The default implementation can be overloaded for arbitrary index handlers. 
"""
function pairedindices(::AbstractIndexHandler, arr::AbstractArray{T,N}, 
                        shift::CartesianIndex{N})::CartesianIndices{N, NTuple{N, UnitRange{Int}}} where {T,N}
    ranges = tuple((UnitRange(max(first(ax), first(ax)+shift[i]), 
                              min(last(ax), last(ax)+shift[i])) 
                    for (i, ax) in enumerate(axes(arr)))...)
    
    zip(CartesianIndices(ranges .- tuple(shift)), CartesianIndices(ranges))
end

# """ 
#     get_elided_species(cons_laws)::Vector

# Returns a list of species ``[ s_1, ... ]`` which can be removed from the reaction system
# description using the provided matrix of conservation laws. It is important that this
# function is fully deterministic (for now, otherwise bugs might occur)
# """
# function get_elided_species(cons_laws::AbstractMatrix{Int})::Vector{Int}
#     n_cons, n_spec = size(cons_laws)
#     ret = zeros(Int, n_cons)
#     for i in 1:n_cons
#         cons_law = cons_laws[i,:]
#         cons_law_specs = filter(j -> cons_law[j] != 0, 1:n_spec)
#         possible_specs = filter(j -> !(j in ret[1:i-1]), cons_law_specs)
#         idx = possible_specs[end]
#         ret[i] = idx
#     end
    
#     @assert length(ret) == length(unique(ret))
    
#     return ret
# end

# """
#     reduced_species(sys::FSPSystem)

# Return indices of reduced species for the reaction system.

# See also: [`elided_species`](@ref)
# """
# reduced_species(sys::FSPSystem) = sys.specs_sep[1:size(sys.cons_laws,2) - size(sys.cons_laws,1)]

# """
#     elided_species(sys::FSPSystem)

# Return indices of elided (removed) species for the reaction system.

# See also: [`reduced_species`](@ref)
# """
# elided_species(sys::FSPSystem) = isempty(sys.cons_laws) ? [] : sys.specs_sep[end-size(sys.cons_laws,1)+1,end]

# function FSPSystem(rs::ReactionSystem)
#     cons_laws = get_conservation_laws(netstoichmat(rs))
    
#     n_specs = length(species(rs))
#     n_cons = size(cons_laws, 1)
#     n_red_specs = n_specs - n_cons
#     specs_sep = zeros(Int, n_specs)
#     specs_elided = @view specs_sep[n_red_specs+1:end]
#     specs_elided .= get_elided_species(cons_laws)
#     specs_sep[1:n_red_specs] .= [ i for i in 1:n_specs if !(i in specs_elided) ]
    
#     cons_syms = [ gensym("c$i") for i in 1:n_cons ]
    
#     return FSPSystem(rs, cons_laws, specs_sep, cons_syms)
# end

# """ 
#     get_subs_dict(sys)::Dict

# Replaces the symbols ``A(t)``, ``B(t)``, ... of elided species by
# ``N_1(t) - X(t) - Y(t)``, ``N_2(t) - U(t) - V(t)``, ..., where ``N_i(t)``
# are the conserved quantities of the system.

# See also: [`build_ratefuncs`](@ref)
# """
# function getsubs(::NaiveIndexHandler, sys::FSPSystem)::Dict{Any,Any}
#     ret = Dict()
#     spec_syms = species(sys.rs)
    
#     for (i, spec) in enumerate(elided_species(sys))
#         sym = spec_syms[spec]
#         cons_law = sys.cons_laws[i,:]
        
#         # Does this always work? What if some of the species on the RHS
#         # also end up getting elided at some point?
#         rhs = (Variable(sys.cons_syms[i]) - sum(cons_law[j] * spec_syms[j] for j in 1:length(spec_syms) if j != spec))
#         rhs /= cons_law[i]
        
#         ret[sym] = rhs
#     end
    
#     ret
# end

function getsubstitutions(idxhandler::NaiveIndexHandler, sys::FSPSystem; state_sym::Symbol)
    Dict(symbol => Term(Base.getindex, (state_sym, i)) - idxhandler.offset for (i, symbol) in enumerate(states(sys.rs)))
end

# """
#     build_rhs_header(sys::FSPSystem)::Expr

# Return initialisation code for the RHS function. Unpacks 
# conserved quantities and parameters.

# See also: [`build_rhs`](@ref)
# """
# function build_rhs_header(sys::FSPSystem)::Expr
#     cons_names = Expr(:tuple, sys.cons_syms...)
    
#     quote 
#         (ps::AbstractVector{Float64}, cons::AbstractVector{Int64}) = p
#         $(unpackparams(sys, :ps))
#         $(cons_names) = cons
#     end
# end
