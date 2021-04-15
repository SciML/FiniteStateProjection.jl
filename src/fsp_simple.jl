import AbstractAlgebra

netstoichmat(rs::ReactionSystem) = prodstoichmat(rs) - substoichmat(rs)

""" 
    shifted_indices(arr, shift::CartesianIndex)::CartesianIndices

Returns all Cartesian indices `I` such that `I` and `I .+ shift` are in `arr`
"""
function shifted_indices(arr::AbstractArray{T,N}, shift::CartesianIndex{N})::CartesianIndices{N, NTuple{N, UnitRange{Int}}} where {T,N}
    ranges = tuple((UnitRange(max(first(ax), first(ax)+shift[i]), 
                              min(last(ax), last(ax)+shift[i])) 
                    for (i, ax) in enumerate(axes(arr)))...)
    
    CartesianIndices(ranges)
end

""" 
    get_conservation_laws(netstoichmat::AbstractMatrix{Int})::Matrix{Int}

Given the net stoichiometry matrix of a reaction system, computes a matrix
of conservation laws. Each row contains the stoichiometric coefficients
of a different conserved quantity.
"""
function get_conservation_laws(nsm::AbstractMatrix{Int})::Matrix{Int}
    n_reac, n_spec = size(nsm)
    
    # We basically have to compute the left null space of the matrix
    # over the integers; this is best done using its Smith Normal Form.
    nsm_conv = AbstractAlgebra.matrix(AbstractAlgebra.ZZ, nsm)
    S, T, U = AbstractAlgebra.snf_with_transform(nsm_conv)
    
    # Zero columns of S (which occur after nonzero columns in SNF)
    # correspond to conserved quantities
    n = findfirst(i -> all(S[:,i] .== 0), 1:n_spec)
    if n === nothing
        return zeros(Int, 0, n_spec)
    end
    
    return Matrix(U[:,n:end]')
end

""" 
    get_elided_species(cons_laws)::Vector

Returns a list of species ``[ s_1, ... ]`` which can be removed from the reaction system
description using the provided matrix of conservation laws. It is important that this
function is fully deterministic (for now, otherwise bugs might occur)
"""
function get_elided_species(cons_laws::AbstractMatrix{Int})::Vector{Int}
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

struct FSPSystem
    rs::ReactionSystem
    cons_laws::Matrix{Int}
    specs_sep::Vector{Int}
    cons_syms::Vector{Symbol}
end

"""
    get_conserved_quantities(state, sys::FSPSystem)

Compute conserved quantities for the system at the given state.
"""
get_conserved_quantities(state::AbstractVector{Int}, sys::FSPSystem)::AbstractVector{Int} = 
    sys.cons_laws * state

"""
    reduced_species(sys::FSPSystem)

Return indices of reduced species for the reaction system.

See also: [`elided_species`](@ref)
"""
reduced_species(sys::FSPSystem) = sys.specs_sep[1:size(sys.cons_laws,2) - size(sys.cons_laws,1)]

"""
    elided_species(sys::FSPSystem)

Return indices of elided (removed) species for the reaction system.

See also: [`reduced_species`](@ref)
"""
elided_species(sys::FSPSystem) = isempty(sys.cons_laws) ? [] : sys.specs_sep[end-size(sys.cons_laws,1)+1,end]

function FSPSystem(rs::ReactionSystem)
    cons_laws = get_conservation_laws(netstoichmat(rs))
    
    n_specs = length(species(rs))
    n_cons = size(cons_laws, 1)
    n_red_specs = n_specs - n_cons
    specs_sep = zeros(Int, n_specs)
    specs_elided = @view specs_sep[n_red_specs+1:end]
    specs_elided .= get_elided_species(cons_laws)
    specs_sep[1:n_red_specs] .= [ i for i in 1:n_specs if !(i in specs_elided) ]
    
    cons_syms = [ gensym("c$i") for i in 1:n_cons ]
    
    return FSPSystem(rs, cons_laws, specs_sep, cons_syms)
end

""" 
    get_subs_dict(sys)::Dict

Replaces the symbols ``A(t)``, ``B(t)``, ... of elided species by
``N_1(t) - X(t) - Y(t)``, ``N_2(t) - U(t) - V(t)``, ..., where ``N_i(t)``
are the conserved quantities of the system.

See also: [`build_ratefuncs`](@ref)
"""
function get_subs_dict(sys::FSPSystem)::Dict{Any,Any}
    ret = Dict()
    spec_syms = species(sys.rs)
    
    for (i, spec) in enumerate(elided_species(sys))
        sym = spec_syms[spec]
        cons_law = sys.cons_laws[i,:]
        
        # Does this always work? What if some of the species on the RHS
        # also end up getting elided at some point?
        rhs = (Variable(sys.cons_syms[i]) - sum(cons_law[j] * spec_syms[j] for j in 1:length(spec_syms) if j != spec))
        rhs /= cons_law[i]
        
        ret[sym] = rhs
    end
    
    ret
end

""" 
    build_ratefuncs(sys::FSPSystem)::Vector

Return the rate functions converted to Julia expressions in the reduced variables.
We allow arbitrary offsets to deal with Julia's one-based indexing. Using
`OffsetArrays` (optional) is more natural here.

See also: [`build_rhs`](@ref)
"""
function build_ratefuncs(sys::FSPSystem; state_sym=:idx_in, offset=0)::Vector
    symbols = states(sys.rs)
    subs_dict = get_subs_dict(sys)
    
    ret = []
    idx_subs = Dict(symbols[spec] => Term(:(Base.getindex), (state_sym, i)) - offset for (i, spec) in enumerate(reduced_species(sys)))
    for reac in sys.rs.eqs
        rate_sub = substitute(jumpratelaw(reac), subs_dict)
        ex = toexpr(substitute(rate_sub, idx_subs))
        
        push!(ret, ex)
    end
    
    return ret
end

"""
    build_rhs_header(sys::FSPSystem)::Expr

Return initialisation code for the RHS function. Unpacks 
conserved quantities and parameters.

See also: [`build_rhs`](@ref)
"""
function build_rhs_header(sys::FSPSystem)::Expr
    param_names = Expr(:tuple, map(par -> par.name, params(sys.rs))...)
    cons_names = Expr(:tuple, sys.cons_syms...)
    
    quote 
        (ps::AbstractVector{Float64}, cons::AbstractVector{Int64}) = p
        $(param_names) = ps
        $(cons_names) = cons
    end
end

"""
    build_rhs_firstpass(sys::FSPSystem, rfs)::Expr

Return code for the first pass of the RHS function. Goes through
all reactions and computes the negative part of the CME (probability
flowing out of states). This is a simple array traversal and can be
done in one go for all reactions.

See also: [`build_rhs`](@ref)
"""
function build_rhs_firstpass(sys::FSPSystem, rfs::AbstractVector; jac::Bool=false)::Expr
    # The CME is linear in u; computing the Jacobian is equivalent to
    # dropping the references to u
    jacornot = jac ? 1 : (:(u[idx_in]))
    
    # This would be kinda pointless
    if isempty(rfs) 
        return quote end
    end
        
    first_line = :(du[idx_in] = -$(jacornot) * $(rfs[1]))
    other_lines = (:(du[idx_in] -= $(jacornot) * $rf) for rf in rfs[2:end])
    
    quote
        for idx_in in CartesianIndices(u)
            $first_line
            $(other_lines...)
        end
    end
end

"""
    build_rhs_secondpass(sys::FSPSystem, rfs)::Expr

Return code for the second pass of the RHS function. Goes through
all reactions and computes the positive part of the CME (probability
flowing into states). This requires accessing `du` and `u` at different
locations depending on the net stoichiometries. In order to reduce 
random memory access reactions are processed one by one.

See also: [`build_rhs`](@ref)
"""
function build_rhs_secondpass(sys::FSPSystem, rfs::AbstractVector; jac::Bool=false)::Expr
    S = netstoichmat(sys.rs)
    S = S[:,reduced_species(sys)]
    
    jacornot = jac ? 1 : (:(u[idx_in]))
    
    ret = Expr(:block)
    
    for (i, rf) in enumerate(rfs)
        ex = quote
            for idx_out in shifted_indices(u, $(CartesianIndex(S[i,:]...)))
                idx_in = idx_out - $(CartesianIndex(S[i,:]...))
                du[idx_out] += $(jacornot) * $(rf)
            end
        end
        
        append!(ret.args, ex.args)
    end
    
    return ret
end

"""
    build_rhs(sys::FSPSystem)

Builds the function `f(du,u,p,t)` that defines the right-hand side of the CME, 
for use in the ODE solver. If `expression` is true, returns an expression, else
compiles the function. If `jac` is true, returns the Jacobian function
instead (I haven't actually tested this).
"""
function build_rhs(sys::FSPSystem; expression::Bool=true, jac::Bool=false, offset=0) 
    rfs = build_ratefuncs(sys, offset=offset)
    header = build_rhs_header(sys)

    first_pass = build_rhs_firstpass(sys, rfs, jac=jac)
    second_pass = build_rhs_secondpass(sys, rfs, jac=jac)
    
    args = Expr(:tuple, :du, :u, :p, :t)
    body = Expr(:block, header, first_pass, second_pass)
    
    ex = Expr(:function, args, body) |> MacroTools.striplines |> 
                                        MacroTools.flatten |> 
                                        MacroTools.alias_gensyms
    
    if expression
        return ex
    else
        return @RuntimeGeneratedFunction(ex)
    end
end

"""
    build_ode_func(sys::FSPSystem)

Return an ODEFunction defining the right-hand side of the CME.

Combines the RHS func and its Jacobian to define an `ODEFunction` for 
use with `DifferentialEquations`.

See also: [`build_ode_prob`](@ref)
"""
function build_ode_func(sys::FSPSystem; offset=0)::ODEFunction
    rhs = build_rhs(sys, expression=false, jac=false, offset=offset)
    rhs_jac = build_rhs(sys, expression=false, jac=true, offset=offset)
    
    ODEFunction{true}(rhs, jac=rhs_jac)
end

"""
    build_ode_problem(sys::FSPSystem, u0::AbstractArray, t, p)

Return an `ODEProblem` for use in `DifferentialEquations. 

`u0` is a multidimensional array denoting initial values for the reduced species.
FSP.jl automatically reduces the dimensionality of the ODE system where possible
and elides species that can be expressed in terms of other species using conservation
laws. The reduced species for a given system can be found with [`reduced_species`].
"""
function build_ode_prob(sys::FSPSystem, u0::AbstractArray, tmax::Float64, ps, cons::AbstractVector{Int64}=Int64[])::ODEProblem
    ode_func = build_ode_func(sys, offset=firstindex(u0, 1))
    
    ODEProblem(ode_func, u0, tmax, (ps, cons))
end