function getsubstitutions end

""" 
    build_ratefuncs(idxhandler::AbstractIndexHandler, sys::FSPSystem; state_sym::Symbol)::Vector

Return the rate functions converted to Julia expressions in the reduced variables.
We allow arbitrary offsets to deal with Julia's one-based indexing. Using
`OffsetArrays` (optional) is more natural here.

See also: [`build_rhs`](@ref)
"""
function build_ratefuncs(idxhandler::AbstractIndexHandler, sys::FSPSystem; state_sym::Symbol)::Vector
    #symbols = states(sys.rs)
    substitutions = getsubstitutions(idxhandler, sys, state_sym=state_sym)
    
    return [ toexpr(substitute(jumpratelaw(reac), substitutions)) for reac in sys.rs.eqs ]
    
    # Needs to be implemented!
    #idx_subs = Dict(symbols[spec] => Term(Base.getindex, (state_sym, i)) - offset for (i, spec) in enumerate(reduced_species(sys)))
    #for reac in sys.rs.eqs
        #rate_sub = substitute(jumpratelaw(reac), subs_dict)
        #ex = toexpr(substitute(rate_sub, idx_subs))
    #    ex = toexpr(substitute(jumpratelaw(reac), subs))
        
    #    push!(ret, ex)
    #end
    
    #return ret
end

function unpackparams(sys::FSPSystem, psym::Symbol)::Expr
    param_names = Expr(:tuple, map(par -> par.name, params(sys.rs))...)
     
    quote 
        $(param_names) = ps
    end
end

"""
    build_rhs_header(idxhandler::AbstractIndexHandler, sys::FSPSystem)::Expr

Return initialisation code for the RHS function. Unpacks parameters.

See also: [`build_rhs`](@ref)
"""
function build_rhs_header(::AbstractIndexHandler, sys::FSPSystem)::Expr
    #cons_names = Expr(:tuple, sys.cons_syms...)
    
    # Needs modified
    quote 
        #(ps::AbstractVector{Float64}, cons::AbstractVector{Int64}) = p
        ps::AbstractVector{Float64} = p
        $(unpackparams(sys, :ps))
        #$(cons_names) = cons
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
function build_rhs_firstpass(idxhandler::AbstractIndexHandler, sys::FSPSystem, rfs::AbstractVector)::Expr
    isempty(rfs) && return quote end
        
    first_line = :(du[idx_in] = -u[idx_in] * $(rfs[1]))
    other_lines = (:(du[idx_in] -= u[idx_in] * $(rf)) for rf in rfs[2:end])
    
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
function build_rhs_secondpass(idxhandler::AbstractIndexHandler, sys::FSPSystem, rfs::AbstractVector)::Expr
    isempty(rfs) && return quote end
    
    S = netstoichmat(sys.rs)
    ret = Expr(:block)
    
    for (i, rf) in enumerate(rfs)
        ex = quote
            for (idx_in, idx_out) in pairedindices($(idxhandler), u, $(CartesianIndex(S[i,:]...)))
                du[idx_out] += u[idx_in] * $(rf)
            end
        end
        
        append!(ret.args, ex.args)
    end
    
    return ret
end

"""
    build_rhs(idxhandler::AbstractIndexHandler, sys::FSPSystem)

Builds the function `f(du,u,p,t)` that defines the right-hand side of the CME, 
for use in the ODE solver. If `expression` is true, returns an expression, else
compiles the function. If `jac` is true, returns the Jacobian function
instead (I haven't actually tested this).
"""
function build_rhs(idxhandler::AbstractIndexHandler, sys::FSPSystem; expression::Bool=true) 
    rfs = build_ratefuncs(idxhandler, sys, state_sym=:idx_in)
    header = build_rhs_header(idxhandler, sys)

    first_pass = build_rhs_firstpass(idxhandler, sys, rfs)
    second_pass = build_rhs_secondpass(idxhandler, sys, rfs)
    
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
function build_ode_func(idxhandler::AbstractIndexHandler, sys::FSPSystem)::ODEFunction
    rhs = build_rhs(idxhandler, sys, expression=false)
    
    ODEFunction{true}(rhs)
end

"""
    build_ode_problem(sys::FSPSystem, u0::AbstractArray, t, p)

Return an `ODEProblem` for use in `DifferentialEquations. 

`u0` is a multidimensional array denoting initial values for the reduced species.
FSP.jl automatically reduces the dimensionality of the ODE system where possible
and elides species that can be expressed in terms of other species using conservation
laws. The reduced species for a given system can be found with [`reduced_species`].
"""
function build_ode_prob(idxhandler::AbstractIndexHandler, sys::FSPSystem, u0::AbstractArray, tmax::Float64, p)::ODEProblem
    ode_func = build_ode_func(idxhandler, sys)
    
    ODEProblem(ode_func, u0, tmax, p)
end