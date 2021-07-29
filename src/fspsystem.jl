import AbstractAlgebra

netstoichmat(rs::ReactionSystem) = prodstoichmat(rs) - substoichmat(rs)

""" 
    struct FSPSystem{IHT <: AbstractIndexHandler}
        rs::ReactionSystem
        ih::IHT
        cons_laws::Matrix{Int}
    end

Thin wrapper around `ModelingToolkit.ReactionSystem` for use with this package. 
Automatically computes a matrix of conservation laws (see 
[`conservationlaws`](@ref)). 

Constructor: `FSPSystem(rs::ReactionSystem)`
"""
struct FSPSystem{IHT <: AbstractIndexHandler, RT}
    rs::ReactionSystem
    ih::IHT
    cons_laws::Matrix{Int}
    rfs::RT
end

function FSPSystem(rs::ReactionSystem, ih=NaiveIndexHandler(); combinatoric_ratelaw::Bool=true)
    rfs = create_ratefuncs(rs, ih; combinatoric_ratelaw=combinatoric_ratelaw)
    FSPSystem(rs, ih, conservationlaws(rs), rfs)
end

""" 
    build_ratefuncs(rs, ih; 
                    state_sym::Symbol, combinatoric_ratelaw::Bool)::Vector

Return the rate functions converted to Julia expressions in the state variable 
`state_sym`. Abundances of the species are computed using `getsubstitutions`.

See also: [`getsubstitutions`](@ref), [`build_rhs`](@ref)
"""
function build_ratefuncs(rs::ReactionSystem, ih::AbstractIndexHandler; state_sym::Symbol, combinatoric_ratelaw::Bool=true)
    substitutions = getsubstitutions(ih, rs, state_sym=state_sym)
    
    return map(Catalyst.get_eqs(rs)) do reac
        jrl = jumpratelaw(reac; combinatoric_ratelaw)
        jrl_s = substitute(jrl, substitutions)
        toexpr(jrl_s)
    end
end

function create_ratefuncs(rs::ReactionSystem, ih::AbstractIndexHandler; combinatoric_ratelaw::Bool=true)
    paramsyms = Symbol.(Catalyst.params(rs))
    
    return tuple(map(ex -> compile_ratefunc(ex, paramsyms), 
                     build_ratefuncs(rs, ih; state_sym=:idx_in, combinatoric_ratelaw))...)
end 

function compile_ratefunc(ex_rf, params) 
    # Make this nicer in the future
    ex = :((idx_in, t, $(params...)) -> $(ex_rf)) |> MacroTools.flatten
    @RuntimeGeneratedFunction(ex)
end


""" 
    conservationlaws(rs::ReactionSystem)::Matrix{Int}
    conservationlaws(netstoichmat::AbstractMatrix{Int})::Matrix{Int}

Given the net stoichiometry matrix of a reaction system, computes a matrix
of conservation laws (each represented as a row in the output). 
"""
function conservationlaws(nsm::AbstractMatrix{Int})::Matrix{Int}
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
    
    ret = Matrix(U[:,n:end]')
    
    # If all coefficients for a conservation law are negative
    # we might as well flip them to become positive
    for i in 1:size(ret,1)
        all(ret[i,:] .<= 0) && (ret[i,:] .*= -1)
    end
    
    ret
end

conservationlaws(rs::ReactionSystem) = conservationlaws(netstoichmat(rs))

"""
    conservationlaws(sys::FSPSystem)::Matrix{Int}

Returns matrix of conservation laws associated with the system. Each row of 
the matrix contains the stoichiometric coefficients of a different conserved quantity.
"""
conservationlaws(sys::FSPSystem) = sys.cons_laws

"""
    conservedquantities(state, sys)

Compute conserved quantities for the system at the given state.
"""
conservedquantities(state::AbstractVector, sys) = conservationlaws(sys) * state

