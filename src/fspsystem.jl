import AbstractAlgebra

netstoichmat(rs::ReactionSystem) = prodstoichmat(rs) - substoichmat(rs)

struct FSPSystem
    rs::ReactionSystem
    cons_laws::Matrix{Int}
end

function FSPSystem(rs::ReactionSystem)
    FSPSystem(rs, conservationlaws(netstoichmat(rs)))
end

""" 
    conservationlaws(netstoichmat::AbstractMatrix{Int})::Matrix{Int}
    conservationlaws(sys::FSPSystem)::Matrix{Int}

Given the net stoichiometry matrix of a reaction system, computes a matrix
of conservation laws. Each row contains the stoichiometric coefficients
of a different conserved quantity.
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

conservationlaws(sys::FSPSystem) = sys.cons_laws

"""
    conservedquantities(state, sys::FSPSystem)

Compute conserved quantities for the system at the given state.
"""
conservedquantities(state::AbstractVector{Int}, sys::FSPSystem) = conservationlaws(sys) * state

