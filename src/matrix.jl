function create_sparsematrix(sys::FSPSystem, dims::NTuple, ps, t)
    Ntot = prod(dims)
    lind = LinearIndices(sys.ih, dims)

    I = Int[]
    J = Int[]
    V = Float64[]

    predsize = Ntot * (length(Catalyst.get_eqs(sys.rs)) + 1)

    sizehint!(I, predsize)
    sizehint!(J, predsize)
    sizehint!(V, predsize)

    V = zeros(predsize)

    for idx_cart in singleindices(sys.ih, dims)
        idx_lin = lind[idx_cart]
        push!(I, idx_lin)
        push!(J, idx_lin)

        k = length(I)
        for rf in sys.rfs
            V[k] -= rf(idx_cart, t, ps...)
        end
    end

    S::Matrix{Int64} = netstoichmat(sys.rs)
    for (i, rf) in enumerate(sys.rfs)
        for (idx_cin, idx_cout) in pairedindices(sys.ih, dims, CartesianIndex(S[i,:]...))
            idx_lin = lind[idx_cin]
            idx_lout = lind[idx_cout]
            push!(I, lind[idx_cout])
            push!(J, lind[idx_cin])

            k = length(I)
            V[k] = rf(idx_cin, t, ps...)
        end
    end

    resize!(V, length(I))
    sparse(I, J, V)
end

"""
    Base.convert(::Type{SparseMatrixCSC}, ::FSPSystem, dims::NTuple, ps, t)

Convert the reaction system into a sparse matrix defining the right-hand side of the 
Chemical Master Equation. `dims` is a tuple denoting the dimensions of the FSP and 
`ps` is the tuple of parameters. The sparse matrix works on the flattened version
of the state obtained using `vec`.
"""
function Base.convert(::Type{SparseMatrixCSC}, sys::FSPSystem, dims::NTuple, ps, t)
    return create_sparsematrix(sys, dims, ps, t)
end
