function create_sparsematrix(sys::FSPSystem, dims::NTuple, ps, t)
    Ntot = prod(dims)
    lind = LinearIndices(sys.ih, dims)

    I = Int[]
    J = Int[]
    V = Float64[]

    predsize = Ntot * (length(Catalyst.equations(sys.rs)) + 1)

    sizehint!(I, predsize)
    sizehint!(J, predsize)
    sizehint!(V, predsize)

    for idx_cart in singleindices(sys.ih, dims)
        idx_lin = lind[idx_cart]
        push!(I, idx_lin)
        push!(J, idx_lin)

        rate = 0.0
        for rf in sys.rfs
            rate -= rf(idx_cart, t, ps...)
        end

        push!(V, rate)
    end

    S::Matrix{Int64} = netstoichmat(sys.rs)
    for (i, rf) in enumerate(sys.rfs)
        for (idx_cin, idx_cout) in pairedindices(sys.ih, dims, CartesianIndex(S[:,i]...))
            idx_lin = lind[idx_cin]
            idx_lout = lind[idx_cout]
            push!(I, lind[idx_cout])
            push!(J, lind[idx_cin])

            rate = rf(idx_cin, t, ps...)
            push!(V, rate)
        end
    end

    sparse(I, J, V)
end

function create_sparsematrix_ss(sys::FSPSystem, dims::NTuple, ps)
    Ntot = prod(dims)
    lind = LinearIndices(sys.ih, dims)

    I = Int[]
    J = Int[]
    V = Float64[]

    predsize = 2 * Ntot * length(Catalyst.equations(sys.rs))

    sizehint!(I, predsize)
    sizehint!(J, predsize)
    sizehint!(V, predsize)

    S::Matrix{Int64} = netstoichmat(sys.rs)
    for (i, rf) in enumerate(sys.rfs)
        for (idx_cin, idx_cout) in pairedindices(sys.ih, dims, CartesianIndex(S[:,i]...))
            idx_lin = lind[idx_cin]
            idx_lout = lind[idx_cout]
            push!(I, lind[idx_cout])
            push!(J, lind[idx_cin])

            rate = rf(idx_cin, 0, ps...)
            push!(V, rate)

            push!(I, lind[idx_cin])
            push!(J, lind[idx_cin])
            push!(V, -rate)
        end
    end

    sparse(I, J, V)
end


"""
    Base.convert(::Type{SparseMatrixCSC}, sys::FSPSystem, dims::NTuple, ps, t::Real)

Convert the reaction system into a sparse matrix defining the right-hand side of the
Chemical Master Equation. `dims` is a tuple denoting the dimensions of the FSP and
`ps` is the tuple of parameters. The sparse matrix works on the flattened version
of the state obtained using `vec`.
"""
function Base.convert(::Type{SparseMatrixCSC}, sys::FSPSystem, dims::NTuple, ps, t::Real)
    create_sparsematrix(sys, dims, ps, t)
end

"""
    Base.convert(::Type{SparseMatrixCSC}, sys::FSPSystem, dims::NTuple, ps, ::SteadyState)

Convert the reaction system into a sparse matrix defining the right-hand side of the
Chemical Master Equation, steady-state version.
"""
function Base.convert(::Type{SparseMatrixCSC}, sys::FSPSystem, dims::NTuple, ps, ::SteadyState)
    create_sparsematrix_ss(sys, dims, ps)
end
