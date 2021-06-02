compile_ratefunc(ex, params) = @RuntimeGeneratedFunction(Expr(:->, Expr(:tuple, :idx_in, params...), Expr(:block, ex)))

function create_sparsematrix(ih::AbstractIndexHandler, sys::FSPSystem, dims::NTuple, ps;
                             combinatoric_ratelaw::Bool=true)
    Ntot = prod(dims)
    lind = LinearIndices(dims)
        
    paramsyms::Vector{Symbol} = Symbol.(Catalyst.params(sys.rs))
    ratefuncs::Vector{Function} = map(ex -> compile_ratefunc(ex, paramsyms), 
                                      build_ratefuncs(ih, sys; state_sym=:idx_in, combinatoric_ratelaw))

    I = Int[]
    J = Int[]
    V = Float64[]

    predsize = Ntot * (length(Catalyst.get_eqs(sys.rs)) + 1)

    sizehint!(I, predsize)
    sizehint!(J, predsize)
    sizehint!(V, predsize)

    V = zeros(predsize)

    for idx_cart in singleindices(ih, dims)
        idx_lin = lind[idx_cart]
        push!(I, idx_lin)
        push!(J, idx_lin)

        k = length(I)
        for rf in ratefuncs
            V[k] -= rf(idx_cart, ps...)
        end
    end

    S::Matrix{Int64} = netstoichmat(sys.rs)
    for (i, reaction) in enumerate(Catalyst.get_eqs(sys.rs))
        for (idx_cin, idx_cout) in pairedindices(ih, dims, CartesianIndex(S[i,:]...))
            idx_lin = lind[idx_cin]
            idx_lout = lind[idx_cout]
            push!(I, lind[idx_cout])
            push!(J, lind[idx_cin])

            k = length(I)
            V[k] = ratefuncs[i](idx_cin, ps...)
        end
    end

    resize!(V, length(I))
    sparse(I, J, V)
end

function overwrite_sparsematrix!(out::SparseMatrixCSC, ih::AbstractIndexHandler, 
                                 sys::FSPSystem, dims::NTuple, ps;
                                 combinatoric_ratelaw::Bool=true)
    Ntot = prod(dims)
    lind = LinearIndices(dims)
        
    fill!(out.nzval, 0.0)
    
    paramsyms::Vector{Symbol} = Symbol.(Catalyst.params(sys.rs))
    ratefuncs::Vector{Function} = map(ex -> compile_ratefunc(ex, paramsyms), 
                                      build_ratefuncs(ih, sys; state_sym=:idx_in, combinatoric_ratelaw))
    
    for idx_cart in singleindices(ih, dims)
        idx_lin = lind[idx_cart]
        val = 0.0
        for rf in ratefuncs
            val += rf(idx_cart, ps...)
        end
        
        out[idx_lin, idx_lin] = -val
    end

    S::Matrix{Int64} = netstoichmat(sys.rs)
    for (i, reaction) in enumerate(Catalyst.get_eqs(sys.rs))
        for (idx_cin, idx_cout) in pairedindices(ih, dims, CartesianIndex(S[i,:]...))
            idx_lin = lind[idx_cin]
            idx_lout = lind[idx_cout]
            out[idx_lout, idx_lin] += ratefuncs[i](idx_cin, ps...)
        end
    end

    out
end

function Base.convert(::Type{SparseMatrixCSC}, ih::AbstractIndexHandler, sys::FSPSystem, dims::NTuple, ps; combinatoric_ratelaw::Bool=true)
    return create_sparsematrix(ih, sys, dims, ps; combinatoric_ratelaw)
end

function Base.fill!(out::SparseMatrixCSC, ih::AbstractIndexHandler, sys::FSPSystem, dims::NTuple, ps; combinatoric_ratelaw::Bool=true)
    overwrite_sparsematrix!(out, ih, sys, dims::NTuple, ps; combinatoric_ratelaw)
end
