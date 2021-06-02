compile_ratefunc(ex, params) = @RuntimeGeneratedFunction(Expr(:->, Expr(:tuple, :idx_in, params...), Expr(:block, ex)))

function create_sparsematrix(sys::FSPSystem, ih::AbstractIndexHandler, dims::NTuple, ps;
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

function Base.convert(::Type{SparseMatrixCSC}, sys::FSPSystem, ih::AbstractIndexHandler, 
                      dims::NTuple, ps; combinatoric_ratelaw::Bool=true)
    return create_sparsematrix(sys, ih, dims, ps; combinatoric_ratelaw)
end
