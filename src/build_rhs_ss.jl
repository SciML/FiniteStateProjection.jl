struct SteadyState end

"""
    build_rhs_singlepass_ss(sys::FSPSystem)

Return code for the RHS function in a single pass, for use with
steady state problems. Transitions out of the truncated state space
are ignored.

See also: [`build_rhs_ss`](@ref)
"""
function build_rhs_singlepass_ss(sys::FSPSystem)
    isempty(sys.rfs) && return quote end

    S = netstoichmat(sys.rs)
    ret = Expr(:block, :(fill!(du, 0)))

    for (i, rf) in enumerate(sys.rfs)
        ex = quote
            for (idx_in, idx_out) in pairedindices($(sys.ih), u, $(CartesianIndex(S[:,i]...)))
                rate = u[idx_in] * $(rf.body)
                du[idx_in] -= rate
                du[idx_out] += rate
            end
        end

        append!(ret.args, ex.args)
    end

    return ret
end

##

function build_rhs_ex_ss(sys::FSPSystem; striplines::Bool=false)
    header = build_rhs_header(sys)
    single_pass = build_rhs_singlepass_ss(sys)

    body = Expr(:block, header, single_pass)

    ex = :((du, u, p, t) -> $(body))

    striplines && (ex = MacroTools.striplines(ex))

    ex = ex |> MacroTools.flatten |> MacroTools.prettify

    ex
end

"""
    build_rhs_ss(sys::FSPSystem)

Builds the function `f(du,u,p,t)` that defines the right-hand side of the CME
for use with `SteadyStateProblem`s.
"""
function build_rhs_ss(sys::FSPSystem)
    @RuntimeGeneratedFunction(build_rhs_ex_ss(sys; striplines=false))
end

##

"""
    Base.convert(::Type{ODEFunction}, sys::FSPSystem, ::SteadyState)

Return an `ODEFunction` defining the right-hand side of the CME, for use
with `SteadyStateProblem`s.
"""
Base.convert(::Type{ODEFunction}, sys::FSPSystem, ::SteadyState) = ODEFunction{true}(build_rhs_ss(sys))

"""
    Base.convert(::Type{SteadyStateProblem}, sys::FSPSystem, u0[, p])

Return a `SteadyStateProblem` for use in `DifferentialEquations.
"""
function Base.convert(::Type{SteadyStateProblem}, sys::FSPSystem, u0, p=NullParameters())
    SteadyStateProblem(convert(ODEFunction, sys, SteadyState()), u0, p)
end
