"""
    struct SteadyState end

Marker type used to select the steady-state formulation of the Chemical
Master Equation. Passing a `SteadyState()` instance to the conversion and
matrix-building methods requests the right-hand side in which transitions
out of the truncated state space are dropped (their propensity is set to
zero), as required for solving for the stationary distribution rather than
the time-dependent one.

It is used as a dispatch tag by
[`Base.convert(::Type{ODEFunction}, ::FSPSystem, ::SteadyState)`](@ref),
[`SparseArrays.SparseMatrixCSC(::FSPSystem, ::NTuple, ps, ::SteadyState)`](@ref)
and [`SciMLBase.SteadyStateProblem(::FSPSystem, u0, p)`](@ref).

# Examples
```julia
julia > SteadyState()
SteadyState()
```
"""
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
            for (idx_in, idx_out) in
                pairedindices($(sys.ih), u, $(CartesianIndex(S[:, i]...)))
                rate = u[idx_in] * $(rf.expression)
                du[idx_in] -= rate
                du[idx_out] += rate
            end
        end

        append!(ret.args, ex.args)
    end

    return ret
end

##

function build_rhs_ex_ss(sys::FSPSystem; striplines::Bool = false)
    header = build_rhs_header(sys)
    single_pass = build_rhs_singlepass_ss(sys)

    body = Expr(:block, header, single_pass)

    ex = :((du, u, p, t) -> $(body))

    striplines && (ex = _striplines(ex))

    ex = _prettify(_flatten(ex))

    return ex
end

"""
    build_rhs_ss(sys::FSPSystem)

Builds the function `f(du,u,p,t)` that defines the right-hand side of the CME
for use with `SteadyStateProblem`s.
"""
function build_rhs_ss(sys::FSPSystem)
    return @RuntimeGeneratedFunction(build_rhs_ex_ss(sys; striplines = false))
end

##

"""
    Base.convert(::Type{ODEFunction}, sys::FSPSystem, ::SteadyState)

Return an `ODEFunction` defining the right-hand side of the CME, for use
with `SteadyStateProblem`s.
"""
function Base.convert(::Type{ODEFunction}, sys::FSPSystem, ::SteadyState)
    return ODEFunction{true, SciMLBase.FullSpecialize}(build_rhs_ss(sys))
end

"""
    SciMLBase.SteadyStateProblem(sys::FSPSystem, u0[, p])

Return a `SteadyStateProblem` for use in `DifferentialEquations.
"""
function SciMLBase.SteadyStateProblem(sys::FSPSystem, u0, pmap = SciMLBase.NullParameters())
    return SteadyStateProblem(convert(ODEFunction, sys, SteadyState()), u0, pmap_to_p(sys, pmap))
end
