"""
    unpackparams(sys::FSPSystem, psym::Symbol)

Returns code unpacking the parameters of the system from the symbol
`psym` in the form `(p1, p2, ...) = psym`. This should be called in
all overloads of [`build_rhs_header`](@ref). It is assumed that
the variable `psym` is an `AbstractVector`.

See also: [`build_rhs_header`](@ref), [`build_rhs`](@ref)
"""
function unpackparams(sys::FSPSystem, psym::Symbol)::Expr
    param_names = Expr(:tuple, map(par -> par.name, params(sys.rs))...)
     
    quote 
        $(param_names) = $(psym)
    end
end

"""
    build_rhs_header(idxhandler::AbstractIndexHandler, sys::FSPSystem)::Expr

Return initialisation code for the RHS function, unpacking the parameters
`p` supplied by `DifferentialEquations`. The default implementation
just unpacks parameters from `p`.

See also: [`unpackparams`](@ref), [`build_rhs`](@ref)
"""
function build_rhs_header(sys::FSPSystem)::Expr
    quote 
        ps = p
        $(unpackparams(sys, :ps))
    end
end

##

"""
    build_rhs_firstpass(sys::FSPSystem, rfs)::Expr

Return code for the first pass of the RHS function, for the time-dependent
FSP. Goes through all reactions and computes the negative part of the CME
(probability flowing out of states). This is a simple array traversal and
can be done in one go for all reactions.

See also: [`build_rhs`](@ref)
"""
function build_rhs_firstpass(sys::FSPSystem)
    isempty(sys.rfs) && return quote end

    first_line = :(du[idx_in] = -u[idx_in] * $(sys.rfs[1].body))
    other_lines = (:(du[idx_in] -= u[idx_in] * $(rf.body)) for rf in sys.rfs[2:end])

    quote
        for idx_in in singleindices($(sys.ih), u)
            $first_line
            $(other_lines...)
        end
    end
end

##

"""
    build_rhs_secondpass(sys::FSPSystem)::Expr

Return code for the second pass of the RHS function. Goes through
all reactions and computes the positive part of the CME (probability
flowing into states). This requires accessing `du` and `u` at different
locations depending on the net stoichiometries. In order to reduce
random memory access reactions are processed one by one.

See also: [`build_rhs`](@ref)
"""
function build_rhs_secondpass(sys::FSPSystem)::Expr
    isempty(sys.rfs) && return quote end

    S = netstoichmat(sys.rs)
    ret = Expr(:block)

    for (i, rf) in enumerate(sys.rfs)
        ex = quote
            for (idx_in, idx_out) in pairedindices($(sys.ih), u, $(CartesianIndex(S[i,:]...)))
                du[idx_out] += u[idx_in] * $(rf.body)
            end
        end

        append!(ret.args, ex.args)
    end

    return ret
end

##

function build_rhs_ex(sys::FSPSystem; striplines::Bool=false)
    header = build_rhs_header(sys)
    first_pass = build_rhs_firstpass(sys)
    second_pass = build_rhs_secondpass(sys)

    body = Expr(:block, header, first_pass, second_pass)

    ex = :((du, u, p, t) -> $(body))

    striplines && (ex = MacroTools.striplines(ex))

    ex = ex |> MacroTools.flatten |> MacroTools.prettify

    ex
end

"""
    build_rhs(sys::FSPSystem)

Builds the function `f(du,u,p,t)` that defines the right-hand side of the CME
for use with DifferentialEquations.jl.
"""
function build_rhs(sys::FSPSystem)
    @RuntimeGeneratedFunction(build_rhs_ex(sys; striplines=false))
end

##

"""
    convert(::Type{ODEFunction}, idxhandler::AbstractIndexHandler, sys::FSPSystem;
            combinatoric_ratelaw::Bool=true)

Return an `ODEFunction` defining the right-hand side of the CME.

Creates an `ODEFunction` for use with `DifferentialEquations`. This is where most of
the work in the package happens; for best performance it is suggested to build an `ODEFunction`
once for a given reaction system and reuse it instead of directly converting
a reaction system to an `ODEProblem` (which implicitly calls this function).
"""
function Base.convert(::Type{ODEFunction}, sys::FSPSystem;
                      combinatoric_ratelaw::Bool=true)::ODEFunction
    rhs = build_rhs(sys)
    ODEFunction{true}(rhs)
end

"""
    convert(::Type{ODEProblem}, idxhandler::AbstractIndexHandler, sys::FSPSystem, u0, tmax, p)

Return an `ODEProblem` for use in `DifferentialEquations. This function implicitly
calls `convert(ODEFunction, indexhandler, sys)`. It is usually more efficient to
create an `ODEFunction` first and then use that to create `ODEProblem`s.
"""
function Base.convert(::Type{ODEProblem}, sys::FSPSystem, u0, tint, p;
                      combinatoric_ratelaw::Bool=true)::ODEProblem
     ODEProblem(convert(ODEFunction, sys), u0, tint, p; combinatoric_ratelaw)
end
