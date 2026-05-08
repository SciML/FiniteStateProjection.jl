module FiniteStateProjection

using Reexport
@reexport using Catalyst

using MacroTools
using SparseArrays
using DiffEqBase

import Base: LinearIndices, vec

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

export FSPSystem, DefaultIndexHandler, SteadyState

# Julia 1.12 removed the `mt` field on `Base.MethodList`, breaking
# `MacroTools.prettify`'s `unresolve` step (FluxML/MacroTools.jl#216).
# Replicate `MacroTools.prettify` locally with a 1.12-safe `unresolve`
# so generated CME RHS expressions can still be tidied for display.
_unresolve1(x) = x
@static if VERSION >= v"1.12-"
    _unresolve1(f::Function) = nameof(f)
else
    _unresolve1(f::Function) = methods(f).mt.name
end
_unresolve(ex) = MacroTools.prewalk(_unresolve1, ex)
function _prettify(ex; lines = false, alias = true)
    ex = lines ? ex : MacroTools.striplines(ex)
    ex = MacroTools.flatten(ex)
    ex = _unresolve(ex)
    ex = MacroTools.resyntax(ex)
    return alias ? MacroTools.alias_gensyms(ex) : ex
end

abstract type AbstractIndexHandler end

include("fspsystem.jl")
include("indexhandlers.jl")
include("build_rhs.jl")
include("build_rhs_ss.jl")
include("matrix.jl")

end
