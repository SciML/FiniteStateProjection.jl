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

abstract type AbstractIndexHandler end

include("fspsystem.jl")
include("indexhandlers.jl")
include("build_rhs.jl")
include("build_rhs_ss.jl")
include("matrix.jl")

end
