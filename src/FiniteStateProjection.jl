module FiniteStateProjection

using Reexport
@reexport using Catalyst

using ModelingToolkit
using MacroTools
using SparseArrays

import Base: LinearIndices, vec

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

export FSPSystem, NaiveIndexHandler, ReducingIndexHandler,
       conservedquantities, conservationlaws,
       reducedspecies, elidedspecies

# Check for Catalyst v9.0
if isdefined(Catalyst, :reactioncomplexmap)
    netstoichmat(rn) = Catalyst.netstoichmat(rn)'
end

abstract type AbstractIndexHandler end

include("fspsystem.jl")
include("indexhandlers.jl")
include("build_rhs.jl")
include("matrix.jl")

end
