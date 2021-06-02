module FiniteStateProjection

using Reexport
@reexport using Catalyst

using ModelingToolkit
using MacroTools
using SparseArrays

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

export FSPSystem, AbstractIndexHandler, NaiveIndexHandler, DefaultIndexHandler,
       conservedquantities, conservationlaws,
       pairedindices, singleindices,
       reducedspecies, elidedspecies

include("fspsystem.jl")
include("indexhandlers.jl")
include("build_rhs.jl")
include("matrix.jl")

end
