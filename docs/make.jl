using Documenter
using FiniteStateProjection
using SparseArrays
using FiniteStateProjection: AbstractIndexHandler

makedocs(sitename="FiniteStateProjection.jl", 
         modules=[FiniteStateProjection],
         format = Documenter.HTML(prettyurls = false),
         pages = [
             "Home" => "index.md",
             "Main API" => "mainapi.md",
             "Index Handlers" => "indexhandlers.md",
             "Matrix Conversions" => "matrix.md",
             "Internal API" => "internal.md",
             "Examples" => "examples.md"
         ])
 

deploydocs(repo = "github.com/kaandocal/FiniteStateProjection.jl.git", devbranch="main")
