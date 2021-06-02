using Documenter, FiniteStateProjection
using FiniteStateProjection: AbstractIndexHandler

makedocs(sitename="FiniteStateProjection.jl", 
         format = Documenter.HTML(prettyurls = false),
         pages = [
             "Home" => "index.md",
             "Main API" => "mainapi.md",
             "Index Handlers" => "indexhandlers.md",
             "Internal API" => "internal.md",
             "Examples" => "examples.md"
         ])
 

deploydocs(repo = "github.com/kaandocal/FiniteStateProjection.jl.git", devbranch="main")
