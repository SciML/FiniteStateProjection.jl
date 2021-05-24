using Documenter, FiniteStateProjection

makedocs(sitename="FiniteStateProjection.jl", format = Documenter.HTML(prettyurls = false))

deploydocs(repo = "github.com/kaandocal/FiniteStateProjection.jl.git", devbranch="main")
