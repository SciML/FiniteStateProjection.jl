using SafeTestsets, Pkg

const GROUP = get(ENV, "GROUP", "All")

@time begin
    if GROUP == "All" || GROUP == "Core"
        @safetestset "Telegraph" begin
            include("telegraph.jl")
        end
        @safetestset "FeedbackLoop" begin
            include("feedbackloop.jl")
        end
        @safetestset "BirthDeath2D" begin
            include("birthdeath2D.jl")
        end
    end

    if GROUP == "QA"
        Pkg.activate(joinpath(@__DIR__, "qa"))
        Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
        Pkg.instantiate()
        include("qa/qa.jl")
    end
end
