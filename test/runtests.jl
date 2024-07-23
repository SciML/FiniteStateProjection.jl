using SafeTestsets

@time begin
    @safetestset "Telegraph" begin include("telegraph.jl") end
    @safetestset "FeedbackLoop" begin include("feedbackloop.jl") end
    @safetestset "BirthDeath2D" begin include("birthdeath2D.jl") end
end
