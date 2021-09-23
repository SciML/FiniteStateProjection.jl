using SafeTestsets

@time begin
    @safetestset "BirthDeath2D" begin include("birthdeath2D.jl") end
    @safetestset "FeedbackLoop" begin include("feedbackloop.jl") end
end
