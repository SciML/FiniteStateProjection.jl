using MomentClosure, FiniteStateProjection, Plots, Catalyst, OrdinaryDiffEq, Sundials

function get_raw_moments_FSP(sol::ODESolution, order::Int)

    state_space = size(sol)[1:end-1]
    N = length(state_space)
    no_t_pts = length(sol.u)

    iter_μ = MomentClosure.construct_iter_all(N, order)[2:end]
    μ = Dict([iter => Array{Float64}(undef, no_t_pts) for iter in iter_μ])
    iter_state = Iterators.product((0:i-1 for i in state_space)...)

    for t_pt in 1:no_t_pts
        tslice = sol[t_pt]
        for μ_ind in iter_μ
            μ[μ_ind][t_pt] = sum(prod(n_state.^μ_ind) * tslice[(n_state.+1)...] for n_state in iter_state)
        end
    end
    return μ

end

rn = @reaction_network begin
    (c₁/Ω), 2X → 3X
    (c₂/Ω^2), 3X → 2X
    (c₃*Ω), 0 → X
    (c₄), X → 0
end c₁ c₂ c₃ c₄ Ω

sys = FSPSystem(rn, combinatoric_ratelaw=false)
#[c₁, c₂, c₃, c₄, Ω]
ps = [1., 0.25, 0.25, 1., 100.]

tspan = (0., 2000.)
ic = [70]
dt = 1.0
ts = tspan[1]:dt:tspan[2]

u0 = zeros(2000)
u0[(ic.+1)...] = 1.0

prob = ODEProblem(sys, u0, tspan, ps)
sol = solve(prob, CVODE_BDF(), saveat=dt)
moments_FSP =@time get_raw_moments_FSP(sol, 2)
μ_FSP = moments_FSP[(1,)]
Σ_FSP = moments_FSP[(2,)] .- moments_FSP[(1,)].^2;
plot(μ_FSP, Σ_FSP, legend=:bottomright, xlabel="μ", ylabel="Σ", fmt="svg")

# --------------------------------------------------------------

@parameters ρ σ_on σ_off d
rn = @reaction_network begin
    ρ, G_on --> G_on + M
    (σ_on, σ_off), G_off <--> G_on
    d, M --> 0
end ρ σ_on σ_off d

# This automatically reduces the dimensionality of the
# network by exploiting conservation laws
ih = ReducingIndexHandler(rn)
sys = FSPSystem(rn, ih)

# There is one conserved quantity: G_on + G_off
cons = conservedquantities([1, 0, 0], sys)

# Parameters for our system
ps = [ 15.0, 0.25, 0.15, 1.0 ]

# In the reduced model, G_off = 1 - G_on does not have to be tracked
u0 = zeros(2, 50)
u0[1,1] = 1.0

prob = ODEProblem(sys, u0, (0, 10.0), (ps, cons))
sol = solve(prob, CVODE_BDF(), saveat=0.1)
sol

plot(sol.u[end][1,:] + sol.u[end][2,:])
