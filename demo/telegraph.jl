using FSP
using DifferentialEquations
using PyPlot

@parameters r1 r2 r3 r4

rs = @reaction_network begin
    r1, G_on --> G_on + M
    (r2, r3), G_on <--> G_off
    r4, M --> 0
end r1 r2 r3 r4

sys = FSPSystem(rs)

# There is one conserved quantity: G_on + G_off
cons = get_conserved_quantities([1,0,0], sys)

# Parameters for our system
ps = [ 15.0, 0.25, 0.15, 1.0 ]

# Initial values
# Since G_on + G_off = const. we do not have to model the two
# separately. Use reduced_species(sys) to get the list of 
# species we actually have to model:
# 
# julia> reduced_species(sys)
# 2-element Vector{Int64}:
#  1
#  2
#
u0 = zeros(2, 50)
u0[1,1] = 1.0

prob = build_ode_prob(sys, u0, 10.0, ps, cons)

sol = solve(prob, Vern7(), dense=false, save_everystep=false, atol=1e-6)

plt.suptitle(L"Distribution at $t = 10$")
plt.bar(0:49, sol.u[end][1,:] + sol.u[end][2,:], width=1);
plt.xlabel("# of mRNA")
plt.ylabel("Probability")
plt.xlim(-0.5, 49.5)

plt.show()
