using FiniteStateProjection
using DifferentialEquations
using PyPlot

##

rs = @reaction_network begin
    r1, 0 --> A
    r2, A --> 0
end

##

sys = FSPSystem(rs)

# Parameters for our system
ps = [ :r1 => 10.0, :r2 => 1.0 ]

# Initial values
u0 = zeros(50)
u0[1] = 1.0

##

prob = ODEProblem(sys, u0, 10.0, ps)

##

sol = solve(prob, Vern7(), dense=false, save_everystep=false, abstol=1e-6)

##

plt.suptitle(L"Distribution at $t = 10$")
plt.bar(0:49, sol.u[end], width=1);
plt.xlabel("# of Molecules")
plt.ylabel("Probability")
plt.xlim(-0.5, 49.5)

plt.show()
