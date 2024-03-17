using DifferentialEquations, GroupNumbers, Statistics

# define the differential equation
f(x, r, t) = r * x + x^3 - x^5


# r < -1 is uninteresting as it gives only one equilibrium
# r > 1 results in 3 equilibrium points which properties are all the same
# Therefore only -1 < r < 1 is considered
r_rng = range(-1, 1, length=100)

# x(0) ranging from -1.5 to 1.5 will cover all interesting cases
# which follows from the plots of f(x, -1) and f(x, 1)
x0_rng = range(-1.5, 1.5, length=100)

# uniform norm of f(x, r) on the interesting set is less then 10
# thus choosing the time step around 10^-4 will be sufficient to approximate the solution
# which fastest decay is slower then exp(-10 * x) (characteristic time scale 0.1)
sup_f = abs(f(1.5, -1, 0)) # 5.71875
dt = 1e-2

function get_sol(f, r, x0, dt, tmax)
    prob = ODEProblem(f, x0, (0, tmax), r)
    sol = solve(prob, RK4(), dt=dt, adaptive=false, maxiters = 1e7)
    return sol
end

function final_state(f, r, x0, dt, tmax)
    prob = ODEProblem(f, x0, (0, tmax), r)
    sol = solve(prob, RK4(), dt=dt, adaptive=false, save_everystep=false, maxiters = 1e7)
    return sol.u[end]
end

function find_steady_states(x0_rng, xf_rng, atol)
    unst_index = findall([!isapprox(xf_rng[i], xf_rng[i+1], atol=atol) for i in 1:length(x0_rng)-1])
    unstable = [0.5 * (x0_rng[i] + x0_rng[i+1]) for i in unst_index]

    unst_index = [1, unst_index..., length(xf_rng)]
    stable = [mean(xf_rng[unst_index[i]:unst_index[i+1]]) for i in 1:length(unst_index)-1]
    return stable, unstable
end

function find_stable_points(r, f, x0_rng, dt, tmax, atol)
    # solves the equation to get final positions after a long time
    xfs = [final_state(f, r, x0, dt, tmax) for x0 in x0_rng]
    # groups the solutions and calculates a mean in each bucket
    stable = [mean(x) for x in groupby_numbers(xfs; atol=atol)]
    return sort(stable)
end


# sol = get_sol(f, 0.0, 0.1, dt, 1e3)
# all_stable = [find_stable_points(r, f, x0_rng, dt, 1e2, 1e-2) for r in r_rng]


using CairoMakie

fig = Figure()
ax = Axis(fig[1,1])

for x in eachindex(all_stable)
    scatter!(ax, fill(r_rng[x], length(all_stable[x])), all_stable[x], color=:black, markersize=3)
end


