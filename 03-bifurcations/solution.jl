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

function find_steady_states(r, f, x0_rng, dt, tmax, atol)
    # solves the equation to get final positions after a long time
    xf_rng = [final_state(f, r, x0, dt, tmax) for x0 in x0_rng]

    # finds neighbouring solutions that converge to different attractors
    unst_index = findall([!isapprox(xf_rng[i], xf_rng[i+1], atol=atol) for i in 1:length(x0_rng)-1])
    unstable = [0.5 * (x0_rng[i] + x0_rng[i+1]) for i in unst_index]

    border_index = [0, unst_index..., length(xf_rng)]
    stable = [mean(xf_rng[border_index[i]+1:border_index[i+1]]) for i in 1:length(border_index)-1]
    return stable, unstable
end


# sol = get_sol(f, 0.0, 0.1, dt, 1e3)
all_ss = [find_steady_states(r, f, x0_rng, dt, 2e2, 1e-2) for r in r_rng]


using CairoMakie

fig = Figure()
ax = Axis(fig[1,1])

for i in eachindex(all_ss)
    stable, unstable = all_ss[i]
    scatter!(ax, fill(r_rng[i], length(stable)), stable, color=Makie.wong_colors()[5])
    scatter!(ax, fill(r_rng[i], length(unstable)), unstable,
             color=Makie.wong_colors()[2], marker=:utriangle)
end

lines!(ax, -1..0, x -> 0, color=:black)
lines!(ax, 0..1, x -> 0, color=:red, linestyle=:dash)

lines!(ax, -0.25..0, r -> sqrt(0.5 - sqrt(0.25 + r)), color=:red, linestyle=:dash)
lines!(ax, -0.25..0, r -> -sqrt(0.5 - sqrt(0.25 + r)), color=:red, linestyle=:dash)

lines!(ax, -0.25..1, r -> sqrt(0.5 + sqrt(0.25 + r)), color=:black)
lines!(ax, -0.25..1, r -> -sqrt(0.5 + sqrt(0.25 + r)), color=:black)

ax.xlabel = L"$r$ - parameter"
ax.ylabel = L"$x$ - stable point position"

elem_1 = [LineElement(color = :black, linestyle = nothing)]
elem_2 = [LineElement(color = :red, linestyle = :dash)]
elem_3 = [MarkerElement(color = Makie.wong_colors()[5], marker = :circle)]
elem_4 = [MarkerElement(color = Makie.wong_colors()[2], marker = :utriangle)]
axislegend(ax, [elem_1, elem_2, elem_3, elem_4],
           ["theoretical stable", "theoretical unstable", "numerical stable", "numerical unstable"],
          position=:lt)

save("bifurc-diagram.pdf", fig)

