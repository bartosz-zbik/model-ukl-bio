# imports https://docs.sciml.ai/DiffEqDocs/stable/
using DifferentialEquations
using CairoMakie # for plotting
using LinearAlgebra


##################
# Explicit methods
##################

# define the differential equation
# X - current system state
# p - parameter vector
# t - time
f(X, p, t) = [1023.0 2023.0; (-1024.0) (-2024.0)] * X

# initial conditions
X0 = [1.0, 0.0]
# time domain on which to solve the problem
tspan = (0.0, 0.125)

# specifies the "problem", it's specific for this ecosystem
prob = ODEProblem(f, X0, tspan)

# Solve using different algorithms: https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/
# Euler(), RK4() - specifies the algorithm
sol_a = solve(prob, Euler(), dt=1/2048)
sol_b = solve(prob, Euler(), dt=1/512)
sol_c = solve(prob, Euler(), dt=1/256)

sol_f = solve(prob, RK4(), dt=1/512, adaptive=false)
sol_g = solve(prob, RK4(), dt=1/256, adaptive=false)




##################
# Non Explicit methods
##################

# Non explicit Euler method for a linear equation
function BackwardEulerStep(x, A, dt)
    # x_n+1 = x_n + dt * A x_n+1
    # (1 - A dt) x_n+1 = x_n
    # x_n+1 = (1 - A dt)^-1 x_n
    return ([1 0; 0 1] - A * dt)^-1 * x
end

A = [1023.0 2023.0; (-1024.0) (-2024.0)]

nsteps_d = floor(Int, 1 + tspan[2] / (1/512)) 
sol_d = zeros(1 + nsteps_d, 2) 
sol_d[1, :] .= X0
for i in 2:nsteps_d+1
    sol_d[i, :] .= BackwardEulerStep(sol_d[i-1, :], A, 1/512)
end
time_d = collect(0:1/512:0.125 + 1/512)


nsteps_e = floor(Int, 1 + tspan[2] / (1/256)) 
sol_e = zeros(1 + nsteps_e, 2) 
sol_e[1, :] .= X0
for i in 2:nsteps_e+1
    sol_e[i, :] .= BackwardEulerStep(sol_e[i-1, :], A, 1/256)
end
time_e = collect(0:1/256:0.125 + 1/256)


function exact_sol(A, t, X0)
    return exp(A*t) * X0
end

# Euler methods
fig = Figure()
ax = Axis(fig[1,1], yscale = Makie.pseudolog10)
ax.limits=(nothing, (0, 100))
ax.xlabel = "t"
ax.ylabel = "x"

lines!(ax, sol_a.t, getindex.(sol_a.u, 1), label=L"Euler, dt=\frac{1}{2048}")
lines!(ax, sol_b.t, getindex.(sol_b.u, 1), label=L"Euler, dt=\frac{1}{516}")
lines!(ax, sol_c.t, getindex.(sol_c.u, 1), label=L"Euler, dt=\frac{1}{256}")
lines!(ax, 0..0.125, t -> getindex(exp(A*t) * X0, 1), label="exact", color=:black, linestyle=:dash)

axislegend(ax)
save("euler-methods.pdf", fig)



# Runge–Kutta methods
fig = Figure()
ax = Axis(fig[1,1], yscale = Makie.Symlog10(10.0))
ax.limits=(nothing, (-1e4, 1e4))
ax.yticks = [-1e4, -1e3, -1e2, -10, -5, 0, 5, 10, 1e2, 1e3, 1e4]
ax.xlabel = "t"
ax.ylabel = "x"

lines!(ax, sol_f.t, getindex.(sol_f.u, 1), label=L"Runge–Kutta, $dt=\frac{1}{516}$")
lines!(ax, sol_g.t[1:30], getindex.(sol_g.u, 1)[1:30], label=L"Runge–Kutta, $dt=\frac{1}{256}$")
lines!(ax, 0..0.125, t -> getindex(exp(A*t) * X0, 1), label="exact", color=:black, linestyle=:dash)

axislegend(ax)
save("RK-methods.pdf", fig)



# Backward Euler methods
fig = Figure()
ax = Axis(fig[1,1])
ax.xlabel = "t"
ax.ylabel = "x"

lines!(ax, time_d, sol_d[:, 1], label=L"Backward Euler, $dt=\frac{1}{516}$")
lines!(ax, time_e, sol_e[:, 1], label=L"Backward Euler, $dt=\frac{1}{256}$")
lines!(ax, 0..0.125, t -> getindex(exp(A*t) * X0, 1), label="exact", color=:black, linestyle=:dash)

axislegend(ax, position=:rb)
save("BE-methods.pdf", fig)



# Errors
fig = Figure()
ax = Axis(fig[1,1], yscale = Makie.pseudolog10)
ax.limits = (nothing, (0, 0.3))
ax.yticks = vcat([0], [10.0^i for i in -2:-1], [0.3])
ax.xlabel = "t"
ax.ylabel = "error (euclidian metric)"

scatter!(ax, sol_a.t, norm.([sol_a.u[i] - exact_sol(A, sol_a.t[i], X0) for i in eachindex(sol_a.t)]),
         label=L"Euler, $dt=\frac{1}{2048}$")
scatter!(ax, sol_b.t, norm.([sol_b.u[i] - exact_sol(A, sol_b.t[i], X0) for i in eachindex(sol_b.t)]),
         label=L"Euler, $dt=\frac{1}{512}$")

scatter!(ax, sol_f.t, norm.([sol_f.u[i] - exact_sol(A, sol_f.t[i], X0) for i in eachindex(sol_f.t)]),
         label=L"RK4, $dt=\frac{1}{512}$")

scatter!(ax, time_d, norm.([sol_d[i, :] - exact_sol(A, time_d[i], X0) for i in eachindex(time_d)]),
         label=L"Backward Euler, $dt=\frac{1}{516}$")
scatter!(ax, time_e, norm.([sol_e[i, :] - exact_sol(A, time_e[i], X0) for i in eachindex(time_e)]),
         label=L"Backward Euler, $dt=\frac{1}{516}$")


axislegend(ax)
save("method-errors.pdf", fig)




