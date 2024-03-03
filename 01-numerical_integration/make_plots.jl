using CairoMakie

set_my_theme!()

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
ax = Axis(fig[1,1], yscale = Makie.Symlog10(-1e-3, 3.0))
ax.limits=(nothing, (-1e4, 3))
ax.yticks = [-1e4, -1e3, -1e2, -10, 0, 2, 3]
ax.xlabel = "t"
ax.ylabel = "x"

lines!(ax, sol_f.t, getindex.(sol_f.u, 1), label=L"Runge–Kutta, $dt=\frac{1}{516}$")
lines!(ax, sol_g.t[1:30], getindex.(sol_g.u, 1)[1:30], label=L"Runge–Kutta, $dt=\frac{1}{256}$")
lines!(ax, 0..0.125, t -> getindex(exp(A*t) * X0, 1), label="exact", color=:black, linestyle=:dash)

axislegend(ax, position=:rb)
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
ax = Axis(fig[1,1], yscale = log10)
ax.xlabel = "t"
ax.ylabel = "error (euclidian metric)"

scatter!(ax, sol_a.t[2:end], norm.([sol_a.u[i] - exact_sol(A, sol_a.t[i], X0) for i in 2:length(sol_a.t)]),
         label=L"Euler, $dt=\frac{1}{2048}$")
scatter!(ax, sol_b.t[2:end], norm.([sol_b.u[i] - exact_sol(A, sol_b.t[i], X0) for i in 2:length(sol_b.t)]),
         label=L"Euler, $dt=\frac{1}{512}$")

scatter!(ax, sol_f.t[2:end], norm.([sol_f.u[i] - exact_sol(A, sol_f.t[i], X0) for i in 2:length(sol_f.t)]),
         label=L"RK4, $dt=\frac{1}{512}$")

scatter!(ax, time_d[2:end], norm.([sol_d[i, :] - exact_sol(A, time_d[i], X0) for i in 2:length(time_d)]),
         label=L"Back. Euler, $dt=\frac{1}{516}$")
scatter!(ax, time_e[2:end], norm.([sol_e[i, :] - exact_sol(A, time_e[i], X0) for i in 2:length(time_e)]),
         label=L"Back. Euler, $dt=\frac{1}{256}$")

fig[1, 2] = Legend(fig, ax, "",)
save("method-errors.pdf", fig)

