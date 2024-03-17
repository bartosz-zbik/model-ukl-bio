using CairoMakie
set_theme!(theme_latexfonts())

fig = Figure()
ax = Axis(fig[1,1])

for i in eachindex(all_ss)
    stable, unstable = all_ss[i]
    scatter!(ax, fill(r_rng[i], length(stable)), stable, color=Makie.wong_colors()[5], alpha=0.6)
    scatter!(ax, fill(r_rng[i], length(unstable)), unstable,
             color=Makie.wong_colors()[2], marker=:utriangle, alpha=0.6)
end

lines!(ax, -1..0, x -> 0, color=:black)
lines!(ax, 0..1, x -> 0, color=:black, linestyle=:dash)

lines!(ax, -0.25..0, r -> sqrt(0.5 - sqrt(0.25 + r)), color=:black, linestyle=:dash)
lines!(ax, -0.25..0, r -> -sqrt(0.5 - sqrt(0.25 + r)), color=:black, linestyle=:dash)

lines!(ax, -0.25..1, r -> sqrt(0.5 + sqrt(0.25 + r)), color=:black)
lines!(ax, -0.25..1, r -> -sqrt(0.5 + sqrt(0.25 + r)), color=:black)

ax.xlabel = L"$r$ - parameter"
ax.ylabel = L"$x$ - steady state position"

elem_1 = [LineElement(color = :black, linestyle = nothing)]
elem_2 = [LineElement(color = :black, linestyle = :dash)]
elem_3 = [MarkerElement(color = Makie.wong_colors()[5], marker = :circle, alpha=0.6)]
elem_4 = [MarkerElement(color = Makie.wong_colors()[2], marker = :utriangle, alpha=0.6)]
axislegend(ax, [elem_1, elem_2, elem_3, elem_4],
           ["theoretical stable", "theoretical unstable", "numerical stable", "numerical unstable"],
          position=:lt)

save("bifurc-diagram.pdf", fig)


fig = Figure()
ax = Axis(fig[1,1])
ax.limits = (nothing, (-10, 10))
ax.xlabel = L"x"
ax.ylabel = L"\dot{x} = f(x, r)"

lines!(ax, -2..2, x -> f(x, -1, 0), label=L"r = -1", color=Makie.wong_colors()[1])
lines!(ax, -2..2, x -> f(x, 1, 0), label=L"r = 1", color=Makie.wong_colors()[2], linestyle=:dash)
vlines!(ax, [-1.5, 1.5], label=nothing, color=:magenta, linestyle=:dashdot)
axislegend(ax)

save("f-func.pdf", fig)

