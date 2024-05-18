using DifferentialEquations
using CairoMakie
set_theme!(theme_latexfonts())


# p = (k, r)
function model(n, p, t)
    return p[2] * n * (1.0 - n / p[1]) - n^2 / (1 + n^2)
end


function findsol(k, r, x0)
    tspan = (0.0, 50.0)
    prob = ODEProblem(model, x0, tspan, (k, r))
    sol = solve(prob)
    return sol
end

function make_plot(k, r, x0)
    sol = findsol(k, r, x0)
    tspan = 0..50

    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, tspan, t -> sol(t))
    ax.title = L"k = %$k,~~ r = %$r,~~ x_0 = %$(x0)"
    ax.xlabel = "time"
    ax.ylabel = "population density"
    ax.limits = (nothing, (0.0, nothing))

   return fig 
end

# k = 8
k8_x0_rng = [0.5, 1.0, 3.0, 5.0]
for i in eachindex(k8_x0_rng)
    save("out/model_k8_$i.pdf", make_plot(8.0, 0.5, k8_x0_rng[i]))
end

# k = 5
k5_x0_rng = [0.5, 1.0]
for i in eachindex(k5_x0_rng)
    save("out/model_k5_$i.pdf", make_plot(5.0, 0.5, k5_x0_rng[i]))
end

