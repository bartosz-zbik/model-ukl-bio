using DifferentialEquations, Statistics, Distributions, Printf
using CairoMakie
set_theme!(theme_latexfonts())

# p[1] - Vector of base frequencies
# p[2] - interactions strength parameter K
function kuramoto!(du, u, p, t)
    r, psi = kuramotoorder(u)
    du .= p[1] .- (p[2] * r) * sin.(u .- psi)
end

function kuramotoorder(u)
    # cis(x) = cos(x) + im * sin(x)
    c = sum(cis.(u)) / length(u)
    r = abs(c) # amplitude
    psi = angle(c) # phase
    return r, psi
end

function anliticalcauchyorder(k)
    k_c = 2.0
    if k < k_c
        return 0.0
    else
        return sqrt((k - 2) / k)
    end
end


# finds a single solution for the Kuramoto model
# with random initial conditions and \omega ~ N(\mu=0, \sigma=1)
function findsol(K, tmax, N)
    u0 = 2pi * rand(N)
    p = [randn(N), K]

    prob = ODEProblem(kuramoto!, u0, (0, tmax), p)
    return solve(prob)
end

# finds a single solution for the Kuramoto model
# with random initial conditions and \omega ~ Cauchy(\mu=0, \sigma=1)
function findcauchysol(K, tmax, N)
    # solves the Kuramoto model
    # with random initial conditions and \omega ~ Cauchy(\mu=0, \sigma=1)
    u0 = 2pi * rand(N)
    p = [rand(Cauchy(0, 1), N), K]

    prob = ODEProblem(kuramoto!, u0, (0, tmax), p)
    return solve(prob)
end

# ####################
# Main, making plots
# ####################

# finds nrep solutions
# plots how r evolves in time
function asymptoticnormal(K, tmax, nrep, N=128)
    s = [findsol(K, tmax, N) for i in 1:nrep]
    trng = range(tmax / 2, tmax, length=30)
    m = mean([kuramotoorder(x(t))[1] for x in s, t in trng])
    fig = Figure()
    ax = Axis(fig[1, 1])
    for x in s
        lines!(ax, 0 .. tmax, t -> kuramotoorder(x(t))[1], color=:black, alpha=0.1, label=nothing)
    end
    lines!(ax, 0 .. tmax, t -> mean([kuramotoorder(x(t))[1] for x in s]), label="<r>")
    ax.xlabel = "time"
    ax.ylabel = L"order parameter $r$"
    ax.title = L"""\omega \sim \mathcal{N}(0, 1),~~ K/K_c = %$(@sprintf "%.2f" K/sqrt(8/pi)),~~ <r> = %$(@sprintf "%.3f" m)"""
    ax.limits = (nothing, nothing, -0.05, 1.05) 
    axislegend(ax, position=:lt)
    return fig
end

# finds nrep solutions
# plots how r evolves in time and compares it to the analytical solution
function asymptoticcauchy(K, tmax, nrep, N=128)
    s = [findcauchysol(K, tmax, N) for i in 1:nrep]
    trng = range(tmax / 2, tmax, length=30)
    m = mean([kuramotoorder(x(t))[1] for x in s, t in trng])
    fig = Figure()
    ax = Axis(fig[1, 1])
    for x in s
        lines!(ax, 0 .. tmax, t -> kuramotoorder(x(t))[1], color=:black, alpha=0.1, label=nothing)
    end
    lines!(ax, 0 .. tmax, t -> mean([kuramotoorder(x(t))[1] for x in s]), label="<r>")
    hlines!(ax, [anliticalcauchyorder(K)], label="theory", linestyle=:dash, color=:magenta)
    ax.xlabel = "time"
    ax.ylabel = L"order parameter $r$"
    ax.title = L"""\omega \sim \mathcal{Cauchy}(0, 1),~~ K/K_c=%$(@sprintf "%.2f" K/2),~~ r=%$(@sprintf "%.3f" anliticalcauchyorder(K)),~~ <r>=%$(@sprintf "%.3f" m)"""
    ax.limits = (nothing, nothing, -0.05, 1.05) 
    axislegend(ax, position=:lt)
    return fig
end


K_c_normal = sqrt(8 / pi)
K_rng_normal = [0, 1 / 4, 1 / 2, 0.97, 1.0, 1.03, 3 / 2, 2.0, 10] * K_c_normal

for i in eachindex(K_rng_normal)
    tmax = 200
    nrep = 30
    N = 512
    f = asymptoticnormal(K_rng_normal[i], tmax, nrep, N)
    save("out/kuramoto_normal_$i.pdf", f)
end


K_c_cauchy = 2.0
K_rng_cauchy = [0, 1 / 4, 1 / 2, 0.97, 1.0, 1.03, 3 / 2, 2.0, 10] * K_c_cauchy

for i in eachindex(K_rng_cauchy)
    tmax = 200
    nrep = 30
    N = 512
    f = asymptoticcauchy(K_rng_cauchy[i], tmax, nrep, N)
    save("out/kuramoto_cauchy_$i.pdf", f)
end

