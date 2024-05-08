using DifferentialEquations, Statistics, Distributions
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

function anliticalorder(k)
    k_c = 2.0
    if k < k_c
        return 0.0
    else
        return sqrt((k - 2) / k)
    end
end


# N = 128
# K = 2.0
# tmax = 50.0
# u0 = 2pi * rand(N)
# p = [randn(N), K]
# 
# prob = ODEProblem(kuramoto!, u0, (0, tmax), p)
# sol = solve(prob)


# finds a single solution
function findsol(K, tmax, N=128)
    # solves the Kuramoto model
    # with random initial conditions and \omega ~ Cauchy(\mu=0, \sigma=1)
    u0 = 2pi * rand(N)
    p = [rand(Cauchy(0, 1), N), K]

    prob = ODEProblem(kuramoto!, u0, (0, tmax), p)
    return solve(prob)
end

# finds nrep solutions
# plots how r evolves in time and compares it to the analytical solution
function asymptoticorrder(K, tmax, nrep, N=128)
    s = [findsol(K, tmax, N) for i in 1:nrep]
    trng = range(tmax/2, tmax, length=30)
    m = mean([kuramotoorder(x(t))[1] for x in s, t in trng])
    fig = Figure()
    ax = Axis(fig[1, 1])
    for x in s
        lines!(ax, 0..tmax, t -> kuramotoorder(x(t))[1], color=:black, alpha=0.1, label=nothing)
    end
    lines!(ax, 0..tmax, t -> mean([kuramotoorder(x(t))[1] for x in s]))
    hlines!(ax, [anliticalorder(K)], label="theory")
    ax.xlabel = "time"
    ax.ylabel = L"order parameter $r$" 
    ax.title = "K=$K, r = $(anliticalorder(K)), <r> = $m"
    return fig
end

K_rng = [1/4, 1/2, 0.97, 1.0, 1.03, 3/2, 2.0] * 2.0

for i in eachindex(K_rng)
    tmax = 200
    nrep = 20
    N = 256
    f = asymptoticorrder(K_rng[i],tmax, nrep, N)
    save("kuramoto_$i.pdf", f)
end



function visualizepolar(sol, n, t)
    fig = Figure()
    ax = PolarAxis(fig[1, 1], rlimits=(0, 1.1))
    scatter!(ax, sol(t), ones(n))
    return fig
end

function visualizeorder(sol, tmax)
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, 0 .. tmax, t -> (kuramotoorder(sol(t))[1]))
    return fig
end

