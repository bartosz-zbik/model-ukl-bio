using Statistics
using CairoMakie

function logistic(x, mu)
    return 4 * mu * x * ( 1 - x)
end

function difflogistic(x, mu)
    return 4 * mu * (1 - 2 * x) 
end


function cstable(x, mu)
    z = 2 * x + 1//2
    return (z) - floor(z)
end

function lapunov(f, df, x0::T, p, n::Integer) where T<:Number
    x = x0
    l = Vector{T}(undef, n)
    for i in eachindex(l)
        x = f(x, p)
        l[i] = df(x, p)
    end
    return mean(log.(abs.(l)))
end

function gentr(f, x0::T, p, n::Integer) where T<:Number
    tr = Vector{T}(undef, n)
    tr[1] = x0
    for i in 2:n
        tr[i] = f(tr[i-1], p)
    end
    return tr
end

num_of_mu = 10
mu_rng = [BigInt(i) // num_of_mu for i in 1:num_of_mu]
num_of_x0 = 20
x0_rng = [ BigInt(i) // (num_of_x0 + 1) for i in 1:num_of_x0]

#scatter(mu_rng, lapunov.(logistic, difflogistic, BigInt(1)//3, mu_rng, 100)) 
