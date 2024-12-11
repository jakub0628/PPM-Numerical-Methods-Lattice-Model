include("nvector.jl")

Temps = 1e-3:1e-1:4
betas = 1 ./(Temps)
n = length(betas)
E = Vector{Float64}(undef, n)
M = Vector{Float64}(undef, n)
C = Vector{Float64}(undef, n)
X = Vector{Float64}(undef, n)

@threads for idx in 1:n
    i = betas[idx]
    g = Grid(m=100, n=100, d=3, β=i)
    evolve(g, niters=1_000_000)
    E[idx] = av_energy(g)
    M[idx] = av_mag(g)
    C[idx] = CV(g)
    X[idx] = χ(g)
end

scatter(Temps[10:length(Temps)], C[10:length(C)], xlabel="T/k_B",ylabel="CV/K_B",legend=false)
scatter(Temps, C, xlabel="T/k_B",ylabel="CV/K_B",legend=false)
vline!([2.3, 2.3/2, 2.3/3])