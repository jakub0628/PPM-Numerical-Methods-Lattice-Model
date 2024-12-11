using Plots
using LinearAlgebra
using Distributions

N = 100
Niters = 100000
β = 1e2

s = Array{Float64}(undef, 2)
λ = Array{typeof(s)}(undef, (N, N))
Λ = Array{typeof(λ)}(undef, Niters)

function get_spin()
    spin = rand(Uniform(-1, 1), 2)
    spin = spin ./ norm(spin)
    return(spin)
end

for i ∈ (1:N)
    for j ∈ (1:N)
        λ[i, j] = get_spin()
    end
end

Λ[1] = λ

function energy(i, j, sC)
    sN = λ[mod1(i+1, N), j           ]
    sS = λ[mod1(i-1, N), j           ]
    sW = λ[i,            mod1(j-1, N)]
    sE = λ[i,            mod1(j+1, N)]

    E = - dot(sC, (sN + sS + sW + sE))
    return(E)
end

function flip(λ, i, j)
    sF = get_spin()
    ΔE = energy(i, j, sF) -  energy(i, j, λ[i, j])

    p = rand(Float64)
    pʙ = exp(-β * ΔE)

    if ΔE < 0 || pʙ > p
        λ[i, j] = sF
    end
    
    return(λ)
end

function get_theta(λ)
    θ = Array{Float64}(undef, N, N)
    
    for i ∈ (1:N)
        for j ∈ (1:N)
            θ[i, j] = atan(λ[i, j][2] / λ[i, j][1] )
        end
    end
    
    return(θ)
end

for iter in 1:Niters-1
    for i in rand(1:N) 
        for j in rand(1:N)
            Λ[iter+1] = flip(Λ[iter], i, j)
        end 
    end
end

print(Λ[10][1,1])

A = get_theta(Λ[100])

heatmap((1:N), (1:N), A, c=:lightrainbow, aspect_ratio=:equal, axis=([], false), right_margin = 10Plots.mm)
