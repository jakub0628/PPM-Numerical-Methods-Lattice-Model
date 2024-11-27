using Plots
using LinearAlgebra
using Distributions

N = 100
β = 1e2

function get_spin()
    spin = rand(Uniform(-1, 1), 2)
    spin = spin ./ norm(spin)
end

Λ = Array{Vector}(undef, N, N)

for i ∈ (1:N)
    for j ∈ (1:N)
        Λ[i, j] = get_spin()
    end
end

function energy(i, j, sC)
    sN = Λ[mod1(i+1, N), j           ]
    sS = Λ[mod1(i-1, N), j           ]
    sW = Λ[i,            mod1(j-1, N)]
    sE = Λ[i,            mod1(j+1, N)]

    E = - dot(sC, (sN + sS + sW + sE))
    return(E)
end

function flip(i, j)
    sF = get_spin()
    ΔE = energy(i, j, sF) -  energy(i, j, Λ[i, j])

    p = rand(Float64)
    pʙ = exp(-β * ΔE)

    if ΔE < 0 || pʙ > p
        Λ[i, j] = sF
    else
        end
end

for iter in 1:1e7
    for i in rand(1:N) 
        for j in rand(1:N) 
            flip(i, j)
        end
    end
end

A = Array{Float64}(undef, N, N)

for i ∈ (1:N)
    for j ∈ (1:N)
        A[i, j] = atan(Λ[i, j][2] / Λ[i, j][1] )
    end
end

heatmap((1:N), (1:N), A, c=:lightrainbow, aspect_ratio=:equal, axis=([], false), right_margin = 10Plots.mm)