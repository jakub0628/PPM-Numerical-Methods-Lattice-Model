using Plots
using Ising2D

N = 256
Λ = rand([1, -1], (N, N))
β = 1e10

function energy(i, j)
    sN = Λ[mod1(i+1, N), j           ]
    sS = Λ[mod1(i-1, N), j           ]
    sW = Λ[i,            mod1(j-1, N)]
    sE = Λ[i,            mod1(j+1, N)]

    E = - Λ[i, j] * (sN + sS + sW + sE)
    return(E)
end

function flip(i, j)
    ΔE = -2 * energy(i, j)
    p = rand(Float64)
    pʙ = exp(-β * ΔE)

    if ΔE < 0 || pʙ > p
        Λ[i, j] *= -1
    else
        end
end

for iter in 1:1e8
    for i in rand(1:N) 
        for j in rand(1:N) 
            flip(i, j)
        end
    end
end

heatmap(Λ, aspect_ratio=:equal, axis=([], false), cbar=false)