using Plots

N = 32
Λ = rand([1, -1], (N, N))
β = 1e-10

function getZ()
    Σ = [-1, 1]
    Z = 0
    for spin ∈ Σ, N ∈ Σ, E ∈ Σ, S ∈ Σ, W ∈ Σ
        ΔE = 2 * spin * (N + E + S + W)
        Z += exp(-β * ΔE)
    end
    return(Z)
end

Z = getZ()

function get_neigbours(x)
    x₋ = x == 1 ? N : x-1
    x₊ = x == N ? 1 : x+1
    return(x₋, x₊)
end

function energy(i, j)
    i₋, i₊ = get_neigbours(i)
    j₋, j₊ = get_neigbours(j)
    E = - Λ[i, j] * (Λ[i₋, j] + Λ[i₊, j] + Λ[i, j₋] + Λ[i, j₊])
    return(E)
end

function flip(i, j)
    ΔE = -2 * energy(i, j)
    p = rand(Float64)
    pʙ = exp(-β * ΔE) / Z

    if ΔE < 0 || pʙ > p
        print("ΔE=$(ΔE), p=$(p), pʙ = $(pʙ)\n")
        Λ[i, j] *= -1
    else
        end
end

@gif for iter in 1:1e4
    for i in rand(1:N) 
        for j in rand(1:N) 
            flip(i, j)
            heatmap(Λ)
        end
    end
end
