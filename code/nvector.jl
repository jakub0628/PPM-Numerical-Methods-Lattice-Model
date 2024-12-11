using Plots, LinearAlgebra, Distributions, Parameters, BenchmarkTools, PlutoUI, Statistics, Distributed, Base.Threads

k_B = 1.38e-23 # define the Boltzmann constant

function spin(d::Int)
        # generate and normalize a d-dimensional vector
        vec = rand(Uniform(-1, 1), d)
        return vec ./ norm(vec)
    end

mutable struct Grid
    # create a struct to generate and store the m*n grid of d-dimensional spins
    m::Int # grid height
    n::Int # grid length
    d::Int # spin dimension
    β::Float64 # thermodynamic beta
    Glauber::Bool # Use Glauber or Boltzmann distribution
    J::Float64 # coupling strength
    shape::String # shape of the lattice
    s::Matrix{Vector{Float64}} # grid of the spins
    sum_vec::Vector{Float64} # preallocate sum of neighbours

    function Grid(; m=100, n=100, d=2, β=1e2, Glauber=false, J=1, shape="Square") # creates and fills the m*n grid of d-dimensional spins using the default values
        new(m, n, d, β, Glauber, J, shape, [spin(d) for i in 1:m, j in 1:n], zeros(d)) 
    end
end

function neighbours!(S::Grid, i::Int, j::Int) # obtains the sum of the neighbours for the energy calculation
    # Get the sum of neighbours of s[i, j], periodic boundary conditions
    m, n = S.m, S.n # get the dimensions of the lattice
    s = S.s # get the spins
    
    NN = s[mod1(i-1, m), j] # North neighbour 
    SS = s[mod1(i+1, m), j] # South neighbour
    WW = s[i, mod1(j-1, n)] # West neighbour
    EE = s[i, mod1(j+1, n)] # East neighbour
    NW, SW = 0, 0 
    
    if S.shape == "Triangular"  # if triangular lattice add the extra neighbours
        NW = s[mod1(i-1, m), mod1(j-1, n)] # North-west neighbour
        SW = s[mod1(i+1, m), mod1(j-1, n)] # South-west neighbour
    end
    
    S.sum_vec .= NN .+ SS .+ WW .+ EE .+ NW .+ SW # return sum of the neighbours
end

function evolve(S::Grid; niters::Int=1_000_000, store::Bool=false, red::Int=10_000)
    # evolves the grid for niters iterations using the Glauber or Boltzmann prob and stores them with the appropriate point reduction if needed

    if store
        T = Vector{Matrix{Vector{Float64}}}(undef, div(niters, red))
        # preallocate history vector
    else
        T = nothing
    end
    
    @threads for iter in 1:niters
        i, j = rand(1:S.m), rand(1:S.n) # choose spin
        sf = spin(S.d) # choose random orientation
        ΔE = S.J * dot(neighbours!(S, i, j), (S.s[i, j] - sf)) # calculate the energy difference
        
        S.s[i, j] = ifelse(S.Glauber, # check if the spin should be flipped
            ifelse(
                (rand() < 1/(1 + exp(-S.β * ΔE))), 
                sf, S.s[i, j]
            ), 
            ifelse(
                (rand() < exp(-S.β * ΔE)) || (ΔE<0), 
                sf, S.s[i, j]
            )
        )
        
        if store && iter % red == 0
            T[div(iter, red)] = copy(S.s) # store the current state
        end
    
    end
    
    return store ? (T, red) : nothing
end

function vis(S::Grid ; anim::Bool=false, hist::Tuple{Vector{Matrix{Vector{Float64}}}, Int}=([[[[0.0]] [[0.0]]]], 0)) 
    # plots the function for 1 or 2 dimensional spins with the option to animate if the evolution was stored
    T, red = hist
    
    if S.d == 1 # if 1D spin
        
        if anim
            p = [[j[1] for j in i] for i in T]
            
            @gif for i in 1:length(T)
                heatmap(
                    p[i], 
                    # title="$(i*red) of $(length(T)*red) for $(ifelse(S.Glauber, "Glauber", "Boltzmann"))", 
                    c=:binary, 
                    showaxis=false, 
                    cbar=false,
                    size=(500, 500),
                    aspect_ratio=:equal,
                    title="↑ / ↓ evolution for T / kʙ = $(1/(k_B * S.β))"
                    )
                    end
        else # if only plot of last state
            p = [i[1] for i in S.s]
            heatmap(
                p, 
                # title="$(ifelse(S.Glauber, "Glauber", "Boltzmann"))", 
                c=:binary, 
                showaxis=false,
                cbar=false, 
                size=(500, 500),
                aspect_ratio=:equal,
                )
        end
    
    elseif S.d == 2 # if 2D spin
        
        if anim
            Angs = [atan.(getindex.(i, 2) ./ getindex.(i, 1)) for i in T]
            @gif for i in 1:length(T)
                heatmap(
                    Angs[i], 
                    clims=(-pi/2, pi/2), 
                    c=:lightrainbow, 
                    axis=([], false),
                    cbar=false,
                    size=(500, 500),
                    aspect_ratio=:equal,
                    # title="$(i*red) of $(length(T)*red) for $(ifelse(S.Glauber, "Glauber", "Boltzmann"))"
                    title="θ evolution for T / kʙ = $(1/(k_B * S.β))"
                    )
                    end
        else # if only plot of last state
            Ang = atan.(getindex.(S.s, 2) ./ getindex.(S.s, 1))
            heatmap(
                Ang, 
                clims=(-pi/2, pi/2), 
                c=:lightrainbow, 
                axis=([], false),
                cbar=false,
                size=(500, 500),
                aspect_ratio=:equal,
                # title="$(ifelse(S.Glauber, "Glauber", "Boltzmann"))"
                )
        end
    
    end
end

function av_energy(S::Grid; vis::Bool=false, hist::Tuple{Vector{Matrix{Vector{Float64}}}, Int}=([[[[0.0]] [[0.0]]]], 0)) # returns the average energy of the grid with the option to display its evolution over the iterations
    T, red = hist
    if vis
        energies = [1/2 * mean([-S.J * dot(neighbours!(S, i, j), t[i, j]) for i in 1:S.m, j in 1:S.n]) for t in T] # calculate total energy of each grid
        plot(energies, xlabel="N of iterations / $red", ylabel="Av. E", legend=false, title="Average Energy over each iteration")
    else
        return 1/2 * mean([-S.J * dot(neighbours!(S, i, j), S.s[i, j]) for i in 1:S.m, j in 1:S.n])
    end
end

function av_mag(S::Grid ;vis::Bool=false, hist::Tuple{Vector{Matrix{Vector{Float64}}}, Int}=([[[[0.0]] [[0.0]]]], 0)) # returns the average magnetization of the grid
    T, red = hist
    if vis
        mags = [norm(mean(t)) for t in T]
        plot(mags, xlabel="N of iterations/$(red)", ylabel="Av. |M|", legend=false, title="Average Magnetisation over each iteration")
    else
        norm(mean(S.s))
    end
end

function CV(S::Grid ;vis::Bool=false, hist::Tuple{Vector{Matrix{Vector{Float64}}}, Int}=([[[[0.0]] [[0.0]]]], 0)) # returns the specific heat capacity of the grid
    T, red = hist
    if vis
        E(i, j, t) = -S.J * dot(neighbours!(S, i, j), t[i, j])
        varE(t) = var([E(i, j, t) for i in 1:S.m, j in 1:S.n])
        cs = [S.β^2 * varE(t) for t in T]
        
        plot(cs, xlabel="N of iterations/$(red)", ylabel="CV/k_B", legend=false, title="Heat capacity over each iteration")
    
    else
        E(i, j) = -S.J * dot(neighbours!(S, i, j), S.s[i, j])
        varE = var([E(i, j) for i in 1:S.m, j in 1:S.n])
        S.β^2 * varE
    end
end

function χ(S::Grid ;vis::Bool=false, hist::Tuple{Vector{Matrix{Vector{Float64}}}, Int}=([[[[0.0]] [[0.0]]]], 0)) # returns the magnetic susceptibility of the grid
    T, red = hist
    if vis
        xs = [S.β * (1 - norm(mean(t))^2) for t in T] # calculate X of each grid
        plot(xs, xlabel="N of iterations/$(red)", ylabel="χ", legend=false, title="Magnetic susceptibility over each iteration")
    else
        S.β * (1 - norm(mean(S.s))^2) 
    end
end

function T_sweep(; ) # performs a temperature sweep 
    return nothing
end