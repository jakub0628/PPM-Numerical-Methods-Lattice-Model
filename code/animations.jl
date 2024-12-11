include("nvector.jl")

g1ta = Grid(d=1, β=1e3)
gf1ta = evolve(g1ta, store=true, niters=10_000_000)
vis(g1ta, anim=true,hist=gf1ta)