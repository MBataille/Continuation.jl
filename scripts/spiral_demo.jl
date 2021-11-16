using DrWatson
@quickactivate :Continuation # activate package +  using

using Plots
hmap(A::Matrix) = heatmap(1:size(A, 1), 1:size(A, 2), abs.(A) )

N = 128
σ = 0.4
α = 0.8832
γ = 0.01

params = SpiralParams(N, σ, α, γ)
#S = read_spiral_state(params)
S = read_spiral_x(params, sym=true)
U2 = extend_symmetric(S.u, dim=1)

S₁ = State(gp=S.gp, u=U2)
S₀_conj = conj(S₁)
g_ft = get_g_ft(S.gp)

V, (g, h, i) = RHS_Spiral(S₁, S₀_conj, g_ft)

hmap(V)

extend_symmetric([1 2; 3 4; 5 6; 7 8], dim=2)

U2 = vcat(S.u, reverse(S.u, dims=1))

X = vec(S.u)
#X2 = append!(X, reverse(X))

xs = range(-π, π, length=S.gp.N)

hmap(U2)


import PyPlot

PyPlot.imshow(abs.(U2))
PyPlot.savefig(plotsdir("spiral", "mods.png"))
