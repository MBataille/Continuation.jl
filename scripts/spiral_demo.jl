using DrWatson
@quickactivate :Continuation # activate package +  using

N = 128
σ = 0.4
α = 0.8832
γ = 0.01

params = SpiralParams(N, σ, α, γ)
#S = read_spiral_state(params)
S = read_spiral_x(params, sym=true)

xs = range(-π, π, length=S.gp.N)


g_ft = get_g_ft(S.gp)

using Plots
heatmap(1:size(S.u, 1), 1:size(S.u, 2), abs.(S.u) )#, c=:hsv)
heatmap(xs, xs, abs.(g_ft))