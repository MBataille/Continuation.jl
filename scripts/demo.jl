using DrWatson
@quickactivate("Continuation")
using Continuation
using Plots

params = @dict ϵ=1.165

S = read_state(params)

xs = getXs(S)
heatmap(xs, xs, S.u)

V = RHS_SHE(S)
heatmap(xs, xs, V)
