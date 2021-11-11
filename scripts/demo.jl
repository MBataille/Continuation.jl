using DrWatson
@quickactivate :Continuation
using Plots

N = 256
Δx = 0.4
p = @dict ϵ=1.16 ν=1
params = GenericParams(N, Δx, p)

S = read_state(params)

xs = getXs(S)
heatmap(xs, xs, S.u)

V = RHS_SHE(S);
heatmap(xs, xs, V)
L = sum(abs.(V)) / params.N^2
m = maximum(abs.(V))