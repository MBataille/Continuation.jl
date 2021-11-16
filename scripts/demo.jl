using DrWatson
@quickactivate :Continuation
using Plots
using Setfield

N = 256
Δx = 0.4
p = @dict ϵ=1.16 ν=1
params = GenericParams(N, Δx, p)

S = read_state(params)

xs = getXs(S)
heatmap(xs, xs, S.u)

S = @set S.flags.useFFT = false

V = RHS_SHE(S);
heatmap(xs, xs, V)
L = sum(abs.(V)) / params.N^2
m = maximum(abs.(V))

X = Array(vec(transpose(S.u)))
Y = RHS_SHE(X, S.gp.p, N, Δx)

J = Jacobian_SHE(X, S.gp.p, N, Δx)

X2 = newton(X2, RHS_SHE, Jacobian_SHE, S.gp.p, N, Δx)

S = @set S.gp.p[:ϵ] = 1.175
heatmap(xs, xs, transpose(reshape(X2, (N, N))))