using DrWatson
@quickactivate :Continuation

using Revise
using SparseArrays, LinearAlgebra, DiffEqOperators, Setfield, Parameters
using BifurcationKit
using Plots
const BK = BifurcationKit

N = 1000
Δx = 0.1
p = @dict ϵ=1.16 ν=1
params = GenericParams(N, Δx, p)


# initial condition
x = collect((1:N) * Δx); 
x0 = x[500]
σ = 1
sol0 = @. cos(x) * exp(-((x - x0) / (2 * σ))^2)

plot(x, sol0)
# ok

# newton!
p = @set p[:ϵ] -= 0.1
X2 = newton(X2, RHS_SHE1D, Jacobian_SHE1D, p, N, Δx)
plot(x, X2)


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