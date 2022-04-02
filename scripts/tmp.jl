# used as a notebook
using Plots, GLM, DataFrames
using Polynomials

xs = collect(5:5:50)
vs = [3.98, 6.3, 7.73, 8.66, 9.23, 9.46, 9.65, 9.86, 9.88, 9.89]

V₀ = 10
yp = log.(V₀ .- vs)

plot(xs, yp)

p = Polynomials.fit(xs, yp, 1) 

tau = -1 / p[1]

b(n) = V₀ / (n * π) * (n % 2 == 0 ? 0 : 2)

T = 100
v_ft(t, N) = V₀ / 2 + sum([b(n) * sin(2π * n * t / T) for n ∈ 1:N])

v15 = [v_ft(t, 15) for t ∈ 5:5:40]
plot(collect(1:100), v15)

v100 = [v_ft(t, 100) for t ∈ 5:5:40]
plot!(collect(1:100), v100)

vs = vs[1:8]
yp = log.(v100 .- vs)
xs = collect(5:5:40)
