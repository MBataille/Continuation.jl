#=
    CODE FOR SIMULATION AND CONTINUATION OF 
    SWIFT-HOHENBERG with Non-local RAMAN response:

        ∂u/∂t = η + μu - u³ + β ∂²u/∂τ² - 4/3 ∂⁴u/∂τ⁴ + ∫χ(τ-τ')u(τ')dτ'

    Author: Martin Bataille-Gonzalez.


    -------------------------------
    Where,
        χ(τ) = √3/2 a fᵣ * exp(-τ₀τ/τ₂) * sin(τ₀τ/τ₁)
    and,
        a = τ₀(τ₁² + τ₂²) / (τ₁τ₂²)

    Parameters will be,
        τ₀ = 14 fs
        fᵣ = 0.18
        β = -1.5
        y = 0.12 = η (?)    
        σ = -0.35 ⟹ μ = σ - √3/2fᵣ 
    From, Parra-Rivas et al (2021) PRA.
        fᵣ = 0.18
        τ₁ = 12.2 fs
        τ₂ = 32 fs
=#
using DifferentialEquations

using Plots

function χ(τ::Number, p)
    return sqrt(3) / 2 * p[:a] * p[:fᵣ] * exp(-p[:τ₀]/p[:τ₂] * τ) * sin.(p[:τ₀]/p[:τ₁] * τ)
end

function χ(τ::AbstractArray, p)
    [χ(τᵢ, p) for τᵢ=τ]
end

function spectral_derivative(u, Δx, order)
    N = length(u)

    k = fftfreq(N) * N
    L = (N-1) * Δx

    return ifft( (2π / L * im) ^ order * k.^ order .* fft(u))
end

∫χu(u) = ifft(ft_coupling * fft(u))

∂ₓₓ(u, Δx) = spectral_derivative(u, Δx, 2)
∂ₓₓₓₓ(u, Δx) = spectral_derivative(u, Δx, 4)

# Parameters
L = 300
Δx = 0.1
τ = 0:Δx:L
τ₀ = 14
fᵣ = 0.18
β = -1
η = 0.0
σ = -0.35
μ = σ - sqrt(3)/2 * fᵣ
τ₁ = 12.2
τ₂ = 32
a = τ₀ * (τ₁ ^ 2 + τ₂ ^ 2) / (τ₁ * τ₂ ^ 2)

p = (τ₀ = τ₀, fᵣ = fᵣ, β = β, η = η, μ = μ, τ₁ = τ₁, τ₂ = τ₂, a=a, L=L, Δx = Δx)

coupling = χ(τ, p)

plot(τ, coupling)

using FFTW

FFTW.set_provider!("mkl")
ft_coupling = fft(coupling)



coupling[200]
decay = τ₂ / τ₀

function SHRaman(u, p, t)
    return p[:η] .+ p[:μ] * u .- u .^ 3 .+ p[:β] * ∂ₓₓ(u, p[:Δx]) 
            .- 4 / 3 * ∂ₓₓₓₓ(u, p[:Δx]) .+ ∫χu(u)
end

u = im * zeros(length(τ))

tspan = (0.0, 1.0)
prob = ODEProblem(SHRaman, u, tspan, p)
sol = solve(prob, RK4())

plot(τ, real.(sol(1.0)))

u = SHRaman(u, p, 0)