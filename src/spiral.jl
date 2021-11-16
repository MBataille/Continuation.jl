using Trapz

export RHS_Spiral, get_g_ft, extend_symmetric, conj

import Base:conj
conj(S::State) = State(gp=S.gp, u=conj.(S.u))

function coupling_kernel(x::Real, y::Real, σ::Real)
    x = if (x > π)  x - 2π else x end
    y = if (y > π)  y - 2π else y end
    if (x^2 + y^2 ≤ (π * σ)^2)  1 / (π^3 * σ^2) else 0 end
end

function get_g_ft(gp::GenericParams)
    xs = getXs(gp)
    fft([coupling_kernel(x, y, gp.p[:σ]) for y in xs, x in xs])
end

get_g_ft(S::State) = get_g_ft(S.gp)

function extend_symmetric(A::Matrix; dim::Int=1)
    cat(A, reverse(A, dims=dim), dims=dim)
end

function RHS_Spiral(S::State, S₀_conj::State, g_ft::Matrix)
    @unpack α, γ, Ω, vx, vy = S.gp.p

    Gz = ifft(g_ft .* fft(S.u)) * (exp(-im * α) * (2π / S.gp.N) ^ 2)
    z2_G_conjz = S.u .^ 2 .* ifft(g_ft .* fft(conj.(S.u))) * (exp(im * α) * (2π / S.gp.N) ^ 2)

    V = vx * ∂ₓ(S) .+ vy * ∂ₓ(S, axis='y') .- (γ + im * Ω) * S.u .+ 0.5 * Gz .- 0.5 * z2_G_conjz
    
    xs = range(-π, π, length=S.gp.N)
    g = real(trapz((xs, xs), S.u .* ∂ₓ(S₀_conj)))
    h = real(trapz((xs, xs), S.u .* ∂ₓ(S₀_conj, axis='y')))
    i = real(trapz((xs, xs), im * S.u .* S₀_conj.u))

    V, [g, h, i]
end

function RHS_Spiral(V::Array, V₀_conj::Array, g_ft::Array, gp::GenericParams)
    S = State(gp=gp, u=reshape(V, (gp.N, gp.N)))
    S₀_conj = State(gp=gp, u=reshape(V₀_conj, (gp.N, gp.N)))
    V, pinning = RHS_Spiral(S, S₀_conj, g_ft)
    append!(vec(V), pinning)
end

function RHS_Spiral_Sym()
    # a
end