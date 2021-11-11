using Trapz

export RHS_Spiral, get_g_ft

conj(S::State) = State(gp=S.gp, u=conj.(S.u))

function coupling_kernel(x::Float64, y::Float64, σ::Float64)
    x = if (x > π)  x - 2π else x end
    y = if (y > π)  y - 2π else y end
    if (x^2 + y^2 ≤ (π * σ)^2)  1 / (π^3 * σ^2) else 0 end
end

function get_g_ft(gp::GenericParams)
    xs = getXs(gp)
    fft([coupling_kernel(x, y, gp.p[:σ]) for y in xs, x in xs])
end

get_g_ft(S::State) = get_g_ft(S.gp)

function extend_symmetric()
    # a
end

function RHS_Spiral(S::State, S₀_conj::State, g_ft::Matrix)
    Gz = ifft(g_ft .* fft(S.u)) * (exp(-im * α) * (2π / S.N) ^ 2)
    z2_G_conjz = S.u ^ 2 .* ifft(g_ft .* fft(conj.(S.u))) * (exp(im * α) * (2π / S.N) ^ 2)

    @unpack α, γ, Ω, vx, vy = S.gp.p

    V = vx * ∂ₓ(S) .+ vy * ∂ₓ(S, axis='y') .- (γ + im * Ω) * S.u .+ 0.5 * Gz .- 0.5 * z2_G_conjz
    
    xs = range(-π, π, length=S.gp.N)
    g = real(trapz((xs, xs), S.u .* ∂ₓ(S₀_conj)))
    h = real(trapz((xs, xs), S.u .* ∂ₓ(S₀_conj, axis='y')))
    i = real(trapz((xs, xs), im * S.u .* S₀_conj.u))

    V, [g, h, i]
end

# function RHS_Spiral_Sym()