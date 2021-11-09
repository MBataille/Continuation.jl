#using FFTW
function RHS_Spiral(S::State, g_ft::Matrix)
    A_ft = fft(S.u)
    Gz = ifft(g_ft .* fft(S.u)) * (exp(-im * α) * (2π / S.N) ^ 2)
    z2_G_conjz = S.u ^ 2 .* ifft(g_ft .* fft(conj.(S.u))) * (exp(im * α) * (2π / S.N) ^ 2)

    γ = S.gp.p[:γ]
    Ω = S.gp.p[:Ω]
    α = S.gp.p[:α]
    vx = S.gp.p[:vx]
    vy = S.gp.p[:vy]

    vx * ∂ₓ(S) .+ vy * ∂ₓ(S, axis='y') .- (γ + im * Ω) * S.u .+ 0.5 * Gz .- 0.5 * z2_G_conjz
end