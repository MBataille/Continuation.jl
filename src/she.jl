#============= SHE related functions ============#

export RHS_SHE, Jacobian_SHE


function RHS_SHE(S::State)
    @unpack ϵ, ν = S.gp.p
    u = S.u
    ϵ * u .- u .^ 3 .- (ν * ∇²(S)) .- ∇⁴(S)
end

# Vectorized (flattened) version of the RHS 
function RHS_SHE(V::Array, p::Dict, N::Int, Δx::Real)
    @unpack ϵ, ν = p
    ϵ * V .- V .^3 .- (ν * ∇²(V, N, Δx)) .- ∇⁴(V, N, Δx)
end

# Jacobian will go here
function Jacobian_SHE(V::Array, p::Dict, N::Int, Δx::Real)
    @unpack ϵ, ν = p
    spdiagm(0 => ϵ .- 3 .* V .^ 2 ) .- ν * ∇²(N, Δx) .- ∇⁴(N, Δx)
end
