#============= SHE related functions ============#

export RHS_SHE, Jacobian_SHE, RHS_SHE1D, Jacobian_SHE1D


function RHS_SHE(S::State)
    @unpack ϵ, ν = S.gp.p
    u = S.u
    ϵ * u .- u .^ 3 .- (ν * ∇²(S)) .- ∇⁴(S)
end

# Vectorized (flattened) version of the RHS 
function RHS_SHE(V::Array, p::Dict, N::Int, Δx::Real; dim::Int=2)
    @unpack ϵ, ν = p
    ϵ * V .- V .^3 .- (ν * ∇²(V, N, Δx; dim=dim)) .- ∇⁴(V, N, Δx; dim=dim)
end

# Jacobian will go here
function Jacobian_SHE(V::Array, p::Dict, N::Int, Δx::Real; dim::Int=2)
    @unpack ϵ, ν = p
    spdiagm(0 => ϵ .- 3 .* V .^ 2 ) .- ν * ∇²(N, Δx; dim=dim) .- ∇⁴(N, Δx; dim=dim)
end

## One dimensional case ##

RHS_SHE1D(V::Array, p::Dict, N::Int, Δx::Real) = RHS_SHE(V, p, N, Δx; dim=1)
Jacobian_SHE1D(V::Array, p::Dict, N::Int, Δx::Real) = Jacobian_SHE(V, p, N, Δx; dim=1)



