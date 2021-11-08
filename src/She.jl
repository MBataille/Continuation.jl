using Core
#============= SHE related functions ============#

function RHS_SHE(S::State)
    ϵ = S.p.ϵ
    ν = S.p.ν
    u = S.u
    ϵ * u .- u .^ 3 .- (ν * ∇²(S)) .- ∇⁴(S)
end

# Jacobian will go here


export RHS_SHE
