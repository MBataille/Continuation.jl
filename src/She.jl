#============= SHE related functions ============#

export RHS_SHE


function RHS_SHE(S::State)
    ϵ = S.gp.p[:ϵ]
    ν = S.gp.p[:ν]
    u = S.u
    ϵ * u .- u .^ 3 .- (ν * ∇²(S)) .- ∇⁴(S)
end

# Vectorized (flattened) version of the RHS 
function RHS_SHE(V::Array, gp::GenericParams)
    u = reshape(V, (gp.N, gp.N))
    vec(RHS_SHE(State(u=u, gp=gp)))
end
# Jacobian will go here

function jacobian_element()
 # sth
end
