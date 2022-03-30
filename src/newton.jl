using Pardiso

export newton, L₁

L₁(V::Array) = sum(abs.(V)) / length(V)
L₂(V::Array, xs::Array) = trapz(xs, V.^2) 

function newton(V₀::Array, rhs::Function, jac::Function, p::Dict, N::Int, Δx::Real)
    max_iter = 100
    tol = 1e-8

    V = copy(V₀)

    ps = MKLPardisoSolver()

    for i in 1:max_iter
        Y = rhs(V, p, N, Δx)
        error = L₁(Y)
        
        println("Running iteration n $i, error is $error")
        if error < tol
            break
        end

        J = jac(V, p, N, Δx)
        Nx, Ny = size(J)
        Nv = length(V)
        V += solve(ps, J, -Y)
    end
    V
end