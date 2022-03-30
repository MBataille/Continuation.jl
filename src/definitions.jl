using FFTW, Parameters, DiffEqOperators, Setfield, Parameters
FFTW.set_provider!("mkl")
using LinearAlgebra, Plots, SparseArrays


#============== Exports =============#
export GenericParams, Params, Flags, State, getXs
export ∂ₓ, ∂ₓₓ, ∂ₓₓₓₓ, ∇², ∇⁴, LaplaceOperator

#================ Classes ===============#

#  PARAMS contain all relevant parameters
# of the state. They will later be stored
# in State.
struct GenericParams
    N::Int
    Δx::Real
    p::Dict
end

#  FLAGS contain all relevant flags of the
# state, i.e. wether to compute derivatives
# using FFT or FD, etc.
@with_kw struct Flags
    useFFT::Bool = true
    sym::Bool = false
    switch::Bool = false
end

#  STATE is the main class, it will contain
# the solution u, and relevant parameters
# in the form of a dictionary, it MUST contain
# N and Δx
@with_kw struct State
    gp::GenericParams
    u::Matrix
    flags::Flags = Flags()
end

#============== Base functions =============#

function LaplaceOperator(N, Δx)
    D2x = CenteredDifference(2, 2, Δx, N)
    Qx = PeriodicBC(Float64)
    sp = sparse(D2x * Qx)[1]
    kron(sparse(I,  N, N), sp) + kron(sp, sparse(I, N, N))
end

function LaplaceOperator1D(N, Δx)
    D2x = CenteredDifference(2, 2, Δx, N)
    Qx = Dirichlet0BC(typeof(Δx))
    sparse(D2x * Qx)[1]
end

#  Base function that computes matrix derivatives
# along x or y, for arbitrary order.
function deriv_mat(A::Matrix, Δx::Real; axis::Char='x', order::Int=1, useFFT::Bool=true)
    
    N = size(A)[1];    
    if useFFT
        freqs = fftfreq(N) * N;
        L = (N-1) * Δx;

        if axis == 'y'
            dim = 1;
            k = [freqs[i] for i in 1:N, j in 1:N];
        elseif axis == 'x'
            dim = 2;
            k = [freqs[j] for i in 1:N, j in 1:N]; 
        else
            println("Axis not understood");
        end
        real.(ifft( (2π / L * im) ^ order * k .^ order .* fft(A) ))
    else
        B = similar(A) # assuming order = 2
        for i in 1:N, j in 1:N
            j₋ = if j > 1 j - 1 else N end
            j₊ = if j < N j + 1 else 1 end
            i₊ = if i < N i + 1 else 1 end
            i₋ = if i > 1 i - 1 else N end
            B[j, i] = A[j₋, i] + A[j₊, i] + A[j, i₋] + A[j, i₊] - 4 * A[j, i] 
        end
        B
    end
end

function getXs(gp::GenericParams)
    (0:gp.N-1) * gp.Δx
end

function getXs(S::State)
    getXs(S.gp)
end

### Definitions to simplify the code

∂ₓ(S::State; axis='x') = deriv_mat(S.u, S.gp.Δx, axis=axis, order=1, useFFT=S.flags.useFFT);
∂ₓₓ(S::State; axis='x') = deriv_mat(S.u, S.gp.Δx, axis=axis, order=2, useFFT=S.flags.useFFT);
∂ₓₓₓₓ(S::State; axis='x') = deriv_mat(S.u, S.gp.Δx, axis=axis, order=4, useFFT=S.flags.useFFT);

∇²(S::State) = ∂ₓₓ(S) + ∂ₓₓ(S, axis='y');
∇⁴(S::State) = ∂ₓₓₓₓ(S) + ∂ₓₓₓₓ(S, axis='y');

function ∇²(N::Int, Δx::Real; dim::Int=2)
    if dim == 2
        LaplaceOperator(N, Δx)
    elseif dim == 1
        LaplaceOperator1D(N, Δx)
    end
end

∇⁴(N::Int, Δx::Real; dim::Int=2) = ∇²(N, Δx; dim)^2

∇²(V::Array, N::Int, Δx::Real; dim::Int=2) = ∇²(N, Δx; dim) * V
∇⁴(V::Array, N::Int, Δx::Real; dim::Int=2) = ∇⁴(N, Δx; dim) * V


