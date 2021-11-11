using FFTW
using Parameters
FFTW.set_provider!("mkl")


#============== Exports =============#
export GenericParams, Params, Flags, State, getXs
export ∂ₓ, ∂ₓₓ, ∂ₓₓₓₓ, ∇², ∇⁴

#================ Classes ===============#

#  PARAMS contain all relevant parameters
# of the state. They will later be stored
# in State.
struct GenericParams
    N::Int
    Δx::Float64
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

#  Base function that computes matrix derivatives
# along x or y, for arbitrary order.
function deriv_mat(A::Matrix, Δx::Float64; axis::Char='x', order::Int64=1)
    
    #factor = (2π * im / L);
    N = size(A)[1];
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
end

function getXs(gp::GenericParams)
    (0:gp.N-1) * gp.Δx
end

function getXs(S::State)
    getXs(S.gp)
end

### Definitions to simplify the code

∂ₓ(S::State; axis='x') = deriv_mat(S.u, S.gp.Δx, axis=axis, order=1);
∂ₓₓ(S::State; axis='x') = deriv_mat(S.u, S.gp.Δx, axis=axis, order=2);
∂ₓₓₓₓ(S::State; axis='x') = deriv_mat(S.u, S.gp.Δx, axis=axis, order=4);

∇²(S::State) = ∂ₓₓ(S) + ∂ₓₓ(S, axis='y');
∇⁴(S::State) = ∂ₓₓₓₓ(S) + ∂ₓₓₓₓ(S, axis='y');


