using FFTW
FFTW.set_provider!("mkl")

#================ Classes ===============#

#  PARAMS contain all relevant parameters
# of the state. They will later be stored
# in State.
struct Params
    N::Int
    L::Float64  
    ϵ::Float64
    ν::Float64
end

#  FLAGS contain all relevant flags of the
# state, i.e. wether to compute derivatives
# using FFT or FD, etc.
struct Flags
    useFFT::Bool
end

#  STATE is the main class, it will contain
# the solution u, and relevant parameters
struct State
    p::Params
    u::Matrix
end

#============== Base functions =============#

#  Base function that computes matrix derivatives
# along x or y, for arbitrary order.
function deriv_mat(A::Matrix{Float64}, L::Float64; axis::Char='x', order::Int64=1)
    
    #factor = (2π * im / L);
    N = size(A)[1];
    freqs = fftfreq(N) * N;
    
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

### Definitions to simplify the code

∂ₓ(S::State; axis='x') = deriv_mat(S.u, S.p.L, axis=axis, order=1);
∂ₓₓ(S::State; axis='x') = deriv_mat(S.u, S.p.L, axis=axis, order=2);
∂ₓₓₓₓ(S::State; axis='x') = deriv_mat(S.u, S.p.L, axis=axis, order=4);

∇²(S::State) = ∂ₓₓ(S) + ∂ₓₓ(S, axis='y');
∇⁴(S::State) = ∂ₓₓₓₓ(S) + ∂ₓₓₓₓ(S, axis='y');


#============== Exports =============#
export Params, Flags, State
export ∂ₓ, ∂ₓₓ, ∂ₓₓₓₓ, ∇², ∇⁴



