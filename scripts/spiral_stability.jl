using DrWatson
@quickactivate :Continuation

using SparseArrays, Pardiso, LinearAlgebra
using PyCall, Printf

ddir = "/media/hdd/chimera/pathfollow_fftg001/"
#ddir = "../../chimera/pathfollow_fftg001/"


function readJac(filename)
    py"""
    import numpy as np
    import scipy.sparse as sp
    def read_jac(filename):
        mat = sp.coo_matrix(sp.load_npz(filename))
        return mat.row + 1, mat.col + 1, mat.data
    """
    sparse(py"read_jac"(filename)...)
end

function readJac(N, σ, α, γ)
    wdir = params_to_fname(N, σ, α, γ)
    if isfile(wdir * "/jac2.npz")
        readJac(wdir * "/jac2.npz")
    else
        readJac(wdir * "/jac1.npz")
    end
end

function float2str(fl::Real)
    replace(string(fl), "." => "")
end

function a2str(α)
    replace((@sprintf "a%.4f" α), "." => "")
end

function params_to_fname(N, σ, α, γ)
    ddir * "N$N" * "s040" * a2str(α) * "g" * float2str(γ)
end

N, σ, α, γ = 128, 0.4, 0.8073, 0.01
jac = readJac(N, σ, α, γ)

using PyPlot
pygui(true)
spy(jac)

diff = jac - transpose(jac)

maximum(diff)

ps = PardisoSolver()

val, vec = eigen(jac)

using BenchmarkTools

# 2.86 GB for N=64 ~ 7s
@btime lu(jac);
F = lu(jac);

detlog = sum(log10.(Complex.(diag(F.L)))) + sum(log10.(Complex.(diag(F.U))))

det(F)
