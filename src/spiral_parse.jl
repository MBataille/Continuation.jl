using NPZ
using PyCall

export read_spiral_state, read_spiral_x, SpiralParams

function __init__()
    
end

function SpiralParams(N::Int, σ::Real, α::Real, γ::Real)
    Δx = 2π / (N - 1)
    _p = @dict σ α γ
    GenericParams(N, Δx, _p)
end

function SpiralParams(N::Int, σ::Real, α::Real, γ::Real, vx::Real, vy::Real, Ω::Real)
    Δx = 2π / (N - 1)
    _p = @dict σ α γ vx vy Ω
    GenericParams(N, Δx, _p)
end

function SpiralParams(gp::GenericParams, vx::Real, vy::Real, Ω::Real)
    SpiralParams(gp.N, gp.p[:σ], gp.p[:α], gp.p[:γ], vx, vy, Ω)
end

function spiral_p_to_name(gp::GenericParams, ext::String)
    _p = @strdict N=gp.N sigma=gp.p[:σ] alpha=gp.p[:α] gamma=gp.p[:γ]
    datadir("sims", "spirals", savename(_p, ext, digits=4))
end

function spiral_p_to_name(pre::String, gp::GenericParams, ext::String)
    _p = @strdict N=gp.N sigma=gp.p[:σ] alpha=gp.p[:α] gamma=gp.p[:γ]
    datadir("sims", "spirals", savename(pre, _p, ext, digits=4))
end

function read_scalars(filename)
    lines = readlines(filename)
    splitted = split(lines[2], " ")
    [parse(Float64, splitted[1]), parse(Float64, splitted[3])]
end

function read_spiral_x(gp::GenericParams; sym::Bool=false, switch::Bool=false, vx::Real=0.0)
    filename = spiral_p_to_name("x", gp, "npz")

    py"""
    import numpy as np

    def readX(filename):
        data = np.load(filename)
        return data['X1']
    """

    X = py"readX"(filename)
    N_rows = if sym Int(gp.N/2) else gp.N end
    N_cols = gp.N
    println(N_rows)

    N2 = N_cols * N_rows
    #u = permutedims(reshape(X[1:N2] + im * X[N2+1:2N2], (N_cols, N_rows)), [2, 1])
    u = Array(transpose(reshape(X[1:N2] + im * X[N2+1:2N2], (N_cols, N_rows))))
    if switch
        if sym
            vy = 0
            Ω, α = X[end-1: end]
        else
            Ω, α, vy = X[end-2:end] # last 3 elements
        end
    else
        if sym
            Ω, vx = X[end-1: end]
            α = gp.p[:α]
            vy = 0
        else
            Ω, vx, vy = X[end-2:end]
            α = gp.p[:α]
        end
    end
    gp = SpiralParams(gp.N, gp.p[:σ], α, gp.p[:γ], vx, vy, Ω)
    flags = Flags(sym=sym, switch=switch)
    State(gp=gp, u=u, flags=flags)
end

function read_spiral_state(gp::GenericParams)
    filename_npz = spiral_p_to_name(gp, "npz")
    filename_dat = spiral_p_to_name(gp, "dat")
    u = py"readState"(filename_npz)
    vx, vy = read_scalars(filename_dat) * 2π / (10 * gp.N)
    gp = SpiralParams(gp, vx, vy, 0) # change this
    State(gp=gp, u=u)
end