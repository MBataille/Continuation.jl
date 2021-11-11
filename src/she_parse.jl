using DrWatson
using DelimitedFiles, Formatting

export read_state

#============== Parsing functions =============#

# rips blank space before new lines in Seba's data
function clean_txt_file(filename)
    open(filename, "r") do f
        s = read(f, String)
    end
    s = replace(s, " \r\n" => "\n")
    open(filename, "w") do f
        write(f, s)
    end
end

function she_p_to_name(gp::GenericParams)
    _p = @strdict N=gp.N dx=gp.Δx eps=gp.p[:ϵ]
    datadir("sims", "she", savename(_p, "values", digits=4))
end

# loads state for given params
function read_state(gp::GenericParams)
    filename = she_p_to_name(gp)
    clean_txt_file(filename)
    u = readdlm(filename, ' ', Float64, '\n')

    State(gp=gp, u=u) 
end