using Core
using DelimitedFiles, Formatting

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

# converts float to string replacing the 
# decimal point with a given separator
function float_to_str(fl::Float64; sep="_")
    s = format(fl)
    replace(s, "." => sep)
end

function params_to_filename(p::Params)
    folder = format("data/N{}_dx{}/", p.N, float_to_str(p.L / (p.N-1)))
    file = format("epsilon{}.values", float_to_str(p.ϵ, sep="."))
    folder * file
end

# loads state for given params
function read_state(p::Params)
    filename = params_to_filename(p)
    clean_txt_file(filename)
    u = readdlm(filename, ' ', Float64, '\n')

    S = State(p, u) 
end

export read_state