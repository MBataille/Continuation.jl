using Continuation
using Documenter

DocMeta.setdocmeta!(Continuation, :DocTestSetup, :(using Continuation); recursive=true)

makedocs(;
    modules=[Continuation],
    authors="Martin Bataille",
    repo="https://github.com/mbataille/Continuation.jl/blob/{commit}{path}#{line}",
    sitename="Continuation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mbataille.github.io/Continuation.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mbataille/Continuation.jl",
    devbranch="main",
)
