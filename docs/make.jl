using GraphsOptim
using Documenter

DocMeta.setdocmeta!(GraphsOptim, :DocTestSetup, :(using GraphsOptim); recursive=true)

makedocs(;
    modules=[GraphsOptim],
    authors="Guillaume Dalle <22795598+gdalle@users.noreply.github.com> and contributors",
    repo="https://github.com/gdalle/GraphsOptim.jl/blob/{commit}{path}#{line}",
    sitename="GraphsOptim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://gdalle.github.io/GraphsOptim.jl",
        assets=String[],
        edit_branch="main",
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/gdalle/GraphsOptim.jl",
    devbranch="main",
)
