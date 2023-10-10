using GraphsOptim
using Documenter

DocMeta.setdocmeta!(GraphsOptim, :DocTestSetup, :(using GraphsOptim); recursive=true)

base_url = "https://github.com/JuliaGraphs/GraphsOptim.jl/blob/main/"
index_path = joinpath(@__DIR__, "src", "index.md")
readme_path = joinpath(dirname(@__DIR__), "README.md")

open(index_path, "w") do io
    println(
        io,
        """
        ```@meta
        EditURL = "$(base_url)README.md"
        ```
        """,
    )
    for line in eachline(readme_path)
        println(io, line)
    end
end

makedocs(;
    modules=[GraphsOptim],
    authors="Guillaume Dalle, Aurora Rossi and contributors",
    sitename="GraphsOptim.jl",
    format=Documenter.HTML(;
        repolink="https://github.com/JuliaGraphs/GraphsOptim.jl",
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliagraphs.org/GraphsOptim.jl",
    ),
    pages=["Home" => "index.md", "Algorithms" => "algorithms.md"],
)

deploydocs(;
    repo="github.com/JuliaGraphs/GraphsOptim.jl", devbranch="main", push_preview=true
)
