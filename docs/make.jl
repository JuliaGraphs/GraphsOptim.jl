using GraphsOptim
using Documenter

DocMeta.setdocmeta!(GraphsOptim, :DocTestSetup, :(using GraphsOptim); recursive=true)

generated_path = joinpath(@__DIR__, "src")
base_url = "https://github.com/JuliaGraphs/GraphsOptim.jl/blob/main/"
isdir(generated_path) || mkdir(generated_path)

open(joinpath(generated_path, "index.md"), "w") do io
    # Point to source license file
    println(
        io,
        """
        ```@meta
        EditURL = "$(base_url)README.md"
        ```
        """,
    )
    # Write the contents out below the meta block
    for line in eachline(joinpath(dirname(@__DIR__), "README.md"))
        println(io, line)
    end
end

makedocs(;
    modules=[GraphsOptim],
    authors="Guillaume Dalle, Aurora Rossi and contributors",
    sitename="GraphsOptim.jl",
    format=Documenter.HTML(;
        repolink="https://github.com/JuliaGraphs/GraphsOptim.jl/",
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGraphs.github.io/GraphsOptim.jl",
        assets=String[],
        edit_link=:commit,
    ),
    pages=["Home" => "index.md", "Algorithms" => "algorithms.md"],
)

deploydocs(; repo="github.com/JuliaGraphs/GraphsOptim.jl", devbranch="main")
