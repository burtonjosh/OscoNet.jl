using OscoNet
using Documenter

DocMeta.setdocmeta!(OscoNet, :DocTestSetup, :(using OscoNet); recursive = true)

makedocs(;
    modules = [OscoNet],
    authors = "Joshua Burton",
    repo = "https://github.com/burtonjosh/OscoNet.jl/blob/{commit}{path}#{line}",
    sitename = "OscoNet.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://burtonjosh.github.io/OscoNet.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => ["Tutorial" => "examples/tutorial.md", "Application to real data" => "examples/real-data.md"],
        "Library" => ["Public" => "lib/public.md", "Internals" => "lib/internals.md"]
        ],
)

deploydocs(; repo = "github.com/burtonjosh/OscoNet.jl", devbranch = "main")
