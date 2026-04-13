using Documenter, NeutralAtoms

push!(LOAD_PATH,"../src/")

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    edit_link = "main",
    assets = [joinpath("assets", "logo.ico")],
    size_threshold_ignore = ["library.md"]
)

makedocs(
    sitename="NeutralAtoms.jl",
    modules=[NeutralAtoms],
    checkdocs=:all,
    format=format,
    authors="M.Y. Goloshchapov",
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API" => "library.md"
    ]
    )

deploydocs(
    repo = "https://github.com/mgoloshchapov/NeutralAtoms.jl",
    devbranch = "main",
)
