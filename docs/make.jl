using Documenter, MetidaBioeq
#using DocumenterLaTeX

makedocs(
    modules = [MetidaBioeq],
    sitename = "MetidaBioeq.jl",
    authors = "Vladimir Arnautov",
    pages = [
        "Home" => "index.md",
        "Home" => "api.md",
    ],
)

deploydocs(repo = "github.com/PharmCat/MetidaBioeq.jl.git", push_preview = true,
)
