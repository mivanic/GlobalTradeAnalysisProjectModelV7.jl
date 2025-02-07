using Documenter, GlobalTradeAnalysisProjectModelV7

DocMeta.setdocmeta!(GlobalTradeAnalysisProjectModelV7, :DocTestSetup, :(using GlobalTradeAnalysisProjectModelV7); recursive=true)


const _PAGES = [
    "Introduction" => ["index.md"],
]


makedocs(
    sitename="GlobalTradeAnalysisProjectModelV7.jl",
    format = Documenter.HTML(),
    modules = [GlobalTradeAnalysisProjectModelV7],
    pages = _PAGES
)



deploydocs(
    repo = "github.com/mivanic/GlobalTradeAnalysisProjectModelV7.jl",
    branch = "gh-pages",
    push_preview = true
    #versions = ["stable" => "v^", "v#.#" ],
)
