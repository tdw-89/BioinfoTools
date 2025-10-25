using BioinfoTools
using Documenter

DocMeta.setdocmeta!(BioinfoTools, :DocTestSetup, :(using BioinfoTools); recursive=true)

const pages = [
    "Overview" => "index.md",
    "API reference" => "api.md",
]

html_config = Documenter.HTML(
    prettyurls = get(ENV, "CI", "false") == "true",
    edit_link = "main",
    highlightjs = "github",
    assets = ["assets/extra.css"],
)

makedocs(; 
    modules = [BioinfoTools],
    authors = "thomaswolfe <thomas_wolfe@student.uml.edu> and contributors",
    sitename = "BioinfoTools.jl",
    format = html_config,
    pages = pages,
    warnonly = Documenter.except(:missing_docs),
)
