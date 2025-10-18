using BioinfoTools
using Documenter

DocMeta.setdocmeta!(BioinfoTools, :DocTestSetup, :(using BioinfoTools); recursive=true)

makedocs(;
    modules=[BioinfoTools],
    authors="thomaswolfe <thomas_wolfe@student.uml.edu> and contributors",
    sitename="BioinfoTools.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
