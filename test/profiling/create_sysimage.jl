using Pkg
Pkg.activate(joinpath(@__DIR__, "../../"))
using PackageCompiler

using BioinfoTools

PackageCompiler.create_sysimage(
    :BioinfoTools, 
    sysimage_path=joinpath(@__DIR__, "BioinfoTools_sysimage.so"),
    precompile_statements_file=joinpath(@__DIR__, "compile_trace.jl")
)