using BioinfoTools
using Test

@testset "BioinfoTools.jl" begin
    include("Types_tests.jl")
    include("Data_tests.jl")
    include("Enrichment_tests.jl")
    include("GFF_tests.jl")
    include("Paralogs_tests.jl")
end