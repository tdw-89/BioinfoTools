using BioinfoTools
using Test

@testset "BioinfoTools.jl" begin
    include("GenomeTypes_tests.jl")
    include("GenomicData_tests.jl")
    include("EnrichmentUtils_tests.jl")
    include("LoadGFF_tests.jl")
    include("ParalogUtils_tests.jl")
end