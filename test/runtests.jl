using BioinfoTools
using Test

@testset "BioinfoTools.jl" begin
    include("GenomeTypes_tests.jl")
    include("GenomicData_tests.jl")
end
