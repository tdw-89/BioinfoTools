include("../src/genomic_data.jl")
using .GenomicData
using Test
@testset "GenomicData Types" begin
    chrom = ChromData("chr1", UInt16[1, 2, 3, 4])
    @test chrom.name == "chr1"
    @test chrom.signal isa Vector{UInt16}
    sample = SampleData("sampleA", [ChromData("chr1", UInt16[1, 1, 1, 1])], 10)
    @test sample.n_reads == 10
    @test length(sample.chroms) == 1
    experiment = Experiment("exp", [SampleData("sampleA", [ChromData("chr1", UInt16[1, 1, 1, 1])], 10)])
    @test experiment.name == "exp"
    @test length(experiment.samples) == 1
end
@testset "Peak binning" begin
    peak_file = joinpath(@__DIR__, "data", "micro.narrowPeak")
    chrom_lengths_path = tempname()
    open(chrom_lengths_path, "w") do io
        write(io, "DDB0215018\t50000\n")
    end
    experiment = binpeaks(peak_file, chrom_lengths_path)
    @test length(experiment.samples) == 1
    sample = experiment.samples[1]
    chrom = only(sample.chroms)
    @test sum(chrom.signal) > 0
    @test chrom.signal[5]
    @test chrom.signal[1609]
    @test !chrom.signal[1610]
    rm(chrom_lengths_path; force=true)
end
@testset "Replicate averaging" begin
    chrom1 = ChromData("chr1", UInt16[1, 2, 3, 4])
    chrom2 = ChromData("chr1", UInt16[2, 4, 6, 8])
    sample1 = SampleData("rep1", [chrom1], 10)
    sample2 = SampleData("rep2", [chrom2], 20)
    avg = average_bam_replicates([sample1, sample2])
    @test avg.n_reads ≈ 15
    @test all(isapprox.(avg.chroms[1].signal, [1.5, 3.0, 4.5, 6.0]))
    peak_chrom1 = ChromData("chrA", BitVector([true, false, true, false]))
    peak_chrom2 = ChromData("chrA", BitVector([false, true, true, false]))
    peak_sample1 = SampleData("p1", [peak_chrom1], -1)
    peak_sample2 = SampleData("p2", [peak_chrom2], -1)
    avg_peak = average_peak_replicates([peak_sample1, peak_sample2])
    @test avg_peak.chroms[1].signal ≈ [0.5, 0.5, 1.0, 0.0]
end
@testset "Helpers" begin
    warnings = PeakWarnings()
    chrom_dict = Dict("chr1" => falses(10))
    addpeak!(chrom_dict, "chr1", 2, 4, warnings)
    @test chrom_dict["chr1"][2:4] == trues(3)
    @test contiguous_values([0.0, 0.0, 1.0, 0.0], 0.0) == [1:2, 4:4]
    contiguous = contiguous_values(Float64[0, 2, 2, 0, 3])
    @test contiguous == [(2:3, 2.0), (5:5, 3.0)]
    peak_sample = SampleData("floatSample", [ChromData("chr1", [0.0, 1.0, 1.0, 0.0])], -1)
    bed_df = sample_to_bed_df(peak_sample)
    @test size(bed_df, 1) == 1
    @test bed_df.Score[1] == 1000
end
@testset "Gene helpers" begin
    struct FakeScaffold
        name::String
    end
    mutable struct FakeGene
        scaffold::FakeScaffold
        gene_start::Int
        gene_end::Int
        id::String
        binsignals::Vector{BitVector}
        signals::Vector{Vector{UInt16}}
        regions::Vector{Any}
        samples::Vector{String}
        rnas::Vector{Any}
    end
    mutable struct FakeGenome
        genes::Tuple{Vector{String}, Vector{FakeGene}}
    end
    gene = FakeGene(FakeScaffold("chr1"), 2, 4, "geneA", BitVector[], Vector{Vector{UInt16}}[], Any[], String[], Any[])
    genome = FakeGenome((String["geneA"], [gene]))
    chrom = ChromData("chr1", BitVector([false, true, true, false, false]))
    sample = SampleData("binSample", [chrom], -1)
    experiment = Experiment("peak", [sample])
    addtogenes!(genome, experiment; peak_data=true, regions=false)
    @test !isempty(genome.genes[2][1].binsignals)
    fake_gene = (; rnas = Any[(; expression = [1.0, 2.0])])
    @test has_expr(fake_gene)
end
