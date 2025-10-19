using Test
using BioinfoTools
using Interpolations
import Interpolations: AbstractInterpolation

Base.getindex(itp::AbstractInterpolation, r::AbstractVector{<:Real}) = [itp[x] for x in r]
Base.getindex(itp::AbstractInterpolation, r::AbstractRange{<:Real}) = [itp[x] for x in r]
Base.getindex(itp::AbstractInterpolation, x::Real) = itp(x)

const GT = BioinfoTools.GenomeTypes
const EU = BioinfoTools.EnrichmentUtils

function build_gene(; strand='+', region_len=800, gene_start=400, gene_end=520, binsignal::BitVector=BitVector(fill(false, region_len)))
    scaffold = GT.Scaffold("chr1", GT.Feature[], GT.Feature[], GT.Feature[], 1, region_len, "chromosome")
    signal_collection = nothing
    bin_collection = binsignal === nothing ? nothing : [binsignal]
    region = GT.Region(scaffold, missing, 1, region_len, missing, signal_collection, bin_collection, ["sample"])
    GT.Gene(
        scaffold,
        missing,
        "gene",
        "gene",
        strand,
        missing,
        gene_start,
        gene_end,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        [region],
        missing,
        gene_start,
        gene_end,
        nothing,
        nothing,
        nothing,
        ["sample"],
    )
end

@testset "EnrichmentUtils.getrange" begin
    @testset "Upstream positive strand" begin
        gene = build_gene(strand='+', gene_start=400, binsignal=BitVector(fill(true, 800)))
        range = EU.getrange(gene, "upstream")
        @test range == 1:149
    end

    @testset "Downstream negative strand" begin
        gene = build_gene(strand='-', gene_start=250, gene_end=360, binsignal=BitVector(fill(true, 800)))
        range = EU.getrange(gene, "downstream")
        @test range == 1:(gene.gene_start - gene.regions[1].region_start)
    end

    @testset "Promoter window" begin
        gene = build_gene(strand='+', gene_start=450, binsignal=BitVector(fill(true, 800)))
        range = EU.getrange(gene, "promoter")
        @test minimum(range) == (gene.gene_start - gene.regions[1].region_start) - 250
        @test maximum(range) == (gene.gene_start - gene.regions[1].region_start) + 150
    end

    @testset "Invalid range type" begin
        gene = build_gene(binsignal=BitVector(fill(true, 800)))
        @test_throws ErrorException EU.getrange(gene, "invalid")
    end
end

@testset "EnrichmentUtils.to_percent" begin
    @testset "Resample linear signal" begin
        signal = collect(0.0:99.0)
        percent = EU.to_percent(signal)
        @test length(percent) == 100
        @test percent == signal
    end

    @testset "Constant signal remains constant" begin
        signal = fill(5.5, 100)
        percent = EU.to_percent(signal)
        @test percent == signal
    end
end

@testset "EnrichmentUtils.getsiginrange" begin
    @testset "Positive strand slice" begin
        bins = BitVector([i <= 60 for i in 1:800])
        gene = build_gene(strand='+', gene_start=300, gene_end=420, binsignal=bins)
        range = GT.GeneRange(GT.TSS(), GT.TSS(), 0, 20, GT.INTEGER(), GT.INTEGER())
        result = EU.getsiginrange(gene, range, 1; peak_data=true, clamped=false)
        start = gene.gene_start - gene.regions[1].region_start + 1
        @test result == bins[start:start + 20]
    end

    @testset "Negative strand reversal" begin
        bins = BitVector([i <= 40 for i in 1:800])
        gene = build_gene(strand='-', gene_start=320, gene_end=440, binsignal=bins)
        range = GT.GeneRange(GT.TSS(), GT.TSS(), 0, 15, GT.INTEGER(), GT.INTEGER())
        result = EU.getsiginrange(gene, range, 1; peak_data=true, clamped=false)
        tss = (gene.regions[1].region_end - gene.gene_end) + 1
        expected = reverse(bins[tss:tss + 15])
        @test result == expected
    end
    
    @testset "Out-of-bounds returns missing" begin
        bins = BitVector(fill(true, 40))
        gene = build_gene(strand='+', region_len=40, gene_start=15, gene_end=30, binsignal=bins)
        range = GT.GeneRange(GT.TES(), GT.TES(), 0, 50, GT.INTEGER(), GT.INTEGER())
        @test EU.getsiginrange(gene, range, 1; peak_data=true, clamped=false) === missing
    end

    @testset "Percentage offsets" begin
        bins = BitVector([mod(i, 3) == 0 for i in 1:500])
        gene = build_gene(strand='+', region_len=500, gene_start=120, gene_end=260, binsignal=bins)
        gene_len = gene.gene_end - gene.gene_start + 1
        tss = gene.gene_start - gene.regions[1].region_start + 1

        range = GT.GeneRange(GT.TSS(), GT.TSS(), 25, 40, GT.PERCENTAGE(), GT.INTEGER())
        result = EU.getsiginrange(gene, range, 1; peak_data=true, clamped=false)

        start_offset = min(round(Int, (gene_len / 100) * 25), gene_len)
        expected_start = tss + start_offset
        expected_stop = tss + 40
        @test result == bins[expected_start:expected_stop]
    end
end

@testset "EnrichmentUtils.siginrange" begin
    @testset "Range within bounds" begin
        bins = BitVector(fill(true, 200))
        gene = build_gene(strand='+', gene_start=100, gene_end=180, binsignal=bins)
        range = GT.GeneRange(GT.TSS(), GT.TSS(), -5, 5, GT.INTEGER(), GT.INTEGER())
        @test EU.siginrange(gene, range, 1; peak_data=true)
    end

    @testset "Range exceeding bounds" begin
        bins = BitVector(fill(true, 60))
        gene = build_gene(strand='+', region_len=60, gene_start=20, gene_end=40, binsignal=bins)
        range = GT.GeneRange(GT.REGION(), GT.REGION(), -20, 200, GT.INTEGER(), GT.INTEGER())
        @test !EU.siginrange(gene, range, 1; peak_data=true)
    end

    @testset "Percentage offsets" begin
        bins = BitVector([i <= 120 for i in 1:400])
        gene = build_gene(strand='+', region_len=400, gene_start=50, gene_end=170, binsignal=bins)
        range_valid = GT.GeneRange(GT.TSS(), GT.TES(), 10, 10, GT.PERCENTAGE(), GT.INTEGER())
        @test EU.siginrange(gene, range_valid, 1; peak_data=true)

        range_invalid = GT.GeneRange(GT.TSS(), GT.TES(), 90, 500, GT.PERCENTAGE(), GT.INTEGER())
        @test !EU.siginrange(gene, range_invalid, 1; peak_data=true)
    end
end
