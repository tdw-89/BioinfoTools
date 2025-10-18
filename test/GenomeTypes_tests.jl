using Test
using BioSequences
using BioinfoTools

const GT = BioinfoTools.GenomeTypes

@testset "GenomeTypes" begin
    @testset "Type Hierarchy" begin
        @test isabstracttype(GT.Feature)
        @test isabstracttype(GT.RegElement)
        @test GT.RegElement <: GT.Feature
    end

    @testset "Scaffold and Contig Construction" begin
        scaffold = GT.Scaffold("scaf1", Vector{GT.Feature}(), missing, missing, 1, 1_000, "chromosome")
        @test scaffold.name == "scaf1"
        @test scaffold.scaffold_start == 1
        @test scaffold.level == "chromosome"

        contig = GT.Contig(scaffold, missing, missing, missing, 10, 900)
        @test contig.scaffold === scaffold
        @test contig.contig_start == 10
        @test contig.contig_end == 900

        push!(scaffold.contigs, contig)
        @test length(scaffold.contigs) == 1

        @test_throws MethodError GT.Scaffold(42, missing, missing, missing, missing, missing, missing)
    end

    @testset "Gene and Regulatory Elements" begin
        scaffold = GT.Scaffold("chrX", missing, missing, missing, missing, missing, missing)
        contig = GT.Contig(scaffold, missing, missing, missing, missing, missing)
        dna = dna"ACGTACGT"
        signals = [UInt16[1, 2, 3]]
        binsignals = [BitVector([true, false, true])]
        samples = ["sample1"]
        annotation = GT.Annotation(missing, "hypothetical")
        annotations = Dict("note" => Vector{GT.Feature}([annotation]))

        gene = GT.Gene(
            scaffold,
            contig,
            "Gene1",
            "GENE0001",
            '+',
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            annotations,
            100,
            900,
            dna,
            signals,
            binsignals,
            samples,
        )

        @test gene isa GT.Feature
        @test gene.sequence === dna
        @test gene.annotations === annotations
        @test gene.signals === signals
        @test gene.samples == samples

        promoter = GT.Promoter([gene], '+', missing, 80, 120, dna, signals, binsignals, samples)
        @test promoter.genes[1] === gene
        @test promoter.promoter_start == 80

        enhancer = GT.Enhancer(contig)
        @test enhancer.contig === contig
    end

    @testset "Transcript Features" begin
        scaffold = GT.Scaffold("chr1", missing, missing, missing, missing, missing, missing)
        contig = GT.Contig(scaffold, missing, missing, missing, 1, 1_000)
        gene = GT.Gene(
            scaffold,
            contig,
            "Gene2",
            "GENE0002",
            '-',
            missing,
            10,
            500,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            missing,
            10,
            500,
            nothing,
            nothing,
            nothing,
            ["sampleA"],
        )

        exon = GT.Exon(gene, missing, 10, 100)
        intron = GT.Intron(gene, missing, 101, 200)
        rna = GT.RNA(
            "RNA1",
            "mRNA",
            gene,
            [exon],
            [intron],
            ["sampleA"],
            10,
            500,
            [0.5, 1.0],
        )

        @test exon.gene === gene
        @test intron.intron_end == 200
        @test rna.exons[1] === exon
        @test rna.expression == [0.5, 1.0]
    end

    @testset "Intervals and Repeats" begin
        scaffold = GT.Scaffold("chr2", missing, missing, missing, missing, missing, missing)
        contig = GT.Contig(scaffold, missing, missing, missing, missing, missing)
        annotation = GT.Annotation(missing, "regulatory hotspot")
        region_annotations = Dict("default" => GT.Annotation[annotation])
        signals = [Float64[0.1, 0.2]]
        binsignals = [BitVector([false, true])]
        samples = ["sampleZ"]

        region = GT.Region(
            scaffold,
            contig,
            1,
            200,
            region_annotations,
            signals,
            binsignals,
            samples,
        )

        repeat = GT.Repeat(
            scaffold,
            contig,
            [region],
            "LINE",
            "L1",
            5,
            150,
            missing,
            signals,
            binsignals,
            samples,
        )

        segment = GT.Segment([region], 1, 200)

        @test region.annotations == region_annotations
        @test repeat.family == "LINE"
        @test repeat.regions[1] === region
        @test segment.prev_segment[1] === region
    end
end
