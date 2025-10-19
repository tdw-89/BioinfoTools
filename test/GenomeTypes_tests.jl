using Test
using BioSequences
using BioinfoTools

const GT = BioinfoTools.GenomeTypes

@testset "GenomeTypes" begin

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

    @testset "Reference Genome" begin
        scaffold = GT.Scaffold("chrR", GT.Feature[], missing, missing, 1, 1_000, "chromosome")
        gene = GT.Gene(scaffold, missing, "GeneRef", "GENE_REF", '+')
        genome = GT.RefGenome(
            Dict("chrR" => scaffold),
            Dict{String, GT.Contig}(),
            GT.Enhancer[],
            GT.Promoter[],
            (String["GENE_REF"], GT.Gene[gene]),
            GT.Intron[],
            GT.Exon[],
            GT.RNA[],
            Dict{String, Vector{GT.Segment}}(),
            GT.Repeat[],
            GT.Region[],
            GT.Annotation[],
        )

        @test haskey(genome.scaffolds, "chrR")
        @test genome.genes[1][1] == "GENE_REF"
        @test genome.genes[2][1] === gene
    end

    @testset "Gene Range Utilities" begin
        default_range = GT.GeneRange(GT.TSS(), GT.TES())
        @test default_range.start_offset == 0
        @test default_range.stop_offset == 0
        @test default_range.start_offset_type isa GT.INTEGER
        @test default_range.stop_offset_type isa GT.INTEGER

        percent_range = GT.GeneRange(GT.TSS(), GT.TSS(), 25, 10, GT.PERCENTAGE(), GT.INTEGER())
        @test percent_range.start_offset == 25
        @test percent_range.start_offset_type isa GT.PERCENTAGE

        @test_throws ErrorException GT.GeneRange(GT.REGION(), GT.TSS(), 10, 5, GT.PERCENTAGE(), GT.INTEGER())
    end

    @testset "Helper Methods" begin
        range = GT.GeneRange(GT.TSS(), GT.TSS(), -2, 4, GT.INTEGER(), GT.INTEGER())
        @test GT.to_vector(range) == -2:4
        @test_throws ErrorException GT.to_vector(GT.GeneRange(GT.TSS(), GT.TES()))

        scaffold = GT.Scaffold(
            "chrH",
            GT.Feature[],
            GT.Gene[],
            GT.Feature[],
            missing,
            missing,
            missing,
        )

        gene_a = GT.Gene(scaffold, missing, "A", "A", '+', missing, missing, missing,
                         missing, missing, missing, missing, missing, missing, missing,
                         missing, 200, 400)
        gene_b = GT.Gene(scaffold, missing, "B", "B", '-', missing, missing, missing,
                         missing, missing, missing, missing, missing, missing, missing,
                         missing, 150, 280)
        gene_c = GT.Gene(scaffold, missing, "C", "C", '+', missing, missing, missing,
                         missing, missing, missing, missing, missing, missing, missing,
                         missing, missing, 500)

        append!(scaffold.genes, [gene_a, gene_b, gene_c])

        genome = GT.RefGenome(
            Dict("chrH" => scaffold),
            Dict{String, GT.Contig}(),
            GT.Enhancer[],
            GT.Promoter[],
            (String["B", "A", "C"], GT.Gene[gene_b, gene_a, gene_c]),
            GT.Intron[],
            GT.Exon[],
            GT.RNA[],
            Dict{String, Vector{GT.Segment}}(),
            GT.Repeat[],
            GT.Region[],
            GT.Annotation[],
        )

        GT.sortgenes!(genome)
        sorted_genes = genome.scaffolds["chrH"].genes
        @test sorted_genes[1] === gene_b
        @test sorted_genes[2] === gene_a
        @test sorted_genes[end] === gene_c

        @test get(genome, "A") === gene_a
        @test get(genome, "missing-gene") === missing
        results = get(genome, ["A", "missing-gene", "B"])
        @test results[1] === gene_a
        @test results[2] === missing
        @test results[3] === gene_b

        @test GT.hasoverlap(0, 12, 10, 20) === false
        @test GT.hasoverlap(0, 5, 20, 25) === true
        @test GT.hasoverlap(gene_a, gene_b) === true

        @test GT.overlaplength(0, 5, 20, 25) == 16
        @test GT.overlaplength(0, 100, 10, 90) == 0
        @test GT.overlaplength(gene_a, gene_b) == 81
    end
end
