using Test
using BioinfoTools
using DataFrames
using GFF3

const LG = BioinfoTools.LoadGFF
const GT = BioinfoTools.GenomeTypes

@testset "LoadGFF" begin
    
    @testset "get_id - String parsing" begin
        # Test ID with gene: prefix
        @test LG.get_id("gene:ENSG00012345") == "ENSG00012345"
        
        # Test ID with Gene: prefix  
        @test LG.get_id("Gene:WBGene00001234") == "WBGene00001234"
        
        # Test ID with transcript: prefix
        @test LG.get_id("transcript:ENST00012345") == "ENST00012345"
        
        # Test ID with hyphen separator
        @test LG.get_id("gene-CELE_2L52.1") == "CELE_2L52.1"
        
        # Test ID without prefix (should return as-is)
        @test LG.get_id("GENE0001") == "GENE0001"
        
        # Test ID with multiple separators
        @test LG.get_id("gene:ENSG-12345-AB") == "ENSG-12345-AB"
    end
    
    @testset "get_id - Dictionary lookup" begin
        # Test successful ID retrieval with default "ID" field
        dict = Dict("ID" => ["gene:TEST001"], "Name" => ["TestGene"])
        @test LG.get_id("ID", dict) == "TEST001"
        
        # Test successful ID retrieval with custom field
        dict = Dict("gene_id" => ["gene-CUSTOM123"], "Name" => ["CustomGene"])
        @test LG.get_id("gene_id", dict) == "CUSTOM123"
        
        # Test missing ID field returns nothing
        dict = Dict("Name" => ["GeneWithoutID"])
        @test_logs (:warn, r"ID not found") LG.get_id("ID", dict) === nothing
        
        # Test missing custom field returns nothing
        dict = Dict("ID" => ["gene:TEST002"])
        @test_logs (:warn, r"'custom_field' not found") LG.get_id("custom_field", dict) === nothing
        
        # Test ID field with complex prefix
        dict = Dict("ID" => ["transcript:ENST-00012345-1"])
        @test LG.get_id("ID", dict) == "ENST-00012345-1"
    end
    
    @testset "parseattributes" begin
        # Test basic attribute parsing
        attrs = [
            "ID" => ["gene-123"],
            "Name" => ["TestGene"],
            "Dbxref" => ["GeneID:456", "WormBase:WB789"]
        ]
        result = LG.parseattributes(attrs)
        @test result isa Dict{String, Vector{String}}
        @test result["ID"] == ["gene-123"]
        @test result["Name"] == ["TestGene"]
        @test result["Dbxref"] == ["GeneID:456", "WormBase:WB789"]
        
        # Test empty attributes
        empty_attrs = Pair{String, Vector{String}}[]
        result_empty = LG.parseattributes(empty_attrs)
        @test result_empty isa Dict{String, Vector{String}}
        @test isempty(result_empty)
        
        # Test single attribute
        single_attr = ["ID" => ["single-gene"]]
        result_single = LG.parseattributes(single_attr)
        @test length(result_single) == 1
        @test result_single["ID"] == ["single-gene"]
    end
    
    @testset "create_region" begin
        # Test positive strand region creation
        scaffold = GT.Scaffold("chr1", missing, GT.Gene[], GT.Repeat[], missing, missing, missing)
        region_pos = LG.create_region('+', 1000, 2000, 500, 300, 10000, scaffold, missing)
        @test region_pos isa GT.Region
        @test region_pos.region_start == 500  # 1000 - 500
        @test region_pos.region_end == 2300   # 2000 + 300
        @test region_pos.scaffold === scaffold
        
        # Test negative strand region creation
        region_neg = LG.create_region('-', 1000, 2000, 500, 300, 10000, scaffold, missing)
        @test region_neg.region_start == 700   # 1000 - 300
        @test region_neg.region_end == 2500    # 2000 + 500
        
        # Test boundary conditions - start at chromosome beginning
        region_start = LG.create_region('+', 100, 500, 500, 100, 10000, scaffold, missing)
        @test region_start.region_start == 1   # max(100-500, 1) = 1
        
        # Test boundary conditions - end at chromosome end
        region_end = LG.create_region('+', 9000, 9500, 100, 1000, 10000, scaffold, missing)
        @test region_end.region_end == 10000   # min(9500+1000, 10000) = 10000
        
        # Test with missing scaffold
        region_missing = LG.create_region('+', 1000, 2000, 100, 100, 10000, missing, missing)
        @test ismissing(region_missing.scaffold)
    end
    
    @testset "getchromlengths" begin
        gff_file = joinpath(@__DIR__, "data", "NC_003280.10_copy.gff")
        
        # Test successful chromosome length extraction
        chrom_lengths = LG.getchromlengths(gff_file)
        @test chrom_lengths isa DataFrame
        @test hasproperty(chrom_lengths, :chrom)
        @test hasproperty(chrom_lengths, :length)
        @test nrow(chrom_lengths) == 1
        @test chrom_lengths[1, :chrom] == "NC_003280.10"
        @test chrom_lengths[1, :length] > 0
        
        # Test with gzipped file
        gff_gz_file = joinpath(@__DIR__, "data", "NC_003280.10.gff.gz")
        if isfile(gff_gz_file)
            chrom_lengths_gz = LG.getchromlengths(gff_gz_file)
            @test chrom_lengths_gz isa DataFrame
            @test nrow(chrom_lengths_gz) >= 1
        end
    end
    
    @testset "averagereplicates!" begin
        # Test basic replicate averaging
        expr_data = DataFrame(
            GeneID = ["GENE1", "GENE2", "GENE3"],
            Length = [1000, 1500, 2000],
            Rep1_A = [10.0, 20.0, 30.0],
            Rep1_B = [12.0, 22.0, 32.0],
            Rep2_A = [15.0, 25.0, 35.0],
            Rep2_B = [13.0, 23.0, 33.0]
        )
        
        LG.averagereplicates!(expr_data)
        
        # Check structure after averaging
        @test ncol(expr_data) == 4  # GeneID, Length, Rep1_A_Rep1_B, Rep2_A_Rep2_B
        @test "GeneID" in names(expr_data)
        @test "Length" in names(expr_data)
        @test expr_data[1, 3] == 11.0  # mean of 10 and 12
        @test expr_data[2, 3] == 21.0  # mean of 20 and 22
        
        # Test error with odd number of replicate columns
        expr_data_odd = DataFrame(
            GeneID = ["GENE1"],
            Length = [1000],
            Rep1 = [10.0],
            Rep2 = [20.0],
            Rep3 = [30.0]
        )
        @test_throws Exception LG.averagereplicates!(expr_data_odd)
    end
    
    @testset "parsegene!" begin
        # Setup
        gff_file = joinpath(@__DIR__, "data", "NC_003280.10_copy.gff")
        reader = open(GFF3.Reader, gff_file)
        chrom_lengths = DataFrame(:chrom => ["NC_003280.10"], :length => [15279421])
        refs = GT.RefGenome()
        
        # Find a gene record
        gene_record = nothing
        for record in reader
            if GFF3.featuretype(record) == "gene"
                gene_record = record
                break
            end
        end
        close(reader)
        
        @test !isnothing(gene_record)
        
        # Test successful gene parsing
        LG.parsegene!(gene_record, refs, chrom_lengths, false)
        
        @test length(refs.genes[1]) == 1
        @test length(refs.genes[2]) == 1
        gene = refs.genes[2][1]
        @test gene isa GT.Gene
        @test gene.gene_start > 0
        @test gene.gene_end > gene.gene_start
        @test gene.strand in ['+', '-']
        @test !isempty(gene.regions)
        @test haskey(refs.scaffolds, "NC_003280.10")
        
        # Test with alternative ID field
        refs2 = GT.RefGenome()
        reader2 = open(GFF3.Reader, gff_file)
        gene_record2 = nothing
        for record in reader2
            if GFF3.featuretype(record) == "gene"
                gene_record2 = record
                break
            end
        end
        close(reader2)
        
        LG.parsegene!(gene_record2, refs2, chrom_lengths, false; alt_id_field="Name")
        @test length(refs2.genes[1]) == 1
        
        # Test error with chromosome not in chrom_lengths
        bad_chrom_lengths = DataFrame(:chrom => ["chr99"], :length => [1000])
        refs3 = GT.RefGenome()
        @test_throws ErrorException LG.parsegene!(gene_record, refs3, bad_chrom_lengths, false)
    end
    
    @testset "parserna!" begin
        # Setup - load a gene first
        gff_file = joinpath(@__DIR__, "data", "NC_003280.10_copy.gff")
        reader = open(GFF3.Reader, gff_file)
        chrom_lengths = DataFrame(:chrom => ["NC_003280.10"], :length => [15279421])
        refs = GT.RefGenome()
        
        gene_record = nothing
        rna_record = nothing
        
        for record in reader
            if GFF3.featuretype(record) == "gene" && isnothing(gene_record)
                gene_record = record
                LG.parsegene!(gene_record, refs, chrom_lengths, false)
            elseif contains(GFF3.featuretype(record), "RNA") && isnothing(rna_record)
                rna_record = record
                break
            end
        end
        close(reader)
        
        @test !isnothing(rna_record)
        @test length(refs.genes[2]) == 1
        
        prev_gene = refs.genes[2][1]
        
        # Test successful RNA parsing with valid parent
        LG.parserna!(rna_record, refs, prev_gene)
        
        @test length(refs.rnas) == 1
        rna = refs.rnas[1]
        @test rna isa GT.RNA
        @test !ismissing(rna.id)
        @test rna.rna_start > 0
        @test rna.rna_end >= rna.rna_start
        @test length(prev_gene.rnas) == 1
        @test prev_gene.rnas[1] === rna
        
        # Test with missing parent gene (warning case)
        refs2 = GT.RefGenome()
        @test_logs (:warn, r"Parent gene.*not found") LG.parserna!(rna_record, refs2, nothing)
        @test length(refs2.rnas) == 1
        @test ismissing(refs2.rnas[1].gene)
    end
    
    @testset "parseexon!" begin
        # Setup - load gene and RNA first
        gff_file = joinpath(@__DIR__, "data", "NC_003280.10_copy.gff")
        reader = open(GFF3.Reader, gff_file)
        chrom_lengths = DataFrame(:chrom => ["NC_003280.10"], :length => [15279421])
        refs = GT.RefGenome()
        
        gene_record = nothing
        rna_record = nothing
        exon_record = nothing
        
        for record in reader
            if GFF3.featuretype(record) == "gene" && isnothing(gene_record)
                gene_record = record
                LG.parsegene!(gene_record, refs, chrom_lengths, false)
            elseif contains(GFF3.featuretype(record), "RNA") && isnothing(rna_record)
                rna_record = record
                LG.parserna!(rna_record, refs, refs.genes[2][1])
            elseif GFF3.featuretype(record) == "exon" && isnothing(exon_record)
                exon_record = record
                break
            end
        end
        close(reader)
        
        @test !isnothing(exon_record)
        prev_gene = refs.genes[2][1]
        prev_rna = refs.rnas[1]
        
        # Test successful exon parsing
        LG.parseexon!(exon_record, refs, prev_rna, prev_gene)
        
        @test length(refs.exons) == 1
        exon = refs.exons[1]
        @test exon isa GT.Exon
        @test exon.gene === prev_gene
        @test exon.rna === prev_rna
        @test exon.exon_start > 0
        @test exon.exon_end >= exon.exon_start
        @test length(prev_gene.exons) == 1
        @test length(prev_rna.exons) == 1
        
        # Test with missing parent RNA (warning case)
        refs2 = GT.RefGenome()
        @test_logs (:warn, r"could not find parent") LG.parseexon!(exon_record, refs2, nothing, prev_gene) === nothing
        @test length(refs2.exons) == 0
        
        # Test with missing parent gene (warning case)
        refs3 = GT.RefGenome()
        @test_logs (:warn, r"could not find parent") LG.parseexon!(exon_record, refs3, prev_rna, nothing) === nothing
        @test length(refs3.exons) == 0
    end
    
    @testset "loadgff" begin
        gff_file = joinpath(@__DIR__, "data", "NC_003280.10_copy.gff")
        
        # Test loading genes only
        genome = LG.loadgff(gff_file; feature_type="gene")
        @test genome isa GT.RefGenome
        @test length(genome.genes[1]) > 0
        @test length(genome.genes[2]) > 0
        @test haskey(genome.scaffolds, "NC_003280.10")
        
        # Test loading all features
        genome_all = LG.loadgff(gff_file; feature_type="all")
        @test genome_all isa GT.RefGenome
        @test length(genome_all.genes[1]) > 0
        @test length(genome_all.rnas) > 0
        @test length(genome_all.exons) > 0
        
        # Test with existing genome (should return nothing)
        existing_genome = GT.RefGenome()
        result = LG.loadgff(gff_file; feature_type="gene", genome=existing_genome)
        @test result === nothing
        @test length(existing_genome.genes[1]) > 0
        
        # Test with chromosome lengths file
        chrom_file = joinpath(@__DIR__, "data", "test_chrom_lengths.txt")
        open(chrom_file, "w") do io
            println(io, "NC_003280.10\t15279421")
        end
        
        genome_with_chrom = LG.loadgff(gff_file, chrom_file; feature_type="gene")
        @test genome_with_chrom isa GT.RefGenome
        @test length(genome_with_chrom.genes[1]) > 0
        
        rm(chrom_file)
    end
    
    @testset "loadgenome" begin
        gff_file = joinpath(@__DIR__, "data", "NC_003280.10_copy.gff")
        
        # Test with single file path
        genome = LG.loadgenome(gff_file; feature_type="gene")
        @test genome isa GT.RefGenome
        @test length(genome.genes[1]) > 0
        
        # Test with vector of file paths
        genome_vec = LG.loadgenome([gff_file]; feature_type="gene")
        @test genome_vec isa GT.RefGenome
        @test length(genome_vec.genes[1]) > 0
        
        # Test with directory path
        data_dir = joinpath(@__DIR__, "data")
        genome_dir = LG.loadgenome(data_dir; feature_type="gene")
        @test genome_dir isa GT.RefGenome
        @test length(genome_dir.genes[1]) > 0
        
        # Test with alternative ID field
        genome_alt = LG.loadgenome(gff_file; feature_type="gene", alt_id_field="Name")
        @test genome_alt isa GT.RefGenome
        @test length(genome_alt.genes[1]) > 0
    end
    
    @testset "addsequence!" begin
        using BioSequences
        
        # Create a test genome and repeat
        genome = GT.RefGenome()
        scaffold = GT.Scaffold("chr1", missing, GT.Gene[], GT.Repeat[], missing, missing, missing)
        genome.scaffolds["chr1"] = scaffold
        
        repeat_elem = GT.Repeat(
            scaffold,
            missing,
            GT.Region[],
            "LINE",
            "L1",
            100,
            200,
            missing,
            nothing,
            nothing,
            nothing
        )
        push!(genome.repeats, repeat_elem)
        
        # Add sequence to repeat
        seq = dna"ACGTACGTACGT"
        LG.addsequence!(genome, repeat_elem, seq)
        
        # Check that the repeat was replaced with updated sequence
        @test length(genome.repeats) == 1
        @test !ismissing(genome.repeats[1].sequence)
        @test genome.repeats[1].sequence == seq
        @test genome.repeats[1].family == "LINE"
        @test genome.repeats[1].type == "L1"
    end
    
    @testset "Integration test - Full workflow" begin
        # Test the complete workflow of loading a GFF file
        gff_file = joinpath(@__DIR__, "data", "NC_003280.10_copy.gff")
        
        # Load genome with all features
        genome = LG.loadgenome(gff_file; feature_type="all")
        
        # Verify genome structure
        @test genome isa GT.RefGenome
        @test length(genome.genes[1]) > 0
        @test length(genome.genes[2]) > 0
        @test length(genome.rnas) > 0
        @test length(genome.exons) > 0
        
        # Verify gene structure
        first_gene = genome.genes[2][1]
        @test first_gene.gene_start > 0
        @test first_gene.gene_end > first_gene.gene_start
        @test first_gene.strand in ['+', '-']
        
        # Verify RNA-gene relationships
        if length(first_gene.rnas) > 0
            first_rna = first_gene.rnas[1]
            @test first_rna.gene === first_gene
        end
        
        # Verify exon-RNA-gene relationships
        if length(genome.exons) > 0
            first_exon = genome.exons[1]
            @test !ismissing(first_exon.gene)
            @test !ismissing(first_exon.rna)
        end
        
        # Verify scaffold structure
        @test haskey(genome.scaffolds, "NC_003280.10")
        scaffold = genome.scaffolds["NC_003280.10"]
        @test length(scaffold.genes) > 0
        @test scaffold.genes[1] === first_gene
    end
end
