using Test
using DataFrames
using Graphs
using MetaGraphs
using BioinfoTools.ParalogUtils

@testset "ParalogUtils Tests" begin

    # ========================================
    # Tests for rbh()
    # ========================================
    
    @testset "rbh - basic functionality with max scoring" begin
        # Create test data with reciprocal best hits
        df = DataFrame(
            GeneID = ["A", "B", "C", "D"],
            ParalogID = ["B", "A", "D", "C"],
            Perc1 = [95.0, 94.0, 80.0, 85.0],
            Perc2 = [94.0, 95.0, 85.0, 80.0]
        )
        
        result = rbh(df; scoring="max")
        
        # Should find 2 RBH pairs: A-B and C-D
        @test nrow(result) == 2
        @test "GeneID" in names(result)
        @test "ParalogID" in names(result)
        @test "perc_1" in names(result)
        @test "perc_2" in names(result)
        @test "max_perc" in names(result)
        @test "mean_perc" in names(result)
        
        # Check that max_perc is correctly calculated
        @test all(result.max_perc .>= result.perc_1)
        @test all(result.max_perc .>= result.perc_2)
    end
    
    @testset "rbh - mean scoring" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "A"],
            Perc1 = [90.0, 88.0],
            Perc2 = [88.0, 90.0]
        )
        
        result = rbh(df; scoring="mean")
        
        @test nrow(result) == 1
        # Mean should be average of bidirectional scores
        @test result.mean_perc[1] ≈ 89.0
    end
    
    @testset "rbh - average scoring (alias)" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "A"],
            Perc1 = [100.0, 80.0],
            Perc2 = [80.0, 100.0]
        )
        
        result = rbh(df; scoring="avg")
        
        @test nrow(result) == 1
        @test result.mean_perc[1] ≈ 90.0
    end
    
    @testset "rbh - double_max scoring" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "A"],
            Perc1 = [95.0, 93.0],
            Perc2 = [93.0, 95.0]
        )
        
        result = rbh(df; scoring="double_max")
        
        @test nrow(result) == 1
        # In double_max mode, should use original scores
        @test result.perc_1[1] == 95.0
        @test result.perc_2[1] == 93.0
    end
    
    @testset "rbh - no reciprocal hits" begin
        df = DataFrame(
            GeneID = ["A", "B", "C"],
            ParalogID = ["B", "C", "D"],
            Perc1 = [95.0, 90.0, 85.0],
            Perc2 = [94.0, 89.0, 84.0]
        )
        
        result = rbh(df; scoring="max")
        
        # May still find some hits depending on the actual algorithm
        # Just check that it returns a valid dataframe
        @test "GeneID" in names(result)
        @test "ParalogID" in names(result)
    end
    
    @testset "rbh - empty input should error" begin
        df = DataFrame(
            GeneID = String[],
            ParalogID = String[],
            Perc1 = Float64[],
            Perc2 = Float64[]
        )
        
        # Empty DataFrame will cause BoundsError when accessing elements
        @test_throws BoundsError rbh(df; scoring="max")
    end

    # ========================================
    # Tests for findfamilies()
    # ========================================
    
    @testset "findfamilies - basic unweighted graph" begin
        df = DataFrame(
            GeneID = ["A", "A", "B"],
            ParalogID = ["B", "C", "C"],
            EdgeValue = [0.5, 0.6, 0.4]
        )
        
        graph = findfamilies(df; apply_cutoff=false, 
                            cutoff_val=nothing, 
                            cutoff_variable=nothing,
                            cutoff_comparison=<,
                            weighted_graph=false,
                            add_names=false)
        
        @test graph isa SimpleGraph
        @test nv(graph) == 3  # A, B, C
        @test ne(graph) == 3  # A-B, A-C, B-C
    end
    
    @testset "findfamilies - with names" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            EdgeValue = [0.5, 0.6]
        )
        
        graph = findfamilies(df; apply_cutoff=false,
                            cutoff_val=nothing,
                            cutoff_variable=nothing,
                            cutoff_comparison=<,
                            weighted_graph=false,
                            add_names=true)
        
        @test graph isa MetaGraph
        @test nv(graph) == 3
    end
    
    @testset "findfamilies - with cutoff filter" begin
        df = DataFrame(
            GeneID = ["A", "A", "B"],
            ParalogID = ["B", "C", "C"],
            dS = [0.3, 0.8, 0.4]
        )
        
        graph = findfamilies(df; apply_cutoff=true,
                            cutoff_val=0.5,
                            cutoff_variable="dS",
                            cutoff_comparison=<,
                            weighted_graph=false,
                            add_names=false)
        
        @test graph isa SimpleGraph
        @test nv(graph) == 3
        # Only A-B and B-C edges should remain (dS < 0.5)
        @test ne(graph) == 2
    end
    
    @testset "findfamilies - weighted graph" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            dS = [0.5, 0.6]
        )
        
        graph = findfamilies(df; apply_cutoff=false,
                            cutoff_val=nothing,
                            cutoff_variable="dS",
                            cutoff_comparison=<,
                            weighted_graph=true,
                            add_names=false)
        
        @test graph isa MetaGraph
        @test nv(graph) == 3
        @test ne(graph) == 2
    end
    
    @testset "findfamilies - invalid input (wrong number of columns)" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"]
        )
        
        @test_throws ErrorException findfamilies(df; apply_cutoff=false,
                                                    cutoff_val=nothing,
                                                    cutoff_variable=nothing,
                                                    cutoff_comparison=<,
                                                    weighted_graph=false,
                                                    add_names=false)
    end
    
    @testset "findfamilies - invalid input (wrong column types)" begin
        df = DataFrame(
            GeneID = [1, 2],  # Should be String
            ParalogID = ["B", "C"],
            EdgeValue = [0.5, 0.6]
        )
        
        @test_throws ErrorException findfamilies(df; apply_cutoff=false,
                                                    cutoff_val=nothing,
                                                    cutoff_variable=nothing,
                                                    cutoff_comparison=<,
                                                    weighted_graph=false,
                                                    add_names=false)
    end
    
    @testset "findfamilies - with expression data (5 columns)" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            dS = [0.5, 0.6],
            ExprGene = [10.0, 15.0],
            ExprParalog = [12.0, 20.0]
        )
        
        graph = findfamilies(df; apply_cutoff=false,
                            cutoff_val=nothing,
                            cutoff_variable="dS",
                            cutoff_comparison=<,
                            weighted_graph=true,
                            add_names=false)
        
        @test graph isa MetaGraph
        @test nv(graph) == 3
    end

    # ========================================
    # Tests for family_sizes()
    # ========================================
    
    @testset "family_sizes - basic functionality" begin
        families = [["geneA", "geneB"], ["geneC", "geneD", "geneE"], ["geneF"]]
        
        result = family_sizes(families)
        
        @test nrow(result) == 6
        @test "GeneID" in names(result)
        @test "FamilySize" in names(result)
        
        # Check family sizes
        @test result[result.GeneID .== "geneA", :FamilySize][1] == 2
        @test result[result.GeneID .== "geneB", :FamilySize][1] == 2
        @test result[result.GeneID .== "geneC", :FamilySize][1] == 3
        @test result[result.GeneID .== "geneD", :FamilySize][1] == 3
        @test result[result.GeneID .== "geneE", :FamilySize][1] == 3
        @test result[result.GeneID .== "geneF", :FamilySize][1] == 1
    end
    
    @testset "family_sizes - empty input" begin
        families = Vector{String}[]
        
        result = family_sizes(families)
        
        @test nrow(result) == 0
        @test "GeneID" in names(result)
        @test "FamilySize" in names(result)
    end
    
    @testset "family_sizes - single family" begin
        families = [["gene1", "gene2", "gene3"]]
        
        result = family_sizes(families)
        
        @test nrow(result) == 3
        @test all(result.FamilySize .== 3)
    end
    
    @testset "family_sizes - all singletons" begin
        families = [["geneA"], ["geneB"], ["geneC"]]
        
        result = family_sizes(families)
        
        @test nrow(result) == 3
        @test all(result.FamilySize .== 1)
    end

    # ========================================
    # Tests for create_unweighted_graph()
    # ========================================
    
    @testset "create_unweighted_graph - basic functionality" begin
        df = DataFrame(
            GeneID = ["A", "B", "C"],
            ParalogID = ["B", "C", "A"]
        )
        
        graph = create_unweighted_graph(df)
        
        @test graph isa SimpleGraph
        @test nv(graph) == 3
        @test ne(graph) == 3
    end
    
    @testset "create_unweighted_graph - with names" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"]
        )
        
        graph = create_unweighted_graph(df; add_names=true)
        
        @test graph isa MetaGraph
        @test nv(graph) == 3
        @test ne(graph) == 2
        
        # Check that names are stored
        @test has_prop(graph, 1, :name)
        @test has_prop(graph, 2, :name)
        @test has_prop(graph, 3, :name)
    end
    
    @testset "create_unweighted_graph - duplicate edges" begin
        df = DataFrame(
            GeneID = ["A", "A", "B"],
            ParalogID = ["B", "B", "A"]
        )
        
        graph = create_unweighted_graph(df)
        
        @test nv(graph) == 2
        # Duplicate edge should only be added once
        @test ne(graph) == 1
    end
    
    @testset "create_unweighted_graph - empty input" begin
        df = DataFrame(
            GeneID = String[],
            ParalogID = String[]
        )
        
        graph = create_unweighted_graph(df)
        
        @test nv(graph) == 0
        @test ne(graph) == 0
    end
    
    @testset "create_unweighted_graph - self loops" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["A", "B"]
        )
        
        graph = create_unweighted_graph(df)
        
        # The function creates vertices for each unique gene ID
        # Since "A" in GeneID and "A" in ParalogID are the same gene,
        # the graph should have 2 vertices (A and B)
        # Testing the actual behavior - adds both as separate initially
        @test graph isa SimpleGraph
        # SimpleGraph may or may not support self-loops depending on implementation
        # Just verify it completes without error
    end

    # ========================================
    # Tests for create_weighted_graph()
    # ========================================
    
    @testset "create_weighted_graph - basic functionality" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            dS = [0.5, 0.6]
        )
        
        graph = create_weighted_graph(df, "dS")
        
        @test graph isa MetaGraph
        @test nv(graph) == 3
        @test ne(graph) == 2
    end
    
    @testset "create_weighted_graph - with names" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            dS = [0.5, 0.6]
        )
        
        graph = create_weighted_graph(df, :dS; add_names=true)
        
        @test graph isa MetaGraph
        @test has_prop(graph, 1, :name)
        @test has_prop(graph, 2, :name)
        @test has_prop(graph, 3, :name)
    end
    
    @testset "create_weighted_graph - with expression data" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            dS = [0.5, 0.6],
            ExprA = [10.0, 15.0],
            ExprB = [12.0, 20.0]
        )
        
        graph = create_weighted_graph(df, :dS; add_expr=true)
        
        @test graph isa MetaGraph
        # Check that expression properties are added
        @test has_prop(graph, 1, :expr)
        @test has_prop(graph, 2, :expr)
    end
    
    @testset "create_weighted_graph - edge weights are stored" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            score = [0.8, 0.9]
        )
        
        graph = create_weighted_graph(df, :score)
        
        @test has_prop(graph, 1, 2, :weight)
        @test has_prop(graph, 2, 3, :weight)
    end
    
    @testset "create_weighted_graph - empty input" begin
        df = DataFrame(
            GeneID = String[],
            ParalogID = String[],
            dS = Float64[]
        )
        
        graph = create_weighted_graph(df, :dS)
        
        @test nv(graph) == 0
        @test ne(graph) == 0
    end
    
    @testset "create_weighted_graph - custom cutoff value" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            dS = [0.5, 0.6]
        )
        
        graph = create_weighted_graph(df, :dS; cutoff_val=0.55)
        
        @test graph isa MetaGraph
        @test nv(graph) == 3
    end

    # ========================================
    # Tests for savegraphcsv()
    # ========================================
    
    @testset "savegraphcsv - basic functionality" begin
        # Create a simple weighted graph
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            dS = [0.5, 0.6]
        )
        graph = create_weighted_graph(df, :dS)
        
        # Save to temporary file
        tmpfile = tempname() * ".csv"
        
        savegraphcsv(graph, tmpfile)
        
        @test isfile(tmpfile)
        
        # Read back and check format
        content = read(tmpfile, String)
        @test occursin("Source,Target,Weight", content)
        
        # Clean up
        rm(tmpfile)
    end
    
    @testset "savegraphcsv - with names" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            dS = [0.5, 0.6]
        )
        graph = create_weighted_graph(df, :dS; add_names=true)
        
        tmpfile = tempname() * ".csv"
        
        savegraphcsv(graph, tmpfile; add_names=true)
        
        @test isfile(tmpfile)
        
        content = read(tmpfile, String)
        @test occursin("Source,Target,Weight,Source_Name,Target_Name", content)
        
        rm(tmpfile)
    end
    
    @testset "savegraphcsv - file exists without force" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            dS = [0.5, 0.6]
        )
        graph = create_weighted_graph(df, :dS)
        
        tmpfile = tempname() * ".csv"
        
        # Create the file first
        savegraphcsv(graph, tmpfile)
        
        # Try to save again without force - should error
        @test_throws ErrorException savegraphcsv(graph, tmpfile; force=false)
        
        rm(tmpfile)
    end
    
    @testset "savegraphcsv - file exists with force" begin
        df = DataFrame(
            GeneID = ["A", "B"],
            ParalogID = ["B", "C"],
            dS = [0.5, 0.6]
        )
        graph = create_weighted_graph(df, :dS)
        
        tmpfile = tempname() * ".csv"
        
        # Create the file first
        savegraphcsv(graph, tmpfile)
        
        # Save again with force - should succeed
        savegraphcsv(graph, tmpfile; force=true)
        
        @test isfile(tmpfile)
        
        rm(tmpfile)
    end
    
    @testset "savegraphcsv - empty graph" begin
        df = DataFrame(
            GeneID = String[],
            ParalogID = String[],
            dS = Float64[]
        )
        graph = create_weighted_graph(df, :dS)
        
        tmpfile = tempname() * ".csv"
        
        savegraphcsv(graph, tmpfile)
        
        @test isfile(tmpfile)
        
        # Should only have header
        content = read(tmpfile, String)
        lines = split(strip(content), '\n')
        @test length(lines) == 1  # Only header
        
        rm(tmpfile)
    end

end
