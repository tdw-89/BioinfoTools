module ParalogUtils

"""
Utilities for analyzing paralogous gene relationships.

This module provides functions for:
- Reciprocal best hit (RBH) analysis
- Gene family identification using graph-based methods
- Graph creation and manipulation for paralog relationships
- Export utilities for network analysis

# TODO
- Add a tie-breaking scheme to each scoring method
"""

using DataFrames
using Graphs
using MetaGraphs
using EzXML

function rbh_ds(paralog_df::DataFrame)
    @assert typeof(paralog_df[1,1]) <: AbstractString
    @assert typeof(paralog_df[1,2]) <: AbstractString
    @assert typeof(paralog_df[1,3]) <: AbstractFloat

    unique_ids = unique(vcat(paralog_df[:,1], paralog_df[:,2]))
    ids_to_ind_dict = Dict(unique_ids[i] => i for i in eachindex(unique_ids))
    ind_to_ids_dict = Dict(i => unique_ids[i] for i in eachindex(unique_ids))

    # Create a matrix of zeros
    rbh_matrix = zeros(Float64, length(unique_ids), length(unique_ids))
    orig_mat = zeros(Float64, length(unique_ids), length(unique_ids))

    # Fill in the matrix, treating gene 'i' as the query and gene 'j' as the subject
    for row in eachrow(paralog_df)
        

        i = ids_to_ind_dict[row[1]]
        j = ids_to_ind_dict[row[2]]

        orig_mat[i,j] = row[3]
        orig_mat[j,i] = row[3]
        rbh_matrix[i,j] = row[3]
        rbh_matrix[j,i] = row[3]
    end

    rbh_gene, rbh_paralog = String[], String[]
    matched_inds = Int[]
    score_i_origs, score_j_origs, max_scores, mean_scores = Float64[], Float64[], Float64[], Float64[]
    for i in 1:size(rbh_matrix)[1]

        score_i, max_i = findmax(rbh_matrix[i,:])[1:2]
        score_j, max_j = findmax(rbh_matrix[:,max_i])[1:2]
        score_i_orig = orig_mat[i,max_i]
        score_j_orig = orig_mat[max_i,i]
        max_score = max(score_i_orig, score_j_orig)
        mean_score = (score_i_orig + score_j_orig) / 2
        
        if i == max_j && max_i ∉ matched_inds && max_j ∉ matched_inds
            
            id_i = ind_to_ids_dict[max_i]
            id_j = ind_to_ids_dict[max_j]
            
            push!(rbh_gene, id_i)
            push!(rbh_paralog, id_j)
            push!(matched_inds, max_i)
            push!(matched_inds, max_j)
            push!(score_i_origs, score_i_orig)
            push!(score_j_origs, score_j_orig)
            push!(max_scores, max_score)
            push!(mean_scores, mean_score)
        end
    end

    return DataFrame(
        "GeneID" => rbh_gene, 
        "ParalogID" => rbh_paralog, 
        "ds" => score_i_origs
        )
end

"""
    rbh(paralog_df; scoring="max")

Identify reciprocal best hits (RBH) between paralogs based on similarity scores.

# Arguments
- `paralog_df::DataFrame`: DataFrame with at least 4 columns:
  1. GeneID (String)
  2. ParalogID (String)
  3. Percent identity from gene to paralog (Float), or dS value if using `scoring="ds"`
  4. Percent identity from paralog to gene (Float), or ignored if using `scoring="ds"`
- `scoring::String="max"`: Scoring method for determining best hits
  - `"max"` or `"maximum"`: Use maximum of bidirectional scores
  - `"mean"`, `"avg"`, or `"average"`: Use mean of bidirectional scores
  - `"double_max"`: Use original bidirectional scores
  - `"ds"`: Use dS values (column 3 should contain dS values)

# Returns
- `DataFrame`: Contains reciprocal best hit pairs with columns:
  - `GeneID`: First gene in RBH pair
  - `ParalogID`: Second gene in RBH pair
  - `perc_1`: Original percent identity (gene → paralog)
  - `perc_2`: Original percent identity (paralog → gene)
  - `max_perc`: Maximum of the two scores
  - `mean_perc`: Mean of the two scores

# Examples
```julia
df = DataFrame(
    GeneID = ["A", "B", "C"],
    ParalogID = ["B", "A", "D"],
    Perc1 = [95.0, 94.0, 80.0],
    Perc2 = [94.0, 95.0, 85.0]
)
rbh_pairs = rbh(df; scoring="max")
```
"""
function rbh(paralog_df::DataFrame; scoring::String="max")

    scoring = lowercase(scoring)

    if !(scoring in ["max", "maximum", "double_max", "mean", "average", "avg", "ds"])

        error("Invalid scoring method. Must be 'ds', 'max', 'maximum', 'double_max', 'mean', 'avg', or 'average'.")
    end

    scoring = scoring in ["max", "maximum"] ? "max" : scoring in ["mean", "avg", "average"] ? "mean" : scoring

    if scoring == "ds"
        return rbh_ds(paralog_df)
    end

    @assert typeof(paralog_df[1,1]) <: AbstractString
    @assert typeof(paralog_df[1,2]) <: AbstractString
    @assert typeof(paralog_df[1,3]) <: AbstractFloat
    @assert typeof(paralog_df[1,4]) <: AbstractFloat

    unique_ids = unique(vcat(paralog_df[:,1], paralog_df[:,2]))
    ids_to_ind_dict = Dict(unique_ids[i] => i for i in eachindex(unique_ids))
    ind_to_ids_dict = Dict(i => unique_ids[i] for i in eachindex(unique_ids))

    # Create a matrix of zeros
    rbh_matrix = zeros(Float64, length(unique_ids), length(unique_ids))
    orig_mat = zeros(Float64, length(unique_ids), length(unique_ids))

    # Fill in the matrix, treating gene 'i' as the query and gene 'j' as the subject
    for row in eachrow(paralog_df)

        i = ids_to_ind_dict[row[1]]
        j = ids_to_ind_dict[row[2]]
        
        orig_mat[i,j] = row[3]
        orig_mat[j,i] = row[4]

        if scoring == "max"

            max_perc_temp = max(row[3], row[4])
            rbh_matrix[i,j] = max_perc_temp
            rbh_matrix[j,i] = max_perc_temp
        elseif scoring == "mean"

            mean_perc_temp = (row[3] + row[4]) / 2
            rbh_matrix[i,j] = mean_perc_temp
            rbh_matrix[j,i] = mean_perc_temp
        else

            rbh_matrix[i,j] = orig_mat[i,j]
            rbh_matrix[j,i] = orig_mat[j,i]
        end
    end

    rbh_gene, rbh_paralog = String[], String[]
    matched_inds = Int[]
    perc_i_origs, perc_j_origs, max_percs, mean_percs = Float64[], Float64[], Float64[], Float64[]
        
    for i in 1:size(rbh_matrix)[1]

        perc_i, max_i = findmax(rbh_matrix[i,:])[1:2]
        perc_j, max_j = findmax(rbh_matrix[:,max_i])[1:2]
        perc_i_orig = orig_mat[i,max_i]
        perc_j_orig = orig_mat[max_i,i]
        max_perc = max(perc_i_orig, perc_j_orig)
        mean_perc = (perc_i_orig + perc_j_orig) / 2
        
        if i == max_j && max_i ∉ matched_inds && max_j ∉ matched_inds
            
            id_i = ind_to_ids_dict[max_i]
            id_j = ind_to_ids_dict[max_j]
            
            push!(rbh_gene, id_i)
            push!(rbh_paralog, id_j)
            push!(matched_inds, max_i)
            push!(matched_inds, max_j)
            push!(perc_i_origs, perc_i_orig)
            push!(perc_j_origs, perc_j_orig)
            push!(max_percs, max_perc)
            push!(mean_percs, mean_perc)
        end
    end

    return DataFrame("GeneID" => rbh_gene, "ParalogID" => rbh_paralog, "perc_1" => perc_i_origs, "perc_2" => perc_j_origs, "max_perc" => max_percs, "mean_perc" => mean_percs)
end

"""
    findfamilies(paralog_df; apply_cutoff=false, cutoff_val=nothing, cutoff_variable=nothing, cutoff_comparison, weighted_graph=false, add_names=false)

Identify gene families from paralog relationships using graph-based clustering.

This function creates a graph from paralog pairs and identifies connected components
as gene families. Optionally applies cutoffs and creates weighted graphs.

# Arguments
- `paralog_df::DataFrame`: DataFrame with at least 3 columns:
  1. GeneID (String)
  2. ParalogID (String)
  3. Edge value (Float) - e.g., dS, dN, average % identity - used for weighting if `weighted_graph=true`
  4. (Optional) Expression value for GeneID (Float) - used if expression data is available
  5. (Optional) Expression value for ParalogID (Float) - used if expression data is available
- `apply_cutoff::Bool=false`: Whether to apply filtering cutoff
- `cutoff_val::Union{Float64, Nothing}=nothing`: Cutoff threshold value (required if `apply_cutoff=true`)
- `cutoff_variable::Union{String, Symbol, Nothing}=nothing`: Column name for cutoff comparison (required if `apply_cutoff=true`)
- `cutoff_comparison::Function`: Comparison function (e.g., `<`, `>`, `<=`) (required if `apply_cutoff=true`)
- `weighted_graph::Bool=false`: Create weighted MetaGraph instead of unweighted SimpleGraph
- `add_names::Bool=false`: Add gene names as vertex properties

# Returns
- `Union{SimpleGraph, MetaGraph}`: Graph representing paralog relationships
  - Unweighted: `SimpleGraph` with vertices as gene indices
  - Weighted: `MetaGraph` with edge weights and optional metadata

# Notes
- Assumes pairwise symmetry: edge value from a→b equals b→a
- Expression values (columns 4-5) are optional
- Connected components in the graph represent gene families

# Examples
```julia
df = DataFrame(
    GeneID = ["A", "A", "B"],
    ParalogID = ["B", "C", "C"],
    dS = [0.5, 0.6, 0.4]
)
graph = findfamilies(df; apply_cutoff=true, cutoff_val=0.55, 
                     cutoff_variable=:dS, cutoff_comparison=<)
```
"""
function findfamilies(paralog_df::DataFrame; apply_cutoff::Bool=false,
                                            cutoff_val::Union{Float64, Nothing}, 
                                            cutoff_variable::Union{String, Symbol, Nothing}=nothing, 
                                            cutoff_comparison::F, 
                                            weighted_graph::Bool=false, 
                                            add_names::Bool=false) where F <: Function

    if ncol(paralog_df) != 5 && ncol(paralog_df) != 3

        error("Invalid input DataFrame. Must have 3, or 5 columns.")
    elseif ncol(paralog_df) == 5 && !(typeof(paralog_df[1,4]) <: AbstractFloat && typeof(paralog_df[1,5]) <: AbstractFloat)

        error("Invalid input DataFrame. Columns 4 and 5 must be <: AbstractFloat.")
    elseif typeof(paralog_df[1,1]) <: AbstractString && typeof(paralog_df[1,2]) <: AbstractString && typeof(paralog_df[1,3]) <: AbstractFloat

        nothing
    else

        error("Invalid input DataFrame. Items in columns 1 and 2 must be <: AbstractString, and column 3 must be <: AbstractFloat.")
    end

    if apply_cutoff
        
        filtered_df = filter(row -> cutoff_comparison(row[cutoff_variable], cutoff_val), paralog_df)
    else
        filtered_df = paralog_df
    end

    if weighted_graph

        return create_weighted_graph(filtered_df, cutoff_variable, cutoff_val=isnothing(cutoff_val) ? Inf : cutoff_val, add_names=add_names, add_expr=ncol(paralog_df) == 5) # Change this so adding exprsesion is an explicit parameter, not assumed (also change so that exprssion values can be in other columns)
    else

        return create_unweighted_graph(filtered_df, add_names=add_names)
    end
end

"""
    family_sizes(name_vecs)

Calculate family sizes for genes grouped by their gene families.

# Arguments
- `name_vecs::Vector{Vector{String}}`: Vector of vectors, where each inner vector 
  contains gene names belonging to the same family

# Returns
- `DataFrame`: Contains columns:
  - `GeneID`: Gene identifier
  - `FamilySize`: Number of genes in the family (including this gene)

# Examples
```julia
families = [["geneA", "geneB"], ["geneC", "geneD", "geneE"], ["geneF"]]
sizes_df = family_sizes(families)
# Returns: geneA and geneB have size 2, geneC/D/E have size 3, geneF has size 1
```
"""
function family_sizes(name_vecs::Vector{Vector{String}})

    family_sizes = length.(name_vecs)
    gene_names, family_size = String[], Int[]
    
    for (i, vec) in enumerate(name_vecs)

        for gene_name in vec

            push!(gene_names, gene_name)
            push!(family_size, family_sizes[i])
        end
    end

    return DataFrame("GeneID" => gene_names, "FamilySize" => family_size)
end

"""
    create_unweighted_graph(paralog_df; add_names=false)

Create an unweighted simple graph from paralog relationships.

This function constructs a graph where vertices represent genes and edges represent
paralog relationships. Each unique gene becomes a vertex, and paralog pairs are 
connected by edges.

# Arguments
- `paralog_df::DataFrame`: DataFrame with at least 2 columns:
  1. GeneID (String)
  2. ParalogID (String)
- `add_names::Bool=false`: If true, returns a MetaGraph with gene names stored as 
  vertex properties; if false, returns a SimpleGraph

# Returns
- `Union{SimpleGraph, MetaGraph}`: 
  - `SimpleGraph`: Basic graph with numbered vertices (if `add_names=false`)
  - `MetaGraph`: Graph with `:name` property storing gene IDs (if `add_names=true`)

# Examples
```julia
df = DataFrame(
    GeneID = ["A", "B", "C"],
    ParalogID = ["B", "C", "A"]
)
graph = create_unweighted_graph(df)
graph_named = create_unweighted_graph(df; add_names=true)
```
"""
function create_unweighted_graph(paralog_df::DataFrame; add_names::Bool=false)

    ids = String[]
    add_vert_1, add_vert_2 = true, true
    dup_graph = SimpleGraph()
    
    for i in 1:size(paralog_df)[1]

        id_1 = paralog_df[i,1]
        id_2 = paralog_df[i,2]

        id_1_vert = findfirst(ids .== id_1)
        id_2_vert = findfirst(ids .== id_2)

        if isnothing(id_1_vert)

            push!(ids, id_1)
            id_1_vert = length(ids)
            add_vert_1 = true
        else

            add_vert_1 = false
        end

        if isnothing(id_2_vert)

            push!(ids, id_2)
            id_2_vert = length(ids) 
            add_vert_2 = true
        else

            add_vert_2 = false
        end

        if add_vert_1
                
            add_vertex!(dup_graph)
        end

        if add_vert_2

            add_vertex!(dup_graph)
        end

        if has_edge(dup_graph, id_1_vert, id_2_vert)

            nothing
        else

            add_edge!(dup_graph, id_1_vert, id_2_vert)
        end
    end

    if add_names

        dup_graph = MetaGraph(dup_graph, 1)
        
        for i in eachindex(ids)

            set_prop!(dup_graph, i, :name, ids[i])
        end
    end

    return dup_graph
end

"""
    create_weighted_graph(paralog_df, weight_variable; cutoff_val=Inf, add_names=false, add_expr=false)

Create a weighted MetaGraph from paralog relationships with edge weights.

This function constructs a MetaGraph where edges are weighted by a specified variable
(e.g., dS values, percent identity). Optionally includes gene names and expression data
as vertex properties.

# Arguments
- `paralog_df::DataFrame`: DataFrame with at least 3 columns:
  1. GeneID (String)
  2. ParalogID (String)
  3. Weight variable (Float) - the column specified by `weight_variable`
  4. (Optional) Expression value for GeneID (Float) - required if `add_expr=true` and columns 4-5 contain expression data
  5. (Optional) Expression value for ParalogID (Float) - required if `add_expr=true` and columns 4-5 contain expression data
- `weight_variable::Union{String, Symbol}`: Column name in `paralog_df` to use for edge weights
- `cutoff_val::Float64=Inf`: Default weight value for graph initialization
- `add_names::Bool=false`: Add gene names as `:name` vertex property
- `add_expr::Bool=false`: Add expression values as `:expr` vertex property (requires columns 4-5 to contain expression data)

# Returns
- `MetaGraph`: Weighted graph with:
  - Edge property `:weight` containing -log10 of the weight variable
  - Vertex property `:name` if `add_names=true`
  - Vertex property `:expr` if `add_expr=true`

# Notes
- Edge weights are stored as -log10(weight_variable) for proper scaling
- Expression values from columns 4 and 5 are averaged for each gene

# Examples
```julia
df = DataFrame(
    GeneID = ["A", "B"],
    ParalogID = ["B", "C"],
    dS = [0.5, 0.3],
    expr1 = [10.0, 15.0],
    expr2 = [12.0, 20.0]
)
graph = create_weighted_graph(df, :dS; add_names=true, add_expr=true)
```
"""
function create_weighted_graph(paralog_df::DataFrame, weight_variable::Union{String, Symbol}; cutoff_val::Float64=Inf, add_names::Bool=false, add_expr::Bool=false)

    ids = String[]
    add_vert_1, add_vert_2 = true, true
    dup_graph = MetaGraph(SimpleGraph(), cutoff_val)

    for i in 1:size(paralog_df)[1]

        id_1 = paralog_df[i,1]
        id_2 = paralog_df[i,2]

        id_1_vert = findfirst(ids .== id_1)
        id_2_vert = findfirst(ids .== id_2)

        if isnothing(id_1_vert)

            push!(ids, id_1)
            id_1_vert = length(ids)
            add_vert_1 = true
        else

            add_vert_1 = false
        end

        if isnothing(id_2_vert)

            push!(ids, id_2)
            id_2_vert = length(ids) 
            add_vert_2 = true
        else

            add_vert_2 = false
        end

        if add_vert_1

            add_vertex!(dup_graph)
        end

        if add_vert_2

            add_vertex!(dup_graph)
        end

        if has_edge(dup_graph, id_1_vert, id_2_vert)
            
            nothing
        else

            add_edge!(dup_graph, id_1_vert, id_2_vert)
            set_prop!(dup_graph, id_1_vert, id_2_vert, :weight, paralog_df[i, weight_variable])

            if add_expr
                    
                set_prop!(dup_graph, id_1_vert, :expr, paralog_df[i,4])
                set_prop!(dup_graph, id_2_vert, :expr, paralog_df[i,5])
            end
        end
    end

    if add_names

        for i in 1:length(ids)

            set_prop!(dup_graph, i, :name, ids[i])
        end
    end

    return dup_graph
end

"""
    savegraphcsv(graph, filename; force=false, add_names=false)

Save a MetaGraph to a CSV file in edge list format.

This function exports a weighted graph to a CSV file containing edge information
with source, target, and weight columns. Optionally includes gene names.

# Arguments
- `graph::MetaGraph`: Weighted graph with `:weight` property on edges
- `filename::String`: Output CSV filename (should include .csv extension)
- `force::Bool=false`: Overwrite existing file if true, error if false
- `add_names::Bool=false`: Include gene names in output (requires `:name` vertex property)

# Returns
- `Nothing`: Writes output to file

# Output Format
Without names:
```
Source,Target,Weight
1,2,0.5
2,3,0.3
```

With names:
```
Source,Target,Weight,Source_Name,Target_Name
1,2,0.5,geneA,geneB
2,3,0.3,geneB,geneC
```

# Notes
- Weights are converted from -log10 back to original scale (10^weight)
- Raises error if file exists and `force=false`
- Requires edges to have `:weight` property
- If `add_names=true`, vertices must have `:name` property

# Examples
```julia
savegraphcsv(graph, "paralog_network.csv")
savegraphcsv(graph, "paralog_network.csv"; force=true, add_names=true)
```
"""
function savegraphcsv(graph::MG, filename::String; force::Bool=false, add_names::Bool=false) where MG <: MetaGraph

    if isfile(filename) && !force

        error("File already exists. Use 'force=true' to overwrite.")
    end

    open(filename, "w") do io

        if add_names

            println(io, "Source,Target,Weight,Source_Name,Target_Name")
        else
           
            println(io, "Source,Target,Weight")
        end

        for (edge,weight_dict) in graph.eprops

            temp_str_vec = split(string(edge), " => ")
            v1 = match(r"[0-9]+", temp_str_vec[1]).match
            v2 = match(r"[0-9]+", temp_str_vec[2]).match

            if add_names
                
                v1_name = props(graph, parse(Int, v1))[:name]
                v2_name = props(graph, parse(Int, v2))[:name]

                println(io, "$v1_name,$v2_name,$(10 ^ weight_dict[:weight]),$v1_name,$v2_name")
            else
                
                println(io, "$v1,$v2,$(10 ^ weight_dict[:weight])")
            end
        end
    end
end

export rbh, 
    findfamilies, 
    family_sizes, 
    create_unweighted_graph, 
    create_weighted_graph, 
    savegraphcsv

end # module ParalogUtils