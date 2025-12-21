
module MiscUtils
using StatsBase
using DataFrames
using YeoJohnsonTrans
"""
    add_expression_to_paralogs(paralog_df, expression_df)

Add expression values to a paralog DataFrame by matching gene IDs.

# Arguments
- `paralog_df::DataFrame`: DataFrame with at least 2 columns in the following order:
  1. GeneID (String or AbstractString)
  2. ParalogID (String or AbstractString)
- `expression_df::DataFrame`: DataFrame with at least 2 columns in the following order:
  1. GeneID (String or AbstractString)
  2. Expression value (Float or AbstractFloat)

# Returns
- `DataFrame`: The original paralog DataFrame with two additional columns appended:
  - `GeneExpr`: Expression value for the gene (missing if not found)
  - `ParalogExpr`: Expression value for the paralog (missing if not found)
"""
function add_expression_to_paralogs(paralog_df::DataFrame, expression_df::DataFrame)
    gene_expr_vals = Vector{Union{Float64, Missing}}(undef, size(paralog_df)[1])
    paralog_expr_vals = Vector{Union{Float64, Missing}}(undef, size(paralog_df)[1])
    for i in 1:size(paralog_df)[1]
        gene_id = paralog_df[i, 1]
        paralog_id = paralog_df[i, 2]
        gene_expr_ind = findfirst(expression_df[:, 1] .== gene_id)
        paralog_expr_ind = findfirst(expression_df[:, 1] .== paralog_id)
        gene_expr_val = isnothing(gene_expr_ind) ? missing : expression_df[gene_expr_ind, 2]
        paralog_expr_val = isnothing(paralog_expr_ind) ? missing : expression_df[paralog_expr_ind, 2]
        gene_expr_vals[i] = gene_expr_val
        paralog_expr_vals[i] = paralog_expr_val
    end
    return hcat(paralog_df, DataFrame("GeneExpr" => gene_expr_vals, "ParalogExpr" => paralog_expr_vals))
end

# helper function
function normalize_yj(u::Vector{Float64})
    v = YeoJohnsonTrans.transform(u)
    return zscore(v)
end
end # module