if !@isdefined Feature
    include("genome_types.jl")

end


"""
    This function assumes that the 'paralog' DataFrame has (at least) the following columns in the
    following order, with the following types:
        * 1: GeneID <: AbstractString
        * 2: ParalogID <: AbstractString
    This function also assumes that the 'expression' DataFrame has (at least) the following columns in the
    following order, with the following types:
        * 1: GeneID <: AbstractString
        * 2: Expression value <: AbstractFloat
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
