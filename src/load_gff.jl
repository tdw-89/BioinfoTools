module LoadGFF
using GFF3
using DataFrames
using GenomicFeatures
using InteractiveUtils
using Statistics
using CSV
using CodecZlib
using ..GenomeTypes
import Base.replace!
import BioSequences: BioSequence
const default_upstream = 2000
const default_downstream = 2000
# Misplaced function
"""
    averagereplicates!(expr_data)
Average expression data across paired replicates in a DataFrame.

# Arguments
- `expr_data::DataFrame`: Expression data with replicates in adjacent columns. Must have at least 3 columns. The first two columns are assumed to be gene ID and length, and columns 3 onwards contain replicate expression values in adjacent pairs.

# Returns
- `Nothing` (modifies `expr_data` in place)

# Notes
This function assumes that there are two replicates for each sample, that each pair 
are in adjacent columns starting from column 3, and that the first two columns are gene ID and length respectively.
"""
function averagereplicates!(expr_data::DataFrame)
    new_names = String[]
    for i in 3:2:ncol(expr_data)
        expr_data[:,i] = mean(eachcol(expr_data[:,i:(i + 1)]))
        temp_name = names(expr_data)[i] * "_" * names(expr_data)[i + 1]
        push!(new_names, temp_name)
    end
    original_names = names(expr_data)[1:2]
    select!(expr_data, Not(4:2:ncol(expr_data)))
    rename!(expr_data, reduce(vcat, [original_names, new_names]))
end
# Utility functions
# Parse the attributes of a record, returning a Dict(attribute name => value)
"""
    parseattributes(attrs)
Parse GFF3 attribute vector into attribute names and a dictionary.

# Arguments
- `attrs::Vector{Pair{String, Vector{String}}}`: Vector of attribute key-value pairs from a GFF3 record

# Returns
- `Dict{String, Vector{String}}`: A dictionary mapping attribute names to their values
"""
function parseattributes(attrs::Vector{Pair{String, Vector{String}}})
    attr_names = [vec[1] for vec in attrs]
    attr_values = [vec[2] for vec in attrs]
    # Create a dictionary mapping attribute names to their values
    dict = Dict(attr_names .=> attr_values)
    return dict
end
"""
    get_id(id_string)
"""
function get_id(id_string::String)
    id_split = split(id_string, r"[:\-]")
    return length(id_split) > 1 ? join(id_split[2:end], "-") : first(id_split)
end
"""
    get_id(id_field, dict)
Retrieve the gene ID from the attributes dictionary using the specified ID field.
# Arguments
- `id_field::String`: The ID field to retrieve (e.g., "ID", "gene_id")
- `dict::Dict{String, Vector{String}}`: The attributes dictionary
# Returns
- `Union{String, Nothing}`: The gene ID if found, or `nothing` if not found
"""
function get_id(id_field::String, dict::Dict{String, Vector{String}})
    if id_field in keys(dict)
        id_attrs = dict[id_field]
        return get_id(first(id_attrs))
    else
        if id_field == "ID"
            @warn "ID not found in attributes for feature, check GFF file formatting. Skipping..."
        else
            @warn "'$id_field' not found in attributes for feature, check GFF file formatting. Skipping..."
        end
        return nothing
    end
end
"""
    create_region(strand, gene_start, gene_end, upstream, downstream, chrom_length[, scaffold, contig])
Create a `Region` object around a gene with specified upstream/downstream extensions.
# Arguments
- `strand::Char`: Gene strand ('+' or '-')
- `gene_start::Int`: Gene start position
- `gene_end::Int`: Gene end position
- `upstream::Int`: Number of bases upstream to include
- `downstream::Int`: Number of bases downstream to include
- `chrom_length::Int`: Length of the chromosome
- `scaffold::Union{Scaffold,Missing}=missing`: Associated scaffold
- `contig::Union{Contig,Missing}=missing`: Associated contig
# Returns
- `Region`: A Region object spanning the specified interval
"""
function create_region(strand::Char,
                       gene_start::Int, 
                       gene_end::Int, 
                       upstream::Int, 
                       downstream::Int, 
                       chrom_length::Int,
                       scaffold::Union{Scaffold,Missing}=missing,
                       contig::Union{Contig,Missing}=missing)
    return Region(scaffold,
                    contig,
                    strand == '+' ? max(gene_start - upstream, 1) : max(gene_start - downstream, 1),
                    strand == '+' ? min(gene_end + downstream, chrom_length) : min(gene_end + upstream, chrom_length),
                    Dict{String, Vector{Annotation}}(),
                    Vector{Union{Vector{UInt8}, Vector{UInt16}, Vector{UInt32}, Vector{UInt64}}}(),
                    BitVector[],
                    String[])
end
"""
    getchromlengths(gff_file)
Extract chromosome names and lengths from a GFF file.
# Arguments
- `gff_file::String`: Path to a GFF file (can be gzipped)
# Returns
- `DataFrame`: A DataFrame with columns `:chrom` and `:length`
"""
function getchromlengths(gff_file::String)
    if endswith(gff_file, ".gz")
        reader = GFF3.Reader(GzipDecompressorStream(open(gff_file, "r")))
    else
        reader = open(GFF3.Reader, gff_file)
    end
    chroms = String[]
    lengths = Int[]
    for record in reader
        if GFF3.seqid(record) ∉ chroms
            push!(chroms, GFF3.seqid(record))
            push!(lengths, GFF3.seqend(record))
        end
    end
    close(reader)
    return DataFrame(:chrom=>chroms, :length=>lengths)
end

"""
    create_chromlengths_file(gff_file, output_file)
Create a chromosome lengths file from a GFF file.
# Arguments
- `gff_file::String`: Path to the input GFF file
- `output_file::String`: Path to the output chromosome lengths file
# Returns
- `Nothing`
"""
function create_chromlengths_file(gff_file::String, output_file::String)
    chrom_lengths = getchromlengths(gff_file)
    CSV.write(output_file, chrom_lengths, writeheader=false)
end
"""
    loadgenome(gff_files, chrom_lengths_file=nothing; feature_type="gene", alt_id_field=nothing)
Load one or more GFF files into a `RefGenome` object.
# Arguments
- `gff_files::Union{Vector{String}, String}`: Path to a GFF file, directory containing GFF files, or vector of GFF file paths
- `chrom_lengths_file::Union{String, Nothing}=nothing`: Optional path to chromosome lengths file
- `feature_type::String="gene"`: Type of feature to load (e.g., "gene", "all")
- `alt_id_field::Union{String, Nothing}=nothing`: Alternative ID field to use instead of "ID"
# Returns
- `RefGenome`: A reference genome object populated with features from the GFF file(s)
"""
function loadgenome(gff_files::Union{Vector{String}, String}, chrom_lengths_file::Union{String, Nothing}=nothing; feature_type::String="gene", alt_id_field::Union{String, Nothing}=nothing)
    if typeof(gff_files) == String
        if isdir(gff_files)
            gff_files = readdir(gff_files, join=true)
            gff_files = gff_files[endswith.(gff_files, ".gff") .|| endswith.(gff_files, ".gff3")]
        else
            gff_files = [gff_files]
        end
    end
    ref_genome = nothing
    for i in eachindex(gff_files)
        if i == 1
            ref_genome = loadgff(gff_files[i], chrom_lengths_file; feature_type=feature_type, alt_id_field=alt_id_field)
        else
            loadgff(gff_files[i], chrom_lengths_file; feature_type=feature_type, genome=ref_genome, alt_id_field=alt_id_field)
        end
    end
    return ref_genome
end
"""
    loadgff(gff_file, chrom_lengths_file=nothing; feature_type="gene", genome=nothing, alt_id_field=nothing)
Parse a single GFF file and populate gene/RNA/exon information into a reference genome.
# Arguments
- `gff_file::String`: Path to a GFF file (can be gzipped)
- `chrom_lengths_file::Union{String, Nothing}=nothing`: Optional path to chromosome lengths file
- `feature_type::String="gene"`: Type of feature to load ("gene", "all", etc.)
- `genome::Union{RefGenome, Nothing}=nothing`: Existing RefGenome to add features to
- `alt_id_field::Union{String, Nothing}=nothing`: Alternative ID field name
# Returns
- `Union{RefGenome, Nothing}`: Returns RefGenome if `genome` is nothing, otherwise nothing
"""
function loadgff(gff_file::String, chrom_lengths_file::Union{String, Nothing}=nothing; 
                                    feature_type::String="gene", 
                                    genome::Union{RefGenome, Nothing}=nothing,
                                    alt_id_field::Union{String, Nothing}=nothing)
    # Load the chromosome lengths files, and make sure the chromosome names are strings:
    using_chrom_file = false
    if isnothing(chrom_lengths_file)
        chrom_lengths = getchromlengths(gff_file)
    else
        using_chrom_file = true
        chrom_lengths = CSV.read(chrom_lengths_file, DataFrame, header=false)
        chrom_lengths[!, 1] = string.(chrom_lengths[!, 1])
    end
    if endswith(gff_file, ".gz")
        reader = GFF3.Reader(GzipDecompressorStream(open(gff_file, "r")))
    else
        reader = open(GFF3.Reader, gff_file)
    end
    ref_genome = nothing
    ret_ref = false
    if isnothing(genome)
        ref_genome = RefGenome()
        ret_ref = true
    else
        ref_genome = genome
    end
    # Create place holders for the prev_gene and prev_rna variables
    prev_gene::Union{Nothing, Gene} = nothing
    prev_rna::Union{Nothing, RNA} = nothing
    if lowercase(feature_type) == "all"
        for record in reader
            if GFF3.featuretype(record) == "gene"
                parsegene!(record, ref_genome, chrom_lengths, using_chrom_file; alt_id_field=alt_id_field)
                prev_gene = last(ref_genome.genes[2])
            elseif contains(GFF3.featuretype(record), "RNA")
                parserna!(record, ref_genome, prev_gene)
                prev_rna = last(ref_genome.rnas)
            elseif GFF3.featuretype(record) == "exon"
                parseexon!(record, ref_genome, prev_rna, prev_gene)
            else
                # UNDER CONSTRUCTION
            end
        end
        if ret_ref
            return ref_genome
        else
            return nothing
        end
    elseif lowercase(feature_type) == "gene"
        for record in reader
            if GFF3.featuretype(record) == "gene"
                parsegene!(record, ref_genome, chrom_lengths, using_chrom_file; alt_id_field=alt_id_field)
            end
        end
        if ret_ref
            return ref_genome
        else
            return nothing
        end
    else
        # UNDER CONSTRUCTION
    end
end
"""
    parsegene!(record, refs, chrom_lengths, using_chrom_file; alt_id_field=nothing)
Parse a gene record from a GFF file and add it to the reference genome.

# Arguments
- `record::GFF3.Record`: GFF3 record to parse
- `refs::RefGenome`: Reference genome to add the gene to
- `chrom_lengths::DataFrame`: DataFrame containing chromosome lengths (must have chromosome names in column 1 and lengths in column 2)
- `using_chrom_file::Bool`: Whether a chromosome lengths file was used
- `alt_id_field::Union{String, Nothing}=nothing`: Alternative ID field name (if `nothing`, uses "ID")

# Returns
- `Nothing`
"""
function parsegene!(record::GFF3.Record, refs::RefGenome, chrom_lengths::DataFrame, using_chrom_file::Bool; alt_id_field::Union{String, Nothing}=nothing)
    scaffold = nothing
    contig = missing
    source = ""
    gene_start = -1
    gene_end = -1
    score = -1
    strand = '.'
    phase = 0
    id = ""
    name = ""
    region = missing
    chrom_length = 0
    seqid = GFF3.seqid(record)
    if seqid ∉ chrom_lengths[!, 1]
        error("chromosome ID '$seqid' not found in chromosome lengths")
    end
    chrom_ind = findfirst(chrom_lengths[:,1] .== seqid)
    chrom_length = chrom_lengths[chrom_ind,2]
    if haskey(refs.scaffolds, seqid)
        scaffold = refs.scaffolds[seqid]
    else 
        refs.scaffolds[seqid] = Scaffold(seqid,
                                        missing, 
                                        Gene[],
                                        Repeat[], 
                                        missing,
                                        missing,
                                        missing)
        scaffold = refs.scaffolds[seqid]
    end
    try
        source = GFF3.source(record)
    catch
        nothing
    end
    gene_start = GFF3.seqstart(record) 
    gene_end = GFF3.seqend(record)
    try
        score = GFF3.score(record)
    catch
        nothing
    end
    strand = GFF3.strand(record)
    if strand == STRAND_POS 
        strand = '+'
    elseif strand == STRAND_NEG
        strand = '-'
    else
        error("unrecognized/missing strand info: $strand for gene with attributes: $(GFF3.attributes(record)))")
    end
    try
        phase = GFF3.phase(record)
    catch
        nothing
    end
    attrs = GFF3.attributes(record)
    dict = parseattributes(attrs)
    id_field = isnothing(alt_id_field) ? "ID" : alt_id_field
    id = get_id(id_field, dict)
    if isnothing(id)
        return nothing
    end
    if  "Name" in keys(dict)
        name = dict["Name"][1]
    end
    # Add default region, which stretches 'default_upstream' bases down
    # from 'gene_start', and 'default_downstream' bases up from 'gene_end'
    # (for (+)-sense genes).
    #
    if !(ismissing(gene_start)) # DEBUG CHANGE?
        region = create_region(strand, gene_start, gene_end, default_upstream, default_downstream, chrom_length, scaffold, contig)
    else
        region = missing
    end
    temp_gene = Gene(scaffold,
                        contig,
                        name,
                        id,
                        strand,
                        RegElement[],
                        missing,
                        missing,
                        Intron[],
                        Exon[],
                        missing,
                        missing,
                        RNA[],
                        Segment[],
                        Region[region],
                        Dict{String, Vector{Annotation}}(),
                        gene_start,
                        gene_end,
                        nothing,
                        Vector{Union{Vector{UInt8}, Vector{UInt16}, Vector{UInt32}, Vector{UInt64}}}(),
                        BitVector[],
                        String[])
    if !(ismissing(scaffold))
        push!(scaffold.genes, temp_gene)
    end
    push!(refs.regions, region)
    push!(refs.genes[1], id)
    push!(refs.genes[2], temp_gene)
end
"""
    parserna!(record, refs, prev_gene)
Parse an RNA record from a GFF file and add it to the reference genome.
# Arguments
- `record::GFF3.Record`: GFF3 record to parse
- `refs::RefGenome`: Reference genome to add the RNA to
- `prev_gene::Union{Gene, Nothing}`: Previously parsed gene (parent of this RNA)
# Returns
- `Nothing`
"""
function parserna!(record::GFF3.Record, refs::RefGenome, prev_gene::Union{Gene, Nothing})
    rna_type = GFF3.featuretype(record)
    exon_vec = Exon[]
    intron_vec = Intron[]
    rna_start = GFF3.seqstart(record)
    rna_end = GFF3.seqend(record)
    attrs = Dict(GFF3.attributes(record))
    rna_id = get_id("ID", attrs)
    parent_gene = haskey(attrs, "Parent") ? get_id(attrs["Parent"][1]) : nothing
    if isnothing(prev_gene) || isnothing(parent_gene) || (parent_gene != prev_gene.id && parent_gene != prev_gene.name)
        @warn "Parent gene for RNA: $rna_id of type: $rna_type not found"
        if !isnothing(parent_gene)
            println("Expected parent gene: $parent_gene")
            if !isnothing(prev_gene)
                println("Got: $(prev_gene.id)")
            end
        end
        push!(refs.rnas, RNA(rna_id, rna_type, missing, exon_vec, intron_vec, String[], rna_start, rna_end, Float64[]))
    else
        new_rna = RNA(rna_id, rna_type, prev_gene, exon_vec, intron_vec, String[], rna_start, rna_end, Float64[])
        push!(refs.rnas, new_rna)
        push!(prev_gene.rnas, new_rna)
    end
end
"""
    parseexon!(record, refs, prev_rna, prev_gene)
Parse an exon record from a GFF file and add it to the reference genome.
# Arguments
- `record::GFF3.Record`: GFF3 record to parse
- `refs::RefGenome`: Reference genome to add the exon to
- `prev_rna::Union{RNA, Nothing}`: Previously parsed RNA (parent of this exon)
- `prev_gene::Union{Gene, Nothing}`: Previously parsed gene (grandparent of this exon)
# Returns
- `Nothing`
"""
function parseexon!(record::GFF3.Record, refs::RefGenome, prev_rna::Union{RNA, Nothing}, prev_gene::Union{Gene, Nothing})
    parent_id = get_id(Dict(GFF3.attributes(record))["Parent"][1])
    parent_rna = isnothing(prev_rna) ? nothing : prev_rna.id
    parent_gene = isnothing(prev_gene) ? nothing : prev_gene.id
    if isnothing(parent_rna) || isnothing(parent_gene) || parent_id != parent_rna
        @warn "could not find parent gene/rna for exon with stated parent $parent_id, check GFF file formatting. Skipping..."
        return nothing
    elseif ismissing(prev_rna.gene) || prev_rna.gene.id != parent_gene
        @warn "could not find parent gene/rna for exon with stated parent $parent_id, check GFF file formatting. Skipping..."
        return nothing
    end
    new_exon = Exon(prev_gene, prev_rna, GFF3.seqstart(record), GFF3.seqend(record))
    push!(prev_gene.exons, new_exon)
    push!(prev_rna.exons, new_exon)
    push!(refs.exons, new_exon)
end
"""
    addsequence!(ref_genome, repeat_elem, sequence)
Add a sequence to a specific repeat element in the reference genome.
# Arguments
- `ref_genome`: Reference genome containing the repeat element
- `repeat_elem::Repeat`: Repeat element to add sequence to
- `sequence::BioSequence`: Biological sequence to add
# Returns
- `Nothing`
"""
function addsequence!(ref_genome, repeat_elem::Repeat, sequence::B) where B<:BioSequence
    new_repeat = Repeat(repeat_elem.scaffold, 
                        repeat_elem.contig, 
                        repeat_elem.regions, 
                        repeat_elem.family, 
                        repeat_elem.type, 
                        repeat_elem.repeat_start,
                        repeat_elem.repeat_end,
                        sequence,
                        repeat_elem.signals,
                        repeat_elem.binsignals,
                        repeat_elem.samples)
    replace!(ref_genome.repeats, repeat_elem => new_repeat)
end
# Base show overloading
"""
    Base.show(io, x::RefGenome)
Display a summary of a `RefGenome` object showing present and missing fields.
"""
Base.show(io::IO, x::RefGenome) = begin
    type = typeof(x)
    fields = fieldnames(type)
    fields_missing = [isempty(getfield(x, i)) for i in fields]
    fields_present = fields[.!(fields_missing)]
    fields_missing = fields[fields_missing]
    if any(:id .== fields_present)
        id = x.id
        print(io, "$type: '$id', with fields: ")
        for field in fields_present
            print(io, "'$field' ")
        end
    else
        print(io, "$type with fields: ")
        for field in fields_present
            print(io, "'$field' ")
        end
    end
    print(io, "\n          empty fields: ")
    for field in fields_missing
        print(io, "'$field' ")
    end
end
"""
    Base.show(io, x::Tuple{Vector{String}, Vector{Gene}})
Display a summary of a gene-name/gene-data vector pair.
"""
Base.show(io::IO, x::Tuple{Vector{String}, Vector{Gene}}) = begin
    len = length(x[1])
    print(io, "Gene-name/Gene-data vector pair with $len genes")
end
export loadgenome,
       getchromlengths,
       loadgff,
       create_region,
       parsegene!,
       parserna!,
       parseexon!,
       parseattributes,
       averagereplicates!,
       addpromoters!,
       addsequences!,
       addsequence!
end # module LoadGFF