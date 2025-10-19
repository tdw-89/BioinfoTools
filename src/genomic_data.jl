module GenomicData

using StatsBase
using DataFrames
using XAM
using CSV
using CodecZlib: GzipDecompressorStream
using BED

export ChromData,
       SampleData,
       Experiment,
       PeakWarnings,
       getallrecords,
       getallcountvectors,
       average_bam_replicate_groups,
       average_peak_replicate_groups,
       average_bam_replicates,
       average_peak_replicates,
       binpeaks,
       binpeakshomer,
       addtogenes!,
       addexpression!,
       addtorepeats!,
       has_expr,
       sample_to_bed_df,
       contiguous_values,
       addpeak!

"""
    ChromData(name, signal)

Container pairing a chromosome identifier with a numeric signal vector.
"""
mutable struct ChromData{T<:Real,V<:AbstractVector{T}}
    name::String
    signal::V

    function ChromData{T,V}(name::String, signal::V) where {T<:Real,V<:AbstractVector{T}}
        return new{T,V}(name, signal)
    end
end

ChromData(name::AbstractString, signal::V) where {T<:Real,V<:AbstractVector{T}} =
    ChromData{T,V}(String(name), signal)

"""
    SampleData(name, chroms[, n_reads])

Bundle of chromosome signals for a single experimental sample.
"""
mutable struct SampleData{T<:Real,V<:AbstractVector{T},R<:Real}
    name::String
    chroms::Vector{ChromData{T,V}}
    n_reads::R

    function SampleData{T,V,R}(name::String, chroms::Vector{ChromData{T,V}}, n_reads::R) where {T<:Real,V<:AbstractVector{T},R<:Real}
        return new{T,V,R}(name, chroms, n_reads)
    end
end

SampleData(name::AbstractString, chroms::Vector{ChromData{T,V}}, n_reads::R) where {T<:Real,V<:AbstractVector{T},R<:Real} =
    SampleData{T,V,R}(String(name), chroms, n_reads)
SampleData(name::AbstractString, chroms::Vector{ChromData{T,V}}) where {T<:Real,V<:AbstractVector{T}} =
    SampleData{T,V,Int}(String(name), chroms, -1)

"""
    Experiment(name, samples)

Collection of `SampleData` objects forming a cohesive experiment.
"""
mutable struct Experiment{T<:Real,V<:AbstractVector{T},R<:Real}
    name::String
    samples::Vector{SampleData{T,V,R}}

    function Experiment{T,V,R}(name::String, samples::Vector{SampleData{T,V,R}}) where {T<:Real,V<:AbstractVector{T},R<:Real}
        return new{T,V,R}(name, samples)
    end
end

Experiment(name::AbstractString, samples::Vector{SampleData{T,V,R}}) where {T<:Real,V<:AbstractVector{T},R<:Real} =
    Experiment{T,V,R}(String(name), samples)

"""
    PeakWarnings()

Counters and toggles for suppressing repeated warnings while processing peak files.
"""
mutable struct PeakWarnings
    counts::Vector{Int}
    switches::Vector{Bool}
end

PeakWarnings() = PeakWarnings([0, 0], [true, true])

"""
    _foreach_bed_record(path, fn)

Internal helper that iterates over every BED record contained in `path`, calling `fn` for each record.
Handles gzipped inputs transparently.
"""
function _foreach_bed_record(path::AbstractString, fn::Function)
    lower_path = lowercase(path)

    if endswith(lower_path, ".gz") && !endswith(lower_path, ".bgz")
        open(path, "r") do fh
            stream = GzipDecompressorStream(fh)
            reader = BED.Reader(stream)
            try
                for record in reader
                    fn(record)
                end
            finally
                close(stream)
            end
        end
    else
        reader = BED.Reader(path; index=nothing)
        try
            for record in reader
                fn(record)
            end
        finally
            close(BED.BioGenerics.IO.stream(reader))
        end
    end

    return nothing
end

_foreach_bed_record(fn::Function, path::AbstractString) = _foreach_bed_record(path, fn)

"""
    getallrecords(bam_file_path; mapq_signal=false)

Load a BAM file and accumulate per-base signal. Returns `(chrom_data, total_read_count)`.
"""
function getallrecords(bam_file_path::String; mapq_signal::Bool=false)
    return open(BAM.Reader, bam_file_path) do bam_reader
        chroms = bam_reader.refseqnames
        chrom_lengths = bam_reader.refseqlens
        chrom_count = length(chroms)
        signal_eltype = mapq_signal ? Float64 : UInt16
        count_vec_list = Vector{ChromData{signal_eltype, Vector{signal_eltype}}}(undef, chrom_count)
        chrom_indices = Dict{String, Int}(chroms[i] => i for i in eachindex(chroms))

        for i in 1:chrom_count
            chr_len = chrom_lengths[i]
            buffer = mapq_signal ? zeros(Float64, chr_len) : zeros(UInt16, chr_len)
            count_vec_list[i] = ChromData(chroms[i], buffer)
        end

        temp_record = BAM.Record()
        total_read_count = 0

        while !eof(bam_reader)
            empty!(temp_record)
            read!(bam_reader, temp_record)
            total_read_count += 1

            if BAM.ismapped(temp_record)
                start_pos = BAM.position(temp_record)
                end_pos = BAM.rightposition(temp_record)
                chrom_name = BAM.refname(temp_record)
                ind = chrom_indices[chrom_name]
                signal = count_vec_list[ind].signal

                if mapq_signal
                    mapq = mean(BAM.mappingquality(temp_record))
                    @inbounds signal[start_pos:end_pos] .+= mapq
                else
                    @inbounds signal[start_pos:end_pos] .+= one(eltype(signal))
                end
            end
        end

        return count_vec_list, total_read_count
    end
end

"""
    getallcountvectors(bam_file_list; exp_name="seq_exp", mapq_signal=false)

Apply [`getallrecords`](@ref) to each BAM in `bam_file_list` and aggregate into an `Experiment`.
"""
function getallcountvectors(bam_file_list::Vector{String}; exp_name::String="seq_exp", mapq_signal::Bool=false)
    sample_names = basename.(bam_file_list)
    ret_chroms, total_reads = getallrecords(bam_file_list[1]; mapq_signal=mapq_signal)
    signal_eltype = eltype(ret_chroms[1].signal)
    signal_type = typeof(ret_chroms[1].signal)
    read_type = typeof(total_reads)
    samples = Vector{SampleData{signal_eltype, signal_type, read_type}}(undef, length(bam_file_list))
    samples[1] = SampleData(sample_names[1], ret_chroms, total_reads)

    for (idx, bam_path) in enumerate(bam_file_list[2:end], 2)
        chroms, n_reads = getallrecords(bam_path; mapq_signal=mapq_signal)
        samples[idx] = SampleData(sample_names[idx], chroms, n_reads)
    end

    return Experiment(exp_name, samples)
end

"""
    average_bam_replicate_groups(bam_data; replicate_list=nothing)

Average BAM replicate groups according to the provided replicate specification.
"""
function average_bam_replicate_groups(bam_data::Experiment{T,V,R}; replicate_list::Union{String, DataFrame, Nothing}=nothing) where {T<:Real,V<:AbstractVector{T},R<:Real}
    n_chroms = length(bam_data.samples[1].chroms)
    @assert all(length(sample.chroms) == n_chroms for sample in bam_data.samples)

    if isa(replicate_list, String)
        replicate_list = CSV.read(replicate_list, DataFrame; header=false)
    elseif isnothing(replicate_list)
        sample_names = getfield.(bam_data.samples, :name)
        replicate_list = DataFrame(["$i" => sample_name for (i, sample_name) in enumerate(sample_names)])
    end

    avg_samples = SampleData{Float64,Vector{Float64},Float64}[]
    for i in 1:nrow(replicate_list)
        replicate_names = replicate_list[i, :]
        group = [sample for sample in bam_data.samples if sample.name in replicate_names]
        push!(avg_samples, average_bam_replicates(group))
    end

    return Experiment("avg_" * bam_data.name, avg_samples)
end

"""
    average_peak_replicate_groups(peak_data; replicate_list=nothing)

Average binary peak replicate groups analogous to [`average_bam_replicate_groups`](@ref).
"""
function average_peak_replicate_groups(peak_data::Experiment{Bool,BitVector,R}; replicate_list::Union{String, DataFrame, Nothing}=nothing) where {R<:Real}
    if isa(replicate_list, String)
        replicate_list = CSV.read(replicate_list, DataFrame; header=false)
    elseif isnothing(replicate_list)
        sample_names = getfield.(peak_data.samples, :name)
        replicate_list = DataFrame(["$i" => sample_name for (i, sample_name) in enumerate(sample_names)])
    end

    avg_samples = SampleData{Float64,Vector{Float64},R}[]
    for i in 1:nrow(replicate_list)
        replicate_names = replicate_list[i, :]
        group = [sample for sample in peak_data.samples if sample.name in replicate_names]
        push!(avg_samples, average_peak_replicates(group))
    end

    return Experiment("avg_" * peak_data.name, avg_samples)
end

"""
    average_bam_replicates(replicate_vec)

Average a set of BAM replicates after normalizing by total read count.
"""
function average_bam_replicates(replicate_vec::Vector{SampleData{T,V,R}}) where {T<:Real,V<:AbstractVector{T},R<:Real}
    count_mean = mean(getfield.(replicate_vec, :n_reads))
    count_ratios = count_mean ./ getfield.(replicate_vec, :n_reads)
    chrom_count = length(replicate_vec[1].chroms)
    return_chroms = [ChromData(replicate_vec[1].chroms[i].name, zeros(Float64, length(replicate_vec[1].chroms[i].signal))) for i in 1:chrom_count]
    inv_rep_count = inv(length(replicate_vec))

    for chrom_idx in 1:chrom_count
        target = return_chroms[chrom_idx].signal
        for (rep_idx, replicate) in enumerate(replicate_vec)
            source = replicate.chroms[chrom_idx].signal
            weight = count_ratios[rep_idx]
            @inbounds @simd for pos in eachindex(target)
                target[pos] += source[pos] * weight
            end
        end
        target .*= inv_rep_count
    end

    return SampleData("avg_" * replicate_vec[1].name, return_chroms, count_mean)
end

"""
    average_peak_replicates(replicate_vec)

Average binary peak replicates by computing per-position mean occupancy.
"""
function average_peak_replicates(replicate_vec::Vector{SampleData{Bool,BitVector,R}}) where {R<:Real}
    n_chroms = length(replicate_vec[1].chroms)
    @assert all(length(replicate.chroms) == n_chroms for replicate in replicate_vec) "All samples must have the same number of chromosomes"
    @assert all(all(length(replicate_vec[1].chroms[j].signal) == length(replicate.chroms[j].signal) for j in 1:n_chroms) for replicate in replicate_vec) "All samples must have the same chromosome lengths"

    return_chroms = [ChromData(replicate_vec[1].chroms[i].name, zeros(Float64, length(replicate_vec[1].chroms[i].signal))) for i in 1:n_chroms]
    inv_rep_count = inv(length(replicate_vec))

    for chrom_idx in 1:n_chroms
        target = return_chroms[chrom_idx].signal
        for replicate in replicate_vec
            @inbounds target .+= replicate.chroms[chrom_idx].signal
        end
        target .*= inv_rep_count
    end

    return SampleData("avg_" * replicate_vec[1].name, return_chroms, one(R))
end

"""
    addpeak!(chrom_dict, chrom, peak_start, peak_end, warnings)

Flag the specified range within `chrom_dict` as occupied.
"""
function addpeak!(chrom_dict::Dict{String, BitVector}, chrom::AbstractString, peak_start::Int, peak_end::Int, warnings::PeakWarnings)
    if haskey(chrom_dict, chrom)
        signal = chrom_dict[chrom]
        _set_peak!(signal, chrom, peak_start, peak_end, warnings)
    elseif length(chrom) > 3 && haskey(chrom_dict, chrom[4:end])
        trimmed = chrom[4:end]
        if warnings.switches[2]
            @warn "Chromosome '$(chrom)' not found in chromosome lengths file, using '$(trimmed)' instead."
            warnings.counts[2] += 1
            if warnings.counts[2] == 10
                @warn "Too many chromosome name warnings were issued. Suppressing further warnings..."
                warnings.switches[2] = false
            end
        end
        signal = chrom_dict[trimmed]
        _set_peak!(signal, chrom, peak_start, peak_end, warnings)
    else
        error("Chromosome $(chrom) not found in chromosome lengths file")
    end
end

function _set_peak!(signal::BitVector, chrom::AbstractString, peak_start::Int, peak_end::Int, warnings::PeakWarnings)
    try
        @inbounds signal[peak_start:peak_end] .= true
    catch BoundsError
        if warnings.switches[1]
            @warn "Bounds error for peak on chrom $(chrom)"
            warnings.counts[1] += 1
            if warnings.counts[1] == 10
                @warn "Too many bounds error warnings were issued. Suppressing further warnings..."
                warnings.switches[1] = false
            end
        end
        len = length(signal)
        @inbounds signal[clamp(peak_start, 1, len):clamp(peak_end, 1, len)] .= true
    end
end

"""
    binpeaks(peak_files, chrom_lengths_file)

Convert one or more narrowPeak/BED files into an `Experiment` with binary occupancy signals.
"""
function binpeaks(peak_files::Union{String, Vector{String}}, chrom_lengths_file::Union{String, Nothing}=nothing)
    if isa(peak_files, String)
        peak_files = [peak_files]
    end
    chrom_lengths_file === nothing && throw(ArgumentError("`chrom_lengths_file` must be provided"))

    chrom_lengths_df = CSV.read(chrom_lengths_file, DataFrame; header=false)
    chrom_lengths_df[!, 1] = string.(chrom_lengths_df[!, 1])
    rename!(chrom_lengths_df, ["chrom", "length"])

    template = Dict(String(row.chrom) => falses(row.length) for row in eachrow(chrom_lengths_df))
    samples = Vector{SampleData{Bool,BitVector,Int}}(undef, length(peak_files))
    warnings = PeakWarnings()

    for (idx, peak_path) in enumerate(peak_files)
        for signal_data in values(template)
            fill!(signal_data, false)
        end

        _foreach_bed_record(peak_path) do record
            addpeak!(template, BED.chrom(record), BED.chromstart(record), BED.chromend(record), warnings)
        end

        chroms = [ChromData(chrom, copy(signal)) for (chrom, signal) in template]
        samples[idx] = SampleData(String(split(basename(peak_path), ".")[1]), chroms, -1)
    end

    return Experiment("peak_data", samples)
end

"""
    binpeaks(narrow_peak_dir, chrom_lengths_file)

Directory wrapper for [`binpeaks(::Union{String,Vector{String}}, ...)`](@ref).
"""
function binpeaks(narrow_peak_dir::String, chrom_lengths_file::Union{String, Nothing}=nothing)
    if !isdir(narrow_peak_dir)
        return binpeaks([narrow_peak_dir], chrom_lengths_file)
    end
    peak_files = readdir(narrow_peak_dir, join=true)
    filter!(path -> occursin(r"\.(?:narrowPeak|bed)(?:\.gz)?$", path), peak_files)
    return binpeaks(peak_files, chrom_lengths_file)
end

"""
    binpeakshomer(peak_files, chrom_lengths_file; gff_files=nothing)

Construct binary peak tracks from HOMER peak output.
"""
function binpeakshomer(peak_files::Vector{String}, chrom_lengths_file::Union{String, Nothing}=nothing; gff_files::Union{String, Nothing, Vector{String}}=nothing)
    chrom_lengths_df = DataFrame()
    if chrom_lengths_file === nothing
        isnothing(gff_files) && error("Must provide either a file containing chromosome lengths or a GFF file.")
        gff_iter = isa(gff_files, String) ? (gff_files,) : gff_files
        for gff_file in gff_iter
            chrom_lengths_df = vcat(chrom_lengths_df, getchromlengths(gff_file))
        end
        rename!(chrom_lengths_df, ["chrom", "length"])
    else
        chrom_lengths_df = CSV.read(chrom_lengths_file, DataFrame; header=false)
        chrom_lengths_df[!, 1] = string.(chrom_lengths_df[!, 1])
        rename!(chrom_lengths_df, ["chrom", "length"])
    end

    template = Dict(String(row.chrom) => falses(row.length) for row in eachrow(chrom_lengths_df))
    samples = Vector{SampleData{Bool,BitVector,Int}}(undef, length(peak_files))
    warnings = PeakWarnings()

    for (idx, peak_file) in enumerate(peak_files)
        for signal_data in values(template)
            fill!(signal_data, false)
        end

        temp_df = CSV.read(peak_file, DataFrame; header=7)
        filter!(row -> row."p-value vs Control" <= 0.01, temp_df)
        temp_df = temp_df[:, ["#PeakID", "chr", "start", "end"]]
        rename!(temp_df, ["peakName", "chrom", "start", "end"])

        for peak in eachrow(temp_df)
            addpeak!(template, String(peak.chrom), peak.start + 1, peak.end, warnings)
        end

        chroms = [ChromData(chrom, copy(signal)) for (chrom, signal) in template]
        samples[idx] = SampleData(split(basename(peak_file), ".")[1], chroms, -1)
    end

    return Experiment("peak_data", samples)
end

"""
    addtogenes!(genome, experiment; peak_data=true, regions=true)

Attach experiment signal vectors to genes stored within `genome`.
"""
function addtogenes!(genome, experiment::Experiment{T,V,R}; peak_data::Bool=true, regions::Bool=true) where {T<:Real,V<:AbstractVector{T},R<:Real}
    if peak_data
        experiment.samples[1].chroms[1].signal isa BitVector ||
            error("If 'peak_data' is set to true, then the signal must be a BitVector")
    end

    for gene in genome.genes[2]
        gene_start = gene.gene_start
        gene_end = gene.gene_end
        ismissing(gene.scaffold) && error("'Gene' $(gene.id) is missing chromosome information (a 'Scaffold' object)")
        chrom_name = gene.scaffold.name

        for sample_data in experiment.samples
            for chrom_data in sample_data.chroms
                chrom_data.name == chrom_name || continue

                if peak_data
                    signal_vec = chrom_data.signal[gene_start:gene_end]
                    push!(gene.binsignals, signal_vec)

                    if regions
                        isempty(gene.regions) && error("adding regions inside 'addtogenes!' not implemented yet")
                        region_start = gene.regions[1].region_start
                        region_end = gene.regions[1].region_end
                        region_signal_vec = chrom_data.signal[region_start:region_end]
                        push!(gene.regions[1].binsignals, region_signal_vec)
                    end
                else
                    signal_vec = _convert_signal_slice(chrom_data.signal[gene_start:gene_end])
                    push!(gene.signals, signal_vec)

                    if regions
                        isempty(gene.regions) && error("adding regions inside 'addtogenes!' not implemented yet")
                        region_start = gene.regions[1].region_start
                        region_end = gene.regions[1].region_end
                        region_signal_vec = _convert_signal_slice(chrom_data.signal[region_start:region_end])
                        push!(gene.regions[1].signals, region_signal_vec)
                    end
                end

                break
            end

            push!(gene.samples, sample_data.name)
            if regions
                push!(gene.regions[1].samples, sample_data.name)
            end
        end
    end
end

"""
    addexpression!(ref_genome, expr_data; total_expr=true)

Add expression measurements from `expr_data` into `ref_genome`.
"""
function addexpression!(ref_genome, expr_data::DataFrame; total_expr::Bool=true)
    for row in eachrow(expr_data)
        gene_id = row.GeneID
        gene_idx = findfirst(gene_id .== ref_genome.genes[1])
        isnothing(gene_idx) && continue
        new_rna = RNA("all_transcripts", "mRNA", ref_genome.genes[2][gene_idx], Exon[], Intron[], names(row[3:end]), missing, missing, Float64.(Vector(row[2:end])))
        pushfirst!(ref_genome.genes[2][gene_idx].rnas, new_rna)
    end
end

"""
    addtorepeats!(genome, experiment; peak_data=false, add_to_region=false)

Attach experimental signals to repeat annotations.
"""
function addtorepeats!(genome, experiment::Experiment{T,V,R}; peak_data::Bool=false, add_to_region::Bool=false) where {T<:Real,V<:AbstractVector{T},R<:Real}
    if peak_data
        experiment.samples[1].chroms[1].signal isa BitVector ||
            error("If 'peak_data' is set to true, then the signal must be a BitVector")
    elseif !(experiment.samples[1].chroms[1].signal isa Vector{UInt16})
        error("If 'peak_data' is set to false (default), then the signal must be a UInt16 vector")
    end

    for repeat_elem in genome.repeats
        ismissing(repeat_elem.scaffold) && error("Repeat element is missing associated scaffold information.")
        chrom_name = repeat_elem.scaffold.name
        repeat_start = repeat_elem.repeat_start
        repeat_end = repeat_elem.repeat_end

        for sample_data in experiment.samples
            for chrom_data in sample_data.chroms
                chrom_data.name == chrom_name || continue

                if peak_data
                    if add_to_region
                        isempty(repeat_elem.regions) && error("No method for adding default region to repeat element")
                        region_start = repeat_elem.regions[1].region_start
                        region_end = repeat_elem.regions[1].region_end
                        region_signal_vec = chrom_data.signal[region_start:region_end]
                        push!(repeat_elem.regions[1].binsignals, region_signal_vec)
                    end
                    repeat_signal_vec = chrom_data.signal[repeat_start:repeat_end]
                    push!(repeat_elem.binsignals, repeat_signal_vec)
                else
                    isempty(repeat_elem.regions) && error("No method for adding default region to repeat element")
                    if add_to_region
                        region_start = repeat_elem.regions[1].region_start
                        region_end = repeat_elem.regions[1].region_end
                        region_signal_vec = _convert_signal_slice(chrom_data.signal[region_start:region_end])
                        push!(repeat_elem.regions[1].signals, region_signal_vec)
                    end
                    repeat_signal_vec = _convert_signal_slice(chrom_data.signal[repeat_start:repeat_end])
                    push!(repeat_elem.signals, repeat_signal_vec)
                end

                break
            end

            push!(repeat_elem.samples, sample_data.name)
            if add_to_region
                push!(repeat_elem.regions[1].samples, sample_data.name)
            end
        end
    end
end

"""
    has_expr(gene)

Return true if `gene` carries expression data.
"""
function has_expr(gene)
    return !isempty(gene.rnas) && !isempty(gene.rnas[1].expression)
end

"""
    sample_to_bed_df(sample)

Convert a `SampleData` into a BED-like `DataFrame`.
"""
function sample_to_bed_df(sample::SampleData{T,V,R}) where {T<:Real,V<:AbstractVector{T},R<:Real}
    chroms = String[]
    starts = Int[]
    ends = Int[]
    scores = Int[]

    for chrom in sample.chroms
        for (range, value) in contiguous_values(chrom.signal)
            push!(chroms, chrom.name)
            push!(starts, first(range))
            push!(ends, last(range))
            push!(scores, round(Int, value * 1000))
        end
    end

    df = DataFrame("Chrom" => chroms, "Start" => starts, "End" => ends, "Score" => scores)
    insertcols!(df, 4, :Name => string.(1:nrow(df)))
    return df
end

"""
    contiguous_values(signal, val[, val_complement])

Return contiguous index ranges where `signal` equals (`val_complement=false`) or differs from `val`.
"""
function contiguous_values(signal::AbstractVector{T}, val::T, val_complement::Bool=false) where {T<:Real}
    ranges = UnitRange{Int}[]
    start_idx = 0
    in_range = false

    for i in eachindex(signal)
        matches = val_complement ? signal[i] != val : signal[i] == val
        if matches
            if !in_range
                start_idx = i
                in_range = true
            end
        elseif in_range
            push!(ranges, start_idx:(i - 1))
            in_range = false
        end
    end

    if in_range
        push!(ranges, start_idx:length(signal))
    end

    return ranges
end

"""
    contiguous_values(signal)

Identify contiguous stretches of non-zero signal values and return `(range, value)` pairs.
"""
function contiguous_values(signal::AbstractVector{T}) where {T<:Real}
    ranges = Tuple{UnitRange{Int}, T}[]
    start_idx = 0
    in_range = false
    current_val = zero(T)

    for i in eachindex(signal)
        if signal[i] != zero(T)
            if !in_range
                start_idx = i
                current_val = signal[i]
                in_range = true
            elseif signal[i] != current_val
                push!(ranges, (start_idx:(i - 1), current_val))
                start_idx = i
                current_val = signal[i]
            end
        elseif in_range
            push!(ranges, (start_idx:(i - 1), current_val))
            in_range = false
        end
    end

    if in_range
        push!(ranges, (start_idx:length(signal), current_val))
    end

    return ranges
end

"""
    Base.show(io, chrom::ChromData)

Pretty-print a summary description for `ChromData`.
"""
function Base.show(io::IO, chrom::ChromData)
    print(io, "Chromosome $(chrom.name) with a $(length(chrom.signal)) bp signal")
end

"""
    Base.show(io, experiment::Experiment)
"""
function Base.show(io::IO, experiment::Experiment)
    print(io, "Experiment '$(experiment.name)' with $(length(experiment.samples)) samples")
end

"""
    Base.show(io, sample::SampleData)
"""
function Base.show(io::IO, sample::SampleData)
    reads = sample.n_reads == -1 ? "" : " and $(sample.n_reads) reads"
    print(io, "Sample '$(sample.name)' with $(length(sample.chroms)) sequences$reads")
end

"""
    _convert_signal_slice(signal_slice)

Attempt to convert numeric slices into the smallest unsigned integer vector that can represent them.
"""
function _convert_signal_slice(signal_slice::AbstractVector)
    signal_slice isa AbstractVector{Float64} && return signal_slice
    for T in (UInt8, UInt16, UInt32)
        try
            return convert(Vector{T}, signal_slice)
        catch InexactError
        end
    end
    return convert(Vector{UInt64}, signal_slice)
end

end # module GenomicData
