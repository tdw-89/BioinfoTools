module AlignmentUtils

"""
    pair_aln_perc(seq1, seq2, gap_char='-')

Calculate pairwise alignment percent identity between two aligned sequences.

# Arguments
- `seq1::AbstractString`: First aligned sequence
- `seq2::AbstractString`: Second aligned sequence (must be same length as `seq1`)
- `gap_char::Char`: Character representing gaps in the alignment

# Returns
- `Tuple{Int, Float64, Int, Float64}`: A tuple containing:
  1. Length of seq1 (excluding gaps)
  2. Percent identity from seq1 to seq2
  3. Length of seq2 (excluding gaps)
  4. Percent identity from seq2 to seq1
"""
function pair_aln_perc(seq1::AbstractString, seq2::AbstractString, gap_char::Char='-')
    if length(seq1) != length(seq2)
        error("Sequences must be of equal length for pairwise alignment percentage calculation.")
    end

    seq1_len = sum(c != gap_char for c in seq1)
    seq2_len = sum(c != gap_char for c in seq2)

    seq1_to_seq2 = sum(c1 == c2 && c1 != gap_char for (c1, c2) in zip(seq1, seq2)) / seq1_len * 100
    seq2_to_seq1 = sum(c1 == c2 && c2 != gap_char for (c1, c2) in zip(seq1, seq2)) / seq2_len * 100

    return (
        seq1_len, 
        seq1_to_seq2, 
        seq2_len, 
        seq2_to_seq1
        )
end

"""
    read_fasta_aln(file_path)

Read a FASTA alignment file and return sequences as a dictionary.

# Arguments
- `file_path::AbstractString`: Path to the FASTA alignment file

# Returns
- `Dict{String, String}`: Dictionary mapping sequence headers to sequences
"""
function read_fasta_aln(file_path::AbstractString)
    open(file_path, "r") do io
        records = read(io, String)
        entries = split(records, '>')[2:end]
        seq_dict = Dict{String, String}()
        for entry in entries
            lines = split(entry, '\n')
            header = lines[1]
            seq = join(lines[2:end], "")
            seq_dict[header] = seq
        end
        return seq_dict
    end
end

export pair_aln_perc, read_fasta_aln

end # module AlignmentUtils