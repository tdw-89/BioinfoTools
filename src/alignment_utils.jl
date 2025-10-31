module AlignmentUtils

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
end # module AlignmentUtils