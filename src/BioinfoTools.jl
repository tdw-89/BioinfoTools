module BioinfoTools
include("genome_types.jl")
include("genomic_data.jl")
include("enrichment_utils.jl")
include("load_gff.jl")
using .GenomeTypes
using .GenomicData
using .EnrichmentUtils
using .LoadGFF
export GenomeTypes, GenomicData, EnrichmentUtils, LoadGFF
end
