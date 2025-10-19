module BioinfoTools

include("genome_types.jl")
include("genomic_data.jl")
include("enrichment_utils.jl")

using .GenomeTypes
using .GenomicData
using .EnrichmentUtils

export GenomeTypes, GenomicData, EnrichmentUtils

end
