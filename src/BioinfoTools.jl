module BioinfoTools
include("genome_types.jl")
include("genomic_data.jl")
include("enrichment_utils.jl")
include("load_gff.jl")
include("paralog_utils.jl")
using .GenomeTypes
using .GenomicData
using .EnrichmentUtils
using .LoadGFF
using .ParalogUtils
export GenomeTypes, GenomicData, EnrichmentUtils, LoadGFF, ParalogUtils
end
