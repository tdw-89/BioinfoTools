```@meta
CurrentModule = BioinfoTools
DocTestSetup = quote
    using BioinfoTools
    using BioinfoTools.LoadGFF
    using BioinfoTools.GenomicData
end
```

# BioinfoTools.jl

This package exposes a small set of submodules that cover genome types, data ingestion for BAM/BED
signals, enrichment helpers, GFF loaders, and paralog analysis. Load `BioinfoTools` and bring the
submodules you need into scope; the rest of this manual enumerates their exported APIs.

## Minimal workflow

```julia
using BioinfoTools
using BioinfoTools.LoadGFF
using BioinfoTools.GenomicData

ref = loadgenome("data/annotations.gff3"; feature_type = "gene")
experiment = getallcountvectors(["data/alignments/sample.bam"], exp_name = "demo")
addtogenes!(ref, experiment; peak_data = false)
```

```julia
using BioinfoTools.EnrichmentUtils
using BioinfoTools.ParalogUtils

promoter = getrange(first(last(ref.genes)), "promoter")
promoter_signal = getsiginrange(first(last(ref.genes)), promoter; peak_data = false)
paired = rbh(paralog_df; scoring = "max")
```
