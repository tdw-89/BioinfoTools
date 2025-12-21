"""
BioinfoTools.jl bundles strongly typed genome models, high-throughput signal and peak
processing utilities, enrichment helpers, and paralog analysis tools into a single
package. Load the top-level module to access the submodules listed below or import the
ones you need individually.

# Submodules

- `GenomeTypes` - core data structures for scaffolds, contigs, genes, regions, and
  associated signal tracks.
- `GenomicData` - parsers and utilities for working with BAM/BED signals, replicate
  aggregation, and experiment management.
- `EnrichmentUtils` - helpers for aligning signal windows to genomic features and computing
  enrichment statistics.
- `LoadGFF` - functions for loading GFF/GFF3 annotations into the in-memory genome types.
- `ParalogUtils` - reciprocal-best-hit and gene family utilities powered by Graphs.jl.
"""
module BioinfoTools

# Submodules
include("genome_types.jl")
include("genomic_data.jl")
include("enrichment_utils.jl")
include("load_gff.jl")
include("paralog_utils.jl")
include("alignment_utils.jl")
include("misc_utils.jl")
using .GenomeTypes
using .GenomicData
using .EnrichmentUtils
using .LoadGFF
using .ParalogUtils
using .AlignmentUtils
using .MiscUtils

@doc raw"""
    GenomeTypes

Data structures representing scaffolds, contigs, genes, regulatory regions, and the
signals attached to them. These typed containers keep coordinates, annotations, and
signals bundled for downstream analyses.
""" GenomeTypes
@doc raw"""
    GenomicData

Tools for streaming BAM/BED data, aggregating replicates, binning peaks, and attaching
signals back to annotated genome objects.
""" GenomicData
@doc raw"""
    EnrichmentUtils

Helper routines for carving promoter/TSS/TES windows, validating requested regions, and
extracting signal slices suitable for downstream enrichment or plotting workflows.
""" EnrichmentUtils
@doc raw"""
    LoadGFF

Readers that convert one or more GFF/GFF3 files (plain or gzipped) into `RefGenome`
objects populated with genes, transcripts, and optional repeats.
""" LoadGFF
@doc raw"""
    ParalogUtils

Utilities for discovering reciprocal best hits, building paralog graphs, and exporting
gene family relationships for further network analysis.
""" ParalogUtils
export GenomeTypes, 
        GenomicData, 
        EnrichmentUtils, 
        LoadGFF, 
        ParalogUtils
end