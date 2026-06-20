"""
BioinfoTools.jl bundles strongly typed genome models, high-throughput signal and peak
processing utilities, enrichment helpers, and paralog analysis tools into a single
package. Load the top-level module to access the submodules listed below or import the
ones you need individually.

# Submodules

- `Types` - core data structures for scaffolds, contigs, genes, regions, and
  associated signal tracks.
- `Data` - parsers and utilities for working with BAM/BED signals, replicate
  aggregation, and experiment management.
- `Enrichment` - helpers for aligning signal windows to genomic features and computing
  enrichment statistics.
- `GFF` - functions for loading GFF/GFF3 annotations into the in-memory genome types.
- `Paralogs` - reciprocal-best-hit and gene family utilities powered by Graphs.jl.
"""
module BioinfoTools

# Submodules
include("types.jl")
include("data.jl")
include("enrichment.jl")
include("gff.jl")
include("paralogs.jl")
include("alignment.jl")
include("repeats.jl")
include("misc.jl")
using .Types
using .Data
using .Enrichment
using .GFF
using .Paralogs
using .AlignmentUtils
using .RepeatUtils
using .MiscUtils

@doc raw"""
    Types

Data structures representing scaffolds, contigs, genes, regulatory regions, and the
signals attached to them. These typed containers keep coordinates, annotations, and
signals bundled for downstream analyses.
""" Types
@doc raw"""
    Data

Tools for streaming BAM/BED data, aggregating replicates, binning peaks, and attaching
signals back to annotated genome objects.
""" Data
@doc raw"""
    Enrichment

Helper routines for carving promoter/TSS/TES windows, validating requested regions, and
extracting signal slices suitable for downstream enrichment or plotting workflows.
""" Enrichment
@doc raw"""
    GFF

Readers that convert one or more GFF/GFF3 files (plain or gzipped) into `RefGenome`
objects populated with genes, transcripts, and optional repeats.
""" GFF
@doc raw"""
    Paralogs

Utilities for discovering reciprocal best hits, building paralog graphs, and exporting
gene family relationships for further network analysis.
""" Paralogs
export Types, 
        Data, 
        Enrichment, 
        GFF, 
        Paralogs,
        AlignmentUtils,
        RepeatUtils,
        MiscUtils
end