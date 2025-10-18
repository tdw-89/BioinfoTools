module GenomeTypes

using BioSequences

const Optional{T} = Union{Missing, T}
const OptionalOrNothing{T} = Union{Nothing, T}
const SignalElement = Union{UInt8, UInt16, UInt32, UInt64, Float64}

"""
    Feature

Abstract supertype for all genomic entities represented in this package.
Concrete feature types capture different biological concepts such as
scaffolds, genes, and regulatory elements.
"""
abstract type Feature end

"""
    RegElement

Abstract subtype of [`Feature`](@ref) for regulatory elements such as
enhancers and promoters.
"""
abstract type RegElement <: Feature end

const OptionalFeatureCollection = Union{Missing, AbstractVector{<:Feature}}
const OptionalRegElementCollection = Union{Missing, AbstractVector{<:RegElement}}
const OptionalSignalCollection =
    OptionalOrNothing{AbstractVector{<:AbstractVector{<:SignalElement}}}
const OptionalBitVectorCollection = OptionalOrNothing{AbstractVector{<:BitVector}}
const OptionalSampleCollection = OptionalOrNothing{AbstractVector{<:AbstractString}}
const OptionalAnnotationDict = Union{Missing, Dict{String, Vector{Feature}}}

"""
    Scaffold{N,C,G,R,S,E,L} <: Feature

Container for scaffold-level annotations and relationships between features.

# Type Parameters
- `N`: identifier type, typically `AbstractString` or `Char`.
- `C`: collection type used for contigs (allows `missing`).
- `G`: collection type used for genes (allows `missing`).
- `R`: collection type used for repeats (allows `missing`).
- `S`: coordinate type for scaffold start (allows `missing`).
- `E`: coordinate type for scaffold end (allows `missing`).
- `L`: descriptive level label type (allows `missing`).

# Fields
- `name::N`: scaffold identifier.
- `contigs::C`: contigs associated with the scaffold.
- `genes::G`: genes contained within the scaffold.
- `repeats::R`: repetitive elements located on the scaffold.
- `scaffold_start::S`: genomic start coordinate.
- `scaffold_end::E`: genomic end coordinate.
- `level::L`: assembly level or status label.
"""
mutable struct Scaffold{
    N<:Union{AbstractString,Char},
    C<:OptionalFeatureCollection,
    G<:OptionalFeatureCollection,
    R<:OptionalFeatureCollection,
    S<:Optional{Integer},
    E<:Optional{Integer},
    L<:Optional{AbstractString},
} <: Feature
    name::N
    contigs::C
    genes::G
    repeats::R
    scaffold_start::S
    scaffold_end::E
    level::L
end

"""
    Contig{S,J,E,G,T,U} <: Feature

Represents a contig within a scaffold and stores related genomic features.

# Type Parameters
- `S`: parent scaffold reference type (allows `missing`).
- `J`: collection of miscellaneous child features.
- `E`: collection of enhancers.
- `G`: collection of genes.
- `T`: coordinate type for contig start.
- `U`: coordinate type for contig end.
"""
mutable struct Contig{
    S<:Optional{Scaffold},
    J<:OptionalFeatureCollection,
    E<:OptionalRegElementCollection,
    G<:OptionalFeatureCollection,
    T<:Optional{Integer},
    U<:Optional{Integer},
} <: Feature
    scaffold::S
    junk::J
    enhancers::E
    genes::G
    contig_start::T
    contig_end::U
end

"""
    Enhancer{C} <: RegElement

Regulatory enhancer linked to a parent contig.

# Type Parameters
- `C`: parent contig reference type (allows `missing`).
"""
struct Enhancer{C<:Optional{Contig}} <: RegElement
    contig::C
end

"""
    Gene{S,C,N,I,Str,Cres,TSS,TES,Intr,Exon,CdsS,CdsE,Rna,Seg,Reg,Anno,GStart,GEnd,Seq,Sig,BinSig,Samples} <: Feature

Encapsulates gene-level annotations, relationships, and experimental signals.

# Fields
- `scaffold::S`: scaffold identifier or link.
- `contig::C`: contig identifier or link.
- `name::N`: gene symbol or name.
- `id::I`: stable gene identifier.
- `strand::Str`: genomic strand (`'+'` / `'-'`).
- `cres::Cres`: regulatory elements associated with the gene.
- `tss::TSS`: transcription start site coordinate.
- `tes::TES`: transcription end site coordinate.
- `introns::Intr`: intronic segments.
- `exons::Exon`: exonic segments.
- `cds_start::CdsS`: coding sequence start.
- `cds_end::CdsE`: coding sequence end.
- `rnas::Rna`: RNA isoforms transcribed from the gene.
- `segments::Seg`: other genomic segments linked to the gene.
- `regions::Reg`: genomic regions overlapping the gene.
- `annotations::Anno`: arbitrary annotation dictionary.
- `gene_start::GStart`: genomic start coordinate.
- `gene_end::GEnd`: genomic end coordinate.
- `sequence::Seq`: primary sequence representation.
- `signals::Sig`: quantitative signal tracks.
- `binsignals::BinSig`: binned binary signal tracks.
- `samples::Samples`: sample or condition labels.
"""
struct Gene{
    S<:Union{Missing,Scaffold,AbstractString},
    C<:Union{Missing,Contig,AbstractString},
    N<:Optional{AbstractString},
    I<:Optional{AbstractString},
    Str<:Optional{Char},
    Cres<:OptionalRegElementCollection,
    TSS<:Optional{Integer},
    TES<:Optional{Integer},
    Intr<:OptionalFeatureCollection,
    Exon<:OptionalFeatureCollection,
    CdsS<:Optional{Integer},
    CdsE<:Optional{Integer},
    Rna<:OptionalFeatureCollection,
    Seg<:OptionalFeatureCollection,
    Reg<:OptionalFeatureCollection,
    Anno<:OptionalAnnotationDict,
    GStart<:Optional{Integer},
    GEnd<:Optional{Integer},
    Seq<:OptionalOrNothing{BioSequences.BioSequence},
    Sig<:OptionalSignalCollection,
    BinSig<:OptionalBitVectorCollection,
    Samples<:OptionalSampleCollection,
} <: Feature
    scaffold::S
    contig::C
    name::N
    id::I
    strand::Str
    cres::Cres
    tss::TSS
    tes::TES
    introns::Intr
    exons::Exon
    cds_start::CdsS
    cds_end::CdsE
    rnas::Rna
    segments::Seg
    regions::Reg
    annotations::Anno
    gene_start::GStart
    gene_end::GEnd
    sequence::Seq
    signals::Sig
    binsignals::BinSig
    samples::Samples
end

"""
    Promoter{G,Str,Seg,Start,End,Seq,Sig,BinSig,Samples} <: RegElement

Promoter region tied to one or more genes with optional experimental signals.
"""
struct Promoter{
    G<:Optional{AbstractVector{<:Gene}},
    Str<:Optional{Char},
    Seg<:OptionalFeatureCollection,
    Start<:Optional{Integer},
    End<:Optional{Integer},
    Seq<:OptionalOrNothing{BioSequences.BioSequence},
    Sig<:OptionalSignalCollection,
    BinSig<:OptionalBitVectorCollection,
    Samples<:OptionalSampleCollection,
} <: RegElement
    genes::G
    strand::Str
    segments::Seg
    promoter_start::Start
    promoter_end::End
    sequence::Seq
    signals::Sig
    binsignals::BinSig
    samples::Samples
end

"""
    Intron{G,R,S,E} <: Feature

Intron segment belonging to a gene or RNA transcript.
"""
struct Intron{
    G<:Optional{Gene},
    R<:Optional{Feature},
    S<:Optional{Integer},
    E<:Optional{Integer},
} <: Feature
    gene::G
    rna::R
    intron_start::S
    intron_end::E
end

"""
    Exon{G,R,S,E} <: Feature

Exon segment belonging to a gene or RNA transcript.
"""
struct Exon{
    G<:Optional{Gene},
    R<:Optional{Feature},
    S<:Optional{Integer},
    E<:Optional{Integer},
} <: Feature
    gene::G
    rna::R
    exon_start::S
    exon_end::E
end

"""
    RNA{I,T,G,E,Intr,Samples,Start,End,Expr} <: Feature

Represents an RNA transcript derived from a gene.
"""
struct RNA{
    I<:Optional{AbstractString},
    T<:Optional{AbstractString},
    G<:Optional{Gene},
    E<:Union{Missing, AbstractVector{<:Exon}},
    Intr<:Union{Missing, AbstractVector{<:Intron}},
    Samples<:Optional{AbstractVector{<:AbstractString}},
    Start<:Optional{Integer},
    End<:Optional{Integer},
    Expr<:Optional{AbstractVector{<:AbstractFloat}},
} <: Feature
    id::I
    type::T
    gene::G
    exons::E
    introns::Intr
    samples::Samples
    rna_start::Start
    rna_end::End
    expression::Expr
end

"""
    Annotation{E,A} <: Feature

Generic annotation container for grouping features by label.
"""
struct Annotation{
    E<:OptionalAnnotationDict,
    A<:AbstractString,
} <: Feature
    elements::E
    annotation::A
end

const OptionalAnnotationCollection = Union{Missing, Dict{String, Vector{Annotation}}}

"""
    Region{S,C,RStart,REnd,Anno,Sig,BinSig,Samples} <: Feature

Arbitrary genomic region with optional annotations and signal data.
"""
struct Region{
    S<:Optional{Scaffold},
    C<:Optional{Contig},
    RStart<:Optional{Integer},
    REnd<:Optional{Integer},
    Anno<:OptionalAnnotationCollection,
    Sig<:OptionalSignalCollection,
    BinSig<:OptionalBitVectorCollection,
    Samples<:OptionalSampleCollection,
} <: Feature
    scaffold::S
    contig::C
    region_start::RStart
    region_end::REnd
    annotations::Anno
    signals::Sig
    binsignals::BinSig
    samples::Samples
end

"""
    Repeat{S,C,R,F,T,Start,End,Seq,Sig,BinSig,Samples} <: Feature

Repetitive genomic element with optional sequence and signal annotations.
"""
struct Repeat{
    S<:Optional{Scaffold},
    C<:Optional{Contig},
    R<:Optional{AbstractVector{<:Region}},
    F<:Optional{AbstractString},
    T<:Optional{AbstractString},
    Start<:Optional{Integer},
    End<:Optional{Integer},
    Seq<:Optional{BioSequences.BioSequence},
    Sig<:OptionalSignalCollection,
    BinSig<:OptionalBitVectorCollection,
    Samples<:OptionalSampleCollection,
} <: Feature
    scaffold::S
    contig::C
    regions::R
    family::F
    type::T
    repeat_start::Start
    repeat_end::End
    sequence::Seq
    signals::Sig
    binsignals::BinSig
    samples::Samples
end

"""
    Segment{Prev,Start,End} <: Feature

Segment of genomic sequence, optionally linked to preceding segments.
"""
struct Segment{
    Prev<:OptionalFeatureCollection,
    Start<:Optional{Integer},
    End<:Optional{Integer},
} <: Feature
    prev_segment::Prev
    segment_start::Start
    segment_end::End
end

export Feature, RegElement, Scaffold, Contig, Enhancer, Gene, Promoter,
       Intron, Exon, RNA, Annotation, Region, Repeat, Segment

end # module
