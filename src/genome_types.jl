module GenomeTypes

using BioSequences

abstract type Feature end
abstract type RegElement <: Feature end

const _ScaffoldNameUnion = Union{String, Char}
const _MaybeInt = Union{Int, Missing}
const _MaybeString = Union{String, Missing}
const _MaybeChar = Union{Char, Missing}
const _MaybeFeatureVector = Union{AbstractVector{<:Feature}, Missing}
const _MaybeRegElementVector = Union{AbstractVector{<:RegElement}, Missing}
const _SignalCollection = Union{
    Vector{Union{Vector{UInt8}, Vector{UInt16}, Vector{UInt32}, Vector{UInt64}, Vector{Float64}}},
    Nothing,
}
const _BitSignalCollection = Union{Vector{BitVector}, Nothing}
const _SampleCollection = Union{Vector{String}, Nothing}
const _StringVectorMissing = Union{Vector{String}, Missing}
const _FloatVectorMissing = Union{Vector{Float64}, Missing}
const _SequenceOrNothing = Union{BioSequences.BioSequence, Nothing}
const _SequenceOrMissing = Union{BioSequences.BioSequence, Missing}
const _FeatureAnnotationDict = Union{Dict{String, Vector{Feature}}, Missing}

mutable struct Scaffold{NameType, ContigType, GeneType, RepeatType, StartType, EndType, LevelType} <: Feature
    name::NameType
    contigs::ContigType
    genes::GeneType
    repeats::RepeatType
    scaffold_start::StartType
    scaffold_end::EndType
    level::LevelType
end

function Scaffold(
    name::Name,
    contigs::Contigs=missing,
    genes::Genes=missing,
    repeats::Repeats=missing,
    scaffold_start::Start=missing,
    scaffold_end::End=missing,
    level::Level=missing,
) where {
    Name<:_ScaffoldNameUnion,
    Contigs<:_MaybeFeatureVector,
    Genes<:_MaybeFeatureVector,
    Repeats<:_MaybeFeatureVector,
    Start<:_MaybeInt,
    End<:_MaybeInt,
    Level<:_MaybeString,
}
    Scaffold{Name, Contigs, Genes, Repeats, Start, End, Level}(
        name,
        contigs,
        genes,
        repeats,
        scaffold_start,
        scaffold_end,
        level,
    )
end

const _MaybeScaffold = Union{Scaffold, Missing}

mutable struct Contig{
    ScaffoldType,
    JunkType,
    EnhancerType,
    GeneType,
    StartType,
    EndType,
} <: Feature
    scaffold::ScaffoldType
    junk::JunkType
    enhancers::EnhancerType
    genes::GeneType
    contig_start::StartType
    contig_end::EndType
end

function Contig(
    scaffold::ScaffoldType=missing,
    junk::JunkType=missing,
    enhancers::EnhancerType=missing,
    genes::GeneType=missing,
    contig_start::Start=missing,
    contig_end::End=missing,
) where {
    ScaffoldType<:_MaybeScaffold,
    JunkType<:_MaybeFeatureVector,
    EnhancerType<:_MaybeFeatureVector,
    GeneType<:_MaybeFeatureVector,
    Start<:_MaybeInt,
    End<:_MaybeInt,
}
    Contig{ScaffoldType, JunkType, EnhancerType, GeneType, Start, End}(
        scaffold,
        junk,
        enhancers,
        genes,
        contig_start,
        contig_end,
    )
end

const _MaybeContig = Union{Contig, Missing}

struct Enhancer{ContigType} <: RegElement
    contig::ContigType
end

function Enhancer(contig::C=missing) where {C<:_MaybeContig}
    Enhancer{C}(contig)
end

const _MaybeScaffoldOrString = Union{Scaffold, String, Missing}
const _MaybeContigOrString = Union{Contig, String, Missing}

struct Gene{
    ScaffoldType,
    ContigType,
    NameType,
    IdType,
    StrandType,
    CresType,
    TssType,
    TesType,
    IntronsType,
    ExonsType,
    CdsStartType,
    CdsEndType,
    RnaType,
    SegmentsType,
    RegionsType,
    AnnotationsType,
    GeneStartType,
    GeneEndType,
    SequenceType,
    SignalsType,
    BinSignalsType,
    SamplesType,
} <: Feature
    scaffold::ScaffoldType
    contig::ContigType
    name::NameType
    id::IdType
    strand::StrandType
    cres::CresType
    tss::TssType
    tes::TesType
    introns::IntronsType
    exons::ExonsType
    cds_start::CdsStartType
    cds_end::CdsEndType
    rnas::RnaType
    segments::SegmentsType
    regions::RegionsType
    annotations::AnnotationsType
    gene_start::GeneStartType
    gene_end::GeneEndType
    sequence::SequenceType
    signals::SignalsType
    binsignals::BinSignalsType
    samples::SamplesType
end

function Gene(
    scaffold::ScaffoldType=missing,
    contig::ContigType=missing,
    name::NameType=missing,
    id::IdType=missing,
    strand::StrandType=missing,
    cres::CresType=missing,
    tss::TssType=missing,
    tes::TesType=missing,
    introns::IntronsType=missing,
    exons::ExonsType=missing,
    cds_start::CdsStartType=missing,
    cds_end::CdsEndType=missing,
    rnas::RnaType=missing,
    segments::SegmentsType=missing,
    regions::RegionsType=missing,
    annotations::AnnotationsType=missing,
    gene_start::GeneStartType=missing,
    gene_end::GeneEndType=missing,
    sequence::SequenceType=nothing,
    signals::SignalsType=nothing,
    binsignals::BinSignalsType=nothing,
    samples::SamplesType=nothing,
) where {
    ScaffoldType<:_MaybeScaffoldOrString,
    ContigType<:_MaybeContigOrString,
    NameType<:_MaybeString,
    IdType<:_MaybeString,
    StrandType<:_MaybeChar,
    CresType<:_MaybeRegElementVector,
    TssType<:_MaybeInt,
    TesType<:_MaybeInt,
    IntronsType<:_MaybeFeatureVector,
    ExonsType<:_MaybeFeatureVector,
    CdsStartType<:_MaybeInt,
    CdsEndType<:_MaybeInt,
    RnaType<:_MaybeFeatureVector,
    SegmentsType<:_MaybeFeatureVector,
    RegionsType<:_MaybeFeatureVector,
    AnnotationsType<:_FeatureAnnotationDict,
    GeneStartType<:_MaybeInt,
    GeneEndType<:_MaybeInt,
    SequenceType<:_SequenceOrNothing,
    SignalsType<:_SignalCollection,
    BinSignalsType<:_BitSignalCollection,
    SamplesType<:_SampleCollection,
}
    Gene{
        ScaffoldType,
        ContigType,
        NameType,
        IdType,
        StrandType,
        CresType,
        TssType,
        TesType,
        IntronsType,
        ExonsType,
        CdsStartType,
        CdsEndType,
        RnaType,
        SegmentsType,
        RegionsType,
        AnnotationsType,
        GeneStartType,
        GeneEndType,
        SequenceType,
        SignalsType,
        BinSignalsType,
        SamplesType,
    }(
        scaffold,
        contig,
        name,
        id,
        strand,
        cres,
        tss,
        tes,
        introns,
        exons,
        cds_start,
        cds_end,
        rnas,
        segments,
        regions,
        annotations,
        gene_start,
        gene_end,
        sequence,
        signals,
        binsignals,
        samples,
    )
end

const _MaybeGeneVector = Union{AbstractVector{<:Gene}, Missing}

struct Promoter{
    GeneType,
    StrandType,
    SegmentsType,
    StartType,
    EndType,
    SequenceType,
    SignalsType,
    BinSignalsType,
    SamplesType,
} <: RegElement
    genes::GeneType
    strand::StrandType
    segments::SegmentsType
    promoter_start::StartType
    promoter_end::EndType
    sequence::SequenceType
    signals::SignalsType
    binsignals::BinSignalsType
    samples::SamplesType
end

function Promoter(
    genes::Genes=missing,
    strand::Strand=missing,
    segments::Segments=missing,
    promoter_start::Start=missing,
    promoter_end::End=missing,
    sequence::Sequence=nothing,
    signals::Signals=nothing,
    binsignals::BinSignals=nothing,
    samples::Samples=nothing,
) where {
    Genes<:_MaybeGeneVector,
    Strand<:_MaybeChar,
    Segments<:_MaybeFeatureVector,
    Start<:_MaybeInt,
    End<:_MaybeInt,
    Sequence<:_SequenceOrNothing,
    Signals<:_SignalCollection,
    BinSignals<:_BitSignalCollection,
    Samples<:_SampleCollection,
}
    Promoter{
        Genes,
        Strand,
        Segments,
        Start,
        End,
        Sequence,
        Signals,
        BinSignals,
        Samples,
    }(
        genes,
        strand,
        segments,
        promoter_start,
        promoter_end,
        sequence,
        signals,
        binsignals,
        samples,
    )
end

struct Intron{
    GeneType,
    RnaType,
    StartType,
    EndType,
} <: Feature
    gene::GeneType
    rna::RnaType
    intron_start::StartType
    intron_end::EndType
end

function Intron(
    gene::GeneType=missing,
    rna::RnaType=missing,
    intron_start::Start=missing,
    intron_end::End=missing,
) where {
    GeneType<:Union{Gene, Missing},
    RnaType<:Union{Feature, Missing},
    Start<:_MaybeInt,
    End<:_MaybeInt,
}
    Intron{GeneType, RnaType, Start, End}(gene, rna, intron_start, intron_end)
end

struct Exon{
    GeneType,
    RnaType,
    StartType,
    EndType,
} <: Feature
    gene::GeneType
    rna::RnaType
    exon_start::StartType
    exon_end::EndType
end

function Exon(
    gene::GeneType=missing,
    rna::RnaType=missing,
    exon_start::Start=missing,
    exon_end::End=missing,
) where {
    GeneType<:Union{Gene, Missing},
    RnaType<:Union{Feature, Missing},
    Start<:_MaybeInt,
    End<:_MaybeInt,
}
    Exon{GeneType, RnaType, Start, End}(gene, rna, exon_start, exon_end)
end

struct RNA{
    IdType,
    TypeType,
    GeneType,
    ExonType,
    IntronType,
    SamplesType,
    StartType,
    EndType,
    ExpressionType,
} <: Feature
    id::IdType
    type::TypeType
    gene::GeneType
    exons::ExonType
    introns::IntronType
    samples::SamplesType
    rna_start::StartType
    rna_end::EndType
    expression::ExpressionType
end

function RNA(
    id::IdType=missing,
    type::TypeType=missing,
    gene::GeneType=missing,
    exons::ExonType=missing,
    introns::IntronType=missing,
    samples::SamplesType=missing,
    rna_start::StartType=missing,
    rna_end::EndType=missing,
    expression::ExpressionType=missing,
) where {
    IdType<:_MaybeString,
    TypeType<:_MaybeString,
    GeneType<:Union{Gene, Missing},
    ExonType<:Union{AbstractVector{<:Exon}, Missing},
    IntronType<:Union{AbstractVector{<:Intron}, Missing},
    SamplesType<:_StringVectorMissing,
    StartType<:_MaybeInt,
    EndType<:_MaybeInt,
    ExpressionType<:_FloatVectorMissing,
}
    RNA{
        IdType,
        TypeType,
        GeneType,
        ExonType,
        IntronType,
        SamplesType,
        StartType,
        EndType,
        ExpressionType,
    }(
        id,
        type,
        gene,
        exons,
        introns,
        samples,
        rna_start,
        rna_end,
        expression,
    )
end

struct Annotation{
    ElementType,
    AnnotationType,
} <: Feature
    elements::ElementType
    annotation::AnnotationType
end

function Annotation(
    elements::ElementsType=missing,
    annotation::AnnotationType="",
) where {
    ElementsType<:_FeatureAnnotationDict,
    AnnotationType<:AbstractString,
}
    Annotation{ElementsType, AnnotationType}(elements, annotation)
end

const _AnnotationDict = Union{Dict{String, Vector{Annotation}}, Missing}

struct Region{
    ScaffoldType,
    ContigType,
    StartType,
    EndType,
    AnnotationType,
    SignalsType,
    BinSignalsType,
    SamplesType,
} <: Feature
    scaffold::ScaffoldType
    contig::ContigType
    region_start::StartType
    region_end::EndType
    annotations::AnnotationType
    signals::SignalsType
    binsignals::BinSignalsType
    samples::SamplesType
end

function Region(
    scaffold::ScaffoldType=missing,
    contig::ContigType=missing,
    region_start::StartType=missing,
    region_end::EndType=missing,
    annotations::AnnotationType=missing,
    signals::SignalsType=nothing,
    binsignals::BinSignalsType=nothing,
    samples::SamplesType=nothing,
) where {
    ScaffoldType<:_MaybeScaffold,
    ContigType<:_MaybeContig,
    StartType<:_MaybeInt,
    EndType<:_MaybeInt,
    AnnotationType<:_AnnotationDict,
    SignalsType<:_SignalCollection,
    BinSignalsType<:_BitSignalCollection,
    SamplesType<:_SampleCollection,
}
    Region{
        ScaffoldType,
        ContigType,
        StartType,
        EndType,
        AnnotationType,
        SignalsType,
        BinSignalsType,
        SamplesType,
    }(
        scaffold,
        contig,
        region_start,
        region_end,
        annotations,
        signals,
        binsignals,
        samples,
    )
end

const _RegionVector = Union{AbstractVector{<:Region}, Missing}
const _SequenceOrMissingRepeat = Union{BioSequences.BioSequence, Missing}

struct Repeat{
    ScaffoldType,
    ContigType,
    RegionType,
    FamilyType,
    TypeType,
    StartType,
    EndType,
    SequenceType,
    SignalsType,
    BinSignalsType,
    SamplesType,
} <: Feature
    scaffold::ScaffoldType
    contig::ContigType
    regions::RegionType
    family::FamilyType
    type::TypeType
    repeat_start::StartType
    repeat_end::EndType
    sequence::SequenceType
    signals::SignalsType
    binsignals::BinSignalsType
    samples::SamplesType
end

function Repeat(
    scaffold::ScaffoldType=missing,
    contig::ContigType=missing,
    regions::RegionType=missing,
    family::FamilyType=missing,
    type::TypeType=missing,
    repeat_start::StartType=missing,
    repeat_end::EndType=missing,
    sequence::SequenceType=missing,
    signals::SignalsType=nothing,
    binsignals::BinSignalsType=nothing,
    samples::SamplesType=nothing,
) where {
    ScaffoldType<:_MaybeScaffold,
    ContigType<:_MaybeContig,
    RegionType<:_RegionVector,
    FamilyType<:_MaybeString,
    TypeType<:_MaybeString,
    StartType<:_MaybeInt,
    EndType<:_MaybeInt,
    SequenceType<:_SequenceOrMissingRepeat,
    SignalsType<:_SignalCollection,
    BinSignalsType<:_BitSignalCollection,
    SamplesType<:_SampleCollection,
}
    Repeat{
        ScaffoldType,
        ContigType,
        RegionType,
        FamilyType,
        TypeType,
        StartType,
        EndType,
        SequenceType,
        SignalsType,
        BinSignalsType,
        SamplesType,
    }(
        scaffold,
        contig,
        regions,
        family,
        type,
        repeat_start,
        repeat_end,
        sequence,
        signals,
        binsignals,
        samples,
    )
end

struct Segment{
    PrevType,
    StartType,
    EndType,
    AnnotationType,
    NextType,
    SignalsType,
    BinSignalsType,
    SamplesType,
} <: Feature
    prev_segment::PrevType
    segment_start::StartType
    segment_end::EndType
    annotations::AnnotationType
    next_segment::NextType
    signals::SignalsType
    binsignals::BinSignalsType
    samples::SamplesType
end

function Segment(
    prev_segment::PrevType=missing,
    segment_start::StartType=missing,
    segment_end::EndType=missing,
    annotations::AnnotationType=missing,
    next_segment::NextType=missing,
    signals::SignalsType=nothing,
    binsignals::BinSignalsType=nothing,
    samples::SamplesType=nothing,
) where {
    PrevType<:_MaybeFeatureVector,
    StartType<:_MaybeInt,
    EndType<:_MaybeInt,
    AnnotationType<:_AnnotationDict,
    NextType<:_MaybeFeatureVector,
    SignalsType<:_SignalCollection,
    BinSignalsType<:_BitSignalCollection,
    SamplesType<:_SampleCollection,
}
    Segment{
        PrevType,
        StartType,
        EndType,
        AnnotationType,
        NextType,
        SignalsType,
        BinSignalsType,
        SamplesType,
    }(
        prev_segment,
        segment_start,
        segment_end,
        annotations,
        next_segment,
        signals,
        binsignals,
        samples,
    )
end

mutable struct RefGenome{
    ScaffoldDictType,
    ContigDictType,
    EnhancerVectorType,
    PromoterVectorType,
    GeneTupleType,
    IntronVectorType,
    ExonVectorType,
    RnaVectorType,
    SegmentDictType,
    RepeatVectorType,
    RegionVectorType,
    AnnotationVectorType,
}
    scaffolds::ScaffoldDictType
    contigs::ContigDictType
    enhancers::EnhancerVectorType
    promoters::PromoterVectorType
    genes::GeneTupleType
    introns::IntronVectorType
    exons::ExonVectorType
    rnas::RnaVectorType
    segments::SegmentDictType
    repeats::RepeatVectorType
    regions::RegionVectorType
    annotations::AnnotationVectorType
end

function RefGenome(
    scaffolds::Dict{String, Scaffold},
    contigs::Dict{String, Contig},
    enhancers::Vector{Enhancer},
    promoters::Vector{Promoter},
    genes::Tuple{Vector{String}, Vector{Gene}},
    introns::Vector{Intron},
    exons::Vector{Exon},
    rnas::Vector{RNA},
    segments::Dict{String, Vector{Segment}},
    repeats::Vector{Repeat},
    regions::Vector{Region},
    annotations::Vector{Annotation},
)
    RefGenome{
        typeof(scaffolds),
        typeof(contigs),
        typeof(enhancers),
        typeof(promoters),
        typeof(genes),
        typeof(introns),
        typeof(exons),
        typeof(rnas),
        typeof(segments),
        typeof(repeats),
        typeof(regions),
        typeof(annotations),
    }(
        scaffolds,
        contigs,
        enhancers,
        promoters,
        genes,
        introns,
        exons,
        rnas,
        segments,
        repeats,
        regions,
        annotations,
    )
end

RefGenome() = RefGenome(
    Dict{String, Scaffold}(),
    Dict{String, Contig}(),
    Vector{Enhancer}(),
    Vector{Promoter}(),
    (Vector{String}(), Vector{Gene}()),
    Vector{Intron}(),
    Vector{Exon}(),
    Vector{RNA}(),
    Dict{String, Vector{Segment}}(),
    Vector{Repeat}(),
    Vector{Region}(),
    Vector{Annotation}(),
)

abstract type RangeAnchor end
struct TSS <: RangeAnchor end
struct REGION <: RangeAnchor end
struct TES <: RangeAnchor end

abstract type OffSetType end
struct INTEGER <: OffSetType end
struct PERCENTAGE <: OffSetType end

struct GeneRange
    range_start::RangeAnchor
    range_stop::RangeAnchor
    start_offset::Int
    stop_offset::Int
    start_offset_type::OffSetType
    stop_offset_type::OffSetType

    function GeneRange(
        start::RangeAnchor,
        stop::RangeAnchor,
        start_offset::Int,
        stop_offset::Int,
        start_offset_type::OffSetType,
        stop_offset_type::OffSetType,
    )
        if isa(start_offset_type, PERCENTAGE) && !isa(start, TSS)
            error()
        end

        if isa(stop_offset_type, PERCENTAGE) && !isa(stop, TSS)
            error()
        end

        new(start, stop, start_offset, stop_offset, start_offset_type, stop_offset_type)
    end
end

GeneRange(start::RangeAnchor, stop::RangeAnchor) =
    GeneRange(start, stop, 0, 0, INTEGER(), INTEGER())

GeneRange(start::RangeAnchor, stop::RangeAnchor, start_offset::Int, stop_offset::Int) =
    GeneRange(start, stop, start_offset, stop_offset, INTEGER(), INTEGER())

const Junk = Feature
const mRNA = RNA

Base.show(io::IO, x::Feature) = begin
    type = typeof(x)
    fields = fieldnames(type)
    fields_missing = [ismissing(getfield(x, i)) for i in fields]
    fields_present = fields[.!fields_missing]
    fields_missing = fields[fields_missing]

    if any(:id .== fields_present)
        id = getfield(x, :id)
        print(io, "$type: '$id', with fields: ")
        for field in fields_present
            print(io, "'$field' ")
        end
    else
        print(io, "$type with fields: ")
        for field in fields_present
            print(io, "'$field' ")
        end
    end

    print(io, "\n     missing fields: ")
    for field in fields_missing
        print(io, "'$field' ")
    end
end

Base.show(io::IO, x::Vector{Feature}) = begin
    type = typeof(x)
    len = length(x)
    println(io, "$type with $len items")
end

export Feature,
       RegElement,
       Scaffold,
       Contig,
       Junk,
       Enhancer,
       Gene,
       Promoter,
       mRNA,
       Intron,
       Exon,
       Segment,
       Repeat,
       Annotation,
       Region,
       RNA,
       RefGenome,
       RangeAnchor,
       TSS,
       REGION,
       TES,
       OffSetType,
       INTEGER,
       PERCENTAGE,
       GeneRange

end # module
