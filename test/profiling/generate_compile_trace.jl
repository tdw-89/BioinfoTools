using Pkg
Pkg.activate(joinpath(@__DIR__, "../../"))

using BioinfoTools

test_genome = GFF.loadgenome(joinpath(@__DIR__, "../data/genomic.gff.gz"), feature_type="all")
for gene in test_genome.genes[2]
    println("Processing gene: ", gene.id)
end