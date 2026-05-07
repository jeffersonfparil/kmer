module kmer

include("io.jl")
include("graph.jl")
include("search.jl")

export NUCLEOTIDES, Kmer, get_kmers
export MetaGraphIndex, add_sample!, compress
export reconstruct_row, search

end
