# Graph representation for k-mers and their annotations
mutable struct MetaGraphIndex
    k::Int
    samples::Vector{String}
    kmer_to_id::Dict{UInt64, Int} # Map k-mer to row index
    adjacency::Vector{Int}        # Simple successor for RowDiff
    annotations::Vector{BitVector} # The "Colors" (Raw)
    
    MetaGraphIndex(k, samples) = new(k, samples, Dict(), Int[], BitVector[])
end

# Add a sample and its sequence to the graph index
function add_sample!(idx::MetaGraphIndex, sample_name::String, seq::String)
    # k = 5; samples = ["sample-1", "sample-2", "sample-3"]; idx = MetaGraphIndex(k, samples); sample_name = samples[1]; seq = "ATCTCGATCGATCGACTACG"
    # Find the sample ID (index) for the given sample name
    sample_id = findfirst(idx.samples .== sample_name)
    kmers = get_kmers(seq, idx.k)
    for (i, km) in enumerate(kmers)
        # i = 1; km = kmers[i]
        # Check if the k-mer is already in the index; if not, add it and initialize its annotation
        if !haskey(idx.kmer_to_id, km.data)
            push!(idx.annotations, BitVector(zeros(length(idx.samples))))
            idx.kmer_to_id[km.data] = length(idx.annotations)
            # Topology: track successor for RowDiff
            next_id = (i < length(kmers)) ? 0 : 1 # Simplified
            push!(idx.adjacency, next_id) 
        end
        # Update the annotation for the current k-mer and sample
        row_idx = idx.kmer_to_id[km.data]
        # Mark the presence of the k-mer in the sample by setting the corresponding bit
        idx.annotations[row_idx][sample_id] = true
    end
    # idx
end

# TODO: Implement graph traversal and compression logic (RowDiff) for the MetaGraphIndex

struct CompressedMetaGraph
    kmer_to_id::Dict{UInt64, Int}
    # RowDiff transformed matrix: stores XOR deltas
    diff_matrix::Vector{BitVector}
    anchors::BitVector # Which rows are stored raw (sinks or random samplings)
end

function compress(idx::MetaGraphIndex)
    n_rows = length(idx.annotations)
    diff_matrix = Vector{BitVector}(undef, n_rows)
    anchors = BitVector(zeros(n_rows))
    
    for i in 1:n_rows
        # If node has a successor in the graph, store the XOR delta
        successor_idx = (i < n_rows) ? i + 1 : 0 # In reality, use graph traversal
        
        if successor_idx != 0
            diff_matrix[i] = idx.annotations[i] .⊻ idx.annotations[successor_idx]
        else
            diff_matrix[i] = idx.annotations[i]
            anchors[i] = true
        end
    end
    return CompressedMetaGraph(idx.kmer_to_id, diff_matrix, anchors)
end