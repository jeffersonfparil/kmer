# Graph representation for k-mers and their annotations
"""
MetaGraphIndex(k, samples)

A mutable graph index for k-mers and their annotations.

Fields
- `k`: k-mer size.
- `samples`: list of sample names.
- `kmer_to_id`: mapping from k-mer hashes to row indices.
- `adjacency`: simple successor indices for RowDiff compression.
- `annotations`: bit vectors representing sample presence per k-mer.
"""
mutable struct MetaGraphIndex
    k::Int64
    samples::Vector{String}
    kmer_to_id::Dict{UInt64, Int64} # Map k-mer to row index
    adjacency::Vector{Int64}        # Simple successor for RowDiff
    annotations::Vector{BitVector} # The "Colors" (Raw); each row corresponds to a k-mer, each column to a sample

    MetaGraphIndex(k::Int64, samples::Vector{String}) = new(k, samples, Dict(), Int64[], BitVector[])
end

"""
add_sample!(idx, sample_name, seq)

Add a sample and its sequence to the graph index.
Updates the k-mer mapping, adjacency, and annotations for the given sample.
"""
# Add a sample and its sequence to the graph index
function add_sample!(idx::MetaGraphIndex, sample_name::String, seq::String)
    # samples = ["sample_1", "sample_2", "sample_3"]; idx = MetaGraphIndex(5, samples); 
    # sample_name = "sample_1"; seq = "ATCTCGATCGATCGACTACG"
    # sample_name = "sample_2"; seq = "ATCCCGATCGAACGACTTCG"
    # sample_name = "sample_3"; seq = "ATCTTGATCCATCGACTTCG"

    # Find the sample ID (index) for the given sample name
    sample_id = findfirst(idx.samples .== sample_name)
    kmers = get_kmers(seq, idx.k)
    for (i, km) in enumerate(kmers)
        # i = 1; km = kmers[i]
        @show i
        # Check if the k-mer is already in the index; if not, add it and initialize its annotation
        if !haskey(idx.kmer_to_id, km.data)
            push!(idx.annotations, BitVector(zeros(length(idx.samples))))
            idx.kmer_to_id[km.data] = length(idx.annotations)
            # Topology: track successor for RowDiff graph traversal
            next_kmer_data = (i < length(kmers)) ? kmers[i + 1].data : 0
            push!(idx.adjacency, next_kmer_data) 
        end
        # Update the annotation for the current k-mer and sample
        row_idx = idx.kmer_to_id[km.data]
        # Mark the presence of the k-mer in the sample by setting the corresponding bit
        idx.annotations[row_idx][sample_id] = true
    end
    # idx
end

"""
    build_traversal_order(idx::MetaGraphIndex)

Build a traversal order for k-mers based on graph adjacency.
Returns a vector of row indices in traversal order and a mapping of kmers to their traversal positions.
"""
function build_traversal_order(idx::MetaGraphIndex)
    # Determine the number of k-mer rows in the annotations matrix.
    # This is the size of the graph in terms of distinct k-mers stored.
    n_kmers = length(idx.annotations)
    # Track which k-mer rows have already been visited in the traversal.
    # We use a BitVector for compactness and constant-time membership checks.
    visited = BitVector(zeros(n_kmers))
    # The output order of row indices as we follow adjacency chains.
    traversal_order = Int64[]
    # Iterate over every possible row index so that disconnected chains
    # and isolated k-mers are also included in the final order.
    for start_id in 1:n_kmers
        # start_id = 2
        # If this row was already included in a previous chain, skip it.
        if visited[start_id]
            continue
        end
        # Begin traversing from an unvisited starting row.
        current_id = start_id
        while (current_id != 0) && !visited[current_id]
            # Append the current row to the traversal order.
            push!(traversal_order, current_id)

            # Mark the current row as visited to avoid revisiting it in later chains.
            visited[current_id] = true

            # Follow the adjacency link stored for this k-mer row.
            # The adjacency vector stores the next k-mer's raw hash value.
            next_kmer_data = idx.adjacency[current_id]

            # Convert the raw hash to the corresponding row index if it exists.
            # If the adjacency entry is zero or not present in the index, end the chain.
            current_id = if next_kmer_data != 0 && haskey(idx.kmer_to_id, next_kmer_data)
                idx.kmer_to_id[next_kmer_data]
            else
                0
            end
        end
    end

    # Return the final traversal order for subsequent compression.
    return traversal_order
end

struct CompressedMetaGraph
    kmer_to_id::Dict{UInt64, Int}
    # RowDiff transformed matrix: stores XOR deltas
    diff_matrix::Vector{BitVector}
    anchors::BitVector # Which rows are stored raw (sinks or random samplings)
    traversal_order::Vector{Int64} # Order of k-mers for traversal
end

"""
compress(idx::MetaGraphIndex)

Compress the MetaGraphIndex using RowDiff compression.
Stores XOR deltas between consecutive k-mers in the traversal order.
"""
function compress(idx::MetaGraphIndex)
    n_rows = length(idx.annotations)
    traversal_order = build_traversal_order(idx)
    diff_matrix = Vector{BitVector}(undef, n_rows)
    anchors = BitVector(zeros(n_rows))
    
    # Initialize all diff_matrix entries
    for i in 1:n_rows
        diff_matrix[i] = BitVector(zeros(length(idx.samples)))
    end
    
    # Apply RowDiff compression along the traversal order
    for i in 1:length(traversal_order)
        current_row = traversal_order[i]
        
        if i < length(traversal_order)
            # Store XOR delta between current and next in traversal order
            next_row = traversal_order[i + 1]
            diff_matrix[current_row] = idx.annotations[current_row] .⊻ idx.annotations[next_row]
        else
            # Last in chain: store raw value and mark as anchor
            diff_matrix[current_row] = idx.annotations[current_row]
            anchors[current_row] = true
        end
    end
    
    return CompressedMetaGraph(idx.kmer_to_id, diff_matrix, anchors, traversal_order)
end