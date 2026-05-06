# TODO: Test, add more comments and documentation

function reconstruct_row(c_idx::CompressedMetaGraph, row_id::Int)
    current_row = copy(c_idx.diff_matrix[row_id])
    curr = row_id
    
    # Cascade through the deltas until we hit an anchor
    while !c_idx.anchors[curr] && curr < length(c_idx.diff_matrix)
        curr += 1
        current_row .⊻= c_idx.diff_matrix[curr]
    end
    return current_row
end

function search(c_idx::CompressedMetaGraph, query::String, k::Int)
    query_kmers = get_kmers(query, k)
    if isempty(query_kmers) return BitVector() end
    
    # Start with all bits set (Universal set)
    result = nothing
    
    for km in query_kmers
        if haskey(c_idx.kmer_to_id, km.data)
            row_id = c_idx.kmer_to_id[km.data]
            row = reconstruct_row(c_idx, row_id)
            
            if isnothing(result)
                result = row
            else
                result .&= row # Intersection: sequence must exist in the sample
            end
        else
            return BitVector() # K-mer missing, whole sequence missing
        end
    end
    return result
end

# Initialize with 3 samples
index = MetaGraphIndex(31, ["Human_Gut", "Soil_Alpha", "Ocean_Surface"])

# Index some mock data
add_sample!(index, "Human_Gut", "ATGC...") 
add_sample!(index, "Ocean_Surface", "ATGC...")

# Compress and Search
compressed = compress(index)
hits = search(compressed, "ATGC", 31)

println("Found in samples: ", findall(hits))