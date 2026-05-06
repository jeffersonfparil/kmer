# Define a mapping from nucleotides to their 2-bit representations
const NUCLEOTIDES = Dict('A' => 0x0, 'C' => 0x1, 'G' => 0x2, 'T' => 0x3)

# Define a struct to represent a k-mer, which consists of the encoded data and the length of the k-mer
struct Kmer
    data::UInt64
    k::Int
end

function Base.show(io::IO, kmer::Kmer)
    data = string(kmer.data, base=2, pad=kmer.k * 2)
    print(io, "Kmer(data=$data, k=$(kmer.k))")
end

# Extract k-mers using a sliding window and bit-shifting
function get_kmers(seq::String, k::Int; verbose::Bool=false)
    # seq = "ATCTCGATCGATCGACTACG"; k=5; verbose=false
    # Initialize an array to hold the k-mers
    kmers = Kmer[]
    # Create a bitmask to keep only the last k nucleotides (2 bits per nucleotide)
    mask = (UInt64(1) << (2 * k)) - 1
    # Use a rolling hash approach to compute the k-mers efficiently
    current = UInt64(0)
    # Iterate through the sequence and update the current k-mer using bit-shifting
    for (i, char) in enumerate(seq)
        # i = 1; char = seq[i]
        # Shift the current k-mer left by 2 bits and add the new nucleotide, then apply the mask
        current = ((current << 2) | NUCLEOTIDES[char]) & mask
        verbose ? println("i=$i; kmer=$(Kmer(current, k))") : nothing
        # Once we have processed at least k nucleotides, we can start adding k-mers to the array
        if i >= k
            push!(kmers, Kmer(current, k))
        end
    end
    verbose ? println.(kmers) : nothing
    return kmers
end