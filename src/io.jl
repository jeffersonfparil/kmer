# Define a mapping from nucleotides to their 2-bit representations
const NUCLEOTIDES = Dict('A' => 0x0, 'C' => 0x1, 'G' => 0x2, 'T' => 0x3)

"""
A struct representing a k-mer, a substring of length `k` from a DNA sequence.

The `data` field is a `UInt64` encoding the nucleotides using 2 bits per nucleotide:
- A: 00
- C: 01
- G: 10
- T: 11

The `k` field specifies the length of the k-mer (number of nucleotides).
"""
struct Kmer
    data::UInt64
    k::Int64
end

"""
    Base.show(io::IO, kmer::Kmer)

Custom display for `Kmer` instances.

This method formats the internal `UInt64` representation of the k-mer as a
binary string padded to `2 * k` bits and prints it in the form
`Kmer(data=..., k=...)`.

# Examples
```jldoctest
julia> x = Kmer(0b0001101100, 5);

julia> show(x)
Kmer(data=0001101100, k=5)

julia> print(x)
Kmer(data=0001101100, k=5)
```
"""
function Base.show(io::IO, kmer::Kmer)
    data = string(kmer.data, base=2, pad=kmer.k * 2)
    print(io, "Kmer(data=$data, k=$(kmer.k))")
end

"""
    get_kmers(seq::String, k::Int64; verbose::Bool=false)

Extract k-mers from a DNA sequence using a rolling bit-shift approach.

# Arguments
- `seq::String`: the DNA sequence containing characters `A`, `C`, `G`, and `T`.
- `k::Int64`: the length of each k-mer.
- `verbose::Bool=false`: if `true`, print intermediate k-mer construction steps.

# Returns
- `Vector{Kmer}`: the encoded k-mers for the input sequence.

# Examples
```jldoctest
julia> get_kmers("ATCTCGATCGATCGACTACG", 5)
16-element Vector{Kmer}:
 Kmer(data=0011011101, k=5)
 Kmer(data=1101110110, k=5)
 Kmer(data=0111011000, k=5)
 Kmer(data=1101100011, k=5)
 Kmer(data=0110001101, k=5)
 Kmer(data=1000110110, k=5)
 Kmer(data=0011011000, k=5)
 Kmer(data=1101100011, k=5)
 Kmer(data=0110001101, k=5)
 Kmer(data=1000110110, k=5)
 Kmer(data=0011011000, k=5)
 Kmer(data=1101100001, k=5)
 Kmer(data=0110000111, k=5)
 Kmer(data=1000011100, k=5)
 Kmer(data=0001110001, k=5)
 Kmer(data=0111000110, k=5)
```
"""
function get_kmers(seq::String, k::Int64; verbose::Bool=false)
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