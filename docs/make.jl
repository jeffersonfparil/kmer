using kmer
using Documenter

DocMeta.setdocmeta!(kmer, :DocTestSetup, :(using kmer); recursive=true)

makedocs(;
    modules=[kmer],
    authors="jeffersonparil@gmail.com",
    sitename="kmer.jl",
    format=Documenter.HTML(;
        canonical="https://jeffersonfparil.github.io/kmer.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jeffersonfparil/kmer.jl",
    devbranch="main",
)
