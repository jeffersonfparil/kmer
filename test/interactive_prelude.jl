using Pkg
Pkg.activate(".")
try
    Pkg.update()
catch
    nothing
end
using kmer

# using PkgTemplates
# t = Template(;
#     user="jeffersonfparil",
#     authors=["jeffersonparil@gmail.com"],
#     dir="./",
#     julia=v"1.12",
#     plugins=[
#         License(; name="GPL-3.0+", path=nothing, destination="LICENSE.md"),
#         CompatHelper(),
#         GitHubActions(;
#         osx=false,
#         windows=false,
#         ),
#         Documenter{GitHubActions}(),
#         Git(;
#             ignore=[
#                 "*.code-workspace",
#                 "*.jl.*.cov",
#                 "*.jl.cov",
#                 "*.jl.mem",
#                 ".DS_Store",
#                 "/docs/Manifest.toml",
#                 "/docs/build/",
#                 "Manifest.toml",
#                 "docs/build/",
#                 "tmp/",
#                 "*.svg",
#                 "*.jld2",
#                 "*.tsv",
#                 "*.csv",
#                 "*.txt"
#             ],
#             ssh=true
#         ),
#     ],
# )
# t("kmer")