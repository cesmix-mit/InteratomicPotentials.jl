pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add InteratomicPotentials to environment stack

using InteratomicPotentials
using Documenter
using DocumenterCitations
using Literate

DocMeta.setdocmeta!(
    InteratomicPotentials,
    :DocTestSetup,
    :(using InteratomicPotentials);
    recursive = true,
)

bib = CitationBibliography(joinpath(@__DIR__, "citation.bib"))

# Generate examples

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR = joinpath(@__DIR__, "src/generated")

examples = [
    "Lennard Jones Cluster" => "LJCluster/feature_exploration.jl"
]

for (_, example_path) in examples
    s = split(example_path, "/")
    sub_path, file_name = string(s[1:end-1]...), s[end]
    example_filepath = joinpath(EXAMPLES_DIR, example_path)
    Literate.markdown(example_filepath,
                      joinpath(OUTPUT_DIR, sub_path),
                      documenter = true)
end

examples = [title => joinpath("generated", replace(example_path, ".jl" => ".md"))
            for (title, example_path) in examples]

makedocs(
      root    =  joinpath(dirname(pathof(InteratomicPotentials)), "..", "docs"),
      source  = "src",
      build   = "build",
      clean   = true,
      doctest = true,
      modules = [InteratomicPotentials],
      repo    = "https://github.com/cesmix-mit/InteratomicPotentials.jl/blob/{commit}{path}#{line}",
      highlightsig = true,
      sitename = "InteratomicPotentials.jl",
      expandfirst = [],
      draft = false,
      pages = [
        "Home" => "index.md",
        "InteratomicPotentials Interface" => "interface.md",
        "Examples" => examples,
        "API Reference" => "api.md",
        "References" => "bibliography.md"
      ],
      format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://cesmix-mit.github.io/InteratomicPotentials.jl",
        assets = String[],
      ),
      plugins=[bib]
)

deploydocs(;
    repo = "github.com/cesmix-mit/InteratomicPotentials.jl",
    devbranch = "main",
    push_preview = true
)
