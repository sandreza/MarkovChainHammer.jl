using Documenter
using MarkovChainHammer

mch_methods = Any[
    "Overview"=>"mch_methods.md",
    "Basics"=> "basics.md",
]

module_overview = Any[
    "Overview"=>"module_overview.md",
    "TransitionMatrix" => "transition_matrix.md",
    "Trajectory" => "trajectory.md",
    "Clustering" => "clustering.md",
    "Utils" => "utils.md",
]

makedocs(
    sitename="Markov Chain Hammer",
    format=Documenter.HTML(collapselevel=1),
    pages=[
        "Home" => "index.md",
        "Markov Chains" => mch_methods,
        "Modules" => module_overview,
        "Function Index" => "function_index.md",
    ],
    modules=[MarkovChainHammer]
)

deploydocs(repo="github.com/sandreza/MarkovChainHammer.jl.git")
