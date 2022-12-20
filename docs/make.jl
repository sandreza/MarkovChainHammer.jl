using Documenter
using MarkovChainHammer

api_dir = "API/"
dt_dir = "Discrete Time/"
ct_dir = "Continuous Time/"

api = Any[
    "Overview"=> api_dir * "overview.md",
    "Discrete Time"=> (api_dir * dt_dir) .* ["transfer_operators.md", "empirical_transfer_operators.md", "convergence.md"],
    "Continuous Time" => (api_dir * ct_dir) .* ["generators.md", "holding_times.md", "empirical_generator.md",  "convergence.md"],
]


mod_dir = "Modules/"
modules = Any[
    "Overview"=> mod_dir * "module_overview.md",
    "TransitionMatrix" => mod_dir * "transition_matrix.md",
    "Trajectory" => mod_dir * "trajectory.md",
    "Clustering" => mod_dir * "clustering.md",
    "Utils" => mod_dir * "utils.md",
]

makedocs(
    sitename="Markov Chain Hammer",
    format=Documenter.HTML(collapselevel=1),
    pages=[
        "Home" => "index.md",
        "API" => api,
        "Modules" => modules,
        "Function Index" => "function_index.md",
    ],
    modules=[MarkovChainHammer]
)

deploydocs(repo="github.com/sandreza/MarkovChainHammer.jl.git")
