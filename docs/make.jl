using Documenter
using MarkovChainHammer

mch_methods = Any[
    "Home"=>"mch_methods.md",
    "Markov Chains"=>"basics.md",
]

makedocs(
    sitename="Markov Chain Hammer",
    format=Documenter.HTML(collapselevel=1),
    pages=[
        "Home" => "index.md",
        "Markov Chains" => mch_methods,
        "Function Index" => "function_index.md",
    ],
    modules=[MarkovChainHammer]
)

deploydocs(repo="github.com/sandreza/MarkovChainHammer.jl.git")
