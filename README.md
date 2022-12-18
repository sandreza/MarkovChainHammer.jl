<!-- Title -->
<h1 align="center">
  MarkovChainHammer.jl
</h1>

A toolkit for analyzing, generating, and constructing both continuous and discrete markov chains with a finite state space. See the [documentation](https://sandreza.github.io/MarkovChainHammer.jl/dev/) for the mathematical underpinnings, repository functionality, and examples.

 <a href="https://mit-license.org">
    <img alt="MIT license" src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square">
  </a>
 <a href="https://sandreza.github.io/MarkovChainHammer.jl/dev">
    <img alt="Documentation" src="https://img.shields.io/badge/documentation-stable%20release-red?style=flat-square">
  </a>
 <a href="https://github.com/SciML/ColPrac">
    <img alt="ColPrac: Contributor's Guide on Collaborative Practices for Community Packages" src="https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet?style=flat-square">

## Contents
* [Tools](#tools)
* [Installation instructions](#installation-instructions)
* [Contributing](#contributing)

## Tools

By default MarkovChainHammer exports no functions and has modules that
1. Construct [transfer operators](https://en.wikipedia.org/wiki/Transfer_operator) from data
2. Construct [Markov chains](https://en.wikipedia.org/wiki/Markov_chain) from transfer operators
3. Detects [communities](https://en.wikipedia.org/wiki/Community_structure) from transfer operators
4. Exports utilities for plotting with [Makie](https://github.com/MakieOrg/Makie.jl)

See the [examples](https://github.com/sandreza/MarkovChainHammer.jl/tree/main/examples) for inspiration on how the utilities can be used.
## Installation instructions

MarkovChainHammer is an ***unregistered*** Julia package that requires Julia 1.8+. To install it,

1. [Download Julia](https://julialang.org/downloads/).
2. Launch Julia and type 
```julia
julia> using Pkg

julia> Pkg.add("https://github.com/sandreza/MarkovChainHammer.jl.git")
```

## Contributing 

We follow [Julia conventions](https://docs.julialang.org/en/v1/manual/style-guide/) and recommend reading through [ColPrac](https://docs.sciml.ai/ColPrac/stable/) as a standard guide contributing to Julia Software. New issues and pull requests are welcome!


