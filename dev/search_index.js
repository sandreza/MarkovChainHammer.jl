var documenterSearchIndex = {"docs":
[{"location":"mch_methods/#Markov-Chains","page":"Home","title":"Markov Chains","text":"","category":"section"},{"location":"mch_methods/","page":"Home","title":"Home","text":"This repository gathers various convenience tools for analyzing and generating Markov Chains with a finite state space. ","category":"page"},{"location":"mch_methods/","page":"Home","title":"Home","text":"The following section contains a review of Markov Chains","category":"page"},{"location":"mch_methods/","page":"Home","title":"Home","text":"Basics","category":"page"},{"location":"basics/#sec:basics","page":"Markov Chains","title":"Markov Chains","text":"","category":"section"},{"location":"basics/","page":"Markov Chains","title":"Markov Chains","text":"under construction","category":"page"},{"location":"function_index/#List-of-functions-in-MarkovChainHammer","page":"Function Index","title":"List of functions in MarkovChainHammer","text":"","category":"section"},{"location":"function_index/","page":"Function Index","title":"Function Index","text":"Modules = [ MarkovChainHammer.TransitionMatrix, MarkovChainHammer.Trajectory, MarkovChainHammer.Utils]","category":"page"},{"location":"function_index/#MarkovChainHammer.TransitionMatrix.generator-Tuple{Any}","page":"Function Index","title":"MarkovChainHammer.TransitionMatrix.generator","text":"generator(markov_chain; dt=1)\n\nDescription\n\nCalculate the generator matrix from a markov chain.\n\nArguments\n\nmarkov_chain::AbstractVector: A vector of integers representing the state of a markov chain at each time step.\ndt::Real: The time step between each state.\n\nReturns\n\ngenerator_matrix::Matrix: The generator matrix of the markov chain.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#MarkovChainHammer.TransitionMatrix.holding_times-Tuple{Any, Any}","page":"Function Index","title":"MarkovChainHammer.TransitionMatrix.holding_times","text":"holding_times(markov_chain, number_of_states; dt=1)\n\nDescription\n\nCalculate the holding times of a markov chain.\n\nArguments\n\nmarkov_chain::AbstractVector: A vector of integers representing the state of a markov chain at each time step.\nnumber_of_states::Integer: The number of states in the markov chain.\ndt::Real: The time step of the markov chain.\n\nReturns\n\nholding_times::Vector{Vector{Real}}: A vector of vectors of holding times for each state.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#MarkovChainHammer.TransitionMatrix.perron_frobenius-Tuple{Any}","page":"Function Index","title":"MarkovChainHammer.TransitionMatrix.perron_frobenius","text":"perron_frobenius(markov_chain)\n\nDescription\n\nCalculate the perron-frobenius matrix from a markov chain.\n\nArguments\n\nmarkov_chain::AbstractVector: A vector of integers representing the state of a markov chain at each time step.\n\nReturns\n\nperron_frobenius_matrix::Matrix: The perron-frobenius matrix of the markov chain.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#MarkovChainHammer.TransitionMatrix.symmetric_generator-Tuple{Any, Any}","page":"Function Index","title":"MarkovChainHammer.TransitionMatrix.symmetric_generator","text":"symmetric_generator(markov_chain, symmetries; dt=1)\n\nDescription\n\nCalculate the generator matrix from a markov chain with symmetries.\n\nArguments\n\nmarkov_chain::AbstractVector: A vector of integers representing the state of a markov chain at each time step.\nsymmetries::AbstractVector: A vector of functions that are symmetries of the markov chain.\ndt::Real: The time step between each state.\n\nReturns\n\ngenerator_matrix::Matrix: The generator matrix of the markov chain.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#MarkovChainHammer.TransitionMatrix.symmetric_perron_frobenius-Tuple{Any, Any}","page":"Function Index","title":"MarkovChainHammer.TransitionMatrix.symmetric_perron_frobenius","text":"symmetric_perron_frobenius(markov_chain, symmetries)\n\nDescription\n\nCalculate the perron-frobenius matrix from a markov chain with symmetries.\n\nArguments\n\nmarkov_chain::AbstractVector: A vector of integers representing the state of a markov chain at each time step.\nsymmetries::AbstractVector: A vector of functions that are symmetries of the markov chain.\n\nReturns\n\nperron_frobenius_matrix::Matrix: The perron-frobenius matrix of the markov chain.\n\n\n\n\n\n","category":"method"},{"location":"function_index/#MarkovChainHammer.Utils.histogram-Tuple{Any}","page":"Function Index","title":"MarkovChainHammer.Utils.histogram","text":"function histogram(     array;     bins=minimum([100, length(array)]),     normalization=:uniform,     custom_range=false )\n\nDescription\n\nUtility function for barplot in GLMakie\n\nArguments\n\narray; one dimensional sequence of numbers \n\nKeyword Arguments\n\nbins; how many buckets to use, always uniform\nnormalization; how much to weight each value of array (default 1/length(array))\ncustom_range; range for uniform bucket\n\n\n\n\n\n","category":"method"},{"location":"#MarkovChainHammer.jl","page":"Home","title":"MarkovChainHammer.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for MarkovChainHammer.jl","category":"page"}]
}
