module MetidaBioeq

using Metida, MetidaBase, MetidaNCA, GLM, MixedModels, Distributions, Tables, 
DataFrames, CategoricalArrays

import Base: show

include("types.jl")
include("bioequivalence.jl")

end 
