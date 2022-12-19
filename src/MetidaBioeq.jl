module MetidaBioeq

using Metida, MetidaBase, MetidaNCA, GLM, MixedModels, Distributions, Tables, 
DataFrames, CategoricalArrays

using MetidaBase.StatsModels

import Base: show

include("types.jl")
include("bioequivalence.jl")

end 
