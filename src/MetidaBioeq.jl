module MetidaBioeq

using Metida, MetidaBase, MetidaNCA, GLM, MixedModels, Tables, 
DataFrames, CategoricalArrays


using MetidaBase.StatsModels, MetidaBase.Distributions

import Base: show

include("types.jl")
include("bioequivalence.jl")

end 
