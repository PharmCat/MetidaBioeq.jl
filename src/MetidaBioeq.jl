module MetidaBioeq

using Metida, MetidaBase, MetidaNCA, GLM, MixedModels,
DataFrames, CategoricalArrays


using MetidaBase.StatsModels, MetidaBase.Distributions, MetidaBase.Tables

import Base: show

include("types.jl")
include("bioequivalence.jl")

end 
