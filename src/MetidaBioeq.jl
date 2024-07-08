module MetidaBioeq

using  MetidaNCA, Metida, GLM, MixedModels, DataFrames, CategoricalArrays, Distributions, StatsBase

import Base: show

export result, bioquivalence, estimate, makeseq

include("types.jl")
include("bioequivalence.jl")

end 
