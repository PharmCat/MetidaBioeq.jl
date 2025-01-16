module MetidaBioeq

    using  MetidaNCA, Metida, GLM, MixedModels, DataFrames, CategoricalArrays, Distributions, StatsBase

    import Base: show

    import MetidaBase: cvfromsd, cvfromvar, sdfromcv

    export result, bioquivalence, estimate, makeseq, cvfromsd, sdfromcv

    include("types.jl")
    include("bioequivalence.jl")

end 
