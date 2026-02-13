# MetidaBioeq
# Copyright © 2022-2025 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# SPDX-License-Identifier: MIT
module MetidaBioeq

    using  Metida, GLM, MixedModels, DataFrames, CategoricalArrays, Distributions, StatsBase

    import Base: show

    import MetidaBase: cvfromsd, cvfromvar, sdfromcv

    export result, bioquivalence, estimate, makeseq, cvfromsd, sdfromcv

    include("types.jl")
    include("bioequivalence.jl")

end 
