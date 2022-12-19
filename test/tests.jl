using Test
using DataFrames, CSV, CategoricalArrays, StatsModels

path     = dirname(@__FILE__)
io       = IOBuffer();

bedf2x2 = CSV.File(joinpath(path, "csv", "2x2rds1.csv")) |> DataFrame
transform!(bedf2x2, :Subject => categorical, renamecols = false)
transform!(bedf2x2, :Period => categorical, renamecols = false)
bedf2x2.logVar = log.(bedf2x2.Var)

bedf2x2x4 = CSV.File(joinpath(path, "csv", "2x2x4rds1.csv")) |> DataFrame 
transform!(bedf2x2x4, :Subject => categorical, renamecols = false)
transform!(bedf2x2x4, :Period => categorical, renamecols = false)

@testset "  Test" begin

    @test MetidaBioeq.nomissing(bedf2x2x4, :logVar) == false
    @test MetidaBioeq.nomissing(bedf2x2x4, [:logVar]) == false
    @test MetidaBioeq.nomissing(bedf2x2, :logVar) == true
    @test MetidaBioeq.nomissing(bedf2x2, [:logVar]) == true

    be = MetidaBioeq.bioequivalence(bedf2x2x4, 
    vars = :logVar, 
    subject = :Subject, 
    formulation = :Formulation, 
    period = :Period,
    sequence = :Sequence, 
    autoseq = true)

    berBmet  = MetidaBioeq.result(be;  estimator = "met", method = "B")
    berCmet  = MetidaBioeq.result(be;  estimator = "met", method = "C")
    berBmm   = MetidaBioeq.result(be;  estimator = "mm", method = "B")
    berAglm  = MetidaBioeq.result(be;  estimator = "glm", method = "A")


    be2 = MetidaBioeq.bioequivalence(bedf2x2, 
    vars = :logVar, 
    subject = :Subject, 
    formulation = :Formulation, 
    period = :Period,
    sequence = :Sequence, 
    autoseq = true)

    bedf2x2seq = MetidaBioeq.makeseq(bedf2x2, 
    subject = :Subject, 
    formulation = :Formulation, 
    period = :Period)

    @test bedf2x2.Sequence == bedf2x2seq

    beAglm  = MetidaBioeq.result(be2;  estimator = "glm", method = "A")
    beBmm   = MetidaBioeq.result(be2;  estimator = "mm", method = "B")
    beBmm   = MetidaBioeq.result(be2;  estimator = "met", method = "B")

    bedf2x2.psbj = collect(1:size(bedf2x2, 1))

    be3 = MetidaBioeq.bioequivalence(bedf2x2, 
    vars = :logVar, 
    subject = :psbj, 
    formulation = :Formulation, 
    autoseq = true)

    berPmm  = MetidaBioeq.result(be3;  estimator = "glm")
end


# Examples
#=
using MixedModels
mm = fit(MixedModel, 
@formula(Var ~ Formulation + Period + Sequence + (1| Subject)), 
bedf2x2x4;
REML=true
)
coef(mm)
stderror(mm)[2]
nobs(mm) - rank(hcat(mm.X, mm.reterms[1]))

using GLM
ols = fit(LinearModel, 
@formula(Var ~ Formulation + Period + Sequence +  Subject), 
bedf2x2)

nobs(ols) - dof(ols) + 1
nobs(ols) - rank(ols.mm.m)
=#

#=
# Deprecated
m1 = @eval Metida.@covstr(1| $(be.subject))
m2 = Metida.@covstr(1|subject)
lmm = Metida.LMM(@formula(logPK~sequence+period+treatment), ds;
    random = Metida.VarEffect(m1, Metida.SI),
    )
lmm = Metida.fit(Metida.LMM, @formula(logPK~sequence+period+treatment), ds;
    random = Metida.VarEffect(m1, Metida.SI),
    )
=#