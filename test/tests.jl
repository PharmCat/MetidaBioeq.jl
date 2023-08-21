using Test, DataFrames, CSV, CategoricalArrays, StatsModels

path     = dirname(@__FILE__)
io       = IOBuffer();

rds1 = bedf2x2 = CSV.File(joinpath(path, "csv", "2x2rds1.csv")) |> DataFrame
transform!(bedf2x2, :Subject => categorical, renamecols = false)
transform!(bedf2x2, :Period => categorical, renamecols = false)
bedf2x2.logVar = log.(bedf2x2.Var)

bedf2x2x4 = CSV.File(joinpath(path, "csv", "2x2x4rds1.csv")) |> DataFrame 
transform!(bedf2x2x4, :Subject => categorical, renamecols = false)
transform!(bedf2x2x4, :Period => categorical, renamecols = false)

@testset "  Basic test" begin

    @test MetidaBioeq.nomissing(bedf2x2x4, :logVar) == false
    @test MetidaBioeq.nomissing(bedf2x2x4, [:logVar]) == false
    @test MetidaBioeq.nomissing(bedf2x2, :logVar) == true
    @test MetidaBioeq.nomissing(bedf2x2, [:logVar]) == true

    # replicate(repeated) design 
    # 2x2x4
    be = MetidaBioeq.bioequivalence(bedf2x2x4, 
    vars = :logVar, 
    subject = :Subject, 
    formulation = :Formulation, 
    reference = "R",
    period = :Period,
    sequence = :Sequence, 
    autoseq = true)
    berBmet  = MetidaBioeq.estimate(be;  estimator = "met", method = "B")
    @test berBmet.df[1,1] == "Formulation: T - R"
    @test berBmet.df[1,2] == "logVar"
    @test berBmet.df[1,:level] == 90.0
    berCmet  = MetidaBioeq.estimate(be;  estimator = "met", method = "C")
    berBmm   = MetidaBioeq.estimate(be;  estimator = "mm", method = "B")
    berAglm  = MetidaBioeq.estimate(be;  estimator = "glm", method = "A")
    @test berBmet.method  == "B"
    @test berBmet.estimator== "met"
    @test berCmet.method == "C"
    @test berCmet.estimator == "met"
    @test berBmm.method == "B"
    @test berBmm.estimator== "mm"
    @test berAglm.method == "A"
    @test berAglm.estimator== "glm"

    @test_nowarn MetidaBioeq.result(berCmet)
    # Inappropriate 
    # Try "glm" && "B" "C"
    beres = MetidaBioeq.estimate(be;  estimator = "glm", method = "B")
    @test beres.method  == "B"
    @test beres.estimator== "mm"
    beres = MetidaBioeq.estimate(be;  estimator = "glm", method = "C")
    @test beres.method  == "C"
    @test beres.estimator== "met"
    beres = MetidaBioeq.estimate(be;  estimator = "mm", method = "C")
    @test beres.method  == "C"
    @test beres.estimator== "met"
    beres = MetidaBioeq.estimate(be;  estimator = "met", method = "P")
    @test beres.method  == "B"
    @test beres.estimator== "met"
    beres = MetidaBioeq.estimate(be;  estimator = "mm", method = "P")
    @test beres.method  == "B"
    @test beres.estimator== "mm"



    # Crossover design 2X2
    #
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

    beAglm  = MetidaBioeq.estimate(be2;  estimator = "glm", method = "A")
    beBmm   = MetidaBioeq.estimate(be2;  estimator = "mm", method = "B")
    beBmet   = MetidaBioeq.estimate(be2;  estimator = "met", method = "B")
    @test beAglm.method  == "A"
    @test beAglm.estimator == "glm"
    @test beBmm.method  == "B"
    @test beBmm.estimator== "mm"
    @test beBmet.method == "B"
    @test beBmet.estimator== "met"

    # PARALLEL DESIGN 
    #
    bedf2x2.psbj = collect(1:size(bedf2x2, 1))
    transform!(bedf2x2, :psbj => categorical, renamecols = false)
    be3 = MetidaBioeq.bioequivalence(bedf2x2, 
    vars = :logVar, 
    subject = :psbj, 
    formulation = :Formulation, 
    autoseq = true)

    berPglm  = MetidaBioeq.estimate(be3;  estimator = "glm")
    @test berPglm.method    == "P" 
    @test berPglm.estimator == "glm"
    
    # Inappropriate
    # Try mm and met
    berPglm  = MetidaBioeq.estimate(be3;  estimator = "glm", method = "B")
    @test berPglm.method    == "P" 
    @test berPglm.estimator == "glm"
    berPglm  = MetidaBioeq.estimate(be3;  estimator = "mm")
    @test berPglm.method    == "P" 
    @test berPglm.estimator == "glm"
    berPglm  = MetidaBioeq.estimate(be3;  estimator = "met")
    @test berPglm.method    == "P" 
    @test berPglm.estimator == "glm"

    # 2x2 not LOG-TRANSFORMED var (logt = false)
    #
    be2 = MetidaBioeq.bioequivalence(bedf2x2, 
    vars = :Var, 
    subject = :Subject, 
    formulation = :Formulation, 
    period = :Period,
    sequence = :Sequence, 
    design = "2x2",
    autoseq = true,
    logt = false)

    # 2x2
    # 
    beAglm  = MetidaBioeq.estimate(be2;  estimator = "glm", method = "A")
    @test beAglm.df[1,2] == "log(Var)"
    beBmm   = MetidaBioeq.estimate(be2;  estimator = "mm", method = "B")
    beBmet  = MetidaBioeq.estimate(be2;  estimator = "met", method = "B")
    @test beAglm.method == "A"
    @test beAglm.estimator == "glm"
    @test beBmm.method == "B"
    @test beBmm.estimator == "mm"
    @test beBmet.method == "B"
    @test beBmet.estimator == "met"

    # Inappropriate
    # Try C
    beBmet   = MetidaBioeq.estimate(be2;  estimator = "met", method = "C")
    @test beBmet.method == "B"
    @test beBmet.estimator == "met"
    beAglm   = MetidaBioeq.estimate(be2;  estimator = "glm", method = "P")
    @test beAglm.method == "A"
    @test beAglm.estimator == "glm"

    # More than 1 var; reverce reference
    be2 = MetidaBioeq.bioequivalence(bedf2x2, 
    vars = [:Var, :logVar], 
    subject = :Subject, 
    formulation = :Formulation, 
    reference = "T",
    period = :Period,
    sequence = :Sequence)
    beres =  MetidaBioeq.estimate(be2;  estimator = "glm", method = "A")
    @test beres.df[1,1] == "Formulation: R - T"
    @test_nowarn  MetidaBioeq.estimate(be2;  estimator = "mm", method = "B")
    @test_nowarn  MetidaBioeq.estimate(be2;  estimator = "met", method = "B")

end

rds2 = CSV.File(joinpath(path, "csv", "2x2rds2.csv")) |> DataFrame
transform!(rds2, :Subj => categorical, renamecols = false)
transform!(rds2, :Per => categorical, renamecols = false)
rds2.logVar = log.(rds2.Var)

rds3 = CSV.File(joinpath(path, "csv", "2x2rds3.csv")) |> DataFrame
transform!(rds3, :Subj => categorical, renamecols = false)
transform!(rds3, :Per => categorical, renamecols = false)
rds3.logVar = log.(rds3.Var)

rds4 = CSV.File(joinpath(path, "csv", "2x2rds4.csv")) |> DataFrame
transform!(rds4, :Subj => categorical, renamecols = false)
transform!(rds4, :Per => categorical, renamecols = false)
rds4.logVar = log.(rds4.Var)

rds5 = CSV.File(joinpath(path, "csv", "2x2rds5.csv")) |> DataFrame
transform!(rds5, :Subj => categorical, renamecols = false)
transform!(rds5, :Per => categorical, renamecols = false)
rds5.logVar = log.(rds5.Var)

rds6 = CSV.File(joinpath(path, "csv", "2x2rds6.csv")) |> DataFrame
transform!(rds6, :Subj => categorical, renamecols = false)
transform!(rds6, :Per => categorical, renamecols = false)
rds6.logVar = log.(rds6.Var)

rds7 = CSV.File(joinpath(path, "csv", "2x2rds7.csv")) |> DataFrame
transform!(rds7, :Subj => categorical, renamecols = false)
transform!(rds7, :Per => categorical, renamecols = false)
rds7.logVar = log.(rds7.Var)

rds8 = CSV.File(joinpath(path, "csv", "2x2rds8.csv")) |> DataFrame
transform!(rds8, :Subj => categorical, renamecols = false)
transform!(rds8, :Per => categorical, renamecols = false)
rds8.logVar = log.(rds8.Var)

@testset "  Validation" begin
    # RDS 1
    be = MetidaBioeq.bioequivalence(rds1, 
    vars = :Var, 
    subject = :Subject, 
    formulation = :Formulation, 
    period = :Period,
    sequence = :Sequence, 
    design = "2x2",
    autoseq = true,
    logt = false)
    beres  = MetidaBioeq.estimate(be;  estimator = "glm", method = "A")
    df     = MetidaBioeq.result(beres)

    @test isapprox(df.GMR[1], 95.09, atol = 0.01)
    @test isapprox(df.LCI[1], 90.76, atol = 0.01) 
    @test isapprox(df.UCI[1], 99.62, atol = 0.01)

    # RDS 2
    be = MetidaBioeq.bioequivalence(rds2, 
    vars = :Var, 
    subject = :Subj, 
    formulation = :Trt, 
    period = :Per,
    sequence = :Seq, 
    design = "2x2",
    autoseq = true,
    logt = false)
    beres  = MetidaBioeq.estimate(be;  estimator = "glm", method = "A")
    df     = MetidaBioeq.result(beres)

    @test isapprox(df.GMR[1], 71.10, atol = 0.01)
    @test isapprox(df.LCI[1], 51.45, atol = 0.01) 
    @test isapprox(df.UCI[1], 98.26, atol = 0.01)

    # RDS 3
    be = MetidaBioeq.bioequivalence(rds3, 
    vars = :Var, 
    subject = :Subj, 
    formulation = :Trt, 
    period = :Per,
    sequence = :Seq, 
    design = "2x2",
    autoseq = true,
    logt = false)
    beres  = MetidaBioeq.estimate(be;  estimator = "glm", method = "A")
    df     = MetidaBioeq.result(beres)

    @test isapprox(df.GMR[1], 58.56, atol = 0.01)
    @test isapprox(df.LCI[1], 39.41, atol = 0.01) 
    @test isapprox(df.UCI[1], 87.03, atol = 0.01)

    # RDS 4
    be = MetidaBioeq.bioequivalence(rds4, 
    vars = :Var, 
    subject = :Subj, 
    formulation = :Trt, 
    period = :Per,
    sequence = :Seq, 
    design = "2x2",
    autoseq = true,
    logt = false)
    beres  = MetidaBioeq.estimate(be;  estimator = "glm", method = "A")
    df     = MetidaBioeq.result(beres)

    @test isapprox(df.GMR[1], 71.10, atol = 0.01)
    @test isapprox(df.LCI[1], 51.45, atol = 0.01) 
    @test isapprox(df.UCI[1], 98.26, atol = 0.01)


    # RDS 5
    be = MetidaBioeq.bioequivalence(rds5, 
    vars = :Var, 
    subject = :Subj, 
    formulation = :Trt, 
    period = :Per,
    sequence = :Seq, 
    design = "2x2",
    autoseq = true,
    logt = false)
    beres  = MetidaBioeq.estimate(be;  estimator = "glm", method = "A")
    df     = MetidaBioeq.result(beres)

    @test isapprox(df.GMR[1], 91.83, atol = 0.01)
    @test isapprox(df.LCI[1], 55.71, atol = 0.01) 
    @test isapprox(df.UCI[1], 151.37, atol = 0.01)

    # RDS 6
    be = MetidaBioeq.bioequivalence(rds6, 
    vars = :Var, 
    subject = :Subj, 
    formulation = :Trt, 
    period = :Per,
    sequence = :Seq, 
    design = "2x2",
    autoseq = true,
    logt = false)
    beres  = MetidaBioeq.estimate(be;  estimator = "glm", method = "A")
    df     = MetidaBioeq.result(beres)

    @test isapprox(df.GMR[1], 99.89, atol = 0.01)
    @test isapprox(df.LCI[1], 93.37, atol = 0.01) 
    @test isapprox(df.UCI[1], 106.86, atol = 0.01)

    # RDS 7
    be = MetidaBioeq.bioequivalence(rds7, 
    vars = :Var, 
    subject = :Subj, 
    formulation = :Trt, 
    period = :Per,
    sequence = :Seq, 
    design = "2x2",
    autoseq = true,
    logt = false)
    beres  = MetidaBioeq.estimate(be;  estimator = "glm", method = "A")
    df     = MetidaBioeq.result(beres)

    @test isapprox(df.GMR[1], 92.15, atol = 0.01)
    @test isapprox(df.LCI[1], 88.46, atol = 0.01) 
    @test isapprox(df.UCI[1], 95.99, atol = 0.01)

    # RDS 8
    be = MetidaBioeq.bioequivalence(rds8, 
    vars = :Var, 
    subject = :Subj, 
    formulation = :Trt, 
    period = :Per,
    sequence = :Seq, 
    design = "2x2",
    autoseq = true,
    logt = false)
    beres  = MetidaBioeq.estimate(be;  estimator = "glm", method = "A")
    df     = MetidaBioeq.result(beres)

    @test isapprox(df.GMR[1], 93.42, atol = 0.01)
    @test isapprox(df.LCI[1], 86.81, atol = 0.01) 
    @test isapprox(df.UCI[1], 100.55, atol = 0.01)

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