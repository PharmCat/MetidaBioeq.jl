using Test, DataFrames, CSV, CategoricalArrays, StatsModels
using MetidaNCA

path     = dirname(@__FILE__)
io       = IOBuffer();

rdsdict    = Dict()
rdsdict[1] = rds1 = bedf2x2 = CSV.File(joinpath(path, "csv", "2x2rds1.csv")) |> DataFrame
transform!(bedf2x2, :Subj => categorical, renamecols = false)
transform!(bedf2x2, :Per => categorical, renamecols = false)
bedf2x2.logVar = log.(bedf2x2.Var)

bedf2x2x4 = CSV.File(joinpath(path, "csv", "2x2x4rds1.csv")) |> DataFrame 
transform!(bedf2x2x4, :Subject => categorical, renamecols = false)
transform!(bedf2x2x4, :Period => categorical, renamecols = false)

refvals = CSV.File(joinpath(path, "csv", "ciref.csv")) |> DataFrame

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
    @test MetidaBioeq.result(berBmet)[1, 1] == "Formulation: T - R"
    @test MetidaBioeq.result(berBmet)[1, 2] == "logVar"
    @test MetidaBioeq.result(berBmet)[1,:level] == 90.0
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
    subject = :Subj, 
    formulation = :Trt, 
    period = :Per,
    sequence = :Seq, 
    autoseq = true)
    bedf2x2seq = MetidaBioeq.makeseq(bedf2x2, 
    subject = :Subj, 
    formulation = :Trt, 
    period = :Per)
    @test bedf2x2.Seq == bedf2x2seq

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
    formulation = :Trt, 
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
    subject = :Subj, 
    formulation = :Trt, 
    period = :Per,
    sequence = :Seq, 
    design = "2x2",
    autoseq = true,
    logt = false)

    # 2x2
    # 
    beAglm  = MetidaBioeq.estimate(be2;  estimator = "glm", method = "A")
    @test MetidaBioeq.result(beAglm)[1,2] == "log(Var)"
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
    subject = :Subj, 
    formulation = :Trt, 
    reference = "T",
    period = :Per,
    sequence = :Seq)
    beres =  MetidaBioeq.estimate(be2;  estimator = "glm", method = "A")
    @test MetidaBioeq.result(beres)[1,1] == "Trt: R - T"
    @test_nowarn  MetidaBioeq.estimate(be2;  estimator = "mm", method = "B")
    @test_nowarn  MetidaBioeq.estimate(be2;  estimator = "met", method = "B")

end

rdsdict[2] = CSV.File(joinpath(path, "csv", "2x2rds2.csv")) |> DataFrame
rdsdict[3] = CSV.File(joinpath(path, "csv", "2x2rds3.csv")) |> DataFrame
rdsdict[4] = CSV.File(joinpath(path, "csv", "2x2rds4.csv")) |> DataFrame
rdsdict[5] = CSV.File(joinpath(path, "csv", "2x2rds5.csv")) |> DataFrame
rdsdict[6] = CSV.File(joinpath(path, "csv", "2x2rds6.csv")) |> DataFrame
rdsdict[7] = CSV.File(joinpath(path, "csv", "2x2rds7.csv")) |> DataFrame
rdsdict[8] = CSV.File(joinpath(path, "csv", "2x2rds8.csv")) |> DataFrame


for (k,v) in rdsdict
    transform!(v, :Subj => categorical, renamecols = false)
    transform!(v, :Per => categorical, renamecols = false)
    v.logVar = log.(v.Var)
end

# Schütz, H., Labes, D., & Fuglsang, A. (2014). Reference datasets for 2-treatment, 2-sequence, 2-period bioequivalence studies. The AAPS journal, 16(6), 1292–1297. https://doi.org/10.1208/s12248-014-9661-0

@testset "  Validation 2x2" begin
    cidf = refvals[1:8, 4:6]

    for (k, v) in rdsdict
        be = MetidaBioeq.bioequivalence(v, 
        vars = :Var, 
        subject = :Subj, 
        formulation = :Trt, 
        period = :Per,
        sequence = :Seq, 
        design = "2x2",
        autoseq = true,
        logt = false,
        info = false)
        beres  = MetidaBioeq.estimate(be;  estimator = "glm", method = "A")
        df     = MetidaBioeq.result(beres)

        @test isapprox(df.GMR[1], cidf[k, "GMR"], atol = 0.01)
        @test isapprox(df.LCI[1], cidf[k, "LCI"], atol = 0.01) 
        @test isapprox(df.UCI[1], cidf[k, "UCI"], atol = 0.01)
    end
    for (k, v) in rdsdict
        be = MetidaBioeq.bioequivalence(v, 
        vars = :logVar, 
        subject = :Subj, 
        formulation = :Trt, 
        period = :Per,
        sequence = :Seq, 
        design = "2x2",
        autoseq = true,
        logt = true,
        info = false)
        beres  = MetidaBioeq.estimate(be;  estimator = "glm", method = "A")
        df     = MetidaBioeq.result(beres)

        @test isapprox(df.GMR[1], cidf[k, "GMR"], atol = 0.01)
        @test isapprox(df.LCI[1], cidf[k, "LCI"], atol = 0.01) 
        @test isapprox(df.UCI[1], cidf[k, "UCI"], atol = 0.01)
    end
end


# Reference Datasets for Studies in a Replicate Design Intended for Average Bioequivalence with Expanding Limits

rdsdict2 = Dict()
for i = 1:30
    rdsdict2[i] = CSV.File(joinpath(path, "csv", "rds$i.csv")) |> DataFrame
    transform!(rdsdict2[i], :subject => categorical, renamecols = false)
    transform!(rdsdict2[i], :period => categorical, renamecols = false)
    rdsdict2[i].logVar = log.(rdsdict2[i].PK)
end

# Schütz, H., Labes, D., Tomashevskiy, M., la Parra, M. G., Shitova, A., & Fuglsang, A. (2020). Reference Datasets for Studies in a Replicate Design Intended for Average Bioequivalence with Expanding Limits. The AAPS journal, 22(2), 44. https://doi.org/10.1208/s12248-020-0427-6

@testset "  Validation replicate (model B)" begin
    cidf = refvals[61:90, :]

    for i = 1:30
        @testset "  RDS $i" begin
            be = MetidaBioeq.bioequivalence(rdsdict2[i], 
            vars = :PK, 
            subject = :subject, 
            formulation = :treatment, 
            period = :period,
            sequence = :sequence, 
            autoseq = false,
            seqcheck = false,
            dropcheck = false,
            logt = false,
            info = false)

            beres = MetidaBioeq.estimate(be;  estimator = "mm", method = "B")
            df     = MetidaBioeq.result(beres)

            @test isapprox(df.LCI[1], cidf[i, "LCI"], atol = 0.01) 
            @test isapprox(df.UCI[1], cidf[i, "UCI"], atol = 0.01)
            @test isapprox(df.DF[1],  cidf[i, "DF"], atol = 0.1)

            beres = MetidaBioeq.estimate(be;  estimator = "met", method = "B")
            df     = MetidaBioeq.result(beres)

            @test isapprox(df.LCI[1], cidf[i, "LCISPSS"], atol = 0.01) 
            @test isapprox(df.UCI[1], cidf[i, "UCISPSS"], atol = 0.01)
            if !isapprox(df.DF[1], cidf[i, "DFSPSS"], atol = 0.1)
                @info "DF for Metida method B = $(df.DF[1]), reference value SAS = $(cidf[i, "DF"]), reference value SPSS = $(cidf[i, "DFSPSS"]). CI check passed."
            end
            @test isapprox(beres.models[1].result.reml, cidf[i, "REMLSPSS"], atol = 0.01)

        end
    end

    @testset "  Already log-transformed" begin
        for i = 1:30
        
            be = MetidaBioeq.bioequivalence(rdsdict2[i], 
            vars = :logVar, 
            subject = :subject, 
            formulation = :treatment, 
            period = :period,
            sequence = :sequence, 
            autoseq = false,
            seqcheck = false,
            dropcheck = false,
            logt = true,
            info = false)

            beres = MetidaBioeq.estimate(be;  estimator = "mm", method = "B")
            df     = MetidaBioeq.result(beres)
            @test isapprox(df.LCI[1], cidf[i, "LCI"], atol = 0.01) 
            @test isapprox(df.UCI[1], cidf[i, "UCI"], atol = 0.01)
            beres = MetidaBioeq.estimate(be;  estimator = "met", method = "B")
            df     = MetidaBioeq.result(beres)
            @test isapprox(df.LCI[1], cidf[i, "LCISPSS"], atol = 0.01) 
            @test isapprox(df.UCI[1], cidf[i, "UCISPSS"], atol = 0.01)
        end
    end

end

spssinfmsg = Dict()
spssinfmsg["*"] = """
 - The final Hessian matrix is not positive definite although all convergence criteria
are satisfied. The MIXED procedure continues despite this warning. Validity of subsequent
results cannot be ascertained."""

spssinfmsg["**"] = """
 - Iteration was terminated but convergence has not been achieved. The MIXED procedure
continues despite this warning. Subsequent results produced are based on the last iteration.
Validity of the model fit is uncertain."""
spssinfmsg[missing] = "" 


@testset "  Validation replicate (model C)" begin
    cidf = refvals[91:120, :]

    for i = 1:29 # RDS 30 not estimable with SPSS 
        @testset "  RDS $i" begin
            be = MetidaBioeq.bioequivalence(rdsdict2[i], 
            vars = :PK, 
            subject = :subject, 
            formulation = :treatment, 
            period = :period,
            sequence = :sequence, 
            autoseq = false,
            seqcheck = false,
            dropcheck = false,
            logt = false,
            info = false)

            beres = MetidaBioeq.estimate(be;  estimator = "met", method = "C")
            df     = MetidaBioeq.result(beres)

            if !isapprox(df.LCI[1], cidf[i, "LCISPSS"], atol = 0.1)
                @info "RDS$i LCI for Metida method C = $(df.LCI[1]), reference value SPSS = $(cidf[i, "LCISPSS"]). "

            end
            if !isapprox(df.UCI[1], cidf[i, "UCISPSS"], atol = 0.1)
                @info "RDS$i UCI for Metida method C = $(df.UCI[1]), reference value SPSS = $(cidf[i, "UCISPSS"]). "

            end
  
            @test isapprox(beres.models[1].result.reml, cidf[i, "REMLSPSS"], atol = 0.01)

        end
    end
end

rdsdict3 = Dict()
for i = 1:11
    rdsdict3[i] = CSV.File(joinpath(path, "csv", "pds$i.txt")) |> DataFrame
    transform!(rdsdict3[i], :Subj => categorical, renamecols = false)
    rdsdict3[i].logVar = log.(rdsdict3[i].Var)
end

# Fuglsang, A., Schütz, H., & Labes, D. (2015). Reference datasets for bioequivalence trials in a two-group parallel design. The AAPS journal, 17(2), 400–404. https://doi.org/10.1208/s12248-014-9704-6

@testset "  Validation parallel, no welch correction" begin
    cidf = refvals[20:30, :]

    for i = 1:11
        @testset "  RDS $i" begin
            be = MetidaBioeq.bioequivalence(rdsdict3[i], 
            vars = :Var, 
            subject = :Subj, 
            formulation = :Treat, 
            autoseq = false,
            seqcheck = false,
            dropcheck = false,
            logt = false,
            info = true)

            beres = MetidaBioeq.estimate(be;  estimator = "glm")
            df     = MetidaBioeq.result(beres)
            
            @test isapprox(df.LCI[1], cidf[i, "LCI"], atol = 0.01) 
            @test isapprox(df.UCI[1], cidf[i, "UCI"], atol = 0.01)
        end
    end

    for i = 1:11
        @testset "  RDS $i log-transformed" begin
            be = MetidaBioeq.bioequivalence(rdsdict3[i], 
            vars = :logVar, 
            subject = :Subj, 
            formulation = :Treat, 
            autoseq = false,
            seqcheck = false,
            dropcheck = false,
            logt = true,
            info = true)

            beres = MetidaBioeq.estimate(be;  estimator = "glm")
            df     = MetidaBioeq.result(beres)

            @test isapprox(df.LCI[1], cidf[i, "LCI"], atol = 0.01) 
            @test isapprox(df.UCI[1], cidf[i, "UCI"], atol = 0.01)
        end
    end
end


@testset "  NCA parallel design" begin
    ncads  = CSV.File(joinpath(dirname(pathof(MetidaNCA)), "..", "test", "csv",  "pkdata2.csv")) |> DataFrame
    ncares = nca(ncads, :Time, :Concentration, [:Subject, :Formulation])
    ncadf = DataFrame(ncares)
    transform!(ncadf, :Subject => categorical, renamecols = false)
    be = MetidaBioeq.bioequivalence(ncadf, 
            vars = [:Cmax, :AUClast], 
            subject = :Subject, 
            formulation = :Formulation, 
            autoseq = false,
            seqcheck = false,
            dropcheck = false,
            logt = false,
            info = false)

            beres = MetidaBioeq.estimate(be;  estimator = "glm")
            df    = MetidaBioeq.result(beres)

            @test df.Parameter[1] == "Formulation: T - R"
            @test df.Metric[1] == "log(Cmax)"
            @test isapprox(df.DF[1], 8.0, atol = 0.01)
            @test isapprox(df.LCI[1], 77.98, atol = 0.01) 
            @test isapprox(df.UCI[1], 149.52, atol = 0.01)

            @test df.Parameter[2] == "Formulation: T - R"
            @test df.Metric[2] == "log(AUClast)"
            @test isapprox(df.DF[2], 8.0, atol = 0.01)
            @test isapprox(df.LCI[2], 90.41, atol = 0.01) 
            @test isapprox(df.UCI[2], 154.16, atol = 0.01)
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