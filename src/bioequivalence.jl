#

getcol = Tables.getcolumn

function nomissing(data, col)
    c = getcol(data, col)
    !any(ismissing, c)
end

function nomissing(data, cols::AbstractVector)
    for col in cols
        if !nomissing(data, col)  return false end
    end
    true
end

function functional_term(f, arg_expr...)
	expr = Expr(:call, Symbol(f), arg_expr...)
	eval(:(@formula 0 ~ $expr)).rhs
end

"""
    bioquivalence(data;
    vars = nothing,
    subject = :subject,
    period = :period,
    formulation = :formulation,
    sequence = :sequence,
    stage = nothing,
    reference = nothing,
    design = nothing,
    io::IO = stdout,
    seqcheck = true,
    dropcheck = true,
    logt = true)

* `vars` - variabel's column(s);
* `subject` - subject's column;
* `period` - period's column;
* `formulation` - formulation's column;
* `sequence` -sequence's column;
* `stage` - stage's column;
* `reference` - reference value for `formulation` column;
* `design` - design: "parallel", "2X2", "2X2X2", "2X2X4", ets. (formulations X sequences X periods);
* `seqcheck` - check sequencs;     
* `dropcheck` - dropuot check;
* `logt` - if `true` (default) data is already log-transformed, else `log()` will be used.
"""
function bioequivalence(data;
    vars = nothing,
    subject::Union{String, Symbol},
    period::Union{String, Symbol, Nothing} = nothing,
    formulation::Union{String, Symbol},
    sequence::Union{String, Symbol, Nothing} = nothing,
    stage::Union{String, Symbol, Nothing} = nothing,
    reference::Union{String, Symbol, Nothing} = nothing,
    design::Union{String, Symbol, Nothing} = nothing,
    io::IO = stdout,
    seqcheck::Bool = true,
    dropcheck::Bool = true,
    info::Bool = true,
    warns::Bool = true,
    autoseq::Bool = false,
    logt::Bool = true)

    if !isa(vars, Vector{<:Symbol}) 
        if isa(vars, Symbol) vars = [vars] 
        elseif isa(vars, String) vars = [Symbol(vars)] 
        elseif isa(vars, Vector{<:String}) vars = Symbol.(vars) 
        else
            error("Not supported vars type")
        end
    end

    if isa(design, Symbol) design = string(design) end
    if isa(design, String) && design != "parallel" design = uppercase(design) end

    if isa(subject, String) subject = Symbol(subject) end
    if isa(formulation, String) formulation = Symbol(formulation) end
    if isa(period, String) period = Symbol(period) end
    if isa(sequence, String) sequence = Symbol(sequence) end

    dfnames = Symbol.(Tables.columnnames(data))

    fac = [subject, formulation]

    fac ⊆ dfnames || error("Subject or formulation column not found in dataframe!")
    # Check subject and formulation column is categorical
    if isa(Tables.getcolumn(data, subject), AbstractVector{<:Real}) && !isa(Tables.getcolumn(data, subject), AbstractCategoricalArray)
        warns && @warn "Seems subject column is not categorical."
    end
    if isa(Tables.getcolumn(data, formulation), AbstractVector{<:Real}) && !isa(Tables.getcolumn(data, formulation), AbstractCategoricalArray)
        warns && @warn "Seems formulation column is not categorical."
    end

    # Subject column can't have missing data
    nomissing(data, subject) || error("Subject column have missing data")
    # Formulation column can't have missing data
    nomissing(data, formulation) || error("Formulation column have missing data")

    # Unique subjects and formulations
    subjects     = unique(getcol(data, subject))
    formulations = sort!(unique(getcol(data, formulation)))

    subjnum = length(subjects)
    obsnum  = size(data, 1)
    # If reference not defined - first level used as base
    if isnothing(reference)
        info && @info "Reference formulation not specified. First used: \"$(first(formulations))\"."
        reference = first(formulations)
    else
        reference ∈ formulations || error("Reference formulation \"$(reference)\" not found in dataframe.")
    end

    dropout   = nothing
    periods   = nothing
    sequences = nothing

    # For parallel design period and sequence are nothing
    if isnothing(period) && isnothing(sequence) && isnothing(design)
        subjnum == length(Tables.getcolumn(data, subject)) || error("Trial design seems parallel, but subjects not unique!")
        design = "parallel"
        info && @info "Parallel desigh used."
    end

    # check if design is not parallel
    if isnothing(design) || design != "parallel"
        # Period should be defined
        !isnothing(period) || error("Trial design seems NOT parallel, but period is nothing")
        # Sequence should be defined
        autoseq || !isnothing(sequence) || error("Trial design seems NOT parallel, but sequence is `nothing`, autoseq is `false`")
        # 
        period ∈ dfnames || error("Period not found in dataframe!")
        # Check period column is categorical
        if isa(Tables.getcolumn(data, period), AbstractVector{<:Real}) && !isa(Tables.getcolumn(data, period), AbstractCategoricalArray)
            @warn "Seems period column is not categorical."
        end
        # If sequence defined it should be in table
        if !isnothing(sequence)
            sequence ∈ dfnames || error("Sequence not found in dataframe!")
        else
            # Check sequence column is categorical
            if isa(Tables.getcolumn(data, sequence), AbstractVector{<:Real}) && !isa(Tables.getcolumn(data, sequence), AbstractCategoricalArray)
                @warn "Seems sequence column is not categorical."
            end
        end

        periods = sort!(unique(Tables.getcolumn(data, period)))
        push!(fac, period)
        # Period column can't have missing data
        nomissing(data, period) || error("Period column have missing data")
        # Check sequences
        if autoseq || seqcheck
            subjdict = Dict()
            for p in periods
                for i = 1:obsnum
                    if getcol(data, period)[i] == p
                        subj = getcol(data, subject)[i]
                        if haskey(subjdict, subj)
                            subjdict[subj] *= string(getcol(data, formulation)[i])
                        else
                            subjdict[subj] = string(getcol(data, formulation)[i])
                        end
                    end
                end
            end
        end

        if isnothing(sequence) && autoseq
            sequences = unique(values(subjdict))
        elseif isnothing(sequence)
            error("Sequence is nothing, but autoseq is false")
        else
            info && autoseq && @info "autoseq is `true`, but sequence defined - sequence column used"
            sequences = unique(getcol(data, sequence))
            nomissing(data, sequence) || error("Sequence column have missing data")
            push!(fac, sequence)
        end


        if dropcheck
            if !isnothing(vars) && !nomissing(data, vars)
                info && @info "Dropuot(s) found in dataframe!"
                dropout = true
            elseif !isnothing(vars)
                info && @info "No dropuot(s) found in dataframe!"
                dropout = false
            end
        end

        if seqcheck && !isnothing(sequence)
            for i = 1:obsnum
                if getcol(data, sequence)[i] != subjdict[getcol(data, subject)[i]]
                    error("Sequence error or data is incomplete! \n Subject: $(getcol(data, subject)[i]), Sequence: $(getcol(data, sequence)[i]), auto: $(subjdict[getcol(data, subject)[i]])")
                end
            end
            if length(unique(length.(sequences))) > 1
                error("Some sequence have different length!")
            end
            info && @info "Sequences looks correct..."
        end

        if isnothing(design)
            info && @info "Trying to find out the design..."
            design = "$(length(formulations))X$(length(sequences))X$(length(periods))"
            info && @info  "Seems design type is: $design"
        else
            if design == "2X2" design = "2X2X2" end
            spldes = split(design, "X")
            if length(spldes) != 3 error("Unknown design type. Use fXsXp format or \"2Х2\".") end
            if length(formulations) != parse(Int, spldes[1]) error("Design error: formulations count wrong!") end
            if length(sequences) != parse(Int, spldes[2]) error("Design error: sequences count wrong!") end
            if length(periods) != parse(Int, spldes[3]) error("Design error: periods count wrong!") end
            info && @info "Design type seems fine..."
        end
    end

    if !isnothing(stage)
        stage ⊆ dfnames || error("Stage column not found in dataframe!")
        if !(design in ["parallel", "2X2", "2X2X2"]) @warn "Stage is defined but design not equal \"parallel\", \"2X2\" or \"2X2X2\"." end
    end

    Bioequivalence(
        vars,
        data,
        design,
        dropout,
        subject,
        period,
        formulation,
        sequence,
        stage,
        reference,
        subjects,
        periods,
        formulations,
        sequences,
        logt)
end

"""
    makeseq(data;
        subject = :subject,
        period = :period,
        formulation = :formulation)

Make sequence vector from `data` and `subject`, `period`, `formulation` columns.
"""
function makeseq(data;
    subject = :subject,
    period = :period,
    formulation = :formulation)

    dfnames = Symbol.(Tables.columnnames(data))

    subject ∈ dfnames     || error("Subject column not found in dataframe!")
    period ∈ dfnames      || error("Period column not found in dataframe!")
    formulation ∈ dfnames || error("Formulation column not found in dataframe!")
    # Subject column can't have missing data
    nomissing(data, subject)     || error("Subject column have missing data")
    # Formulation column can't have missing data
    nomissing(data, formulation) || error("Formulation column have missing data")
    # Period column can't have missing data
    nomissing(data, period)      || error("Period column have missing data")

    # Unique periods
    periods      = sort!(unique(getcol(data, period)))

    obsnum  = size(data, 1)

    subjdict = Dict()
    for p in periods
        for i = 1:obsnum
            if getcol(data, period)[i] == p
                subj = getcol(data, subject)[i]
                if haskey(subjdict, subj)
                    subjdict[subj] *= string(getcol(data, formulation)[i])
                else
                    subjdict[subj] = string(getcol(data, formulation)[i])
                end
            end
        end
    end
    map(x-> subjdict[x], getcol(data, subject))
end

"""
    result(be; estimator = "auto", method = "auto", supresswarn = false)

`method` - Model settings.

* if `method == "auto"` than method `A` used for "2X2" and "2X2X2"  designes, method `P` for "parallel" design and method `B` for any other.

*Methods:*

* `A` using GLM and model `@formula(var ~ formulation + period + sequence + subject)`
* `B` using MixedModels and model `@formula(var ~ formulation + period + sequence + (1|subject))` or Metida and model `@lmmformula(v ~ formulation + period + sequence, random = 1|subject:SI)`
* `C` using Metida and model `@lmmformula(v ~ formulation + period + sequence, random = formulation|subject:CSH, repeated = formulation|subject:DIAG)`
* `P` using GLM and model `@formula(var ~ formulation)`

`estimator` - Estimator settings.

* if `estimator == "auto"` than GLM used for "parallel" design; for "2X2" design used GLM if no droputs and MixedModels if `dropout == true`; for other designes with method `C` Metida used and MixedModel for other cases.

*Estimators:*

* "glm" for GLM (https://juliastats.org/GLM.jl/stable/)
* "mm" for MixedModels (https://juliastats.org/MixedModels.jl/stable/)
* "met" for Metida (https://pharmcat.github.io/Metida.jl/stable/)

*Other autosettings:*

If design is "parallel" `estimator` set as "glm" and `method` as "P".


If design is "2X2" and method is "P" or "C" than if `estimator` == "glm" method set as "A" and "B" for other estimators. 


If design not "parallel" or "2X2": 

if method not "A", "B" or "C" than set as "A" for "glm" ann as B for other estimators;

if `estimator` == "glm" and `method` == "B" than `estimator` set as "mm", if `estimator` == "glm" or "mm" and `method` == "C" than `estimator` set as "met".


Reference:

EMA: [GUIDELINE ON THE INVESTIGATION OF BIOEQUIVALENCE](https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-bioequivalence-rev1_en.pdf)

EMA: [GUIDELINE ON THE INVESTIGATION OF BIOEQUIVALENCE, Annex I](https://www.ema.europa.eu/en/documents/other/31-annex-i-statistical-analysis-methods-compatible-ema-bioequivalence-guideline_en.pdf)

"""
function result(be; estimator = "auto", method = "auto", supresswarn = false, alpha = 0.05)

    length(be.formulations) > 2 &&  error("More than 2 formulations not supported yet")
    design = be.design
    if method == "auto"
        if design in ("2X2", "2X2X2") 
            method = "A" 
        elseif design == "parallel" 
            method = "P"
        else
            method = "B"  
        end
    end
    # Define estimator: MixedModel / Metida / GLM
    if estimator == "auto" 
        if design == "parallel" 
            estimator = "glm"
        else
            # Check DS correctness, mixed modela used for incomplete data
            if  design in ("2X2", "2X2X2") 
                if be.dropout 
                    estimator = "mm"
                else
                    estimator = "glm"
                end
            else
                if method == "C"
                    estimator = "met"
                else
                    estimator = "mm"
                end
            end

        end
    end
    
    # MODEL SELECTION
    if design == "parallel"

        if estimator != "glm" && !supresswarn @warn("Design is parallel, but estimator not GLM, GLM will be used!") end 
        if method != "P" && !supresswarn @warn("Method not P (parallel), for parallel simple GLM model will be used!") end
        estimator = "glm"
        method = "P"
        if be.logt
            models = [@eval @formula($i ~ $(be.formulation)) for i in be.vars]
        else
            models = [@eval @formula(log(Term($i)) ~ $(be.formulation)) for i in be.vars]
        end

    elseif design in ("2X2", "2X2X2") 

        if method in ("P", "C") 
            method == "C" && !supresswarn && @warn("Method C can't be used with 2X2 design!")
            method == "P" && !supresswarn && @warn("Method for parallel design can't be used with 2X2 design!")
            if estimator == "glm"
                !supresswarn && @warn("Method A will be used!")
                method = "A"
            else
                !supresswarn && @warn("Method B will be used!")
                method = "B"
            end
        end
  
        if estimator == "glm"
            if be.logt
                models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence) + $(be.subject)) for i in be.vars]
            else
                models = [begin 
                    rfo =  @eval  @formula(0 ~ $(be.formulation) + $(be.period) + $(be.sequence) + $(be.subject)) 
                    lhs = functional_term(log, i)
                    FormulaTerm(lhs, rfo.rhs) 
                end for i in be.vars]
            end
        elseif estimator == "met"
            if be.logt
                models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence)) for i in be.vars]
            else
                models = [begin 
                rfo = @eval @formula(0 ~ $(be.formulation) + $(be.period) + $(be.sequence))
                lhs = functional_term(log, i)
                FormulaTerm(lhs, rfo.rhs) 
            end for i in be.vars]
            end
        elseif estimator == "mm"
            if be.logt
                models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence) + (1|  $(be.subject) )) for i in be.vars]
            else
                models = [begin 
                rfo = @eval @formula(0 ~ $(be.formulation) + $(be.period) + $(be.sequence) + (1| $(be.subject) ))
                lhs = functional_term(log, i)
                FormulaTerm(lhs, rfo.rhs) 
                end for i in be.vars]
            end
        else
            error("Unknown estimator!")
        end
        
    else
        if !(method in ("A", "B", "C"))
            if estimator == "glm" 
                !supresswarn && @warn("Method P or unknown, method changed to \"A\"!")
                method = "A"
            else
                !supresswarn && @warn("Method P or unknown, method changed to \"B\"!")
                method = "B"
            end
        end 
        if estimator == "glm" && method == "B"
            !supresswarn && @warn("Method B used, estimator changed to MixedModels.jl!")
            estimator = "mm"
        elseif estimator == "glm" && method == "C"
            !supresswarn && @warn("Method C used, estimator changed to Metida.jl!")
            estimator = "met"
        elseif estimator == "mm" && method == "C"
            !supresswarn && @warn("Method C used, estimator changed to Metida.jl!")
            estimator = "met"
        end

        if estimator == "glm"
            if be.logt
                models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence) + $(be.subject)) for i in be.vars]
            else
                models = [begin
                rfo = @eval @formula(0 ~ $(be.formulation) + $(be.period) + $(be.sequence) + $(be.subject)) 
                lhs = functional_term(log, i)
                FormulaTerm(lhs, rfo.rhs)
            end for i in be.vars]
            end
        elseif estimator == "met"
            if be.logt
                models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence)) for i in be.vars]
            else
                models = [begin 
                rfo = @eval @formula(0 ~ $(be.formulation) + $(be.period) + $(be.sequence))
                lhs = functional_term(log, i)
                FormulaTerm(lhs, rfo.rhs) 
                end for i in be.vars]
            end
        elseif estimator == "mm"
            if be.logt
                models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence) + (1| $(be.subject) ))  for i in be.vars]
            else
                models = [begin 
                rfo = @eval @formula(0 ~ $(be.formulation) + $(be.period) + $(be.sequence) + (1| $(be.subject) ))
                lhs = functional_term(log, i)
                #lhs = term(i)
                FormulaTerm(lhs, rfo.rhs) 
                end for i in be.vars]
                
            end
        else
            error("Unknown estimator!")
        end
    end
    
    # If GLM used 
    if estimator == "glm"

        results = [fit(LinearModel, m, be.data; contrasts = Dict(be.formulation => DummyCoding(base = be.reference))) for m in models]

        df = DataFrame(Parameter = String[], Metric = String[], PE = Float64[], SE = Float64[], DF = Float64[], lnLCI = Float64[], lnUCI = Float64[], GMR = Float64[], LCI = Float64[], UCI = Float64[], level = Float64[])
            for i in results
                DF = dof_residual(i)
                CI = confint(i, 1-2alpha)[2,:]
                PE = coef(i)[2]
                push!(df, (string(coefnames(i)[2], " - ", be.reference),
                    coefnames(i.mf.f.lhs),
                    PE,
                    stderror(i)[2],
                    DF,
                    CI[1],
                    CI[2],
                    exp(PE),
                    exp(CI[1]),
                    exp(CI[2]),
                    (1-2alpha)*100
                    ))
            end

    # If Metida Used
    elseif estimator == "met"

        if method == "B"
           
            results = [fit!(LMM(m, be.data;
            random = Metida.VarEffect(@eval(Metida.@covstr(1| $(be.subject))), Metida.SI), 
            contrasts = Dict(be.formulation => DummyCoding(base = be.reference)))) for m in models]
            
        elseif method == "C"
    
            results = [fit!(LMM(m, be.data; 
            random = Metida.VarEffect(@eval(Metida.@covstr($(be.formulation)|$(be.subject))), Metida.CSH),
            repeated = Metida.VarEffect(@eval(Metida.@covstr($(be.formulation)|$(be.subject))), Metida.DIAG),
            contrasts = Dict(be.formulation => DummyCoding(base = be.reference)))) for m in models]
            
        else
            error("Method A used or unknown method!")
        end
        # Take resulst from models
        df = DataFrame(Parameter = String[], Metric = String[], PE = Float64[], SE = Float64[], DF = Float64[], lnLCI = Float64[], lnUCI = Float64[], GMR = Float64[], LCI = Float64[], UCI = Float64[], level = Float64[])
        for i in results
            DF = dof_satter(i, 2)
            dist = TDist(DF)
            PE = coef(i)[2]
            SE = stderror(i)[2]
            lnLCI = PE - SE * quantile(dist, 1-alpha)
            lnUCI = PE + SE * quantile(dist, 1-alpha)
            push!(df, (string(coefnames(i)[2], " - ", be.reference),
                coefnames(i.mf.f.lhs),
                PE,
                SE,
                DF,
                lnLCI,
                lnUCI,
                exp(PE),
                exp(lnLCI),
                exp(lnUCI),
                (1-2alpha)*100
                ))
        end

    # If MixedModels Used
    elseif estimator == "mm"

        results = [fit(MixedModel, m, be.data; 
        contrasts = Dict(be.formulation => DummyCoding(base = be.reference)),
        REML=true
        ) for m in models]

        df = DataFrame(Parameter = String[], Metric = String[], PE = Float64[], SE = Float64[], DF = Float64[], lnLCI = Float64[], lnUCI = Float64[], GMR = Float64[], LCI = Float64[], UCI = Float64[], level = Float64[])
        for i in results
            DF = nobs(i) - rank(hcat(i.X, i.reterms[1]))
            dist = TDist(DF)
            PE = coef(i)[2]
            SE = stderror(i)[2]
            lnLCI = PE - SE * quantile(dist, 1-alpha)
            lnUCI = PE + SE * quantile(dist, 1-alpha)

            push!(df, (string(coefnames(i)[2], " - ", be.reference),
                responsename(i),
                PE,
                SE,
                DF,
                lnLCI,
                lnUCI,
                exp(PE),
                exp(lnLCI),
                exp(lnUCI),
                (1-2alpha)*100
                ))
        end


    end

    BEResults(results, df, estimator, method)
    
end