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

struct Bioequivalence
    vars
    data
    design
    dropout
    subject
    period
    formulation
    sequence
    reference
    subjects
    periods
    formulations
    sequences
end

struct BEResults
    models
    df
end

"""
    bioquivalence(data;
    vars = nothing,
    subject = :subject,
    period = :period,
    formulation = :formulation,
    sequence = :sequence,
    reference = nothing,
    design = nothing,
    io::IO = stdout,
    seqcheck = true,
    dropcheck = true)


"""
function bioequivalence(data;
    vars = nothing,
    subject = :subject,
    period = :period,
    formulation = :formulation,
    sequence = nothing,
    reference = nothing,
    design = nothing,
    io::IO = stdout,
    seqcheck = true,
    dropcheck = true,
    info = true,
    autoseq = false)

    if isa(vars, Symbol) vars = [vars] 
    elseif isa(vars, String) vars = [Symbol(vars)] 
    elseif isa(vars, Vector{<:String}) vars = Symbol.(vars) 
    else
        error("Not supported vars type")
    end

    if isa(design, Symbol) design = string(design) end
    if isa(design, String) && design != "parallel" design = uppercase(design) end

    dfnames = Symbol.(Tables.columnnames(data))

    fac = [subject, formulation]

    fac ⊆ dfnames || error("Subject or formulation column not found in dataframe!")

    nomissing(data, subject) || error("Subject column have missing data")

    nomissing(data, formulation) || error("Formulation column have missing data")

    obsnum = size(data, 1)

    subjects = unique(Tables.getcolumn(data, subject))

    subjnum = length(subjects)

    formulations = sort!(unique(Tables.getcolumn(data, formulation)))

    if isnothing(reference)
        @info "Reference formulation not specified. First used: \"$(first(formulations))\"."
        reference = first(formulations)
    else
        reference ∈ formulations || error("Reference formulation \"$(reference)\" not found in dataframe.")
    end

    dropout = nothing
    periods = nothing
    sequences = nothing

    if isnothing(period) && isnothing(sequence) && isnothing(design)
        length(subjects) == length(Tables.getcolumn(data, subject)) || error("Trial design seems parallel, but subjects not unique!")
        design = "parallel"
        println(io, "Parallel desigh used.")
    end


    if isnothing(design) || design != "parallel"

        !isnothing(period) || error("Trial design seems NOT parallel, but period is nothing")

        autoseq || !isnothing(sequence) || error("Trial design seems NOT parallel, but sequence is nothing, autoseq is `false`")

        period ∈ dfnames || error("Period not found in dataframe!")

        if !isnothing(sequence)
            sequence ∈ dfnames || error("Sequence not found in dataframe!")
        end

        periods = sort!(unique(Tables.getcolumn(data, period)))

        push!(fac, period)

        nomissing(data, period) || error("Period column have missing data")

        if autoseq || seqcheck
            subjdict = Dict()
            for p in periods
                for i = 1:obsnum
                    if Tables.getcolumn(data, period)[i] == p
                        subj = Tables.getcolumn(data, subject)[i]
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
            if info && autoseq @info "autoseq is true, but sequence defined - sequence column used" end
            sequences = unique(getcol(data, sequence))
            push!(fac, sequence)

            nomissing(data, sequence) || error("Sequence column have missing data")
        end


        if dropcheck
            if !isnothing(vars) && !nomissing(data, vars)
                dropout = true
                @info "Dropuot(s) found in dataframe!"
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
            design = Symbol("$(length(formulations))X$(length(sequences))X$(length(periods))")
            info && @info  "Seems design type is: $design"
        else
            spldes = split(design, "X")
            if length(spldes) != 3 && design != "2X2" error("Unknown design type. Use fXsXp format or \"2Х2\".") end
            if length(formulations) != parse(Int, spldes[1]) error("Design error: formulations count wrong!") end
            if length(sequences) != parse(Int, spldes[2]) error("Design error: sequences count wrong!") end
            if length(periods) != parse(Int, spldes[3]) error("Design error: periods count wrong!") end
            info && @info "Design type seems fine..."
        end
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
        reference,
        subjects,
        periods,
        formulations,
        sequences)
end

function result(be; estimator = "auto", method = "auto", supresswarn = false)

    length(be.formulations) > 2 &&  error("More than 2 formulations not supported yet")
    
    if method == "auto"
        if disgn in ("parallel", "2X2") 
            method == "A" 
        else
            method == "B"  
        end
    end
    # Define estimator: MixedModel / Metida / GLM
    if estimator == "auto" 
        if design == "parallel" 
            estimator = "glm"
        else
            # Check DS correctness, mixed modela used for incomplete data
            if  design == "2X2"
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
    if be.design == "parallel"

        if esimator != "glm" && supresswarn @warn("Design is parallel, but estimator not GLM, GLM will be used!") end 
        if method != "A" && supresswarn @warn("Method not A, for parallel design GLM (Method A) will be used!") end
        esimator = "glm"
        method = "A"
        models = [@eval @formula($i ~ $(be.formulation)) for i in be.vars]

    elseif be.design == "2X2"

        if method != "A" && supresswarn @warn("Method not A, Method A will be used!") end
        if estimator == "glm"
            models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence) + $(be.subject)) for i in be.vars]
        elseif estimator == "met"
            models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence)) for i in be.vars]
        elseif estimator == "mm"
            models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence) + (1|  $(be.subject) )) for i in be.vars]
        else
            error("Unknown estimator!")
        end
        
    else
        if estimator == "glm" && method == "B"
            supresswarn || @warn("Method B used, estimator changed to MixedModels.jl!")
            estimator == "mm"
        elseif estimator == "glm" && method == "C"
            supresswarn || @warn("Method C used, estimator changed to Metida.jl!")
            estimator == "met"
        elseif estimator == "mm" && method == "C"
            supresswarn || @warn("Method C used, estimator changed to Metida.jl!")
            estimator == "met"
        end

        if estimator == "glm"
            models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence) + $(be.subject)) for i in be.vars]
        elseif estimator == "met"
            models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence)) for i in be.vars]
        elseif estimator == "mm"
            models = [@eval @formula($i ~ $(be.formulation) + $(be.period) + $(be.sequence) + (1| $(be.subject) ))  for i in be.vars]
        else
            error("Unknown estimator!")
        end
    end
    
    # If GLM used 
    if estimator == "glm"

        results = [fit(LinearModel, m, be.data; contrasts = Dict(be.formulation => DummyCoding(base = be.reference))) for m in models]

        df = DataFrame(Parameter = [], PE = [], SE = [], DF = [], lnLCI = [], lnUCI = [], GMR = [], LCI = [], UCI = [])
            for i in results
                DF = dof_residual(i)
                CI = confint(i, 0.9)[2,:]
                PE = coef(i)[2]
                push!(df, (string(coefnames(i)[2], " - ", be.reference),
                    PE,
                    stderror(i)[2],
                    DF,
                    CI[1],
                    CI[2],
                    exp(PE),
                    exp(CI[1]),
                    exp(CI[2])))
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
        df = DataFrame(Parameter = [], PE = [], SE = [], DF = [], lnLCI = [], lnUCI = [], GMR = [], LCI = [], UCI = [])
        for i in results
            DF = dof_satter(i, 2)
            dist = TDist(DF)
            PE = coef(i)[2]
            SE = stderror(i)[2]
            lnLCI = PE - SE * quantile(dist, 0.95)
            lnUCI = PE + SE * quantile(dist, 0.95)
            push!(df, (string(coefnames(i)[2], " - ", be.reference),
                PE,
                SE,
                DF,
                lnLCI,
                lnUCI,
                exp(PE),
                exp(lnLCI),
                exp(lnUCI)))
        end

    # If MixedModels Used
    elseif estimator == "mm"

        results = [fit(MixedModel, m, be.data; 
        contrasts = Dict(be.formulation => DummyCoding(base = be.reference)),
        REML=true
        ) for m in models]

        df = DataFrame(Parameter = [], PE = [], SE = [], DF = [], lnLCI = [], lnUCI = [], GMR = [], LCI = [], UCI = [])
        for i in results
            DF = nobs(i) - rank(hcat(i.X, i.reterms[1]))
            dist = TDist(DF)
            PE = coef(i)[2]
            SE = stderror(i)[2]
            lnLCI = PE - SE * quantile(dist, 0.95)
            lnUCI = PE + SE * quantile(dist, 0.95)

            push!(df, (string(coefnames(i)[2], " - ", be.reference),
                PE,
                SE,
                DF,
                lnLCI,
                lnUCI,
                exp(PE),
                exp(lnLCI),
                exp(lnUCI)))
        end


    end

    BEResults(results, df)
    
end