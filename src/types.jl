struct Bioequivalence
    vars::AbstractVector
    data
    design::String
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
    logt
end

struct BEResults
    models
    df
    estimator
    method
end

function Base.show(io::IO, obj::Bioequivalence)
    println(io, "Bioequivalence object:")
    println(io, "Design: $(obj.design)")
    println(io, "Dropouts: $(obj.dropout)")
    println(io, "Subject: $(obj.subject) ($(length(obj.subjects)))")
    print(io, "Formulation: $(obj.formulation)")
    if isnothing(obj.formulations)
        println(io, ", Ref: $(obj.reference)")
    else
        print(io, " ($(obj.formulations[1])")
        if length(obj.formulations) > 1
            for i = 2:length(obj.formulations)
                print(io, ", $(obj.formulations[i])")
            end
        end
        print(io, ") Ref: $(obj.reference)")
    end
    if !isnothing(obj.period) print(io, "\nPeriod: $(obj.period)") end
    if !isnothing(obj.sequence) print(io, "\nSequence: $(obj.sequence)") end
end

function Base.show(io::IO, obj::BEResults)
    println(io, obj.df)
    print(io, "Estimator: $(obj.estimator)")
end