struct Bioequivalence{D}
    vars::AbstractVector
    data::D
    design::String
    dropout::Union{Bool, Nothing}
    subject
    period
    formulation
    sequence
    stage
    reference
    subjects::Vector
    periods::Vector
    formulations::Vector
    sequences::Vector
    logt::Bool
end

struct BEResults{D} 
    be::Bioequivalence{D}
    models::Vector
    df::Dict
    estimator::String
    method::String
end

function Base.show(io::IO, obj::Bioequivalence)
    println(io, "Bioequivalence object:")
    println(io, "Design: $(obj.design)")
    println(io, "Dropouts: $(obj.dropout)")
    println(io, "Subject: $(obj.subject) ($(length(obj.subjects)))")
    print(io, "Formulation: $(obj.formulation) ($(obj.formulations))")
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
    if !isnothing(obj.period) print(io, "\nPeriod: $(obj.period) ($(obj.periods))") end
    if !isnothing(obj.sequence) 
        print(io, "\nSequence: $(obj.sequence) ($(obj.sequences))") 
    end
    if !isnothing(obj.stage) print(io, "\nStage: $(obj.stage)") end
end

function Base.show(io::IO, obj::BEResults)
    println(io, "  Result:")
    println(io, obj.df[:result])
    print(io, "Estimator: $(obj.estimator) ($(obj.method))")
end