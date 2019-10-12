using Parameters


@with_kw struct SafetyParams
    maximumjointvelocity::Float64 = 40.0  # [rad/s]
    maximumjointangle::Float64 = π/4  # [rad]
    maximumCOMdeviation::Float64 = 0.2 # [m]
    nominalCOMposition::Vector{Float64} = [0.0, 0.0, 0.34]  # [m]
end

@with_kw mutable struct SafetyStatus
    safe::Bool = true
end

function zeroouttorques!(jointtorques::Vector)
    jointtorques .= zeros(size(jointtorques))
end

function safejointangles(jointpos::Vector, safetyparams::SafetyParams)
    return all(abs.(jointpos) .< safetyparams.maximumjointangle)
end

function safejointvelocities(jointvel::Vector, safetyparams::SafetyParams)
    return all(abs.(jointvel) .< safetyparams.maximumjointvelocity)
end

function safeCOMposition(x::Vector, safetyparams::SafetyParams)
    return all(abs.(x - safetyparams.nominalCOMposition) .< safetyparams.maximumCOMdeviation)
end

function issafe(x, jointpos::Vector{T}, jointvel::Vector{T}, safetyparams::SafetyParams) where {T<:Number}
    if ~safejointangles(jointpos, safetyparams)
        println("Aborting: A joint exceeded the maximum allowable angle of $(safetyparams.maximumjointangle) rads")
        return false
    end
    if ~safejointvelocities(jointvel, safetyparams)
        println("Aborting: a joint velocity exceeded the maximum allowable velocity of $(safetyparams.maximumjointvelocity) rad/s")
        return false
    end
    if ~safeCOMposition(x, safetyparams)
        println("Aborting: COM position deviated beyond safe threshold of $(safetyparams.maximumCOMdeviation)")
        return false
    end
    return true
end