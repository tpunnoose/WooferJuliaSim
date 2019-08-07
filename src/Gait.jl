@with_kw struct GaitParams
	# Default == Trot

	# add in a cyclic array for the phases?
	num_phases::Int64 = 4

	# 4xnum_phase array of contacts for each gait phase
	contact_phases::Array{Int64} = [1 1 1 0; 1 0 1 1; 1 0 1 1; 1 1 1 0]

	phase_times::Vector{Float64} = [0.1, 0.3, 0.1, 0.3]

	phase_length::Float64 = sum(phase_times)
	alpha::Float64 = 0.5
end

function getPhase(t::AbstractFloat, gait_params::GaitParams)
	phase_time = t % gait_params.phase_length

	for i in 1:gait_params.num_phases
		if phase_time < sum(gait_params.phase_times[1:i])
			return i
		end
	end
end

function nextPhase(phase::Integer, gait_params::GaitParams)
	if (phase == gait_params.num_phases)
		return 1
	else
		return phase+1
	end
end