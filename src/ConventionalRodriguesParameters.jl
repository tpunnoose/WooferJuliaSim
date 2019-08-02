# for a description of CRPs and stereographic projections of quaternions
# see: Schaub, Hanspeter, and John L. Junkins. "Stereographic orientation
# 		parameters for attitude dynamics: A generalization of the Rodrigues
#		parameters." Journal of the Astronautical Sciences 44.1 (1996): 1-19.

function quaternion2CRP!(g::Vector{Float64}, q::Vector{Float64})
	# g is the Conventional Rodriuges Parameters
	g[1] = q[2]/q[1]
	g[2] = q[3]/q[1]
	g[3] = q[4]/q[1]
end

function CRP2quaternion!(q::Vector{Float64}, g::Vector{Float64})
	# g is the Conventional Rodriuges Parameters
	q[1] = 1/sqrt(1 + g'*g)
	q[2] = g[1]*q[1]
	q[3] = g[2]*q[1]
	q[4] = g[3]*q[1]
end

function CRP2DCM!(R::Array{Float64, 2}, g::Vector{Float64})
	R[1,1] = 1 + g[1]^2 - g[2]^2 - g[3]^2
	R[1,2] = 2*(g[1]*g[2] + g[3])
	R[1,3] = 2*(g[1]*g[3] - g[2])
	R[2,1] = 2*(g[1]*g[2] - g[3])
	R[2,2] = 1 - g[1]^2 + g[2]^2 - g[3]^2
	R[2,3] = 2*(g[2]*g[3] + g[1])
	R[3,1] = 2*(g[3]*g[1] + g[2])
	R[3,2] = 2*(g[3]*g[2] - g[1])
	R[3,3] = 1 - g[1]^2 - g[2]^2 + g[3]^2
end
