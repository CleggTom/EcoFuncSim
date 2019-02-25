module EcoFuncSim

using OrdinaryDiffEq, DiffEqCallbacks, Distributions, LinearAlgebra

include("Boltzmann.jl")
include("Community.jl")
include("Parameters.jl")
include("DiffEqs.jl")
include("Simulate.jl")

export Boltz_params,boltz
export Community,rand_matrix,mod_matrix,growth_vectors
export make_params, simulate,mean_simulate

end # module
