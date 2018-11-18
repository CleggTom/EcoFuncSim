module EcoFuncSim

using OrdinaryDiffEq, DiffEqCallbacks, StaticArrays

include("Boltzmann.jl")
include("Parameters.jl")
include("DiffEqs.jl")
include("Simulate.jl")

export Boltz_params, boltz, Community, make_params, simulate

end # module
