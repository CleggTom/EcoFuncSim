const k = 8.6173303e-5

"""
    Boltz_params(B0::Float64,EB::Float64,Tref::Float64)

Container type for the parameters of the Boltzmann Arrhenius equation.
"""
struct Boltz_params
    B0::Float64
    EB::Float64
    Tref::Float64
end

"""
    boltz(B::Boltz_params,T::Float64)

Function to calculate rate based on the Boltzmann parameters. Takes a `Boltz_params`
object to paramerise the equation.
"""
function boltz(B::Boltz_params,T::Float64)
    return( B.B0::Float64 * exp((-B.EB::Float64/k) * ((1/T) - (1/B.Tref::Float64))) )
end
