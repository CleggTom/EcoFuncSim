#Code to define and generate the Parameters object that is used in the simuations
"""
    Parameters(Com::Community,n_sp::Int64)

Type containing the parameters for a community. Includes rates (`Com`) and the number of
species (`n_sp`).
"""
struct Parameters
    Com::Community
    n_sp::Int64
end

"""
    make_params(U::Vector{Boltz_params},R::Vector{Boltz_params},a::Array{Boltz_params,2},T::Float64;intra::Bool = true)

function that makes a `Parameter` object given arrays of `Boltz_params`.
"""
function make_params(U::Vector{Boltz_params},R::Vector{Boltz_params},
    a::Array{Boltz_params,2},T::Float64)
    #check array sizes are correct
    @assert length(U) == length(R) == size(a)[1] == size(a)[2] "Parameter arrays not the same length"

    u = boltz.(U,T)
    r = boltz.(R,T)
    a = boltz.(a,T)

    Com = Community(u,r,a)

    return(Parameters(Com,length(u)))
end

"""
    make_params(U::Vector{Float64},R::Vector{Float64},a::Array{Float64,2}; intra::Bool = true)

function that makes a `Parameter` object given arrays of `Float64`.
"""
function make_params(U::Vector{Float64},R::Vector{Float64},a::Array{Float64,2};
    intra::Bool = true)
    #check array sizes are correct
    @assert length(U) == length(R) == size(a)[1] == size(a)[2] "Parameter arrays not the same length"

    if !intra
        for sp = 1:length(U)
            a[sp,sp] = 1.0
        end
    end

    Com = Community(U,R,a)

    return(Parameters(Com,length(U)))
end
