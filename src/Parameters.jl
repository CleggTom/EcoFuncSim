"""
    Community(U::Vector{Float64},R::Vector{Float64},a::Array{Float64,2})

Container type for rate parameters in the simulations. Note that these are the actual
rates (as opposed to their temperature dependancys) that have been calculated from the
parameters already.
"""
struct Community
    U::Vector{Float64}
    R::Vector{Float64}
    a::Array{Float64,2}
end

"""
    Parameters(Com::Community,n_sp::Int64)

Type containing the parameters for a community. Includes rates (`Com`) and the number of
species (`n_sp`).
"""
struct Parameters
    Com::Community
    n_sp::Int64
end


function make_params(U::Vector{Boltz_params},R::Vector{Boltz_params},
    a::Array{Boltz_params,2},T::Float64;intra = true)
    #check array sizes are correct
    @assert length(U) == length(R) == size(a)[1] == size(a)[2] "Parameter arrays not the same length"

    u = boltz.(U,T)
    r = boltz.(R,T)
    a = boltz.(a,T)

    if !intra
        for sp = 1:length(u)
            a[sp,sp] = 1.0
        end
    end

    Com = Community(u,r,a)

    return(Parameters(Com,length(u)))
end

function make_params(U::Vector{Float64},R::Vector{Float64},a::Array{Float64,2})
    #check array sizes are correct
    @assert length(U) == length(R) == size(a)[1] == size(a)[2] "Parameter arrays not the same length"

    Com = Community(U,R,a)

    return(Parameters(Com,length(U)))
end
