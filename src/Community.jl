#Code for functions that generate communities
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
    key

an array containing the required parameters to fully describe a random community model.
`μX` refers to the mean of parameters X whilst `σX` refers to the standard deviation.
"""
const key = [:μER, :μEU, :μEa_ij, :μEa_ii,
             :σER, :σEU, :σEa_ij, :σEa_ii,
             :μR0, :μU0, :μa0_ij, :μa0_ii,
             :σR0, :σU0, :σa0_ij, :σa0_ii]

#generate growth parameters
"""
    rand_growth(p::Dict{Symbol,Float64},nSp::Int; Tref::Float64 = 285.0,rand_U0::Bool = true, rand_R0::Bool = true,rand_EU::Bool = true, rand_ER::Bool = true)

Generates the normalisation constants for uptake and respiration parameters. Returns an array of `Boltz_params` for both rates (i.e. `[[U],[R]]`).

# Arguments
- `p::Dict{Symbol,Float64}`: parameter dictionary with elements contained `key`
- `nSp::Int`:  number of species in the community.
- `Tref::Float64 = 285.0`: refernce temperature to use
- `rand_U0::Bool = true`: whether U0 values should be drawn randomly. If not they are set as `p[:μU0]`.
- `rand_R0::Bool = true`: whether R0 values should be drawn randomly. If not they are set as `p[:μR0]`.
- `rand_EU::Bool = true`: whether EU values should be drawn randomly. If not they are set as `p[:μEU]`.
- `rand_ER::Bool = true`: whether ER values should be drawn randomly. If not they are set as `p[:μER]`.


"""
function rand_growth(p::Dict{Symbol,Float64},nSp::Int; Tref::Float64 = 285.0,
                        rand_U0::Bool = true, rand_R0::Bool = true,
                        rand_EU::Bool = true, rand_ER::Bool = true)

    @assert all(in.(key,Ref(keys(p)))) "Make sure all parameters are provided: \n $key"

    U0 = rand_U0 ? rand(TruncatedNormal(p[:μU0],p[:σU0],0,Inf),nSp) : fill(p[:μU0],nSp)
    R0 = rand_R0 ? rand(TruncatedNormal(p[:μR0],p[:σR0],0,Inf),nSp) : fill(p[:μR0],nSp)

    EU = rand_EU ? rand(TruncatedNormal(p[:μEU],p[:σEU],0,Inf),nSp) : fill(p[:μEU],nSp)
    ER = rand_ER ? rand(TruncatedNormal(p[:μER],p[:σER],0,Inf),nSp) : fill(p[:μER],nSp)

    U = Boltz_params.(U0,EU,Tref)
    R = Boltz_params.(R0,ER,Tref)

    return(hcat([U,R]...))
end

##interaction matrix generation
"""
    rand_matrix(p::Dict{Symbol,Float64},nSp::Int,C::Float64; Tref::Float64 = 285.0,rand_a0_ij::Bool = true, rand_a0_ii::Bool = true,rand_Ea_ij::Bool = true, rand_Ea_ii::Bool = true,pNeg::Float64 = 0.0, sym::Bool = true,intra::Bool = true)

Function that generates a random interaction matrix with interaction weights. This function returns a NxN matrix containing the `a0` values.

# Arguments
- `p::Dict{Symbol,Float64}`: parameter dictionary with elements contained `key`.
- `nSp::Int`: number of species in the community.
- `C::Float64`: Connectance of the interaction matrix.
- `Tref::Float64 = 285.0`: Reference temperature to use
- `rand_a0_ij::Bool = true`: whether interspecific interaction constants should be randomly sampled. If not they are set as `p[:μa0_ij]`.
- `rand_a0_ii::Bool = true`: whether intraspecific interaction constants should be randomly sampled. If not they are set as `p[:μa0_ii]`. Only used if `intra == false`.
- `rand_Ea_ij::Bool = true`: whether interspecific E values should be random. If not they are set as `p[:μEa_ij]`.
- `rand_Ea_ii::Bool = true`: whether intraspecific E values should be random. If not they are set as `p[:μEa_ii]`.
- `pNeg::Float64 = 0.0`: proportion of links which are assigned a negative value
- `sym::Bool = true`: whether the matrix should be symmetric
- `intra::Bool = true`: whether intraspecific interactions should be generated to ensure stability. This is done by ensuring diagonal domminance of the matrix.
"""
function rand_matrix(p::Dict{Symbol,Float64},nSp::Int,C::Float64; Tref::Float64 = 285.0,
                     rand_a0_ij::Bool = true, rand_a0_ii::Bool = true,
                     rand_Ea_ij::Bool = true, rand_Ea_ii::Bool = true,
                     pNeg::Float64 = 0.0, sym::Bool = true,
                     intra::Bool = true)

    @assert all(in.(key,Ref(keys(p)))) "Make sure all parameters are provided: \n $key"

    #make a0
    a0 = rand_a0_ij ? rand(TruncatedNormal(p[:μa0_ij],p[:σa0_ij],0,Inf),nSp,nSp) : fill(p[:μa0_ij],nSp,nSp)


    #set negative interactions
    neg_index = sample(1:nSp^2,Int(floor((nSp^2) * pNeg)),replace=false)
    a0[neg_index] .= -a0[neg_index]

    #set connectance
    C_index = sample(1:nSp^2,Int(floor( (1-C)*nSp^2 ) ),replace = false)
    a0[C_index] .= 0.0

    if sym
        a0 = Array(Symmetric(a0))   #set symmetry
    end

    #set intra
    if intra
        a0[diagind(a0)] .= 0.0
        a0[diagind(a0)] .= sum(abs.(a0),dims=2)[:,1]
    else
        a0[diagind(a0)] .= rand_a0_ii ? rand(TruncatedNormal(p[:μa0_ii],p[:σa0_ii],0,Inf),nSp) : p[:μa0_ii]
    end

    #get all Ea
    Ea = rand_Ea_ij ? rand(TruncatedNormal(p[:μEa_ij],p[:σEa_ij],0,Inf),nSp,nSp) : fill(p[:μEa_ij],nSp,nSp)
    #set connectance
    Ea[C_index] .= 0.0

    if sym
        Ea = Array(Symmetric(Ea))   #set symmetry
    end
    #set Ea_ii
    Ea[diagind(Ea)] .= rand_Ea_ii ? rand(TruncatedNormal(p[:μEa_ii],p[:σEa_ii],0,Inf),nSp) : p[:μEa_ii]

    return(Boltz_params.(a0,Ea,Tref))
end

"""
    function mod_matrix(p::Dict{Symbol,Float64},nSp::Int,C::Float64;rand_ij::Bool = true, rand_ii::Bool = true,pNeg::Float64 = 0.0, sym::Bool = true, intra::Bool = true, nGrp::Int = 1, pMod::Float64 = 1.0)

Function that generates a modular interaction matrix with interaction weights. This function returns a NxN matrix containing the `a0` values.

# Arguments
- `p::Dict{Symbol,Float64}`: parameter dictionary with elements contained `key`.
- `nSp::Int`: number of species in the community.
- `C::Float64`: Connectance of the interaction matrix.
- `rand_ij::Bool = true`: whether interspecific interaction constants should be randomly sampled. If not they are set as `p[:μa0_ij]`.
- `rand_ii::Bool = true`: whether intraspecific interaction constants should be randomly sampled. If not they are set as `p[:μa0_ii]`. Only used if `intra == false`.
- `pNeg::Float64 = 0.0`: proportion of links which are assigned a negative value
- `sym::Bool = true`: whether the matrix should be symmetric
- `intra::Bool = true`: whether intraspecific interactions should be generated to ensure stability. This is done by ensuring diagonal domminance of the matrix.
- `nGrp::Int = 1`: The number of modules to create in the matrix
- `pMod::Float64 = 1.0`: The relative connectance between inter- and intra-modular links. Values closer to 0 indicate more links within the module.
"""
function mod_matrix(p::Dict{Symbol,Float64},nSp::Int,C::Float64; Tref::Float64 = 285.0,
                    rand_a0_ij::Bool = true, rand_a0_ii::Bool = true,
                    rand_Ea_ij::Bool = true, rand_Ea_ii::Bool = true,
                    pNeg::Float64 = 0.0, sym::Bool = true, intra::Bool = true,
                    nGrp::Int = 1, pMod::Float64 = 1.0)

    #check input
    @assert all(in.(key,Ref(keys(p)))) "Make sure all parameters are provided: \n $key"
    @assert nSp % nGrp == 0 "Number of species must be divisible by group size."
    @assert true "pMod and nGrp and C are not compatable"

    #set up a0 structure
    a0 = zeros(nSp,nSp)
    #work out group size and arrangement
    index = collect(1:nSp)
    grp_size = Int(nSp/nGrp)
    grps = [index[((grp_size*(i-1))+1) : (grp_size*(i)) ] for i in 1:nGrp]
    #work out inter- and intra-group connectance
    pGroup = grp_size / nSp #get the group size as a propotion
    aMin = -(pGroup - C) / (1-pGroup) # work out the smallest value C_exp = a * C_int allows
    a = aMin + ((1.0 - aMin) * pMod) #Get a via linear interp
    if isnan(a) a = 1.0 end #if grp size = 1 we need to manualy set a
    C_int = C / (a + pGroup*(1-a))
    C_ext = C_int * a

    for i = 1:nSp
        i_group = findall(in.(i,grps))
        for j = 1:nSp
            j_group = findall(in.(j,grps))
            if i_group == j_group
                if rand() < C_int
                    a0[i,j] = rand_a0_ij ? rand(TruncatedNormal(p[:μa0_ij],p[:σa0_ij],0,Inf)) : p[:μa0_ij]
                end
            else
                if rand() < C_ext
                    a0[i,j] = rand_a0_ij ? rand(TruncatedNormal(p[:μa0_ij],p[:σa0_ij],0,Inf)) : p[:μa0_ij]
                end
            end
        end
    end

    #set negative interactions
    neg_index = sample(1:nSp^2,Int(floor((nSp^2) * pNeg)))
    a0[neg_index] .= -a0[neg_index]

    if sym
    a0 = Array(Symmetric(a0))  #set symmetry
    end

    #set intra
    if intra
        a0[diagind(a0)] .= 0.0
        a0[diagind(a0)] .= sum(abs.(a0),dims=2)[:,1]
    elseif rand_a0_ii
        a0[diagind(a0)] .= rand(TruncatedNormal(p[:μa0_ii],p[:σa0_ii],0,Inf),nSp)
    else
        a0[diagind(a0)] .= p[:μa0_ii]
    end

    #E values
    #get all Ea
    Ea = rand_Ea_ij ? rand(TruncatedNormal(p[:μEa_ij],p[:σEa_ij],0,Inf),nSp,nSp) : fill(p[:μEa_ij],nSp,nSp)
    #set Ea_ii
    Ea[diagind(Ea)] .= rand_Ea_ii ? rand(TruncatedNormal(p[:μEa_ii],p[:σEa_ii],0,Inf),nSp) : p[:μEa_ii]
    #set connectance
    Ea[findall(a0 .== 0)] .= 0.0

    if sym
        Ea = Array(Symmetric(Ea))   #set symmetry
    end


    return(Boltz_params.(a0,Ea,Tref))
end





# function connectance(x,nSp)
#     sum(x .!= 0.0) / nSp^2
# end
