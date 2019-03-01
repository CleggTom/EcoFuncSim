##First We Load the EcoFuncSim Package
using EcoFuncSim
eco = EcoFuncSim

#we then set up some simulation parameters
n_Sp = 100 #number of species
Temp = 285.0 #simulation temperature in Kelvin
max_t = 10.0 #Time to run the simulation too.
C = 0.5 #network connectance
C0 = ones(n_Sp) #vector of starting biomasses


#parameter distributions see ?eco.key for more information.
params = Dict(:μU0=>100.0,:σU0=>10.0,:μEU=>0.32,:σEU=>0.1,
             :μR0=>50.0 ,:σR0=>10.0,:μER=>0.65,:σER=>0.1,
             :μa0_ii=>1.0,:σa0_ii=>0.3,:μEa_ii=>0.5,:σEa_ii=>0.1,
             :μa0_ij=>1.0 / (n_Sp-1),:σa0_ij=>0.3/(n_Sp-1),:μEa_ij=>0.5,:σEa_ij=>0.1)

#we can now use the functions in the Package to generate the actual community parameters
#For more information look at the documentation (?eco.f(...))

#first the intrinsic growth vectors
growth = eco.rand_growth(params,n_Sp)

#and the interaction matrix (note there is also a mod_matrix function to generate modular interaction matricies)
a = eco.rand_matrix(params,n_Sp,C)

#Then we convert this into a Parameter object at the given temperature
p = make_params(growth[:,1],growth[:,2],a,Temp)

#which is then used to run the simulation
#This returns the solution object as documented in the DiffEqs parameters
sol = eco.simulate(C0,p;t_start=0.0,t_end=max_t)

#We can use this to get the solution at a given time
sol(1.0)

#or the whole soution array
sol.u
