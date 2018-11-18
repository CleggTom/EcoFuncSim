function simulate(c0::Vector{Float64},p::Parameters;t_start=0.0,t_end=1.0)

    prob = ODEProblem(Lotka,c0,(t_start,t_end),p)
    sol = solve(prob,Tsit5(),callback = PositiveDomain())

    return(sol)
end
