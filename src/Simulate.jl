function simulate(c0::Vector{Float64},p::Parameters;t_start=0.0,t_end=1.0)

    prob = ODEProblem(lotka_volterra,c0,(t_start,t_end),p)
    sol = solve(prob,Rosenbrock23())

    return(sol)
end



function mean_simulate(c0::Vector{Float64},p::Parameters;t_start=0.0,t_end=1.0)
    prob = ODEProblem(mean_lotka_volterra,c0,(t_start,t_end),p)
    sol = solve(prob,Rosenbrock23(),callback = PositiveDomain())

    return(sol)
end
