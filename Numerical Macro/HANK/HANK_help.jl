#=  filename:       HANK_help.jl
    description:    auxiliary functions to solve HANK model
    author:         Daniele Caratelli (danicaratelli@gmail.com)
=#

function solveHANK(K0::Number,tau0::Number,Tran::Array{Float64,2},params:Dict)

    #loading parameter values
    β = params["β"];
    α = params["α"];
    δ = params["δ"];
    η = params["η"];

    #step 1:
        #compute stationary employment L
    Tstationary = Tran^200;
    L = Tstationary[end,end];

    #step 2:
        #compute the wage w and the interest rate r based on guesses for capital
        #stock K
    w = (1-α)*(K0/L)^(α);
    R = α*(N/K0)^(1-α) - δ;

    #step 3:
        #solving HH's problem
    cpol, apol = HHegm(params,tau0,w,R-1,Tran);

    #step 4:
        #computing the stationary distribution of assets via simulation


end


function HHegm(params::Dict,τ::Number,w::Number,r::Number,Tran::Array{Float64,2},tol=1e-10,maxiter=5000)
    #Inputs:
        #params:    parameter dictionary
        #τ:         tax rate
        #w:         wage rate
        #r:         interest rate
        #Tran:      transition matrix
        #e_status:  employment status (1-> employed, 0->unemployed)

    #Outputs:
        #cpol:      consumption policy function
        #apol:      asset policy function

    println("Starting EGM...")

    #loading parameter values
    β = params["β"];
    α = params["α"];
    δ = params["δ"];
    η = params["η"];
    N = params["N"];
    a_grid = params["a_grid"];

    dist = Inf;
    iter = 0;
    cpol_old = ones(N,2);
    cpol_new = ones(N,2);
    while dist>tol && iter<maxiter
        cpol_tmp = (β*cpol_old.^(-η)*(1+(1-τ)*r)*Tran').^(-1/η);

        #assets in current period
        a_tmp = (a_grid .+ cpol_tmp .- [(1-τ)*w b])./(1+(1-τ)*r);
        for m=1:2

            #binding constraint
            loc_constrained = findall(a_grid.<=a_tmp[1,m]);
            cpol_new[loc_constrained,m] = a_grid[loc_constrained]*(1+(1-τ)*r) .+ (abs(m-2)*(1-τ)*w .+ b*(m-1)) .- a_grid[1];

            #non-binding constraint
            cgrid_end = findall(a_grid.>a_tmp[end,m]);
            #linear interpolation
            loc_nonconstrained = setdiff(1:N,union(loc_constrained,cgrid_end));
            knots = (a_tmp[:,m],);
            itp = interpolate(knots, cpol_tmp[:,m], Gridded(Linear()));
            cpol_new[loc_nonconstrained,m] = itp.(a_grid[loc_nonconstrained]);
            #corner
            cpol_new[cgrid_end,m] .= itp(a_tmp[end,m]);
        end

        dist = maximum(abs.(cpol_new.-cpol_old)[:]);
        if iter%100==0
            println("Iteration #"*string(iter)*", dif = "*string(dist));
        end
        iter+=1;
        cpol_old = copy(cpol_new); #updating policy function
    end

    #imputing savings policy function
    apol_new = a_grid.*(1+(1-τ)*r) .+ [(1-τ)*w b] .- cpol_new;
    println("Ending EGM.");
    println("\n");
    return cpol_new, apol_new;
end


function simul_dist()

end
