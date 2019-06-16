#=  filename:       HANK_help.jl
    description:    auxiliary functions to solve HANK model
    author:         Daniele Caratelli (danicaratelli@gmail.com)
=#

function solveHANK(K0::Number,tau0::Number,Tran::Array{Float64,2},params::Dict)

    #loading parameter values
    α = param["α"];
    δ = param["δ"];

    #step 1:
        #compute stationary employment L
    Tstationary = Tran^200;
    L = Tstationary[end,end];

    #step 2:
        #compute the wage w and the interest rate r based on guesses for capital
        #stock K
    K0 = 15
    w = (1-α)*(K0/L)^(α);
    r = α*(L/K0)^(1-α) - δ;

    #step 3:
        #solving HH's problem
    cpol, apol = HHegm(param,tau0,w,r,Tran);

    #step 4:
        #computing the stationary distribution of assets via simulation
    Fdist = simul_dist(apol,N,a_grid,Tran);

    #step 5:
        #aggregate capital and taxes
    Kstar = sum(Fdist.*[a_grid a_grid]);
    println(Kstar);
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
    β = param["β"];
    α = param["α"];
    δ = param["δ"];
    η = param["η"];
    N = param["N"];
    b = param["b"];
    a_grid = param["a_grid"];

    dist = Inf;
    iter = 0;
    cpol_old = ones(N,2);
    cpol_new = ones(N,2);
    while dist>tol && iter<maxiter
        cpol_tmp = (β*(1+(1-τ)*r)*cpol_old.^(-η)*Tran').^(-1/η);

        #assets in current period
        a_tmp = (a_grid .+ cpol_tmp .- [b (1-τ)*w])./(1+(1-τ)*r);

        for m=1:2

            #binding constraint
            loc_constrained = findall(a_grid.<=a_tmp[1,m]);
            cpol_new[loc_constrained,m] = a_grid[loc_constrained]*(1+(1-τ)*r) .+ (b*abs(m-2) + (m-1)*(1-τ)*w) .- a_grid[1];

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
    apol_new = a_grid.*(1+(1-τ)*r) .+ [b (1-τ)*w] .- cpol_new;
    println("Ending EGM.");
    println("\n");
    return cpol_new, apol_new;
end




function simul_dist(apol::Array{Float64,2},N::Int64,a_grid::Array{Float64,1},Tran::Array{Float64,2},tol=1e-10,maxiter=1000)

    #finding the policy for assets in index form
    apol_idx = Int.(zeros(apol));
    #first finding apol_idx for those at constraint
    locsi = findall(apol.==a_grid[1]);
    apol_idx[locsi] .= 1;

    for m=1:2
        if m==1
            st = locsi[end][1] + 1;
        else
            st = 1;
        end
        for n=st:N
            locsi = findlast(apol[n,m].>=a_grid)
            apol_idx[n,m] = locsi;
        end
    end

    #K = maxiter;
    #eps_rand = rand(K,N,2);
    Fdist = ones(N,2)./(2*N);
    dist = Inf;
    it = 0;
    while dist>tol && it<maxiter
        Fdist_new = zeros(N,2);
        for m=1:2
            us = unique(apol_idx[:,m],dims=1)
            for j=1:length(us)
               sumj = sum(Fdist[apol_idx[:,m].==us[j],m]);
                for mm = 1:2
                    Fdist_new[us[j],mm] = Fdist_new[us[j],mm] + Tran[m,mm]*sumj;#sum(Φold[locsj,m]);
                end
            end
        end
        dist = maximum(abs.(Fdist_new.-Fdist));
        Fdist = copy(Fdist_new);
        it +=1;
    end
    return Fdist;
end





#=
function invapol(apol,a_grid,N)
    amin = a_grid[1]; #borrowing limit

    #constructing finer version of grid
    #a_grid_mid = mean([a_grid[2:end] a_grid[1:end-1]],dims=2);
    #a_grid_fine = sort([a_grid;a_grid_mid],dims=1);
    #NN = length(a_grid_fine);

    apol_inv = zeros(N,2);
    apol_inv_idx = Int.(zeros(N,2));
    for m=1:2
        #filling in for those at borrowing constraint
        last_con = findlast(apol[:,m].<=amin);
        if last_con==nothing
            apol_inv[1,m] = a_grid[last_con];
            apol_inv_idx[1,m] = last_con;
        else
            last_con = 0;
        end

        #interpolating over monotone region
        for n=2:N
            ll = findlast(apol[:,m].<=a_grid[n]);
            if ll<N
                ww = (apol[ll+1,m] - a_grid[n])/(apol[ll+1,m]-apol[ll,m]);
                apol_inv[n,m] = ww*a_grid[ll] + (1-ww)*a_grid[ll+1];
            else
                apol_inv[n,m] = a_grid[ll];
            end
        end
    end
end
=#
