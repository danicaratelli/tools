## author: daniele caratelli
## description: mainKS.jl solves the basic Krussel-Smith (1998) model as in
## Heer-Maussner example 8.3

using Dierckx, Statistics, PyPlot

##  Parameters
#grids
    #assets
amin = 0;   #min assets
amax = 12;  #max assets
na = 101;   #number of gridpoints for assets
    #capital
Kmin = 1.5;   #min assets
Kmax = 7.5;  #max assets
nk = 15;    #number of gridpoints for capital
#preference
β = 0.96;   #discount factor
η = 1.5;    #risk aversion
ζ = 0.15;   #share of wages to unemployed
#production
δ = 0.1;    #depreciation factor
α = 0.36;   #capital share
#states
Zb = 0.97;  #bad productivity level
Zg = 1.03;  #good productivity level
NN = 5000;   #number of households
#simulation details
Nsimul = 3000;  #number of simulations
Ngarbage = 500; #number of initial draws to throw out

include("helper.jl")


## Step 1 (assuming t= 0 is good state)
N0 = 0.95;
K0 = ((1/β - δ)/α)^(1/(α-1))*N0;
#computing the next period employment
N1 = 1 .- us;

## Step 2-3
Hmat = [0.2 0.9; 0.2 0.9];  #parameters for the law of motion of capital

## Step 4
Zs = [Zg;Zb];
cpol, apol, tol_res, iter_res = solveHH(agrid,Kgrid,Zs,Hmat,0,0);

plot_policies(agrid,Kgrid,cpol)


## Step 5
    #simulating aggregate states (1=good, 2=bad)
simulZ = simul_states(Nsimul);
    #simulating individual employment/unemployment states
dist_employment = simul_employment(NN,Nsimul,simulZ);
    #simulating distribution of HHs over time
@time A_sim, Kpath = simulate_HH_assets(NN,Nsimul,3,simulZ,dist_employment,apol,Kgrid,agrid);

## Step 6
simulZ_clean = simulZ[Ngarbage:Nsimul-2];
Kpath_clean = Kpath[Ngarbage:Nsimul-1];
#LOM parameters for capital
    #good state
ig = findall(simulZ_clean.==1);
ib = findall(simulZ_clean.==2);
lom_params_good = [ones(length(ig)) log.(Kpath_clean[ig])]\log.(Kpath_clean[ig.+1]);
lom_params_bad = [ones(length(ib)) log.(Kpath_clean[ib])]\log.(Kpath_clean[ib.+1]);
Hmat_new = [lom_params_good'; lom_params_bad'];

## Step 7
@time Hmat_star, dstar, itstar = estimate_LOM(agrid,Kgrid,Hmat,Zs,simulZ,dist_employment,Ngarbage,cpol,1e-6,20);

## Final Results
cpol_final, apol_final, tol_res, iter_res = solveHH(agrid,Kgrid,Zs,Hmat_star,0,0,0);
A_sim_final, Kpath_final = simulate_HH_assets(NN,Nsimul,5,simulZ,dist_employment,apol_final,Kgrid,agrid);

    #plot for capital
close()
plot(1:(Nsimul-Ngarbage),Kpath_final[Ngarbage:end-1]);
xlabel("time")
title("Capital over time",fontsize=14)
    #plot distribution of agents
close()
hist(mean(A_sim_final[:,1:(Nsimul-Ngarbage)],dims=2),bins=50)
xlabel("assets")
title("Distribution of asset by people (average over time)",fontsize=14)
