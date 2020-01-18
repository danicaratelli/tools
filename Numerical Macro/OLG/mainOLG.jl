## author: daniele caratelli
## description: mainOLG.jl solves an overlapping generations code as in
## Heer-Maussner example 10

using Dierckx, Statistics, PyPlot, Distributions
import Distributions: pdf, Normal, quantile
##  Parameters
    #years
T = 40;
TR = 20;
    #survival probabilities (taken from Life Tables for the United States Social Security)
Ss  = 1 .- [0.00495 0.00082 0.00592 0.00993 0.01572 0.2234];
x = [0, 30, 60, 65, 70, 100];
spl = Spline1D(x, Ss[:]);
Ss = spl(collect(range(20,stop=100,step=1)));
Ss = min.(Ss,1);
Ss = [Ss[1:59];0];
    #productivity profile
Es = [0.6 1 1.075 1.15 1.08];
x = [20, 30, 40, 50, 60];
spl = Spline1D(x, Es[:]; w=ones(length(x)), k=1, bc="nearest", s=0.0)
Es = spl(collect(range(20,stop=59,step=1)));
    #assets
kmin = 0;   #min assets
kmax = 12;  #max assets
na = 101;   #number of gridpoints for assets
    #capital
Kmin = 1.5;   #min assets
Kmax = 7.5;  #max assets
nk = 15;    #number of gridpoints for capital
    #preference
β = 1.011;      #discount factor
η = 1.5;        #risk aversion
hbar = 0.3;     #hours worked
    #production
δ = 0.06;    #depreciation factor
α = 0.36;   #capital share
    #productivity
ρ = 0.96;       #productivity persistence
σz = 0.045;     #productivity variance
nz = 9;         #number of Tauchen discretization
σ0 = 0.38;      #income variance at age 20
#simulation details
Nsimul = 3000;  #number of simulations
Ngarbage = 500; #number of initial draws to throw out

include("helper.jl")

#constructing productivity process
M,zs = tauchen(nz, ρ, σz); zs = collect(zs);
#constructing income distribution at t=20
zs_mid = mean([zs[1:end-1] zs[2:end]],dims=2);  #mid points for consecutive z in zs
cdf_mass = cdf.(Normal(0,σ0),[-Inf;zs_mid]);
y0_mass = diff(cdf.(Normal(0,σ0),[-Inf;zs_mid]),dims=1);
y0_mass = [y0_mass;1-sum(y0_mass)];             #weight on each element of zgrid
                                                #for income at t=20

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
