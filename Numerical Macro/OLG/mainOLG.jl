## author: daniele caratelli
## description: mainOLG.jl solves an overlapping generations code as in
## Heer-Maussner example 10

using Dierckx, Statistics, PyPlot, Distributions, Roots
#import Distributions: pdf, Normal, quantile
##  Parameters
    #years
T = 40;
TR = 20;
TT = T + TR;
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
    #capital
Kmin = 1e-6;   #min assets
Kmax = 40;  #max assets
nk = 601;    #number of gridpoints for capital
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

## Step 1 (computing aggregate employment)
μs = [1;map(x->prod(Ss[1:x]),1:TT-1)];
μs = μs.*(1/sum(μs));
N = hbar*sum(μs[1:T]);


## Step 2 (guesses)
Kguess = 2.2;       #guessing aggregate capital (= guessing the interest rate)
τguess = 0.1;       #guessing labor income tax
## Steo 3
wguess = (1-α)*(Kguess/N)^(α);
bguess = b(Kguess,τguess);



## Solving HH problem once
@time Cs_sample, As_sample = solveHH(kgrid,Kguess,τguess,zs,M,T,TT);

## Solving forward for distributions of consumption and asseets given initial
#  distribution over productivity states y0_mass
@time K_dist, Ks, Cs = distribution_forward(kgrid,Kguess,τguess,zs,y0_mass,k0,M,T,TT);
