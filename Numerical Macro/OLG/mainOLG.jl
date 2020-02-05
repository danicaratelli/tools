## author: daniele caratelli
## description: mainOLG.jl solves an overlapping generations code as in
## Heer-Maussner example 10

using Dierckx, Statistics, PyPlot, Distributions, Roots, StatsBase

include("helper.jl")
include("paramsOLG.jl")

## Step 1 (computing aggregate employment)
μs = [1;map(x->prod(Ss[1:x]),1:TT-1)];
μs = μs.*(1/sum(μs));
N = sum(exp.(Es).*Ss[1:T].*μs[1:T]*hbar); #sum(hbar*exp.(Es).*μs[1:T]);
#sumc(sumc(muy.*exp(ye').*ef))*hbar/sumc(mass[1:t])

## Step 2 (guesses)
Kguess = 5;        #guessing aggregate capital (= guessing the interest rate)
τguess = 0.12;           #guessing labor income tax
## Steo 3
wguess = (1-α)*(Kguess/N)^(α);
bguess = b(Kguess,τguess);


## Solving HH problem once
k0 = 0;     #initial capital for new-born groups
## Find right aggregate capital
Kstar = clear_markets(kgrid,Kguess,τguess,zs,y0_mass,μs,k0,M,T,TT);
K_dist, Ks, Cs = distribution_forward(kgrid,Kstar,τguess,zs,y0_mass,k0,M,T,TT);
#assets age profile
plot_profile(Ks,"Assets",T,TT);
#consumption age profile
plot_profile(Cs,"Consumption",T,TT);

## Plotting Lorenz curve
Kall = sum(Ks);                     #total capital in economy
qs = [collect(0:0.1:0.9);0.99;1];   #quantiles of interest
plot_gini(K_dist,Kall,kgrid,qs)
