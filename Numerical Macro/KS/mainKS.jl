## author: daniele caratelli
## description: mainKS.jl solves the basic Krussel-Smith (1998) model as in
## Heer-Maussner example 8.3

using Dierckx #Interpolations


##  Parameters
#grids
    #assets
amin = 0;   #min assets
amax = 12;  #max assets
na = 101;   #number of gridpoints for assets
    #capital
Kmin = 2.5;   #min assets
Kmax = 5.5;  #max assets
nk = 10;    #number of gridpoints for capital
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



## Step 1 (assuming t= 0 is good state)
N0 = 0.95;
K0 = ((1/β - δ)/α)^(1/(α-1))*N0;
#computing the next period employment
ug = 0.0386;    #probability of unemployment in good state
ub = 0.1073;    #probability of unemployment in bad state
us = [ug;ub];
N1 = 1 .- us;

## Step 2-3
Hmat = [0.2 0.9; 0.2 0.9];  #parameters for the law of motion of capital
