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
#simulation details
Nsimul = 3000;  #number of simulations
Ngarbage = 500; #number of initial draws to throw out



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

## Step 4
@time cpol, Vfun, tol_res, iter_res = solveHH(agrid,Kgrid,Zs,Hmat,1);

#plotting consumption policy function for low and high capital and employed and unemployed
close()
figure()
#low capital
subplot(2,1,1)
plot(agrid,cpol[1,:,1],label=L"$\epsilon = e$")
plot(agrid,cpol[2,:,1],label=L"$\epsilon = u$")
title("K = "*string(Kgrid[1]))
legend()
xticks([])
#high capital
subplot(2,1,2)
plot(agrid,cpol[1,:,end],label=L"$\epsilon = e$")
plot(agrid,cpol[2,:,end],label=L"$\epsilon = u$")
title("K = "*string(Kgrid[end]))
subplots_adjust(hspace=0.5)
xlabel("assets")


## Step 5
# here we simulated the dynamics of the distribution function first by drawing a
# series of states according to the transition function. We initialize everything
# at the good state and cut off the first few observations

    #first getting simulations for aggregate state (starting from good state)
simulZ = Int.(zeros(Nsimul));
simulZ[1] = 1;
rnd_nums = rand(Nsimul);
tmpTmat = cumsum(ΓZ,dims=2);
for i=2:Nsimul
    simulZ[i] = findfirst(rnd_nums[i] .<= tmpTmat[simulZ[i-1],:]);
end

    #employment distribution over time and workers
dist_employment = Int.(zeros(NN,Nsimul));
dist_employment[1:Int(floor(ug*NN)),1] .= 1;    #initilization for employed/unemployed
rnd_nums = rand(NN,Nsimul);
tmpTmat = cumsum(Γ,dims=2);
for i=2:Nsimul
    tmats = tmpTmat[2*(simulZ[i-1])-1 .+ dist_employment[:,i-1],:];
    tmats[i] = map(x->findfirst(rnd_nums[x,i].<tmats[x,:]),1:NN);
    for ii=1:NN
        dist_employment[ii,i] =
    end
end
