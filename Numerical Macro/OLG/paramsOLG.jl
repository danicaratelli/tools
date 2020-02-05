## author: daniele caratelli
## description: paramsOLG.jl containes the parameters and more for mainOLG.jl
## Heer-Maussner example 10


##  Parameters

    #preferences
β = 1.011;      #discount factor
η = 1.5;        #risk aversion
hbar = 0.3;     #hours worked
    #production
δ = 0.06;    #depreciation factor
α = 0.36;   #capital share
    #productivity
ρ = 0.96;       #productivity persistence
σz = 0.045;     #productivity variance
σ0 = 0.38;      #income variance at age 20
    #lifetime
T = 40;         #years of work
TR = 20;        #years of retirement
TT = T + TR;    #total lifespan

## Survival probabilities
    #survival probabilities (taken from Life Tables for the United States Social Security)
Ss  = 1 .- [0.00495 0.00082 0.00592 0.00993 0.01572 0.2234];
x = [0, 30, 60, 65, 70, 100];
spl = Spline1D(x, Ss[:]);
Ss = spl(collect(range(20,stop=100,step=1)));
Ss = min.(Ss,1);
Ss = [Ss[1:59];0];
Ss[1:findmax(Ss)[2]] .= 1;

## Productivity profile
Es = [0.6 1 1.075 1.15 1.08];
x = [20, 30, 40, 50, 60];
spl = Spline1D(x, Es[:]; w=ones(length(x)), k=1, bc="nearest", s=0.0)
Es = spl(collect(range(20,stop=59,step=1)));

#constructing productivity process using Tauchen
nz = 9;         #number of Tauchen discretization
M,zs = tauchen(nz, ρ, σz); zs = collect(zs);

## Initial distribution over productivities
zs_mid = mean([zs[1:end-1] zs[2:end]],dims=2);  #mid points for consecutive z in zs
cdf_mass = cdf.(Normal(0,σ0),[-Inf;zs_mid]);
y0_mass = diff(cdf.(Normal(0,σ0),[-Inf;zs_mid]),dims=1);
y0_mass = [y0_mass;1-sum(y0_mass)];             #weight on each element of zgrid
                                                #for income at t=20

## Grids
    #capital
Kmin = 1e-6;   #min assets
Kmax = 40;  #max assets
nk = 601;    #number of gridpoints for capital
kgrid = collect(range(Kmin,stop=Kmax,length=nk));
