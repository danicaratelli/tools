#=
filename:       HANK_main.jl
description:    Replicated model in Chapter 7.1 of Herr and Maussner, a simple
                heterogeneous agent model. It solves for the steady state
                (distributions and aggregates) in general equilibrium.
author:         Daniele Caratelli (danicaratelli@gmail.com)
=#

#*************** PACKAGES ***************#
using QuantEcon, Interpolations, Distributions, PyPlot
pygui(true)
include("HANK_help.jl")


#*************** INPUTS ***************#

#Employment-Unemployment transition matrix
puu = 0.5000; pue = 0.5000;
peu = 0.0435; pee = 0.9565
Tran = [puu pue;
        peu pee];

#Initial guess for capital stock K and tax rate tau
K0 = 5;
tau0 = 0.3;

#other parameters:
params = Dict();
params["β"] = 0.98;     #discount factor
params["α"] = 1/3;      #capital share
params["δ"] = 0.025;    #depreciation rate
params["η"] = 2;        #risk aversion parameter
params["b"] = 0.5;      #home production

#Asset grid
N = 10000;
params["N"] = N;
amin = -2;
amax = 3000;
a_grid = exp.(collect(range(log(1e-4),stop=log(amax-amin),length=N))) + amin;
params["a_grid"] = a_grid;

K0 = 250; tau0 = 0.3;
