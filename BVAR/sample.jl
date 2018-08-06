using  CSV, DataFrames, Dates, Distributions, PyPlot

## Options
λ = Inf;
ϵ = 0.00001;
D = [1 1 1];
p = 13;

## including files
include("help_funcs.jl");
include("bvar_m.jl");
include("var_m.jl");
include("forecast.jl");
include("irf.jl");


#load data
Data = CSV.read("Data.csv",dateformat="yyyy-mm-dd");
#transforming the data
Data[:Employment] = log.(Data[:Employment])*100;
Data[:Inflation] = log.(Data[:Inflation])*100;

data_full = convert(Array{Float64,2},Data[:,2:end]);

data = data_full[25:528,:];
#running Bayesian VAR

BVAR_model = bvar_m(data,p,λ,D[:]);
#running VAR
VAR_model = var_m(data,p);

K = 300;
irf_, irf_high, irf_low = irf(BVAR_model,48,λ,D[:],1)
irf, irf_high, irf_low, irf_med = irf(VAR_model,48,λ,D[:],1)
