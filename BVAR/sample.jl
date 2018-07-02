using  CSV, DataFrames, Dates

## Options
λ = 10000;
ϵ = 0.00001;
D = [1 1 0];
p = 4;

## including files
include("help_funcs.jl");
include("bvar_m.jl");
include("var_m.jl");
include("forecast.jl");

#load data
Data = CSV.read("Data.csv",dateformat="yyyy-mm-dd");

data = convert(Array{Float64,2},Data[:,2:end]);

#running Bayesian VAR
BVAR_model = bvar_m(data,p,λ,D[:]);
#running VAR
VAR_model = var_m(data,p);

#forecasting models h periods ahead
bYh = forecast(BVAR_model,2);
Yh = forecast(VAR_model,2);