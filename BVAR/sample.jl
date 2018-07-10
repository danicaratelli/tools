using  CSV, DataFrames, Dates, Distributions


## Options
λ = 1/15.4564;
ϵ = 0.00001;
D = [1 1 1];
p = 13;

## including files
include("help_funcs.jl");
include("bvar_m.jl");
include("var_m.jl");
include("forecast.jl");


#load data
Data = CSV.read("Data.csv",dateformat="yyyy-mm-dd");
#transforming the data
Data[:Employment] = log.(Data[:Employment])*100;
Data[:Inflation] = log.(Data[:Inflation])*100;

data_full = convert(Array{Float64,2},Data[:,2:end]);

data = data_full[38:end-12,:];
#running Bayesian VAR
BVAR_model = bvar_m(data,p,λ,D[:]); 
#running VAR
VAR_model = var_m(data,p);

#forecasting models h periods ahead
bYh = forecast(BVAR_model,1);
Yh = forecast(VAR_model,2);



K = 10_000;
