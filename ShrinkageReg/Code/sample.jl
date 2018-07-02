using Dates, FredData, Distributions

FRED_API_key = "4f9f790c74d8b1e1f808400bc343fbe3"; #St. Louis Fed Alfred jey
start_date = "2000-01-01";
end_date = "2016-12-01";
f = Fred(FRED_API_key);

GDP = get_data(f,"GDPC1",vintage_dates=today(),observation_start=start_date,observation_end=end_date,units="pca");

Employment = get_data(f,"PAYEMS",vintage_dates=today(),observation_start=start_date,observation_end=end_date,units="pca",frequency="q");

Consumption = get_data(f,"PCEC96",vintage_dates=today(),observation_start=start_date,observation_end=end_date,units="pca",frequency="q");


Y = convert(Array{Float64},GDP.df[:value]); Y=Y[:,:];
X = convert(Array{Float64},[Employment.df[:value] Consumption.df[:value]]);

include("b_reg.jl");

B = b_reg(X,Y,1000,100,"normal","normal",0.01);
