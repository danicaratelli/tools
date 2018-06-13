using FredData
n = length(Var_names);
start_date = "1960-01-01";
end_date = "2009-01-01";
#fetching data
f = Fred("4f9f790c74d8b1e1f808400bc343fbe3")
Data = DataFrame();
for i=1:n
    df_tmp = get_data(f, Mnems[i];frequency="q",observation_start=start_date,observation_end=end_date,units=Transformations[i])
    if i==1
        Data[:date] = df_tmp.data[:date];
    end
    Data[Symbol(Var_names[i])] = df_tmp.data[:value];
end