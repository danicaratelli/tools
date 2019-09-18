function Data_organizer(Data, Var_names, st_, en_, Transformations,T,n)
    #this function takes different series (the y's)
    #and organizes them in a matrix that will be used
    #for an OLS regression

    #The inputs are:
    #1) Data            ---> the data from the .csv file as is
    #2) Mnems           ---> the series needed in the VAR
    #3) st_             ---> the starting point
    #4) end_date        ---> the end point
    #5) Transformations ---> the transformations to apply to each series
    #6) T               ---> the size of the data
    #7) n               ---> the number of variables

    if st_!=1;
        Data = deleterows!(Data,1:st_-1);
    end
    if en_!=size(Data,1);
       Data = deleterows!(Data,en_+1:size(Data,1));
    end

    Y = repmat([NaN],length(Var_names),T);

    Data_cols = map(x->string(x),names(Data));

    if sum(Transformations.==intersect(Transformations,["lin","ln"]))!==n
        Y = Y[:,2:end];
    end
    for k=1:length(Var_names)
        loc_k = find(Data_cols.==Var_names[k])[1];
        tmp_data = float.(Data[:,loc_k]);

        if Transformations[k]=="lin";
            yy = lin(tmp_data);
        elseif Transformations[k]=="log";
            yy = ln(tmp_data);
        elseif Transformations[k]=="diff";
            yy = diff(tmp_data);
        elseif Transformations[k]=="pc";
            yy = pc(tmp_data);
        elseif Transformations[k]=="pcq";
            yy = pcq(tmp_data);
        end
        Y[k,:] = yy[end-size(Y,2)+1:end];
    end
    return Y;
end
