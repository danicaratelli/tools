function RMSE(var_num,Actual,Predicted);
    Actual = Actual[:,var_num];
    Predicted = Predicted[:,var_num];
    nan_list1 = find(!isnan(Actual));
    nan_list2 = find(!isnan(Predicted));
    nan_list = intersect(nan_list1,nan_list2);
    Actual = Actual[nan_list];
    Predicted = Predicted[nan_list];
    if size(Actual,1)>size(Predicted,2);
        Actual = Actual[end-size(Actual,1)+1:end];
    end
    deviation = Actual.-Predicted;
    deviation = deviation[~isnan(deviation)];
    RMSE = sqrt(sum(deviation.^2));
    return RMSE;
end
