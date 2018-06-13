function var_reg(Data,p,con=1)
    #this function runs the VAR and outputs a β and Σ
    #the inputs are: i.   Y, the actual Data
    #                ii.  p, the number of lags
    #                ii.  con, 1 if constant included, 0 otherwise

   
    T = size(Data,2);
    n = size(Data,1);
    #organize Y (LHS of regression)
    Y = Data[:,p+1:end]';

    #organize X (RHS of regression)
    X = zeros(T-p,n*p);
    for i=1:p;
        X[:,(i-1)*n+1:i*n] = Data[:,end-T+p-i+1:end-i]';
    end

    ## Running regression (normal OLS) ##

    #adding constant term
    if con==1;
        X = [ones(size(X,1),1) X];
    end

    β = (X'*X)^-1*X'*Y; #coefficients of regression
    Σ = (Y-X*β)'*(Y-X*β)/(T-p); #estimated variance of ϵ
    Σ = floor(Σ*1e12)/1e12;
    Fits = X*β;
    
    dict = Dict("coeff" => β, "var" => Σ, "dep" => Y, "indep" => X, "fits" => Fits);
    return dict;
end
    
