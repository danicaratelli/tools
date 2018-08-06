function var_m(data::Array{Float64,2},p::Int64)
    n = size(data,2);
    T = size(data,1);
    #organize Y (LHS of regression)
    Y = data[p+1:end,:];

    #organize X (RHS of regression)
    X = zeros(T-p,n*p);
    for i=1:p;
        X[:,(i-1)*n+1:i*n] = data[end-T+p-i+1:end-i,:];
    end
    X = [X ones(size(X,1),1)];

    ## VAR output
    #estimated parameters of VAR
    B = inv(X'*X)*X'*Y;
    #estimated variance-covariance matrix
    Σ = (Y-X*B)'*(Y-X*B)/size(Y,1);;
    #fitted values
    fits = X*B;

    M = model(B,Σ,fits,p,Y,X,data,T,n);
    return M;
end
