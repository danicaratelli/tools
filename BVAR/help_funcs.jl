type model
    coeff::Array{Float64,2};
    variance::Array{Float64,2};
    fits::Array{Float64,2};
    lags::Int64;
    Y::Array{Float64,2};
    X::Array{Float64,2};
    T::Int64;
    n::Int64;
end

#1st estimate the variance of residuals from AR(p) for each variable
function ar_reg(Y,p)    
    T, n = size(Y);
    S = zeros(n);
    for i=1:n
        y = Y[p+1:end,i];
        yy = zeros(size(y,1),p);
        for j=1:p
            yy[:,j] = Y[p+1-j:end-j,i];
        end
        b = yy\y;
        y_hat = yy*b;
        eps = y.-y_hat;
        S[i] = sqrt.((eps'*eps)/(T-n-p));
    end
    return S;
end
