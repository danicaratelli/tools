type model
    coeff::Array{Float64,2};
    variance::Array{Float64,2};
    fits::Array{Float64,2};
    lags::Int64;
    Y::Array{Float64,2};
    X::Array{Float64,2};
    data::Array{Float64,2};
    T::Int64;
    n::Int64;
end
import Base.copy
Base.copy(m::model) = model(copy(m.coeff), copy(m.variance), copy(m.fits),
                            copy(m.lags), copy(m.Y), copy(m.X), copy(m.data),
                            copy(m.T), copy(m.n))

type irf_st
    median::DataFrame;
    qnt_high::DataFrame;
    qnt_low::DataFrame;
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
        yy = [yy ones(size(yy,1),1)]
        b = yy\y;
        y_hat = yy*b;
        eps = y.-y_hat;
        S[i] = sqrt.((eps'*eps)/(T-2p-1));
    end
    return S;
end
