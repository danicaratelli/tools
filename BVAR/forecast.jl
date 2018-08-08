function forecast(M::model,h::Int64)
    n = M.n; #number of variables
    p = M.lags;
    Y = M.Y;
    B = M.coeff; #model coefficients
    C = [B[end,:]' zeros(1,n*(p-1))]'; #constant
    A = zeros(n*p,n*p);
    A[1:n,:] = B[1:end-1,:]';
    A[n+1:end,1:n*(p-1)] = eye(n*(p-1));

    A_tmp = eye(size(A,1),size(A,2));
    for i=1:h-1
        A_tmp = A_tmp+(A)^i;
    end

    fore = (A^h)*reshape(flipdim(Y[end-p+1:end,:],1)',n*p,1) + A_tmp*C;

    return fore[1:n];
end
