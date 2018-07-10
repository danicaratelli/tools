function irf_jl(M::model,h::Int64,K::Int64)
    T = M.T;
    n = M.n;
    p = M.lags;
    Y = M.Y[1:T-p,:];
    X = M.X[1:T-p,:];
    Y_star = M.Y;
    X_star = M.X;
    
    B = M.coeff; #model coefficients
    C = [B[end,:]' zeros(1,n*(p-1))]'; #constant
    A = zeros(n*p,n*p);
    A[1:n,:] = B[1:end-1,:]';
    A[n+1:end,1:n*(p-1)] = eye(n*(p-1));


    eps = M.Y.-M.fits;
    Su = (eps'*eps)/size(Y_star,1);

    CC = zeros(n*p,n);
    CC[1:n,1:n] = chol(Su)';
    JJ = [eye(n) zeros(n,n*(p-1))];
    Aj = eye(n*p);

    IRF = zeros(h+1,n,n);

    for j=0:h
        IRF[j+1,:,:] = JJ*Aj*CC*diagm(1./diag(CC));
        Aj = Aj*A;
    end

    # continue from here....
    
    iW = InverseWishart(n+n*(p-1)+2+2+T-(n*p+1),M.variance);
    PSIs = rand(iW,K);
    
    for k=1:K
        PSI = PSIs[k];
        B = chol(PSI)';
        D = diagm(diag(PSI));
        C = B*inv(D^(1/2));
        Sig = Symmetric(kron(PSI,inv(BVAR_model.X'*BVAR_model.X)));
        knorm = MvNormal(vec(BVAR_model.coeff),convert(Array{Float64,2},Sig));
    end
end
