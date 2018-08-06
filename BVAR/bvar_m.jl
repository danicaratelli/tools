function bvar_m(data::Array{Float64,2},p::Int64,λ::Number,D::Array{Int64,1})
    T, n = size(data);

    S = ar_reg(data,p);
    #organize Y (LHS of regression)
    Y = data[p+1:end,:];

    #organize Y (LHS of regression)
    Y = data[p+1:end,:];

    #organize X (RHS of regression)
    X = zeros(T-p,n*p);
    for i=1:p;
        X[:,(i-1)*n+1:i*n] = data[end-T+p-i+1:end-i,:];
    end
    X = [X ones(size(X,1),1)];

        #creting dummy data
    Y_d = [diagm((D.*S)[:])./λ;
           zeros(n*(p-1),n);
           diagm(S);
           zeros(1,n)];
    X_d = [kron(diagm(1:p),diagm(S)./λ) zeros(n*p,1);
           zeros(n,n*p) zeros(n,1);
           zeros(1,n*p) ϵ];

    Y_star = [Y' Y_d']';
    X_star = [X' X_d']';

    ## BVAR output
    #estimated parameters of BVAR
    B = inv(X_star'*X_star)*X_star'*Y_star;
    #estimated variance-covariance matrix
    Σ = (Y_star-X_star*B)'*(Y_star-X_star*B)/size(Y_star,1);
    Σ = UpperTriangular(Σ)+ UpperTriangular(Σ)'.-diagm(diag(Σ)); #ensuring symmetry
    #fitted values
    fits = X_star*B;

    BVAR_model = model(B,Σ,fits,p,Y_star,X_star,data,T,n);

    return BVAR_model;
end
