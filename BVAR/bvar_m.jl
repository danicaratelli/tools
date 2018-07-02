function bvar_m(data::Array{Float64,2},p::Int64,λ::Number,D::Array{Int64,1})
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
    
    ## Augmenting data with dummy data (to match prior)

    S = ar_reg(Y,p);
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
    Σ = (Y_star-X_star*B)'*(Y_star-X_star*B); 
    #fitted values
    fits = X_star*B;

    BVAR_model = model(B,Σ,fits,p,Y_star,X_star);
    
    return BVAR_model;
end
