function imp_res(β,Σₑ,S,K,con);
    #inputs:
    #β = coefficients coming from the VAR (including the constant);
    #Σₑ = variance-covariance matrix of residuals
    #S = selection matrix
    #K = periods of impulse
    #con = 0 if the VAR was run with no regression, 1 if it
         # it was run with a constant term
    
    n = size(Σₑ,1);
    p = Int8((size(β,1)-con)/n);
    
    β¹ = β[1+con:end,:]; #selecting only the A's

    #constructing the A matrix
    #taking out first row from β since we already have the 
    #constant term
    #finding A of interest in β¹
    A = zeros(n*p,n*p);
    A[1:n,:] = β¹';
    A[n+1:end,1:n*(p-1)] = eye(n*(p-1));
    B₀ = chol(Σₑ)';

    output = zeros(K,n,n);
    #we apply a unitary shock to the orthogonalized system
    for k=1:K;
        tmp = (S*A^(k-1)*S'*B₀);
        for j=1:n;
            output[k,:,j] = tmp[:,j]'; #multiply by D¹
        end
    end
    
    #normalizing shocks to 1
    #for i=1:n;
    #    output[:,:,i] = (1./output[1,i,i])*output[:,:,i];
    #end

    return output;
end
