function dummy_fores(Y,A,C,Σ,p,T,S);
    #description: this file computes a dummy forecast for Yₜ from p+1 to T adding a stochastic term 
    #"yₜ = c + Φ₁ yₜ₋₁ + ... + Φₚ yₜ₋ₚ + ϵₜ "
    #where ϵₜ∼𝒩(0,Σ)


     #inputs:
    #A = coefficients coming from the VAR (excluding the constant);
    #C = constant term from VAR;
    #Σₑ = variance-covariance matrix of residuals
    #S = selection matrix
    #K = periods of impulse
    #con = 0 if the VAR was run with no regression, 1 if it
         # it was run with a constant term
    d=MvNormal(zeros(size(Σ,1)),Σ); #normal distribution from which we pick ϵ's
    
    Y_out = zeros(n,T);
    
    Y_out[:,1:p] = Y;

    

    for k=p+1:T;
        E = rand(d,1);
        Y_in = reshape(flipdim(Y_out[:,k-p:k-1],2),n*p,1);#data we are 
        #forecasting with, contains p lags for n variables
        Y_out[:,k] = S*(forecast(A,C,Y_in,n,p,1))+E;
    end
    
    return Y_out; #this is the pretend 
end
