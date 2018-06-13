function b_reg(X::Array{Float64,2},Y::Array{Float64,2},M=1000,cut=100,prior_mean="flat",prior_var="flat",λ=1);

    ## filename: shrink_reg
    ## author: Daniele Caratelli (danicaratelli@gmail.com)
    ## date = 03/07/2017

    ## file description this function takes in

    ## notes: running this function requires
    # Distributions


    ## Defining the distributions we are drawing from:

    # For flat prior
    # β |σ²,X,Y ∼ N(β_hat,(X'⋅X)¹⋅σ²), where β_hat = (X'⋅X)¹⋅X⋅Y
    # σ²|β,X,Y ∼ IG(n/2,-C/2), where C = νs² + (β-β_hat)'⋅X'⋅X⋅(β-β_hat)

    # For normal prior on β
    # β |σ²,X,Y ∼ N(β_bar,σ²Ωₙ), where β_bar = Ωₙ⋅(Ω₀⁻¹⋅β₀ +X'⋅X⋅β_hat)
    # β₀ is the mean of the prior distribution of β
    # Ω₀=λ⋅diag(cov(X)), Ωₙ=(Ω₀⁻¹ + X'⋅X)⁻¹
    # σ²|β,X,Y ∼ IG(n+1/2,C/2), where C = νs² + (β-β₀)'⋅X'⋅X⋅(β-β₀)
   
    
    #dimensions and degrees of freedom
    n = size(Y,1);
    k = size(X,2);
    ν = n-k;

    #estimated parameters
    β₀ = zeros(k,1);
    β_hat = ((X'*X)^-1)*X'*Y;
    if prior_mean=="normal";
        Ω₀ = λ*diagm(diag(cov(X)));
        Ωₙ = (Ω₀^-1 + X'*X)^-1;
        β_bar = Ωₙ*((Ω₀^-1)*β₀ + X'*X*β_hat);
    end

    #extra useful varaibles
    s2 = ((Y-X*β_hat)'*(Y-X*β_hat)/ν)[1];

    
    df = 10;
    B = zeros(M,k); Var = zeros(M);

    for j=1:M;
        if prior_mean=="flat" && prior_var=="flat"
            #drawing variance
            var = rand(InverseGamma(n/2,((ν*s2 + (β₀-β_hat)'*
                                          X'*X*(β₀-β_hat))/2)[1]),1)[1];
            #drawing coefficients
            β₀ = rand(MvNormal(β_hat[:],((X'*X)^-1)*var),1);
        elseif prior_mean=="normal" && prior_var=="flat"
            #drawing variance
            var = rand(InverseGamma((n+1)/2,
                                    ((ν*s2 +
                                      (β₀-β_bar)'*(Ωₙ^-1)*(β₀-β_bar) +
                                      (β_hat-β_bar)'*X'*X*(β_hat-β_bar) +
                                      (β_hat-β_bar)'*(Ω₀^-1)*(β_hat-β_bar))/2)[1]),1)[1];
            #drawing coefficients
            β₀ = rand(MvNormal(β_bar[:],var*Ωₙ),1);
        elseif prior_mean=="normal" && prior_var=="normal"
            #drawing variance
            var = rand(InverseGamma((n+df+1)/2,
                                    ((ν*s2 + 1 +
                                      (β₀-β_bar)'*(Ωₙ^-1)*(β₀-β_bar) +
                                      (β_hat-β_bar)'*X'*X*(β_hat-β_bar) +
                                      (β_hat-β_bar)'*(Ω₀^-1)*(β_hat-β_bar))/2)[1]),1)[1];
            #drawing coefficients
            β₀ = rand(MvNormal(β_bar[:],var*Ωₙ),1);
        else
            error("function cannot handle one of your prior distributions -- line 73");
        end
        Var[j] = var;
        B[j,:] = β₀;
    end

    return B[cut+1:end,:];
end
