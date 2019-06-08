#=  filename:       helper.jl
    description:    auxiliary functions to solve HANK model
    author:         Daniele Caratelli (danicaratelli@gmail.com)
=#

function solveHANK(K0::Number,tau0::Number,Tran::Array{Float64,2},params:Dict)

    #loading parameter values
    β = params["β"];
    α = params["α"];
    δ = params["δ"];
    η = params["η"];

    #computing N from stationary distribution Tran
    Tstat = Tran^200;
    N = Tstat[end,end];

    #imply wage and interest rate from guesses
    w = (1-α)*(K0/N)^(α);
    r = α*(N/K0)^(1-α) - δ;

end


function HHegm(params::Dict,τ::Number,Tran::Array{Float64,2},tol=1e-10,maxiter=5000)
    #loading parameter values
    β = params["β"];
    α = params["α"];
    δ = params["δ"];
    η = params["η"];
    N = params["N"];

    dist = Inf;
    it = 0;
    cpol_old = ones(N,2);
    while dist>tol && it<maxiter
        
    end
end
