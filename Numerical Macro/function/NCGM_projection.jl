## Simple Neoclassical Growth Model using Projection method and Chebyshev polynomials

using Statistics, PyPlot, Distributions, NLsolve
pygui(true);

#Projection options:
include("chebyshev.jl");
ncheby = 6;     #number of Chebyshev nodes
ndim = 2;       #number of state variables (capital and productivity)

#Parameters
β = 0.98;   #discount factor
α = 1/3;    #capital share
δ = 0.1;    #depreciation rate
ν = 0.5;    #labor-consumption elasticity
γ = 2;      #inverse EIS
#tauchen productivity shocks
σ = 1;      #variance of productivity shocks
λ = 0.9;    #persistence of productivity shocks
Π, zgrid = tauchen(5, λ, σ); zgrid = collect(zgrid);
#grid boundaries
kmin = 1e-3; kmax = 20; kgrid = collect(range(kmin,stop=kmax,length=100));
zmin = zgrid[1]; zmax = zgrid[end];

#Step 1: guess labor policy function by guessing coefficients
knodes = chebynodes(ncheby); #nodes for capital
ks = Z(knodes,kmin,kmax);
znodes = chebynodes(ncheby); #nodes for productivity
zs = Z(znodes,zmin,zmax);
coeff_guess = zeros(ncheby,ncheby);
coeff_guess[1] = 0.5;

#Step 2: solve for consumption policy using static FOC
lfun(coeff,ks,zs) = Cheby_eval(coeff,T_deg,(X(ks,kmin,kmax),X(zs,zmin,zmax)));
function cfun(coeff,ks,zs)
    #vectorized states (k and z)
    v_states = combine_vecs((ks,zs));
    cres = ((1-α).*exp.(v_states[:,2]).*v_states[:,1].^(α).*(lfun(coeff,ks,zs)[:]).^(-α)).*(1 .- lfun(coeff,ks,zs)[:]).*(ν/(1-ν));
    #bring back to un-vectorized shape
    cres = reshape(cres,length(ks),length(zs));
    return cres;
end

#Step 3: define residual (dynamic FOC)
function res(coeff,cpol,lpol,ks,zs,ntauchen)
    nk = length(ks);
    nz = length(zs);

    #RHS of Euler Equation
    v_vars = [cpol[:] lpol[:]];
    v_states = combine_vecs((ks,zs));
    RHS = ν*(1 .- v_vars[:,2]).^((1-ν)*(1-γ)) .* (v_vars[:,1]).^(ν*(1-γ)+1);
    RHS = reshape(RHS,nk,nz);
    #LHS of Euler Equation
    knext = exp.(v_states[:,2]) .* v_states[:,1].^(α) .* v_vars[:,2].^(1-α) .+ (1-δ).*v_states[:,1] .- v_vars[:,1];
    lnext = lfun(coeff,knext,zs);
    cnext = cfun(coeff,knext,zs);
    v_vars_next = [cnext[:] lnext[:]];
    v_states_next = combine_vecs((knext,zs));
    LHS1 = ((1-δ).+ α*exp.(v_states_next[:,2]).*(v_states_next[:,1]).^(α-1).*(v_vars_next[:,2]).^(α))
    LHS2 = ν*(1 .- v_vars_next[:,2]).^((1-ν)*(1-γ)) .* (v_vars_next[:,1]).^(ν*(1-γ)+1);
    LHS = LHS2.*LHS1;
    LHS = reshape(LHS,(nk,nz,ntauchen));
    #computing E[LHS]
    E_LHS = zeros(nk,nz);
    for i=1:ntauchen
        E_LHS[:,i] .+= LHS[:,:,i]*Π[1,:];
    end

    resid = ((RHS .- E_LHS).^2)[:];
end

T_deg = combine_vecs((collect(1:ncheby),collect(1:ncheby)));


function solve_projection(Π,knodes,znodes)
    nk = length(knodes);
    nz = length(znodes);
    coeff_guess = zeros(nk,nz); coeff_guess[1] = 0.2;

    Fres = [zeros(1)];
    function objfun!(Fres,x)

    end

end



function clear_markets(agrid,thetas,y0_mass,a0,pars)
    #defining my objective function
    F = zeros(4);
    function objfun!(F,x)
        K = x[1];
        L = x[2];
        τL = x[3];
        b = x[4];
        if (K<0) | (L<0) | (τL>1) | (τL<0) | (b<0)
            F[1] = Inf;
            F[2] = Inf;
            F[3] = Inf;
            F[4] = Inf;
        else
            minW = w(K,L)*thetas[1];
            pars["L"] = L;
            pars["K"] = K;
            pars["τL"] = τL;
            pars["b"] = b;
            K_E_dist, K_U_dist, Ks, Cs, Ls =  distribution_forward(agrid,thetas,y0_mass,a0,pars);
            Y = (K^pars["α"])*(L^(1-pars["α"]));
            F[1] = abs(K - sum(Ks));                    #capital market clearing
            F[2] = abs(L - sum(Ls));                    #labor market clearing
            F[3] = abs(Y - (sum(Cs) + pars["d"]*K));      #goods market clearing
            F[4] = abs(b*sum(K_U_dist) - w(K,L)*τL*L); #government budget constraint
        end
    end

    xsol = nlsolve(objfun!, [pars["K"];pars["L"];pars["τL"];pars["b"]],iterations=100,show_trace=true,xtol=1e-6,ftol=1e-6);
    return xsol;
end

function tauchen(N::Integer, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::Integer=3) where {T1 <: Real, T2 <: Real}
    # Get discretized space
    a_bar = n_std * sqrt(σ^2 / (1 - ρ^2));
    y = range(-a_bar, stop=a_bar, length=N);
    d = y[2] - y[1]

    # Get transition probabilities
    Π = zeros(promote_type(T1, T2), N, N)
    dst = Normal(0,1);
    std_norm_cdf(x) = cdf(dst,x);
    for row = 1:N
        # Do end points first
        Π[row, 1] = std_norm_cdf((y[1] - ρ*y[row] + d/2) / σ)
        Π[row, N] = 1 - std_norm_cdf((y[N] - ρ*y[row] - d/2) / σ)

        # fill in the middle columns
        for col = 2:N-1
            Π[row, col] = (std_norm_cdf((y[col] - ρ*y[row] + d/2) / σ) -
                           std_norm_cdf((y[col] - ρ*y[row] - d/2) / σ))
        end
    end

    yy = y .+ μ / (1 - ρ) # center process around its mean (wbar / (1 - rho)) in new variable

    # renormalize. In some test cases the rows sum to something that is 2e-15
    # away from 1.0, which caused problems in the MarkovChain constructor
    Π = Π./sum(Π, dims = 2)

    return Π, yy;
end
