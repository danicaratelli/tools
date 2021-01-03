## Simple Neoclassical Growth Model using Projection method and Chebyshev polynomials

using Statistics, PyPlot, Distributions, NLsolve
pygui(true);

#Projection options:
include("chebyshev.jl");
ncheby = 5;     #number of Chebyshev nodes
ndim = 2;       #number of state variables (capital and productivity)

#Parameters
β = 0.98;   #discount factor
α = 1/3;    #capital share
δ = 0.1;    #depreciation rate
ν = 0.5;    #labor-consumption elasticity
γ = 2;      #inverse EIS
#tauchen productivity shocks
σ = 0.007;      #variance of productivity shocks
λ = 0.9;    #persistence of productivity shocks

#grid boundaries
kmin = 1e-3; kmax = 20; kgrid = collect(range(kmin,stop=kmax,length=100));
zmin = -0.8; zmax = 0.8;

#Step 1: guess labor policy function by guessing coefficients
knodes = chebynodes(ncheby); #nodes for capital
ks = Z(knodes,kmin,kmax);
znodes = chebynodes(ncheby); #nodes for productivity
zs = Z(znodes,zmin,zmax);
zgrid = copy(zs);
#tauchen -> transition matrix
Π, zs = tauchen_mine(zs, λ, σ);

#labor policy function
#guess a constant labor = 0.33
lguess = 0.5*ones(ncheby,ncheby);

coeff_guess, vec_locs = chebycoeffs((ncheby,ncheby),lguess,(knodes,znodes));
lfun(coeff,ks,zs) = Cheby_eval(coeff,T_deg,(X(ks,kmin,kmax),X(zs,zmin,zmax)));

#Step 2: solve for consumption policy using static FOC
function cfun(coeff,ks,zs)
    #vectorized states (k and z)
    v_states = combine_vecs((ks,zs));
    cres = ((1-α).*exp.(v_states[:,2]).*v_states[:,1].^(α).*(lfun(coeff,ks,zs)[:]).^(-α)).*(1 .- lfun(coeff,ks,zs)[:]).*(ν/(1-ν));
    #bring back to un-vectorized shape
    cres = reshape(cres,length(ks),length(zs));
    return cres;
end

#Step 3: define residual (dynamic FOC)
lpol = lfun(coeff_guess,ks,zs);
cpol = cfun(coeff_guess,ks,zs);
function res(coeff,cpol,lpol,ks,zs,Π)
    nk = length(ks);
    nz = length(zs);

    #RHS of Euler Equation
    v_vars = [cpol[:] lpol[:]];
    v_states = combine_vecs((ks,zs));
    RHS = ν*(1 .- v_vars[:,2]).^((1-ν)*(1-γ)) .* (v_vars[:,1]).^(ν*(1-γ)-1);
    RHS = reshape(RHS,nk,nz);
    #LHS of Euler Equation
    knext = exp.(v_states[:,2]) .* v_states[:,1].^(α) .* v_vars[:,2].^(1-α) .+ (1-δ).*v_states[:,1] .- v_vars[:,1];
    #@assert(sum(knext.<0)==0,"Negative capital: problem with either coefficient guess or boundaries.")
    if sum(knext.<0)==0
        lnext = lfun(coeff,knext,zs);
        cnext = cfun(coeff,knext,zs);
        v_vars_next = [cnext[:] lnext[:]];
        v_states_next = combine_vecs((knext,zs));
        LHS1 = ((1-δ) .+ α*exp.(v_states_next[:,2]).*(v_states_next[:,1]).^(α).*(v_vars_next[:,2]).^(-α));
        LHS2 = ν*(1 .- v_vars_next[:,2]).^((1-ν)*(1-γ)) .* (v_vars_next[:,1]).^(ν*(1-γ)-1);
        LHS = LHS2.*LHS1;
        LHS = reshape(LHS,(nk,nz,nz));
        #computing E[LHS]
        E_LHS = zeros(nk,nz);
        for i=1:nz
            E_LHS[:,i] .+= LHS[:,:,i]*Π[i,:];
        end
        resid = ((RHS .- β*E_LHS).^2)[:];
    else
        resid = [1e6];
    end
end

T_deg = combine_vecs((collect(1:ncheby),collect(1:ncheby)));


function solve_projection(coeff_guess,Π,knodes,znodes,kgrid,zs)
    nk = length(knodes);
    nz = length(znodes);
    #coeff_guess = zeros(nk,nz); coeff_guess[1] = 0.4;

    Fres = [100.];
    function objfun!(Fres,x)
        coff_it = reshape(x,(nk,nz));
        #println("here")
        #println(coff_it[1])
        lpol = lfun(coff_it,kgrid,zs);
        if sum(lpol.<=0)==0
            cpol = cfun(coff_it,kgrid,zs);
            res_it = res(coff_it,cpol,lpol,kgrid,zs,Π);
            Fres[1] = maximum(res_it);
        else
            Fres[1] = 1e6;
        end
        #println(maximum(Fres))
    end
    xsol = nlsolve(objfun!,coeff_guess[:],iterations=2,show_trace=true,xtol=1e-6,ftol=1e-10);
    return xsol;
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

#NEED TO REWRITE THIS TO ACCOUNT FOR specific nodes...
function tauchen_mine(y::Array{Float64,1}, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::Integer=3) where {T1 <: Real, T2 <: Real}
    # Get discretized space
    N = length(y);
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
