## author: daniele caratelli
## description: helper.jl contained the functions used in mainOLG.jl
## Heer-Maussner example 10

#useful functions
function u(c,η)
    if η==1
        return ln.(c);
    else
        return (c.^(1-η))/(1-η);
    end
end

function du(c,η)
    if η==1
        return 1 ./c;
    else
        return c.^(-η);
    end
end

function duinv(x,η)
    if η==1
        return 1 ./x;
    else
        return x.^(-1/η);
    end
end

w(K) = (1-α)*(K/N)^(α);
r(K) = α*(K/N)^(α-1) - δ;
b(τ,K) = (τ*w(K)*N)/sum(μs[T+1:TT]);

std_norm_cdf(x::Number) = cdf(Normal(0,1),x);

##-----     HOUSEHOLD PROBLEM     -----##
"""
    solveHH(kgrid,K,zs,M,T,TR,TT)

    Solving household problem via Endogenous Grid Method.

   Inputs:
   -------
   ``kgrid``    :       array (nk),  grid for assets\n
   ``K``        :       Number,  aggregate capital\n
   ``τ``        :       Number,  income tax rate\n
   ``zs``       :       array (nZ),  aggregate productivities\n
   ``M``        :       array (nZ,nZ),   transition of individual productivity\n
   ``T``        :       Number,  years of working life\n
   ``TT``       :       Number,  total life = T+TR\n
   -------
   ```
"""
function solveHH(kgrid,K,τ,zs,M,T,TT)
    nk = length(kgrid);
    nz = length(zs);

    #initiliazing matrices
        #consumption
    Cs = zeros(nk,nz,TT);     #1st dim -> capital
                                #2nd dim -> productivity shock
                                #3rd dim -> age
        #assets
    As = zeros(size(Cs));

    #last period of life consumption
    Cs[:,:,end] .= (1+r(K)).*kgrid .+ b(τ,K);
    #structures for interpolation
    xqi = zeros(nk);
    xqia = zeros(nk);
    xqpi = zeros(nk);

    for t=TT-1:-1:T+1
        Exp_V = du(Cs[:,:,t+1],η)*M';
        c_prev = duinv(β*Ss[t]*(1+r(K))*Exp_V,η);

        k_prev = (kgrid .+ c_prev .- b(τ,K))/(1+r(K));
        #because in retirement productivity does not matter, I just use iz = 1
        kprev_low_int, kprev_high_int, w_low_int = interpolate_coord(k_prev[:,1],kgrid,xqi,xqia,xqpi);
        Cs[:,:,t] .= c_prev[kprev_low_int,:].*w_low_int .+ c_prev[kprev_high_int,:].*(1 .- w_low_int);

        #binding borrowing constraint
        iconstained = (kgrid .<= k_prev[1,1]); #points which are at the borrowing constraint
        if sum(iconstained)>0
            Cs[iconstained,:,t] .= (1 + r(K))*kgrid[iconstained] .+ b(τ,K) .-kgrid[1];
        end
    end
    As[:,:,T+1:TT] = (1 + r(K)).*kgrid .+ b(τ,K) .- Cs[:,:,T+1:TT];

    #working life
    for t=T:-1:1
        Exp_V = du(Cs[:,:,t+1],η)*M';
        c_prev = duinv(β*Ss[t]*(1+r(K))*Exp_V,η);
        k_prev = (kgrid .+ c_prev .- w(K)*hbar*(1-τ)*exp.(Es[t].+repeat(zs',nk,1)))/(1+r(K));
        #looping over productivity state
        for iz=1:nz
           kprev_low_int, kprev_high_int, w_low_int = interpolate_coord(k_prev[:,iz],kgrid,xqi,xqia,xqpi);
           Cs[:,iz,t] = c_prev[kprev_low_int,iz].*w_low_int .+ c_prev[kprev_high_int,iz].*(1 .- w_low_int);

           #binding borrowing constraint
           iconstained = (kgrid .<= k_prev[1,iz]); #points which are at the borrowing constraint
           if sum(iconstained)>0
              Cs[iconstained,iz,t] = (1 + r(K))*kgrid[iconstained] .+ w(K)*hbar*(1-τ)*exp(Es[t]*zs[iz]) .- kgrid[1];
          end
        end
       As[:,:,t] = (1 + r(K)).*kgrid .+ w(K)*hbar*(1-τ)*exp.(Es[t].+repeat(zs',nk,1)) .- Cs[:,:,t];
    end

    return Cs, As;
end


##-----     DISTRIBUTION FORWARD     -----##
"""
    distribution_forward(kgrid,K,τ,zs,y0_mass,k0,M,T,TT)

    Updating capital distribution over time.

   Inputs:
   -------
   ``kgrid``    :       array (nk),  grid for assets\n
   ``K``        :       Number,  aggregate capital\n
   ``τ``        :       Number,  income tax rate\n
   ``zs``       :       array (nz),  aggregate productivities\n
   ``y0_mass``  :       array (n),  initial distribution over zs\n
   ``k0``       :       Number,  initial asset holdings\n
   ``M``        :       array (nz,nz),   transition of individual productivity\n
   ``T``        :       Number,  years of working life\n
   ``TT``       :       Number,  total life = T+TR\n
   -------
   ```
"""
function distribution_forward(kgrid,K,τ,zs,y0_mass,k0,M,T,TT)
    nk = length(kgrid);
    nz = length(zs);

    #initializing distribution
    K_dist = zeros(nk,nz,TT);   #assets
    C_dist = zeros(nk,nz,TT);   #consumption

    #aggregate variables
    Ks = zeros(TT);
    Cs = zeros(TT);

    #structures for interpolation
    xqi = zeros(nk);
    xqia = zeros(nk);
    xqpi = zeros(nk);

    #HH problem to guide update in distribution
    Cs_sample, As_sample = solveHH(kgrid,Kguess,τguess,zs,M,T,TT);

    #period t=1 distribution

    #allocating agents according to initial capital holding k0
    loc0, w0 = interpolate1D(k0,kgrid);
    K_dist[loc0,:,1] .= w0*y0_mass[:]*μs[1];
    K_dist[loc0+1,:,1] .= (1 - w0)*y0_mass[:]*μs[1];
    Ks[1] = sum(K_dist[:,:,1].*kgrid);
    Cs[1] = sum(K_dist[:,:,1].*Cs_sample[:,:,1]);

    #updating distribution
    for t=2:TT
        for iz=1:nz
            locs_low_t, locs_high_t, wt = interpolate_coord(kgrid,As_sample[:,iz,t-1],xqi,xqia,xqpi);
                #points on lower bound of interpolation
            dups_low = findall(diff(locs_low_t).==0) .+ 1;
            nodup = setdiff(1:nk,dups_low);
            locs_low_t_nodup = locs_low_t[nodup];
            K_dist[locs_low_t_nodup,:,t] = K_dist[locs_low_t_nodup,:,t] .+ Ss[t-1]*K_dist[nodup,iz,t-1].*wt[nodup].*repeat(M[iz,:]',length(nodup),1);
                #taking care of doubles
            for id in dups_low
                K_dist[locs_low_t[id],:,t] = K_dist[locs_low_t[id],:,t] .+ Ss[t-1]*K_dist[id,iz,t-1].*wt[id].*M[iz,:];
            end
                #points on upper bound of interpolation
            dups_high = findall(diff(locs_high_t).==0) .+ 1;
            nodup = setdiff(1:nk,dups_high);
            locs_high_t_nodup = locs_high_t[nodup];
            K_dist[locs_high_t_nodup,:,t] = K_dist[locs_high_t_nodup,:,t] .+ Ss[t-1]*K_dist[nodup,iz,t-1].*(1 .- wt[nodup]).*repeat(M[iz,:]',length(nodup),1);
                #taking care of doubles
            for id in dups_high
                K_dist[locs_high_t[id],:,t] = K_dist[locs_high_t[id],:,t] .+ Ss[t-1]*K_dist[id,iz,t-1].*(1 - wt[id]).*M[iz,:];
            end
        end
        Ks[t] = sum(K_dist[:,:,t].*kgrid);
        Cs[t] = sum(K_dist[:,:,t].*Cs_sample[:,:,t]);
    end
    return K_dist, Ks, Cs
end


##-----     MARKET CLEARING     -----##
"""
    clear_markets(kgrid,K,τ,zs,y0_mass,k0,M,T,TT)

    Finding market-clearing prices for GE.

   Inputs:
   -------
   ``kgrid``    :       array (nk),  grid for assets\n
   ``K``        :       Number,  guess of aggregate capital\n
   ``τ``        :       Number,  guess of income tax rate\n
   ``zs``       :       array (nz),  aggregate productivities\n
   ``y0_mass``  :       array (n),  initial distribution over zs\n
   ``y0_mass``  :       array (n),  initial distribution over zs\n
   ``k0``       :       Number,  initial asset holdings\n
   ``M``        :       array (nz,nz),   transition of individual productivity\n
   ``T``        :       Number,  years of working life\n
   ``TT``       :       Number,  total life = T+TR\n
   -------
   ```
"""
function clear_markets(kgrid,K,τ,zs,y0_mass,μs,k0,M,T,TT)
    function obj_fun(x)
        K_dist, Ks, Cs = distribution_forward(kgrid,x[1],τguess,zs,y0_mass,k0,M,T,TT);
        #sol = abs(Ks[:]'*μs[:] - x);
        sol = abs(sum(K_dist.*kgrid) - x);
        return sol;
    end

    xsol = find_zero(obj_fun,K);
    return xsol;
end


##-----     INTERPOLATIONS FUNCTIONS     -----##
"""
    interpolate1D(xq, xgrid)

    Linear interpolate in 1 dimension.

   Inputs:
   -------
   ``x``        :       Number,     number for placement on the grid\n
   ``xgrid``    :       array (n),  grid on which x is to be placed\n
   -------
   ```
"""
function interpolate1D(x,xgrid)
    nx = length(xgrid);
    if x<xgrid[1]
        xloc_min = 1;
        wx = 1;
    elseif x>xgrid[end]
        xloc_min = nx-1;
        wx = 0;
    else
        xloc_min = findlast(x.>=xgrid);
        wx = 1 - (x - xgrid[xloc_min])/(xgrid[xloc_min + 1] - xgrid[xloc_min]);
    end
    return xloc_min, wx;
end

"""
    interpolate_coord(x, x1, xqi, xqia ,xqpi)

    Linear interpolate:  `xq = xqpi * x[xqi] + (1-xqpi)*x[xqia]`.

    Code converted from Using the Sequence-Space Jacobian to Solve and Estimate Heterogeneous-Agent Models

   Inputs:
   -------
   ``x``        :       array (n), ascending data points\n
   ``xq``       :       array (nq), query points\n
   ``xqi``      :       array (nq), empty (to be filled with indices of lower bracketing gridpoints)\n
   ``xqia``     :       array (nq), empty (to be filled with indices of upper bracketing gridpoints)\n
   ``xqpi``     :       array (nq), empty (to be filled with weights of lower bracketing gridpoints)\n
   -------
   ```
"""
function interpolate_coord(x,xq,xqi,xqia,xqpi)
    #size of arrays
    nxq, nx = size(xq,1), size(x,1);

    #sort and keep track of initial order
    ind_new = sortperm(xq);
    ind_init = sortperm(ind_new);
    xq = xq[ind_new];

    #take care of value below and above minimum
    id = findall((x[1] .<= xq) .& (xq .< x[end]));
    xqi[(xq .< x[1])] .= 1;
    xqpi[(xq .< x[1])] .= 1;
    xqi[(xq .> x[nx])] .= nx;
    xqpi[(xq .> x[nx])] .= 1;

    #interpolation
    xi = 1;
    x_low = x[1];
    x_high = x[2];

    for xqi_cur in id
        xq_cur = xq[xqi_cur];
        while xi < (nx - 1)
           if x_high>=xq_cur
                break
            end
            xi += 1
            x_low = x_high;
            x_high = x[xi + 1];
        end
        xqpi[xqi_cur] = (x_high - xq_cur)/(x_high - x_low);
        xqi[xqi_cur] = xi;
    end

    # revert back to initial order
    xqpi[:] = xqpi[ind_init];
    xqi[:] = xqi[ind_init];

    # Compute index of point above, or same if last on the list
    xqia[:] = xqi[:] .+ 1
    xqia[(xqia .>= nx)] .= nx;
    xqia[(xq .< x[1])] .= xqi[(xq .< x[1])];

    return Int.(xqi), Int.(xqia), xqpi;
end


##-----     AGE PROFILE     -----##
function plot_profile(X,xname,T,TT)
    figure()
    plot(20:20+TT-1,X)
    ys = ylim();
    vlines(20+T,ys[1],ys[2],linestyle="dotted")
    xlim(20,20+TT-1);
    ylim(ys[1],ys[2])
    xlabel("Age");
    ylabel(xname);
end


##-----     LORENZ CURVE     -----##
function plot_gini(K_dist,Kall,kgrid,qs)
    Kdist_all = sum(K_dist,dims=(2,3))[:];
    Kdist_all = Kdist_all./sum(Kdist_all);
    #quantiles for assets
    worker_qnts = map(x->findfirst(cumsum(Kdist_all).>=x),qs);

    figure()
    plot(qs,map(x->sum(Kdist_all[1:x].*kgrid[1:x])/Kall,worker_qnts),marker="*",label="model")
    plot(qs,qs,marker="s",label="equal distribution")
    xlabel("proportion of Workers");
    ylabel("proportion of Asset");
    legend()
end


##-----     TAUCHEN DISCRETIZATION (taken from Quantecon)     -----##
"""
    tauchen(N::Integer, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::Integer=3)

    Tauchen discretization taken form Quantecon

   Inputs:
   -------
   - `N::Integer`: Number of points in markov process
   - `ρ::Real` : Persistence parameter in AR(1) process
   - `σ::Real` : Standard deviation of random component of AR(1) process
   - `μ::Real(0.0)` : Mean of AR(1) process
   - `n_std::Integer(3)` : The number of standard deviations to each side the process
   -------
"""
function tauchen(N::Integer, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::Integer=3) where {T1 <: Real, T2 <: Real}
    # Get discretized space
    a_bar = n_std * sqrt(σ^2 / (1 - ρ^2))
    y = range(-a_bar, stop=a_bar, length=N)
    d = y[2] - y[1]

    # Get transition probabilities
    Π = zeros(promote_type(T1, T2), N, N)
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
