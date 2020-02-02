##

# constructing grids
kgrid = collect(range(Kmin,stop=Kmax,length=nk));


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
#=function b(K)
    wtmp = w(K)
    f(b) = b - 0.3*wtmp*hbar*(1 - (b/(wtmp*N))*sum(μs[T+1:TT]))*sum(μs[1:T].*Es[1:T]);
    b = find_zero(f, (0, wtmp), Bisection(), atol=1e-8);
    return b;
end
function τ(K)
    wtmp = w(K);
    btmp = b(K);
    return (btmp*sum(μs[T+1:TT]))/(wtmp*N)
end=#

std_norm_cdf(x::Number) = cdf(Normal(0,1),x);

##-----     HOUSEHOLD PROBLEM     -----##
"""
    solveHH(kgrid,K,zs,M,T,TR,TT)

    Solving household problem via Endogenous Grid Method.

   Inputs:
   -------
   ``kgrid``    :       array (na),  grid for assets\n
   ``K``        :       Number,  aggregate capital\n
   ``zs``       :       array (nZ),  aggregate productivities\n
   ``M``        :       array (nZ,nZ),   transition of individual productivity\n
   ``T``        :       Number,  years of working life\n
   ``TR``       :       Number,  years of retirement\n
   ``TT``       :       Number,  total life = T+TR\n
   -------
   ```
"""
function solveHH(kgrid,K,τ,zs,M,T,TR,TT)
    nk = length(kgrid);
    nz = length(zs);

    #initiliazing matrices
    Cs = zeros(nk,nz,TT);       #1st dim -> capital
                                #2nd dim -> productivity shock
                                #3rd dim -> age
    Vs = zeros(nk,nz,TT);
    Ks = zeros(nk,nz,TT);
    #last period of life consumption and last period utility
    Cs[:,:,end] .= (1+r(K)).*kgrid .+ b(τ,K);
    Vs[:,:,end] .= u(Cs[:,:,end],η);

    #structures for interpolation
    xqi = zeros(nk);
    xqia = zeros(nk);
    xqpi = zeros(nk);

    #useful derivatives of kgrid
    kstep = kgrid[2] - kgrid[1];
    kmidpoint = (kgrid[2:end] .+ kgrid[1:end-1])./2;


    #Value function Iteration
    cgrid = (1+r(K))*repeat(kgrid,1,nk) .- repeat(kgrid',nk,1) .+ b(τ,K);
    c0s = (cgrid .<0);
    cgrid[c0s] .= 100;
    Us = u(cgrid,η);
    Us[c0s] .= -Inf;

    #retirement
    for t=TT-1:-1:T+1
        Vnext = β*Ss[t].*(Vs[:,:,t+1]*M');
        for iz=1:nz
            vs, locs = findmax(Us .+ repeat(Vnext[:,iz]',nk,1),dims=2);
            locs = map(x->x[2],locs[:]);
            Ks[:,iz,t] = kgrid[locs];
            Vs[:,iz,t] = vs;
        end
        Cs[:,:,t] = (1+r(K))*kgrid .- Ks[:,:,t] .+ b(τ,K);
    end

    #working life
    for t=T:-1:1
        Vnext = β*Ss[t].*(Vs[:,:,t+1]*M');
        for iz=1:nz
            cgrid = (1+r(K))*repeat(kgrid,1,nk) .- repeat(kgrid',nk,1) .+ (1-τ)*w(K)*hbar*exp(Es[t] + zs[iz]);
            c0s = (cgrid .<0);
            cgrid[c0s] .= 100;
            Us = u(cgrid,η);
            Us[c0s] .= -Inf;

            vs, locs = findmax(Us .+ repeat(Vnext[:,iz]',nk,1),dims=2);
            locs = map(x->x[2],locs[:]);
            Ks[:,iz,t] = kgrid[locs];
            Vs[:,iz,t] = vs;
            Cs[:,iz,t] = (1+r(K))*kgrid .- Ks[:,iz,t] .+ (1-τ)*w(K)*hbar*exp(Es[t] + zs[iz]);
        end
    end

    #=@time for t=TT-1:-1:T+1 #;-1:1 #TT-1:-1:T+1
        if t>T      #retirement
            #next period value function
            spl = Spline1D(kgrid,Vs[:,1,t+1]);
            #numerical derivative of nexst period value function
            dVt_num = diff(Vs[:,1,t+1],dims=1)./kstep;
            dspl = Spline1D(kmidpoint,dVt_num;bc="extrapolate");

            dVR(k,c) = evaluate(dspl, (1+r(K))*k .+ b(τ,K) .- c);
            #computing consumption at t using FOC of HH problem
            function solveC_ret(k)
                objC(c) = du(c,η) - β*Ss[t]*dVR(k,c)[1];
                solC = find_zero(objC,k/((TT-t)+1));
                return solC;
            end

            #computing consumption c_t
            Cs[:,:,t] .= map(x->solveC_ret(x),kgrid);
            #computing assets k_{t+1}
            Ks[:,:,t] .= (1+r(K))*kgrid .+ b(τ,K) .- Cs[:,1,t];
            iconstrained = findall(Ks[:,1,t] .< Kmin);
            if length(iconstrained)>0
                Ks[iconstrained,:,t] .= kgrid[1];
                Cs[iconstrained,:,t] .= (1 + r(K)).*kgrid[iconstrained] .+ b(τ,K) .- Ks[iconstrained,1,t];
            end
            Vs[:,:,t] .= u(Cs[:,1,t],η) .+ (β*Ss[t])*evaluate(spl, Ks[:,1,t]);
        else        #working life
            for iz=1:nz
                #next period value function
                spl = Spline2D(kgrid,zs,Vs[:,:,t+1]);
                #numerical derivative of next period value function
                dVt_num = diff(Vs[:,:,t+1],dims=1)./kstep;
                dspl = Spline2D(kmidpoint,zs,dVt_num);

                dVW(k,c,z) = evalgrid(dspl, (1+r(K))*[k] .+ (1-τ)*w(K)*hbar*exp.(Es[t]+z) .- c, zs);
                #computing consumption at t using FOC of HH problem
                function solveC_work(k)
                    objC(c) = du(c,η) - (β*Ss[t]*dVW(k,c,zs[iz])*M[iz,:])[1];
                    solC = find_zero(objC,0.1);#k/((TT-t)+1));
                    return solC;
                end
                #computing consumption c_t
                Cs[:,iz,t] = map(x->solveC_work(x),kgrid);
                #linear extrapolate if necessary
                cdoubles = findall(Cs[:,iz,t].==Cs[end,iz,t]);
                if length(cdoubles)>1
                    m = (Cs[cdoubles[1]-1,iz,t] - Cs[cdoubles[1] - 2,iz,t])/(kgrid[cdoubles[1]-1] - kgrid[cdoubles[1] - 2]);
                    for ii in cdoubles[2:end]
                        Cs[ii,iz,t] = Cs[ii-1,iz,t] + m*(kgrid[ii] - kgrid[ii-1]);
                    end
                end
                #computing assets k_{t+1}
                Ks[:,iz,t] .= (1+r(K))*kgrid .+ (1-τ)*w(K)*hbar*exp.(Es[t].+zs[iz]) .- Cs[:,iz,t];
                iconstrained = (Ks[:,iz,t] .< Kmin);
                if sum(iconstrained)>0
                    Cs[iconstrained,iz,t] = (1 + r(K)).*kgrid[iconstrained] .+ (1-τ)*w(K)*hbar*exp.(Es[t]*zs[iz])  .-kgrid[1];
                    Ks[iconstrained,iz,t] .= kgrid[1];
                end
                Vs[:,iz,t] = u(Cs[:,iz,t],η) .+ (β*Ss[t])*evalgrid(spl, Ks[:,iz,t], zs)*M[iz,:];
            end
        end=#
    end
end

#=

            end
        end
    end



   if t>T
        k_prev = (kgrid .+ c_prev .- b(K))/(1+r(K));
    else t<=T
        k_prev = (kgrid .+ c_prev .- w(K)*hbar*(1-τ(K))*exp.(Es[t].+repeat(zs',nk,1)))/(1+r(K));
    end
    for iz=1:nz
        kprev_low_int, kprev_high_int, w_low_int = interpolate_coord(k_prev[:,iz],kgrid,xqi,xqia,xqpi);
        Cs[:,iz,t] = c_prev[kprev_low_int,iz].*w_low_int .+ c_prev[kprev_high_int,iz].*(1 .- w_low_int);

        #binding borrowing constraint
        iconstained = (kgrid .<= k_prev[1,iz]); #points which are at the borrowing constraint
        if t>T && sum(iconstained)>0
            Cs[iconstained,iz,t] = (1 + r(K))*kgrid[iconstained] .+ b(K) .-kgrid[1];
        elseif t<=T && sum(iconstained)>0
            Cs[iconstained,iz,t] = (1 + r(K))*kgrid[iconstained] .+ w(K)*hbar*(1-τ(K))*exp(Es[t]*zs[iz]) .- kgrid[1];
        end
        #itp = Spline1D(k_prev[:,iz],c_prev[:,iz],k=1,bc="extrapolate");
        #Cs[.!iconstained,iz,t] = itp(kgrid[.!iconstained]);
    end

    ##computing asset policy corresponding to consumption policy
    As = zeros(size(Cs));
        #retirement
    As[:,:,T+1:TT] = (1 + r(K)).*kgrid .+ b(K) .- Cs[:,:,T+1:TT];
        #working life
    ass_t(t) = (1 + r(K)).*kgrid .+ w(K)*hbar*(1-τ(K))*exp.(Es[t].+repeat(zs',nk,1)) .- Cs[:,:,t];
    for t=1:T
        As[:,:,t] = ass_t(t);
    end
    return Cs, As;
end

=#
##-----     DISTRIBUTION FORWARD     -----##
"""
    distribution_forward(kgrid,K,zs,M,T,TR,TT)

    Solving household problem via Endogenous Grid Method.

   Inputs:
   -------
   ``kgrid``    :       array (na),  grid for assets\n
   ``K``        :       Number,  aggregate capital\n
   ``zs``       :       array (nZ),  aggregate productivities\n
   ``M``        :       array (nZ,nZ),   transition of individual productivity\n
   ``T``        :       Number,  years of working life\n
   ``TR``       :       Number,  years of retirement\n
   ``TT``       :       Number,  total life = T+TR\n
   -------
   ```
"""
function distribution_forward(kgrid,K,zs,y0_mass,M,T,TR,TT)
    k0 = 0;     #capital in first period of life
    nk = length(kgrid);
    nz = length(zs);

    #initializing capital distribution
    K_dist = zeros(nk,nz,TT);
    #structures for interpolation
    xqi = zeros(nk);
    xqia = zeros(nk);
    xqpi = zeros(nk);

        #period t=1: agents have k=0 and the initial distribution over productivity z
    K_dist[1,:,1] .= y0_mass[:];
    for t=2:TT
        for iz=1:nz
            kl_int, kh_int, wl_int = interpolate_coord(kgrid, As[:,iz,t-1], xqi, xqia ,xqpi);
            for ik=1:nk
                K_dist[kl_int[ik],:,t] .= K_dist[kl_int[ik],:,t] .+ wl_int[kl_int[ik]].*K_dist[ik,iz,t-1].*M[iz,:]
                K_dist[kh_int[ik],:,t] .= K_dist[kh_int[ik],:,t] .+ (1 - wl_int[kh_int[ik]]).*K_dist[ik,iz,t-1].*M[iz,:]
            end
        end
    end

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

#=
# Plotting policy functions from the household problem
function plot_policies(agrid,Kgrid,cpol)
    #plotting consumption policy function for low and high capital and employed and unemployed
    close()
    figure()
    #low capital
    subplot(2,1,1)
    plot(agrid,cpol[1,:,1],label=L"$\epsilon = e$")
    plot(agrid,cpol[2,:,1],label=L"$\epsilon = u$")
    title("K = "*string(Kgrid[1]))
    legend()
    xticks([])
    #high capital
    subplot(2,1,2)
    plot(agrid,cpol[1,:,end],label=L"$\epsilon = e$")
    plot(agrid,cpol[2,:,end],label=L"$\epsilon = u$")
    title("K = "*string(Kgrid[end]))
    subplots_adjust(hspace=0.5)
    xlabel("assets")
end


##-----     SIMULATION of STATES     -----##
"""
    simul_states(Nsimul)

    Simulate aggregate states (good = 1, bad = 2).

   Inputs:
   -------
   ``Nsimul``    :       Integer,   number of simulations\n
   -------
"""
function simul_states(Nsimul)
    dist_aggregate = Int.(zeros(Nsimul));
    dist_aggregate[1] = 1;
    rnd_nums = rand(Nsimul);
    tmpTmat = cumsum(ΓZ,dims=2);
    for i=2:Nsimul
        dist_aggregate[i] = findfirst(rnd_nums[i] .<= tmpTmat[dist_aggregate[i-1],:]);
    end
    return dist_aggregate;
end

"""
    simul_employment(NN,Nsimul,dist_aggregate)

    Simulate idiosyncratic (i.e. employment) states (employed = 0, unemployed = 1).

   Inputs:
   -------
   ``NN``        :       Integer,   number of agents\n
   ``Nsimul``    :       Integer,   number of simulations\n
   ``dist_aggregate``    :       array,     simulated aggregate states\n
   -------
"""
function simul_employment(NN,Nsimul,dist_aggregate)
    #employment distribution over time and workers
    dist_employment = Int.(zeros(NN,Nsimul));
    dist_employment[1:Int(floor(ug*NN)),1] .= 1;    #initilization for employed/unemployed
    rnd_nums = rand(NN,Nsimul);
    #creating the overall
    tmpTmats = zeros(2,2,4);
    tmpTmats[:,:,1] = cumsum(Γgg,dims=2);
    tmpTmats[:,:,2] = cumsum(Γbb,dims=2);
    tmpTmats[:,:,3] = cumsum(Γgb,dims=2);
    tmpTmats[:,:,4] = cumsum(Γbg,dims=2);
    for i=2:Nsimul
        if dist_aggregate[i-1]==1 && dist_aggregate[i]==1
            tmpTmat = tmpTmats[:,:,1];
        elseif dist_aggregate[i-1]==2 && dist_aggregate[i]==2
            tmpTmat = tmpTmats[:,:,2];
        elseif dist_aggregate[i-1]==1 && dist_aggregate[i]==2
            tmpTmat = tmpTmats[:,:,3];
        elseif dist_aggregate[i-1]==2 && dist_aggregate[i]==1
            tmpTmat = tmpTmats[:,:,4];
        end

        tmats = tmpTmat[1 .+ dist_employment[:,i-1],:];
        dist_employment[:,i] = map(x->findfirst(rnd_nums[x,i].<tmats[x,:]),1:NN) .- 1;
    end
    return dist_employment;
end


##-----     SIMULATE ASSETS     -----##
"""
    simulate_HH_assets(NN,Nsimul,astart,dist_aggregate,dist_employment,apol,Kgrid,agrid)

    Simulate individual agent's distribution of assets and the corresponding
    aggregate capital path.

   Inputs:
   -------
   ``NN``               :       Integer,   number of agents\n
   ``Nsimul``           :       Integer,   number of simulations\n
   ``astart``           :       Number,    initial asset holdings of each agent\n
   ``dist_aggregate``   :       array,     simulated aggregate states\n
   ``dist_employment``  :       array,     simulated idiosyncratic states\n
   ``aold``             :       array (2*nZ,na,nk),  asset policy fnct\n
   ``Kgrid``            :       array (nk),  grid for capital\n
   ``agrid``            :       array (na),  grid for assets\n
   -------
"""
function simulate_HH_assets(NN,Nsimul,astart,dist_aggregate,dist_employment,apol,Kgrid,agrid)
    #simulating using the HH policies to get aggregate capital in each period
    Kpath = zeros(Nsimul);
    asset_sim = zeros(NN,Nsimul);
    xfill1 = Int.(zeros(NN));
    xfill2 = Int.(zeros(NN));
    xfill3 = zeros(NN);
    asset_sim[:,1] .= astart;


    for i=2:Nsimul
        Klast = mean(asset_sim[:,i-1]);
        Kpath[i-1] = Klast;
        locK, wK = interpolate1D(Klast,Kgrid);
        #locas, was = interpolate2D(asset_sim[:,i-1],agrid);
        locas, locasH, was = interpolate_coord(agrid,asset_sim[:,i-1],xfill1,xfill2,xfill3);

        x_indxs = 2*(dist_aggregate[i-1]-1) + 1;

        #getting the individual components
        #low asset low K
        A_alow_Klow = dist_employment[:,i-1].*apol[x_indxs+1,locas,locK] .+ (1 .- dist_employment[:,i-1]).*apol[x_indxs,locas,locK];

        #high assets and low capital
        A_ahigh_Klow = dist_employment[:,i-1].*apol[x_indxs+1,locasH,locK] .+ (1 .- dist_employment[:,i-1]).*apol[x_indxs,locasH,locK];

        #low assets and high capital
        A_alow_Khigh = dist_employment[:,i-1].*apol[x_indxs+1,locas,locK+1] .+ (1 .- dist_employment[:,i-1]).*apol[x_indxs,locas,locK+1];

        #high assets and high capital
        A_ahigh_Khigh = dist_employment[:,i-1].*apol[x_indxs+1,locasH,locK+1] .+ (1 .- dist_employment[:,i-1]).*apol[x_indxs,locasH,locK+1];

        asset_sim[:,i] = (was*wK).*A_alow_Klow .+ (was*(1-wK)).*A_alow_Khigh .+ ((1 .- was)*wK).*A_ahigh_Klow .+ ((1 .- was)*(1-wK)).*A_ahigh_Khigh;
    end
    return asset_sim, Kpath;
end


##-----     ESTIMATE Law of Motion for CAPITAL     -----##
"""
    estimate_LOM(agrid,Kgrid,Hmat,Zs,dist_aggregate,dist_idiosyncratic,Ngarbage,cpolguess,tol=1e-3,maxiter=20)

    Estimate coefficient of capital law of motion by iteration.

   Inputs:
   -------
   ``agrid``             :       array (na),  grid for assets\n
   ``Kgrid``             :       array (nk),  grid for capital\n
   ``Hmat``              :       array (nZ,nZ),   law of motion parameters\n
   ``Zs``                :       array (nZ),  aggregate productivities\n
   ``dist_aggregate``    :       array,     simulated aggregate states\n
   ``dist_idiosyncratic``:       array,     simulated idiosyncratic states\n
   ``Ngarbage``          :       Integer,   number of period to discard\n
   ``cpolguess``         :       array (2*nZ,na,nk),  guess of consumption policy fnct\n
       --- optional inputs:\n
    ``tol``      :       Number,  algorithm tolerance\n
    ``maxiter``  :       Integer,  maximum number of iterations\n
   -------
"""
function estimate_LOM(agrid,Kgrid,Hmat,Zs,dist_aggregate,dist_idiosyncratic,Ngarbage,cpolguess,tol=1e-3,maxiter=20)

    NN, Nsimul = size(dist_idiosyncratic);
    itnum = 0;
    dist = Inf;
    dist_aggregate_clean = dist_aggregate[Ngarbage:Nsimul-2];
        #good state
    ig = findall(dist_aggregate_clean.==1);
        #bad state
    ib = findall(dist_aggregate_clean.==2);

    while dist>tol && itnum<maxiter
        #policy function for household
        cpol, apol, tol_res, iter_res = solveHH(agrid,Kgrid,Zs,Hmat,cpolguess);
        cpolguess = copy(cpol);
        #simulation of distributions
        A_sim, Kpath = simulate_HH_assets(NN,Nsimul,5,dist_aggregate,dist_idiosyncratic,apol,Kgrid,agrid);
        Kpath_clean = Kpath[Ngarbage:Nsimul-1]; #selecting relevant periods
        #updating LOM of capital
        lom_params_good = [ones(length(ig)) log.(Kpath_clean[ig])]\log.(Kpath_clean[ig.+1]);
        lom_params_bad = [ones(length(ib)) log.(Kpath_clean[ib])]\log.(Kpath_clean[ib.+1]);
        Hmat_new = [lom_params_good'; lom_params_bad'];

        dist = maximum(abs.(Hmat_new .- Hmat));
        if (itnum%1)==0
            println("Iteration = "*string(itnum)*";      dist = "*string(dist));
        end
        itnum += 1;
        Hmat = 0.5*Hmat_new + 0.5*Hmat;
    end

    return Hmat, dist, itnum;
end

=#

##-----     TAUCHEN DISCRETIZATION     -----##
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
