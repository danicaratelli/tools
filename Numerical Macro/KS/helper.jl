##

## Transition function (Z,ϵ) where Z is aggregate productivity and ϵ is
## employment status. For details on how to construct the transition matrix see
## Heer-Maussner Ch. 8.3

ug = 0.0386;    #probability of unemployment in good state
ub = 0.1073;    #probability of unemployment in bad state
us = [ug;ub];
ΓZ = [0.8 0.2; 0.2 0.8];    #transition matrix for aggregate state
Γgg = [0.9615 0.0385; 0.9581 0.0419];
Γbb = [0.9525 0.0475; 0.3952 0.6048];
Γgb = [(1-ub)/(1-ug) (ub-ug)/(1-ug); 0 1];
Γbg = [1 0; (ub-ug)/ub ug/ub];
Γ = [Γgg*0.8 Γgb*0.2; Γbg*0.2 Γbb*0.8];

# constructing grids
agrid = collect(range(amin,stop=amax,length=na));
Kgrid = collect(range(Kmin,stop=Kmax,length=nk));


#useful functions
N(Z) = floor.(Z)*(1-ug) .+ (1 .- floor.(Z))*(1-ub);
u(Z) = 1 .- N(Z);
w(K,Z) = (1-α)*Z.*(K./N(Z)).^(α);
r(K,Z) = α*Z.*(N(Z)/K).^(1-α) .- δ;
τ(K,Z) = 1 ./((w(K,Z).*N(Z) .+ r(K,Z)*K)./(ζ*w(K,Z).*u(Z)) .+ 1);


# endogenous grid method for HH problem
function solveHH_assets(agrid,Kgrid,Zs,Hmat,iV=1,cold=0,iprint=0,tol=1e-6,maxiter=1000)
    na = length(agrid);
    nK = length(Kgrid);
    nZ = length(Zs);

    #initiliazing matrices
    if cold==0
        cold = repeat(0.2*agrid,1,nZ*2,nK);
        cold = permutedims(cold,[2,1,3]) .+ 0.01;
    end
    cnew = zeros(nZ*2,na,nk);

    itnum = 0;
    dist = Inf;

    while dist>tol && itnum<maxiter
        for ik = 1:nk
            k = Kgrid[ik];
            for iz = 1:nZ
                z = Zs[iz];
                knew = exp(Hmat[iz,:]'*[1;log(k)]);
                klow_int, wint = interpolateK(knew,Kgrid);
                #consumption guess given k' instead of k
                cold_int = wint*cold[:,:,klow_int] .+ (1-wint)*cold[:,:,klow_int+1];
                #computing aggregate
                rtoday = r(k,z);   #today's interest rate
                wtoday = w(k,z);   #today's interest rate
                τtoday = τ(k,z);   #today's interest rate
                rnew = r(knew,Zs);   #next period interest rate
                wnew = w(knew,Zs);   #next period interest rate
                τnew = τ(knew,Zs);   #next period interest rate
                #computing consumption from the Euler equation (z)
                ctoday = (β*Γ[2*(iz-1)+1:2*iz,:]*((repeat((1 .+ (1 .- τnew).*rnew),1,2)'[:]).*(cold_int.^(-η)))).^(-1/η);
                #looping over employed (ie=0) and unemployed (ie=1)
                for ie = 0:1
                    if ie==0
                        atoday = (agrid .+ ctoday[ie+1,:] .- (1-τtoday)*wtoday)./(1 + (1 - τtoday)*rtoday);
                    elseif ie==1
                        atoday = (agrid .+ ctoday[ie+1,:] .- ζ*(1-τtoday)*wtoday)./(1 + (1 - τtoday)*rtoday);
                    end
                    #getting consumption today over agrid rather than atoday
                    ctoday_new = zeros(na);
                    #binding borrowing constraint
                    iconstained = agrid .<= atoday[1]; #points which are at the borrowing constraint
                    if ie==0 && sum(iconstained)>0
                        ctoday_new[iconstained] = (1 + (1-τtoday)*rtoday)*agrid[iconstained] .+ (1-τtoday)*wtoday .- agrid[1];
                    elseif ie==1 && sum(iconstained)>0
                        ctoday_new[iconstained] = (1 + (1-τtoday)*rtoday)*agrid[iconstained] .+ ζ*(1-τtoday)*wtoday .- agrid[1];
                    end
                    #interpolating today's assets
                    itp = Spline1D(atoday,ctoday[ie+1,:],k=1,bc="extrapolate")
                    #itp = interpolate((atoday,), ctoday[ie+1,:], Gridded(Linear()));
                    ctoday_new[.!iconstained] = itp(agrid[.!iconstained]);
                    #inserting updated ctoday_new in cnew
                    cnew[2*(iz-1)+1 + ie,:,ik] = ctoday_new;
                end
            end
        end
        dist = maximum(abs.(cnew .- cold));
        if (iprint==1) && (itnum%20)==0
            println("Iteration = "*string(itnum)*";      dist = "*string(dist));
        end
        itnum += 1;
        cold = copy(cnew);
    end

    ##computing asset policy corresponding to consumption policy
    apol = zeros(size(cold));
    for ik=1:nk
        for iz = 1:nZ
            apol[2*(iz-1)+1,:,ik] = (1 + (1-τ(Kgrid[ik],Zs[iz]))*r(Kgrid[ik],Zs[iz]))*agrid .+ (1-τ(Kgrid[ik],Zs[iz]))*w(Kgrid[ik],Zs[iz]) .- cold[2*(iz-1)+1,:,ik];
            apol[2*(iz-1)+2,:,ik] = (1 + (1-τ(Kgrid[ik],Zs[iz]))*r(Kgrid[ik],Zs[iz]))*agrid .+ ζ*(1-τ(Kgrid[ik],Zs[iz]))*w(Kgrid[ik],Zs[iz]) .- cold[2*(iz-1)+2,:,ik];
        end
    end

    return cold, apol, dist, itnum;
end


# linear interpolation routine (for singleton or vector)
function interpolateK(k,Kgrid)
    nK = length(Kgrid);
    if length(k)==1
            nK = length(Kgrid);
        if k<Kgrid[1]
            kloc_min = 1;
            wk = 1;
        elseif k>Kgrid[end]
            kloc_min = nK-1;
            wk = 0;
        else
            kloc_min = findlast(k.>=Kgrid);
            wk = 1 - (k - Kgrid[kloc_min])/(Kgrid[kloc_min + 1] - Kgrid[kloc_min]);
        end
    else
        kloc_min = Int.(zeros(length(k)));
        wk = zeros(length(k));
            #below minimum
        kloc_min[k.<Kgrid[1]] .= 1;
        wk[k.<Kgrid[1]] .= 1;
        kloc_min[k.>Kgrid[end]] .= nk - 1;
        wk[k.>Kgrid[1]] .= 0;
        kintra = (k.>=Kgrid[1]) .* (k.<=Kgrid[end]);
        kloc_min[kintra] = map(x->findlast(x.>=Kgrid),k[kintra]);
        wk[kintra] = 1 .- (k[kintra] .- Kgrid[kloc_min[kintra]])./(Kgrid[kloc_min[kintra] .+ 1] .- Kgrid[kloc_min[kintra]]);
    end
    return kloc_min, wk;
end


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


## Simuation of aggregate states
function simul_states(Nsimul)
    simulZ = Int.(zeros(Nsimul));
    simulZ[1] = 1;
    rnd_nums = rand(Nsimul);
    tmpTmat = cumsum(ΓZ,dims=2);
    for i=2:Nsimul
        simulZ[i] = findfirst(rnd_nums[i] .<= tmpTmat[simulZ[i-1],:]);
    end
    return simulZ;
end

function simul_employment(NN,Nsimul,simulZ)
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
        if simulZ[i-1]==1 && simulZ[i]==1
            tmpTmat = tmpTmats[:,:,1];
        elseif simulZ[i-1]==2 && simulZ[i]==2
            tmpTmat = tmpTmats[:,:,2];
        elseif simulZ[i-1]==1 && simulZ[i]==2
            tmpTmat = tmpTmats[:,:,3];
        elseif simulZ[i-1]==2 && simulZ[i]==1
            tmpTmat = tmpTmats[:,:,4];
        end

        tmats = tmpTmat[1 .+ dist_employment[:,i-1],:];
        dist_employment[:,i] = map(x->findfirst(rnd_nums[x,i].<tmats[x,:]),1:NN) .- 1;
    end
    return dist_employment;
end


function simulate_HH_assets(NN,Nsimul,astart,simulZ,dist_employment,apol,Kgrid,agrid)
    #simulating using the HH policies to get aggregate capital in each period
    Kpath = zeros(Nsimul);
    asset_sim = zeros(NN,Nsimul);
    asset_sim[:,1] .= astart;


    for i=2:Nsimul
        Klast = mean(asset_sim[:,i-1]);
        Kpath[i-1] = Klast;
        locK, wK = interpolateK(Klast,Kgrid);
        locas, was = interpolateK(asset_sim[:,i-1],agrid);
        ind_extract = zeros(Nsimul)
        x_indxs = 2*(simulZ[i-1]-1) .+ 1 .+ dist_employment[:,i-1];

        #getting the individual components
        #low asset low K
        A_alow_Klow = dist_employment[:,i-1].*apol[2*(simulZ[i-1]-1)+2,locas,locK] .+ (1 .- dist_employment[:,i-1]).*apol[2*(simulZ[i-1]-1)+1,locas,locK];

        #high assets and low capital
        A_ahigh_Klow = dist_employment[:,i-1].*apol[2*(simulZ[i-1]-1)+2,locas.+1,locK] .+ (1 .- dist_employment[:,i-1]).*apol[2*(simulZ[i-1]-1)+1,locas.+1,locK];

        #low assets and high capital
        A_alow_Khigh = dist_employment[:,i-1].*apol[2*(simulZ[i-1]-1)+2,locas,locK+1] .+ (1 .- dist_employment[:,i-1]).*apol[2*(simulZ[i-1]-1)+1,locas,locK+1];

        #high assets and high capital
        A_ahigh_Khigh = dist_employment[:,i-1].*apol[2*(simulZ[i-1]-1)+2,locas.+ 1,locK+1] .+ (1 .- dist_employment[:,i-1]).*apol[2*(simulZ[i-1]-1)+1,locas.+ 1,locK+1];

        asset_sim[:,i] = (was*wK).*A_alow_Klow .+ (was*(1-wK)).*A_alow_Khigh .+ ((1 .- was)*wK).*A_ahigh_Klow .+ ((1 .- was)*(1-wK)).*A_ahigh_Khigh;
    end
    return asset_sim, Kpath;
end


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
        cpol, apol, tol_res, iter_res = solveHH_assets(agrid,Kgrid,Zs,Hmat,0,cpolguess);
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
