##

## Transition function (Z,ϵ) where Z is aggregate productivity and ϵ is
## employment status. For details on how to construct the transition matrix see
## Heer-Maussner Ch. 8.3

ΓZ = [0.8 0.2; 0.2 0.8];    #transition matrix for aggregate state
Γgg = [0.9615 0.0385; 0.9581 0.0419];
Γbb = [0.9525 0.0475; 0.3952 0.6048];
Γgb = [(1-ub)/(1-ug) (ub-ug)/(1-ug); 0 1];
Γbg = [1 0; (ub-ug)/ub ug/ub];
Γ = [Γgg*0.8 Γgb*0.2; Γbg*0.2 Γbb*0.8];

## Grids
    #asset holdings
agrid = collect(range(amin,stop=amax,length=na));
Kgrid = collect(range(Kmin,stop=Kmax,length=nk));


## Household
    #useful functions
N(Z) = floor.(Z)*(1-ug) .+ (1 .- floor.(Z))*(1-ub);
u(Z) = 1 .- N(Z);
w(K,Z) = (1-α)*Z.*(K./N(Z)).^(α);
r(K,Z) = α*Z.*(N(Z)/K).^(1-α) .- δ;
τ(K,Z) = 1 ./((w(K,Z).*N(Z) .+ r(K,Z)*K)./(ζ*w(K,Z).*u(Z)) .+ 1);


function solveHH(agrid,Kgrid,Zs,Hmat,iV=1,tol=1e-6,maxiter=1000)
    na = length(agrid);
    nK = length(Kgrid);
    nZ = length(Zs);

    #initiliazing matrices
    cold = repeat(0.2*agrid,1,nZ*2,nK);
    cold = permutedims(cold,[2,1,3]) .+ 0.01;
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
        if (itnum%20)==0
            println("Iteration = "*string(itnum)*";      dist = "*string(dist));
        end
        itnum += 1;
        cold = copy(cnew);
    end

    if iV==1
        #computing the corresponding Value function
        Vold = (cold.^(1-η))/(1-η);
        Vnew = zeros(2*nZ,na,nK);

        itnum = 0;
        dist = Inf;
        while dist>tol && itnum<maxiter
            for ik=1:nK
                Vnew[:,:,ik] = (cpol[:,:,ik].^(1-η))/(1-η) .+ β*Γ*Vold[:,:,ik];
            end
            dist = maximum(abs.(Vnew .- Vold));
            Vold = copy(Vnew);
            itnum += 1;
        end
        return cold, Vold, dist, itnum
    else
        return cold, dist, itnum
    end
end


function interpolateK(k,Kgrid)
    nK = length(Kgrid);
    if k<Kgrid[1]
        kloc_min = 1;
        w = 1;
    elseif k>Kgrid[end]
        kloc_min = nK-1;
        w = 0;
    else
        kloc_min = findlast(k.>=Kgrid);
        w = 1 - (k - Kgrid[kloc_min])/(Kgrid[kloc_min + 1] - Kgrid[kloc_min]);
    end
    return kloc_min, w;
end
