##

## Transition function (Z,ϵ) where Z is aggregate productivity and ϵ is
## employment status. For details on how to construct the transition matrix see
## Heer-Maussner Ch. 8.3

ug = 0.0386;   #probability of unemployment in good state
ub = 0.1073;    #probability of unemployment in bad state
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


function solveHH(agrid,Kgrid,Zs,Hmat)
    na = length(agrid);
    nK = length(Kgrid);
    nZ = length(Zs);

    #initiliazing matrices
    cguess = repeat(0.5*agrid,1,nZ*2,nK);
    cguess = permutedims(cguess,[2,1,3]);

    for ik = 1:nk
        k = Kgrid[ik];
        for iz = 1:nZ
            z = Zs[iz];
            knew = exp(Hmat[iz,:]'*[1;log(k)]);
            klow_int, wint = interpolateK(knew,Kgrid);
            #consumption guess given k' instead of k
            cguess_int = wint*cguess[:,:,klow_int] .+ (1-wint)*cguess[:,:,klow_int+1];
            #computing aggregate
            rnew = r(knew,Zs);   #next period interest rate
            wnew = w(knew,Zs);   #next period interest rate
            τnew = τ(knew,Zs);   #next period interest rate
            #computing consumption from the Euler equation (z)
            ctoday = β*Γ[2*(iz-1)+1:2*iz,:]*repeat(repeat((1 .+ (1 .- τnew).*rnew),1,2)'[:],1,na).*(cguess_int.^(η));
            for ie = 0:1

            end
        end
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
