## Author: Daniele Caratelli
#email: danicaratelli@gmail.com (don't hesitate to contact if you are
#confused, at worst I will not answer you)

#description: We estimate a VAR for " yₜ = c + Ψ₁yₜ₋₁ + ... + Ψₚyₜ₋ₚ + ϵ "
#             where ϵ ∼ N(0,Σ)
#             We spit out Impulse Response functions for a 1σ shock to a variable



#------------- Bureau ---------------
workspace();

## Path setting##
code_path = pwd();
fun_path =  joinpath(code_path,"functions");
cd("..");
root = pwd();
out_path = joinpath(root,"Outdata");
data_path = joinpath(root,"Data");

## Packages used ##
using DataFrames, PyPlot, Dates, Distributions, FredData

cd(code_path);
include("options.jl");
include(joinpath(fun_path,"Data_organizer.jl"));
include(joinpath(fun_path,"var_reg.jl"));
include(joinpath(fun_path,"forecast.jl"));
include(joinpath(fun_path,"RMSE.jl"));
include(joinpath(fun_path,"imp_res.jl"));
include(joinpath(fun_path,"dummy_fores.jl"));
include(joinpath(fun_path,"transformations.jl"));
#------------- Bureau ---------------





#------------- Loading Data ---------------
include(joinpath(data_path,"data_build.jl"));

#constructing date vector
dates = Data[Symbol(dates_mnems)]
dates = eval(parse("dates = Data[:"*dates_mnems*"]"));
#st_ = find(dates.==Date(start_date))[1]; en_ = find(dates.==Date(end_date))[1];
st_ = 1; en_ = length(dates);
dates = dates[st_:en_];


Y_stacked = Data_organizer(Data,Var_names, st_, en_, Transformations, size(dates,1),n); #prepare the data

if !isequal(intersect(Transformations,["lin","ln"]),Transformations);
    y_stacked = zeros(Y_stacked[:,2:end]);
     dates = dates[2:end];
else
    y_stacked = zeros(Y_stacked);
end

#------------- Loading Data ---------------






#------------- Running VAR ---------------
T = size(dates,1); #number of periods

var_res = var_reg(Y_stacked,p); #running regression

β = var_res["coeff"]; Σ = var_res["var"]; Fits = var_res["fits"]; Y = var_res["dep"];

#Constructing the coefficient matrix in companion form:
C = [β[1,:]' zeros(1,n*(p-1))]'; #constant term
β¹ = β[2:end,:];
A = zeros(n*p,n*p);
A[1:n,:] = β¹';
A[n+1:end,1:n*(p-1)] = eye(n*(p-1));

#Selection matrix
S = zeros(n,n*p);
S[:,1:n] = eye(n,n);
#------------- Running VAR ---------------


#------------- Impulse Response Functions ---------------


YMAT = zeros(T,n,ndraws);

for k=1:ndraws;
    YMAT[:,:,k] = dummy_fores(Y_stacked[:,1:p],A,C,Σ,p,T,S)'; #dummy forecast with sampled errors
end

IRFs = zeros(nperiods,n,n,ndraws);

B = zeros(size(β,1),n,ndraws);
V = zeros(n,n,ndraws);

for kk=1:ndraws;
    #running the VAR estiation
    var_res = var_reg(YMAT[:,:,kk]',p); #running regression
    β_temp = var_res["coeff"]; Var_temp = var_res["var"];
    B[:,:,kk] = β_temp; V[:,:,kk] = Var_temp;
    IRFs[:,:,:,kk] = imp_res(β_temp,Var_temp,S,nperiods,1);
    IRFs[:,:,:,kk] = IRFs[:,:,:,kk];
end

for s = 1:n;
    for r = 1:n;
        #percentiles of the computed IRFs
        IRF_conf = zeros(nperiods,3); #the standard errors of the bootstrapping

        for k=1:nperiods;
            set = reshape(IRFs[k,r,s,:],ndraws);
            set_sort = sort(set,1);
            IRF_conf[k,1] = set_sort[Integer(((50-conf_int/2)/100)*ndraws)];
            IRF_conf[k,2] = set_sort[Integer((50/100)*ndraws)];
            IRF_conf[k,3] = set_sort[Integer(((50+conf_int/2)/100)*ndraws)];
        end
        fig = figure(figsize=(15,10));
        fill_between(0:nperiods-1,IRF_conf[:,3],IRF_conf[:,1],color="gray",alpha = 0.4,
                     label=string(conf_int)*" % confidence bands");
        p_imp = plot(0:nperiods-1,IRF_conf[:,2],color="black",linewidth=2, label="median IRF");
        xlabel("Horizon",fontsize=15);
        legend();
        PyPlot.title("Response of "*Var_names[r]*" to "*Var_names[s]*" shock",fontsize=17);
        savefig(joinpath(out_path,irf_fold,"IRF_r"*Vars[s]*"_s"*Vars[r]*".png"))
        close("all");
    end
end
#------------- Impulse Response Functions ---------------
