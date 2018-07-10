function irf = irfBVAR(model,Par,FIT);
OPTS.disp = 0;
OPTS.tol = 1e-2;
OPTS.disp = 0;

% MODELS

% We assume a recursive identification scheme, therefore
% slow variables should be listed before the interest rate which should be 
% followed by the fast variables  

SMALL  =     [33 115 87];  %  EMP, FFR, CPI
CEE    =     [33 115 113 87 73 77 78]; % EMP, FFR, CPI, PCOM, M2, NBR, TR
%Main aggregate variables
MEDIUM =     [33 115 113 2 3 6 20 25 51 109 125  129 87 72 73 77 78 83 93 104]; 
% Whole panel
LARGE  =     [33 115 113 1:32 34:70  109:112 114 116:132 87 72:86 88:95 96:108]; 


%==========================================================================
%  Parameters 
%==========================================================================

Var_fit = Par.Var_fit;     % variables on which you want to evaluate the fit

start_train_y  =  Par.start_train_y; % starting dates of the presample
start_train_mq =  Par.start_train_mq; 
end_train_y    =  Par.end_train_y;    % end dates of the presample
end_train_mq   =  Par.end_train_mq;   

Var_VD = Par.Var_VD; % Key variables for Variance Decomposition
Var_MP = Par.Var_MP; % index of the monetary policy variable

start_irf_y  = Par.start_irf_y; % starting dates for the computation of the IRF 
start_irf_mq = Par.start_irf_mq;
end_irf_y    = Par.end_irf_y;   % end dates for the computation of the IRF   
end_irf_mq   = Par.end_irf_mq;

lagIRF = Par.lag;      % lags for the responses
hor_VD = Par.hor_VD;   % forecast horizons for the variance decomposition

KK = Par.KK;

Gibbs = Par.Gibbs; % number of draws for the confidence bands

p = Par.p;   % Number of lags in the VAR

                     
GRID = [0:.025:5 50]; % Grid to look for the tightness hyperparamter


%==========================================================================
% Preparing the panel
%==========================================================================
%% Load the data from hof.xls (Stock and Watson (2005))
[A,B] = xlsread('hof.xls');
matlabDates = datenum('30-Dec-1899') + A(2:end,1);
dates = matlabDates(1:end,1);
time = datevec(dates); 
time = time(:,1:2);

DATA   = A(2:end,2:end);
transf = A(1,2:end); %% Vector of transformations
series = B(1,2:end); %% The mnemonic for the variables in the panel

VFit = zeros(size(transf)); VFit(Var_fit) = 1; %% indicator for Var_fit variables
VMP  = zeros(size(transf)); VMP(Var_MP)   = 1; %% indicator for monetary policy variable
VVD  = zeros(size(transf)); VVD(Var_VD)   = 1; %% indicator for variance decomposition variables

% choosing the required variables in the VAR
eval(['VarList =', model,';']); %% vector of series in the VAR model

DATA      = DATA(:,VarList);
transf    = transf(VarList);
seriesVAR = series(VarList);

% variables for the fit
VFit = VFit(VarList);
VFit = find(VFit==1);
seriesFit = seriesVAR(VFit);

% Monetary policy variable
VMP = VMP(VarList);
VMP = find(VMP==1);
seriesMP = seriesVAR(VMP);

% Variables of interest for Variance Decomposition
VVD = VVD(VarList);
VVD = find(VVD==1);
seriesVD = seriesVAR(VVD);


%--------------------------------------------------------------------------
% Transforming the data

X = zeros(size(DATA));
X(:,ismember(transf,[1 2]))   = DATA(:,ismember(transf,[1 2]));
X(:,ismember(transf,[4 5 6])) = log(DATA(:,ismember(transf,[4 5 6])))*100;

% indicator variable for non-stationary (1) and stationary (0) variables
iRW = ones(size(transf)); iRW(ismember(transf,[1,4])) = 0; 

%--------------------------------------------------------------------------
N = length(iRW);

%==========================================================================
% Choosing the overall tightness
%==========================================================================
% overall tightness is governed by hyperparameter pi (which is the inverse
% of lambda from the paper)

% pi = 0: prior variance =infty -> OLS
% pi = infty: prior variance = 0, random walk model

% Setting the training sample
j0 = find((time(:,1)==start_train_y) & (time(:,2)==start_train_mq));
j  = find((time(:,1)==end_train_y)   & (time(:,2)==end_train_mq));
x  = X(j0:j,:);                                      %  The training sample

% computing fit for the grid of pi's
GRIDsearch = GRID*sqrt(N*p);

for jpi = 1:length(GRIDsearch)
    pi = GRIDsearch(jpi);   % The overall tightness prior
    mu = pi*KK;             % Set the prior on the sum of the coefficients
    % Estimate the Bayesian VAR to get the fit corresponding to the prior pi
    [Xpr,fit(:,jpi)] = bvarLitt(x,p,pi,mu,-1,iRW);
end;

% REMARK, take the average fit across all variables 
% AVGfit = mean(fit,1);

% Take the average fit across Var_fit variables
AVGfit = mean(fit(VFit,:),1); 

% search for pi yielding the right fit
[temp,Jstar] = min(abs(AVGfit-FIT));
FITstar = AVGfit(Jstar);
pi      = GRIDsearch(Jstar);

% corresponding hyperparameter for the prior on sum of coefficients
% (mu is the inverse of tau from the paper)
mu = KK*pi;

%==========================================================================
% MAtrix of regressors
%==========================================================================
% Indices for the start and end of the evaluation period
start_irf = find((time(:,1)==start_irf_y) & (time(:,2)==start_irf_mq));
end_irf   = find((time(:,1)==end_irf_y)   & (time(:,2)==end_irf_mq));


XI = X(start_irf:end_irf,:);          % The data for the irf
[T,N] = size(XI);                     % The size of the panel

Z = [];
for jp = 1:p
    Z = [Z XI(p+1-jp:end-jp,:)];
end
Z = [Z ones(size(Z(:,1)))];
Y = XI(p+1:end,:);

for jn = 1:N
    [war,Aar,Car] = arfit(XI(:,jn),p,p);
    SS0(jn) = sqrt(Car);
end;
MM0 = mean(XI);


%==========================================================================
% Dummies
%==========================================================================
%%% Construct dummy for litterman prior
Yrw = pi*[diag(SS0.*iRW);zeros(N*(p-1),N)];
Zrw = pi*[diag(kron(1:p,SS0)) zeros(N*p,1)];

%%% Construct dummy for the constant
Ycs = 1e-5*[zeros(1,N)];
Zcs = 1e-5*[zeros(1,N*p) 1];

%%% Construct dummy for the sum of the coefficients
Ylr = mu*diag(MM0.*iRW);
Zlr = mu*[kron(ones(1,p),diag(MM0.*iRW)) zeros(N,1)];

%% Construct dummies for prior on covariance matrix of residual;
Ycv = diag(SS0);
Zcv = zeros(N,N*p+1);

%put together all the information
Ypr = [Yrw;Ylr; Ycv; Ycs]; 
Zpr = [Zrw; Zlr; Zcv; Zcs];

Tpr = size(Ypr,1);

%==========================================================================
% Posterior
%==========================================================================

ZZinv = inv(Zpr'*Zpr + Z'*Z);
ZY    = Zpr'*Ypr + Z'*Y;
YY    = Ypr'*Ypr + Y'*Y;

beta = ZZinv*ZY;

%==========================================================================
% IRF 
%==========================================================================

e  = [Y; Ypr]-[Z; Zpr]*beta;
Su = 1/(T+Tpr)*e'*e;

k = size(beta,1);

AA = zeros(N*p);
AA(1:N,:) = beta(1:end-1,:)';
AA(N+1:end,1:N*(p-1)) = eye(N*(p-1));

CC = zeros(N*p,N); 
CC(1:N,1:N) = chol(Su)';
JJ = [eye(N) zeros(N,N*(p-1))];
AAj = eye(N*p);

lag = max([lagIRF,hor_VD]);
for j = 0:lag; 
    IRF(:,:,j+1) = JJ*AAj*CC*diag(1./diag(CC));
    AAj = AAj*AA;
end;

%==========================================================================
% Confidence bands for the IRF
%==========================================================================

CSuInv = chol(inv(Su*(T+Tpr)));
CZZinv = chol(ZZinv);

jg = 1;

while jg <= Gibbs
    
    Z    = randn(T+Tpr-k+2,N)*CSuInv;
    SuG  = inv(Z'*Z);
    CsuG = chol(SuG);

    temp  = randn(size(beta));
    betaG = beta + CZZinv'*temp*CsuG;

    AAG = zeros(N*p);
    AAG(1:N,:) = betaG(1:end-1,:)';
    AAG(N+1:end,1:N*(p-1)) = eye(N*(p-1));

    CCG = zeros(N*p,N); 
    CCG(1:N,1:N) = CsuG';
    JJ = [eye(N) zeros(N,N*(p-1))];

    AAj = eye(N*p);
    for j = 0:lag; IRFG(:,:,j+1,jg) = JJ*AAj*CCG*diag(1./diag(CCG));AAj = AAj*AAG;end;

    jg = jg+1;
    
end;

%==========================================================================
% Variance decomposition
%==========================================================================

D = diag(diag(Su));
for h = 1:lag+1
    V_MP(:,h)  = diag(squeeze(IRF(:,VMP,h))*D(VMP,VMP)*squeeze(IRF(:,VMP,h))');
    V_ALL(:,h) =  diag(squeeze(IRF(:,:,h))*D*squeeze(IRF(:,:,h))');
end;
VarDecomp = cumsum(V_MP,2)./cumsum(V_ALL,2)*100;

%==========================================================================
% Loading the structure
%==========================================================================

irf.IRF  = IRF; clear IRF;
irf.IRFG = IRFG; clear IRFG;
irf.VarNames = seriesVAR;
irf.VarIdx   = VarList;

irf.VarFit = seriesFit;

irf.VarMP  = seriesMP;
irf.VDkey  = VarDecomp(VVD,hor_VD)';
irf.VD_all = VarDecomp(:,hor_VD);
irf.IRF_lag = squeeze(irf.IRF(:,VMP,lagIRF+1));

irf.fit     = FITstar;
fcst.lambda = 1./pi;
fcst.tau    = 1./(KK*pi);



