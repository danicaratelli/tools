function fcst = fcstBFAVAR(Par,nFact,FIT)
OPTS.disp = 0;

% Variables in the VAR
SMALL  =     [33 115 87];             % EMP, FFR, CPI
% Variables to calculate the factors
LARGE  =     [33 115 113 1:32 34:70  109:112 114 116:132 87 72:86 88:95 96:108]; 

%==========================================================================
%  Parameters 
%==========================================================================
Var_fit = Par.Var_fit;     % variables on which you want to evaluate the fit

start_train_y  =  Par.start_train_y;  % starting dates of the presample
start_train_mq =  Par.start_train_mq; 
end_train_y    =  Par.end_train_y;    % end dates of the presample
end_train_mq   =  Par.end_train_mq;   

start_eval_y  = Par.start_eval_y;  % starting dates for the out-of-sample evaluation
start_eval_mq = Par.start_eval_mq;
end_eval_y    = Par.end_eval_y;    % end dates for the out-of-sample evaluation   
end_eval_mq   = Par.end_eval_mq;

Var_eval = Par.Var_eval; % Variables for frecast evaluation 
hor_eval = Par.hor_eval; % Forecast horizons for evaluation

KK = Par.KK;


p = Par.p;         % Number of lags in the VAR
Jwind = Par.Jwind; % number of the observations used each time for estimation in the rolling scheme
                   % number longer than the sample size will result in a recursive scheme             

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

%  Data for the estimations of the factors
DATA   = DATA(:,LARGE);
transf = transf(LARGE);

X(:,ismember(transf,1))     = DATA(2:end,ismember(transf,1));
X(:,ismember(transf,2))     = DATA(2:end,ismember(transf,2))-...
    DATA(1:end-1,ismember(transf,2));
X(:,ismember(transf,4))     = log(DATA(2:end,ismember(transf,4)))*100;
X(:,ismember(transf,[5 6])) = (log(DATA(2:end,ismember(transf,[5 6])))-...
    log(DATA(1:end-1,ismember(transf,[5 6]))))*100;


Xall = X; clear X;

% Data fo VAR
VFit = zeros(size(transf));  VFit(Var_fit)   = 1; %% indicator for Var_fit variables
VEval = zeros(size(transf)); VEval(Var_eval) = 1; %% indicator for evaluation variables

% choosing the required variables in the VAR
DATA   = A(2:end,2:end);
transf = A(1,2:end); 
DATA   = DATA(:,SMALL);
transf = transf(SMALL);
seriesVAR = series(SMALL);

% variables for the fit
VFit = VFit(SMALL);
VFit = find(VFit==1);
seriesFit = seriesVAR(VFit);
% variables for the forecast evaluation
VEval = VEval(SMALL);
VEval = find(VEval==1);
seriesEval = seriesVAR(VEval);

%--------------------------------------------------------------------------
% Transforming the data

X(:,ismember(transf,[1 2]))   = DATA(2:end,ismember(transf,[1 2]));
X(:,ismember(transf,[4 5 6])) = log(DATA(2:end,ismember(transf,[4 5 6])))*100;
time = time(2:end,:);

% indicator variable for non-stationary (1) and stationary (0) variables
iRW = ones(size(transf)); iRW(ismember(transf,[1,4])) = 0; 
iRW = [iRW zeros(1,nFact)];

%--------------------------------------------------------------------------

[TT,NN] = size(X);                   % The size of the panel

%==========================================================================
% Choosing the overall tightness
%==========================================================================
% overall tightness is governed by hyperparameter pi(which is the inverse
% of lambda from the paper)
% pi = 0: prior variance =infty -> OLS
% pi = infty: prior variance = 0, random walk model

% Setting the training sample
j0 = find((time(:,1)==start_train_y) & (time(:,2)==start_train_mq));
j  = find((time(:,1)==end_train_y)   & (time(:,2)==end_train_mq));

x   = X(j0:j,:);                                     %  The training sample
xpc = Xall(j0:j,:);

% computing fit for the grid of pi's
GRIDsearch = GRID*sqrt(NN*p);

for jpi = 1:length(GRIDsearch)
    pi = GRIDsearch(jpi);    % The overall tightness prior
    mu = pi*KK;              % Set the prior on the sum of the coefficients
    % Factors
    xstd = (xpc-ones(size(xpc(:,1)))*mean(xpc))*diag(1./std(xpc));
    [F,temp] = eigs(xstd*xstd',nFact,'lm',OPTS);
    % Estimate the Bayesian VAR to get the fit corresponding to the prior pi
    [Xpr,fit(:,jpi)] = bvarLitt([x F],p,pi,mu,-1,iRW);
    
end;

% REMARK, take the average fit across all variables 
% AVGfit = mean(fit,1);

%%Take the average fit across KEY variables
AVGfit = mean(fit(VFit,:),1); 

% search for pi yielding the right fit
[temp,Jstar] = min(abs(AVGfit-FIT));
FITstar      = AVGfit(Jstar);
pi           = GRIDsearch(Jstar);

% corresponding hyperparameter for the prior on sum of coefficients
% (mu is the inverse of tau from the paper)
mu = KK*pi;

%==========================================================================
% Out-of-sample exercise
%==========================================================================
hor = max(hor_eval); % Define the maximum horizon to evaluate the forecasts

% Indices for the start and end of the evaluation period
start_eval = find((time(:,1)==start_eval_y) & (time(:,2)==start_eval_mq));
end_eval   = find((time(:,1)==end_eval_y)   & (time(:,2)==end_eval_mq));


% all the variables below are of the size nobs x nvar x nhor
% for a given horizon h, non-zero entries are from 
% start_eval-hor+h till  end_eval-1+h
% the other entries are 0

tru  = zeros(TT,NN,hor);
pred = zeros(TT,NN,hor);
rw   = zeros(TT,NN,hor);

for j = start_eval-hor:end_eval-1

    if j-Jwind+1>0       %% If rolling scheme
        j0 = j-Jwind+1;
    else                 %% If recursive scheme
        j0 = 1;
    end;

    x = X(j0:j,:); %% The available data at the beginning of the evaluation period
    % Factors
    xpc  = Xall(j0:j,:);
    xstd = (xpc-ones(size(xpc(:,1)))*mean(xpc))*diag(1./std(xpc));
    [F,temp] = eigs(xstd*xstd',nFact,'lm',OPTS);
    %% Estimate the Bayesian VAR to get the forecasts
    Xpr = bvarLitt([x F],p,pi,mu,hor,iRW);
        
    for h = 1:hor                                             %% loop across horizon
        tru(j+h,:,h) = X(j+h,:)-X(j,:);                       %% number to be forecast: h month difference h steps ahead
        pred(j+h,:,h)= Xpr(end-hor+h,1:NN)-Xpr(end-hor,1:NN); %% forecast h steps ahead of the h mnths differences
        rw(j+h,:,h)  = mean(x(h+1:end,:)-x(1:end-h,:));       %% prediction for constant growth model, to be used for IP
    end;

end;

%==========================================================================
% MSFE for chosen variables
%==========================================================================

%the msfe will be given in the order as in VarList!!!!

%MSFE's
% From the model
MSFEvar = squeeze(mean((tru(start_eval:end_eval,VEval,hor_eval)...
    -pred(start_eval:end_eval,VEval,hor_eval)).^2,1));
% From the benchmark model
MSFErw  = squeeze(mean((tru(start_eval:end_eval,VEval,hor_eval)...
    -rw(start_eval:end_eval,VEval,hor_eval)).^2,1));

% From the model relative to the benchmark
rMSFE = MSFEvar./MSFErw;

%==========================================================================
% Loading the structure
%==========================================================================

fcst.fcst = pred;
fcst.tru  = tru;
fcst.ben  = rw;
fcst.time = time;

fcst.VarNames = seriesVAR;
fcst.VarIdx   = SMALL;

fcst.VarFit   = seriesFit;
fcst.VarEval  = seriesEval;

fcst.MSFE    = MSFEvar;
fcst.rMSFE   = rMSFE;
fcst.MSFEben = MSFErw;


fcst.fit    = FITstar;
fcst.lambda = 1./pi;
fcst.tau    = 1./(KK*pi);

fcst.nFact = nFact;
