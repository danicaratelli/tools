function fcst = fcstOLSVAR(model,Par)

% MODELS

SMALL  =     [33 115 87];              % EMP, FFR, CPI
CEE    =     [33 115 113 87 73 77 78]; % EMP, FFR, CPI, PCOM, M2, NBR, TR
%Main aggregate variables
MEDIUM =     [33 115 113 2 3 6 20 25 51 109 125  129 87 72 73 77 78 83 93 104]; 
% Whole panel
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


p = Par.p;         % Number of lags in the VAR
Jwind = Par.Jwind; % number of the observations used each time for estimation in the rolling scheme
                   % number longer than the sample size will result in a recursive scheme             


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

VFit  = zeros(size(transf));  VFit(Var_fit)  = 1; %% indicator for Fit variables
VEval = zeros(size(transf)); VEval(Var_eval) = 1; %% indicator for evaluation variables

% choosing the required variables in the VAR
eval(['VarList =', model,';']); %% vector of series in the VAR model

DATA      = DATA(:,VarList);
transf    = transf(VarList);
seriesVAR = series(VarList);

% variables for the fit
VFit = VFit(VarList);
VFit = find(VFit==1);
seriesFit = seriesVAR(VFit);
% variables for the forecast evaluation
VEval = VEval(VarList);
VEval = find(VEval==1);
seriesEval = seriesVAR(VEval);



%--------------------------------------------------------------------------
% Transforming the data

X = zeros(size(DATA));
X(:,ismember(transf,[1 2]))   = DATA(:,ismember(transf,[1 2]));
X(:,ismember(transf,[4 5 6])) = log(DATA(:,ismember(transf,[4 5 6])))*100;

% indicator variable for non-stationary (1) and stationary (0) variables
iRW = ones(size(transf));     iRW(ismember(transf,[1,4])) = 0; 

%--------------------------------------------------------------------------

[TT,NN] = size(X);                   % The size of the panel


%==========================================================================
% Calculating the fit
%==========================================================================
% Setting the training sample
j0 = find((time(:,1)==start_train_y) & (time(:,2)==start_train_mq));
j  = find((time(:,1)==end_train_y) & (time(:,2)==end_train_mq));
x  = X(j0:j,:);                                      %  The training sample

% Estimate the OLS VAR to get the fit on a presample
[Xpr,fit] = varOLS(x,p,-1);


%% REMARK, take the average fit across all variables
% AVGfit = mean(fit,1);

% Take the average fit across Fit variables
fit = mean(fit(VFit,:),1);

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

    x = X(j0:j,:);        %% The available data at each time point of the evaluation exercise

    %% get the predictions from VAR (the original data are augmented by the prediction as additional rows)
    Xpr = varOLS(x,p,hor);

    for h = 1:hor                                       %% loop across horizon
        tru(j+h,:,h) = X(j+h,:)-X(j,:);                 %% number to be forecast: h month difference h steps ahead
        pred(j+h,:,h)= Xpr(end-hor+h,:)-Xpr(end-hor,:); %% forecast h steps ahead of the h mnths differences
        rw(j+h,:,h)  = mean(x(h+1:end,:)-x(1:end-h,:)); %% prediction for constant growth model
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
fcst.VarIdx   = VarList;

fcst.VarFit   = seriesFit;
fcst.VarEval  = seriesEval;

fcst.MSFE    = MSFEvar;
fcst.rMSFE   = rMSFE;
fcst.MSFEben = MSFErw;


fcst.fit    = fit;
fcst.lambda = inf;
