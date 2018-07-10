%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           This is the script to replicate the results from              %
%                "Bayesian VARs with Large Panels"                        %
%                by Banbura, Giannone and Reichlin                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The data are from Stock and Watson (2005)

clear
clc

%==========================================================================
%                      P A R A M E T E R S                                %
%==========================================================================

% THE ONLY VARYING PARAMETER IS KK

% KK is the tightness of the prior on the sum of coeffcient with respect to the Litterman prior
% In the paper: tau = lambda/KK, if KK=0 there is no prior on the sum of coefficients

% KK IS DIFFERENT FOR DIFFERENT TABLES:

P.KK = 0;   % Baseline prior (Tables 1-3, A.1, Figure B.1)
% P.KK = 0.1; % Baseline prior+prior on the sum of coefficients (Tables 4-5, A.2, B.1-B.4, Figures 1-2)
 

%--------------------------------------------------------------------------

% REMAINING PARAMTERS ARE THE SAME FOR ALL THE TABLES IN THE PAPER

%--------------------------------------------------------------------------
% TIGHTNESS

P.Var_fit = [33 115 87];    % variables on which you want to evaluate the fit
                            % EMP, FFR, CPI

% Training sample for finding the tightness
P.start_train_y = 1960; P.start_train_mq = 2;  % starting dates of the presample
P.end_train_y   = 1970;   P.end_train_mq = 1;  % end dates of the presample

%--------------------------------------------------------------------------
% FORECAST EVALUATION

P.start_eval_y = 1971; P.start_eval_mq = 1;   % starting dates for the out-of-sample evaluation
P.end_eval_y   = 2003;   P.end_eval_mq = 1;   % end dates for the out-of-sample evaluation

% Variables of interest for the forecast evaluation (for those MSFE will be
% calculated)
P.Var_eval = [33 115 87];         % EMP, FFR, CPI

P.hor_eval = [1 3 6 12];          % Forecast horizons for evaluation

%--------------------------------------------------------------------------
% IMPULSE RESPONSE FUNCTION

P.start_irf_y = 1961; P.start_irf_mq = 1;   % starting dates for the computation of the IRF
P.end_irf_y   = 2002;   P.end_irf_mq = 12;  % end dates for the computation of the IRF

P.Var_MP = 87;                    %index of the monetary policy variable

% Key variables for Variance Decomposition (e.g. Table 4)
P.Var_VD = [33 115 87];           % EMP, FFR, CPI

P.lag    = [0 3 6 12 24 36 48];   % lags for the responses
P.hor_VD = [1 3 6 12 24 36 48];   % forecast horizons for the variance decomposition

P.Gibbs = 200;                    % number of draws for the confidence bands

% Variables for the plots for the impulse responses for all the models
P.Var_plotAll = [33 115 87];      % EMP, FFR, CPI

P.lag_plot = 48;                  % maximal lag for plotting the IRFs

% Confidence levels for the impulse responses
P.level1 = .90;
P.level2 = .68;


%--------------------------------------------------------------------------
% OTHER PARAMETERS

P.p = 13;      % Number of lags in the VAR
P.Jwind = 120; % number of the observations used each time for estimation in the rolling scheme
               % number longer than the sample size will result in a recursive scheme             

%--------------------------------------------------------------------------
% REMARK
%
% As the computations for the large model require long time and large
% memory the lines corresponding to this model were commented out.
% (They can be restored on a good computer)

 

%==========================================================================
%                         T A B L E S                                     %
%==========================================================================

%--------------------------------------------------------------------------
% Bayesian VAR with fits from the OLS for SMALL
% Tables 1, 4
%--------------------------------------------------------------------------
fcstS   = fcstOLSVAR('SMALL',P);

fcstC   = fcstBVAR('SMALL',  P,fcstS.fit);
fcstM   = fcstBVAR('MEDIUM',P,fcstS.fit);
% fcstL   = fcstBVAR('LARGE', P,fcstS.fit);

Table1 = [fcstS.rMSFE(:),fcstC.rMSFE(:),fcstM.rMSFE(:);...
    fcstS.lambda,fcstC.lambda,fcstM.lambda];
% Table1 = [fcstS.rMSFE(:),fcstC.rMSFE(:),fcstM.rMSFE(:),fcstL.rMSFE(:);...
%     fcstS.lambda,fcstC.lambda,fcstM.lambda,fcstL.lambda];

%--------------------------------------------------------------------------
% OLS VAR with p = 13, BIC, BVAR for SMALL and CEE 
% Table 2 
%--------------------------------------------------------------------------

fcstSols = fcstOLSVAR('SMALL',P);
fcstSbic = fcstOLSVARBIC('SMALL',P);


fcstCols = fcstOLSVAR('CEE',P);
fcstCbic = fcstOLSVARBIC('CEE',P);

Table2   = [fcstSols.rMSFE(:),fcstSbic.rMSFE(:),fcstS.rMSFE(:),...
    fcstCols.rMSFE(:),fcstCbic.rMSFE(:),fcstC.rMSFE(:)];
% 
%--------------------------------------------------------------------------
% FAVAR with p = 13, BIC, BVAR  
% Table 3
%--------------------------------------------------------------------------
% one factor
fcstFols1 = fcstOLSFAVAR(P,1);
fcstFbic1 = fcstOLSFAVARBIC(P,1);
fcstF1    = fcstBFAVAR(P,1,fcstS.fit);
% three factors
fcstFols3 = fcstOLSFAVAR(P,3);
fcstFbic3 = fcstOLSFAVARBIC(P,3);
fcstF3    = fcstBFAVAR(P,3,fcstS.fit);

Table3 = [fcstFols1.rMSFE(:),fcstFbic1.rMSFE(:),fcstF1.rMSFE(:),...
    fcstFols3.rMSFE(:),fcstFbic3.rMSFE(:),fcstF3.rMSFE(:)];


%--------------------------------------------------------------------------
% Bayesian VAR with different fits 
% Table A.1, A.2
%--------------------------------------------------------------------------
fits = [0.25 0.5 0.75];
nFits = length(fits);
TableA1 = zeros((length(P.hor_eval)*length(P.Var_eval)+1)*nFits,3);
% TableB1 = zeros((length(P.hor_eval)*length(P.Var_eval)+1)*nFits,4)
for i = 1:length(fits)
    fcstSf(i)   = fcstBVAR('SMALL', P,fits(i));
    fcstCf(i)   = fcstBVAR('CEE',   P,fits(i));
    fcstMf(i)   = fcstBVAR('MEDIUM',P,fits(i));
%     fcstLf(i)   = fcstBVAR('LARGE', P,fits(i));

    TableA1(i:nFits:end,:) = [fcstSf(i).rMSFE(:),fcstCf(i).rMSFE(:),fcstMf(i).rMSFE(:);...
        fcstSf(i).lambda,fcstCf(i).lambda,fcstMf(i).lambda];
%     TableB1(i:nFits:end,:) = [fcstSf(i).rMSFE(:),fcstCf(i).rMSFE(:),fcstMf(i).rMSFE(:),fcstLf(i).rMSFE(:);...
%         fcstSf(i).lambda,fcstCf(i).lambda,fcstMf(i).lambda,fcstLf(i).lambda];
end

%--------------------------------------------------------------------------
% Variance Decomposition
% Table 5
%--------------------------------------------------------------------------
irfS = irfBVAR('SMALL', P,fcstS.fit);
irfC = irfBVAR('CEE',   P,fcstS.fit);
irfM = irfBVAR('MEDIUM',P,fcstS.fit);
% irfL = irfBVAR('LARGE', P,fcstS.fit);

Table5 = [irfS.VDkey(:),irfC.VDkey(:),irfM.VDkey(:)];
% Table5 = [irfS.VDkey(:),irfC.VDkey(:),irfM.VDkey(:),irfL.VDkey(:)];

%--------------------------------------------------------------------------
% Impulse Response Function + Variance decomposition
% Tables B.1-B.4
%--------------------------------------------------------------------------

TableS = [irfS.IRF_lag,irfS.VD_all];
TableC = [irfC.IRF_lag,irfC.VD_all];
TableM = [irfM.IRF_lag,irfM.VD_all];
% TableL = [irfL.IRF_lag,irfL.VD_all];

%==========================================================================
%                           G R A P H S
%==========================================================================
% Figure 1, B1
P.Models= {'SMALL','CEE','MEDIUM','LARGE'}
IRFplotAllmod1(P,irfS,irfC,irfM);
% IRFplotAllmod(P,irfS,irfC,irfM,irfL);

% Figure 2
IRFplot1mod1(P,irfM);
% IRFplot1mod(P,irfL);

