function [X,fit] = varOLS(X,p,H);


[T,N] = size(X);

% Regressors
Z = [];

for jp = 1:p
    Z = [Z X(p+1-jp:end-jp,:)];
end

Z = [Z ones(size(Z(:,1)))];

Y = X(p+1:end,:);


% Estimate of the parameters
ZZinv = inv(Z'*Z);
ZY    = Z'*Y;

beta = ZZinv*ZY;

% Fit
if nargout>1
    SS = 1/(T-p+1)*(Y-Z*beta)'*(Y-Z*beta);
   fit = 1 - diag(SS)./diag(cov(diff(X)));
end;

% Predictions
pred = [];
for j = 0:H
    
    X = [X;pred];
    Z = [];
    for jp = 0:p-1
        Z = [Z X(end-jp,:)];
    end
    
    Z = [Z 1];
    pred = Z(end,:)*beta;
    
end;
