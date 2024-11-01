% OUTPUTS
% VAR = sturcutre containing
%   VAR.coeff.AR = cell of VAR estimates
%   VAR.res = OLS residuals
%   VAR.cov = covariance matrix of estimates

function [VAR] = VAR_est(X,K)
    
    M = size(X,2);    
    lags = lagmatrix(X,0:K);
    lags = lags(K+1:end,:);
    Y = lags(:,1:M);
    W = lags(:,M+1:end);
    pi=Y'*W*pinv(W'*W);

    for k = 1:K
        VAR.coef.AR{1,k} = pi(:,(k-1)*M+1:k*M);
    end

    VAR.res = Y - W*pi';
    
    VAR.cov = zeros(size(pi,1)*size(pi,2));
    for i = 1:M
        ui = VAR.res(:,i);
        Wu = W.*ui;
        VAR.cov((i-1)*size(pi,2)+1:i*size(pi,2),(i-1)*size(pi,2)+1:i*size(pi,2)) = pinv(W'*W) * (Wu'*Wu) * pinv(W'*W);
    end
end