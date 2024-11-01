% OUTPUT
% C = companion matrix

% INPUT
% A = row cell of VAR coefficients' matrices

function [C] = A2C(A)
    K = size(A,2); % number of lags
    M = size(A{1,1},1); % number of variables
    C = zeros(M*K); % empty companion matrix
    
    for k = 1:1:K
        C(1:M,(k-1)*M+1:k*M) = A{1,k};
    end
    C(M+1:end,1:M*(K-1))=eye(size(C(M+1:end,1:M*(K-1))));
end