% OUTPUTS
% X = final series

% INPUTS
% T = length of final series
% A = cell of VAR parameters
% th = share of T with unitary variance (i.e., from 1 to T*th v(e)=1, then v(e)=4)

function [X] = VARH_sim(T,A,th)

    % paramters
    K = size(A,2); % number of lags
    M = size(A{1,1},1); % number of variables

    % innovations
    Th = round(T*th);
    e = randn(T,M);
    e(Th:end,:) = e(Th:end,:)*2;
    
    % VAR
    X = zeros(T+K,M); % empty matrix
    C = A2C(A); % companion matrix

    W1=reshape(flip(X(1:K,:)',1),M*K,1);
    
    for t=K+1:1:T+K
        W1 = C * W1 + [e(t-K,:)';zeros(M*(K-1),1)];
        X(t,:) = W1(1:M,:)';
    end
    
    X=X(K+1:end,:);

end
