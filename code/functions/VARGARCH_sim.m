% OUTPUTS
% X = final series

% INPUTS
% T = length of final series
% A = cell of VAR parameters
% G = cell of GARCH parameters ({i,1} = constant; {i,2} = row vector of ARCH params; {i,3} = row vector of GARCH params)
% W = weight matrix of GARCH process (P x M)

function [X] = VARGARCH_sim(T,A,G,W)

    % dimension check
    if size(W,2) ~= size(A{1,1},1)
        error('sizes of W and A do not match')
    end

    % paramters
    K = size(A,2); % number of lags
    M = size(A{1,1},1); % number of variables
    P = size(W,1); % number of GARCH processes

    % GARCH
    U=zeros(T,P);
    for i = 1:P
        U(:,i) = GARCH_sim(T,G);
    end
    e = U*W; 
    
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
