% OUTPUTS
% X = final series

% INPUTS
% T = length of final series
% A = cell of VAR parameters
% G = cell of MA parameters
% W = weight matrix of MA process (P x M)

function [X] = VARMA_sim(T,A,G,W)

    % dimension check
    if size(W,1) ~= size(G{1,1},1)
        error('sizes of W and G do not match')
    end

    if size(W,2) ~= size(A{1,1},1)
        error('sizes of W and A do not match')
    end

    % paramters
    K = size(A,2); % number of lags
    M = size(A{1,1},1); % number of variables
    Q = size(G,2)-1; % highest lag of the MA process
    P = size(W,1); % number of MA processes

    % VMA
    Gm = cell2mat(G);
    if rank(G{1,1} - eye(P)) > 0
        error('The first matrix of the VMA is not the identity')
    end

    u = randn(T+Q,P); % iid components
    U = lagmatrix(u,0:Q);
    U = U(Q+1:end,:);
    e = U*Gm'*W; 
    
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
