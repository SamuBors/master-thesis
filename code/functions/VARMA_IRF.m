% OUTPUTS
% IRF = VARMA IRFs (the first element of the 3rd dimension correspond to
% h=0)

% INPUTS
% A = cell of VAR parameters
% G = cell of MA parameters
% Ho = IRF horizons

function [IRF] = VARMA_IRF(A,G,Ho)
    
    M = size(A{1,1},1);
    P = size(A,2);
    Q = size(G,2)-1;
    
    A11 = A2C(A);
    if Q > 0
        A12 = zeros(M*P,M*Q);
        aux = cell2mat(G);
        A12(1:M,:) = aux(:,M+1:end);
        A21 = zeros(M*Q,M*P);
        aux=cell(1,Q);
        for i =1:Q
            aux{1,i} = zeros(M);
        end
        A22 = A2C(aux);
        
        B = [A11,A12;A21,A22];
    else
        B = A11;
    end
    
    J = zeros(M,M*(P+Q));
    J(1:M,1:M) = eye(M);
    
    H1 = zeros(M*P,M);
    H1(1:M,1:M) = eye(M);
    if Q > 0
        H2 = zeros(M*Q,M);
        H2(1:M,1:M) = eye(M);
        H = [H1;H2];
    else
        H = H1;
    end
    
    IRF = zeros(M,M,size(Ho,2));
    ch = 0;
    for h = Ho
        ch = ch+1;
        IRF(:,:,ch)=J*B^h*H;
    end
end

