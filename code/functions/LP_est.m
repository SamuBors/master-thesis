% OUTPUTS
% LP = IRF coefficients
% s = variances of the coefficients

% INPUTS
% X = data
% H = IRF horizons (row vector)
% lags = number of lags
% K_max = max number of lags

function [LP,s] = LP_est(X,H,K,met)
    
    M = size(X,2); 
    if max(H) >= size(X,1)
        error(['Maximum IRF horizon (' num2str(max(H)) ') is higher or equal than T (' num2str(size(X,1)) ')'])
    end

    if strcmp(met,'c')
        
        % LP estimation
        LP = zeros(M,M,size(H,2));
        s = zeros(M,M,size(H,2));
        ch = 0;
        for h = H
            ch = ch + 1;
            if h == 0
                LP(:,:,ch) = eye(M);
                s(:,:,ch) = zeros(M);
            else
                Y = X(K+h:end,:);
                W = lagmatrix(X,0:(K-1));
                W = W(K:end,:);
                W = W(1:size(Y,1),:);
                pi=pinv(W'*W)*W'*Y;
                LP(:,:,ch)=pi(1:M,:)';
                
                Yh = W*pi;
                xi = Y - Yh;
    
                Q = W'*W;

                for i = 1:M
                    u = xi(:,i);
                    s1 = (W.*u)'*(W.*u);
                    s2 = zeros(M*K);
                    aux = W.*u;
                    auxla = lagmatrix(W.*u,1:h);
                    for j = 1:h
                        aux1 = aux(j+1:end,:);
                        auxla1 = auxla(j+1:end,(j-1)*M*K+1:j*M*K);
                        s2 = s2 + (1 - j/(h+1))*(aux1'*auxla1 + auxla1'*aux1);
                    end 
                    V = diag(pinv(Q)*(s1+s2)*pinv(Q));
                    s(i,:,ch) = V(1:M)';
                end

            end
        end

    elseif strcmp(met,'la')
        K = K + 1;

        % LP estimation
        LP = zeros(M,M,size(H,2));
        s = zeros(M,M,size(H,2));
        ch = 0;
        for h = H
            ch = ch + 1;
            if h == 0
                LP(:,:,ch) = eye(M);
                s(:,:,ch) = zeros(M);
            else
                Y = X(K+h:end,:);
                W = lagmatrix(X,0:(K-1));
                W = W(K:end,:);
                W = W(1:size(Y,1),:);
                pi=pinv(W'*W)*W'*Y;
                LP(:,:,ch)=pi(1:M,:)';
        
                Yh = W*pi;
                xi = Y - Yh;
    
                WY = W(:,1:M);
                WX = W(:,M+1:end);
    
                u = WY-WX*pinv(WX'*WX)*WX'*WY;
    
                S = u'*u/size(u,1);
    
                for i = 1:M
                    s(i,:,ch) = diag((1/size(u,1))^2 * (pinv(S)*(u.*xi(:,i))'*(u.*xi(:,i))*pinv(S)));
                end
            end
        end
    end
end
