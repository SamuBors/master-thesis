% OUTPUTS
% X = B bootstrapped series

% INPUTS
% X = original series (for starting values)
% A = cell of VAR estimates
% res = innovations (to be resampled)
% method = bootstrap method (string)
% B = bootstrap iterations

function [Xb] = VAR_boot(X,A,res,method,B)

    K = size(A,2); % number of lags
    M = size(A{1,1},1); % number of variables
    T = size(res,1)+K;

    % centering the residuals
    resc = res - mean(res);
    
    if strcmp(method,'ri') % ri = recursive iid bootstrap
        boot=datasample(resc,T*B);
        e = zeros(T,M*B);
        for b = 1:B
            e(:,M*(b-1)+1:M*b) = boot(T*(b-1)+1:T*b,:);
        end

        % DGP
        Xb = zeros(T+K,M*B); % empty matrix
        Xb(1:K,:) = repmat(X(1:K,:),1,B); % setting the starting values
        C = kron(eye(B),A2C(A)); % companion matrix
        
        % if max(abs(eig(A2C(A)))) > 1
        %     disp(['Non-stationary VAR, max eigenvalue = ' num2str(max(abs(eig(A2C(A)))))])
        % elseif max(abs(eig(A2C(A)))) == 1
        %     ur=sum(abs(eig(A2C(A)))==1);
        %     disp(['Presence of ' num2str(ur) ' unit roots'])
        % end
        
        W=reshape(flip(Xb(1:K,:)',1),M*K*B,1);
        
        for t=K+1:1:T+K
            W = C * W + [e(t-K,:)';zeros(M*B*(K-1),1)];
            Xb(t,:) = W(1:M*B,:)';
        end
        
        Xb=Xb(K+1:end,:);
        
    elseif strcmp(method,'w') % w = wild bootstrap
        aux = repmat(resc,1,B);
        S = kron(eye(B),ones(M));
        w = randn(size(aux))*sqrt(S)'/sqrt(2);
        e = aux .* w;

        % DGP
        Xb = zeros(T+K,M*B); % empty matrix
        Xb(1:K,:) = repmat(X(1:K,:),1,B); % setting the starting values
        C = kron(eye(B),A2C(A)); % companion matrix
        
        % if max(abs(eig(A2C(A)))) > 1
        %     disp(['Non-stationary VAR, max eigenvalue = ' num2str(max(abs(eig(A2C(A)))))])
        % elseif max(abs(eig(A2C(A)))) == 1
        %     ur=sum(abs(eig(A2C(A)))==1);
        %     disp(['Presence of ' num2str(ur) ' unit roots'])
        % end
        
        W=reshape(flip(Xb(1:K,:)',1),M*K*B,1);
        
        for t=K+1:1:T
            W = C * W + [e(t-K,:)';zeros(M*B*(K-1),1)];
            Xb(t,:) = W(1:M*B,:)';
        end

    end
    
    
end
