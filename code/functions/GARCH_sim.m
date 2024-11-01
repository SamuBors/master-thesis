% OUTPUTS
% X = final series

% INPUTS
% T = length of final series
% G = cell of GARCH parameters ({1,1} = constant; {1,2} = row vector of ARCH params; {1,3} = row vector of GARCH params)

function X = GARCH_sim(T,G)

    % params
    o = G{1,1}; %constant
    a = G{1,2}; % ARCH params
    b = G{1,3}; % GARCH params
    p = size(G{1,2},2); % ARCH lags, p in GARCH(p,q)
    q = size(G{1,3},2);  % GARCH lags, q in GARCH(p,q)

    
    % initialization
    z = randn(T,1); % iid component
    s = ones(T+q,1); % sigma
    x = zeros(T+p,1); % GARCH

    % process
    for t = 1:T
        aux = o + a'*(x(t:t+p-1).^2) + b'*(s(t:t+q-1).^2);
        s(t+q) = sqrt(aux);
        x(t+p) = s(t+q)*z(t);
    end

    X = x(p+1:end);

end

% NOTES
% to generate an ARCH(p) set G{1,3} = 0