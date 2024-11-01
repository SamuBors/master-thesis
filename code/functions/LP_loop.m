% OUTPUTS
% LPe = LP estimates
% LBe = asymptotic C.I. LB estimates
% UBe = asymptotic C.I. UB estimates
% LBb = recursive bootstrap C.I. LB estimates
% UBb = recursive bootstrap C.I. UB estimates
% LBw = wild bootstrap C.I. LB estimates
% UBw = wild bootstrap C.I. UB estimates

% INPUTS
% A = VAR parameters cell
% T = time series length
% inn = innovations type
% lp = LP method ('c' = classical, 'la' = lag-augmented)
% H = IRF horizons vector
% K = VAR lag-order
% B = bootstrap replicates
% M = number of variables in the system
% th = share of periods before the structural break in the innovations covariance matrix (only used if inn = 'het')

function [LPe,LBe,UBe,LBb,UBb,LBw,UBw] = LP_loop(A,T,G,inn,lp,H,K,B,M,th)

    % DGP selection
    if strcmp(inn,'iid')
        X = VARMA_sim(T,A,{eye(M)},eye(M));
    elseif strcmp(inn,'miss')
        X = VARMA_sim(T,A,G,eye(M));
    elseif strcmp(inn,'het')
        X = VARH_sim(T,A,th);
    elseif strcmp(inn,'garch')
        X = VARGARCH_sim(T,A,G,eye(M));
    end

    if size(A{1,1},1) > 1 && mean(A{1,1} == eye(M),'all') == 1
        K = 1;
    end
    
    % LP estimation
    [LPe,s] = LP_est(X,H,K,lp);

    % asymptotic CI
    UBe = LPe+s.^0.5*norminv(.95); 
    LBe = LPe+s.^0.5*norminv(.05);
    
    % VAR estimation
    VAR = VAR_est(X,K);
    res = VAR.res;
    Ah = VAR.coef.AR;
    Ch = A2C(Ah);

    % RECURSIVE BOOTSTRAP

    % initializing the bootstrap matrices
    LPTb = zeros(M,M,size(H,2));
    Rb = zeros(M,M,size(H,2),B);
    
    % bootstrap true IRF
    for ch = 1:size(H,2)
        h = H(ch);
        aux = Ch^h;
        LPTb(:,:,ch) = aux(1:M,1:M);
    end
    
    % bootstrapping procedure
    Xb = VAR_boot(X,Ah,res,'ri',B);
    for b = 1:B
        Xbb = Xb(:,M*(b-1)+1:M*b);
        [LPb,sb] = LP_est(Xbb,H,K,lp);
        Rb(:,:,:,b)=(LPb-LPTb)./sb.^0.5;
    end
    
    % bootstrap CI
    LBR = prctile(Rb,5,4);
    UBR = prctile(Rb,95,4);
    LBb = LBR .* s.^0.5 + LPe;
    UBb = UBR .* s.^0.5 + LPe;

    % WILD BOOTSTRAP

    % initializing the bootstrap matrices (LPTb is the same as above)
    Rw = zeros(M,M,size(H,2),B);

    % bootstrap true IRF
    Xbw = VAR_boot(X,Ah,res,'w',B);
    for b = 1:B
        Xbbw = Xbw(:,M*(b-1)+1:M*b);
        [LPw,sw] = LP_est(Xbbw,H,K,lp);
        Rw(:,:,:,b)=(LPw-LPTb)./sw.^0.5;
    end
    
    % bootstrap CI
    LBR = prctile(Rw,5,4);
    UBR = prctile(Rw,95,4);    
    LBw = LBR .* s.^0.5 + LPe;
    UBw = UBR .* s.^0.5 + LPe;
    
end