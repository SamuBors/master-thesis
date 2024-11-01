% OUTPUTS
% IRFe = 5-D array of IRF estimates
% LBe = 5-D array of asymptotic C.I. LB estimates
% UBe = 5-D array of asymptotic C.I. UB estimates
% LBb = 5-D array of recursive bootstrap C.I. LB estimates
% UBb = 5-D array of recursive bootstrap C.I. UB estimates
% LBw = 5-D array of wild bootstrap C.I. LB estimates
% UBw = 5-D array of wild bootstrap C.I. UB estimates

% INPUTS
% T = time series length
% G = cell of MA or GARCH parameters
% inn = innovations type
% H = vector of IRF horizons
% K = lag order of VAR
% B = bootstrap repetitions
% M = number of variables in the system
% rhos = persistence level vector (it is only needed that the length of the array is equal to the number of persistence level, nothing else)
% met = DGP (AR1 or VAR3)
% th = share of periods before the structural break in the innovations covariance matrix (only used if inn = 'het')


function [IRFe,LBe,UBe,LBb,UBb,LBw,UBw] = loop(T,G,inn,H,K,B,M,rhos,met,th)

    % initializing the results matrices
    IRFe = zeros(M,M,size(H,2),size(rhos,1),2);
    LBe = zeros(M,M,size(H,2),size(rhos,1),2);
    UBe = zeros(M,M,size(H,2),size(rhos,1),2);
    LBb = zeros(M,M,size(H,2),size(rhos,1),2);
    UBb = zeros(M,M,size(H,2),size(rhos,1),2);
    LBw = zeros(M,M,size(H,2),size(rhos,1),2);
    UBw = zeros(M,M,size(H,2),size(rhos,1),2);


    for cp = 1:size(rhos,1)

        % 1 = LP
        [IRFe(:,:,:,cp,1),LBe(:,:,:,cp,1),UBe(:,:,:,cp,1),LBb(:,:,:,cp,1),UBb(:,:,:,cp,1),LBw(:,:,:,cp,1),UBw(:,:,:,cp,1)] = LP_loop(coeff_loop(met,cp),T,G,inn,'c',H,K,B,M,th);
        % 2 = LPLA
        [IRFe(:,:,:,cp,2),LBe(:,:,:,cp,2),UBe(:,:,:,cp,2),LBb(:,:,:,cp,2),UBb(:,:,:,cp,2),LBw(:,:,:,cp,2),UBw(:,:,:,cp,2)] = LP_loop(coeff_loop(met,cp),T,G,inn,'la',H,K,B,M,th);
        
    end
end