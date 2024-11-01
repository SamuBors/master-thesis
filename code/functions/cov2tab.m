% OUTPUS
% tab = coverage tables

% INPUTS
% cov = logical array of coverage (1 if true value is inside, 0 otherwise)
% i = responding variable
% j = impulse variable
% rhos = array of persistence levels
% H = array of IRF horizons

function [tab] = cov2tab(cov,i,j,rhos,H)
    bm = mean(cov,6);
    bmm = squeeze(bm(i,j,:,:));
    tab = zeros(size(bmm)+1);
    tab(2:end,2:end) = bmm;
    tab(2:end,1) = H';
    tab(1,2:end) = rhos';
end