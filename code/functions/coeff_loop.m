% OUTPUT
% A = row cell of VAR coefficients' matrices

% INPUT
% met = DGP type (AR1 or VAR3)
% cp = persistence scenario (from 1 to 4)

function [A] = coeff_loop(met,cp)
    if strcmp(met,'AR1')
        % AR(1) params
        if cp == 1
            A={0};
        elseif cp == 2
            A={.5};
        elseif cp == 3
            A={.95};
        elseif cp == 4
            A={1};
        end
    elseif strcmp(met,'VAR3')
        % VAR(3) params (M=2)
        A = cell(1,3);
        if cp == 1
            A{1,1}=zeros(2);
            A{1,2}=zeros(2);
            A{1,3}=zeros(2);
        elseif cp == 2
            A{1,1}=[0.3,0.1;-0.4,0.1];
            A{1,2}=-[0.02,-0.05;0.3,-0.3];
            A{1,3}=[0.1,0;0,0.04];
        elseif cp == 3
            A{1,1}=[0.3,0.1065;-0.4,0.7];
            A{1,2}=-[0.3,-0.7;0.5,-0.1];
            A{1,3}=[0.1,0;0,0.3];
        elseif cp == 4
            A{1,1}=eye(2);
            A{1,2}=zeros(2);
            A{1,3}=zeros(2);
        end
    end
end