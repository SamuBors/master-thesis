% inputs
% C = structure containing:
%   C.T = vector of time series lengths
%   C.inn = vertical cell of innovations types (iid, het, garch, miss)
%   C.inv = vertical cell of innovations invertibility: y = invertible, n =
%   non invertible, empty = non MA innovations
%   C.met = vertical cell of DGP (AR1 or VAR3)
%   C.I = number of simulations in each run
%   C.B = number of bootstrap repetitions in each run

function [] = sim(C)
    
    % number of methods
    RUN = size(C.T,1);

    for run = 1:RUN

        % inputs
        T = C.T(run);
        inn = C.inn{run};
        inv = C.inv{run};
        met = C.met{run};
        I = C.I(run);
        B = C.B(run);
        
        % VAR DGP
        if T >= 100
            H = [1,6,12,36,60]; % Horizons for tables
        elseif T == 25
            H = [1,5,10,15,20]; 
        else
            error('T is neither 25 nor greater than 100, H?')
        end
        
        if strcmp(met,'AR1')
            K = 1;
            M = 1;
        elseif strcmp(met,'VAR3')
            K = 3;
            M = 2;
        end % select K and M
        
        if strcmp(inn,'het')
            th = 1/2;
        else
            th = 1;
        end % select th = 0.5 if heteroskedastic innovations (if not choose 1, yet th is not used in the other cases)
        
        if strcmp(inn,'miss')
            n = 1/3;
        else
            n = inf;
        end % select n
        
        G{1,1} = eye(M);
        if strcmp(inn,'miss')
            if strcmp(inv,'y') % invertible MA coeff
                if strcmp(met,'AR1')
                    G{1,2} = -1;
                    G{1,3} = 0.5;
                    G{1,4} = -0.5;
                    G{1,5} = -0.7;
                    G{1,6} = 0.6;
                    G{1,7} = 0.3;
                    G{1,8} = 0.4;
                    G{1,9} = 0.8;
                    G{1,10} = -0.5;
                elseif strcmp(met,'VAR3')
                    G{1,2} = [0.4,.5;-1,.7];
                    G{1,3} = [.3,-.3;-0.9,1];
                    G{1,4} = [-.5,.1;-0.7,0.8];
                    G{1,5} = [1.2,0;-0.8,0.6];
                    G{1,6} = -[.4,-.3;-0.3,0.5];
                    G{1,7} = [0,.5;1,0];
                    G{1,8} = -[0.1,-0.2;0.2,0.1];
                    G{1,9} = [.1,-.1;.2,0.1];
                    G{1,10} = [.9,.4;-.4,.3];
                end
            else % non-invertible MA coeff
                if strcmp(met,'AR1')
                    G{1,2} = -1;
                    G{1,3} = 1.5;
                    G{1,4} = -0.8;
                    G{1,5} = -1;
                    G{1,6} = 1.2;
                    G{1,7} = 2;
                    G{1,8} = 1.2;
                    G{1,9} = 0.8;
                    G{1,10} = -1;
                elseif strcmp(met,'VAR3')
                    G{1,2} = [1.5,.5;-1,.7];
                    G{1,3} = [.3,-.3;-1,1];
                    G{1,4} = [-.5,.1;-0.7,0.8];
                    G{1,5} = [1.2,0;-1,0.6];
                    G{1,6} = -[.4,-.3;-0.3,0.5];
                    G{1,7} = [1,.5;1,1.5];
                    G{1,8} = -[0.1,-0.2;0.2,0.2];
                    G{1,9} = [.1,-.1;.2,0.1];
                    G{1,10} = [1,1.5;-.4,.3];
                end
            end
            for q = 2:size(G,2)
                G{1,q} = G{1,q}*T^(-n);
            end
        elseif strcmp(inn,'garch')
            G = cell(1,3);
            G{1,1} = 0.1;
            G{1,2} = 0.3;
            G{1,3} = 0.7;
        end
        
        rhos = zeros(1,1);
        tpm = cell(1,size(rhos,1));
        for cp = 1:4
            A = coeff_loop(met,cp);
        
            % creating rhos for tables
            rhos(cp,1) = max(abs(eig(A2C(A))));
            
            % true IRF cell
            if strcmp(inn,'miss')
                G1 = G;
            else
                G1 = {eye(M)};
            end
            tpm{1,cp} = VARMA_IRF(A,G1,H);
        end
        
        % coverage VAR (parallel)
        tic
       
        IRFe = zeros(M,M,size(H,2),size(rhos,1),2,I);
        LBe = zeros(M,M,size(H,2),size(rhos,1),2,I);
        UBe = zeros(M,M,size(H,2),size(rhos,1),2,I);
        LBb = zeros(M,M,size(H,2),size(rhos,1),2,I);
        UBb = zeros(M,M,size(H,2),size(rhos,1),2,I);
        LBw = zeros(M,M,size(H,2),size(rhos,1),2,I);
        UBw = zeros(M,M,size(H,2),size(rhos,1),2,I);
        
        %parfor i = 1:I
        for i = 1:I
            disp(['i=' num2str(i)])
            [IRFe(:,:,:,:,:,i),LBe(:,:,:,:,:,i),UBe(:,:,:,:,:,i),LBb(:,:,:,:,:,i),UBb(:,:,:,:,:,i),LBw(:,:,:,:,:,i),UBw(:,:,:,:,:,i)] = loop(T,G,inn,H,K,B,M,rhos,met,th);
        end
        
        toc

        % mat files
        date = datetime('today');
        if strcmp(inv,'n')
            inn = 'miss_ninv';
        end
        name =  [met '_' inn '_T' num2str(T) '_I' num2str(I) '_B' num2str(B) '_' num2str(day(date)) '_' num2str(month(date))];
        save(['outputs/mat/' name '.mat'],"IRFe","LBe","UBe","LBb","UBb","LBw","UBw");
        disp('Doubles saved as mat')
    end
end