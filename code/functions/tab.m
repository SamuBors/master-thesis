% INPUTS
% LBe = 5-D array of asymptotic C.I. LB estimates
% UBe = 5-D array of asymptotic C.I. UB estimates
% LBb = 5-D array of recursive bootstrap C.I. LB estimates
% UBb = 5-D array of recursive bootstrap C.I. UB estimates
% LBw = 5-D array of wild bootstrap C.I. LB estimates
% UBw = 5-D array of wild bootstrap C.I. UB estimates
% C = 
% i = responding variable
% j = impulse variable
% wirf used only when inn = 'miss' ('y'= w.r.t. VARMA IRF, 'n' = w.r.t. VAR IRF)

function [] = tab(LBe,UBe,LBb,UBb,LBw,UBw,C,i,j,wirf)
    
    % inputs
    T = C.T;
    inn = C.inn;
    inv = C.inv;
    met = C.met;
    I = C.I;
    B = C.B;
    
    % VAR DGP
    if T >= 100
        H = [1,6,12,36,60]; % Horizons for tables
    elseif T == 25
        H = [1,5,10,15,20]; 
    else
        error('T is neither 25 nor greater than 100, H?')
    end
    
    if strcmp(met,'AR1')
        M = 1;
    elseif strcmp(met,'VAR3')
        M = 2;
    end % select K and M
    
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
        if strcmp(wirf,'n') && strcmp(inn,'miss')
            G1 = G;
        else
            G1 = {eye(M)};
        end
        tpm{1,cp} = VARMA_IRF(A,G1,H);
    end

    % from ci to tables arrays
    covb = zeros(M,M,size(H,2),size(rhos,1),2,I);
    cove = zeros(M,M,size(H,2),size(rhos,1),2,I);
    covw = zeros(M,M,size(H,2),size(rhos,1),2,I);
    
    lene = UBe - LBe;
    lenb = UBb - LBb;
    lenw = UBw - LBw;
    
    % storing results for coverage
    for cp = 1:4
        for ch = 1:5
            tp = tpm{1,cp}(:,:,ch);
            covb(:,:,ch,cp,:,:) = LBb(:,:,ch,cp,:,:) <= tp & UBb(:,:,ch,cp,:,:) >= tp;
            cove(:,:,ch,cp,:,:) = LBe(:,:,ch,cp,:,:) <= tp & UBe(:,:,ch,cp,:,:) >= tp;
            covw(:,:,ch,cp,:,:) = LBw(:,:,ch,cp,:,:) <= tp & UBw(:,:,ch,cp,:,:) >= tp;
        end
    end

    % table
    tab.b_cov_c = cov2tab(covb(:,:,:,:,1,:),i,j,rhos,H);
    tab.e_cov_c = cov2tab(cove(:,:,:,:,1,:),i,j,rhos,H);
    tab.w_cov_c = cov2tab(covw(:,:,:,:,1,:),i,j,rhos,H);
    tab.b_len_c = cov2tab(lenb(:,:,:,:,1,:),i,j,rhos,H);
    tab.e_len_c = cov2tab(lene(:,:,:,:,1,:),i,j,rhos,H);
    tab.w_len_c = cov2tab(lenw(:,:,:,:,1,:),i,j,rhos,H);
    
    tab.b_cov_la = cov2tab(covb(:,:,:,:,2,:),i,j,rhos,H);
    tab.e_cov_la = cov2tab(cove(:,:,:,:,2,:),i,j,rhos,H);
    tab.w_cov_la = cov2tab(covw(:,:,:,:,2,:),i,j,rhos,H);
    tab.b_len_la = cov2tab(lenb(:,:,:,:,2,:),i,j,rhos,H);
    tab.e_len_la = cov2tab(lene(:,:,:,:,2,:),i,j,rhos,H);
    tab.w_len_la = cov2tab(lenw(:,:,:,:,2,:),i,j,rhos,H);
    
    tab.fin = zeros((size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),14);
    tab.fin(:,3) = reshape(tab.e_cov_c(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    tab.fin(:,4) = reshape(tab.b_cov_c(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    tab.fin(:,5) = reshape(tab.w_cov_c(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    tab.fin(:,6) = reshape(tab.e_cov_la(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    tab.fin(:,7) = reshape(tab.b_cov_la(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    tab.fin(:,8) = reshape(tab.w_cov_la(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    
    tab.fin(:,9) = reshape(tab.e_len_c(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    tab.fin(:,10) = reshape(tab.b_len_c(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    tab.fin(:,11) = reshape(tab.w_len_c(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    tab.fin(:,12) = reshape(tab.e_len_la(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    tab.fin(:,13) = reshape(tab.b_len_la(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    tab.fin(:,14) = reshape(tab.w_len_la(2:end,2:end),(size(tab.b_cov_c,1)-1)*(size(tab.b_cov_c,2)-1),1);
    
    tab.fin(:,2) = repmat(H',4,1);
    for cp = 1:4
        tab.fin((cp-1)*5+1:cp*5,1) = repmat(rhos(cp),5,1);
    end
    
    tab.fin = array2table(tab.fin);
    tab.fin_og = tab.fin;
    tab.fin = varfun(@(x) round(x,3),tab.fin);
    tab.fin = renamevars(tab.fin,1:width(tab.fin),["rho","h","LP","LP_b","LP_w","LPLA","LPLA_b","LPLA_w","LP_m","LP_b_m","LP_w_m","LPLA_m","LPLA_b_m","LPLA_w_m"]);
    
    % storing
    
    % name
    date = datetime('today');
    if strcmp(inv,'n')
        inn = 'miss_ninv';
    end
    if strcmp(wirf,'n')
        name =  [met '_' inn '_T' num2str(T) '_I' num2str(I) '_B' num2str(B) '_' num2str(day(date)) '_' num2str(month(date))];
    else
        name =  [met '_' inn '_WIRF_T' num2str(T) '_I' num2str(I) '_B' num2str(B) '_' num2str(day(date)) '_' num2str(month(date))];
    end

    % csv
    writetable(tab.fin,['outputs/csv/' name '.csv']);
    disp('Table saved as csv')
    
    % latex
    A = tab.fin;
    A = removevars(A, {'rho'});
    fid = fopen(['outputs/tables/' name '.tex'], 'w');
    fprintf(fid, ['\\begin{tabular}{' repmat('c', 1, width(A)) '}\n']);
    fprintf(fid, '\\toprule\n');
    fprintf(fid, [' & \\multicolumn{' num2str((width(A)-1)/2) '}{c}{Coverage} & \\multicolumn{' num2str((width(A)-1)/2) '}{c}{Average length} \\\\\n']);
    fprintf(fid, ['\\cmidrule(l{2pt}r{2pt}){2-' num2str((width(A)-1)/2 + 1) '} \\cmidrule(l{2pt}r{2pt}){' num2str((width(A)-1)/2 + 2) '-' num2str(width(A)) '}\n']);
    fprintf(fid, '$h$ & LP & $\\text{LP}_b$ & $\\text{LP}_w$ & LPLA & $\\text{LPLA}_b$ & $\\text{LPLA}_w$ & LP & $\\text{LP}_b$ & $\\text{LP}_w$ & LPLA & $\\text{LPLA}_b$ & $\\text{LPLA}_w$ \\\\\n');
    fprintf(fid, '\\midrule\n');
    output = cell(size(A,1), 1);
    formatStr = '%d & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f  \\\\\n' ;
    for j = 1:size(rhos,1)
        for i = 1:size(H,2)
            rowData = num2cell([table2array(A((j-1)*size(H,2)+i,:))]);
            output{i} = sprintf(formatStr, rowData{:});
        end
        if strcmp(met,'AR1')
            fprintf(fid, ['\\multicolumn{' num2str(width(A)) '}{c}{$\\rho=' num2str(rhos(j)) '$} \\\\\n']);
        else
            fprintf(fid, ['\\multicolumn{' num2str(width(A)) '}{c}{$\\max|\\lambda_i|=' num2str(rhos(j)) '$} \\\\\n']);
        end
        fprintf(fid, '%s', output{:});
        if j < size(rhos,1)
            fprintf(fid, '\\midrule\n');
        end
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fclose(fid);
    disp('Table saved as tex')
end