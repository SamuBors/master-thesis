clc
clear

%% directory
addpath functions

%% INPUTS for simulations
% number of different runs
runa = inputdlg({'Number of runs'}, 'Number of runs', [1 30], {'1'}); 
run = str2double(runa{1});

clear runa

C.T = zeros(run,1);
C.inn = cell(run,1);
C.inv = cell(run,1);
C.met = cell(run,1);
C.I = zeros(run,1);
C.B = zeros(run,1);

for i = 1:run
    % time
    Ta = inputdlg({'Time periods'}, 'Time periods', [1 30], {'100'}); 
    C.T(i) = str2double(Ta{1});
    
    % types of errors
    inna = inputdlg({'iid, garch, het or miss?'}, 'Type of innovations', [1 30], {'iid'}); 
    C.inn{i} = inna{1};
    
    % invertibility
    if strcmp(C.inn{i},'miss')
        inva = inputdlg({'Invertible MA? (y/n)'}, 'Invertibility', [1 30], {'y'}); 
        C.inv{i} = inva{1};
    else
        C.inv{i} = '';
    end
    
    % DGP
    meta = inputdlg({'AR1 or VAR3?'}, 'DGP', [1 30], {'VAR3'}); 
    C.met{i} = meta{1};
    
    % simulations
    Ia = inputdlg({'Number of simulations'}, 'Number of simulations', [1 30], {'2500'}); 
    C.I(i) = str2double(Ia{1});
    
    % bootstrap repetition
    Ba = inputdlg({'Number of bootstrap repetitions'}, 'Number of bootstrap repetitions', [1 30], {'1000'}); 
    C.B(i) = str2double(Ba{1});

    clear Ta inna meta Ia Ba inva
end

%% simulations
sim(C)

%% loading bounds
% first mat (cumulative = new)
load('outputs/mat/VAR3_miss_ninv_T100_I500_B1000_26_9.mat')
LBec = LBe;
UBec = UBe;
LBbc = LBb;
UBbc = UBb;
LBwc = LBw;
UBwc = UBw;
clear LBe UBe LBb UBb LBw UBw

% second and following mats (cumulative = cumulative + new)
load('outputs/mat/VAR3_miss_ninv_T100_I1000_B1000_20_9.mat')
LBec = cat(6,LBec,LBe);
UBec = cat(6,UBec,UBe);
LBbc = cat(6,LBbc,LBb);
UBbc = cat(6,UBbc,UBb);
LBwc = cat(6,LBwc,LBw);
UBwc = cat(6,UBwc,UBw);
clear LBe UBe LBb UBb LBw UBw

load('outputs/mat/VAR3_miss_ninv_T100_I1000_B1000_26_9.mat')
LBec = cat(6,LBec,LBe);
UBec = cat(6,UBec,UBe);
LBbc = cat(6,LBbc,LBb);
UBbc = cat(6,UBbc,UBb);
LBwc = cat(6,LBwc,LBw);
UBwc = cat(6,UBwc,UBw);
clear LBe UBe LBb UBb LBw UBw

% renaming and clearing
LBe = LBec;
UBe = UBec;
LBb = LBbc;
UBb = UBbc;
LBw = LBwc;
UBw = UBwc;
clear LBec UBec LBbc UBbc LBwc UBwc

%% input for tables
% time
Ta = inputdlg({'Time periods'}, 'Time periods', [1 30], {'100'}); 
C.T = str2double(Ta{1});

% types of errors
inna = inputdlg({'iid, garch, het or miss?'}, 'Type of innovations', [1 30], {'iid'}); 
C.inn = inna{1};

% invertibility
if strcmp(C.inn,'miss')
    inva = inputdlg({'Invertible MA? (y/n)'}, 'Invertibility', [1 30], {'y'}); 
    C.inv = inva{1};
else
    C.inv = '';
end

% DGP
meta = inputdlg({'AR1 or VAR3?'}, 'DGP', [1 30], {'VAR3'}); 
C.met = meta{1};

% simulations
Ia = inputdlg({'Number of simulations'}, 'Number of simulations', [1 30], {'2500'}); 
C.I = str2double(Ia{1});

% bootstrap repetition
Ba = inputdlg({'Number of bootstrap repetitions'}, 'Number of bootstrap repetitions', [1 30], {'1000'}); 
C.B = str2double(Ba{1});

clear Ta inna meta Ia Ba inva

%% tables
tab(LBe,UBe,LBb,UBb,LBw,UBw,C,1,1,'n')
if strcmp(C.inn,'miss')
    disp('WIRF')
    tab(LBe,UBe,LBb,UBb,LBw,UBw,C,1,1,'y')
end

doub = inputdlg({'Save doubles?'}, 'Doubles', [1 30], {'n'}); 
if strcmp(doub{1},'y')
    date = datetime('today');
    if strcmp(C.inv,'n')
        C.inn = 'miss_ninv';
    end
    name =  [C.met '_' C.inn '_T' num2str(C.T) '_I' num2str(C.I) '_B' num2str(C.B) '_' num2str(day(date)) '_' num2str(month(date))];
    save(['outputs/mat/' name '.mat'],"IRFe","LBe","UBe","LBb","UBb","LBw","UBw");
    disp('Doubles saved as mat')
end

