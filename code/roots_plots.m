%% 
clc
clear

%% directory
addpath functions

%% VAR roots plot
for cp = 2:4
    A = coeff_loop('VAR3',cp);
    r = VAR_roots(A);
    x = -1:0.001:1;
    y = (1-x.^2).^0.5;
    y1 = -y;
    plot(real(r),imag(r),"o",'Color','r')
    hold on
    plot(x,y,'Color','k')
    plot(x,y1,'Color','k')
    xline(0,'LineWidth',0.1)
    yline(0,'LineWidth',0.1)
    grid on
    axis equal
    lim = max([max(real(r)),min(real(r)),max(imag(r)),min(imag(r))])+0.25;
    xlim([-lim,lim])
    ylim([-lim,lim])
    ax = gca;
    exportgraphics(ax,['outputs/img/roots_cp' num2str(cp) '.png'],'Resolution',300) 
    hold off
end

%% invertible MA roots plot (AR)
T = 25;
n = 0.3;
G{1,1} = 1;
G{1,2} = -1;
G{1,3} = 0.5;
G{1,4} = -0.5;
G{1,5} = -0.7;
G{1,6} = 0.6;
G{1,7} = 0.3;
G{1,8} = 0.4;
G{1,9} = 0.8;
G{1,10} = -0.5;

for q = 2:size(G,2)
    G{1,q} = G{1,q}*T^(-n);
end

r = VMA_roots(G);
x = -1:0.001:1;
y = (1-x.^2).^0.5;
y1 = -y;
plot(real(r),imag(r),"o",'Color','r')
hold on
plot(x,y,'Color','k')
plot(x,y1,'Color','k')
xline(0,'LineWidth',0.1)
yline(0,'LineWidth',0.1)
grid on
axis equal
lim = max([max(real(r)),min(real(r)),max(imag(r)),min(imag(r))])+0.25;
xlim([-lim,lim])
ylim([-lim,lim])
ax = gca;
exportgraphics(ax,['outputs/img/roots_MA_AR1_' num2str(T) '.png'],'Resolution',300) 
hold off

%% invertible VMA roots plot (VAR)
T = 25;
n = 0.3;
G{1,1} = eye(2);
G{1,2} = [0.4,.5;-1,.7];
G{1,3} = [.3,-.3;-0.9,1];
G{1,4} = [-.5,.1;-0.7,0.8];
G{1,5} = [1.2,0;-0.8,0.6];
G{1,6} = -[.4,-.3;-0.3,0.5];
G{1,7} = [0,.5;1,0];
G{1,8} = -[0.1,-0.2;0.2,0.1];
G{1,9} = [.1,-.1;.2,0.1];
G{1,10} = [.9,.4;-.4,.3];

for q = 2:size(G,2)
    G{1,q} = G{1,q}*T^(-n);
end

r = VMA_roots(G);
disp(abs(r))
x = -1:0.001:1;
y = (1-x.^2).^0.5;
y1 = -y;
plot(real(r),imag(r),"o",'Color','r')
hold on
plot(x,y,'Color','k')
plot(x,y1,'Color','k')
xline(0,'LineWidth',0.1)
yline(0,'LineWidth',0.1)
grid on
axis equal
lim = max([max(real(r)),min(real(r)),max(imag(r)),min(imag(r))])+0.25;
xlim([-lim,lim])
ylim([-lim,lim])
ax = gca;
exportgraphics(ax,['outputs/img/roots_MA_VAR3_' num2str(T) '.png'],'Resolution',300) 
hold off

%% non-invertible MA roots plot (AR)
T = 25;
n = 0.3;
G{1,1} = 1;
G{1,2} = -1;
G{1,3} = 1.5;
G{1,4} = -0.8;
G{1,5} = -1;
G{1,6} = 1.2;
G{1,7} = 2;
G{1,8} = 1.2;
G{1,9} = 0.8;
G{1,10} = -1;

for q = 2:size(G,2)
    G{1,q} = G{1,q}*T^(-n);
end

r = VMA_roots(G);
disp(abs(r))
x = -1:0.001:1;
y = (1-x.^2).^0.5;
y1 = -y;
plot(real(r),imag(r),"o",'Color','r')
hold on
plot(x,y,'Color','k')
plot(x,y1,'Color','k')
xline(0,'LineWidth',0.1)
yline(0,'LineWidth',0.1)
grid on
axis equal
lim = max([max(real(r)),min(real(r)),max(imag(r)),min(imag(r))])+0.25;
xlim([-lim,lim])
ylim([-lim,lim])
ax = gca;
exportgraphics(ax,['outputs/img/roots_MA_AR1_' num2str(T) '_NINV.png'],'Resolution',300) 
hold off

%% non-invertible VMA roots plot (VAR)
T = 25;
n = 0.3;
G{1,1} = eye(2);
G{1,2} = [1.5,.5;-1,.7];
G{1,3} = [.3,-.3;-1,1];
G{1,4} = [-.5,.1;-0.7,0.8];
G{1,5} = [1.2,0;-1,0.6];
G{1,6} = -[.4,-.3;-0.3,0.5];
G{1,7} = [1,.5;1,1.5];
G{1,8} = -[0.1,-0.2;0.2,0.2];
G{1,9} = [.1,-.1;.2,0.1];
G{1,10} = [1,1.5;-.4,.3];

for q = 2:size(G,2)
    G{1,q} = G{1,q}*T^(-n);
end

r = VMA_roots(G);
disp(abs(r))
x = -1:0.001:1;
y = (1-x.^2).^0.5;
y1 = -y;
plot(real(r),imag(r),"o",'Color','r')
hold on
plot(x,y,'Color','k')
plot(x,y1,'Color','k')
xline(0,'LineWidth',0.1)
yline(0,'LineWidth',0.1)
grid on
axis equal
lim = max([max(real(r)),min(real(r)),max(imag(r)),min(imag(r))])+0.25;
xlim([-lim,lim])
ylim([-lim,lim])
ax = gca;
exportgraphics(ax,['outputs/img/roots_MA_VAR3_' num2str(T) '_NINV.png'],'Resolution',300) 
hold off

%% MA creation check
I = 1000;
T = 1000;
n = inf;
G{1,1} = 1;
G{1,2} = -1;
G{1,3} = 0.5;
G{1,4} = -0.5;
G{1,5} = -0.7;
G{1,6} = 0.6;
G{1,7} = 0.3;
G{1,8} = 0.4;
G{1,9} = 0.8;
G{1,10} = -0.5;

for q = 2:size(G,2)
    G{1,q} = G{1,q}*T^(-n);
end

r = VMA_roots(G);
disp(abs(r))

M = 1;
Q = size(G,2)-1; % highest lag of the MA process
Gm = cell2mat(G);
cor = zeros(21,I);
for i = 1:I
    u = randn(T+Q,1); % iid components
    U = lagmatrix(u,0:Q);
    U = U(Q+1:end,:);
    e = U*Gm'; 
    cor(:,i) = autocorr(e);
end
corm = mean(cor,2);
disp(corm)

%% VARMA function check
I = 1000;
A = coeff_loop('VAR3',2);
e = zeros(21,I);
for i = 1:I
    X = VARMA_sim(T,A,G,W);
    mod = VAR_est(X,2);
    e(:,i) = autocorr(mod.res(:,1));
end
em = mean(e,2);
disp(em)





