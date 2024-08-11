clc;clear;

load('good plots\plots_for_paper\3_l1.mat');
load('good plots\plots_for_paper\3_l1data.mat');
Xl1_dc = X_interval;

load('good plots/plots_for_paper/3_l2.mat');
load('good plots/plots_for_paper/3_l2data.mat');
Xl2_dc = X_interval;

N2 = 35;N=8;m=4;
n = length(Xl1_dc)/(4*N2);

% Europa frame
l2 = 1.0205;
l1 = 0.9798;
plot(l1,0,'*k','MarkerSize',9); hold on;

for i=0:n-1
    colorr = [1-(i/n), 0, (i/n)];
    plot_collected(Xl1_dc(1+(4*N2*i):4*N2+(4*N2*i)),'position-connect',colorr);    % X0
    plot_collected(Xl1_dc(1+(4*N2*i):4*N2+(4*N2*i)),'position',colorr);
    pause(0.1)
end
