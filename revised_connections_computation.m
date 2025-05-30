clc;clear;

%%%% Units used in this script are important to understand --- firstly
%%%% there is a distance scaling were ganymede's semi major is taken as 5
%%%% times this canonical unit, secondly there is a mass scaling involved
%%%% where 10^(20)kg is the mass canonical unit

%%%% Functions in use : crtbp_J, diffcorr_l1_FBP, diffcorr_l2_FBP, PCC4BP_eqn

tStep = 1*10^(-3);

MU = 10^20;         %kg
DU = 1070400000/5;  %m
%TU = 7.155*24*3600; %s
TU = 4.8824e+04; % To make Europa mean motion = 1
%TU = 2.2922e+05; % To make Callisto mean motion = 1

M = 19000000;
m_1 = 893;
m_2 = 480;
m_3 = 1480;
m_4 = 1075.9;
e_1 = 0.004;        
e_2 = 0.009;
e_3 = 0.0013;
e_4 = 0.0074;
a_1 = 1.97; 
a_2 = 3.135; 
a_3 = 5;
a_4 = 8.79;
R_europa = 1560000/DU;
a_syn = a_2;

G = 6.6743*10^(-11)*(MU*TU^2/(DU^3));   % Taking care of rescaling in length, mass and time

laplace_resonance = [4,2,1]; % Io-Europa-Ganymede
T_1 = ((4*pi^2/(G*M))*a_1^3)^(1/2);
T_2 = ((4*pi^2/(G*M))*a_2^3)^(1/2);
T_3 = ((4*pi^2/(G*M))*a_3^3)^(1/2);
T_4 = ((4*pi^2/(G*M))*a_4^3)^(1/2);

% Rotating to Inertial %%%%%
n_1 = sqrt(G*(M+m_1)/(a_1^3));
n_2 = sqrt(G*(M+m_2)/(a_2^3));
n_3 = sqrt(G*(M+m_3)/(a_3^3));
n_4 = sqrt(G*(M+m_4)/(a_4^3));
Rot_mat = @(n,t) [cos(n*t) -sin(n*t); sin(n*t) cos(n*t)];

%%%%%% Setup for the synodic frame %%%%% (m_2:EUROPA | m_4:CALLISTO)
mu1 = m_2/(M+m_2);
mu2 = m_3/(M+m_2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  L1 Computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% L1 (between Jupiter and moons)

f = @(x, mu) (x+mu)^2*x*(x-(1-mu))^2 - (1-mu)*(x-(1-mu))^2 + mu*(x+mu)^2;
g = @(x) f(x, mu1);
l1_pos_io = fzero(g, 0);

A_3d = crtbp_J(l1_pos_io, 0, 0, mu1);
A = [A_3d(1:2,1:2) A_3d(1:2,4:5); A_3d(4:5,1:2) A_3d(4:5,4:5)];

[l1_eigvec_io, l1_eigv_io] = eig(A);
l1_eigvec_io = real(l1_eigvec_io); % Cause these eigenvectors will be used in real space

T_p = 3*pi/abs(l1_eigv_io(3,3)); % Third eigenvalue of the system in the form +-i*nu
del_x = 8e-3;  % 3e-3  5e-3
tspan = 0:0.001:T_p;
x0 = [l1_pos_io;0;0;0] + del_x*l1_eigvec_io(:,3);
x0 = [x0; 0];                                                                  % The fifth state is the angle of third primary with first
[a_l1, b_l1, STM_b] = diffcorr_l1_FBP(x0, mu1, 0, [(a_3/a_syn)+mu1,0], n_3, tStep);      % a --> initial state corrected   2*b --> time period of the periodic orbit  STM_b --> STM matrix from 0 to b
[t_dc, Xl1_dc] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,0,[(a_3/a_syn)+mu1,0], n_3, 3), 0:0.001:2*b_l1, a_l1);
Xl1_dc(:,5) = x0(5)*ones([length(Xl1_dc(:,5)),1]);  % third primary angle for all the entries of periodic orbit should have same value as the initial condition x0
w_l1 = 2*pi/(2*b_l1);
% plot(Xl1_dc(:,1),Xl1_dc(:,2),'-k','LineWidth',2);hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  L2 Computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L2 (beyond Jupiter and moons)
 
% f = @(x, mu) (x+mu)^2*x*(x-(1-mu))^2 - (1-mu)*(x-(1-mu))^2 - mu*(x+mu)^2;
% g = @(x) f(x, mu1);
% l2_pos_io = fzero(g, 0);
% 
% A_3d = crtbp_J(l2_pos_io, 0, 0, mu1);
% A = [A_3d(1:2,1:2) A_3d(1:2,4:5); A_3d(4:5,1:2) A_3d(4:5,4:5)];
% 
% [l2_eigvec_io, l2_eigv_io] = eig(A);
% l2_eigvec_io = real(l2_eigvec_io); % Cause these eigenvectors will be used in real space
% 
% T_p = 3*pi/abs(l2_eigv_io(3,3)); % Third eigenvalue of the system in the form +-i*nu
% %del_x = -1.0345*10^(-1);        % 4x T_3
% %del_x = -7.72084057245*10^(-2);  % ~3.5x T_3
% del_x = -2e-2;
% tspan = 0:0.001:T_p;
% x0 = [l2_pos_io;0;0;0] + del_x*l2_eigvec_io(:,3);
% x0 = [x0; 0];                                                                  % The fifth state is the angle of third primary with first 
% [a_l2, b_l2, STM_b] = diffcorr_l2_FBP(x0, mu1, 0, [(a_3/a_syn)+mu1,0], n_3, tStep);      % a --> initial state corrected   2*b --> time period of the periodic orbit  STM_b --> STM matrix from 0 to b
% [t_dc, Xl2_dc] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,0,[(a_3/a_syn)+mu1,0], n_3, 3), 0:0.001:2*b_l2, a_l2);
% Xl2_dc(:,5) = x0(5)*ones([length(Xl2_dc(:,5)),1]);  % third primary angle for all the entrie  s of periodic orbit should have same value as the initial condition x0 
% w_l2 = 2*pi/(2*b_l2);
% %plot(Xl2_dc(:,1),Xl2_dc(:,2),'-k','LineWidth',2);hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% syms t
% figure(1)
% hold on;
% plot(Xl1_dc(:,1), Xl1_dc(:,2),'-k'); hold on;
% plot(Xl2_dc(:,1), Xl2_dc(:,2),'-k'); hold on;
% fplot((1-mu1)+R_europa*sin(t), R_europa*cos(t));
% %hold on;fplot(((a_3/a_syn)+mu1)*sin(t),((a_3/a_syn)+mu1)*cos(t),'r','Linewidth',4)
% pause(0.000000001)
% hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTATION OF PERIODIC ORBIT OVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_dc = Xl1_dc;
l_pos = l1_pos_io;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTATION OF QUASI-PERIODIC ORBITS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dphii_l1, phii_l1, x_l1] = get_STM(t_dc(end),X_dc(1,:)',4,0);  % Taking out STM for periodic orbit from TBP but now calculating it for FBP
[rvec, rho] = eig(phii_l1);
if angle(rho(3,3)) ~= 0 && angle(rho(3,3)) ~= pi
    rho_l1 = angle(rho(3,3));
else
    rho_l1 = 1;
end

theta_l1 = parameterize(X_dc,l_pos,'2D');
%%%%%%% Taking periodic orbit from 3BP as an initial guess for Newton's method, find invariant cicle for the stroboscopic map

%%%%%%%% DFT of initial states section %%%%%%%%
N2 = 75; % Total number of discrete points on invariant circle (KEEP IT ODD)  %41 (l1)  %39 (l2)    %commons 35  |||| Callisto (25 works for 9e-3)
m2 = (N2-1)/2;  % FS truncation order
k = -(N2-1)/2:(N2-1)/2;
%%% (For now) Use Newton method only for one trajectory (i.e only one vector from collected state vector)
m = 4;  % Order of legendre polynomial (SHOULD BE EVEN)  %4 (l1)
N = 8; %50 % Number of collocation points %12 (l1)  (NEW: N=8(l1 Europa))        %%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT : CHANGING THIS IS NOT CHANGING THE ACCURACY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newton_iter = 25;  % (l1 35| l2 22)  (NEW: 2_l1.mat - 9) l1 - 11   l2 - 16   (new l1 family 37)   25-Callisto

len_circle = 4*N2;   % Total number of discrete invariant circle
len_coll = N*(m+1);   % Total number of collocation + continuity constraint
w = exp(-1i*2*pi/N2);

% Initial approximation of family tangent vector
% X_tang = [];
% for i=1:length(t_dc)
%     [dPh, PHI, xxff] = get_STM(2*b_l1,Xl1_dc(i,:)',4,1);
%     [evec,eval] = eig(PHI);
%     uu = real(exp(1i*theta_l1(i))*evec(:,3));
%     X_tang = [X_tang 10^(-5)*uu/norm(uu)];
% end
U = 0;
U_exp = X_dc(:,1:4) + 1*10^(-5)*U';  % Linear approximation of invariant curve (length(t_dc),5)

% Creating initial Family Tangent vector
% del = 1.5*del_x;
% X_tang = get_nearby_orbit(l1_pos_io,del,l1_eigvec_io,[mu1,a_2,a_3]);

X_coll0 = [];
T_coll0 = [];
X_tang0 = [];
THETA = [];
IDX = [];
for i=1:N2
    angle = ((i-1)/N2)*2*pi;
    [val, idx] = min(abs(theta_l1-angle));
    IDX = [IDX idx];
    THETA = [THETA; theta_l1(idx)];
    X_coll0 = [X_coll0; U_exp(idx,1:4)'];    % Size n*N2
    T_coll0 = [T_coll0; t_dc(idx)];
    %X_tang0 = [X_tang0; X_tang(1:4,idx)];   % Collected equivalent of family tangent vector
end

%%%%%%%%%%%% Newton's methods iterations %%%%%%%%%%%%%%%%%%
DATA = [];
%mu_range = 0:mu2/20:mu2;
theta_range = 0;
for theta_init=theta_range
%X_tang_vec = X_tang0;
X_tang_vec = zeros(4*N2,1);
%X_tang_vec = X_tang_vec/sqrt(norm(X_tang_vec)^2+rho_l1^2);
w_coll0 = rho_l1*(n_3-1);  % For phase condition

% creating DFT and inverse DFT matrix
DFT = zeros(length(X_coll0),length(X_coll0));
Q = zeros(length(X_coll0),length(X_coll0));
for i=1:N2
    Q(1+4*(i-1):4+4*(i-1),1+4*(i-1):4+4*(i-1)) = exp(-1i*(-0.5*(N2-1)+(i-1))*rho_l1)*eye(4);
    for j=1:N2
        DFT(1+4*(i-1):4+4*(i-1),1+4*(j-1):4+4*(j-1)) = w^((j-1)*(-0.5*(N2-1) + (i-1)))*eye(4);
    end
end
KKK = sparse(DFT'*Q*DFT/N2);  % 4*N2 square  (Lujan)
D_rho = DFT'*(kron(1i*diag(-(N2-1)/2:(N2-1)/2),eye(4)))*DFT/N2;  % 4*N2 square  (Lujan)
X_dtheta = real(D_rho*X_coll0);   % 4*N2 vector
LL = sparse(real(D_rho*KKK));     % 4*N2 square
%%%%%%%%%%%%%%%%%%%%%% Gauss Legendre section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmda = 1e-12; % Just an initial guess   1e-12   (Changes required in F_coll for augmented part calculation of dvdtheta)
psued_arc = 1e-18;   % -1e-4 (l1 3e-3)   1e-12 (l2 2e-2)     1e-18

arc_flag = 0;
X_prev = X_coll0;
rho_prev = rho_l1;
rho_tang = 0;
X_interval = kron(ones(N*(m+1)+1,1),X_coll0); % x_interval for each invariant point on orbit
X_interval = [X_interval; rho_l1; lmda]; % THIS VARIABLE NEED TO BE CONVERGED DURING NEWTON'S METHOD (HAS ALL THE INTERVAL) TOTAL SIZE --> N2*(4*N*(m+1)+4)+1 SQUARE
G_mat = [];
J_lmda = [];   % Sets up the jacobian for lmda unfolding parameter
XRHO = [];  % KEEPS TRACK OF RHOs FOR EACH CALCULATED ORBIT
% for h=1:newton_iter
%            J_total = zeros((len_circle*len_coll+len_circle)+2,(len_circle*len_coll+len_circle)+2);    % REFER THE SIZE INFO FROM OLIKARA FIG.2.1 (Added coll+cont+quasiperiod+phase | Xij+Xn0+rho)
%            G_total = [];
%            Rot = real(KKK); 
%            X_init = X_interval(1 : 4*N2);
%            X_ends = X_interval(end-(4*N2)-1 : end-2);
%            x_rot = Rot*X_ends;
% 
%            % Adding collocation + continuity condition
%            for i=1:N
%                tau_interval = [(i-1)/N i/N];
%                X_ = X_interval(1+len_circle*(m+1)*(i-1):(len_circle*(m+2))+len_circle*(m+1)*(i-1));
%                [J_lmda, J, G] = BVP(X_, tau_interval(1), 1/N, m, N2, [X_interval(end-1), X_interval(end), ((a_3/a_syn)+mu1), n_3, mu1, mu2, l_pos, w_coll0], D_rho, theta_init);
%                J_total(1+len_circle*(m+1)*(i-1):(len_circle*(m+1))+len_circle*(m+1)*(i-1), 1+len_circle*(m+1)*(i-1):(len_circle*(m+2))+len_circle*(m+1)*(i-1)) = J;   % N*(4*N2*(m+1)) X N*(4*N2*(m+1)+4)
%                %J_total(1+len_circle*(m+1)*(i-1):(len_circle*m)+len_circle*(m+1)*(i-1), end) = J_lmda;
%                G_total = [G_total; G]; 
%            end
% 
%            % Adding QuasiPeriodic Condition
%            G_total = [G_total; X_init-x_rot];  % 4*N2
%            J_total(end-(4*N2)-1:end-2, 1:4*N2) = eye(4*N2); 
%            J_total(end-(4*N2)-1:end-2, end-(4*N2)-1:end-2) = -Rot;    % d/d(X_N,0)
% 
%            J_total(end-(4*N2)-1:end-2, end-1) = LL*X_ends;
% 
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%            %%%%%%%%% Adding Phase Condition
%            fff = F_coll(X_prev,[0,0],X_interval(end),[(a_3/a_syn)+mu1,0],[n_3, mu1, mu2, 4, l_pos],theta_init);  % For phase condition
%            Xdiff = (X_init-X_prev);
%            X_dtheta1 = (1/w_coll0)*(fff - (n_3-1)*X_dtheta);
%            G_total = [G_total; dot(Xdiff, X_dtheta1)];  
%            J_total(end-1, 1:4*N2) = X_dtheta1';
%            J_total(end-1, end-1) = dot(Xdiff,-(1/w_coll0^2)*(fff - (n_3-1)*X_dtheta)*(n_3-1));   % creates no issues
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            %%%%%%%%% Adding Pseudo-arclength Condition
%            G_total = [G_total; (1/N2)*dot((X_init-X_coll0), X_tang_vec) + (X_interval(end-1)-rho_l1)*rho_tang - psued_arc];
%            J_total(end, 1:4*N2) = X_tang_vec';
%            J_total(end, end-1) = rho_tang;
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%            if arc_flag == 0
%                 X_prev = X_init;
%                 rho_prev = X_interval(end-1);
%            end
% 
%            G_mat = [G_mat norm(G_total)];
%            if norm(G_total)<1e-18  % Check if the solution is very close to the approx periodic orbit then break (l1 - 1e-6  |  l2 - 1e-8) abs(G_total(end))<1e-18 && G_total(end) ~= -psued_arc
%                 fprintf('Break %i',h)
%                 break
%            end
% 
%            X_interval = X_interval - pinv(J_total)*G_total;
%            w_coll0 = X_interval(end-1)*(n_3-1);
%            %%%% Update tangent vector for use in pseudo-arclength and phase condition
%            X_dtheta = real(D_rho*X_init);  % X_dtheta extracted from previous newton stage X_init
%            for i=1:N2
%                 Q(1+4*(i-1):4+4*(i-1),1+4*(i-1):4+4*(i-1)) = exp(-1i*(-0.5*(N2-1)+(i-1))*X_interval(end-1))*eye(4);
%            end
%            KKK = sparse(DFT'*Q*DFT/N2);
%            LL = sparse(real(D_rho*KKK));
% 
%            if arc_flag == 0
%                X_tang_vec = X_interval(1:4*N2) - X_prev;
%                rho_tang = X_interval(end-1) - rho_l1; % rho_l1 or rho_prev (Use rho_l1 for 3e-3)
%                X_tang_vec = X_tang_vec/sqrt(norm(X_tang_vec)^2+(rho_tang)^2);  %Normalization
%                rho_tang = rho_tang/sqrt(norm(X_tang_vec)^2+(rho_tang)^2);   %Normalization
%                arc_flag = 1;
%            end
% 
%            if norm(G_total(1:end-2))<5e-5  % Only checking for collocation + continuity solutions  (5e-5)
%                arc_flag = 1;
%            end
%            if h>0      %9
%                colorr = [rand, 0, rand];
%                XRHO = [XRHO X_interval(end-1)]
%                figure(1);
%                %plot_collected(X_interval(1:4*N2),'position-connect1',colorr);    % X0
%                %plot_collected(X_interval(1:4*N2),'position',colorr);
%                % figure(2);
%                % plot_collected(X_interval((m+1)*4*N*N2 + 1: (m+1)*4*N*N2 + 4*N2),'position-connect2 -',colorr);    % Xends
%                % plot_collected(X_interval((m+1)*4*N*N2 + 1: (m+1)*4*N*N2 + 4*N2),'position',colorr);   % Xends
%                %plot_collected(real(KKK)*X_interval((m+1)*4*N*N2 + 1: (m+1)*4*N*N2 + 4*N2),'position',[0 1 0]);  % Xrot
%                %pause(0.0000001)
%                % % % % % % % calc_jacobi(X_interval(1:4*N2),mu1)
%                % % % % % % % if abs(vpa(calc_jacobi(X_interval(1:4*N2),mu1),6) - (-1.50100)) < 0.0001 %0.00007      % -1.50309 => 1c_l1    -1.50280 => 2c_l1
%                % % % % % % %     save('good plots\plots_for_paper\3c_l2','X_interval')
%                % % % % % % %     save('good plots\plots_for_paper\3c_l2data','J_total','G_total')
%                % % % % % % %     break
%                % % % % % % % end
%            end
% end
% figure(1);
% plot_collected(X_interval(1:4*N2),'position',colorr);
% 
% A_dim = [len_circle*len_coll, len_circle];
% B_dim = [len_circle*len_coll, len_circle*len_coll];
% C_dim = [len_circle*len_coll, 2];  % 2 no. of free parameters rho and lmda
% 
% Al1 = J_total(1:A_dim(1), 1:A_dim(2));
% Bl1 = J_total(1:B_dim(1), 1+A_dim(2):B_dim(2)+A_dim(2));
% Cl1 = J_total(1:C_dim(1), 1+B_dim(2)+A_dim(2):B_dim(2)+A_dim(2)+C_dim(2));
% 
% Phi_l1 = [eye(length(Al1(1,:))); -inv(Bl1)*Al1];
% PPPHHH1 = Phi_l1(end-(4*N2)+1:end,:);
% Gw1 = Rot*PPPHHH1;
% [torus_vec1, torus_val1] = eig(Gw1);
% eval1 = [];
% for i=1:length(torus_val1)
%     eval1 = [eval1 torus_val1(i,i)]; 
% end
% max_v = max(abs(eval1));
% min_v = min(abs(eval1));
% DATA = [DATA; [theta_init max_v min_v]]
end
    %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STABILITY & MANIFOLD CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: RESULTS HIGHLY DEPEND ON VARIABLE "w", FIX IT!!!
clear Q; clear DFT;

flag = 0;  %0 --> Callisto 1 --> Europa
load('good plots\plots_for_paper\3c_l1.mat');
load('good plots\plots_for_paper\3c_l1data.mat');
Xl1_dc = X_interval;
Jtot_l1 = J_total;
Gtot_l1 = G_total;

load('good plots/plots_for_paper/3c_l2.mat');
load('good plots/plots_for_paper/3c_l2data.mat');
Xl2_dc = X_interval;
Jtot_l2 = J_total;
Gtot_l2 = G_total;

% These stability matrices are only taking collocation + continuity constraints
N2_l1 = 25; N2_l2 = 25; %  Europa pairs: [35 29    41 35]  Callisto pairs: [25 25]
len_circle = 4*N2_l1;
A_dim = [len_circle*len_coll, len_circle];
B_dim = [len_circle*len_coll, len_circle*len_coll];
C_dim = [len_circle*len_coll, 2];  % 2 no. of free parameters rho and lmda

Al1 = Jtot_l1(1:A_dim(1), 1:A_dim(2));
Bl1 = Jtot_l1(1:B_dim(1), 1+A_dim(2):B_dim(2)+A_dim(2));
Cl1 = Jtot_l1(1:C_dim(1), 1+B_dim(2)+A_dim(2):B_dim(2)+A_dim(2)+C_dim(2));

len_circle = 4*N2_l2;
A_dim = [len_circle*len_coll, len_circle];
B_dim = [len_circle*len_coll, len_circle*len_coll];
C_dim = [len_circle*len_coll, 2];  % 2 no. of free parameters rho and lmda

Al2 = Jtot_l2(1:A_dim(1), 1:A_dim(2));
Bl2 = Jtot_l2(1:B_dim(1), 1+A_dim(2):B_dim(2)+A_dim(2));
Cl2 = Jtot_l2(1:C_dim(1), 1+B_dim(2)+A_dim(2):B_dim(2)+A_dim(2)+C_dim(2));

Phi_l1 = [eye(length(Al1(1,:))); -inv(Bl1)*Al1]; Psi_l1 = [eye(length(Cl1(1,:))); -inv(Bl1)*Cl1];
for i=1:N2_l1
    Q(1+4*(i-1):4+4*(i-1),1+4*(i-1):4+4*(i-1)) = exp(-1i*(-0.5*(N2_l1-1)+(i-1))*Xl1_dc(end-1))*eye(4);
    for j=1:N2_l1
        DFT(1+4*(i-1):4+4*(i-1),1+4*(j-1):4+4*(j-1)) = w^((j-1)*(-0.5*(N2_l1-1) + (i-1)))*eye(4);
    end
end
Rot_l1 = real(DFT'*Q*DFT/N2_l1);
PPPHHH1 = Phi_l1(end-(4*N2_l1)+1:end,:);
Gw1 = Rot_l1*PPPHHH1;
[torus_vec1, torus_val1] = eig(Gw1);

Phi_l2 = [eye(length(Al2(1,:))); -inv(Bl2)*Al2]; Psi_l2 = [eye(length(Cl2(1,:))); -inv(Bl2)*Cl2];
Q = []; DFT = [];
for i=1:N2_l2
    Q(1+4*(i-1):4+4*(i-1),1+4*(i-1):4+4*(i-1)) = exp(-1i*(-0.5*(N2_l2-1)+(i-1))*Xl2_dc(end-1))*eye(4);
    for j=1:N2_l2
        DFT(1+4*(i-1):4+4*(i-1),1+4*(j-1):4+4*(j-1)) = w^((j-1)*(-0.5*(N2_l2-1) + (i-1)))*eye(4);
    end
end
Rot_l2 = real(DFT'*Q*DFT/N2_l2);
PPPHHH2 = Phi_l2(end-(4*N2_l2)+1:end,:);
Gw2 = Rot_l2*PPPHHH2;
[torus_vec2, torus_val2] = eig(Gw2);

% Plot calculated quasi-periodic orbits
% hold on;
% plot_collected(Xl1_dc(1:4*N2_l1),'position',[1 0 0]) 
% plot_collected(Xl2_dc(1:4*N2_l2),'position',[0 0 1])

if flag == 1
    R = 1560000/DU;syms T;%fplot((1-mu1)+R*sin(T),R*cos(T),'k','LineWidth',2);
    global mu1;
    mu1 = m_2/(M+m_2);
    mu2 = m_3/(M+m_2);
    a_syn = a_2;
    TU = 4.8824e+04;
else
    R = 2403000/DU;syms T;%fplot((1-mu1)+R*sin(T),R*cos(T),'k','LineWidth',2);
    global mu1;
    mu1 = m_4/(M+m_4);
    mu2 = m_3/(M+m_4);
    a_syn = a_4;
    TU = 2.2922e+05;
end
G = 6.6743*10^(-11)*(MU*TU^2/(DU^3));  % This gravitational parameter used for manifold calculation only
n_3 = sqrt(G*(M+m_3)/(a_3^3));

opt = odeset('RelTol',1e-15,'AbsTol',1e-15,'Events',@ps_crossing);

% Stable/unstable manifold calculation L2 stable and L1 unstable manifold
eval1 = []; eval2 = [];
for i=1:length(torus_val1)
    eval1 = [eval1 torus_val1(i,i)]; 
end
for i=1:length(torus_val2)
    eval2 = [eval2 torus_val2(i,i)]; 
end
[max_v, max_ind] = max(abs(eval1)); % Finding maximum eigenvalue for torus1
[min_v, min_ind] = min(abs(eval2)>2e-1); % Finding minimum eigenvalue for torus2
check2 = abs(torus_val2(min_ind,min_ind)); check1 = abs(torus_val1(max_ind,max_ind));
eigvec1 = real(torus_vec1(:,max_ind));
eigvec2 = real(torus_vec2(:,min_ind));

epsilon_st = -5.572e-1;  %(blue)  E - 3e-15           || C - 9e-3       -5.572e-1 gives interesting results    -5.57215e-1 better results
epsilon_unst = 5.66e-1; %(red)  E - -3e-15  -3e-3 || C - -3e-5          5.66e-1 gives Ganymede transition energy for Callisto frame (kk=22 unstable)   5.67e-4 interesting
propag_length = 700*24*3600/TU; %300*pi/abs(n_3-1);  % 300 days, 800 days
tspan_st = 0:-0.01:-propag_length;  % for L2
tspan_unst = 0:0.01:propag_length;  % for L1
Xst0 = Xl2_dc(1:4*N2_l2) + epsilon_st*eigvec2; % Stable direction (mag less than 1)  (end-6)  (106)
Xunst0 = Xl1_dc(1:4*N2_l1) + epsilon_unst*eigvec1;  % Unstable direction (mag more than 1)  (29)
STORE_POINCARE_STABLE = [];
STORE_POINCARE_UNSTABLE = [];
relevants = [];

% for kk=0:N2_l2-1
%     [T_st, M_st, t_st, X_st, i_st] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,mu2,[(a_3/a_syn)+mu1,0], n_3, 3), tspan_st, [Xst0(1+4*kk:4+4*kk);0], opt);
%     if any(sqrt(sum((M_st(:,1:2)-[1-mu1 0]).^2, 2))<=R)  %29 18 12
%         fprintf('removing %d\n',kk)
%     else
%         %if any(sqrt(sum((M_st(:,1:2)).^2, 2))>=(a_3/a_syn)-0.1)    % This is done to find trajectories that are close to Ganymede in Europa frame 
%         if kk==33   
%             relevants = [relevants kk];
%             % X_iner = plot_inertial(T_st, M_st);
%             % figure(3);plot(X_iner(:,1),X_iner(:,2),'-b','LineWidth',1);hold on;fplot(-mu1+((a_3/a_syn)+mu1)*sin(T),((a_3/a_syn)+mu1)*cos(T),'-m','LineWidth',2);hold on;fplot(-mu1+sin(T),cos(T),'-k','LineWidth',2);figure(1)
%             % hold on; plot(M_st(:,1),M_st(:,2),'-b','LineWidth',1)
%             % pause(0.0001)
%             % figure(5);
%             y = transform_frame(M_st, (a_3/a_syn)+mu1, n_3, mu1, m_3/(M+m_3), 0);
%             % for i=1:length(M_st)
%             %     ss = calc_jacobi(y(i,:)',m_3/(M+m_3),m_2/(M+m_3),a_2/a_syn);
%             %     plot(T_st(i)*TU/(3600*24),ss,'.b');hold on;
%             % end 
%             idxx = find(abs(y(:,2))<=9e-3); % Poincare for y=0
%             hold on;plot(y(idxx,1),y(idxx,3),'.b')
%             break;
%         end
%     end
% 
%     if ~isempty(X_st)
%         intersects = X_st(abs(X_st(:,2))<0.015,:);
%     else
%         intersects = [];
%     end
%     STORE_POINCARE_STABLE = [STORE_POINCARE_STABLE intersects'];
% end

for kk=0:N2_l1-1
    [T_unst, M_unst, t_unst, X_unst, i_unst] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,mu2,[(a_3/a_syn)+mu1,0], n_3, 3), tspan_unst, [Xunst0(1+4*kk:4+4*kk);0], opt);
    % if kk==22
    %     kk
    % end
    if any(sqrt(sum((M_unst(:,1:2)-[1-mu1 0]).^2, 2))<=R)
        %fprintf('removing\n')
    else
        %if any(sqrt(sum((M_unst(:,1:2)).^2, 2))<(a_3/a_syn))  % This is done to find trajectories that are close to Ganymede in Callisto frame
        if kk==22
            relevants = [relevants kk];
            % % % % % X_iner = plot_inertial(T_unst, M_unst);
            % % % % % figure(3);plot(X_iner(:,1),X_iner(:,2),'-r');hold on;fplot(-mu1+((a_3/a_syn)+mu1)*sin(T),((a_3/a_syn)+mu1)*cos(T),'-m','LineWidth',2);fplot(-mu1+sin(T),cos(T),'-g','LineWidth',2);figure(1)
            % % % % % hold on; plot(M_unst(:,1),M_unst(:,2),'-r','LineWidth',0.8)
            % % % % % pause(0.0001)
            y = transform_frame(M_unst, (a_3/a_syn)+mu1, n_3, mu1, m_3/(M+m_3), 0);
            % for i=1:length(M_unst)
            %     ss = calc_jacobi(M_unst(i,:)',mu1,mu2,a_3/a_syn);
            %     plot(T_unst(i)*TU/(3600*24),ss,'.r');hold on;+
            % end
            idxx = find(abs(y(:,2))<=9e-3); % Poincare for y=0
            hold on;plot(y(idxx,1),y(idxx,3),'.r')
            break;
        end
    end

    if ~isempty(X_unst)
        intersects = X_unst(abs(X_unst(:,2))<0.015,:);
    else
        intersects = [];
    end
    STORE_POINCARE_UNSTABLE = [STORE_POINCARE_UNSTABLE intersects'];
end

% figure(1);R_gany_orb = (a_3/a_syn)+mu1;syms T;fplot(-mu1+R_gany_orb*sin(T),R_gany_orb*cos(T),'-m','LineWidth',3);

% figure(2);plot(STORE_POINCARE_STABLE(2,:),STORE_POINCARE_STABLE(4,:),'.b','MarkerSize',11);hold on;
% plot(STORE_POINCARE_UNSTABLE(2,:),STORE_POINCARE_UNSTABLE(4,:),'.r','MarkerSize',11)
% xline(R,'--');
% xline(-R,'--')

%%%%%%%%%%%%%% FLI mapping %%%%%%%%%%%%%%%%%%
% % % % % % % % % % [T_unst, M_unst, t_unst, X_unst, i_unst] = ode15s(@(t,x) PCC4BP_eqn(t,x,mu1,mu2,[(a_3/a_syn)+mu1,0], n_3, 3), tspan_unst, [Xunst0(1+4*22:4+4*22);0], opt);
% % % % % % % % % % %[T_st, M_st, t_st, X_st, i_st] = ode15s(@(t,x) PCC4BP_eqn(t,x,mu1,mu2,[(a_3/a_syn)+mu1,0], n_3, 3), tspan_st, [Xst0(1+4*33:4+4*33);0], opt);
% % % % % % % % % % X0 = Xunst0(1+4*22:4+4*22);
% % % % % % % % % % %X0 = Xst0(1+4*33:4+4*33);
% % % % % % % % % % 
% % % % % % % % % % %%% perturbed initial condition on manifolds  (Callisto --> eigvec1, Europa --> eigvec2)
% % % % % % % % % % pert_perpend = cross([0 0 1], [eigvec1(1:2); 0]);
% % % % % % % % % % pert_perpend4D = 1e-2*[pert_perpend(1:2)'; 0; 0];
% % % % % % % % % % X_pert0 = [X0+pert_perpend4D; 0];
% % % % % % % % % % [Tpert_unst, Mpert_unst, tpert_unst, Xpert_unst, ipert_unst] = ode15s(@(t,x) PCC4BP_eqn(t,x,mu1,mu2,[(a_3/a_syn)+mu1,0], n_3, 3), tspan_unst, X_pert0, opt);
% % % % % % % % % % 
% % % % % % % % % % delxt = M_unst-Mpert_unst;
% % % % % % % % % % LCE = [];
% % % % % % % % % % time = [];
% % % % % % % % % % for i=1:length(delxt)
% % % % % % % % % %     if Tpert_unst(i) ~= 0
% % % % % % % % % %         time = [time abs(Tpert_unst(i))];
% % % % % % % % % %         LCE = [LCE log(norm(delxt(i,1:2))/norm(pert_perpend4D(1:2)))/time(end)];
% % % % % % % % % %     end
% % % % % % % % % % end
% % % % % % % % % % [a_err,index] = min(abs(time*TU/(24*3600)-60.9131)); % unst - 60.9131   st - 73.2643
% % % % % % % % % % figure(2);plot(time(1:end)*TU/(24*3600), LCE(1:end),'-r','LineWidth',1.5);hold on;plot(time(index+1:end)*TU/(24*3600), LCE(index+1:end), '--','Color', [0,0.4,0.2],'LineWidth',1.5);  %orange shade-[0.9,0.4,0]  blue shade-[0,0.4,0.2]
% % % % % % % % % % figure(4);plot(M_unst(:,1),M_unst(:,2),'-r',Mpert_unst(:,1),Mpert_unst(:,2),'-b')



















%%%%%%%%%%%%%%%%%%%%%%%%%%% Custom functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J_lmda, J, G] = BVP(x,tau,tau_diff,m,N2,extras,D_rho, theta_init)
% x --> collected state between i to i+1 collocation (4*m); last element x(i+1,0)  (size --> 4*N2*(m+2) vector)
% tau_diff --> tau(i,0) to tau(i+1,0)
% NOTE (imp)  --> 2*pi*tau_j is the angle of ganymede at tau_j parameter value

rho_l1 = extras(1);
lmda = extras(2);
n_3 = extras(4);
mu1 = extras(5);
mu2 = extras(6);
l1_pos_io = extras(7);
w0 = extras(8);
T = 2*pi/abs(n_3-1);

theta_0 = theta_init; %Initial phase of Ganymede (exists in func:BVP and func:F_coll)

J_pert = [zeros(2) eye(2); -eye(2) zeros(2)];
dydx = [1 0 0 1; 0 1 -1 0; 0 0 1 0; 0 0 0 1]';
S_pert = kron(eye(N2), -dydx'*J_pert*dydx);

G = []; % Main constriant vector for Newton iteration
J = zeros(4*N2*(m+1), 4*N2*(m+2));  % Main Jacobian vector for Newton iteration (represents block (b2) in fig 2.1 Olikara)
syms p;
roots = double(vpasolve(legendreP(m,p))); %(tau_caps) (size --> m) (-1,1)
roots = [-1; roots];  % size m+1 [-1,1)
J_lmda = [];
% collocation constraint section 
for j=1:m  % iterating over collocation points btw a single tau interval
    % Constructing collocation constraint vector (1 to m)
    tau_j = tau + 0.5*(roots(j+1)+1)*tau_diff;  % preparing for F_coll (RESCALING)
    x_dyn = x(1+4*N2*j:4*N2+4*N2*j); % This should only be used for J_lmda, F_coll and PCC4BP (4*N2x1) (Contains Ganymede's current position and required for stroboscopic mapping in F_coll)
    f = F_coll(x_dyn,[tau_j,rho_l1*tau_j],lmda,[extras(3),0],[n_3, mu1, mu2, 4, l1_pos_io], theta_0);
    S = 0;
    dvdtheta = (1/w0)*(f-(n_3-1)*real(D_rho*x_dyn));
    %J_lmda = [J_lmda; 0.5*tau_diff*S_pert*vf_pert(x_dyn, 2*pi*tau_j, l1_pos_io)];
    J_lmda = [J_lmda; 0.5*tau_diff*S_pert*T*dvdtheta];
    for k=1:m+1
        % define lagrange basis polynomial
        [diff_lagrange, lagr] = D_lagr(roots(j+1),roots,k);
        S = S + (x(1+4*N2*(k-1):4*N2+4*N2*(k-1))*diff_lagrange);
        
        if (k-1)~=j
            % Jacobian collocation section (includes off diagonals) % (2 to m+1)
            J(1+4*N2*(j-1):4*N2+4*N2*(j-1), 1+4*N2*(k-1):4*N2+4*N2*(k-1)) = diff_lagrange*eye(4*N2);
        else
            % Jacobian collocation section (includes diagonals) % (2 to m+1)    
            JJ_block = [];
            for i=1:N2
                J_block = PCC4BP_J(x_dyn(1+4*(i-1)), x_dyn(2+4*(i-1)), theta_0 + sign(n_3-1)*2*pi*tau_j, mu1, mu2, extras(3));
                if length(JJ_block) ~= 0        
                    JJ_block = [JJ_block zeros(length(JJ_block(:,1)),4); zeros(4,length(JJ_block(1,:))) J_block(1:4, 1:4)];
                else
                    JJ_block = J_block(1:4, 1:4);
                end
            end
            %J(1+4*N2*(j-1):4*N2+4*N2*(j-1), 1+4*N2*j:4*N2+4*N2*j) = diff_lagrange*eye(4*N2) - 0.5*tau_diff*(JJ_block - lmda*S_pert*vf_pert(x(1+4*N2*j:4*N2+4*N2*j), 2*pi*tau_j, l1_pos_io));
            J(1+4*N2*(j-1):4*N2+4*N2*(j-1), 1+4*N2*j:4*N2+4*N2*j) = diff_lagrange*eye(4*N2) - 0.5*tau_diff*(JJ_block - lmda*S_pert*dvdtheta);
        end
    end
    g = S - 0.5*tau_diff*f;  % 4x1
    G = [G; g];
end
% continuity constraint section
S = 0;
for j=0:m
    % Constructing continuty constraint
    [deriv_lagr, lagr] = D_lagr(1,roots,j+1);
    S = S + x(1+4*N2*j:4*N2+4*N2*j)*lagr;
end
G = [G; x(end-(2*4*N2-1):end-(4*N2))-x(end-(4*N2-1):end)];
%G = [G; S-x(end-3:end)];
% Jacobian continuity section % (X(i+1,0))
J(end-(4*N2-1):end, end-(2*4*N2-1):end-(4*N2)) = eye(4*N2);
J(end-(4*N2-1):end, end-(4*N2-1):end) = -eye(4*N2);
end






function [D_L,Q] = D_lagr(tau,roots,k)
m = length(roots);
syms x;
S = 0;

% Computing derivative
for i=1:m
    if i~=k
       L = 1;
       for j=1:m
           if j~=i && j~=k
               L = L*(x-roots(j))/(roots(k)-roots(j));
           end
       end
       S = S + (1/(roots(k)-roots(i)))*L;
    end 
end

% Computing basis function
Q = 1;
for i=1:m
    if i~=k
        Q = Q*(x-roots(i))/(roots(k)-roots(i));
    end
end

x = tau; 
D_L = double(subs(S));
Q = double(subs(Q));
end






function DV = vf_pert(X, theta2, l1_pos_io)
N2 = length(X)/4;
DV = [];
for i=1:N2
    XX = X(1+4*(i-1):4+4*(i-1));
    x = XX(1);
    y = XX(2);
    
    l1 = 2*y*(csc(theta2))^2*cot(theta2)*(y/(x^2 + y^2));
    l2 = (csc(theta2))^2 - 2*y*(csc(theta2))^2*cot(theta2)*(x/(x^2 + y^2));
    l3 = (sec(theta2))^2 - 2*(x-l1_pos_io)*(sec(theta2))^2*tan(theta2)*(y/(x^2 + y^2));
    l4 = 2*(x-l1_pos_io)*(sec(theta2))^2*tan(theta2)*(x/(x^2 + y^2));
    
    DV = [DV; l1; l2; l3; l4];
end
end




function X = get_nearby_orbit(pos,del_x,eigvec,extras)
tStep = 1*10^(-3);
mu1 = extras(1); a_2 = extras(2); a_3 = extras(3);
x0 = [pos;0;0;0] + del_x*eigvec(:,3);
x0 = [x0; 0];                                                                  % The fifth state is the angle of third primary with first 
[a_l1, b_l1, STM_b] = diffcorr_l1_FBP(x0, mu1, 0, [(a_3/a_2)-mu1,0], 0, tStep);      % a --> initial state corrected   2*b --> time period of the periodic orbit  STM_b --> STM matrix from 0 to b
[t, X] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,0,[(a_3/a_2)+mu1,0], 0, 3), 0:0.001:2*b_l1, a_l1);
X(:,5) = x0(5)*ones([length(X(:,5)),1]);  % third primary angle for all the entries of periodic orbit should have same value as the initial condition x0
end

function [position,isterminal,direction] = ps_crossing(t,y)
  global mu1;
  position = y(2); % The value that we want to be zero --- Poincare section existing in the y=0 line and x>0
  isterminal = 0;  % Halt integration 
  direction = 1;   % The zero can be approached from either direction -- moving towards positive y direction
end
