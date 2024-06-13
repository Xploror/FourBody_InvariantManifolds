clc;clear;

%%%% Units used in this script are important to understand --- firstly
%%%% there is a distance scaling were ganymede's semi major is taken as 5
%%%% times this canonical unit, secondly there is a mass scaling involved
%%%% where 10^(20)kg is the mass canonical unit

%%%% Functions in use : crtbp_J, diffcorr_l1_FBP, diffcorr_l2_FBP, PCC4BP_eqn

tStep = 1*10^(-3);

MU = 10^20;         %kg
DU = 1070400000/5;  %m
TU = 7.155*24*3600; %s

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
del_x = 3*10^(-3);
tspan = 0:0.001:T_p;
x0 = [l1_pos_io;0;0;0] + del_x*l1_eigvec_io(:,3);
x0 = [x0; 0];                                                                  % The fifth state is the angle of third primary with first 
[a_l1, b_l1, STM_b] = diffcorr_l1_FBP(x0, mu1, 0, [(a_3/a_2)-mu1,0], n_3, tStep);      % a --> initial state corrected   2*b --> time period of the periodic orbit  STM_b --> STM matrix from 0 to b
[t_dc, Xl1_dc] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,0,[(a_3/a_2)-mu1,0], n_3, 3), 0:0.001:2*b_l1, a_l1);
Xl1_dc(:,5) = x0(5)*ones([length(Xl1_dc(:,5)),1]);  % third primary angle for all the entries of periodic orbit should have same value as the initial condition x0 
w_l1 = 2*pi/(2*b_l1); 
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
% del_x = -1*10^(-3);
% tspan = 0:0.001:T_p;
% x0 = [l2_pos_io;0;0;0] + del_x*l2_eigvec_io(:,3);
% x0 = [x0; 0];                                                                  % The fifth state is the angle of third primary with first 
% [a_l2, b_l2, STM_b] = diffcorr_l2_FBP(x0, mu1, 0, [(a_3/a_2)-mu1,0], n_3, tStep);      % a --> initial state corrected   2*b --> time period of the periodic orbit  STM_b --> STM matrix from 0 to b
% [t_dc, Xl2_dc] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,0,[(a_3/a_2)-mu1,0], n_3, 3), 0:0.001:2*b_l2, a_l2);
% Xl2_dc(:,5) = x0(5)*ones([length(Xl2_dc(:,5)),1]);  % third primary angle for all the entries of periodic orbit should have same value as the initial condition x0 
% w_l2 = 2*pi/(2*b_l2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% syms t
% figure(1)
% hold on;
% plot(Xl1_dc(:,1), Xl1_dc(:,2),'-k'); hold on;
% plot(Xl2_dc(:,1), Xl2_dc(:,2),'-k'); hold on;
% fplot((1-mu1)+R_europa*sin(t), R_europa*cos(t));
% %hold on;fplot(((a_3/a_2)-mu1)*sin(t),((a_3/a_2)-mu1)*cos(t),'r','Linewidth',4)
% pause(0.000000001)
% hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTATION OF PERIODIC ORBIT OVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTATION OF QUASI-PERIODIC ORBITS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dphii_l1, phii_l1, x_l1] = get_STM(t_dc(end),Xl1_dc(1,:)',4,1);  % Taking out STM for periodic orbit from TBP but now calculating it for FBP
[rvec, rho] = eig(phii_l1);
rho_l1 = angle(rho(3,3));
%rho_l2 = w_l2/(n_3-1);

theta_l1 = parameterize(Xl1_dc,l1_pos_io,'2D');
%theta_l2 = parameterize(Xl2_dc,l2_pos_io,'2D');
%%%%%%% Taking periodic orbit from 3BP as an initial guess for Newton's method, find invariant cicle for the stroboscopic map

%%%%%%%% DFT of initial states section %%%%%%%%
N2 = 49; % Total number of discrete points on invariant circle (KEEP IT ODD)
m2 = (N2-1)/2;  % FS truncation order
k = -(N2-1)/2:(N2-1)/2;

% U = [];
% for i=1:length(t_dc)
%     PHI = get_STM(2*b_l1,Xl1_dc(i,:)',4,1);
%     [evec,eval] = eig(PHI);
%     uu = real(exp(1i*theta_l1(i))*evec(:,3));
%     U = [U uu/norm(uu)];
% end
U = 0;
U_exp = Xl1_dc(:,1:4) + 1*10^(-5)*U';  % Linear approximation of invariant curve (length(t_dc),5)

% Creating initial Family Tangent vector
del = 1.5*del_x;
X_tang = get_nearby_orbit(l1_pos_io,del,l1_eigvec_io,[mu1,a_2,a_3]);

X_coll0 = [];
T_coll0 = [];
X_tang0 = [];
THETA = [];
for i=1:N2
    angle = ((i-1)/N2)*2*pi;
    [val, idx] = min(abs(theta_l1-angle));
    THETA = [THETA; theta_l1(idx)];
    X_coll0 = [X_coll0; U_exp(idx,1:4)'];    % Size n*N2
    T_coll0 = [T_coll0; t_dc(idx)];
    X_tang0 = [X_tang0; X_tang(idx,1:4)'];   % Collected equivalent of family tangent vector
end
X_tang_vec = X_coll0 - X_tang0;
w_coll0 = 2*pi/t_dc(end);  % For phase condition

% creating DFT and inverse DFT matrix
DFT = zeros(length(X_coll0),length(X_coll0));
Q = zeros(length(X_coll0),length(X_coll0));
w = exp(-1i*2*pi/N2);
for i=1:N2
    Q(1+4*(i-1):4+4*(i-1),1+4*(i-1):4+4*(i-1)) = exp(-1i*(-0.5*(N2-1)+(i-1))*rho_l1)*eye(4);
    for j=1:N2
        DFT(1+4*(i-1):4+4*(i-1),1+4*(j-1):4+4*(j-1)) = w^(-(j-1)*(0.5*(N2-1) - (i-1)))*eye(4);
    end
end
KKK = sparse(DFT'*Q*DFT/N2);  % 4*N2 square  (Lujan)
D_rho = DFT'*(kron(1i*diag(1:N2),eye(4)))*DFT/N2;  % 4*N2 square  (Lujan)
X_dtheta = real(D_rho*X_coll0);   % 4*N2 vector
LL = sparse(real(D_rho*KKK));     % 4*N2 square
%%%%%%%%%%%%%%%%%%%%%% Gauss Legendre section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmda = 1e-14; % Just an initial guess
psued_arc = 1e-8;

%%% (For now) Use Newton method only for one trajectory (i.e only one vector from collected state vector)
m = 2;  % Order of legendre polynomial (SHOULD BE EVEN)
N = 10; %50 % Number of collocation points
newton_iter = 6;
%%%%%%%%%%%% Newton's methods iterations %%%%%%%%%%%%%%%%%%
X_interval = kron(ones(N*(m+1)+1,1),X_coll0); % x_interval for each invariant point on orbit
X_interval = [X_interval; rho_l1]; % THIS VARIABLE NEED TO BE CONVERGED DURING NEWTON'S METHOD (HAS ALL THE INTERVAL) TOTAL SIZE --> N2*(4*N*(m+1)+4)+1 SQUARE
G_mat = [];
len_circle = 4*N2;   % Total number of discrete invariant circle
len_coll = N*(m+1);   % Total number of collocation + continuity constraint
    for h=1:newton_iter
               J_total = zeros((len_circle*len_coll+len_circle),(len_circle*len_coll+len_circle)+1);    % REFER THE SIZE INFO FROM OLIKARA FIG.2.1 (Added coll+cont+quasiperiod+phase | Xij+Xn0+rho)
               G_total = [];
               Rot = real(KKK);
               X_init = X_interval(1 : 4*N2);
               X_ends = X_interval(end-(4*N2) : end-1);
               x_rot = Rot*X_ends;
               % Adding collocation + continuity condition
               for i=1:N
                   tau_interval = [(i-1)/N i/N];
                   X_ = X_interval(1+len_circle*(m+1)*(i-1):(len_circle*(m+2))+len_circle*(m+1)*(i-1));
                   [J, G] = BVP(X_, tau_interval(1), 1/N, m, N2, [rho_l1, lmda, ((a_3/a_2)-mu1), n_3, mu1, mu2, l1_pos_io, b_l1]);
                   J_total(1+len_circle*(m+1)*(i-1):(len_circle*(m+1))+len_circle*(m+1)*(i-1), 1+len_circle*(m+1)*(i-1):(len_circle*(m+2))+len_circle*(m+1)*(i-1)) = J;   % N*(4*N2*(m+1)) X N*(4*N2*(m+1)+4)
                   G_total = [G_total; G]; 
               end

               % Adding QuasiPeriodic Condition
               G_total = [G_total; X_init-x_rot];  % 4*N2
               J_total(end-(4*N2)+1:end, 1:4*N2) = eye(4*N2); J_total(end-(4*N2)+1:end, end-(4*N2):end-1) = -Rot;    % d/d(X_N,0)
               
               J_total(1:4*N2, end) = LL*X_ends;

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %%%%%%%%% Adding Phase Condition
               % fff = F_coll(X_coll0,[0,0],lmda,[(a_3/a_2)-mu1,0],[n_3, mu1, mu2, 4, l1_pos_io, 2*b_l1]);  % For phase condition
               % Xdiff = (X_init-X_coll0);
               % X_dtheta1 = (1/w_coll0)*(fff - (n_3-1)*X_dtheta);
               % G_total = [G_total; dot(X_init, X_dtheta1)];  
               % J_total(end, 1:4*N2) = X_dtheta1';
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %%%%%%%%% Adding Pseudo-arclength Condition
               % G_total = [G_total; (X_init-X_tang0)'*X_tang_vec - psued_arc];
               % J_total(end, 1:4*N2) = X_tang_vec';
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               if norm(G_total)<1e-12
                    fprintf('Break %i',h)
                    break
               end

               G_mat = [G_mat norm(G_total)];
               X_interval = X_interval - pinv(J_total)*G_total;
               %%%% Update tangent vector for use in pseudo-arclength and phase condition
               %X_tang_vec = X_init - X_tang0;
               X_dtheta = real(D_rho*X_init);  % X_dtheta extracted from previous newton stage X_init
               for i=1:N2
                    Q(1+4*(i-1):4+4*(i-1),1+4*(i-1):4+4*(i-1)) = exp(-1i*(-0.5*(N2-1)+(i-1))*X_interval(end))*eye(4);
               end
               KKK = sparse(DFT'*Q*DFT/N2);
               LL = sparse(real(D_rho*KKK));

               plot_collected(X_interval(1:4*N2),'position',[1 0 0]);
               plot_collected(X_interval(5*4*N2 + 1: 5*4*N2 + 4*N2),'position',[0 0 1]);
    end
    col = rand(1,3);
    plot_collected(X_interval(1:4*N2), 'position', col);
    pause(0.01)
%     plot(X_interval(1),X_interval(2),'.r',X_interval(end-3),X_interval(end-2),'.b')
%     hold on;
%     pause(0.000000001)
    %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%% Custom functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J, G] = BVP(x,tau,tau_diff,m,N2,extras)
% x --> collected state between i to i+1 collocation (4*m); last element x(i+1,0)  (size --> 4*N2*(m+2) ve)
% tau_diff --> tau(i,0) to tau(i+1,0)
% NOTE (imp)  --> 2*pi*tau_j is the angle of ganymede at tau_j parameter value

rho_l1 = extras(1);
lmda = extras(2);
n_3 = extras(4);
mu1 = extras(5);
mu2 = extras(6);
l1_pos_io = extras(7);
b_l1 = extras(8);
T = 2*pi/(n_3-1);

J_pert = [zeros(2) eye(2); -eye(2) zeros(2)];
dydx = [1 0 0 1; 0 1 -1 0; 0 0 1 0; 0 0 0 1]';
S_pert = kron(eye(N2), -dydx'*J_pert*dydx);

G = []; % Main constriant vector for Newton iteration
J = zeros(4*N2*(m+1), 4*N2*(m+2));  % Main Jacobian vector for Newton iteration (represents block (b2) in fig 2.1 Olikara)
syms p;
roots = double(vpasolve(legendreP(m,p))); %(tau_caps) (size --> m) (-1,1)
roots = [-1; roots];  % size m+1 [-1,1)
% collocation constraint section 
for j=1:m  % iterating over collocation points btw a single tau interval
    % Constructing collocation constraint vector (1 to m)
    tau_j = tau + 0.5*(roots(j+1)+1)*tau_diff;  % preparing for F_coll (RESCALING)
    x_dyn = x(1+4*N2*j:4*N2+4*N2*j); % This should only be used for F_coll and PCC4BP (4*N2x1) (Contains Ganymede's current position and required for stroboscopic mapping in F_coll)
    f = F_coll(x_dyn,[tau_j,rho_l1*tau_j],lmda,[extras(3),0],[n_3, mu1, mu2, 4, l1_pos_io, 2*b_l1]);
    S = 0;
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
                J_block = PCC4BP_J(x_dyn(1+4*(i-1)), x_dyn(2+4*(i-1)), 2*pi*tau_j, mu1, mu2, extras(3));
                if length(JJ_block) ~= 0
                    JJ_block = [JJ_block zeros(length(JJ_block(:,1)),4); zeros(4,length(JJ_block(1,:))) J_block(1:4, 1:4)];
                else
                    JJ_block = J_block(1:4, 1:4);
                end
            end
            J(1+4*N2*(j-1):4*N2+4*N2*(j-1), 1+4*N2*j:4*N2+4*N2*j) = diff_lagrange*eye(4*N2) - 0.5*tau_diff*(JJ_block - lmda*S_pert*vf_pert(x(1+4*N2*j:4*N2+4*N2*j), 2*pi*tau_j, l1_pos_io));
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
    x = X(1);
    y = X(2);
    
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
[t, X] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,0,[(a_3/a_2)-mu1,0], 0, 3), 0:0.001:2*b_l1, a_l1);
X(:,5) = x0(5)*ones([length(X(:,5)),1]);  % third primary angle for all the entries of periodic orbit should have same value as the initial condition x0
end
