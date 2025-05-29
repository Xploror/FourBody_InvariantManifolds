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

%%%%%% Loading tori #####################
load('good plots\plots_for_paper\3_l1.mat');
load('good plots\plots_for_paper\3_l1data.mat');
Xl1_dc = X_interval;
Jtot_l1 = J_total;
Gtot_l1 = G_total;

load('good plots/plots_for_paper/3_l2.mat');
load('good plots/plots_for_paper/3_l2data.mat');
Xl2_dc = X_interval;
Jtot_l2 = J_total;
Gtot_l2 = G_total;

N2_l1 = 35; N2_l2 = 35; %35 29
N2 = N2_l1;  % Assuming N2_l1 and N2_l2 are same
len_circle = 4*N2_l1;
len_coll = 8*(4+1); % PLEASE store N and m values when calculating QP tori
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
w = exp(-1i*2*pi/N2);
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

% Stable/unstable manifold calculation
% L2 stable and L1 unstable manifold
epsilon = 1e-6;
Xst0 = Xl2_dc(1:4*N2_l2) + epsilon*real(torus_vec2(:,5)); % Stable direction (mag less than 1)  (end-15)
Xunst0 = Xl1_dc(1:4*N2_l1) + epsilon*real(torus_vec1(:,29));  % Unstable direction (mag more than 1)  (29)

b = 1; %Choosing index of the point on the l1 circle
d = 3; %Choosing index of the point on the l2 circle
w0_u = Xunst0(1+4*b:4+4*b);
w1_s = Xst0(1+4*d:4+4*d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hard coded variabls for Jovian system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_l1 = 0.164;
l1_pos_io = 0.9798; %J-E
l2_pos_io = 1.0205; %J-E
l1_pos_cal = 0.9736; %J-C
l2_pos_cal = 1.0268; %J-C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Initial computations %%%%%%%%%%%%%%%%%%%%%%%%%
D_rho = DFT'*(kron(1i*diag(-(N2-1)/2:(N2-1)/2),eye(4)))*DFT/N2;
Xl1_dtheta = real(D_rho*Xl1_dc(1:4*N2_l1));
Xl2_dtheta = real(D_rho*Xl2_dc(1:4*N2_l2));
wl1 = Xl1_dc(end-1)*(n_3-1);
wl2 = Xl2_dc(end-1)*(n_3-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTATION OF HETEROCLINIC CONNECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Gauss Legendre section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmda = 1e-14; % Just an initial guess
psued_arc = 1e-12;

%%% (For now) Use Newton method only for one trajectory (i.e only one vector from collected state vector)
N2 = 1; % Number of points in the circle for which connections has to be calculated
m = 9;  % Order of legendre polynomial (SHOULD BE EVEN)
N = 21; %50 % Number of collocation points
newton_iter = 6;
%%%%%%%%%%%% Newton's methods iterations %%%%%%%%%%%%%%%%%%

T_span = pi/abs(n_3-1); % T_span --> time taken from L1 to L2
[ttt, Xxx] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,mu2,[(a_3/a_2)-mu1,0], n_3, 3), 0:0.001:T_span, w0_u);      %%% generating initial guess manifold %%%%%
Xxx = reshape(Xxx(1:N*(m+1)+1,:)', [4*(N*(m+1)+1),1]);

fffl1 = F_coll(Xl1_dc(1:4*N2_l1),[0,0],lmda,[(a_3/a_2)-mu1,0],[n_3, mu1, mu2, 4, l1_pos_io]);
fffl2 = F_coll(Xl2_dc(1:4*N2_l2),[0,0],lmda,[(a_3/a_2)-mu1,0],[n_3, mu1, mu2, 4, l2_pos_io]);
ZZZun = (1/wl1)*(fffl1 - (n_3-1)*Xl1_dtheta);
ZZZs = (1/wl2)*(fffl2 - (n_3-1)*Xl2_dtheta);

b_desir = 1;
X_0 = Xl1_dc(1:4); % State vector 4x1 (initial state in the quasi-periodic tori)
u0 = 0; u1= 0; 
X_var = Xxx; % x_var size --> 4*(N*(m+1)+1)
X_var = [X_var; u0; b; u1; d; T_span]; % THIS VARIABLE NEED TO BE CONVERGED DURING NEWTON'S METHOD (HAS ALL THE INTERVAL) TOTAL SIZE --> 4*(N*(m+1)+1)+2  FORMAT -> [Xi,j, u0, b, u1, d, T, lmda]
G_mat = [];
len_J_total = 4*N*(m+1);  % Total number of collocation + continuity constraint
    for h=1:newton_iter
               J_total = zeros((4*N*(m+1)+4)+5,(4*N*(m+1)+4)+5);    % REFER THE SIZE INFO FROM OLIKARA FIG.2.1 (Added coll+cont+boundaryequality+pseudo-arc | Xij+Xn0+T_span+lambda)
               G_total = [];

               % X_ends = []; X_init = [];
               % for j=1:N2
               %     X_init = [X_init; X_var(1 + (4*(m+1)*N+4)*(j-1): 4 + (4*(m+1)*N+4)*(j-1))];
               %     X_ends = [X_ends; X_var(4*(m+1)*N+1 + (4*(m+1)*N+4)*(j-1): 4*(m+1)*N+4 + (4*(m+1)*N+4)*(j-1))];
               % end
               % Adding collocation + continuity condition
               for j=1:N2
                   X_int = X_var(1+(len_J_total+4)*(j-1):(len_J_total+4)+(len_J_total+4)*(j-1));  % Takes entire trajectory of N2 number of points
                   for i=1:N   % iterate over equally spaced collocation times
                       tau_interval = [(i-1)/N i/N];
                       XX = X_int(1+4*(m+1)*(i-1):(4*(m+2))+4*(m+1)*(i-1));
                       [J_lmda, J, G] = BVP(XX, tau_interval(1), 1/N, m, [rho_l1, lmda, ((a_3/a_2)-mu1), n_3, mu1, mu2, l1_pos_io, X_var(end)]);
                       J_total(1+4*(m+1)*(i-1):(4*(m+1))+4*(m+1)*(i-1), 1+4*(m+1)*(i-1):(4*(m+2))+4*(m+1)*(i-1)) = J;   % N2*(4*N*(m+1)) X N2*(4*N*(m+1)+4)
                       J_total(1+4*(m+1)*(i-1):(4*(m+1))+4*(m+1)*(i-1), end-1) = -0.5*(1/N)*F_coll(XX(1:end-3),[tau_interval(1),rho_l1*tau_interval(1)],lmda,[((a_3/a_2)-mu1),0],[n_3, mu1, mu2, 4, l1_pos_io]);
                       %J_total(1+4*(m+1)*(i-1):(4*(m+1))+4*(m+1)*(i-1), end) = J_lmda;
                       G_total = [G_total; G]; 
                   end         
               end

               %%%% Boundary conditions to lie on st/unst manifolds %%%%%%
               G_total = [G_total; X_var(1:4)-w0_u];
               J_total(4*N*(m+1)+1:4*N*(m+1)+4, 1:4) = eye(4);  % 4x4
               J_total(4*N*(m+1)+1:4*N*(m+1)+4, 4*N*(m+1)+4+1) = -(Xl1_dtheta(1+4*X_var(end-3):4+4*X_var(end-3)));  % 4x1
               J_total(4*N*(m+1)+1:4*N*(m+1)+4, 4*N*(m+1)+4+2) = -(ZZZun(1+4*X_var(end-3):4+4*X_var(end-3)));  % 4x1
            
               G_total = [G_total; X_var(end-9:end-6)-w1_s];
               J_total(4*N*(m+1)+5:4*N*(m+1)+8, 4*N*(m+1)+1:4*N*(m+1)+4) = eye(4);  % 4x4
               J_total(4*N*(m+1)+5:4*N*(m+1)+8, 4*N*(m+1)+4+3) = -(Xl2_dtheta(1+4*X_var(end-1):4+4*X_var(end-1)));  % 4x1
               J_total(4*N*(m+1)+5:4*N*(m+1)+8, 4*N*(m+1)+4+4) = -(ZZZs(1+4*X_var(end-1):4+4*X_var(end-1)));  % 4x1
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               %%%%%%%%%%% h(z) = 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               G_total = [G_total; X_var(end-3)-b_desir];
               J_total(4*N*(m+1)+9, 4*N*(m+1)+4+2) = 1;
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               %%%%%%%%% Adding Pseudo-arclength Condition
               % G_total = [G_total; 0];
               % for i=1:N2
               %      J_TOT(end, 1 + (4*(m+1)*N+4)*(i-1):4 + (4*(m+1)*N+4)*(i-1)) = X_tang_vec(1+4*(i-1):4+4*(i-1))';
               % end
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               G_mat = [G_mat norm(G_total)];
               X_var = X_var - pinv(J_total)*G_total;
               
               if norm(G_total)<1e-15
                    fprintf('Break %i',h)
                    break
               end

               w0_u = Xunst0(1+4*int32(X_var(end-3)):4+4*int32(X_var(end-3)));
               w1_s = Xst0(1+4*int32(X_var(end-1)):4+4*int32(X_var(end-1)));

    end
    col = rand(1,3);
    plot_collected(X_var(1:end-5), 'position', col);
    pause(0.01)
    %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%% Custom functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J_lmda, J, G] = BVP(x,tau,tau_diff,m,extras)
% x --> collected state between i to i+1 collocation (4*m); last element x(i+1,0)
% tau_diff --> tau(i,0) to tau(i+1,0)
% NOTE (imp)  --> 2*pi*tau_j is the angle of ganymede at tau_j parameter value

rho_l1 = extras(1);
lmda = extras(2);
n_3 = extras(4);
mu1 = extras(5);
mu2 = extras(6);
l1_pos_io = extras(7);
T = extras(8);

J_pert = [zeros(2) eye(2); -eye(2) zeros(2)];
dydx = [1 0 0 1; 0 1 -1 0; 0 0 1 0; 0 0 0 1]';
S_pert = -dydx'*J_pert*dydx;

G = []; % Main constriant vector for Newton iteration
J = zeros(4*(m+1), 4*(m+2));  % Main Jacobian vector for Newton iteration 
syms p;
roots = double(vpasolve(legendreP(m,p))); %(tau_caps) (size --> m) (-1,1)
roots = [-1; roots];  % size m+1 [-1,1)
J_lmda = -0.5*tau_diff*S_pert*vf_pert(x(1:4), 0, l1_pos_io);
% collocation constraint section 
for j=1:m  % iterating over collocation points btw a single tau interval
    % Constructing collocation constraint vector (1 to m)
    tau_j = tau + 0.5*(roots(j+1)+1)*tau_diff;  % preparing for F_coll (RESCALING)
    x_dyn = [x(1+4*j:4+4*j); tau_j*2*pi]; % This should only be used for F_coll (5x1) (Contains Ganymede's current position and required for stroboscopic mapping in F_coll)
    f = F_coll(x_dyn,[tau_j,rho_l1*tau_j],lmda,[extras(3),0],[n_3, mu1, mu2, 4, l1_pos_io]);
    f = f(1:4);
    S = 0;
    dvdtheta = vf_pert(x(1+4*j:4+4*j), 2*pi*tau_j, l1_pos_io);
    J_lmda = [J_lmda; -0.5*tau_diff*S_pert*dvdtheta];
    for k=1:m+1
        % define lagrange basis polynomial
        [diff_lagrange, lagr] = D_lagr(roots(j+1),roots,k);
        S = S + (x(1+4*(k-1):4+4*(k-1))*diff_lagrange);
        
        if (k-1)~=j
            % Jacobian collocation section (includes off diagonals) % (2 to m+1)
            J(1+4*(j-1):4+4*(j-1), 1+4*(k-1):4+4*(k-1)) = diff_lagrange*eye(4);
        else
            % Jacobian collocation section (includes diagonals) % (2 to m+1)
            J_block = PCC4BP_J(x(1+4*j), x(2+4*j), 2*pi*tau_j, mu1, mu2, extras(3));
            J(1+4*(j-1):4+4*(j-1), 1+4*j:4+4*j) = diff_lagrange*eye(4) - 0.5*tau_diff*(J_block(1:4, 1:4) - lmda*S_pert*vf_pert(x(1+4*j:4+4*j), 2*pi*tau_j, l1_pos_io));
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
    S = S + x(1+4*j:4+4*j)*lagr;
end
G = [G; x(end-7:end-4)-x(end-3:end)];
%G = [G; S-x(end-3:end)];
% Jacobian continuity section % (X(i+1,0))
J(end-3:end, end-(4*2)+1:end-(4*2)+4) = eye(4);
J(end-3:end, end-3:end) = -eye(4);
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
