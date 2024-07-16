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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTATION OF HETEROCLINIC CONNECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Gauss Legendre section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmda = 1e-14; % Just an initial guess
psued_arc = 1e-12;

%%% (For now) Use Newton method only for one trajectory (i.e only one vector from collected state vector)
m = 15;  % Order of legendre polynomial (SHOULD BE EVEN)
N = 51; %50 % Number of collocation points
newton_iter = 6;
%%%%%%%%%%%% Newton's methods iterations %%%%%%%%%%%%%%%%%%
X_0 = 0; % State vector 4x1 (initial state in the quasi-periodic tori)
T_span = 0; % T_span --> time taken from L1 to L2 
X_var = [];
for j=0:N2-1
    X_var = [X_var; [kron(ones(N,1),kron(ones(m+1,1),X_invcir)); X_0]]; % x_interval for each invariant point on orbit
end
X_var = [X_var; T_span, lmda]; % THIS VARIABLE NEED TO BE CONVERGED DURING NEWTON'S METHOD (HAS ALL THE INTERVAL) TOTAL SIZE --> N2*(4*N*(m+1)+4)+1 SQUARE
G_mat = [];
len_J_total = 4*N*(m+1);  % Total number of collocation + continuity constraint
    for h=1:newton_iter
               J_total = zeros((4*N*(m+1)+4),(4*N*(m+1)+4));    % REFER THE SIZE INFO FROM OLIKARA FIG.2.1 (Added coll+cont+quasiperiod+phase | Xij+Xn0+rho)
               G_total = [];

               X_ends = []; X_init = [];
               for j=1:N2
                   X_init = [X_init; X_var(1 + (4*(m+1)*N+4)*(j-1): 4 + (4*(m+1)*N+4)*(j-1))];
                   X_ends = [X_ends; X_var(4*(m+1)*N+1 + (4*(m+1)*N+4)*(j-1): 4*(m+1)*N+4 + (4*(m+1)*N+4)*(j-1))];
               end
               % Adding collocation + continuity condition
               for j=1:N2
                   X_int = X_var(1+(len_J_total+4)*(j-1):(len_J_total+4)+(len_J_total+4)*(j-1));
                   for i=1:N
                       tau_interval = [(i-1)/N i/N];
                       X_ = X_int(1+4*(m+1)*(i-1):(4*(m+2))+4*(m+1)*(i-1));
                       [J, G] = BVP(X_, tau_interval(1), 1/N, m, [rho_l1, lmda, ((a_3/a_2)-mu1), n_3, mu1, mu2, l1_pos_io, b_l1]);
                       %J_total(1+4*(m+1)*(i-1):(4*(m+1))+4*(m+1)*(i-1), 1+4*(m+1)*(i-1):(4*(m+2))+4*(m+1)*(i-1)) = J;                  % d/d(Xij)
                       J_total(1+4*(m+1)*(i-1):(4*(m+1))+4*(m+1)*(i-1), 1+4*(m+1)*(i-1):(4*(m+2))+4*(m+1)*(i-1)) = J;   % N2*(4*N*(m+1)) X N2*(4*N*(m+1)+4)
                       G_total = [G_total; G]; 
                   end         
               end

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %%%%%%%%% Adding Pseudo-arclength Condition
               % G_total = [G_total; (X_init-X_tang0)'*X_tang_vec - psued_arc];
               % for i=1:N2
               %      J_TOT(end, 1 + (4*(m+1)*N+4)*(i-1):4 + (4*(m+1)*N+4)*(i-1)) = X_tang_vec(1+4*(i-1):4+4*(i-1))';
               % end
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               G_mat = [G_mat norm(G_total)];
               X_var = X_var - pinv(J_TOT)*G_total;
               
               if norm(G_total)<1e-15
                    fprintf('Break %i',h)
                    break
               end
    end
    col = rand(1,3);
    plot_collected(X_interval(1:end-1), 'position', col);
    pause(0.01)
%     plot(X_interval(1),X_interval(2),'.r',X_interval(end-3),X_interval(end-2),'.b')
%     hold on;
%     pause(0.000000001)
    %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%% Custom functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J, G] = BVP(x,tau,tau_diff,m,extras)
% x --> collected state between i to i+1 collocation (4*m); last element x(i+1,0)
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
S_pert = -dydx'*J_pert*dydx;

G = []; % Main constriant vector for Newton iteration
J = zeros(4*(m+1), 4*(m+2));  % Main Jacobian vector for Newton iteration 
syms p;
roots = double(vpasolve(legendreP(m,p))); %(tau_caps) (size --> m) (-1,1)
roots = [-1; roots];  % size m+1 [-1,1)
% collocation constraint section 
for j=1:m  % iterating over collocation points btw a single tau interval
    % Constructing collocation constraint vector (1 to m)
    tau_j = tau + 0.5*(roots(j+1)+1)*tau_diff;  % preparing for F_coll (RESCALING)
    x_dyn = [x(1+4*j:4+4*j); tau_j*2*pi]; % This should only be used for F_coll (5x1) (Contains Ganymede's current position and required for stroboscopic mapping in F_coll)
    f = F_coll(x_dyn,[tau_j,rho_l1*tau_j],lmda,[extras(3),0],[n_3, mu1, mu2, 5, l1_pos_io, 2*b_l1]);
    f = f(1:4);
    S = 0;
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
            J(1+4*(j-1):4+4*(j-1), 1+4*j:4+4*j) = diff_lagrange*eye(4) - 0.5*tau_diff*T*(J_block(1:4, 1:4) - lmda*S_pert*vf_pert(x(1+4*j:4+4*j), 2*pi*tau_j, l1_pos_io));
        end
    end
    g = S - 0.5*tau_diff*T*f;  % 4x1
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
x = X(1);
y = X(2);

l1 = 2*y*(csc(theta2))^2*cot(theta2)*(y/(x^2 + y^2));
l2 = (csc(theta2))^2 - 2*y*(csc(theta2))^2*cot(theta2)*(x/(x^2 + y^2));
l3 = (sec(theta2))^2 - 2*(x-l1_pos_io)*(sec(theta2))^2*tan(theta2)*(y/(x^2 + y^2));
l4 = 2*(x-l1_pos_io)*(sec(theta2))^2*tan(theta2)*(x/(x^2 + y^2));

DV = zeros(4,4);
DV(1,1) = l1;DV(1,2) = l2;DV(2,1) = l3;DV(2,2) = l4;
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
