function Y_coll = F_coll(x_coll,parameters,lmda,x1t,extras)

% x_coll --> collected state vector corresponding to perioidic forcing value tau (size-->[n*N,1])
% parameters --> 2 elements --> (angle1 time mapped and scaled to [0,1]) and angle2_map signifying stroboscopic mapping 
% extras contains all the necessary other parameters not directly involved for the functionality of this method
% T1 --> periodic orbit time period    T2 --> periodic forcing time period (Ganymede)
omega = extras(1);mu1 = extras(2);mu2 = extras(3);n = extras(4);l1 = extras(5);
T2 = 2*pi/abs(omega-1);

tau = parameters(1);
angle1 = (omega-1)*tau*T2;

N = length(x_coll)/n;   % DFT points
if n==4
    J = [zeros(2) eye(2); -eye(2) zeros(2)];
    dydx = [1 0 0 1; 0 1 -1 0; 0 0 1 0; 0 0 0 1]';
    S = -dydx'*J*dydx;
elseif n==5
    J = [[zeros(2) eye(2); -eye(2) zeros(2); zeros(1,4)], zeros(5,1)];
    dydx = [1 0 0 1 0; 0 1 -1 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 0]';
    S = -dydx'*J*dydx;
end

Y_coll = [];
for i=1:N   % iterating over DFT points in invariant circle
   x = x_coll(1+n*(i-1):n+n*(i-1));  %5x1
   angle2 = atan2(x(2),x(1)-l1);
   if n==4
       x = [x; angle1];   % If x_coll doesnt have perturber angle included, add it here for each of DFT point (N2)
       dvdx = [x(2)/(sin(angle2))^2 (x(1)-l1)/(cos(angle2))^2 0 0]; %1x4
   elseif n==5
       dvdx = [x(2)/(sin(angle2))^2 (x(1)-l1)/(cos(angle2))^2 0 0 0]; %1x5
   end
   %[tt, yy] = ode45(@(t,x) PCC4BP_eqn(t,x,mu1,mu2,x1t,omega,3), 0:dt:tau, x, opt);  % (strobo section) Note: Integrating 
   yy = PCC4BP_eqn(0,x,mu1,mu2,x1t,omega,3);
   y = yy(1:n) + lmda*S*dvdx';
   Y_coll = [Y_coll; y];
end


end