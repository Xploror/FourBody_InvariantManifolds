function [DPHI,PHI,xf] = get_STM(T, xi, d, flag)

%mu = 2.5263e-05; mu1 = 7.7893e-05; p = 1.5949; n = 6.2864; % for L1 (J-E-G)
mu = 5.6623e-05; mu1 = 7.7890e-05; p = 0.5688; n = 2.3309; % for L1 (J-E-C)
N = length(xi)/d;
if flag ==0
    mu1 = 0;
end

tStep = 1e-3;
time = 0:tStep:T;

f = @(x) PCC4BP_eqn(1, x, mu, mu1, p, n, 3);
phidot = @(A, phi) A*phi;
xf = [];
PHI = 0;
for j=1:N
    x = xi(1+d*(j-1):d+d*(j-1));
    if d==4
        x = [x;0];  % Just add the perturber angle so that PCCFBP can take it 
    end
    phi = eye(length(x));
    PHI_planar = eye(4);
    for i=time
        A = PCC4BP_J(x(1), x(2), x(5), mu, mu1, p);
        A_planar = A(1:4,1:4);

        %%% 6th Order RK method for integration
        % k1 = tStep*phidot(A, phi);
        % k2 = tStep*phidot(A, phi+k1);
        % k3 = tStep*phidot(A, phi+(3*k1+k2)/8);
        % k4 = tStep*phidot(A, phi+(8*k1+2*k2+8*k3)/27);
        % k5 = tStep*phidot(A, phi+((-21+9*21^(1/2))*k1-8*(7-21^(1/2))*k2+48*(7-21^(1/2))*k3-3*(21-21^(1/2))*k4)/392);
        % k6 = tStep*phidot(A, phi+(-5*(231+51*21^0.5)*k1-40*(7+21^0.5)*k2-320*21^0.5*k3+3*(21+121*21^0.5)*k4+392*(6+21^0.5)*k5)/1960);
        % k7 = tStep*phidot(A, phi+(15*(22+7*21^0.5)*k1+120*k2+40*(-5+7*21^0.5)*k3-63*(-2+3*21^0.5)*k4-14*(49+9*21^0.5)*k5+70*(7-21^0.5)*k6)/180);
        % phi = phi + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;

        k1 = tStep*phidot(A_planar, PHI_planar);
        k2 = tStep*phidot(A_planar, PHI_planar+k1);
        k3 = tStep*phidot(A_planar, PHI_planar+(3*k1+k2)/8);
        k4 = tStep*phidot(A_planar, PHI_planar+(8*k1+2*k2+8*k3)/27);
        k5 = tStep*phidot(A_planar, PHI_planar+((-21+9*21^(1/2))*k1-8*(7-21^(1/2))*k2+48*(7-21^(1/2))*k3-3*(21-21^(1/2))*k4)/392);
        k6 = tStep*phidot(A_planar, PHI_planar+(-5*(231+51*21^0.5)*k1-40*(7+21^0.5)*k2-320*21^0.5*k3+3*(21+121*21^0.5)*k4+392*(6+21^0.5)*k5)/1960);
        k7 = tStep*phidot(A_planar, PHI_planar+(15*(22+7*21^0.5)*k1+120*k2+40*(-5+7*21^0.5)*k3-63*(-2+3*21^0.5)*k4-14*(49+9*21^0.5)*k5+70*(7-21^0.5)*k6)/180);
        PHI_planar = PHI_planar + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;
        PHI_dot = phidot(A_planar, PHI_planar);


        k1 = tStep*f(x);
        k2 = tStep*f(x+k1);
        k3 = tStep*f(x+(3*k1+k2)/8);
        k4 = tStep*f(x+(8*k1+2*k2+8*k3)/27);
        k5 = tStep*f(x+((-21+9*21^(1/2))*k1-8*(7-21^(1/2))*k2+48*(7-21^(1/2))*k3-3*(21-21^(1/2))*k4)/392);
        k6 = tStep*f(x+(-5*(231+51*21^0.5)*k1-40*(7+21^0.5)*k2-320*21^0.5*k3+3*(21+121*21^0.5)*k4+392*(6+21^0.5)*k5)/1960);
        k7 = tStep*f(x+(15*(22+7*21^0.5)*k1+120*k2+40*(-5+7*21^0.5)*k3-63*(-2+3*21^0.5)*k4-14*(49+9*21^0.5)*k5+70*(7-21^0.5)*k6)/180);
        x = x + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;  
    end
    xf = [xf; x(1:d)];  % Gives final state after time t of evolution
    PHI(1+4*(j-1):4+4*(j-1), 1+4*(j-1):4+4*(j-1)) = PHI_planar;
    DPHI(1+4*(j-1):4+4*(j-1), 1+4*(j-1):4+4*(j-1)) = PHI_dot;
end

end