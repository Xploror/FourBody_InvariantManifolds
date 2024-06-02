function [PHI,xf] = get_STM(T, xi, d)

mu = 2.5263e-05; mu1 = 7.7893e-05; p = 1.5949; n = 6.2864; % for L1 (J-E-G)
N = length(xi)/d;

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
    for i=time
        A = PCC4BP_J(x(1), x(2), x(5), mu, mu1, p);

        %%% 6th Order RK method for integration
%         k1 = tStep*phidot(A, phi);
%         k2 = tStep*phidot(A, phi+k1);
%         k3 = tStep*phidot(A, phi+(3*k1+k2)/8);
%         k4 = tStep*phidot(A, phi+(8*k1+2*k2+8*k3)/27);
%         k5 = tStep*phidot(A, phi+((-21+9*21^(1/2))*k1-8*(7-21^(1/2))*k2+48*(7-21^(1/2))*k3-3*(21-21^(1/2))*k4)/392);
%         k6 = tStep*phidot(A, phi+(-5*(231+51*21^0.5)*k1-40*(7+21^0.5)*k2-320*21^0.5*k3+3*(21+121*21^0.5)*k4+392*(6+21^0.5)*k5)/1960);
%         k7 = tStep*phidot(A, phi+(15*(22+7*21^0.5)*k1+120*k2+40*(-5+7*21^0.5)*k3-63*(-2+3*21^0.5)*k4-14*(49+9*21^0.5)*k5+70*(7-21^0.5)*k6)/180);
%         phi = phi + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;
        for j=1:length(x)
            pphhii = phi(:,j);
            k1 = tStep*phidot(A, pphhii);
            k2 = tStep*phidot(A, pphhii+k1);
            k3 = tStep*phidot(A, pphhii+(3*k1+k2)/8);
            k4 = tStep*phidot(A, pphhii+(8*k1+2*k2+8*k3)/27);
            k5 = tStep*phidot(A, pphhii+((-21+9*21^(1/2))*k1-8*(7-21^(1/2))*k2+48*(7-21^(1/2))*k3-3*(21-21^(1/2))*k4)/392);
            k6 = tStep*phidot(A, pphhii+(-5*(231+51*21^0.5)*k1-40*(7+21^0.5)*k2-320*21^0.5*k3+3*(21+121*21^0.5)*k4+392*(6+21^0.5)*k5)/1960);
            k7 = tStep*phidot(A, pphhii+(15*(22+7*21^0.5)*k1+120*k2+40*(-5+7*21^0.5)*k3-63*(-2+3*21^0.5)*k4-14*(49+9*21^0.5)*k5+70*(7-21^0.5)*k6)/180);
            pphhii = pphhii + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;
            phi(:,j) = pphhii;
        end


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
    PHI(1+d*(j-1):d+d*(j-1), 1+d*(j-1):d+d*(j-1)) = phi(1:d,1:d);
end

end