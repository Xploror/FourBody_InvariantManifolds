function y = check(x,mu,param,tag)

if length(x) == 6     % 3 dimensional
    r1 = sqrt((x(1)+mu)^2 + x(2)^2 + x(3)^2);
    r2 = sqrt((x(1)-(1-mu))^2 + x(2)^2 + x(3)^2);
    if tag == "energy"
        Ueff = -(1-mu)/r1 - mu/r2 -0.5*(x(1)^2 + x(2)^2) - 0.5*(mu*(1-mu));
        Hamilt = 0.5*(sum(x(4:6).^2)) + Ueff;
        y = Hamilt;
    end
elseif length(x) == 5  % four body problem
    % FORMAT : mu = [mu1 mu2]  | param = [r3 radius of perturber]
    mu1 = mu(1);mu2 = mu(2);
    r1 = sqrt((x(1)+mu1)^2 + x(2)^2);
    r2 = sqrt((x(1)-(1-mu1))^2 + x(2)^2);
    r3 = norm(x(1:2)-param(1)*[cos(x(5));sin(x(5))]);
    r13 = norm([-mu1,0]-param(1)*[cos(x(5));sin(x(5))]);
    if tag == "energy"
        Ueff = -(1-mu1)/r1 - mu1/r2 - mu2/r3 -0.5*(x(1)^2 + x(2)^2) + mu2*x(1)*cos(x(5))/(r13^2) + mu2*x(2)*cos(x(5))/(r13^2);
        Hamilt = 0.5*(sum(x(3:4).^2)) + Ueff;
        y = Hamilt;
    end
else
    r1 = sqrt((x(1)+mu)^2 + x(2)^2);
    r2 = sqrt((x(1)-(1-mu))^2 + x(2)^2);
    if tag == "energy"
        Ueff = -(1-mu)/r1 - mu/r2 -0.5*(x(1)^2 + x(2)^2) - 0.5*(mu*(1-mu));
        Hamilt = 0.5*(sum(x(3:4).^2)) + Ueff;
        y = Hamilt;
    end
end