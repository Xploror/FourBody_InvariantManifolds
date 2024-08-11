function y = calc_jacobi(x, mu, mu2, a_rel)

    r1 = x(1:2)'-[-mu 0];
    r2 = x(1:2)'-[1-mu 0];
    Hamilt = @(x,mu) 0.5*x(3)^2 + 0.5*x(4)^2 - 0.5*x(1)^2 - 0.5*x(2)^2 - (1-mu)/norm(r1) - (mu)/norm(r2);

    if length(x)>4
        r3 = x(1:2)'-[-mu+a_rel*cos(x(5)), a_rel*sin(x(5))];
        r13 = [a_rel*cos(x(5)), a_rel*sin(x(5))];
        Hamilt4bp = @(x,mu) 0.5*x(3)^2 + 0.5*x(4)^2 - 0.5*(1-mu)*norm(r1)^2 - 0.5*mu*norm(r2)^2 - (1-mu)/norm(r1) - (mu)/norm(r2) - (mu2)/norm(r3) + (mu2)*x(1)*cos(x(5))/(norm(r13)^2) + (mu2)*x(2)*sin(x(5))/(norm(r13)^2);
        y = Hamilt4bp(x, mu);  % 4 Body capability 
    else
        y = Hamilt(x, mu);  % 3 Body capability
    end

end