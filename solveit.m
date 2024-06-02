function y = solveit(X,mu,ref,tag)
% Functionality only for 2 dimension right now 

if length(X) == 4 || length(X) == 5
    if tag == "energy"
        Hamilt = @(x) 0.5*X(3)^2 + 0.5*X(4)^2 - 0.5*x^2 - 0.5*X(2)^2 - (1-mu)/sqrt((x+mu)^2 + X(2)^2) - mu/sqrt((x+mu-1)^2 + X(2)^2) - 0.5*(mu*(1-mu)) - ref;
        x_new = fzero(Hamilt, X(1));
        Hamilt(x_new)
        Hamilt(X(1))
    end
elseif length(X) == 6
    if tag == "energy"
        Hamilt = @(x) 0.5*X(4)^2 + 0.5*X(5)^2 + 0.5*X(6)^2 - 0.5*x^2 - 0.5*X(2)^2 - (1-mu)/sqrt((x+mu)^2 + X(2)^2 + X(3)^2) - mu/sqrt((x+mu-1)^2 + X(2)^2 + X(3)^2) - 0.5*(mu*(1-mu)) - ref;
        x_new = fzero(Hamilt, X(1));
    end
end

y = x_new;
end