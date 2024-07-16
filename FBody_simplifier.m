function [G, pos] = FBody_simplifier(x, N2, N, m)
% x --> collected vector | G --> states, pos --> positions
shifts = N*(m+1);

G = [];
for i=1:N
    G = [G x(1 + 4*N2*(m+1)*(i-1) : 4*N2 + 4*N2*(m+1)*(i-1))]; % Just grab all the invariant circle points at Nth timestep for first (perturber) angle
end
G = [G x(1 + 4*N2*(m+1)*N : 4*N2 + 4*N2*(m+1)*N)];

pos = [];
for i=1:N2
    pos = [pos; G(1+4*(i-1):2+4*(i-1), :)];
end

end