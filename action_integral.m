function y = action_integral(x)
% x --> collected vector
pos = []; vel = [];
for i=1:length(x)/4
    pos = [pos x(1+4*(i-1):2+4*(i-1))];
    vel = [vel x(3+4*(i-1):4+4*(i-1))];
end
momenta_x = vel(1,:) - pos(2,:);
momenta_y = vel(2,:) + pos(1,:);
M = [momenta_x; momenta_y];

pos_next = [pos(1,1) pos(1,:); pos(2,1) pos(2,:)];
del_pos = pos - pos_next(:,1:end-1);
S = sum(dot(M,del_pos));

y = S;  % S is the torus action


end