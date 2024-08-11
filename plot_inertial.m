function y = plot_inertial(t, x)
%%%%%% x --> rotating frame coordinates | t --> rescaled time coordinates
Rot = @(t) [cos(t) -sin(t); sin(t) cos(t)];

X = [];
for i=1:length(t)
    X = [X Rot(t(i))*x(i,1:2)'];
end
y = X';
end