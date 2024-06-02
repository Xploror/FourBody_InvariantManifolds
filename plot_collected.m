function p = plot_collected(x, tag, col)
% This script only works for states of size 4 i.e. 2 dimensional
N = length(x)/4;

p = 0;
figure(1);
if tag == "position"
    for i=0:N-1
        plot(x(1+4*i),x(2+4*i),'.','color',col);
        hold on;
        title('Trajectory')
    end
elseif tag == "velocity"
    for i=0:N-1
        quiver(x(1+4*i),x(2+4*i),x(3+4*i),x(4+4*i));
        hold on;
        title('Vector field')
    end
end