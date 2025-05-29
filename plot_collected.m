function p = plot_collected(x, tag, col)
% This script only works for states of size 4 i.e. 2 dimensional
N = length(x)/4;

p = 0;
if tag == "position"
    for i=0:N-1
        plot(x(1+4*i),x(2+4*i),'.','MarkerSize',10,'color',col);
        hold on;
        title('Trajectory')
    end
elseif tag == "velocity"
    for i=0:N-1
        quiver(x(1+4*i),x(2+4*i),x(3+4*i),x(4+4*i));
        hold on;
        title('Vector field')
    end
elseif tag == "position-connect1" || tag == "position-connect2" || tag == "position-connect1 -" || tag == "position-connect2 -"
    X = []; Y = [];
    if tag == "position-connect1" || tag == "position-connect1 -"
        %index = [1 4 7 10 13 16 19 22 25 28 31 34 2 5 8 11 14 17 20 23 26 29 32 35 3 6 9 12 15 18 21 24 27 30 33 1]; %(L1)
        index = [1 25 8 32 15 39 22 5 29 12 36 19 2 26 9 33 16 40 23 6 30 13 37 20 3 27 10 34 17 41 24 7 31 14 38 21 4 28 11 35 18 1];
    elseif tag == "position-connect2" || tag == "position-connect2 -"
        index = [1 5 9 13 17 21 25 29 33 2 6 10 14 18 22 26 30 34 3 7 11 15 19 23 27 31 35 4 8 12 16 20 24 28 32 1];  %(L2)
        %index = [1 12 23 34 4 15 26 37 7 18 29 40 10 21 32 2 13 24 35 5 16 27 38 8 19 30 41 11 22 33 3 14 25 36 6 17 28 39 9 20 31 1];
    end
    X = [X x(1+4*(index-1))]; 
    Y = [Y x(2+4*(index-1))];
    if tag == "position-connect1" || tag == "position-connect2"
        plot(X,Y,'-','LineWidth',1.5,'color',col);
    elseif tag == "position-connect1 -" || tag == "position-connect2 -"
        plot(X,Y,'--','LineWidth',1,'color',col);
    end
    hold on;
    title('Trajectory')
end
end