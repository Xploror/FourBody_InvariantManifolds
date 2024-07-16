function theta = parameterize(x, param, tag)

if tag == "2D"
    % FORMAT: param = [lp]
    lp = param(1);  % libration point location
    theta = atan2(x(:,2), (x(:,1)-lp));
    for i=1:length(theta)
        if theta(i)<0
            theta(i) = 2*pi + theta(i);
        end
    end
end

end