function y = transform_frame(mani, a_rel, n_3, mu, mu_bar, tag)
%%%%%% This funcion assumes that perturber initial angle = 0, if not 0 then
%%%%%% include rotation matrix as stated by kumar.et.al
Rot = @(th) [cos(th) -sin(th); sin(th) cos(th)];
y = zeros(length(mani),4);

for i=1:length(mani)
    pos = Rot(mani(i,5))*[mani(i,1)+mu; mani(i,2)];
    vel = Rot(mani(i,5))*mani(i,3:4)';
    y(i,1) = pos(1)/a_rel - mu_bar;
    y(i,2) = pos(2)/a_rel;
    y(i,3) = (vel(1)+(n_3-1)*pos(2))/(a_rel*n_3);
    y(i,4) = (vel(2)-(n_3-1)*pos(1))/(a_rel*n_3);
end

if tag == 1
    plot(y(:,1),y(:,2),'-k',0,0,'*k');hold on;R = 1;syms T;fplot(-mu_bar+R*sin(T),R*cos(T),'k','LineWidth',2);
end