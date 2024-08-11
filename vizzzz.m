function y = vizzzz(T_unst, M_unst, mu1, mu2, a_rel, a, n_3)

%%%%%%% VIDEO ANIMATION %%%%%%%%%%%
vidObj = VideoWriter('callisto_globaltransfer','MPEG-4');
open(vidObj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Callisto/Europa orbit
figure(1); R = 1;syms T;fplot(-mu1+R*sin(T),R*cos(T),'k','LineWidth',2);grid on;
hold on; fplot(-mu1+(a_rel+mu1)*sin(T),(a_rel+mu1)*cos(T),'-m','LineWidth',2); hold on; fplot(-mu1+(a+mu1)*sin(T),(a+mu1)*cos(T),'-k','LineWidth',2); hold on;

plot(0,0,'ok');hold on;h = plot(0,1,'.');hold on;o = plot(0,1,'.');

pts = gobjects(1393,3);
vanish = 30;

for i=1:length(M_unst)
    zz = plot_inertial(T_unst(i),M_unst(i,:));
    % eeee = calc_jacobi(M_unst(i,:)',mu1, mu2, a_rel);  % Calculating hamiltonian in 4 body system
    % if rem(i,200) == 0
    %     figure(2);plot(i, eeee, '.r');hold on;grid on;
    % end
    delete(h);delete(o);
    figure(1);
    pts(i,1) = plot(zz(1),zz(2),'.r','MarkerSize',15);pause(0.000001);hold on;
    h = plot(-mu1+cos(T_unst(i)),sin(T_unst(i)),'.k','MarkerSize',20);hold on;
    o = plot(-mu1+a_rel*cos(n_3*T_unst(i)),a_rel*sin(n_3*T_unst(i)),'.m','MarkerSize',20);hold on;
    xlabel('X (DU)');ylabel('Y (DU)');
    % h = plot(-mu1+cos(T_unst(i)),sin(T_unst(i)),'.k','MarkerSize',20);hold on;
    % o = plot(-mu1+a_rel*cos(n_3*T_unst(i)),a_rel*sin(n_3*T_unst(i)),'.m','MarkerSize',20);pause(0.000001);
    if i>vanish
        delete(pts(i-vanish+1,:))
    end
    
    fme = getframe(gcf);
    writeVideo(vidObj,fme);

end
close(vidObj)
y = 0;
end