%% Preliminary version of VPQME.

totalT=2;               %Time of evolution in ps
rho=[0.5 0.5; 0.5 0.5]; %Initial density matrix

options=odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4 1e-4]);
[t,Y]=ode45(@VPQME,[0 totalT],[rho(1,1) real(rho(1,2)) imag(rho(1,2)) 0],options);    
%[t,Y]=ode113(@VPQME,[0 totalT],[rho(1,1) real(rho(1,2)) imag(rho(1,2)) 0],options);    

%% Plot of the RDM dynamics

figure(1);hold('on');
plotHandle(1)=plot(t,Y(:,1),'k','LineWidth',5);
plotHandle(2)=plot(t,Y(:,2),'r','LineWidth',5);
plotHandle(3)=plot(t,Y(:,3),'b','LineWidth',5);
plotHandle(4)=plot(t,1-Y(:,1),'color',[0.5 0.5 0.5],'LineWidth',5);

legendHandle=legend(plotHandle,'$\langle1|\rho(t)|1\rangle$','$\Re\langle1|\rho(t)|0\rangle$','$\Im\langle1|\rho(t)|0\rangle$','$\langle0|\rho(t)|0\rangle$');
set(legendHandle,'Interpreter','latex','FontSize',32,'LineWidth',2);

axis([0 totalT min(min(Y)) max(max((Y)))]);
box('on');grid('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',16);
title('Photon Induced Dephasing','interpreter','latex','FontSize',56);
ylabel('Elements of $\rho(t)$','Interpreter','latex','FontSize',40);
xlabel('Time (pico-seconds)','Interpreter','latex','FontSize',40);
yticklabels('auto');

set(gcf, 'Color', 'w');
set(gcf,'renderer','Painters');