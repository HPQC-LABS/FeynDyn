%% Preliminary version of VPQME.

totalT=4/1e11                                                                   %Time of evolution in ps     
rhoi=[0 0; 0 1];                                                         %Initial density matrix

options=odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-2]);
[T,RHO]=ode113(@phononm,[0 totalT],[rhoi(1,1) real(rhoi(1,2)) imag(rhoi(1,2))],options);    %Time and initial conditions  


%Plot of the quantum evolution of the density matrix

figure(1);hold('on')
plotHandle(1)=plot(T,RHO(:,1),'k','LineWidth',5);
plotHandle(2)=plot(T,RHO(:,2),'r','LineWidth',5);
plotHandle(3)=plot(T,RHO(:,3),'b','LineWidth',5);
plotHandle(4)=plot(T,1-RHO(:,1),'color',[0.5 0.5 0.5],'LineWidth',5);

legendHandle=legend(plotHandle,'$\langle1|\rho(t)|1\rangle$','$\Re\langle1|\rho(t)|0\rangle$','$\Im\langle1|\rho(t)|0\rangle$','$\langle0|\rho(t)|0\rangle$');
set(legendHandle,'Interpreter','latex','FontSize',32,'LineWidth',2)

axis([0 totalT min(min(RHO)) 1]);
box('on');grid('on')
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',16);
title('Photon Induced Dephasing at T=50K','interpreter','latex','FontSize',56)
ylabel('Elements of $\rho(t)$','Interpreter','latex','FontSize',40)
xlabel('Time (seconds)','Interpreter','latex','FontSize',40);
yticklabels('auto');

set(gcf, 'Color', 'w');
set(gcf,'renderer','Painters')

%Plot the quantum evolution of the element pxx or p00.


figure(2);hold('on')
plotHandle(1)=plot(T,RHO(:,1),'k','LineWidth',5);

legendHandle=legend(plotHandle,'$\langle X|\rho(t)|X \rangle$');
set(legendHandle,'Interpreter','latex','FontSize',32,'LineWidth',2)

axis([0 totalT min(min(RHO)) 1]);
box('on');grid('on')
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',16);
title('Photon Induced Dephasing at T=50K','interpreter','latex','FontSize',56)
ylabel('Elements of $\rho(t)$','Interpreter','latex','FontSize',40)
xlabel('Time (seconds)','Interpreter','latex','FontSize',40);
yticklabels('auto');

set(gcf, 'Color', 'w');
set(gcf,'renderer','Painters')
