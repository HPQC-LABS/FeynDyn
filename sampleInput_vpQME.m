%%%% Press F5. Figure 1 will appear and show you the dynamics of a quantum dot qubit
%following the model from A.J. Ramsay et al. (2009) "Supplement: Theoretical model of photon induced dephasing" 
%Equations S23 and S24. 

Time=2    %Time of evolution in ps
A=[0.5 0.5; 0.5 0.5]; %Initial density matrix

options=odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4 1e-4]);
[T,Y]=ode113(@phononm,[0 Time],[A(1,1) real(A(1,2)) imag(A(1,2)) 0],options);    %Time and initial conditions  


%Plot of the quantum evolution of the density matrix

figure(1);hold('on')
plotHandle(1)=plot(T,Y(:,1),'k','LineWidth',5);
plotHandle(2)=plot(T,Y(:,2),'r','LineWidth',5);
plotHandle(3)=plot(T,Y(:,3),'b','LineWidth',5);
plotHandle(4)=plot(T,1-Y(:,1),'color',[0.5 0.5 0.5],'LineWidth',5);

legendHandle=legend(plotHandle,'$\langle1|\rho(t)|1\rangle$','$\Re\langle1|\rho(t)|0\rangle$','$\Im\langle1|\rho(t)|0\rangle$','$\langle0|\rho(t)|0\rangle$');
set(legendHandle,'Interpreter','latex','FontSize',32,'LineWidth',2)

axis([0 Time min(min(Y)) max(max((Y)))]);
box('on');grid('on')
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',16);
title('Photon Induced Dephasing','interpreter','latex','FontSize',56)
ylabel('Elements of $\rho(t)$','Interpreter','latex','FontSize',40)
xlabel('Time (pico-seconds)','Interpreter','latex','FontSize',40);
yticklabels('auto');

set(gcf, 'Color', 'w');
set(gcf,'renderer','Painters')

