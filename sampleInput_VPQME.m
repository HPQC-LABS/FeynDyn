% Preliminary version of VPQME.

%% 1) Fundamental constants

kb=1.38064852e-23;              %Boltzmann constant K/K
hbar=1.0545718e-34;             %Planck constant Js

%% 2) Paramaters of the system

Theta=30;                     
tau=14e-12;                     %Full width half maximum in s
T=25;                           %Temperature in K
wc=2.2e12;                      %1/s   
alpha=pi*0.027*1e-24;           % prefactor. In s^2 / rad^2
totalT=4/1e11;                  %Time of evolution in s   

%% 3)Coefficients of the differential equation 

Ot = linspace(0,5/1e11,100);
Omega=(Theta/(2*tau*sqrt(pi)))*exp(-(Ot/(2*tau)).^2);

Kt = linspace(0,4/1e11,100);
K=3e11+Kt*0;            

%% 4) Initial density matrix

rho=[0.5 0.5; 0.5 0.5]; %Initial density matrix

%% 5) ODE configuration
                                    
options=odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4 1e-4]);
[t,Y]=ode45(@(t,rho) VPQME(t,rho,Kt,K,Ot,Omega,kb,hbar,Theta,tau,T,wc,alpha),[0 totalT],[rho(1,1) real(rho(1,2)) imag(rho(1,2)) rho(2,2)],options); 
%[t,Y]=ode113(@VPQME,[0 totalT],[rho(1,1) real(rho(1,2)) imag(rho(1,2))],options);

%% 6) Plot of the RDM dynamics

figure(1);hold('on');
plotHandle(1)=plot(t,Y(:,1),'k','LineWidth',5);
plotHandle(2)=plot(t,Y(:,2),'r','LineWidth',5);
plotHandle(3)=plot(t,Y(:,3),'b','LineWidth',5);
plotHandle(4)=plot(t,Y(:,4),'color',[0.5 0.5 0.5],'LineWidth',5);

legendHandle=legend(plotHandle,'$\langle1|\rho(t)|1\rangle$','$\Re\langle1|\rho(t)|0\rangle$','$\Im\langle1|\rho(t)|0\rangle$','$\langle0|\rho(t)|0\rangle$');
set(legendHandle,'Interpreter','latex','FontSize',32,'LineWidth',2);

axis([0 totalT min(min(Y)) max(max((Y)))]);
box('on');grid('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',16);
title('Photon Induced Dephasing','interpreter','latex','FontSize',56);
ylabel('Elements of $\rho(t)$','Interpreter','latex','FontSize',40);
xlabel('Time (seconds)','Interpreter','latex','FontSize',40);
yticklabels('auto');

set(gcf, 'Color', 'w');
set(gcf,'renderer','Painters');

%% 7) Plot of pxx.

figure(2);hold('on');
plotHandle(1)=plot(t,Y(:,1),'k','LineWidth',5);

legendHandle=legend(plotHandle,'$\langle X|\rho(t)|X \rangle$');
set(legendHandle,'Interpreter','latex','FontSize',32,'LineWidth',2);

axis([0 totalT min(min(Y)) max(max(Y))]);
box('on');grid('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',16);
title('Photon Induced Dephasing at T=25K','interpreter','latex','FontSize',56);
ylabel('Elements of $\rho(t)$','Interpreter','latex','FontSize',40);
xlabel('Time (seconds)','Interpreter','latex','FontSize',40);
yticklabels('auto');

set(gcf, 'Color', 'w');
set(gcf,'renderer','Painters');
set(gcf, 'Color', 'w');
set(gcf,'renderer','Painters');
