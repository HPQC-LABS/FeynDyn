function sampleInput4x4

figure(1);
temperatures=[77 300];
for ii=1:length(temperatures)
     [rho, totalT, finalPoint, splineMesh]=sampleInput4x4engine(temperatures(ii));
    % save(strcat('T=',num2str(temperatures(ii)),'.mat'),'rho');
    
    totalT=10/1e9;
    finalPoint=10000;
    splineMesh=(totalT/finalPoint)/5;
    %load(strcat('T=',num2str(temperatures(ii)),'.mat'));
    
    subplot(length(temperatures),1,ii);hold('on');box('on');
    plotHandle(1)=plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,real(rho(14,:)),0:splineMesh:totalT),'b','LineWidth',3);
    % plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,real(rho(9,:)),0:splineMesh:totalT),'b','LineWidth',3);
    % plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,real(rho(10,:)),0:splineMesh:totalT),'b');
    % plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,imag(rho(10,:)),0:splineMesh:totalT),'LineWidth',3);
    % plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,real(rho(13,:)),0:splineMesh:totalT),'b','LineWidth',3);
    % plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,real(rho(14,:)),0:splineMesh:totalT),'b','LineWidth',3);
    % plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,real(rho(15,:)),0:splineMesh:totalT),'b','LineWidth',3);
    
    plotHandle(2)=plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,imag(rho(14,:)),0:splineMesh:totalT),'r','LineWidth',3);
    % plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,imag(rho(9,:)),0:splineMesh:totalT),'r','LineWidth',3);
    % plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,imag(rho(13,:)),0:splineMesh:totalT),'r','LineWidth',3);
    % plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,imag(rho(14,:)),0:splineMesh:totalT),'r','LineWidth',3);
    % plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,imag(rho(15,:)),0:splineMesh:totalT),'r','LineWidth',3);
    text(3.5, 0.0125, strcat('T=',num2str(temperatures(ii)),' K'), 'HorizontalAlignment', 'left', 'Interpreter','Latex','FontSize',42);
    
    ylabel('$\langle 0$X$|\rho(t)|$XX$\rangle$','Interpreter','Latex','FontSize',36)
    %set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',16,'Position',[0.0958333333333333 0.117718446601942+0.45*(ii-1) 0.895833333333333 0.412621359223301])
    set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',16)
    legendHandle=legend(plotHandle,'$\Re$','$\Im$','FontSize',42,'Interpreter','latex');set(legendHandle,'Interpreter','latex','FontSize',42,'YColor',[1 1 1],'XColor',[1 1 1],'Position',[0.910250208073242 0.142907691451381+0.45*(ii-1) 0.0716292134831461 0.115479115479115])
    axis([0 5 -0.02 0.02])
    if ~(ii==length(temperatures));set(gca,'XTickLabel','');
        % % %%%% INSET %%%% % %
        insetAxesHandle=axes('Position',[0.423809523809524 0.163347022587269 0.533928571428571 0.475]);hold('on')
        plotHandle(1)=plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,real(rho(14,:)),0:splineMesh:totalT),'b','LineWidth',3);
        plotHandle(2)=plot((0:splineMesh:totalT)*1e9,spline(0:totalT/finalPoint:totalT,imag(rho(14,:)),0:splineMesh:totalT),'r','LineWidth',3);
        axis([2 2.1 -0.0015 0.0015])
    end 
end
xlabel('Time (ns)','Interpreter','Latex','FontSize',42);
        
function [rho, totalT, finalPoint, splineMesh]=sampleInput4x4engine(temperature)
%% 1. fundamental constants
kb=1.3806504*10^(-23);        % Joules / Kelvin
hbar=1.054571628*10^(-34);    % Joule * seconds

%% 2. temperature
%temperature=77;
beta=1/(kb*temperature);

%% 3. system hamiltonian
epsilon=hbar*1e12*0.561;
Omega=hbar*1e12*0.015;
J=hbar*1e12*0.111;                     % Rabi frequency in Joules
H=[epsilon Omega Omega 0;
    Omega   0     J     Omega;
    Omega   J     0     Omega;
    0       Omega Omega -epsilon];      % System hamiltonian in Joules

%% 4. initial denisty matrix
rhoInitial=[0 0   0   0;
    0 0.5 0.5 0;
    0 0.5 0.5 0;
    0 0   0   0];

%% 5. system coupling matrix
systemCouplingMatrix=[0 0   0   0;
    0 1   0   0;
    0 0   1   0;
    0 0   0   2];   % The system is coupled to the bath via |1X1|

%% 6. spectral distribution function
Nbath=1;
alpha=0.027*pi;               % prefactor. In ps^2 / rad^2
wc=2.2;                       % cutoff frequency. In rad/ps
dw=0.01;                      % stepsize for w in J(w). In rad/ps
w=dw:dw:14;                   % must start with dw because coth(0)=infinity. In rad/ps
J=alpha*exp(-(w/wc).^2).*w.^3;% spectral density. In rad/ps

J=J*1e12*hbar;                % spectral distribution function. In Joules
%J=J*0;

w=w*1e12;                     % w in s-1
dw=dw*1e12;                   % dw in s-1

%% 7. numerical parameters for Feynman integral
deltaKmax=2;                  % number of time steps before memory kernel dies
totalT=10/1e9;                % total time for the simulation, in seconds
dt=10000/totalT;              % size of time step (delta t)
finalPoint=round(dt*totalT);  % number of timesteps in total
tol=1e-1;
allPointsORjustFinalPoint='allPoints'; % do you just want the density matrix at time=totalT ? or do you want it at every point until then
cpuORgpu='cpu';

%% 8. build input array that represents the densitry matrix (flattened to a column vector) as a function of time (where the columns represent the time steps)
rho=zeros(numel(H),finalPoint+1);
rho(:,1)=reshape(rhoInitial.',[],1); % I'm not sure about the transpose, since I'm not sure about labeling of the states in the column vector for the purposes of the final sum function.

%% 10. run the program and plot all elements of the density matrix as a function of time
wholeDensityMatrixOrJustDiagonals='wholeDensityMatrix';
[rho,elapsedTime]=FeynDyn(Nbath,finalPoint,deltaKmax,totalT,rho,H,systemCouplingMatrix,w,dw,J,temperature,wholeDensityMatrixOrJustDiagonals,allPointsORjustFinalPoint,cpuORgpu);  
splineMesh=(totalT/finalPoint)/5;
%% 11. Make plots

%
% figure(1);hold('on')
% plot(0:totalT/finalPoint:totalT,real(rho(1,:)),'k');
% plot(0:totalT/finalPoint:totalT,real(rho(6,:)),'b','LineWidth',3);
% plot(0:totalT/finalPoint:totalT,real(rho(11,:)),'r');
% plot(0:totalT/finalPoint:totalT,real(rho(16,:)),'g');
% xlabel('time (seconds)');