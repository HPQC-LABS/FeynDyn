%% Press F5. A figure from Nancy Makri's 1995 JCP paper will be plotted.

%% 1. fundamental constants
kb=1.3806504*10^(-23);           % Joules / Kelvin
hbar=1.054571628*10^(-34);       % Joule * seconds

%% 2. temperature
temperature=1/5/kb;
beta=1/(kb*temperature);

%% 3. system hamiltonian
Omega=1/hbar;%hbar*1e12*pi/8;    % Rabi frequency in Joules   %bho
H=hbar*[0 Omega; Omega 0];       % System hamiltonian in Joules
%H=hbar*[Omega Omega; Omega -Omega];

%% 4. initial denisty matrix
rhoInitial=[1 0; 0 0];

%% 5. system coupling matrix
Nbaths=1;
%systemCouplingMatrix=[0 0;0 1]; % The system is coupled to the bath via |1X1|
systemCouplingMatrix=[1 0;0 -1]; % The system is coupled to the bath via sigma_z

%% 6. spectral distribution function
alpha=0.05*pi*hbar;              % prefactor. In ps^2 / rad^2
wc=7.5*Omega;                    % cutoff frequency. In rad/ps
dw=0.01*Omega;                   % stepsize for w in J(w). In rad/ps
w=dw:dw:10*wc;                   % must start with dw because coth(0)=infinity. In rad/ps
J=alpha*exp(-(w/wc)).*w;% spectral density. In rad/ps

%J=J*1e12*hbar;                  % spectral distribution function. In Joules

%w=w*1e12;                       % w in s-1
%dw=dw*1e12;                     % dw in s-1

%% 7. numerical parameters for Feynman integral
deltaKmax=1;                     % number of time steps before memory kernel dies
totalT=25/Omega;%10/1e12;        % total time for the simulation, in seconds
%totalT=15/Omega;
dt=totalT/100;                   % size of time step (delta t)
finalPoint=round(totalT/dt);     % number of timesteps in total
allPointsORjustFinalPoint='allPoints'; % do you just want the density matrix at time=totalT ? or do you want it at every point until then
cpuORgpu='cpu'; 

%% 8. build input array that represents the densitry matrix (flattened to a column vector) as a function of time (where the columns represent the time steps)
rho=zeros(numel(H),finalPoint+1);
rho(:,1)=reshape(rhoInitial.',[],1); % I'm not sure about the transpose, since I'm not sure about labeling of the states in the column vector for the purposes of the final sum function.

%% 9. run the program and plot only <0|rho|0> and <1|rho|1> as a function of time
% wholeDensityMatrixOrJustDiagonals='JustDiagonals';
% 
% [rho_onlyDiagonals,elapsedTime_onlyDiagonals]=FeynDyn(Nbaths,finalPoint,deltaKmax,totalT,rho,H,systemCouplingMatrix,w,dw,J,temperature,wholeDensityMatrixOrJustDiagonals,allPointsORjustFinalPoint,cpuORgpu);
% 
% figure(1);hold('on')
% plot(0:totalT/finalPoint:totalT,real(rho_onlyDiagonals(1,:)));
% plot(0:totalT/finalPoint:totalT,real(rho_onlyDiagonals(4,:)));
% xlabel('time (seconds)');
%% 10. run the program and plot all elements of the density matrix as a function of time
wholeDensityMatrixOrJustDiagonals='wholeDensityMatrix';

[rho_allElements,elapsedTime_allElements]=FeynDyn(Nbaths,finalPoint,deltaKmax,totalT,rho,H,systemCouplingMatrix,w,dw,J,temperature,wholeDensityMatrixOrJustDiagonals,allPointsORjustFinalPoint,cpuORgpu);
totalT=totalT*Omega;
figure(2);hold('on')
plot(0:totalT/finalPoint:totalT,(real(rho_allElements(1,:))-real(rho_allElements(4,:))),'k');
% plot(0:totalT/finalPoint:totalT,real(rho_allElements(1,:)));
% plot(0:totalT/finalPoint:totalT,real(rho_allElements(3,:)),'r');plot(0:totalT/finalPoint:totalT,imag(rho_allElements(3,:)),'--r');
% plot(0:totalT/finalPoint:totalT,real(rho_allElements(3,:)),'r');plot(0:totalT/finalPoint:totalT,-imag(rho_allElements(3,:)),'--r');
% plot(0:totalT/finalPoint:totalT,real(rho_allElements(4,:)));
xlabel('time (seconds)');

