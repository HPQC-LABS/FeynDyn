%% Bearing in mind the effort that went into the development of FeynDyn,
%% the author would indeed appreciate if users can please cite:
%% 1) The paper: N. Dattani (2013) Comp. Phys. Comm. Volume 184, Issue 12, Pg. 2828-2833 , and
%% 2) The code:  N. Dattani (2013) FeynDyn. http://dx.doi.org/10.6084/m9.figshare.823549

%% To make sure your code is the most updated version, please e-mail dattani.nike@gmail.com
%% Bug reports, suggestions, and requests for extensions are more than encouraged: dattani.nike@gmail.com

%% This program is protected by copyright law. Nike S. Dattani , 13/8/12

%% This sample input file, upon pressing F5, plots the real and imaginary parts of <i|rho|j> as a function 
% of time for the hamiltonian presented in Comp. Phys. Comm. Volume 184, Issue 12, Pg. 2828-2833

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. fundamental constants
kb=1.3806504*10^(-23);                 % Joules / Kelvin
hbar=1.054571628*10^(-34);             % Joule * seconds

%% 2. temperature
temperature=1/100/kb;
beta=1/(kb*temperature);

%% 3. system hamiltonian
wc=1/hbar;
Omega=0.1*wc;%hbar*1e12*pi/8;          % Rabi frequency in Joules   %bho
H=0.5*hbar*[0 Omega; Omega 0];         % System hamiltonian in Joules
%H=hbar*[Omega Omega; Omega -Omega];

%% 4. initial denisty matrix
rhoInitial=[1 0; 0 0];

%% 5. system coupling matrix\
%systemCouplingMatrix=[0 0;0 1];       % The system is coupled to the bath via |1X1|
systemCouplingMatrix=[1 0;0 -1];       % The system is coupled to the bath via sigma_z

%% 6. spectral distribution function
Nbaths=1;                              % Number of baths
alpha=0.9*0.022*pi*hbar/2;             % prefactor. In ps^2 / rad^2
%wc=10*Omega;                          % cutoff frequency. In rad/ps
dw=0.01*wc;                            % stepsize for w in J(w). In rad/ps
w=dw:dw:10*wc;                         % must start with dw because coth(0)=infinity. In rad/ps
%J=alpha*exp(-(w/wc)).*w;              % spectral density. In rad/psJ
J=alpha*wc^0.75.*w.^0.25;

%J=J*1e12*hbar;                        % spectral distribution function. In Joules

%w=w*1e12;                             % w in s-1
%dw=dw*1e12;                           % dw in s-1

%% 7. numerical parameters for Feynman integral
deltaKmax=12;                          % number of time steps before memory kernel dies
totalT=200/wc;%10/1e12;                % total time for the simulation, in seconds
%totalT=15/Omega;
dt=totalT/60;                          % size of time step (delta t)
finalPoint=round(totalT/dt);           % number of timesteps in total
allPointsORjustFinalPoint='allPoints'; % do you just want the density matrix at time=totalT ? or do you want it at every point until then
cpuORgpu='cpu'; 

%% 8. build input array that represents the densitry matrix (flattened to a column vector) as a function of time (where the columns represent the time steps)
rho=zeros(numel(H),finalPoint+1);
rho(:,1)=reshape(rhoInitial.',[],1);   % I'm not sure about the transpose, since I'm not sure about labeling of the states in the column vector for the purposes of the final sum function.

%% 9. run the program and plot only <0|rho|0> and <1|rho|1> as a function of time
% wholeDensityMatrixOrJustDiagonals='JustDiagonals';
% 
% [rho_onlyDiagonals,elapsedTime_onlyDiagonals]=FeynDynCode(finalPoint,deltaKmax,totalT,rho,H,systemCouplingMatrix,w,dw,J,temperature,wholeDensityMatrixOrJustDiagonals,allPointsORjustFinalPoint,cpuORgpu);
% 
% figure(1);hold('on')
% plot(0:totalT/finalPoint:totalT,real(rho_onlyDiagonals(1,:)));
% plot(0:totalT/finalPoint:totalT,real(rho_onlyDiagonals(4,:)));
% xlabel('time (seconds)');
%% 10. run the program and plot all elements of the density matrix as a function of time
wholeDensityMatrixOrJustDiagonals='wholeDensityMatrix';

[rho_allElements,elapsedTime_allElements]=FeynDynCode(Nbaths,finalPoint,deltaKmax,totalT,rho,H,systemCouplingMatrix,w,dw,J,temperature,wholeDensityMatrixOrJustDiagonals,allPointsORjustFinalPoint,cpuORgpu);

figure(3);hold('on')
totalT=totalT*wc;
plot(0:totalT/finalPoint*Omega:totalT*Omega,(real(rho_allElements(1,:))-real(rho_allElements(4,:))),'k');
% plot(0:totalT/finalPoint:totalT,real(rho_allElements(1,:)));
% plot(0:totalT/finalPoint:totalT,real(rho_allElements(3,:)),'r');plot(0:totalT/finalPoint:totalT,imag(rho_allElements(3,:)),'--r');
% plot(0:totalT/finalPoint:totalT,real(rho_allElements(3,:)),'r');plot(0:totalT/finalPoint:totalT,-imag(rho_allElements(3,:)),'--r');
 %plot(0:totalT/finalPoint:totalT,real(rho_allElements(4,:)));
xlabel('time (\omega_c^-^1)');

