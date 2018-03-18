%% Press F5. Figures 5 and 6 will be plotted from Eunji Sim (2001) "Quantum dynamics for a system coupled to slow baths:
%  On-the-fly filtered propagator method", J. Chem. Phys. 115, 10, pg. 4450 

%% fundamental constants
kb=1.3806504*10^(-23);            % Joules / Kelvin
hbar=1.054571628*10^(-34);        % Joule * seconds

%% temperature
beta=0.72/hbar;
temperature=1/(kb*beta);

%% system hamiltonian
V12=0.066667;V23=V12;
E2=-1.3333333;E3=2*E2;

H=[0   V12   0; 
   V12 E2  V23;
   0   V23  E3]*hbar;

%% initial denisty matrix
rhoInitial=[1 0 0; 
            0 0 0; 
            0 0 0];

%% system coupling matrix
systemCouplingMatrix=diag([-1 0 1]);   

%% spectral distribution function
alpha=2.094*hbar;              % multiplied by hbar because FeynDyn has an erroneous hbar in MakMak
wc=1.0;                        % cutoff frequency 
dw=0.0001;
w=dw:dw:14;                    % must start with dw because coth(0)=infinity
J=alpha*exp(-(w/wc).^1).*w.^1; % spectral density in ps-1

%% numerical parameters for Feynman integral
finalPoint=3000;               % number of timesteps in total
deltaKmax=6;                   % number of time steps before memory kernel dies
totalT=1000;                   % total time for the simulation, in seconds
allPointsORjustFinalPoint='allPoints'; % do you just want the density matrix at time=totalT ? or do you want it at every point until then
cpuORgpu='cpu'; 

%% build input array that represents the densitry matrix (flattened to a column vector) as a function of time (where the columns represent the time steps)
rho=zeros(numel(H),finalPoint+1);
rho(:,1)=reshape(rhoInitial.',[],1); 

%% run the program to calculate all density matrix elements
wholeDensityMatrixOrJustDiagonals='wholeMatrix';

[rho,elapsedTime]=FeynDyn(finalPoint,deltaKmax,totalT,rho,H,systemCouplingMatrix,w,dw,J,temperature,wholeDensityMatrixOrJustDiagonals,allPointsORjustFinalPoint,cpuORgpu);

%% Plot figure 5 of Eunji Sim's paper: just diagonals
figure(5);hold('on')
plot(0:totalT/finalPoint:totalT,real(rho(1,:)));
plot(0:totalT/finalPoint:totalT,real(rho(5,:)));
plot(0:totalT/finalPoint:totalT,real(rho(9,:)));
xlabel('time');

%% Plot figure 6 of Eunji Sim's paper: just rho(2,3)
figure(6);hold('on')
plot(0:totalT/finalPoint:totalT,real(rho(8,:)));
xlabel('time');