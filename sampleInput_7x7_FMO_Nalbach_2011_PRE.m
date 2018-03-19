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
kb=1.3806504*10^(-23);        % Joules / Kelvin
hbar=1.054571628*10^(-34);    % Joule * seconds

%% 2. temperature
temperature=77;
beta=1/(kb*temperature);

%% 3. system hamiltonian
%Omega=hbar*1e12*pi/8;         % Rabi frequency in Joules
wavenumbers2joules=1.9864475e-23;
Jda=100/5.034e22;
epsilon_d=200/5.034e22;
epsilon_a=300/5.034e22;
%H=[0 Omega/2; Omega/2 0];     % System hamiltonian in Joules
%H=Jda*[0 1;1 0]+[epsilon_d 0;0 epsilon_a];
H=wavenumbers2joules*[240 -87.7 5.5 -5.9 6.7 -13.7 -9.9
   -87.7 315 30.8 8.2 0.7 11.8 4.3
   5.5  30.8 0 -53.5 -2.2 -9.6 6.0
   -5.9 8.2 -53.5 130 -70.7 -17.0 -63.3
   6.7 0.7 -2.2 -70.7 285 81.1 -1.3
   -13.7 11.8 -9.6 -17.0 81.1 435 39.7
   -9.9 4.3 6.0 -63.3 -1.3 39.7 245];

%% 4. initial denisty matrix
Nbaths=7;
rhoInitial=zeros(Nbath);
rhoInitial(1)=1;

%% 5. system coupling matrix
systemCouplingMatrix=diag(ones(7,1));
%systemCouplingMatrix=[0 0;0 1];   % The system is coupled to the bath via |1X1|
%systemCouplingMatrix=[1 0;0 -1]; % The system is coupled to the bath via sigma_z

%% 6. spectral distribution function

wc=1/(50e-15);%53/5.034e22;                       % cutoff frequency. In rad/ps
dw=0.00001*wc;                      % stepsize for w in J(w). In rad/ps
w=dw:dw:50*wc;                   % must start with dw because coth(0)=infinity. In rad/ps
lamda=35*wavenumbers2joules;%2/pi*Jda/5;
%J=pi*lamda/wc*exp(-(w/wc)).*w;% spectral density. In rad/ps

%J=J*1e12*hbar;                % spectral distribution function. In Joules
J=pi*lamda*(w*wc)./((wc^2+w.^2));

%w=w*1e12;                     % w in s-1
%dw=dw*1e12;                   % dw in s-1

%% 7. numerical parameters for Feynman integral
deltaKmax=3;                  % number of time steps before memory kernel dies
%totalT=10/1e12;               % total time for the simulation, in seconds
totalT=900/1e15;
%dt=totalT/200;                 % size of time step (delta t)
finalPoint=84;  % number of timesteps in total
allPointsORjustFinalPoint='allPoints'; % do you just want the density matrix at time=totalT ? or do you want it at every point until then
cpuORgpu='cpu'; 

%% 8. build input array that represents the densitry matrix (flattened to a column vector) as a function of time (where the columns represent the time steps)
rho=zeros(numel(H),finalPoint+1);
rho(:,1)=reshape(rhoInitial.',[],1); % I'm not sure about the transpose, since I'm not sure about labeling of the states in the column vector for the purposes of the final sum function.

%% 9. run the program and plot only <0|rho|0> and <1|rho|1> as a function of time
% wholeDensityMatrixOrJustDiagonals='justDiagonals';
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
figure(2);hold('on')
plot(0:totalT/finalPoint:totalT,real(rho_allElements(1,:)),'r','LineWidth',5);
plot(0:totalT/finalPoint:totalT,real(rho_allElements(9,:)),'g','LineWidth',5);
plot(0:totalT/finalPoint:totalT,real(rho_allElements(17,:)),'b','LineWidth',5);
plot(0:totalT/finalPoint:totalT,real(rho_allElements(25,:)),'m','LineWidth',5);
plot(0:totalT/finalPoint:totalT,real(rho_allElements(33,:)),'c','LineWidth',5);
plot(0:totalT/finalPoint:totalT,real(rho_allElements(41,:)),'y','LineWidth',5);
plot(0:totalT/finalPoint:totalT,real(rho_allElements(49,:)),'Color',[0.5 0.5 0.5],'LineWidth',5);
%plot(0:totalT/finalPoint:totalT,real(rho_allElements(3,:)),'r');plot(0:totalT/finalPoint:totalT,imag(rho_allElements(3,:)),'--r');
%plot(0:totalT/finalPoint:totalT,-real(rho_allElements(3,:)),'r');plot(0:totalT/finalPoint:totalT,imag(rho_allElements(3,:)),'--r');
%plot(0:totalT/finalPoint:totalT,real(rho_allElements(4,:)));
xlabel('time (seconds)');
 %save('rho.dat','rho_allElements','-ascii','-double');