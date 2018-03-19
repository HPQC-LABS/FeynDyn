%% Bearing in mind the effort that went into the development of FeynDyn,
%% the author would indeed appreciate if users can please cite:
%% 1) The paper: N. Dattani (2013) Comp. Phys. Comm. Volume 184, Issue 12, Pg. 2828-2833 , and
%% 2) The code:  N. Dattani (2013) FeynDyn. http://dx.doi.org/10.6084/m9.figshare.823549

%% To make sure your code is the most updated version, please e-mail dattani.nike@gmail.com
%% Bug reports, suggestions, and requests for extensions are more than encouraged: dattani.nike@gmail.com

%% This program is protected by copyright law. Nike S. Dattani , 13/8/12

%% This sample input file, upon pressing F5, plots the 7 diagonals <i|rho|i> as in Fig. 1 of 
% 2015 Wilkins & Dattani, J. Chem. Thy. Comp.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. fundamental constants
kb=1.3806504*10^(-23);                          % Joules / Kelvin
hbar=1.054571628*10^(-34);                      % Joule * seconds

%% 2. temperature
temperature=77;                                 % Kelvin
beta=1/(kb*temperature);

%% 3. system hamiltonian
wavenumbers2joules=1.9864475e-23;

H=[410.0 -87.7   5.5  -5.9   6.7 -13.7  -9.9   % H in wavenumbers (cm-1)
     0.0 530.0  30.8   8.2   0.7  11.8   4.3
     0.0   0.0 210.0 -53.5  -2.2  -9.6   6.0
     0.0   0.0   0.0 320.0 -70.7 -17.0 -63.3
     0.0   0.0   0.0   0.0 480.0  81.1  -1.3
     0.0   0.0   0.0   0.0   0.0 630.0  39.7
     0.0   0.0   0.0   0.0   0.0   0.0 440.0];

 H=H+H'-diag(diag(H));                         % Fill lower triangle, subtract diagonals to avoid doubling
 H=wavenumbers2joules*H;                       % H in Joules

%% 4. initial denisty matrix
rhoInitial=zeros(size(H));                     
rhoInitial(1)=1;                               % Excitation starts localized at site 1

%% 5. system coupling matrix
Nbaths=7;                                      % Number of baths
systemCouplingMatrix=diag(ones(Nbaths,1));     % Each site i couples to its bath via |iXi| coupling

%% 6. spectral distribution function
wc=1/(50e-15);                                 % cutoff frequency. In rad/ps
dw=0.00001*wc;                                 % stepsize for w in J(w). In rad/ps
w=dw:dw:50*wc;                                 % must start with dw because coth(0)=infinity. In rad/ps
lamda=35*wavenumbers2joules;                   % reorganization energy in Joules
J=2*lamda*(w*wc)./((wc^2+w.^2));               % spectral density in Joules
     
%% 7. numerical parameters for Feynman integral
deltaKmax=3;                                   % number of time steps before memory kernel dies
totalT=1000/1e15;                              % total time for the simulation, in seconds
finalPoint=50;                                 % number of timesteps in total
allPointsORjustFinalPoint='allPoints';         % do you just want the density matrix at time=totalT ? or do you want it at every point until then
cpuORgpu='cpu'; 

%% 8. build input array that represents the densitry matrix (flattened to a column vector) as a function of time (where the columns represent the time steps)
rho=zeros(numel(H),finalPoint+1);
rho(:,1)=reshape(rhoInitial.',[],1);          % I'm not sure about the transpose, since I'm not sure about labeling of the states in the column vector for the purposes of the final sum function.

%% 9. run the program and plot all diagoal elements of the density matrix as a function of time
wholeDensityMatrixOrJustDiagonals='justDiagonals';
[rho,elapsedTime]=FeynDynCode(Nbaths,finalPoint,deltaKmax,totalT,rho,H,systemCouplingMatrix,w,dw,J,temperature,wholeDensityMatrixOrJustDiagonals,allPointsORjustFinalPoint,cpuORgpu);

%%
figure(1);hold('on');box('on')
axis([0 totalT 0 1])
mesh=0:totalT/finalPoint:totalT;
splineMesh=0:totalT/finalPoint/100:totalT;
plot(splineMesh,spline(mesh,real(rho(1,:)),splineMesh),'r','LineWidth',5);
plot(splineMesh,spline(mesh,real(rho(9,:)),splineMesh),'g','LineWidth',5);
plot(splineMesh,spline(mesh,real(rho(17,:)),splineMesh),'b','LineWidth',5);
plot(splineMesh,spline(mesh,real(rho(25,:)),splineMesh),'m','LineWidth',5);
plot(splineMesh,spline(mesh,real(rho(33,:)),splineMesh),'c','LineWidth',5);
plot(splineMesh,spline(mesh,real(rho(41,:)),splineMesh),'y','LineWidth',5);
plot(splineMesh,spline(mesh,real(rho(49,:)),splineMesh),'Color',[0.5 0.5 0.5],'LineWidth',5);
%plot(0:totalT/finalPoint:totalT,real(rho(3,:)),'r');plot(0:totalT/finalPoint:totalT,imag(rho(3,:)),'--r');
xlabel('time (seconds)');