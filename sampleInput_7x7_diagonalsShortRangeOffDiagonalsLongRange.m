%% Originally I was having a lot of trouble getting my QUAPI to match the HEOM result of our JCTC. Sites 3-7 were somewhat fine, but Sites 1 and 2 were wrong. The problem was fixed when I discovered that Xu Hui had used Nalbach's hamitlonian rather than the JCTC one (diagonals were different, off-diagonals were still correct). This lead me to believe that diagonals are important for short-time dynamics and off-diagonals were important for the long-time dynamics. Actually I also remember the results being more convincing of this than what I'm getting now. Anyway, it seems that simply fixing the H(1,1) on Nalbach's Hamiltonian fixes everthing. Adding 210 to all diagonals did very little to the dynamics, and this left all diagonals matching the H from JCTC by +/- 15cm-1 except for H(1,1) whcih was off by 40. Setting H(1,1) to be equal and keeping the rest of the diagonals in +/- 15cm-1 discrepancy seemed to be all that was needed to make the dynamics the same.

%% 1. fundamental constants
kb=1.3806504*10^(-23);        % Joules / Kelvin
hbar=1.054571628*10^(-34);    % Joule * seconds

%% 2. temperature
temperature=77;
beta=1/(kb*temperature);

%% 3. system hamiltonian
wavenumbers2joules=1.9864475e-23;

H=[240.0 -87.7   5.5  -5.9   6.7 -13.7  -9.9   % Nalbach
     0.0 315.0  30.8   8.2   0.7  11.8   4.3
     0.0   0.0   0.0 -53.5  -2.2  -9.6   6.0
     0.0   0.0   0.0 130.0 -70.7 -17.0 -63.3
     0.0   0.0   0.0   0.0 285.0  81.1  -1.3
     0.0   0.0   0.0   0.0   0.0 435.0  39.7
     0.0   0.0   0.0   0.0   0.0   0.0 245.0];

H=[410.0 -87.7   5.5  -5.9   6.7 -13.7  -9.9   % Nalbach
     0.0 525.0  30.8   8.2   0.7  11.8   4.3
     0.0   0.0 210.0 -53.5  -2.2  -9.6   6.0
     0.0   0.0   0.0 320.0 -70.7 -17.0 -63.3
     0.0   0.0   0.0   0.0 495.0  81.1  -1.3
     0.0   0.0   0.0   0.0   0.0 645.0  39.7
     0.0   0.0   0.0   0.0   0.0   0.0 455.0]; 
 
% H=[410.0 -87.7   5.5  -5.9   6.7 -13.7  -9.9   % Wilkins-Dattani 2015 JCTC
%      0.0 530.0  30.8   8.2   0.7  11.8   4.3
%      0.0   0.0 210.0 -53.5  -2.2  -9.6   6.0
%      0.0   0.0   0.0 320.0 -70.7 -17.0 -63.3
%      0.0   0.0   0.0   0.0 480.0  81.1  -1.3
%      0.0   0.0   0.0   0.0   0.0 630.0  39.7
%      0.0   0.0   0.0   0.0   0.0   0.0 440.0];

 H=H+H'-diag(diag(H));
 H=wavenumbers2joules*H;

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
J=2*lamda*(w*wc)./((wc^2+w.^2));

%w=w*1e12;                     % w in s-1
%dw=dw*1e12;                   % dw in s-1

%% 7. numerical parameters for Feynman integral
deltaKmax=3;                  % number of time steps before memory kernel dies
totalT=1000/1e15;             % total time for the simulation, in seconds
%totalT=875/1e15;

finalPoint=50;                % number of timesteps in total
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
 %save('rho.dat','rho_allElements','-ascii','-double');