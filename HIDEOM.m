clear('all')
%% 1. fundamental constants
kb=1.3806504*10^(-23);        % Joules / Kelvin
hbar=1.054571628*10^(-34);    % Joule * seconds

%% 2. temperature
temperature=0;
beta=1/(kb*temperature);

%% 3. system hamiltonian
Omega=hbar*1e12*pi/8;         % Rabi frequency in Joules
H=[0 Omega/2; Omega/2 0];     % System hamiltonian in Joules
M=length(H);

V=[1 0 ; 0 -1];

%% 4. spectral distribution function
alpha=0.027*pi;               % prefactor. In ps^2 / rad^2
alpha=0.00054*pi;
wc=2.2;                       % cutoff frequency. In rad/ps
dw=0.1;                       % stepsize for w in J(w). In rad/ps
w_beforeExpand=dw:dw:14;      % must start with dw because coth(0)=infinity. In rad/ps
w_beforeExpand=[-fliplr(w_beforeExpand) w_beforeExpand];   % need to get rid of middle value (w=0 occurs twice)
w=w_beforeExpand;
%w=expand(w_beforeExpand,[M,M]);  % seems it can run out of memory before the HEOM even begins!

J=alpha*exp(-(w/wc).^2).*w.^3;% spectral density. In rad/ps

J=J*1e12*hbar;                % spectral distribution function. In Joules
%J=0*J;

w=w*1e12;                     % w in s-1
dw=dw*1e12;                   % dw in s-1

ThetaRe=(1/hbar)*J.*coth(beta*hbar*w/2); % in s-1
ThetaIm=(1/hbar)*J;                      % in s-1

%H=repmat(H,1,length(w)/M);
%V=repmat(V,1,1,length(w)/M);

%% 5. time mesh
finalPoint=300;           % number of timesteps in total
totalT=10/1e12;            % total time for the simulation, in seconds
dt=totalT/finalPoint;
t=0:dt:totalT;

%% 6. intialize rho and ADOs
rho=zeros(M,M,length(t));rho(1)=1;
rho_w1=zeros(M,M,length(w));
rho_w1_next=rho_w1;
integrand=rho_w1;

rho_w1w2=zeros(M,M,length(w),length(w));
rho_w1w2_next=rho_w1w2;
integrand2=rho_w1w2;

%% 7. propagation
for ii=1:length(t)
    
    for jjj=1:length(w)
        for jj=1:length(w)
            rho_w1w2_next(:,:,jj,jjj)=rho_w1w2(:,:,jj,jjj)-(dt*1i/hbar)*(H*rho_w1w2(:,:,jj,jjj)-rho_w1w2(:,:,jj,jjj)*H+hbar*(w(jj)+w(jjj))*rho_w1w2(:,:,jj,jjj)+hbar*ThetaRe(jjj)*(V*rho_w1(:,:,jj)-rho_w1(:,:,jj)*V)+hbar*ThetaIm(jjj)*(V*rho_w1(:,:,jj)+hbar*ThetaRe(jj)*(V*rho_w1(:,:,jj)-rho_w1(:,:,jj)*V)+hbar*ThetaIm(jj)*(V*rho_w1(:,:,jj)+rho(:,:,jj)*V)));
            integrand2(:,:,jj,jjj)=V*rho_w1w2(:,:,jj,jjj)-rho_w1w2(:,:,jj,jjj)*V;
        end
    end
    
    for jj=1:length(w)
        rho_w1_next(:,:,jj)=rho_w1(:,:,jj)-(dt*1i/hbar)*(hbar*ThetaRe(jj)*(V*rho(:,:,ii)-rho(:,:,ii)*V)+hbar*ThetaIm(jj)*(V*rho(:,:,ii)+rho(:,:,ii)*V)+H*rho_w1(:,:,jj)-rho_w1(:,:,jj)*H+hbar*w(jj)*rho_w1(:,:,jj)+hbar*trapz(w,squeeze(integrand2(:,:,jj,:)),3));
        integrand(:,:,jj)=V*rho_w1(:,:,jj)-rho_w1(:,:,jj)*V;
    end
    
    rho(:,:,ii+1)=rho(:,:,ii)-(dt*1i/hbar)*(H*rho(:,:,ii)-rho(:,:,ii)*H+hbar*(trapz(w,integrand,3)));
    
    rho_w1=rho_w1_next;
    rho_w1w2=rho_w1w2_next;
    disp(ii)
end
plot(t,real(squeeze(rho(1,1,2:end))'));
axis([0 t(end) 0 1])