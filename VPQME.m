function drho = VPQME(t,rho,Kt,K,Ot,Omega,kb,hbar,Theta,tau,T,wc,alpha)
%% Interpolate functions

K = interp1(Kt,K,t); % Interpolate the data set (Kt,K) at time t
Omega = interp1(Ot,Omega,t); % Interpolate the data set (Ot,Omega) at time t

%% Master equation
drho = zeros (4,1);        
drho=[0 0 -Omega 0;0 -K 0 0; Omega 0 -K 0; 0 0 0 0]*rho+[0;0;-Omega/2;-rho(1)];     %Systems of differential equations.
end
