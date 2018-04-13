function drho = VPQME(t,rho)
Omega=1;                                                   % Pulse area      
tau=14;                                                    % Full width half maximum
K=3;                                                       % Spectral density (On construction)
drho(1)=1-((Omega/(2*tau*sqrt(pi)))*exp(-(t./tau)^2))*rho(3);                    % p00 
drho(2)=-K*rho(2);                                                               % Re(p10)  
drho(3)=((Omega/(2*tau*sqrt(pi)))*exp(-(t./tau)^2))*(1/2)*(1-2*rho(1))+K*rho(3); % Re(p10)
drho(4)=1-drho(1);
end