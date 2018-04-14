function drho = phononm(t,rho)
Theta=30                                                   
tau=14e-12                                            %Full width half maximum in s
hbar=1.0545718e-34                                    %Planck constant  J.s
T=25                                                  %Temperature in K
kb=1.38064852e-23                                     %Boltzmann constant J/K
wc=2.2e12                                             %1/s
A=pi*0.027*1e-24                                      %s^2
Omega=(Theta/(2*tau*sqrt(pi)))*exp(-(t./(2*tau))^2)    %1/s
K=coth(hbar*Omega/(2*kb*T))*A*((Omega)^3)*exp(-(Omega/wc)^2)*pi/2   %Spectral density                 


drho = zeros (3,1);        

drho=[0 0 -Omega;0 -K 0; Omega 0 -K]*rho+[0;0;-Omega/2]

end
