function dy = phononm(t,y)
Omega=1     %Pulse area      
tau=14      %Full width half maximum
K=3         %Spectral density (On construction)
dy = zeros (3,1);
dy(1)=-((Omega/(2*tau*sqrt(pi)))*exp(-(t./tau)^2))*y(3);    %p11 
dy(2)=-K*y(2);                                              %Re(p10)  
dy(3)=((Omega/(2*tau*sqrt(pi)))*exp(-(t./tau)^2))*(1/2)*(1-2*y(1))+K*y(3);   %Re(p11)
end

