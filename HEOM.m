function rho=HEOM(rho,rho_now,dt,t,H,V,gamma,ThetaRe,ThetaIm,numberOfLayers,hbar) % does inputting rho_now create another copy of it for this function ? If so you are doubling the memory

for timeStep=2:length(t)
    for k=2+(0:numberOfLayers-1); n=k-2;                                   % written this way becuase 0th layer is at n=2. Only to \cal{N}-1 bcause the last value of the index is a dummy for the n+1
         rho_now(:,:,k)=rho_now(:,:,k)+dt*(-1i/hbar)*(... 
                  -n*(ThetaRe*rho_now(:,:,k-1)-rho_now(:,:,k-1)*ThetaRe... % 0th layer of the hierarchy, done separately since there is not (n-1)th layer of the hierarchy.
                 -1i*(ThetaIm*rho_now(:,:,k-1)+rho_now(:,:,k-1)*ThetaIm))...
                           +H*rho_now(:,:,k)-rho_now(:,:,k)*H + n*gamma*hbar*rho_now(:,:,k)... 
                          -(V*rho_now(:,:,k+1)-rho_now(:,:,k+1)*V)); 
    end
 
rho(:,:,timeStep)=rho_now(:,:,2);
end