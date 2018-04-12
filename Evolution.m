options=odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4]);
[T,Y]=ode113(@phononm,[0 1],[0.5 0.5 0],options);    %Time and initial conditions
plot(T,Y(:,1),'-',T,Y(:,2),'-.',T,Y(:,3),'.')        %Plot
