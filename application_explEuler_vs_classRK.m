
x0= [1.2 0 0 -1.04935750983035];
h= 0.001;

[t,X]= KlassRK4(@Satellit,[0,8],x0,h);
figure(1)
plot(X(:,1),X(:,2))
title('Klassisches RK-Verfahren (Ordnung 4)')

pause

[t,X2]=odeEulerExp(@Satellit,[0,8],x0,h);
figure(2)
plot(X2(:,1),X2(:,2))
title('Explites Euler-Verfahren (Ordnung 1)')


