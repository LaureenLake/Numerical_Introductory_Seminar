x0= [1.2 0 0 -1.04935750983035];
h= 0.002;

[t,X1]= KlassRK4(@Satellit,[0,20],x0,h);
figure(3)
subplot(1,3,1)
plot(X1(:,1),X1(:,2))
title('Klassisches RK-Verfahren (Ordnung 4)')

[t,X2]= expl_RK5(@Satellit,[0,20],x0,h);
subplot(1,3,2)
plot(X2(:,1),X2(:,2))
title('Explizites RK-Verfahren (Ordnung 5)')

%Referenzlösung bestimmen für Bsp. mit h=0.002 und tn=20

options= odeset('RelTol',1e-6,'AbsTol',[1e-10 1e-10 1e-10 1e-10]);
[t3,X3]= ode45('Satellit',[0:h:20-h],x0,options);
subplot(1,3,3)
plot(X3(:,1),X3(:,2))
title('Referenzlösung mit ode45 (adaptive Schrittweitensteuerung)')

pause
errRK4= abs(X3-X1);
errRK5= abs(X3-X2);

figure(4)
subplot(2,1,1)
plot(t,errRK4(:,1),t,errRK4(:,2))
title('globaler Fehler des RK4')
xlabel('t');
ylabel('glob. Fehler'); 
legend('ersteKomponente','zweite Komponente')

subplot(2,1,2)
plot(t,errRK5(:,1),t,errRK5(:,2))
title('globaler Fehler des RK5')
xlabel('t');
ylabel('glob. Fehler'); 
legend('ersteKomponente','zweite Komponente')

