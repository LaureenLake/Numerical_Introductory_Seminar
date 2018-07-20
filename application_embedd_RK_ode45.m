x0= [1.2 0 0 -1.04935750983035];
h= 0.003;
t0=0;
tn=50;

[t,X,err]= expl_RK5bett(@Satellit,[t0,tn],x0,h);

subplot(2,1,1)
plot(t,err(:,1))
xlabel('t')
ylabel('lokaler Fehler')
legend('erste Komponente')

subplot(2,1,2)
plot(t,err(:,2),'r')
xlabel('t')
ylabel('lokaler Fehler')
legend('zweite Komponente')