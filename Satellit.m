function dx = Satellit(t,x)
% DGL der Umlaufbahn eines Satelliten um Erde und Mond

%rel .Mondmasse
mu=1/82.45; 

dx=[x(3);
    x(4);
    x(1)+2*x(4)-(1-mu)*(x(1)+mu)/((x(1)+mu)^2 +x(2)^2)^(3/2)-mu*(x(1)-1+mu)/((x(1)-1+mu)^2+x(2)^2)^(3/2);
    x(2)-2*x(3)-(1-mu)*x(2)/((x(1)+mu)^2 +x(2)^2)^(3/2)-mu*x(2)/((x(1)-1+mu)^2+x(2)^2)^(3/2)];