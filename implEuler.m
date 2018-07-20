function [t,y,anz] = implEuler(f,Jf,tIntervall,y0,h,TOL,IterMax)
% implizites (nichtlineares) Euler-Verfahren mit konstanter 
% Schrittweite h: y_n = y_(n-1) + h*f(t_n,y_n)

%Input: rechte Seite f,Jacobimatrix von f,[t0 tend], Startwert y0,
%       konstante Schrittweite h, sowie Fehlertoleranz und Itermax       
%       für das Newton-Verfahren

%Output: t - Vektor mit den Stellen, an denen ausgewertet wurde
%        y - Lösungsvektor der DGL an den Stellen von t
%        anz - Anzahl an Newtoniterationen für jeden Zeitpunkt

t0=tIntervall(1);
tN=tIntervall(2);

n  = ceil((tN - t0 )/h); % Anzahl Zeitschritte, aufgerundet, damit wir ggf. hinter tn kommen

% Speicherreservierung
[my,ny] = size(y0);
m       = max(my,ny);  % y0 kann als Zeilen- oder Spaltenvektor uebergeben werden
t   = (t0:h:t0+(n-1)*h)';        
y   = zeros(n,m);
anz = zeros(n,1);
I   = eye(length(y0));

% Anfangswerte setzen
y(1,:) = y0;

for k=2:n  
    y0=y(k-1,:);
    y(k,:) = y0;
    fehler = 2*TOL;  %Residuum
    
    while( anz(k) <= IterMax && fehler >= TOL ) 
    %Auswertung der Funktion und der Jacobi-Matrix:
        funy     = y(k,:)'-y0'-h*f(t(k),y(k,:)');
        Jfuny    = I-h*Jf(t(k),y(k,:)');         
        deltay   = GEV_mit_P(Jfuny,-funy);                   
        y(k,:)   = y(k,:)+ deltay';
        fehler   = norm(funy);

        if( anz(k) == IterMax && fehler >= TOL)
            disp('Maximale Anzahl an Iterationen ist erreicht. Das Residuum beträgt: ')
            disp('f(t,y) =')
            disp(funy);
            break
        end
        anz(k)= anz(k)+1;
    end 
end
end

