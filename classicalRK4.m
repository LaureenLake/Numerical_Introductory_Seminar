function [t,y] = KlassRK4(fun,tspan,y0,h)
% [t,y] = KlassRK4(fun,tspan,y0,h)
% Klassisches RK4-Verfahren mit konstanter Schrittweite

%INPUT
% fun              Funktion f der ODE y'=f(t,y) 
% tspan=[t0, tn]   Intervallgrenzen
% y0               Anfangswert (Vektor der Laenge m)
% h                Schrittweite (optional)

%OUTPUT
% t   Spaltenvektor, der die betrachteten Zeitpunkte enthält
% y   Matrix, die in Zeile l die Näherungslösungen zu y(t(l)) enthaelt

% Anfang der Zeitmessung
tic

% Initialisierung
t0 = tspan(1);
tn = tspan(end);

if ( nargin<4 )
    h = (tn - t0)/1000; % ggf. Berechtung von h
end

n  = ceil((tn - t0 )/h); % Anzahl Zeitschritte, aufgerundet, damit wir ggf. hinter tn kommen

% Speicherreservierung
[my,ny] = size(y0);
m       = max(my,ny);  % y0 kann als Zeilen- oder Spaltenvektor uebergeben werden

t = zeros(n,1);        
y = zeros(n,m);

% Anfangswerte setzen
t(1)   = t0;
y(1,:) = y0;

% Zeitschleife des klass. RK4-Verfahrens
for l = 1:n-1
    
    t(l+1)=t(l)+h;
    
    k1=fun(t(l),y(l,:)).';
    
    k2=fun(t(l)+h/2,y(l,:)+h*k1/2).';
    
    k3=fun(t(l)+h/2,y(l,:)+h*k2/2).';
    
    k4=fun(t(l)+h,y(l,:)+h*k3).';
    
    y(l+1,:)=y(l,:) + (h/6)* (k1+2*k2+2*k3+k4);
    
end

% Ende der Zeitmessung
toc