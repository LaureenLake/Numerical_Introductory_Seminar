function [t,y] = expl_RK5(fun,tspan,y0,h)
%explizites RK-Verfahren der Ordung 5 mit konst. Schrittweite h

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

% Zeitschleife des expl. RK-Verfahrens
for l = 1:n-1
    
    t(l+1)=t(l)+h;
    
    k1=fun(t(l),y(l,:)).';
    
    k2=fun(t(l)+h/5,y(l,:)+h*k1/5).';
    
    k3=fun(t(l)+h*3/10,y(l,:)+h*(3/40*k1+9/40*k2)).';
    
    k4=fun(t(l)+h*4/5,y(l,:)+h*(44/45*k1-56/15*k2+32/9*k3)).';
    
    k5=fun(t(l)+h*8/9,y(l,:)+h*(19372/6561*k1-25360/2187*k2+64448/6561*k3-212/729*k4)).';
    
    k6=fun(t(l+1),y(l,:)+h*(9017/3168*k1-355/33*k2+46732/5247*k3+49/176*k4-5103/18656*k5)).';
    
    y(l+1,:)=y(l,:) + h*(35/384*k1+500/1113*k3+125/192*k4-2187/6784*k5+11/84*k6);
    
end

% Ende der Zeitmessung
toc
end

