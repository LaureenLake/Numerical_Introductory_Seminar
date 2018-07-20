function [t,y] = odeEulerExp(f,tIntervall,y0,varargin)
%explizites Euler-Verfahren 
t0=tIntervall(1);
tN=tIntervall(2);

if nargin < 4
    h = 1/1000 * (tN-t0);
    n= 1001;
elseif nargin == 4
    h= varargin{1};
    n=(tN-t0)/h +1;
end

n  = ceil((tN - t0 )/h); % Anzahl Zeitschritte, aufgerundet, damit wir ggf. hinter tn kommen

% Speicherreservierung
[my,ny] = size(y0);
m       = max(my,ny);  % y0 kann als Zeilen- oder Spaltenvektor uebergeben werden

t = zeros(n,1);        
y = zeros(n,m);

% Anfangswerte setzen
t(1)   = t0;
y(1,:) = y0;
% t=(t0:h:tN)';
% y=zeros(n,1);
% y(1)=y0;

for k=1:n-1
    y(k+1,:)=y(k,:)+h*f(t(k),y(k,:))';
end
     
end

