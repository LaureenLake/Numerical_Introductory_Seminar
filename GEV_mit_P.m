function [x] = GEV_mit_P(A,b)
%Gaußsches Eliminationsvefahren für Ax=b mit Spaltenpivotisierung
%   Input: Matrix A, Vektor b
%   Output: Lösungsvektor x

[na,~]= size(A);
[nb,~] = size(b);
if na ~= nb
    error('Die Dimensionen von A und b stimmen nicht überein!')
end

[n,~] = size(b);
h = zeros(n,1);
h = b;

[A,p] = LR_mit_P(A);
 for k=1:n
     % Zeilen von b m.H. des Permutationsvektors vertauschen
     b(k) = h(p(k));
 end
z = LsgGS_vorwaerts(A,b);
x = LsgGLS_rueckwaerts(A,z);

end

