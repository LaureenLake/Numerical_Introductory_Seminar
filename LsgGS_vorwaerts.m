function [z] = LsgGS_vorwaerts(L,b)
%Vorw�rtssubstitution f�r die L�sung von Gleichungssystemen Lz=b
%   Input: Matrix L, Vektor b
%   Output: L�sungsvektor z

[n,~] = size(b);
z = zeros(n,1);
for k=1:n
    if abs(L(k,k)) < eps
        error('Matrix ist singul�r')
    end
    sum = 0;
    for l=1:k-1
        sum = sum + L(k,l)*z(l);
    end
    z(k) = (b(k)-sum);
end
end

