function [z] = LsgGS_vorwaerts(L,b)
%Vorwärtssubstitution für die Lösung von Gleichungssystemen Lz=b
%   Input: Matrix L, Vektor b
%   Output: Lösungsvektor z

[n,~] = size(b);
z = zeros(n,1);
for k=1:n
    if abs(L(k,k)) < eps
        error('Matrix ist singulär')
    end
    sum = 0;
    for l=1:k-1
        sum = sum + L(k,l)*z(l);
    end
    z(k) = (b(k)-sum);
end
end

