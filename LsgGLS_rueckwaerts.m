function [x] = LsgGLS_rueckwaerts(R,z)
%R�ckw�rtssubstitution f�r die L�sung von Gleichungssystemen Rx=z
%   Input: Matrix R, Vektor z
%   Output: L�sungsvektor x

[n,~] = size(z);
x = zeros(n,1);
for k=n:-1:1
    if abs(R(k,k)) < eps
        error('Divison durch 0 nicht moeglich!')
    end
    sum = 0;
    for l=k+1:n
        sum = sum + R(k,l)*x(l);
    end
    x(k) = (z(k)-sum)/R(k,k);
end

end

