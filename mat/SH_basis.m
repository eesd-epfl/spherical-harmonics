% Optimized to work efficiently with Matlab R2020a by Mahmoud S. Shaqfa
function Y_mat = SH_basis(vs, degree)
[phi,theta] = cart2sph(vs(:,1),vs(:,2),vs(:,3));
ind = find(phi<0); phi(ind) = phi(ind)+2*pi; theta = pi/2-theta;

% This can be faster- needs a bit of vectorisation
Y_mat = zeros(size(theta,1), (degree+1)^2);
for n = 0:degree
    Pn = legendre(n,cos(theta))';
    for m = -n:1:n
        norm = sqrt(((2*n+1) * factorial(n-abs(m)))/(4 * pi * factorial(n+abs(m))));
        if m >= 0
            Y_mat(:, n^2 + n + m + 1) = Pn(:, abs(m)+1) .* exp((m*1i).*phi) .* norm;
        else
            Y_mat(:, n^2 + n + m + 1) = (-1)^m .* (Pn(:, abs(m)+1) .* exp((m*1i).*phi) .* norm);
        end
    end
end