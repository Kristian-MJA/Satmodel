% (c) Kristian Asti, July 2020
function coeffs = legpol01(n)
%LEGPOL01(n) - returns the coefficients of the nth Legendre polynomial
% on the interval [0, 1].

coeffs = zeros(1,n+1);

for k = 0:n
    coeffs(n-k+1) = (-1)^(n+k)*nchoosek(n,k)*nchoosek(n+k,k);
end

end