% (c) Kristian Asti, July 2020
function ak = basis_appr_01(fun,phifuns)
%basis_appr(fun,phifuns) - returns the coefficients a_k for a basis
% function approximation. If phi_k, k = 0,...,N-1, are the basis functions,
% then we want that
%   f(x) ~ a_0 * phi_0(x) + a_1 * phi_1(x) + ... + a_{N-1] * phi_{N-1}(x),
% where f is the function to be approximated. This functions computes the
% approximate coefficients by solving a linear system obtained by computing
% values of f at N evenly spaced points on (0,1). This is numerically more
% stable than computing definite integral numerically.
%
% fun: function handle with a single argument, defined on [0,1].
%
% phifuns: cell array containing the basis functions defined on [0,1].

N = length(phifuns);

xk = linspace(0,1,N+2).';
xk = xk(2:end-1);

funvector = fun(xk);
phimat = zeros(N);

for k = 1:N
    phimat(:,k) = polyval(phifuns{k},xk);
end

ak = phimat\funvector;

end