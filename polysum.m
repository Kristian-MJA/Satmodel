% (c) Kristian Asti 2020
function psum = polysum(varargin)
%polysum(P1,P2,...,PN) - returns the sum of N polynomials P1, P2, ... ,PN,
% represented as coefficient (row) vectors compatible with polyval,
% polyder etc.

N = nargin;

maxlength = 0;

% Search for the longest vector in varargin
for k = 1:N
    input_k = varargin{k};
    if length(input_k) > maxlength
        maxlength = length(input_k);
    end
end

polymat = zeros(N,maxlength);

for k = 1:N
    input_k = varargin{k};
    len_k = length(input_k);
    
    polymat(k,:) = [zeros(1,maxlength - len_k) input_k];
end

if N == 1
    psum = polymat;
else
    psum = sum(polymat);
end

% Remove leading zeros
ind = find(psum ~= 0,1,'first');
psum = psum(ind:end);

end