% (c) Kristian Asti, July 2020
function [M,K,B,C,phi1_funs,phi2_funs] = BeamSimLeft2New(gamma,E,I,rho,a,N)
%BeamSimLeft2New(gamma,E,I,rho,a,N) - computes the system matrices for
% simulating an Euler-Bernoulli beam on [-1,0] with a free end at xi = -1
% and a controlled end at xi = 0. Also returns the computed basis functions
% in a cell array. The controlled boundary inputs are:
%   x_1(0,t) and x_1'(0,t), where x_1(xi,t) = (d/dt wL)(xi,t).
% With constant zero inputs, this corresponds to a cantilevered beam
% clamped at xi = 0.
%
% NOTE: this function uses the tools developed for a simulating a beam on
% [0,1]. Thus in this function, xi = 0 acts as xi = -1, similarly xi = 1
% acts as xi = 0, essentially moving the beam to the right by 1.



% Left end basis, phi2(0) = phi2'(0) = 0
phi2 = @(k) polysum(legpol01(k),(2*k+3)/(k+2)*legpol01(k+1),...
                    (k+1)/(k+2)*legpol01(k+2));
        
phi2d = @(k) polyder(phi2(k));
phi2dd = @(k) polyder(polyder(phi2(k)));

% Right end basis, phi1(1) = phi1'(1) = 0    
phi1 = @(k) polysum(legpol01(k),-(2*k+3)/(k+2)*legpol01(k+1),...
                    (k+1)/(k+2)*legpol01(k+2));
        
%phi1d = @(k) polyder(phi1(k));
%phi1dd = @(k) polyder(polyder(phi1(k)));

% Boundary control basis functions at xi = 1
% phi1_u11 = phi1(N-2); % Standard basis function N-2
% phi1_u12 = phi1(N-1); % Standard basis function N-1
phi1_u11 = [-3 6 -2];
phi1_u12 = [1.5 -2 0.5];


% weight function w = 1
%ip_w = @(p1,p2) diff(polyval(polyint(conv(p1,p2)),[0 1]));
% NUMERIC VERSION, WORKS BETTER
ip_w = @(p1,p2) max(abs(p1))*max(abs(p2)) * integral(@(x) ...
            polyval(p1/max(abs(p1)),x).*polyval(p2/max(abs(p2)),x),...
            0,1,'RelTol',0,'AbsTol',1e-14);


M_alpha = zeros(N);
M_beta = zeros(N);

K_alpha1 = zeros(N);
K_alpha2 = zeros(N);
K_beta1 = zeros(N);

B_alpha = zeros(N,2);
B_beta = zeros(N,2);

C_1 = zeros(2,N);
C_2 = zeros(2,N);

phi1_funs = cell(N,1);
phi2_funs = cell(N,1);


for k = 0:N-1
    for m = 0:N-1
        % Basis and test functions corresponding to xi = 1
        
        if k == N - 1
            phi_k1 = phi1_u12;
        elseif k == N - 2
            phi_k1 = phi1_u11;
        else
            phi_k1 = phi1(k);
        end
        
        if m == N - 1
            psi_m1 = phi1_u12;
        elseif m == N - 2
            psi_m1 = phi1_u11;
        else
            psi_m1 = phi1(m);
        end
            
        % Basis and test functions corresponding to xi = 0
        phi_k2 = phi2(k);
        phid_k2 = phi2d(k);
        phidd_k2 = phi2dd(k);
        
        psi_m2 = phi2(m);
        psid_m2 = phi2d(m);
        psidd_m2 = phi2dd(m);
        
        M_alpha(m+1,k+1) = ip_w(phi_k1,psi_m1);
        M_beta(m+1,k+1) = ip_w(phi_k2,psi_m2);
        
        K_alpha1(m+1,k+1) = -gamma/(rho*a)*ip_w(phi_k1,psi_m1);
        K_beta1(m+1,k+1) = -E*I*ip_w(phidd_k2,psi_m1);
        K_alpha2(m+1,k+1) = 1/(rho*a)*ip_w(phi_k1,psidd_m2);
        
        B_beta(m+1,:) = 1/(rho*a)*[polyval(-psid_m2,1) polyval(psi_m2,1)];
    end
    
    phi1_funs{k+1} = phi_k1;
    phi2_funs{k+1} = phi_k2;
    
    C_2(:,k+1) = E*I*[-polyval(phid_k2,1); polyval(phi_k2,1)];
    
end

% M*x' = K*x + B*u
M = [M_alpha zeros(N); zeros(N) M_beta];
K = [K_alpha1 K_beta1; K_alpha2 zeros(N)];
B = [B_alpha; B_beta];

% y = Cx
C = [C_1 C_2];

end