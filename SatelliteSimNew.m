% (c) Kristian Asti, July 2020
% For simulating a felxible satellite with two Euler-Bernoulli beams
% representing the solar panels and a rigid center body connecting them.

%% PARAMETERS AND MATRICES
% Physical constants for the beams: E, I, a, rho, gamma
E = 1;
I = 1;
a = 1;
rho = 2;
% IMPORANT! controls damping strength
gamma = 2;

% Physical constants for the rigid center body: m, Im
m = 1;
Im = 1;

% Number of basis functions used for spectral Galerkin: N
% (NOTE: values around N = 15 should be well enough)
N = 15;

% Left beam matrices: ML, KL, BL, CL in
%   ML*xL' = KL*xL + BL*uL
%       yL = CL*xL
[ML,KL,BL,CL,phi1_funs_L,phi2_funs_L] = BeamSimLeft2New(gamma,E,I,rho,a,N);

% Right beam matrices: MR, KR, BR, CR in
%   MR*xR' = KR*xR + BR*uR
%       yR = CR*xR
[MR,KR,BR,CR,phi1_funs_R,phi2_funs_R] = BeamSimRight2New(gamma,E,I,rho,a,N);

% Rigid center body matrices: Ic, Cc
Mc = eye(2);
Cc = [1/m 0; 0 1/Im];

% Total system matrices: Msys, Ksys in
%   Msys*x' = Ksys*x
Msys = [ML zeros(size(ML)) zeros(length(ML),2); ...
        zeros(size(ML)) MR zeros(length(ML),2); ...
        zeros(2,length(ML)) zeros(2,length(ML)) Mc];

Ksys = [KL zeros(size(KL)) BL*Cc; ...
        zeros(size(KL)) KR BR*Cc; ...
        -CL -CR zeros(2)];
    
Bsys = [zeros(length(ML),2); zeros(length(ML),2); eye(2)];
    
% System output = rigid center body output: ysys = Cc*xc




%% SIMULATION AND INITIAL CONDITIONS
k1 = fzero(@(k) cos(k) + 1./cosh(k),4);
KK1 = (sinh(k1) + sin(k1))/(cosh(k1) + cos(k1));

% --- INITIAL CONDITIONS FOR BEAM ENERGY VARIABLES ---
% Left beam:
%   x1L = rho*a*wL. (. means time derivative)
%   x2L = wL''
% Similarly, right beam:
%   x1R = rho*a*wR.
%   x2R = wR''

% LEFT BEAM INITIAL CONDITION FUNCTIONS. The "working domain" for our
% solvers is [0, 1], so adjust the initial condition accordingly. If your
% original initial condition is defined on [-1, 0], then substitute
% x -> x - 1 for the variable.
x1L_0 = @(x) 0*x;
x2L_0 = @(x) 0.1*(12*(x-1).^2 + 24*(x-1) + 12);

% Initial profile of the left solar panel
wL_0_fun = @(x) 0.1*((x-1).^4 + 4*(x-1).^3 + 6*(x-1).^2);

% RIGHT BEAM INITIAL CONDITION FUNCTIONS. Again, the "working domain"
% for our solvers is [0, 1]. If your original initial condition is defined
% on [0, 1], no need to change anything.
x1R_0 = @(x) 0*x;
x2R_0 = @(x) 0.1*(12*x.^2 - 24*x + 12);

% Initial profile of the right solar panel
wR_0_fun = @(x) 0.1*(x.^4 - 4*x.^3 + 6*x.^2);

% Initial conditions in terms of the basis function coefficients
alphaL_0 = basis_appr_01(x1L_0,phi1_funs_L);
betaL_0 = basis_appr_01(x2L_0,phi2_funs_L);

alphaR_0 = basis_appr_01(x1R_0,phi1_funs_R);
betaR_0 = basis_appr_01(x2R_0,phi2_funs_R);


% --- INITIAL CONDITIONS FOR RIGID CENTER BODY VARIABLES ---
% Variables:
%   x1c =  m*dwc        (linear displacement)
%   x2c = Im*dtheta_c   (angular displacement)

x1c_0 = 0;
x2c_0 = -1;


% --- Initial condition for the whole system variable ---
x_sys_0 = [alphaL_0; betaL_0; alphaR_0; betaR_0; x1c_0; x2c_0];


% --- CONTROL INPUT ---
K = [1 0; 0 1]; % Input adjustment matrix

u1 = @(t) 0*t; %-cos(2*pi*t);
u2 = @(t) 0*t; %0.25*sin(pi*t);
u = @(t) [u1(t); u2(t)];

u_sys = @(t) K*u(t);


% --- ACTUAL SIMULATION ---
% Time span in seconds [t_0 t_end]
%tspan = [0 20];
tspan = 0:0.002:20;
options = odeset('Mass',Msys);
[t, x_sys] = ode15s(@(t, x_sys) Ksys*x_sys + Bsys*u_sys(t), ...
                    tspan, x_sys_0, options);

% --- OUTPUTS ---
% We want to extract yc = [dw_c dtheta_c]^T = Cc * [x1c x2c]^T
%                       = [1/m*x1c 1/Im*x2c]^T,
% i.e. the last two columns of our output matrix x_sys

x1c_sol = 1/m*x_sys(:,end-1);
x2c_sol = 1/Im*x_sys(:,end);


% --- PLOTTING ---
% Edit the plotting parameters according to your needs.

t_plot = t; % Plot every computed time point
%t_plot = t(1:100:end); % Plot every 100th time point
x1c_sol_plot = x1c_sol;
%x1c_sol_plot = x1c_sol(1:100:end);
x2c_sol_plot = x2c_sol;
%x2c_sol_plot = x2c_sol(1:100:end);

figure
plot(t_plot,x1c_sol_plot,'LineWidth',2)
xlabel('Time [s]','FontSize',14)
grid on
hold on

plot(t_plot,x2c_sol_plot,'LineWidth',2)
hold off

legend({'$\frac{dw_c}{dt}(t)$','$\frac{d\theta_c}{dt}(t)$'},...
        'Interpreter','latex','FontSize',14,'Orientation','horizontal')


%% SURFACE PLOT OF SATELLITE PROFILE

xinodes_L = -1:0.01:0;
xinodes_R = 0:0.01:1;

phi1_L_matrix = zeros(length(xinodes_R),N);
phi1_R_matrix = zeros(length(xinodes_R),N);

for ii = 1:N
    phi1_L_matrix(:,ii) = polyval(phi1_funs_L{ii},xinodes_R);
    phi1_R_matrix(:,ii) = polyval(phi1_funs_R{ii},xinodes_R);
end

wL_0 = wL_0_fun(xinodes_R);
wR_0 = wR_0_fun(xinodes_R);

alpha_L_int = zeros(N,length(t));
alpha_R_int = zeros(N,length(t));

for n = 1:N
    alpha_L_n = x_sys(:,n);
    alpha_R_n = x_sys(:,2*N+n);
    alpha_L_int(n,:) = cumtrapz(t,alpha_L_n);
    alpha_R_int(n,:) = cumtrapz(t,alpha_R_n);
end

% Left and right solar panel profiles
wL_profile = 1/(rho*a)*phi1_L_matrix*alpha_L_int + wL_0.';
wR_profile = 1/(rho*a)*phi1_R_matrix*alpha_R_int + wR_0.';


% Graph as surfaces
surfstep = 80;

[T,XI_R] = meshgrid(t(1:surfstep:end),xinodes_R);
surf(T,XI_R,wR_profile(:,1:surfstep:end))
grid off

xlabel('$t$','Interpreter','latex','FontSize',16)
ylabel('$\xi$','Interpreter','latex','FontSize',16)
zlabel('$w(\xi,t)$','Interpreter','latex','FontSize',16)

hold on

[T,XI_L] = meshgrid(t(1:surfstep:end),xinodes_L);
surf(T,XI_L,wL_profile(:,1:surfstep:end))

hold off
