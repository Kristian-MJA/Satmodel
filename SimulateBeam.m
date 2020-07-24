% (c) Kristian Asti, July 2020
% Simulate a single beam on [0,1]
% Cantilevered beam clamped at xi = 0, with optional boundary control
% of (dw/dt)(0,t) = u_1(t) and (dw/(dxi dt))(0,t) = u_2(t).
% Possibility to plot the deflection curve w(xi,t) as a surface or as an
% animated plot.

gamma = 2;
E = 1;
I = 1;
rho = 10;
a = 1;
N = 15;

% BEAM MATRICES: M, K, B, in
%   M*x' = K*x + B*u
[M,K,B,~,phi1_funs,phi2_funs] = BeamSimSingle(gamma,E,I,rho,a,N);

% BEAM VARIABLES
%   x1 = rho*a*w.
%   x2 = w''

k1 = fzero(@(k) cos(k) + 1./cosh(k),2);
KK1 = (sinh(k1) + sin(k1))/(cosh(k1) + cos(k1));

% BEAM INITIAL CONDITION FUNCTIONS. The "working domain"
% for our solvers is [0, 1]. If your original initial condition is defined
% on [0, 1], no need to change anything.
x1_0 = @(x) 0*x;
%x2_0 = @(x) 0.1*k1^2*(-sin(k1*x)-sinh(k1*x)+KK1*(cosh(k1*x)+cos(k1*x)));
x2_0 = @(x) 0.1*(12*x.^2 - 24*x + 12);

% Initial conditions in terms of the basis function coefficients
alpha_0 = basis_appr_01(x1_0,phi1_funs);
beta_0 = basis_appr_01(x2_0,phi2_funs);

x_beam_0 = [alpha_0; beta_0];

% --- CONTROL INPUT ---
Kmat = [1 0; 0 1];  % Input adjustment matrix

u1 = @(t) 0*t;
u2 = @(t) -sin(2*pi*t); %0*t;
u = @(t) [u1(t); u2(t)];

u_beam = @(t) Kmat*u(t);

%% SIMULATION PART
% Time span in seconds [t_0 t_end]
%tspan = [0 20];
tspan = 0:0.1:20;
options = odeset('Mass',M);
[t, x_beam] = ode15s(@(t, x_beam) K*x_beam + B*u_beam(t), ...
                    tspan, x_beam_0, options);

xinodes = 0:0.01:1; % For plotting the solution          
phi1_matrix = zeros(length(xinodes),N);

for ii = 1:N
    phi1_matrix(:,ii) = polyval(phi1_funs{ii},xinodes);
end

% Solving for the deflection of the beam using trapezoidal rule
%w0_fun = @(x) 0.1*(sin(k1*x)-sinh(k1*x)+KK1*(cosh(k1*x)-cos(k1*x)));
w0_fun = @(x) 0.1*(x.^4 - 4*x.^3 + 6*x.^2);
w_beam_0 = w0_fun(xinodes);

alpha_int = zeros(N,length(t));

for n = 1:N
    alpha_n = x_beam(:,n);
    alpha_int(n,:) = cumtrapz(t,alpha_n);
end

w_beam = 1/(rho*a)*phi1_matrix*alpha_int + w_beam_0.';

%% SURFACE PLOT

[T,XI] = meshgrid(t(1:1:end),xinodes);
surf(T,XI,w_beam(:,1:1:end))
grid off

xlabel('$t$','Interpreter','latex','FontSize',16)
ylabel('$\xi$','Interpreter','latex','FontSize',16)
zlabel('$w(\xi,t)$','Interpreter','latex','FontSize',16)
% title([num2str(N) ' basis functions, boundary input $\rho a\frac{\partial w}{\partial\xi \partial t}(0,t) = -0.1\sin(2\pi t)$'],...
%         'Interpreter','latex','FontSize',16)
% title([num2str(N) ' basis functions, boundary input $\frac{1}{\rho a}\frac{\partial w}{\partial\xi \partial t}(0,t) \equiv 0$'],...
%         'Interpreter','latex','FontSize',16)


%% ANIMATION
for tind = 1:1:length(t)
    tnow = t(tind);
    
    if tind < length(t)
        tnext = t(tind+1);
        dt = tnext - tnow;
    end
    
    ysol = w_beam(:,tind);
    
    plot(xinodes,ysol,'LineWidth',2)
    
    axis([min(xinodes) max(xinodes) -0.3 0.3])
    grid on
    title(sprintf('Time = %1.3f',tnow),'FontSize',16)
    xlabel('$\xi$','Interpreter','latex','FontSize',16)
    ylabel('$w(\xi,t)$','Interpreter','latex','FontSize',16)
    
    if(tind == 1)
        pause(3)
    else
        pause(dt)
    end
end         
                