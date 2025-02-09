% Geometric Parameters
L = 300; % Length of macroscale domain

% Mesh size parameters - let each solver decide this independently!
%N = 100^2; % Total number of grid points to use

%Test commit!

% Time/space scale for simulation
T = 1e5; 
tspan = linspace(0,T,1e5); % Interpolate solution on [0,T] with 1e2 points
x = linspace(0,L,N);

% Diffusion parameters
D = 100;% D = 100;
d1=D; d2=1; d3=D;
 
% Kinetic parameters
rho_u  =  0.692; rho_w  =  2.5; gamma_v  =  0.1; gamma_w  =  0.001;
mu_u  =  0.167; mu_w  =  55.56; sigma_u  =  0.01; sigma_w  =  0.5; 
alpha=0.07;

% Kinetic functions
f = @(u,v,w)alpha*v-mu_u*u+rho_u*u.*w./(1+max(w,0))+sigma_u;
g = @(u,v,w)v.*(1-v)-u.*v./(gamma_v+max(v,0));
h = @(u,v,w)rho_w*u.*v./(gamma_w+max(v,0))-mu_w*w+sigma_w;
K = @(t)0.03*t/T;

% Index labels
uN = 1:N; vN = N+1:2*N; wN = 2*N+1:3*N;

% Steady states 
uss = 0.299329353079108; vss = 0.506308497768032; wss = 0.022441474032573;


% Create a normally distributed vector of size 3*N, mean 1, stdev 1e-2.
rng(1);
rand_vec = abs(1+1e-1*randn(3*N,1)); 

% Initial conditions for ODE system
uvH_init = [uss*rand_vec(uN);vss*rand_vec(vN); wss*rand_vec(wN)];
%uvH_init = [U(end,uN)'; U(end,vN)'; U(end,wN)'];