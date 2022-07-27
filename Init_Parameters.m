% Geometric Parameters
L = 1; % Length of macroscale domain

% Mesh size parameters - let each solver decide this independently!
%N = 100^2; % Total number of grid points to use


% Time/space scale for simulation
T = 15000; 
tspan = linspace(0,T,1e3); % Interpolate solution on [0,T] with 1e2 points
x = linspace(0,L,N);

% Diffusion parameters
%d1=1; d2=0.0001; d3=1;
d1=0.001; d2=0.0000199; d3=0.01;
 
% Kinetic parameters
c = 0.01; muu=0.167; pu=0.69167; gu=20;su=0.0250;
rv=1;b=1;pv=0.5555556;gv=0.1;
pw=27.778;gw=0.001;muw=55.55556;sw=46.5327;


% Kinetic functions
f = @(u,v,w)c*v-muu*u+pu*u.*w./(gu+w)+su;
g = @(u,v,w)rv*v.*(1-b*v)-pv*u.*v./(gv+v);
h = @(u,v,w)pw*u.*v./(gw+v)-muw*w+sw;

% Index labels
uN = 1:N; vN = N+1:2*N; wN = 2*N+1:3*N;

% Steady states - TODO NEED TO DISCUSS STEADY STATE ANALYSIS CAREFULLY! FOR
% NOW, PERTURBING AROUND SOL FOR s1=s3=0.
uss = 0.51746; vss = 0.327434; wss = 0.257944;

%Testing stability:
su=0.0250; sw =46.532699999999998; c=0.01;
uss = 0.248075252601345; vss = 0.855808334310556; wss = 0.961482373461736;

uss = 0.180495854926109; vss = 3.061963746749418e-04; wss = 0.858744469117753;

% Create a normally distributed vector of size 3*N, mean 1, stdev 1e-2.
rng(1);
rand_vec = abs(1+1e-4*randn(3*N,1)); 

% Initial conditions for ODE system
uvH_init = [uss*rand_vec(uN);vss*rand_vec(vN); wss*rand_vec(wN)];
