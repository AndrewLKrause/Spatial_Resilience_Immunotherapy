% Geometric Parameters
L = 1; % Length of macroscale domain

% Mesh size parameters
N = 100^2; % Total number of grid points to use


% Time/space scale for simulation
T = 1500; 
tspan = linspace(0,T,1e3); % Interpolate solution on [0,T] with 1e2 points
x = linspace(0,L,N);

% Diffusion parameters
%d1=1; d2=0.0001; d3=1;
d1=0.001; d2=0.0000199; d3=0.01;
 
% Kinetic parameters
c = 0.25; mu1=0.167; p1=0.69167; g1=20;s1=.01;
r2=1;b=1;p2=0.5555556;g2=0.1;
p3=27.778;g3=0.001;mu3=55.55556;s3=.01;

% Kinetic functions
f = @(u,v,w)c*v-mu1*u+p1*u.*w./(g1+w)+s1;
g = @(u,v,w)r2*v.*(1-b*v)-p2*u.*v./(g2+v);
h = @(u,v,w)p3*u.*v./(g3+v)-mu3*w+s3;

% Steady states - TODO NEED TO DISCUSS STEADY STATE ANALYSIS CAREFULLY! FOR
% NOW, PERTURBING AROUND SOL FOR s1=s3=0.
uss = 0.51746; vss = 0.327434; wss = 0.257944;
