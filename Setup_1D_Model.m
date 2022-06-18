

dx = L/(N-1); % Spacing between grid points

% Form the Laplacian
%e = ones(N,1); % Vector of ones to use across the diagonals
%Lap= spdiags([e -2*e e], -1:1, N, N); % Diagonal Laplacian
%Lap(1,1) = -1; Lap(N,N) = -1; % Neumann boundary conditions
%Lap = Lap./dx^2; % Scale the finite-difference operator

% The above is a way of doing this directly for Neumann BCs.
% Instead, let's use the code below to do it with arbitrary BCs:
[~,~,Lap] = laplacian(N,{'NN'}); 
Lap = -Lap./(dx)^2;
% The N is number of gridpoints and NN is Neumann at both boundaries.
% The -is a convention due to treating the Laplacian as a +ive operator.

% Index labels
uN = 1:N; vN = N+1:2*N; wN = 2*N+1:3*N;

% Set up the Jacobian sparsity pattern - important for speed!
Id = eye(N);
JPattern = sparse([Lap, Id, Id; Id, Lap, Id; Id, Id, Lap]);

% Create a normally distributed vector of size 3*N, mean 1, stdev 1e-2.
rng(1);
rand_vec = abs(1+1e-2*randn(3*N,1)); 

% Initial conditions for ODE system
uvH_init = [uss*rand_vec(uN);vss*rand_vec(vN); wss*rand_vec(wN)];