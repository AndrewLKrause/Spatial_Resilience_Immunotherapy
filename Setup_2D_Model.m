%THIS ONLY WORKS IF N=M^2 FOR INTEGER M!
M = sqrt(N);

% Use the code below form square grid Laplacian with arbitrary BCs:
dx = L/(M-1); % Spacing between grid points
[~,~,Lap] = laplacian([M,M],{'NN','NN'}); 
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