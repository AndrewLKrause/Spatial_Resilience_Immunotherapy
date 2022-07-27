
N = 100^2;
Init_Parameters;

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


% The right-hand-side of our discretized ODE system
FH_PDE = @(t, U)[f(U(uN),U(vN),U(wN))+d1*Lap*U(uN);...,
    g(U(uN),U(vN),U(wN))+d2*Lap*U(vN);h(U(uN),U(vN),U(wN))+d3*Lap*U(wN)];

% Solve the PDE - optional could add: 'reltol',1e-9,'AbsTol',1e-9,
opts = odeset('JPattern',JPattern);
[~, U] = ode15s(FH_PDE,tspan,uvH_init,opts);

close all;
figure
imagesc(reshape(U(end,uN),M,M));colorbar;title('$u$','interpreter','latex')
figure
imagesc(reshape(U(end,vN),M,M));colorbar;title('$v$','interpreter','latex')
figure
imagesc(reshape(U(end,wN),M,M));colorbar;title('$w$','interpreter','latex')