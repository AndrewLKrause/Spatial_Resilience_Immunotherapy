
N = 3000;
Init_Parameters;

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

% Set up the Jacobian sparsity pattern - important for speed!
Id = eye(N);
JPattern = sparse([Lap, Id, Id; Id, Lap, Id; Id, Id, Lap]);

% The right-hand-side of our discretized ODE system
FH_PDE = @(t, U)[f(U(uN),U(vN),U(wN))+K(t)+d1*Lap*U(uN);...,
    g(U(uN),U(vN),U(wN))+d2*Lap*U(vN);h(U(uN),U(vN),U(wN))+d3*Lap*U(wN)];

% Solve the PDE - optional could add: 'reltol',1e-9,'AbsTol',1e-9,
opts = odeset('JPattern',JPattern,'reltol',1e-11,'AbsTol',1e-11);
[~, U] = ode15s(FH_PDE,tspan,uvH_init,opts);

% Smooth the initial transient; ideally this can be removed with better
% initial conditions
%U(1:100,:) = repmat(U(101,:),[100,1]);

rows=size(U,1);
for t=1:rows
    tv=sum(U(t,vN))/N;
    tu=sum(U(t,uN))/N;
    tw=sum(U(t,wN))/N;
    vt(t)=tv;
    ut(t)=tu;
    wt(t)=tw;
end



close all;
% plot(x,U(end,uN),'linewidth',2); hold on;
% plot(x,U(end,vN),'linewidth',2);
% plot(x,U(end,wN),'--','linewidth',2);
% 
% legend('$u$','$v$','$w$','interpreter','latex')

f = figure;
%f.Position(3:4) = f.Position(3:4)*1.3; f.Position(1:2) = f.Position(1:2)*0.7;

imagesc(U(:,vN)')
colour_lims = clim;
set(gca,'YDir','normal')
ylabel('$x$','interpreter','latex')
xlabel('$t$','interpreter','latex')
colorbar;
ax = gca;
ax.YTickLabel = {'100','200','300'}
set(gca,'fontsize',22);
colour_limits = clim;

% figure
% gca.Position = gca.Position*1.5;
% colorbar
% axis off;
% set(gca,'fontsize',20);
% clim(colour_limits)

f = figure;
%f.Position(3:4) = f.Position(3:4)*1.3; f.Position(1:2) = f.Position(1:2)*0.7;
plot(ut,'linewidth',2); hold on;
plot(vt,'linewidth',2);
plot(wt,'--','linewidth',2)
legend('$u$','$v$','$w$','interpreter','latex')
xlabel('$t$','interpreter','latex')
set(gca,'fontsize',22);
