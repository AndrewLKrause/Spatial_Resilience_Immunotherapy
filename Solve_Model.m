clear;

dims=1;

m = 10000;
Init_Parameters;



% Form the (Sparse) Laplacian matrix
e=ones(m,1);
Lap=spdiags([e,-2*e,e],[1,0,-1],m,m);

%Neumann BCs
Lap(1,1)=-1;
Lap(end,end)=-1;

if(dims==1)
    %1D Laplacian
    Lap = (1/dx)^2*Lap;
elseif(dims==2)
    %2d Laplacian
    I = speye(m);
    Lap = (1/dx)^2*(kron(Lap,I) + kron(I, Lap));
end

% Set up the Jacobian sparsity pattern - important for speed!
%The second method is much more memory-efficient.
%Id = eye(N);
%JPattern = sparse([Lap, Id, Id; Id, Lap, Id; Id, Id, Lap]);
JPattern = blkdiag(Lap,Lap,Lap)+spdiags(1,[-N,N],3*N,3*N)+spdiags(1,[-2*N,2*N],3*N,3*N);
% The right-hand-side of our discretized ODE system
FH_PDE = @(t, U)[f(U(uN),U(vN),U(wN))+K_u(t)+d1*Lap*U(uN);...,
    g(U(uN),U(vN),U(wN))+d2*Lap*U(vN);h(U(uN),U(vN),U(wN))+K_w(t)+d3*Lap*U(wN)];

% Solve the PDE - optional could add: 'reltol',1e-9,'AbsTol',1e-9,
opts = odeset('JPattern',JPattern,'reltol',1e-11,'AbsTol',1e-11,'MaxStep',T/2000);
[~, U] = ode15s(FH_PDE,tspan,uvH_init,opts);

% Smooth the initial transient; ideally this can be removed with better
% initial conditions
%U(1:100,:) = repmat(U(101,:),[100,1]);

rows=size(U,1);
for t=1:rows
    tv=sum(U(t,vN))/N;
    tu=sum(U(t,uN))/N;
    tw=sum(U(t,wN))/N;
    ut(t)=tu;
    vt(t)=tv;
    wt(t)=tw;
end



close all;

%PlotSol




f = figure;
%f.Position(3:4) = f.Position(3:4)*1.3; f.Position(1:2) = f.Position(1:2)*0.7;

imagesc(U(:,vN)')
colour_lims = clim;
set(gca,'YDir','normal')
ylabel('$x$','interpreter','latex')
xlabel('$t$','interpreter','latex')
colorbar;
ax = gca;
ax.YTick = [1,round(m/3), round(2*m/3), m];
ax.YTickLabel = {'0', num2str(L/3),num2str(2*L/3),num2str(L)};
ax.XTick = [1,rows/2,rows];
ax.XTickLabel = {'0', num2str(round(T/2)),num2str(T)};
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

ax = gca;
ax.XTick = [1,rows/2,rows];
ax.XTickLabel = {'0', num2str(round(T/2)),num2str(T)};
