

uvH_init = [U(end,uN)'; U(end,vN)'; U(end,wN)'];


f = @(u,v,w)c*v-muu*u+pu*u.*w./(gu+w)+su;
g = @(u,v,w)rv*v.*(1-b*v)-pv*u.*v./(gv+v);
h = @(u,v,w)pw*u.*v./(gw+v)-muw*w+sw;

% The right-hand-side of our discretized ODE system
FH_PDE = @(t, U)[f(U(uN),U(vN),U(wN))+d1*Lap*U(uN);...,
    g(U(uN),U(vN),U(wN))+d2*Lap*U(vN);h(U(uN),U(vN),U(wN))+d3*Lap*U(wN)];

% Solve the PDE - optional could add: 'reltol',1e-9,'AbsTol',1e-9,
opts = odeset('JPattern',JPattern,'reltol',1e-9,'AbsTol',1e-9);
[~, U] = ode15s(FH_PDE,tspan,uvH_init,opts);

close all;
plot(x,U(end,uN),'linewidth',2); hold on;
plot(x,U(end,vN),'linewidth',2);
plot(x,U(end,wN),'linewidth',2);

legend('$u$','$v$','$w$','interpreter','latex')

figure

imagesc(flipud(U(:,vN)))
colorbar