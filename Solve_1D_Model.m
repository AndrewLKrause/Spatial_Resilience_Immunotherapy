parameters;
Setup_1D_Model;


% The right-hand-side of our discretized ODE system
FH_PDE = @(t, U)[f(U(uN),U(vN),U(wN))+d1*Lap*U(uN);...,
    g(U(uN),U(vN),U(wN))+d2*Lap*U(vN);h(U(uN),U(vN),U(wN))+d3*Lap*U(wN)];

% Solve the PDE - optional could add: 'reltol',1e-9,'AbsTol',1e-9,
opts = odeset('JPattern',JPattern);
[~, U] = ode15s(FH_PDE,tspan,uvH_init,opts);

close all;
plot(x,U(end,uN),'linewidth',2); hold on;
plot(x,U(end,vN),'linewidth',2);
plot(x,U(end,wN),'linewidth',2);