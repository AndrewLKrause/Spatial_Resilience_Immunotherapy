Init_Parameters;
uN = 1; vN = 2; wN = 3;
% The right-hand-side of our discretized ODE system
FH_ODE = @(t, U)[f(U(uN),U(vN),U(wN));...,
    g(U(uN),U(vN),U(wN));h(U(uN),U(vN),U(wN))];

% Solve the ODE - optional could add: 'reltol',1e-9,'AbsTol',1e-9,
opts = odeset('reltol',1e-12,'AbsTol',1e-12);
uvH_init = [su*(gu*muw+sw)/(muu*(gu*muw+sw)-pu*sw),1e-6,sw/muw]';
[T, U] = ode15s(FH_ODE,tspan,uvH_init,opts);

close all;
plot(T,U(:,uN),'linewidth',2); hold on;
plot(T,U(:,vN),'linewidth',2);
plot(T,U(:,wN),'linewidth',2);

legend('$u$','$v$','$w$','interpreter','latex')