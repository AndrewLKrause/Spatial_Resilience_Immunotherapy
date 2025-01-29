function Gu=sGjac(p,u)
    par = u(p.nu + 1:end); % identify parameters
    n = p.np;
    [f1u, f1v, f1w, f2u, f2v, f2w, f3u, f3v, f3w] = njac(p, u, par); % loading the non-linear jacobian in nodal form, see below
    Fu = [[spdiags(f1u, 0, n, n), spdiags(f1v, 0, n, n), spdiags(f1w, 0, n, n)];
          [spdiags(f2u, 0, n, n), spdiags(f2v, 0, n, n), spdiags(f2w, 0, n, n)];
          [spdiags(f3u, 0, n, n), spdiags(f3v, 0, n, n), spdiags(f3w, 0, n, n)]]; % nonlinear part of the jacobian
    Gu = kron([[par(10), 0, 0]; [0, 1.0, 0]; [0, 0, par(11)]], p.mat.K) - p.mat.M*Fu; % assemble the jacobian
end

function [f1u, f1v, f1w, f2u, f2v, f2w, f3u, f3v, f3w] = njac(p, u, par) % Jacobian for Schnakenberg, nodal version
    u1 = u(1:p.np); % identify solution component 1
    u2 = u(p.np + 1:2*p.np); % identify solution component 2
    u3 = u(2*p.np + 1:3*p.np); % identify solution component 2
    % entries of the jacobian
    f1u = - par(2) + par(3)*u3./(1 + u3);
    f1v = par(1)*ones(size(u2));
    f1w = par(3)*u1./(1 + u3).^2;
    f2u = - u2./(par(5) + u2);
    f2v = 1 - 2*u2 - par(5)*u1./(par(5) + u2).^2;
    f2w = 0.0;
    f3u = par(6)*u2./(par(7) + u2);
    f3v = par(6)*par(7)*u1./(par(7) + u2).^2;
    f3w = - par(8)*ones(size(u3));
end