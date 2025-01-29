function r = sG(p,u)
    u1 = u(1:p.np); % identify solution component 1
    u2 = u(p.np + 1:2*p.np); % identify solution component 2
    u3 = u(2*p.np + 1:3*p.np); % identify solution component 3
    par = u(p.nu + 1:end); % identify parameters
    f1 = par(1)*u2 - par(2)*u1 + par(3)*u1.*u3./(1 + u3) + par(4); % non-linearity for u1 
    f2 = u2.*(1 - u2) - u1.*u2./(par(5) + u2); % non-linearity for u2
    f3 = par(6)*u1.*u2./(par(7) + u2) - par(8)*u3 + par(9); % non-linearity for u3
    f = [f1; f2; f3];
    r = kron([[par(10), 0, 0]; [0, 1.0, 0]; [0, 0, par(11)]], p.mat.K)*u(1:p.nu) - p.mat.M*f; % calculation of the residual
end


