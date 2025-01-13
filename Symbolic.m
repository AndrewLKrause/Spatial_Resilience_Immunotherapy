%This file computes the expressions given in Section 3.1.

clear;
%NB: The commands with 'latex' in them can be uncommented out to produce
%the LaTeX code for the expression.

%Define the symbolic variables and reaction kinetics.
syms alpha rho_u rho_v gamma_w mu_u sigma_u gamma_v rho_w gamma_w mu_w sigma_w f g h u v w lambda;
f = alpha*v -mu_u*u + rho_u*u*w/(gamma_w + w) + sigma_u;
g = v*(1-v) -rho_v*u*v/(gamma_v + v);
h = rho_w*u*v/(gamma_w + v) - mu_w*w  +sigma_w;

%Cancer-free equilibrium (v_0=0)
[u_0, w_0] = solve(subs(f,v,0)==0, subs(g,v,0)==0, subs(h,v,0)==0);
latex(u_0);
latex(w_0);

%Jacobian
J = jacobian([f,g,h],[u,v,w]); 
latex(J);

%Stability of cancer-free state via eigenvalues.
J0 = subs(subs(subs(J,u,u_0),w,w_0),v,0);
latex(eig(J0))

%Cancer coexistence state (NB: Need to clear symbols u_0 and w_0).
clear('w_0'); clear('u_0');
u_0 = solve(g/v==0,u);
w_0 = solve(h==0,w);
latex(w_0);
w_0 = subs(w_0,u,u_0);
latex(simplify(w_0));

[N, D] = numden(collect(subs(subs(f,u,u_0),w,w_0),v));

%Now v_0 satisfies this 5th order polynomial C(1)+C(2)v_0+...+C(6)v_0^5=0.
C = coeffs(collect(N,v),v);
latex(simplify(C));