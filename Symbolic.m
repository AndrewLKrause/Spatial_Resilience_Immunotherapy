%This file computes the expressions given in Section 3.1.

clear;
%NB: The commands with 'latex' in them can be uncommented out to produce
%the LaTeX code for the expression.

%Define the symbolic variables and reaction kinetics.
syms alpha rho_u rho_v gamma_w mu_u sigma_u gamma_v rho_w gamma_w mu_w sigma_w f g h u v w lambda delta_u delta_w k;
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
latex(eig(J0));

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

%Useful for exploring the linearization of the Turing analysis; commented
%out for now as not used.

latex(simplify(C));

%M = J-k^2*diag([delta_u,1,delta_w]);

%C = simplify(charpoly(subs(subs(M,u,u_0),w,w_0)));
%C2 = simplify(charpoly(M));

%Turing analysis of cancer-free state
J = J0;
A = delta_u*delta_w;
B = -(delta_u*J(3,3)+delta_u*delta_w*J(2,2)+delta_w*J(1,1));
C =  (delta_u*J(2,2)*J(3,3)+J(1,1)*J(3,3)+delta_w*J(1,1)*J(2,2)-delta_u*J(2,3)*J(3,2)-J(1,3)*J(3,1)-delta_w*J(1,2)*J(2,1));
D = J(1,1)*J(2,2)*J(3,3)+J(1,1)*J(2,3)*J(3,2)+J(1,2)*J(2,1)*J(3,3)-J(1,2)*J(2,3)*J(3,1)-J(1,3)*J(2,1)*J(3,2)+J(1,3)*J(2,2)*J(3,1);

simplify(B^2-4*A*C)