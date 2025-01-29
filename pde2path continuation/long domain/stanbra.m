function out=stanbra(p,u)
% STANBRA: standard function for output beyond bradat(p) to p.branch.
% here simply: [pars, max(abs(u1)), min(abs(u1))
% highly problem dependent, take this as a template
%
upde=u(1:p.nu); np=p.nu/p.nc.neq; 
out=[u(p.nu+1:end); max(abs(upde(np + 1:2*np))); min(abs(upde(np + 1:2*np)))]; 