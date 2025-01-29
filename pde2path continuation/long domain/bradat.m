function bra=bradat(p)
% BRADAT: set default non-user-defined branch data 
%
%  bra=bradat(p)
%
% bra=[p.file.count;                % 1
%      p.sol.ptype;                 % 2
%      p.sol.ineg;                  % 3
%      getlam(p);                   % 4
%      p.sol.err;                   % 5
%      L2-norm 1st component]       % 6 
ineg=p.sol.ineg; 
if isfield(p,'hopf'); ho=p.hopf; 
    try ineg=ho.ind; end;   
end
ineg=max(ineg); u=p.u; 
try; M=getM(p); n1=floor(p.nu/p.nc.neq); 
    % l2=sqrt(u(n1 + 1:2*n1)'*M(n1 + 1:2*n1,n1 + 1:2*n1)*u(n1 + 1:2*n1)); 
    [po, t, e] = getpte(p);
    l1 = triint(abs(u(n1 + 1:2*n1)), po, t);
catch; l2=0; end % catch an error in Xcont
bra=[p.file.count; p.sol.ptype; ineg'; getlam(p); p.sol.err; l1]; 