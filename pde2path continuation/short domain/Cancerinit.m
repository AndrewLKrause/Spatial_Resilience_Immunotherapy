function p=Cancerinit(p, h, par)
    p = stanparam(p);
    screenlayout(p);
    dir = 'p0';
    p = setfn(p, dir);
    p.nc.lammin = 0.0;
    p.nc.lammax = 0.035;
    p.nc.bisecmax = 50;
    p.nc.neq = 3;
    p.sw.sfem = - 1;
    p.fuha.sG = @sG;
    p.fuha.sGjac = @sGjac; 
    lx = 15.0; % wavenumer of the critical mode
    p.pdeo = stanpdeo1D(lx, h);
    p.np = p.pdeo.grid.nPoints;
    p.nu = p.np*p.nc.neq; 
    p = setfemops(p);
    p.nc.ilam = 4;
    p.sol.xi = 1/p.nu;
    p.sol.ds = 1e-4;
    p.nc.dsmax = 5e-3;
    if strcmp(dir, 'p0')
        v = zeros(p.np, 1);
        u = par(4)*(par(8) + par(9))/(- par(3)*par(9) + par(2)*(par(8) + par(9)))*ones(p.np, 1);
        w = par(9)/par(8)*ones(p.np, 1);
    else
        v = 0.6058511193427998*ones(p.np, 1); % hom.soln 
        u = - (v - 1).*(v + par(5));
        w = - 1 + u*par(3)./(v*par(1) - u*par(2) + u*par(3) + par(4));
    end
    % v = 0.6058511193427998*ones(p.np, 1); % hom.soln 
    % u = - (v - 1).*(v + par(5));
    % w = - 1 + u*par(3)./(v*par(1) - u*par(2) + u*par(3) + par(4));
    p.u = [u; v; w; par'];
    p.plot.pmod = 0;
    p.file.smod = 5;
    p.plot.pcmp = 2;
end