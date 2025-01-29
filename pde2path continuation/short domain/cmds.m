%% 1 - initialising the problem
changedir = 1;

if changedir==1
    sleep()
end

    clc
    close all;
    keep pphome;
    p = [];
    % alpha, mu_u, rho_u, sigma_u, gamma_v, rho_w, gamma_w, mu_w, sigma_w,
    % delta_u, delta_w
    par = [0.07, 0.167, 0.692, 0.0, 0.1, 2.5, 0.001, 55.56, 0.5, 100.0, 100.0]; 
    h = 1e-2;
    p = Cancerinit(p, h, par); % also use nper=80, 120, 160, 200

%% 2 - continue trivial branch to find BP 
    tic;
    p.sw.foldcheck = 1;
    p.sw.bifcheck = 2;
    p = cont(p, 135);
    toc

%% 3 - switch to periodic branch and continue. For comparison of \ and 
    % lssbel, switch off stuff not related to lss 
    p = swibra('p', 'bpt1', 'b1', 0.02);
    p.nc.ds = 1e-2;
    p.nc.dsmax = 1e-2;
    p.sw.spcalc = 1;
    p.sw.foldcheck = 1;
    p.sw.bifcheck = 1;
    p.sw.verb = 2;
    % p.nc.neig = 2;
    t1 = tic;
    p = cont(p, 58);
    t1 = toc(t1); % cont with default settings

%% plots
    
    lw = 4;
    lwst = 4;
    lwun = 2;
    
    ps = 100;
    fs = 60;
    figure(3)
    plotbra('p', 'pt136', 3, 0, 'tyun', ':k', 'tyst', '-k', 'ms', 0, 'fms', 0, 'lwst', lwst, 'lwun', lwun);
    plotbra('p0', 'pt55', 3, 0, 'tyun', ':k', 'tyst', '-k', 'ms', 0, 'fms', 0, 'lwst', lwst, 'lwun', lwun);
    plotbra('b1', 'pt58', 3, 0, 'tyun', ':r', 'tyst', '-r', 'ms', 0, 'fms', 0, 'lwst', lwst, 'lwun', lwun);
    xlabel('$\sigma_u$', 'interpreter', 'latex')
    ylabel('$||v||_{L^1}$', 'interpreter', 'latex')
    set(gca, 'fontsize', fs)