%%updated nondimensional turing code

%diffusion parameters
delta_u = 100;
delta_w = 1;

%Parameters
alp=0.01;
muu=0.167; rhou=0.692;
gamv=0.1;
rhow=2.5 ;gamw=0.001;muw=55.56;

%Number of grid points in each direction of s1-s3:
N=400;
% 
%Ranges of s1, s3
suR = linspace(0,0.06,N);
swR = linspace(0,18,N);
% 
% %reduced range for alp=0.01
% suR = linspace(0.002,0.005,N);
% swR = linspace(13,14,N);
% 
% %reduced range for alp=0.07
% suR = linspace(0.013,0.017,N);
% swR = linspace(0,3,N);


%M will store our stability values - coe is coexistence, cf is cancer-free,
%and coenum is a matrix to store number of feasible states of coexistence.
Mcoe = zeros(N);
Mcf = ones(N);
Mcoenum = zeros(N);
Turing = zeros(N);
Avalues = zeros(N);

for i=1:N
    for j=1:N
        su = suR(i);
        sw = swR(j);

        % Kinetic functions
        f = @(u,v,w)alp.*v-muu.*u+rhou.*u.*w./(1+w)+su;
        g = @(u,v,w)v.*(1-v)-u.*v./(gamv+v);
        h = @(u,v,w)rhow*u.*v./(gamw+v)-muw*w+sw;

        % Jacobian
        J = @(u,v,w)[-muu+rhou.*w./(1+w),     alp,     rhou.*u./(1 + w) - rhou*u.*w./(1 + w).^2;...
            -v./(gamv+v),      1-2.*v - u./(gamv+v) + v.*u./((gamv+v).^2),   0;...
            rhow.*v./(gamw+v), rhow.*u./(gamw + v)-rhow.*v.*u./((gamw+v).^2), -muw];
        Jcf = J(su*(muw+sw)/(muu*(muw+sw)-rhou*sw),0,sw/muw);
        %Compute coefficients for v equation
        c5 = (-muu + rhou)*rhow;
        c4 =-((alp + 2*(-1 + gamv)*(muu - rhou))*rhow);
        c3 = rhow*(alp - alp*gamv + rhou + (-4 + gamv)*gamv*rhou - su) - rhou*sw + muu*(muw - (1 + (-4 + gamv)*gamv)*rhow + sw);
        c2 = -((-1 + gamv)*rhow*(2*gamv*rhou + su)) - (-1 + gamv + gamw)*rhou*sw + alp*(muw + gamv*rhow + sw) + muu*((-1 + gamv + gamw)*muw + 2*(-1 + gamv)*gamv*rhow + (-1 + gamv + gamw)*sw);
        c1 = -gamw*muu*muw + gamv^2*(-muu + rhou)*rhow + muw*su + gamv*rhow*su - gamw*muu*sw + gamw*rhou*sw + su*sw + alp*gamw*(muw + sw) + gamv*(-1 + gamw)*(-rhou*sw + muu*(muw + sw));
        c0= gamw*(gamv*rhou*sw - gamv*muu*(muw + sw) + su*(muw + sw));

       
        %Compute roots of quintic
        v_sols = roots([c5,c4,c3,c2,c1,c0]);
        %Only consider real roots
        v_sols(imag(v_sols) ~= 0) = [];
        %Only consider positive roots smaller than 1/b.
        v_sols(0 > v_sols) = []; v_sols(v_sols > 1+1e-12) = [];
        Mcoenum(i,j) = length(v_sols);


        %Check cancer-free stability
        if((muu*(muw+sw)-rhou*sw)<0)
            Mcf(i,j) = 0;
            %display(["Warning: feasibility error at s_u=",...,
                %num2str(su), ", s_w = ", num2str(sw)]);
        elseif(max(real(eigs(Jcf)))>0)
            Mcf(i,j) = 0;
        else

        end

        %If no roots are found, set M_ij = 0.
        if(isempty(v_sols))
            Mcoe(i,j) = 0; %Interpret this as an infeasible/nonexistence state.
        else
            %Loop over each root found.
            for vss = v_sols'
                %Compute the values of u and w at steady state.
                uss = (1-vss)*(gamv + vss);
                wss = (rhow*uss*vss + gamw*sw + sw*vss)/((gamw + vss)*muw);
                %Check if there is an error in feasibility analysis.

                %Throw out spurious roots?
                if(abs(f(uss,vss,wss)) > 1e-12 || abs(g(uss,vss,wss)) > 1e-12 || abs(h(uss,vss,wss)) > 1e-12 )

                elseif(uss < -1e-14 || wss < -1e-14)
                    display(["Warning: feasibility error at s_u=",...,
                        num2str(su), ", s_w = ", num2str(sw)]);
                    %Mcoe(i,j) = -1;
                else
                    %Compute if any instability occurs - if so, leave value
                    %of M_Ij unchanged
                    if(max(real(eigs(J(uss,vss,wss))))>0)
                        Mcoe(i,j) = Mcoe(i,j) + 0;
                        %If no instability occurs, then this steady state is
                        %stable so add 1 to M_ij
                    else
                        Mcoe(i,j)=Mcoe(i,j) + 1;
                        %in region where coex steady state is stable,
                        %compute condition for Turing instability (note
                        %depends on ss but not explicitly on treatment
                        %parameters
                        J11 = -muu+rhou.*wss./(1+wss);
                        J12 = alp;
                        J13 = rhou.*uss./(1 + wss) - rhou*uss.*wss./(1 + wss).^2;
                        J21 = -vss./(gamv+vss);
                        J22 =  1-2.*vss - uss./(gamv+vss) + vss.*uss./((gamv+vss).^2);
                        J23 = 0;
                        J31 = rhow.*vss./(gamw+vss);
                        J32 =  rhow.*uss./(gamw + vss)-rhow.*vss.*uss./((gamw+vss).^2);
                        J33 = - muw;
                        A = delta_u*delta_w;
                        B = -(delta_u*J33 + delta_u*delta_w*J22 + delta_w*J11);
                        C = delta_u*J22*J33 + J11*J33 + delta_w*J11*J22 - delta_u*J23*J32 - J13*J31 - delta_w*J12*J21;
                        D = - J11*J22*J33 + J11*J23*J32 + J12*J21*J33 - J12*J23*J31 - J13*J21*J32 + J13*J22*J31; 
                        %check to see if it meets Turing instability
                        %criteron
                        mustar = (-2*B + sqrt(4*B^2 - 12*A*C))/(6*A);

                        a0 = A*mustar^3 + B*mustar^2 + C*mustar + D;

                        Avalues(i,j) =a0;
                        if(a0 < -1e-16 && mustar>0)
                            Turing(i,j) = 2; % if hom stables and turing unstable set value of turing matrix to 2
                        elseif(Turing(i,j)==0)
                            Turing(i,j) = 1; %if hom stable and not turing unstable set to 1

                            




                        end
                    end
                end
            end
        end
    end
end

    


close all;

%separate code to plot turing figure
figure("position",[10,10,900,600]);
imagesc(flipud(Turing'))
ax = gca; set(ax,'Fontsize',25);
xlabel('$\sigma_u$','interpreter','latex', FontSize=40)
ylabel('$\sigma_w$','interpreter','latex',FontSize=40)
ax.XTick = [1, ax.XTick];
ax.XTickLabel = round(suR(ax.XTick),3);
ax.YTick = [1, ax.YTick];
ax.YTickLabel = round(flip(swR((ax.YTick))),2);
ax.TickLabelInterpreter = 'latex';
ccmap = [
    0.86, 0.15, 0.50;   
    1.00, 0.69, 0.00;   
    0.39, 0.56, 1.00;   
];
string =["no stable coex","stable", "Turing unstable"];
colormap(flipud(ccmap))
hcb = colorbar('Ticks',linspace(0.5,4-0.5,4),...
        'TickLabels',string(1:3));
clim([0 3]);
hcb.TickLabelInterpreter = 'latex'; hcb.Visible='off';

