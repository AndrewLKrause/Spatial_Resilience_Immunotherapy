%This code will produce two matrices, M and Mn where:
%Mn(i,j) will be the number of feasible coexistence steady states which
%exist, and M(i,j) will be the number of feasible and stable coexistence
%steady states. These correspond to the parameters below with s1 = s1R(i)
%and s3 = s3R(j). TODO: UPDATED NOTATION


%Parameters
c = 0.1;
muu=0.167; pu=0.69167; gu=20;
rv=1;b=1;pv=0.5555556;gv=0.1;
pw=27.778;gw=0.001;muw=55.55556;

%Number of grid points in each direction of s1-s3:
N=400;

%Ranges of s1, s3
suR = linspace(0,0.05,N);
swR = linspace(0,100,N);

%M will store our stability values - coe is coexistence, cf is cancer-free,
%and coenum is a matrix to store number of feasible states of coexistence.
Mcoe = zeros(N);
Mcf = ones(N);
Mcoenum = zeros(N);

for i=1:N
    for j=1:N
        su = suR(i);
        sw = swR(j);

        % Kinetic functions
        f = @(u,v,w)c*v-muu*u+pu*u.*w./(gu+w)+su;
        g = @(u,v,w)rv*v.*(1-b*v)-pv*u.*v./(gv+v);
        h = @(u,v,w)pw*u.*v./(gw+v)-muw*w+sw;

        % Jacobian
        J = @(u,v,w)[-muu+pu.*w./(gu+w),     c,     pu*u.*gu./(gu + w).^2;...
            -pv.*v./(gv+v),      rv.*(1-2.*b.*v)-pv*gv*u./((gv+v).^2),   0;...
            pw*v./(gw+v),pw*gw*u./((gw+v).^2),-muw];
        Jcf = J(su*(gu*muw+sw)/(muu*(gu*muw+sw)-pu*sw),0,sw/muw);
        %Compute coefficients for v equation
        c5 = -b^2*rv^2*pw*(-pu + muu);
        c4 = -2*pw*((-pu + muu)*(b*gv - 1)*rv + pv*c/2)*b*rv;
        c3 = -rv*(pw*(b^2*gv^2 - 4*b*gv + 1)*(-pu + muu)*rv + (b*c*pw*gv + (su*pw + (-gu*muw - sw)*muu + pu*sw)*b - c*pw)*pv);
        c2 = 2*gv*pw*(-pu + muu)*(b*gv - 1)*rv^2 - (((su*pw + (-gu*muw - sw)*muu + pu*sw)*b - c*pw)*gv - gw*((gu*muw + sw)*muu - pu*sw)*b - su*pw + (gu*muw + sw)*muu - pu*sw)*pv*rv + pv^2*c*(gu*muw + sw);
        c1 = -gv^2*pw*(-pu + muu)*rv^2 + ((gw*((gu*muw + sw)*muu - pu*sw)*b + su*pw + (-gu*muw - sw)*muu + pu*sw)*gv - gw*((gu*muw + sw)*muu - pu*sw))*pv*rv + pv^2*(gu*muw + sw)*(c*gw + su);
        c0 = -gw*(gv*((gu*muw + sw)*muu - pu*sw)*rv - pv*su*(gu*muw + sw))*pv;

        %c5 = -b^2*rv^2*pw*(-pu + muu)/pv;
        %c4 = -2*rv*pw*((-pu + muu)*(b*gv - 1)*rv + c*pv/2)*b/pv;
        %c3 = -rv*(pw*(b^2*gv^2 - 4*b*gv + 1)*(-pu + muu)*rv + pv*(b*c*pw*gv + (su*pw + (-gu*muw - sw)*muu + sw*pu)*b - c*pw))/pv;
        %c2 = (2*gv*pw*(-pu + muu)*(b*gv - 1)*rv^2 - (((su*pw + (-gu*muw - sw)*muu + sw*pu)*b - c*pw)*gv - ((gu*muw + sw)*muu - sw*pu)*gw*b - su*pw + (gu*muw + sw)*muu - sw*pu)*pv*rv + c*pv^2*(gu*muw + sw))/pv;
        %c1 = (-gv^2*pw*(-pu + muu)*rv^2 + ((((gu*muw + sw)*muu - sw*pu)*gw*b + su*pw + (-gu*muw - sw)*muu + sw*pu)*gv - ((gu*muw + sw)*muu - sw*pu)*gw)*pv*rv + pv^2*(gu*muw + sw)*(c*gw + su))/pv;
        %c0 = -gw*(((gu*muw + sw)*muu - sw*pu)*gv*rv - su*pv*(gu*muw + sw));

        %Compute roots of quintic
        v_sols = roots([c5,c4,c3,c2,c1,c0]);
        %Only consider real roots
        v_sols(imag(v_sols) ~= 0) = [];
        %Only consider positive roots smaller than 1/b.
        v_sols(0 > v_sols) = []; v_sols(v_sols > 1/b+1e-12) = [];
        Mcoenum(i,j) = length(v_sols);


        %Check cancer-free stability
        if((muu*(gu*muw+sw)-pu*sw)<0)
            %Mcf(i,j) = -1;
            display(["Warning: feasibility error at s_u=",...,
                num2str(su), ", s_w = ", num2str(sw)]);
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
                uss = -rv*(b*vss - 1)*(gv + vss)/pv;
                %wss = (pw*uss*vss + gw*sw + sw*vss)/((gw + vss)*muw);
                wss = (pw*uss*vss + gw*sw + sw*vss)/((gw + vss)*muw);
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
                    end
                end
            end
        end
    end
end

%Bad plotting code
close all;
% figure;
%
% %Number of feasible steady states - x axis is s1R, y axis is s3R
% imagesc(flipud(Mn'))
% colorbar;
% xlabel('$s_u$','interpreter','latex')
% ylabel('$s_u$','interpreter','latex')
% ax = gca; set(ax,'Fontsize',16);
% ax.XTick = [1, ax.XTick];
% ax.XTickLabel = suR(ax.XTick);
% ax.YTick = [1, ax.YTick];
% ax.YTickLabel = flip(swR((ax.YTick)));
% title('# of feasible steady states')
%
%
% figure
%
% %Number of feasible+stable coexistence steady states - x axis is s1R, y axis is s3R
% imagesc(flipud(M'))
% colorbar;
% xlabel('$s_u$','interpreter','latex')
% ylabel('$s_u$','interpreter','latex')
% ax = gca; set(ax,'Fontsize',16);
% ax.XTick = [1, ax.XTick];
% ax.XTickLabel = suR(ax.XTick);
% ax.YTick = [1, ax.YTick];
% ax.YTickLabel = flip(swR((ax.YTick)));
% title('# of feasible+stable COEXISTENCE steady states')
%
% figure
%
% %Number of feasible+stable cancer free steady states - x axis is s1R, y axis is s3R
% imagesc(flipud(Mcf'))
% colorbar;
% xlabel('$s_u$','interpreter','latex')
% ylabel('$s_w$','interpreter','latex')
% ax = gca; set(ax,'Fontsize',16);
% ax.XTick = [1, ax.XTick];
% ax.XTickLabel = suR(ax.XTick);
% ax.YTick = [1, ax.YTick];
% ax.YTickLabel = flip(swR((ax.YTick)));
% title('# of feasible+stable CANCER-FREE steady states')


%IMPORTANT NOTE: This code will yell if it detects multistability in
%coexistence state. It it does so, then plot colours shown are WRONG!
if(max(max(Mcoe))>1)
    'ERROR! Multistability of coexistence equilibria detected!';
end
figure

%Number of feasible+stable cancer free steady states - x axis is s1R, y axis is s3R
imagesc(flipud(Mcoe'+4*Mcf'))
%colorbar;
xlabel('$s_u$','interpreter','latex')
ylabel('$s_w$','interpreter','latex')
ax = gca; set(ax,'Fontsize',14);
ax.XTick = [1, ax.XTick];
ax.XTickLabel = suR(ax.XTick);
ax.YTick = [1, ax.YTick];
ax.YTickLabel = flip(swR((ax.YTick)));
title(['Stable Equilibria for $c=$',num2str(c)],'interpreter','latex')

mx = max(max(Mcoe+4*Mcf));
mn = min(min(Mcoe+4*Mcf));
C=lines(8);
strings = {"None Stable", "1 Coex","2 Coex","3 Coex", "CF", "CF+1Coex",...
    "CF+2Coex","CF+3Coex"};
if(mn <0)
    "Feasibility errors - no plots produced."

else
    
    colormap(flipud(C(mn+1:mx+1,:)))
    colorbar('Ticks',linspace(mn+0.5,mx-0.5,mx-mn+1.5),...
        'TickLabels',strings(1+mn:mx+2));

end


