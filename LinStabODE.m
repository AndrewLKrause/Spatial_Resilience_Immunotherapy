%This code will produce two matrices, M and Mn where:
%Mn(i,j) will be the number of feasible coexistence steady states which
%exist, and M(i,j) will be the number of feasible and stable coexistence
%steady states. These correspond to the parameters below with s1 = s1R(i)
%and s3 = s3R(j). 


%Parameters
c = 0.1;
mu1=0.167; p1=0.69167; g1=20;
r2=1;b=1;p2=0.5555556;g2=0.1;
p3=27.778;g3=0.001;mu3=55.55556;

%Number of grid points in each direction of s1-s3:
N=300;

%Ranges of s1, s3
s1R = linspace(0,0.05,N);
s3R = linspace(0,100,N);

%M will store our stability values
M = zeros(N);
Mn = zeros(N);

for i=1:N
    for j=1:N
        s1 = s1R(i);
        s3 = s3R(j);

        % Kinetic functions
        f = @(u,v,w)c*v-mu1*u+p1*u.*w./(g1+w)+s1;
        g = @(u,v,w)r2*v.*(1-b*v)-p2*u.*v./(g2+v);
        h = @(u,v,w)p3*u.*v./(g3+v)-mu3*w+s3;

        % Jacobian
        J = @(u,v,w)[-mu1+p1.*w./(g1+w),     c,     g1.*p1.*u./((g1+w).^2);...
            -p2.*v./(g2+v),      r2.*(1-2.*b.*v)+p2*g2*u./((g2+v).^2),   0;...
            p3*v./(g3+v),p3*g3*u./((g3+v).^2),-mu3];

        %Compute coefficients for v equation
        c5 = -b^2*r2^2*p3*(-p1 + mu1)/p2;
        c4 = -2*r2*p3*((-p1 + mu1)*(b*g2 - 1)*r2 + c*p2/2)*b/p2;
        c3 = -r2*(p3*(b^2*g2^2 - 4*b*g2 + 1)*(-p1 + mu1)*r2 + p2*(b*c*p3*g2 + (s1*p3 + (-g1*mu3 - s3)*mu1 + s3*p1)*b - c*p3))/p2;
        c2 = (2*g2*p3*(-p1 + mu1)*(b*g2 - 1)*r2^2 - (((s1*p3 + (-g1*mu3 - s3)*mu1 + s3*p1)*b - c*p3)*g2 - ((g1*mu3 + s3)*mu1 - s3*p1)*g3*b - s1*p3 + (g1*mu3 + s3)*mu1 - s3*p1)*p2*r2 + c*p2^2*(g1*mu3 + s3))/p2;
        c1 = (-g2^2*p3*(-p1 + mu1)*r2^2 + ((((g1*mu3 + s3)*mu1 - s3*p1)*g3*b + s1*p3 + (-g1*mu3 - s3)*mu1 + s3*p1)*g2 - ((g1*mu3 + s3)*mu1 - s3*p1)*g3)*p2*r2 + p2^2*(g1*mu3 + s3)*(c*g3 + s1))/p2;
        c0 = -g3*(((g1*mu3 + s3)*mu1 - s3*p1)*g2*r2 - s1*p2*(g1*mu3 + s3));

        %Compute roots of quintic
        v_sols = roots([c5,c4,c3,c2,c1,c0]);
        %Only consider real roots
        v_sols(imag(v_sols) ~= 0) = [];
        %Only consider positive roots smaller than 1/b.
        v_sols(0 > v_sols) = []; v_sols(v_sols > 1/b) = [];
        Mn(i,j) = length(v_sols);

        %If no roots are found, set M_ij = 0.
        if(isempty(v_sols))
            M(i,j) = 0; %Interpret this as an infeasible/nonexistence state.
        else
            %Loop over each root found.
            for vss = v_sols'
                %Compute the values of u and w at steady state.
                uss = -r2*(b*vss - 1)*(g2 + vss)/p2;
                wss = (p3*uss*vss + g3*s3 + s3*vss)/((g3 + vss)*mu3);

                %Check if there is an error in feasibility analysis.
                if(uss < -1e-14 || wss < -1e-14)
                    "Warning: error in feasibility analysis!"
                    M(i,j) = -1;
                    s1
                    s3
                else
                    %Compute if any instability occurs - if so, leave value 
                    %of M_Ij unchanged 
                    if(max(real(eigs(J(uss,vss,wss))))>0)
                        M(i,j) = M(i,j) + 0;

                    %If no instability occurs, then this steady state is
                    %stable so add 1 to M_ij
                    else
                        M(i,j)=M(i,j) + 1;
                    end
                end
            end
        end
    end
end

%Bad plotting code
close all;
figure;

%Number of feasible steady states - x axis is s1R, y axis is s3R
imagesc(flipud(Mn'))
colorbar;
xlabel('$s_1$','interpreter','latex')
ylabel('$s_1$','interpreter','latex')
ax = gca; set(ax,'Fontsize',16);
ax.XTick = [1, ax.XTick];
ax.XTickLabel = s1R(ax.XTick);
ax.YTick = [1, ax.YTick];
ax.YTickLabel = flip(s3R((ax.YTick)));
title('# of feasible steady states')


figure

%Number of feasible+stable steady states - x axis is s1R, y axis is s3R
imagesc(flipud(M'))
colorbar;
xlabel('$s_1$','interpreter','latex')
ylabel('$s_1$','interpreter','latex')
ax = gca; set(ax,'Fontsize',16);
ax.XTick = [1, ax.XTick];
ax.XTickLabel = s1R(ax.XTick);
ax.YTick = [1, ax.YTick];
ax.YTickLabel = flip(s3R((ax.YTick)));
title('# of feasible+stable steady states')


