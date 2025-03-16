figure

if(dims==1)
    plot(x,U(end,uN),'linewidth',2); hold on;
    plot(x,U(end,vN),'linewidth',2);
    plot(x,U(end,wN),'--','linewidth',2);

    legend('$u$','$v$','$w$','interpreter','latex')
else
    figure
    imagesc(reshape(U(end/2,vN),m,m));colorbar;title('$u$','interpreter','latex')
    figure
    imagesc(reshape(U(end,vN),m,m));colorbar;title('$u$','interpreter','latex')
end